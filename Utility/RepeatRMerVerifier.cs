using Tsumiki.Common;
using Tsumiki.IO;

namespace Tsumiki.Utility
{
    /// <summary>
    /// 短い反復配列の経路付け替えに対する拒否権(ABySS RResolver 型)。
    ///
    /// UnitigGraph.V_解決_短い反復 の対応付けは、集計されたペア支持の多数決で
    /// 決めており、その接続を実際に跨いだリードが存在するかを直接は
    /// 確かめていない。生リードから作った厳密な r-mer 集合を使い、
    /// head→repeat→tail という経路の接合点を実際に跨いだリードが
    /// 存在するかを検定する。
    ///
    /// 接合点を跨がない(=どちらか一方の配列の内部だけに収まる)窓は、
    /// 経路の正しさに関わらず必ず真になる(そのunitig自身が既に実在の
    /// 配列である以上、内部のk-merは常に読まれている)ため、判定対象から除く。
    ///
    /// 重要: head/repeat/tail は de Bruijn グラフの辺である以上、隣接する
    /// もの同士は必ずアセンブリの k-1 塩基を共有しており、この共有区間は
    /// repeat 自身の配列の両端にもそのまま現れる。r がこの k 以下だと、
    /// 「接合点を跨ぐ」と判定した窓も共有区間の内側に収まってしまい、
    /// head だけ・repeat だけを独立に読んだリードでも真になる
    /// (=対応付けの正しさを何も検定できない、常に素通りする拒否権になる)。
    /// 呼び出し側は r をこの反復解決で使う k より確実に長く取る責任を持つ
    /// (AssemblyPipeline は k + rMer長のk超過分の既定値 で決めている)。
    ///
    /// r-mer は存在確認だけに使うため、TrustedKmerIndex のような外部ソート付き
    /// カウンタは不要で、単純な HashSet で足りる。r 長は 32 塩基以下
    /// (ulong に 2bit パック可能な上限)。
    /// </summary>
    internal sealed class RepeatRMerVerifier
    {
        private readonly HashSet<ulong> _rMer集合;
        private readonly int _r長;

        private RepeatRMerVerifier(HashSet<ulong> p_rMer集合, int p_r長)
        {
            this._rMer集合 = p_rMer集合;
            this._r長 = p_r長;
        }

        /// <summary>
        /// 生リードファイル群を1回走査し、出現した r-mer(正準形)の集合を作る。
        /// 空・存在しないパスは無視する(片側リードのみの実行に対応するため)。
        /// </summary>
        public static RepeatRMerVerifier V_構築(IEnumerable<string> p_リードパス一覧, int p_r長)
        {
            if (p_r長 is <= 0 or > 32)
            {
                throw new ArgumentException("r-mer length must be between 1 and 32 to pack into a ulong");
            }

            var l_集合 = new HashSet<ulong>();
            foreach (var l_パス in p_リードパス一覧)
            {
                if (string.IsNullOrWhiteSpace(l_パス) || !File.Exists(l_パス))
                {
                    continue;
                }
                using var l_読み込み = new FastqReader(l_パス);
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_リード = l_読み込み.Get_次のリード().A_生リード;
                    if (l_リード is not null)
                    {
                        V_登録_rMer(l_集合, l_リード, p_r長);
                    }
                }
            }
            return new RepeatRMerVerifier(l_集合, p_r長);
        }

        private static void V_登録_rMer(HashSet<ulong> p_集合, string p_リード, int p_r長)
        {
            for (var i = 0; i + p_r長 <= p_リード.Length; i++)
            {
                var l_曖昧か = false;
                for (var j = 0; j < p_r長; j++)
                {
                    if (Util.Get_曖昧塩基か(p_リード[i + j]))
                    {
                        l_曖昧か = true;
                        break;
                    }
                }
                if (!l_曖昧か)
                {
                    _ = p_集合.Add(Get_正準値(p_リード.AsSpan(i, p_r長)));
                }
            }
        }

        /// <summary>
        /// head→repeat→tail の経路上で、head-repeat 接合点と repeat-tail 接合点の
        /// 両方を実際に跨いだ r-mer のうち、生リード由来の集合に見つかった本数。
        /// </summary>
        public int Get_接合点の支持数(string p_head配列, string p_repeat配列, string p_tail配列)
        {
            // head-repeat、repeat-tail は de Bruijn グラフの辺である以上、
            // 必ず(このアセンブリの)k-1塩基を共有しており、そのコピーは
            // repeat 自身の配列の両端にもそのまま現れる。head・tail を
            // そのまま margin に使うと、接合点に重なり区間が二重に並ぶ
            // (head側のコピーの直後に repeat 自身のコピー)テスト配列を
            // 作ってしまい、本物のゲノム配列(重なりは一度しか現れない)には
            // 存在しない配列になる。head・tail からはこの重なりを除いた
            // 固有部分だけを margin に使う。
            var l_重なり長 = Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1);
            var l_head固有 = p_head配列.Length > l_重なり長 ? p_head配列[..^l_重なり長] : string.Empty;
            var l_tail固有 = p_tail配列.Length > l_重なり長 ? p_tail配列[l_重なり長..] : string.Empty;

            var l_margin長 = this._r長 - 1;
            var l_head側 = l_head固有.Length <= l_margin長 ? l_head固有 : l_head固有[^l_margin長..];
            var l_tail側 = l_tail固有.Length <= l_margin長 ? l_tail固有 : l_tail固有[..l_margin長];

            var l_テスト配列 = l_head側 + p_repeat配列 + l_tail側;
            var l_接合点1 = l_head側.Length;
            var l_接合点2 = l_head側.Length + p_repeat配列.Length;

            var l_支持数 = 0;
            for (var i = 0; i + this._r長 <= l_テスト配列.Length; i++)
            {
                var l_窓終端 = i + this._r長; // exclusive
                var l_接合点1を跨ぐ = i < l_接合点1 && l_接合点1 < l_窓終端;
                var l_接合点2を跨ぐ = i < l_接合点2 && l_接合点2 < l_窓終端;
                if (!l_接合点1を跨ぐ && !l_接合点2を跨ぐ)
                {
                    // 接合点を跨がない窓はどちらの経路でも必ず真になるため無視する。
                    continue;
                }

                var l_曖昧か = false;
                for (var j = 0; j < this._r長; j++)
                {
                    if (Util.Get_曖昧塩基か(l_テスト配列[i + j]))
                    {
                        l_曖昧か = true;
                        break;
                    }
                }
                if (!l_曖昧か && this._rMer集合.Contains(Get_正準値(l_テスト配列.AsSpan(i, this._r長))))
                {
                    l_支持数++;
                }
            }
            return l_支持数;
        }

        public bool Get_接合点に支持があるか(string p_head配列, string p_repeat配列, string p_tail配列, int p_閾値)
        {
            return this.Get_接合点の支持数(p_head配列, p_repeat配列, p_tail配列) >= p_閾値;
        }

        /// <summary>
        /// 配列とその逆相補のうち、2bitパック値が小さいほうを返す
        /// (順鎖/逆鎖どちらから読んでも同一のキーに正規化するため)。
        /// KmerKey は現在の実行時引数のk長を前提にするため r-mer(k とは
        /// 別の長さ)には使えず、KmerPacking の長さ非依存な API を使う。
        /// r 長(最大32)は UInt128 の正規化パック値がそのまま ulong に収まる。
        /// </summary>
        private static ulong Get_正準値(ReadOnlySpan<char> p_配列)
        {
            Span<byte> l_塩基ID列 = stackalloc byte[p_配列.Length];
            for (var i = 0; i < p_配列.Length; i++)
            {
                l_塩基ID列[i] = Util.Get_塩基ID(p_配列[i]);
            }
            return (ulong)KmerPacking.Get_正規化パック(l_塩基ID列);
        }
    }
}
