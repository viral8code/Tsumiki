using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Model.Scaffolding;
using Tsumiki.Utility;

namespace Tsumiki.Core.Scaffolding
{
    /// <summary>
    /// GapFiller が埋められなかったスキャフォールドのギャップを、局所アセンブリ
    /// (MEGAHIT/IDBA の localasm 型)で埋める。AssemblyMerger(-mg)の安全な代替。
    ///
    /// -mg は他の k のアセンブリ結果(=既に確定した結論)を持ち込むため、
    /// 同じリードから作った別 k のアセンブリが同じ反復配列で同じ誤りをする
    /// リスクを抱える(統合は誤りを打ち消さず、両方の誤りを取り込む)。
    ///
    /// 局所アセンブリはそれと違い、ギャップの両端に実際にマップされたリードだけを
    /// 集め、そのリードだけからその場でミニアセンブリを組む。これは
    /// 「グローバルなグラフでは低カバレッジに埋もれて削られてしまった領域を、
    /// その領域だけのリードに限定して相対的に高いカバレッジとして扱い直す」
    /// (IDBA-UD の局所カバレッジ閾値と同じ思想)ことで、新しい証拠を持ち込む。
    /// GapFiller が使う信頼できる k-mer 集合はグローバルなカットオフを既に
    /// 適用済みだが、ここではローカルに集めたリードに対してカットオフ1
    /// (=1回でも読まれていれば信頼する)で再構築する。
    /// </summary>
    internal static class LocalAssembler
    {
        /// <summary>
        /// ギャップの両端からアンカーとして使う長さ。
        /// </summary>
        private const int アンカー長 = 300;

        /// <summary>
        /// 1ギャップに集める局所リードの上限。アンカーが反復配列と重なると
        /// 際限なくリードが集まりうるため、暴走を防ぐ。
        /// </summary>
        private const int 局所リード数の上限 = 4000;

        /// <summary>
        /// 局所アセンブリで使う k-mer カットオフ。1回読まれていれば信頼する。
        /// </summary>
        private const ulong 局所カットオフ = 1;

        public static 局所アセンブリ統計 V_充填_ギャップ(
            string p_スキャフォールドパス,
            string p_リード1のパス,
            string p_リード2のパス,
            int p_k長,
            string p_作業ディレクトリ)
        {
            var l_スキャフォールド群 = FastaReader.Get_全エントリ(p_スキャフォールドパス);

            var l_ギャップ一覧 = Get_対象ギャップ一覧(l_スキャフォールド群, p_k長);
            if (l_ギャップ一覧.Count == 0)
            {
                return new 局所アセンブリ統計(0, 0, 0, 0, 0, 0);
            }

            var l_アンカー索引 = Get_アンカー索引(l_ギャップ一覧, p_k長);
            var l_局所リード = Get_局所リード(l_アンカー索引, p_リード1のパス, p_リード2のパス, p_k長, l_ギャップ一覧.Count);

            var l_結果 = new string?[l_ギャップ一覧.Count];
            var l_埋めた数 = 0;
            var l_埋めた塩基数 = 0;
            var l_リード無し数 = 0;
            var l_一意でない数 = 0;
            var l_到達不能数 = 0;

            for (var g = 0; g < l_ギャップ一覧.Count; g++)
            {
                if (l_局所リード[g].Count == 0)
                {
                    l_リード無し数++;
                    continue;
                }

                var l_埋め = Get_局所アセンブリ結果(
                    l_ギャップ一覧[g], l_局所リード[g], p_k長, p_作業ディレクトリ, out var l_判定);
                if (l_埋め != null)
                {
                    l_結果[g] = l_埋め;
                    l_埋めた数++;
                    l_埋めた塩基数 += l_埋め.Length;
                }
                else if (l_判定 == ギャップ充填判定.一意でない)
                {
                    l_一意でない数++;
                }
                else
                {
                    l_到達不能数++;
                }
            }

            V_書き戻し(p_スキャフォールドパス, l_スキャフォールド群, l_ギャップ一覧, l_結果);

            return new 局所アセンブリ統計(
                l_ギャップ一覧.Count, l_埋めた数, l_埋めた塩基数, l_リード無し数, l_一意でない数, l_到達不能数);
        }

        private static List<局所ギャップ> Get_対象ギャップ一覧(
            List<(string A_ID, string A_配列)> p_スキャフォールド群, int p_k長)
        {
            List<局所ギャップ> l_結果 = [];
            for (var s = 0; s < p_スキャフォールド群.Count; s++)
            {
                var l_配列 = p_スキャフォールド群[s].A_配列;
                var i = 0;
                while (i < l_配列.Length)
                {
                    if (l_配列[i] != 'N')
                    {
                        i++;
                        continue;
                    }
                    var l_開始 = i;
                    while (i < l_配列.Length && l_配列[i] == 'N')
                    {
                        i++;
                    }
                    var l_長さ = i - l_開始;
                    // 両端に最低 k 長ぶんの足場が要る(左右のアンカーk-merを取るため)。
                    if (l_長さ <= Consts.ギャップ充填のギャップ長上限 && l_開始 >= p_k長 && i + p_k長 <= l_配列.Length)
                    {
                        var l_左長 = Math.Min(アンカー長, l_開始);
                        var l_右長 = Math.Min(アンカー長, l_配列.Length - i);
                        l_結果.Add(new 局所ギャップ(
                            s, l_開始, l_長さ,
                            l_配列.Substring(l_開始 - l_左長, l_左長),
                            l_配列.Substring(i, l_右長)));
                    }
                }
            }
            return l_結果;
        }

        /// <summary>
        /// 全ギャップの左右アンカーから k-mer 索引を作る。
        /// キーは正規形(順鎖・逆鎖どちらでも同じキーに寄る)。
        /// </summary>
        private static Dictionary<KmerKey, List<int>> Get_アンカー索引(List<局所ギャップ> p_ギャップ一覧, int p_k長)
        {
            Dictionary<KmerKey, List<int>> l_索引 = [];
            for (var g = 0; g < p_ギャップ一覧.Count; g++)
            {
                V_登録_kmer列(l_索引, p_ギャップ一覧[g].A_左アンカー, p_k長, g);
                V_登録_kmer列(l_索引, p_ギャップ一覧[g].A_右アンカー, p_k長, g);
            }
            return l_索引;
        }

        private static void V_登録_kmer列(Dictionary<KmerKey, List<int>> p_索引, string p_配列, int p_k長, int p_ギャップ番号)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (Get_曖昧塩基を含むか(p_配列, i, p_k長))
                {
                    continue;
                }
                var l_鍵 = new KmerKey(p_配列.AsSpan(i, p_k長)).Get_正規形();
                if (!p_索引.TryGetValue(l_鍵, out var l_一覧))
                {
                    l_一覧 = [];
                    p_索引[l_鍵] = l_一覧;
                }
                if (!l_一覧.Contains(p_ギャップ番号))
                {
                    l_一覧.Add(p_ギャップ番号);
                }
            }
        }

        private static bool Get_曖昧塩基を含むか(string p_配列, int p_開始, int p_長さ)
        {
            for (var j = 0; j < p_長さ; j++)
            {
                if (Util.Get_曖昧塩基か(p_配列[p_開始 + j]))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// 生リードを1回走査し、アンカーに触れたリードをギャップごとに集める。
        /// アンカーに触れたペアは両方のリードを局所リード集合に入れる
        /// (相方が、まだ組み込まれていない領域を読んでいる可能性があるため)。
        /// </summary>
        private static List<string>[] Get_局所リード(
            Dictionary<KmerKey, List<int>> p_アンカー索引,
            string p_リード1のパス, string p_リード2のパス, int p_k長, int p_ギャップ数)
        {
            var l_局所リード = new List<string>[p_ギャップ数];
            for (var g = 0; g < p_ギャップ数; g++)
            {
                l_局所リード[g] = [];
            }

            foreach (var l_パス in new[] { p_リード1のパス, p_リード2のパス })
            {
                if (string.IsNullOrWhiteSpace(l_パス) || !File.Exists(l_パス))
                {
                    continue;
                }
                using var l_読み込み = new FastqReader(l_パス);
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_リード = l_読み込み.Get_次のリード().A_生リード;
                    if (l_リード is null || l_リード.Length < p_k長)
                    {
                        continue;
                    }
                    foreach (var g in Get_一致するギャップ(p_アンカー索引, l_リード, p_k長))
                    {
                        if (l_局所リード[g].Count < 局所リード数の上限)
                        {
                            l_局所リード[g].Add(l_リード);
                        }
                    }
                }
            }
            return l_局所リード;
        }

        private static HashSet<int> Get_一致するギャップ(
            Dictionary<KmerKey, List<int>> p_索引, string p_リード, int p_k長)
        {
            HashSet<int>? l_見つかった = null;
            for (var i = 0; i + p_k長 <= p_リード.Length; i++)
            {
                if (Get_曖昧塩基を含むか(p_リード, i, p_k長))
                {
                    continue;
                }
                var l_鍵 = new KmerKey(p_リード.AsSpan(i, p_k長)).Get_正規形();
                if (p_索引.TryGetValue(l_鍵, out var l_一覧))
                {
                    l_見つかった ??= [];
                    foreach (var g in l_一覧)
                    {
                        _ = l_見つかった.Add(g);
                    }
                }
            }
            return l_見つかった ?? [];
        }

        /// <summary>
        /// 1ギャップぶんのミニアセンブリ。左右アンカー配列+局所リードだけから
        /// 使い捨ての TrustedKmerIndex を作り、GapFiller と同じ制約付き探索で
        /// 左アンカー末尾から右アンカー先頭までの経路を探す。
        /// </summary>
        private static string? Get_局所アセンブリ結果(
            局所ギャップ p_ギャップ, List<string> p_局所リード, int p_k長, string p_作業ディレクトリ,
            out ギャップ充填判定 p_判定)
        {
            if (p_ギャップ.A_左アンカー.Length < p_k長 || p_ギャップ.A_右アンカー.Length < p_k長)
            {
                p_判定 = ギャップ充填判定.到達不能;
                return null;
            }

            var l_一時ディレクトリ = Path.Combine(p_作業ディレクトリ, $"localasm_{Guid.NewGuid():N}");
            _ = Directory.CreateDirectory(l_一時ディレクトリ);
            try
            {
                // 索引の構築はギャップの数だけ繰り返される。1件ごとの
                // 統計はログを埋めるだけなので、この区間は記録を止める。
                using var l_休止 = Logger.V_止める_記録();
                using var l_索引 = new TrustedKmerIndex(l_一時ディレクトリ);
                V_登録_全kmer(l_索引, p_ギャップ.A_左アンカー, p_k長);
                V_登録_全kmer(l_索引, p_ギャップ.A_右アンカー, p_k長);
                foreach (var l_リード in p_局所リード)
                {
                    V_登録_全kmer(l_索引, l_リード, p_k長);
                }
                _ = l_索引.V_カットオフ(局所カットオフ);

                var l_左のkmer = Get_kmerバイト列(p_ギャップ.A_左アンカー, p_ギャップ.A_左アンカー.Length - p_k長, p_k長);
                var l_目標kmer = Get_kmerバイト列(p_ギャップ.A_右アンカー, 0, p_k長);
                if (l_左のkmer is null || l_目標kmer is null)
                {
                    p_判定 = ギャップ充填判定.到達不能;
                    return null;
                }

                var l_最小長 = Math.Max(0, p_ギャップ.A_長さ - Consts.ギャップ充填の長さの余裕幅);
                var l_最大長 = p_ギャップ.A_長さ + Consts.ギャップ充填の長さの余裕幅;
                (var l_経路, p_判定) = ConstrainedPathFinder.Get_経路(
                    l_左のkmer, l_目標kmer, l_最小長, l_最大長, l_索引, p_k長);
                return l_経路;
            }
            finally
            {
                try
                {
                    Directory.Delete(l_一時ディレクトリ, recursive: true);
                }
                catch (IOException)
                {
                    // 一時ファイルの掃除に失敗してもアセンブリ結果には影響しない。
                }
            }
        }

        private static void V_登録_全kmer(TrustedKmerIndex p_索引, string p_配列, int p_k長)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (Get_曖昧塩基を含むか(p_配列, i, p_k長))
                {
                    continue;
                }
                var l_kmer = Get_kmerバイト列(p_配列, i, p_k長);
                if (l_kmer is not null)
                {
                    p_索引.V_登録(l_kmer, p_ワーカー番号: 0);
                }
            }
        }

        private static byte[]? Get_kmerバイト列(string p_配列, int p_開始, int p_k長)
        {
            if (p_開始 < 0 || p_開始 + p_k長 > p_配列.Length)
            {
                return null;
            }
            var l_kmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[p_開始 + i]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    return null;
                }
                l_kmer[i] = l_塩基ID;
            }
            return l_kmer;
        }

        private static void V_書き戻し(
            string p_スキャフォールドパス,
            List<(string A_ID, string A_配列)> p_スキャフォールド群,
            List<局所ギャップ> p_ギャップ一覧,
            string?[] p_結果)
        {
            var l_scaffold別ギャップ = p_ギャップ一覧
                .Select((l_ギャップ, l_番号) => (l_ギャップ, l_番号))
                .GroupBy(x => x.l_ギャップ.A_Scaffold番号)
                .ToDictionary(x => x.Key, x => x.ToList());

            using var l_書き込み = new FastaWriter(p_スキャフォールドパス);
            for (var s = 0; s < p_スキャフォールド群.Count; s++)
            {
                var (l_ID, l_配列) = p_スキャフォールド群[s];
                if (!l_scaffold別ギャップ.TryGetValue(s, out var l_該当))
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                    continue;
                }

                var l_出力 = new StringBuilder();
                var l_直前終端 = 0;
                foreach (var (l_ギャップ, l_番号) in l_該当)
                {
                    _ = l_出力.Append(l_配列, l_直前終端, l_ギャップ.A_開始 - l_直前終端);
                    var l_埋め = p_結果[l_番号];
                    _ = l_埋め != null ? l_出力.Append(l_埋め) : l_出力.Append('N', l_ギャップ.A_長さ);
                    l_直前終端 = l_ギャップ.A_開始 + l_ギャップ.A_長さ;
                }
                _ = l_出力.Append(l_配列, l_直前終端, l_配列.Length - l_直前終端);
                l_書き込み.V_書き込み(l_ID, l_出力.ToString());
            }
        }

        public static void V_出力_統計(局所アセンブリ統計 p_統計)
        {
            if (p_統計.A_対象ギャップ数 == 0)
            {
                Logger.V_出力(メッセージID.局所アセンブリ_対象なし);
                return;
            }
            Logger.V_出力(
                メッセージID.局所アセンブリ統計, p_統計.A_埋めたギャップ数, p_統計.A_対象ギャップ数, p_統計.A_埋めた塩基数, p_統計.A_局所リードが集まらなかった数, p_統計.A_一意に定まらなかった数, p_統計.A_到達できなかった数);
        }

        private readonly record struct 局所ギャップ(
            int A_Scaffold番号, int A_開始, int A_長さ, string A_左アンカー, string A_右アンカー);
    }
}
