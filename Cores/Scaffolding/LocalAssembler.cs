using System.Text;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// GapFiller が埋められなかった scaffold のギャップを、局所アセンブリ (MEGAHIT/IDBA の localasm 型) で埋める
    /// </summary>
    internal static class LocalAssembler
    {
        #region 定数

        /// <summary>
        /// ギャップの両端からアンカーとして使う長さ
        /// </summary>
        private const int アンカー長 = 300;

        /// <summary>
        /// 1 ギャップに集める局所リードの上限
        /// </summary>
        /// <remarks>
        /// アンカーが反復配列と重なると際限なくリードが集まりうるため、暴走を防ぐ
        /// </remarks>
        private const int 局所リード数の上限 = 4_000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 残ったギャップを、その両端に付いたリードだけで組み直して埋める
        /// </summary>
        /// <param name="p_スキャフォールドパス">対象の scaffold のパス</param>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>局所アセンブリの集計</returns>
        public static 局所アセンブリ統計 V_充填_ギャップ(string p_スキャフォールドパス, string p_リード1のパス, string p_リード2のパス, int p_k長)
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

                var l_埋め = Get_局所アセンブリ結果(l_ギャップ一覧[g], l_局所リード[g], p_k長, out var l_判定);
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

            return new 局所アセンブリ統計(l_ギャップ一覧.Count, l_埋めた数, l_埋めた塩基数, l_リード無し数, l_一意でない数, l_到達不能数);
        }

        /// <summary>
        /// 局所アセンブリの結果をログへ出力する
        /// </summary>
        /// <param name="p_統計">局所アセンブリの結果</param>
        public static void V_出力_統計(局所アセンブリ統計 p_統計)
        {
            if (p_統計.A_対象ギャップ数 == 0)
            {
                Logger.V_出力(メッセージID.局所アセンブリ_対象なし);
                return;
            }
            Logger.V_出力(メッセージID.局所アセンブリ統計, p_統計.A_埋めたギャップ数, p_統計.A_対象ギャップ数, p_統計.A_埋めた塩基数, p_統計.A_局所リードが集まらなかった数, p_統計.A_一意に定まらなかった数, p_統計.A_到達できなかった数);
        }

        /// <summary>
        /// 指定した文脈長で局所グラフの一意な充填配列を求める
        /// </summary>
        /// <param name="p_ギャップ">充填対象</param>
        /// <param name="p_局所リード">局所に集めたリード</param>
        /// <param name="p_k長">文脈長</param>
        /// <param name="p_判定">探索の結果</param>
        /// <returns>一意な充填配列、未確定なら null</returns>
        internal static string? Get_固定kの局所結果(局所ギャップ p_ギャップ, List<string> p_局所リード, int p_k長, out ギャップ充填判定 p_判定)
        {
            if (p_ギャップ.A_左アンカー.Length < p_k長 || p_ギャップ.A_右アンカー.Length < p_k長)
            {
                p_判定 = ギャップ充填判定.到達不能;
                return null;
            }

            var l_集合 = new LocalKmerSet(p_k長);
            V_登録_全kmer(l_集合, p_ギャップ.A_左アンカー, p_k長);
            V_登録_全kmer(l_集合, p_ギャップ.A_右アンカー, p_k長);
            foreach (var l_リード in p_局所リード)
            {
                V_登録_全kmer(l_集合, l_リード, p_k長);
            }

            var l_左のkmer = Get_kmerバイト列(p_ギャップ.A_左アンカー, p_ギャップ.A_左アンカー.Length - p_k長, p_k長);
            var l_目標kmer = Get_kmerバイト列(p_ギャップ.A_右アンカー, 0, p_k長);
            if (l_左のkmer is null || l_目標kmer is null)
            {
                p_判定 = ギャップ充填判定.到達不能;
                return null;
            }

            var l_最小長 = Math.Max(0, p_ギャップ.A_長さ - Consts.ギャップ充填の長さの余裕幅);
            var l_最大長 = p_ギャップ.A_長さ + Consts.ギャップ充填の長さの余裕幅;
            (var l_経路, p_判定) = ConstrainedPathFinder.Get_経路(l_左のkmer, l_目標kmer, l_最小長, l_最大長, l_集合, p_k長);
            return l_経路;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 埋める対象になるギャップを集めて返す
        /// </summary>
        /// <param name="p_スキャフォールド群">対象の scaffold</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>対象のギャップ</returns>
        private static List<局所ギャップ> Get_対象ギャップ一覧(List<(string A_ID, string A_配列)> p_スキャフォールド群, int p_k長)
        {
            List<局所ギャップ> l_結果 = [];
            for (var s = 0; s < p_スキャフォールド群.Count; s++)
            {
                var l_配列 = p_スキャフォールド群[s].A_配列;
                var l_i = 0;
                while (l_i < l_配列.Length)
                {
                    if (l_配列[l_i] != 'N')
                    {
                        l_i++;
                        continue;
                    }
                    var l_開始 = l_i;
                    while (l_i < l_配列.Length && l_配列[l_i] == 'N')
                    {
                        l_i++;
                    }
                    var l_長さ = l_i - l_開始;

                    // 両端に最低 k 長ぶんの足場が要る (左右のアンカー k-mer を取るため)
                    if (l_長さ <= Consts.ギャップ充填のギャップ長上限 && l_開始 >= p_k長 && l_i + p_k長 <= l_配列.Length)
                    {
                        var l_左長 = Math.Min(アンカー長, l_開始);
                        var l_右長 = Math.Min(アンカー長, l_配列.Length - l_i);
                        l_結果.Add(new 局所ギャップ(s, l_開始, l_長さ, l_配列.Substring(l_開始 - l_左長, l_左長), l_配列.Substring(l_i, l_右長)));
                    }
                }
            }
            return l_結果;
        }

        /// <summary>
        /// 全ギャップの左右アンカーから k-mer 索引を作る
        /// </summary>
        /// <remarks>
        /// キーは正規形 (順鎖・逆鎖どちらでも同じキーに寄る)
        /// </remarks>
        /// <param name="p_ギャップ一覧"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
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

        /// <summary>
        /// 配列の k-mer を、どのギャップの端かと一緒に索引へ登録する
        /// </summary>
        /// <param name="p_索引">登録先の索引</param>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_ギャップ番号">ギャップの番号</param>
        private static void V_登録_kmer列(Dictionary<KmerKey, List<int>> p_索引, string p_配列, int p_k長, int p_ギャップ番号)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (Has曖昧塩基(p_配列, i, p_k長))
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

        /// <summary>
        /// 指定した範囲に曖昧塩基が含まれるか
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始">調べ始める位置</param>
        /// <param name="p_長さ">調べる長さ</param>
        /// <returns>含まれれば true</returns>
        private static bool Has曖昧塩基(string p_配列, int p_開始, int p_長さ)
        {
            for (var j = 0; j < p_長さ; j++)
            {
                if (Util.Is曖昧塩基(p_配列[p_開始 + j]))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// 生リードを 1 回走査し、アンカーに触れたリードをギャップごとに集める
        /// </summary>
        /// <remarks>
        /// アンカーに触れたペアは両方のリードを局所リード集合に入れる (相方が、まだ組み込まれていない領域を読んでいる可能性があるため)
        /// </remarks>
        /// <param name="p_アンカー索引"></param>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_ギャップ数"></param>
        /// <returns></returns>
        private static List<string>[] Get_局所リード(Dictionary<KmerKey, List<int>> p_アンカー索引, string p_リード1のパス, string p_リード2のパス, int p_k長, int p_ギャップ数)
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
                while (l_読み込み.Has続き())
                {
                    var l_リード = l_読み込み.Get_次のリード().A_生リード;
                    if (l_リード is null || l_リード.Length < p_k長)
                    {
                        continue;
                    }
                    foreach (var l_g in Get_一致するギャップ(p_アンカー索引, l_リード, p_k長))
                    {
                        if (l_局所リード[l_g].Count < 局所リード数の上限)
                        {
                            l_局所リード[l_g].Add(l_リード);
                        }
                    }
                }
            }
            return l_局所リード;
        }

        /// <summary>
        /// そのリードが端に当たるギャップの番号を返す
        /// </summary>
        /// <param name="p_索引">ギャップの端の k-mer 索引</param>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>当たったギャップの番号</returns>
        private static HashSet<int> Get_一致するギャップ(Dictionary<KmerKey, List<int>> p_索引, string p_リード, int p_k長)
        {
            HashSet<int>? l_見つかった = null;
            for (var i = 0; i + p_k長 <= p_リード.Length; i++)
            {
                if (Has曖昧塩基(p_リード, i, p_k長))
                {
                    continue;
                }

                var l_鍵 = new KmerKey(p_リード.AsSpan(i, p_k長)).Get_正規形();
                if (p_索引.TryGetValue(l_鍵, out var l_一覧))
                {
                    l_見つかった ??= [];
                    foreach (var l_g in l_一覧)
                    {
                        _ = l_見つかった.Add(l_g);
                    }
                }
            }
            return l_見つかった ?? [];
        }

        /// <summary>
        /// 1 ギャップぶんのミニアセンブリ
        /// </summary>
        /// <remarks>
        /// 左右アンカー配列+局所リードだけから使い捨ての LocalKmerSet を作り、GapFiller と同じ制約付き探索で左アンカー末尾から右アンカー先頭までの経路を探す<br/>
        /// 集めた k-mer は局所リード数の上限ぶんしかなくインメモリで完結するため、TrustedKmerIndex のようなディスク経由のシャード集計は使わない (ギャップの数だけ繰り返すには重すぎる)
        /// </remarks>
        /// <param name="p_ギャップ"></param>
        /// <param name="p_局所リード"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_判定"></param>
        /// <returns></returns>
        private static string? Get_局所アセンブリ結果(局所ギャップ p_ギャップ, List<string> p_局所リード, int p_k長, out ギャップ充填判定 p_判定)
        {
            var l_経路 = Get_固定kの局所結果(p_ギャップ, p_局所リード, p_k長, out p_判定);
            if (l_経路 is not null || p_判定 != ギャップ充填判定.一意でない)
            {
                return l_経路;
            }

            string? l_一致した経路 = null;
            var l_異なるリード = p_局所リード.Distinct(StringComparer.Ordinal).ToList();
            foreach (var l_追加長 in new[] { 10, 20 })
            {
                var l_局所k = p_k長 + l_追加長;
                if (l_異なるリード.Count(x => x.Length >= l_局所k + 1) < 2)
                {
                    continue;
                }

                var l_候補 = Get_固定kの局所結果(p_ギャップ, l_異なるリード, l_局所k, out var l_局所判定);
                if (l_局所判定 == ギャップ充填判定.一意でない)
                {
                    return null;
                }

                if (l_候補 is null)
                {
                    continue;
                }

                if (l_一致した経路 is not null && l_一致した経路 != l_候補)
                {
                    return null;
                }
                l_一致した経路 = l_候補;
            }
            if (l_一致した経路 is not null)
            {
                // 可変文脈は元のグラフで曖昧だった領域にだけ適用する
                // 長い窓が一度も観測されない配列をアンカーの合成で作らない
                var l_接続 = p_ギャップ.A_左アンカー[^p_k長..] + l_一致した経路 + p_ギャップ.A_右アンカー[..p_k長];
                for (var i = 0; i + p_k長 + 1 <= l_接続.Length; i++)
                {
                    var l_窓 = l_接続.Substring(i, p_k長 + 1);
                    var l_逆窓 = Util.V_逆相補(l_窓);
                    if (l_異なるリード.Count(x => x.Contains(l_窓, StringComparison.Ordinal) || x.Contains(l_逆窓, StringComparison.Ordinal)) < 2)
                    {
                        return null;
                    }
                }
                p_判定 = ギャップ充填判定.充填済み;
            }
            return l_一致した経路;
        }

        /// <summary>
        /// 配列のすべての k-mer を索引へ登録する
        /// </summary>
        /// <param name="p_索引">登録先の索引</param>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_k長">k 長</param>
        private static void V_登録_全kmer(LocalKmerSet p_索引, string p_配列, int p_k長)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (Has曖昧塩基(p_配列, i, p_k長))
                {
                    continue;
                }

                var l_kmer = Get_kmerバイト列(p_配列, i, p_k長);
                if (l_kmer is not null)
                {
                    p_索引.V_登録(l_kmer);
                }
            }
        }

        /// <summary>
        /// 指定した位置の k-mer を塩基 ID 列にして返す
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始">k-mer の開始位置</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>塩基 ID 列、曖昧塩基を含む場合は null</returns>
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

        /// <summary>
        /// 埋まったギャップを反映して scaffold を書き直す
        /// </summary>
        /// <param name="p_スキャフォールドパス">書き出し先</param>
        /// <param name="p_スキャフォールド群">対象の scaffold</param>
        /// <param name="p_ギャップ一覧">埋めたギャップ</param>
        /// <param name="p_結果"></param>
        private static void V_書き戻し(string p_スキャフォールドパス, List<(string A_ID, string A_配列)> p_スキャフォールド群, List<局所ギャップ> p_ギャップ一覧, string?[] p_結果)
        {
            var l_scaffold別ギャップ = p_ギャップ一覧
                .Select((l_ギャップ, l_番号) => (l_ギャップ, l_番号))
                .GroupBy(x => x.l_ギャップ.A_足場番号)
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

        #endregion
    }
}
