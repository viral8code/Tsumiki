using Tsumiki.Core.UnitigBuilding;

namespace Tsumiki.Core.Output
{
    /// <summary>
    /// unitig グラフを GFA1 形式で書き出す (SPAdes/Unicycler と同様の診断出力)
    /// </summary>
    /// <remarks>
    /// 決められない分岐は、現状では walk の打ち切り点になるだけで
    /// 「なぜそこで切れたか」の情報が contigs.fasta には残らない<br/>
    /// GFA としてグラフそのものを出力すれば、Bandage 等のビューアで
    /// 「あと何が解ければ閉じるのか」を直接見られる<br/>
    /// 完全長を目指す作業は
    /// 本質的に反復的であり、診断可能性そのものが機能である<br/>
    /// 出力するのは V_結合_コンティグ がバブル除去・反復解決を終えた後の
    /// グラフの状態 (=最終的な walk がどの分岐を残したまま打ち切られたかを
    /// 反映する)<br/>
    /// バブル除去前の生の de Bruijn グラフではない
    /// </remarks>
    internal static class GfaWriter
    {
        /// <summary>
        /// p_ユニティグ配列・p_グラフ の状態を GFA1 として p_パス へ書き出す
        /// </summary>
        /// <remarks>
        /// 頂点は unitig ID(1 始まり) の順鎖/逆鎖のペアで表現されているため、
        /// 各物理的な隣接は双子の辺として 2 回現れる<br/>
        /// 片方だけを 1 本の L 行として出す
        /// </remarks>
        public static void V_出力(
            string p_パス,
            List<string> p_ユニティグ配列,
            UnitigGraph p_グラフ,
            int p_k長,
            IReadOnlyDictionary<int, int>? p_コピー数 = null)
        {
            var l_重なり長 = Math.Max(0, p_k長 - 1);
            var l_ユニティグ数 = (p_ユニティグ配列.Count - 2) / 2;

            using var l_書き込み = new StreamWriter(p_パス);
            l_書き込み.WriteLine("H\tVN:Z:1.0");

            for (var l_ID = 1; l_ID <= l_ユニティグ数; l_ID++)
            {
                var l_配列 = p_ユニティグ配列[l_ID << 1]; // 順鎖側
                var l_深度タグ = p_コピー数 is not null && p_コピー数.TryGetValue(l_ID, out var l_コピー数値)
                    ? $"\tCN:i:{l_コピー数値}"
                    : string.Empty;
                l_書き込み.WriteLine($"S\t{l_ID}\t{l_配列}\tLN:i:{l_配列.Length}{l_深度タグ}");
            }

            // 同じ物理的隣接が v→w と w^1→v^1 の双子として 2 回現れるため、
            // 片方を出したらもう片方は出さない
            HashSet<(int, int)> l_出力済み = [];
            for (var v = 2; v < p_グラフ.A_出辺.Count; v++)
            {
                foreach (var w in p_グラフ.A_出辺[v])
                {
                    if (l_出力済み.Contains((w ^ 1, v ^ 1)))
                    {
                        continue;
                    }
                    _ = l_出力済み.Add((v, w));

                    var l_始点ID = v >> 1;
                    var l_始点向き = (v & 1) == 0 ? "+" : "-";
                    var l_終点ID = w >> 1;
                    var l_終点向き = (w & 1) == 0 ? "+" : "-";
                    l_書き込み.WriteLine($"L\t{l_始点ID}\t{l_始点向き}\t{l_終点ID}\t{l_終点向き}\t{l_重なり長}M");
                }
            }
        }
    }
}
