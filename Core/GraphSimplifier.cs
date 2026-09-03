using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// de Bruijnグラフ(TrustedKmerIndexが保持する厳密なk-mer集合)に対する
    /// グラフ簡略化。現時点では tip clipping(短いdead-end unitigの除去)のみを
    /// 実装する。bubble popping(似た長さの2経路が再収束する構造の解消)は、
    /// 高頻度側を残すためにk-merごとのカバレッジ情報を保持する必要があり
    /// (現状はカットオフ通過/不通過の2値のみ保持)、設計拡張が要るため
    /// 別フェーズへ持ち越す。
    ///
    /// tip clippingの考え方: エラー訂正後も残った少数の未訂正エラーは、
    /// 主経路から分岐する短いdead-end unitig(tip)としてグラフに現れる。
    /// tipを構成するk-merを信頼できる集合から丸ごと除去すると、それが
    /// 分岐点になっていた箇所の分岐が解消され(次数が2以上から1に戻り)、
    /// 本来1本につながるはずだった主経路が正しく1本のunitigとして
    /// 再構築されるようになる。
    /// </summary>
    internal static class GraphSimplifier
    {
        /// <summary>
        /// tipLengthThreshold未満の長さで、かつ少なくとも片端がdead-end
        /// (そちら向きの次数が0)であるunitigを反復的に除去し、簡略化後の
        /// 「unitig開始点」リストを返す。除去のたびにunitigをゼロから
        /// 再構築して次数を再評価するため、最大maxIterations回まで反復する
        /// (それ以上除去されるtipがなくなった時点で早期終了する)。
        ///
        /// 既定の閾値(k*10)は、合成データ(300kbゲノム・1%エラー率・
        /// エラー訂正後)での実測に基づく: k*2(Velvet等でよく使われる値)
        /// では訂正しきれず残った少数のエラー由来の分岐がまだ長めの
        /// tipとして残ってしまい、k*10まで広げてようやく大部分を吸収できた
        /// (収束までに約13反復かかったため既定のmaxIterationsも余裕を見て
        /// 30とした)。bubble popping(未実装、別フェーズ)が入るまでは
        /// tip clippingがエラー由来の断片化を吸収する主な手段になる。
        /// </summary>
        public static List<byte[]> ClipTips(
            TrustedKmerIndex index,
            int kmerLength,
            int? tipLengthThreshold = null,
            int maxIterations = 30)
        {
            var threshold = tipLengthThreshold ?? (10 * kmerLength);
            var unitigMaker = new UnitigMaker(index);
            var firstKmers = index.FindFirstKmers();

            for (var iteration = 1; iteration <= maxIterations; iteration++)
            {
                var unitigs = BuildUnitigs(unitigMaker, firstKmers);

                var tipsRemoved = 0;
                foreach (var unitig in unitigs)
                {
                    if (unitig.Length >= threshold)
                    {
                        continue;
                    }

                    var seqBytes = ToNucleotideIdBytes(unitig);
                    var startDegree = index.CountInEdges(seqBytes.AsSpan(0, kmerLength));
                    var endDegree = index.CountOutEdges(seqBytes.AsSpan(seqBytes.Length - kmerLength, kmerLength));

                    // 片方の端が dead-end(そちら向きに続きがない)であれば tip とみなす。
                    // 両端とも dead-end の場合(=グラフ全体から孤立した短い断片)も対象に含む。
                    if (startDegree == 0 || endDegree == 0)
                    {
                        RemoveUnitigKmers(index, seqBytes, kmerLength);
                        tipsRemoved++;
                    }
                }

                Console.WriteLine($"[TipClipping] Iteration {iteration}: examined {unitigs.Count} unitig(s) (< {threshold}bp threshold check), removed {tipsRemoved} tip(s).");

                if (tipsRemoved == 0)
                {
                    return firstKmers;
                }

                // k-mer集合が縮小されたため、開始点を再検出してから次の反復へ。
                firstKmers = index.FindFirstKmers();
            }

            return firstKmers;
        }

        private static List<string> BuildUnitigs(UnitigMaker unitigMaker, List<byte[]> firstKmers)
        {
            List<string> unitigs = [];
            HashSet<string> seen = [];
            foreach (var kmer in firstKmers)
            {
                var unitig = unitigMaker.MakeUnitig(kmer);
                if (seen.Contains(unitig.Sequence) || seen.Contains(Util.ReverseComprement(unitig.Sequence)))
                {
                    continue;
                }
                _ = seen.Add(unitig.Sequence);
                _ = seen.Add(Util.ReverseComprement(unitig.Sequence));
                unitigs.Add(unitig.Sequence);
            }
            return unitigs;
        }

        private static byte[] ToNucleotideIdBytes(string sequence)
        {
            var bytes = new byte[sequence.Length];
            for (var i = 0; i < sequence.Length; i++)
            {
                bytes[i] = Util.GetSimpleNucleotideID(sequence[i]);
            }
            return bytes;
        }

        private static void RemoveUnitigKmers(TrustedKmerIndex index, byte[] seqBytes, int kmerLength)
        {
            for (var i = 0; i + kmerLength <= seqBytes.Length; i++)
            {
                index.RemoveTrusted(seqBytes.AsSpan(i, kmerLength));
            }
        }
    }
}
