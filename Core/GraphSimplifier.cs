using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// de Bruijnグラフ(TrustedKmerIndexが保持する厳密なk-mer集合)に対する
    /// グラフ簡略化。2種類のアーティファクト除去を行う:
    ///
    /// 1. tip clipping: 短いdead-end unitig(片端の次数が0)を丸ごと除去する。
    ///    dead-endは定義上どこにも合流しないため、丸ごと除去しても他の
    ///    経路が使う配列を壊す心配がない。
    ///
    /// 2. 低カバレッジ端のトリミング: 真のbubble popping(2経路の合流点を
    ///    厳密に特定し、高頻度側だけを残す)ではなく簡易版。合流点を
    ///    明示的には特定せず、代わりに各unitigの両端から「カバレッジが
    ///    baseline(ゲノム全体の典型的な単一コピー相当の深度)に比べて
    ///    著しく低いk-merが続く間」だけを1つずつ剥がしていく。
    ///
    ///    重要な設計上の注意: 当初はunitig全体の平均カバレッジで判定
    ///    しようとしたが、SNP様の短い分岐(数塩基だけ異なりすぐ長い
    ///    共有配列に合流する)では、共有部分(高カバレッジ)に平均が
    ///    引きずられて低カバレッジ側を検出できないことがテストで判明した。
    ///    さらに、検出できてもunitig全体を丸ごと除去すると、合流後の
    ///    共有配列(高カバレッジ側の経路も使っている)まで消してしまい
    ///    別の経路を破壊してしまう。そのため「端から1kmerずつ、カバレッジが
    ///    baseline比で低い間だけ剥がす」方式にした。エラー由来の分岐は
    ///    合流点までの区間(=k-1個のk-mer)が本来低カバレッジのはずなので、
    ///    この区間だけを正確に剥がし、合流後の共有配列には手を付けない。
    /// </summary>
    internal static class GraphSimplifier
    {
        /// <summary>
        /// tip clippingと低カバレッジ端のトリミングを反復的に行い、
        /// 簡略化後の「unitig開始点」リストを返す。除去のたびにunitigを
        /// ゼロから再構築して次数・カバレッジを再評価するため、最大
        /// maxIterations回まで反復する(それ以上変化がなくなった時点で早期終了)。
        ///
        /// 既定の長さ閾値(k*10)は、合成データ(300kbゲノム・1%エラー率・
        /// エラー訂正後)での実測に基づく: k*2(Velvet等でよく使われる値)
        /// では訂正しきれず残った少数のエラー由来の分岐がまだ長めの
        /// tipとして残ってしまい、k*10まで広げてようやく大部分を吸収できた
        /// (収束までに約13反復かかったため既定のmaxIterationsも余裕を見て
        /// 30とした)。この長さ閾値は tip clipping(dead-endの丸ごと除去)
        /// にのみ適用し、低カバレッジ端のトリミングは(合流先の長さに
        /// 依存せず正しく判定できるため)unitigの長さによらず全件に適用する。
        ///
        /// 既定のカバレッジ閾値(baselineの20%)は、SPAdes等が「誤った接続
        /// (erroneous connection)」除去に使う値と同程度の、一般的に
        /// 保守的とされる水準。baselineは全unitigの平均カバレッジの
        /// 長さ加重中央値(短い断片が多数を占めても、ゲノムの大部分を
        /// 占める正しい主経路のカバレッジに引きずられにくくするため)。
        /// </summary>
        public static List<byte[]> ClipTips(
            TrustedKmerIndex index,
            int kmerLength,
            int? tipLengthThreshold = null,
            int maxIterations = 30,
            double lowCoverageFraction = 0.2)
        {
            var threshold = tipLengthThreshold ?? (10 * kmerLength);
            var unitigMaker = new UnitigMaker(index);
            var firstKmers = index.FindFirstKmers();

            for (var iteration = 1; iteration <= maxIterations; iteration++)
            {
                var unitigs = BuildUnitigs(unitigMaker, firstKmers);
                var baseline = WeightedMedianCoverage(index, unitigs, kmerLength);
                var lowCoverageCutoff = baseline * lowCoverageFraction;

                var tipsRemoved = 0;
                var trimmedKmerCount = 0;
                var unitigsTrimmed = 0;
                foreach (var unitig in unitigs)
                {
                    var seqBytes = ToNucleotideIdBytes(unitig);
                    if (seqBytes.Length < kmerLength)
                    {
                        continue;
                    }

                    if (unitig.Length < threshold)
                    {
                        var startDegree = index.CountInEdges(seqBytes.AsSpan(0, kmerLength));
                        var endDegree = index.CountOutEdges(seqBytes.AsSpan(seqBytes.Length - kmerLength, kmerLength));

                        // 片方の端が dead-end(そちら向きに続きがない)であれば tip とみなし、
                        // どこにも合流しないため丸ごと除去してよい。
                        // 両端とも dead-end の場合(=孤立した短い断片)も対象に含む。
                        if (startDegree == 0 || endDegree == 0)
                        {
                            RemoveUnitigKmers(index, seqBytes, kmerLength);
                            tipsRemoved++;
                            continue;
                        }
                    }

                    if (baseline <= 0)
                    {
                        continue;
                    }

                    var trimmed = TrimLowCoverageEdges(index, seqBytes, kmerLength, lowCoverageCutoff);
                    if (trimmed > 0)
                    {
                        trimmedKmerCount += trimmed;
                        unitigsTrimmed++;
                    }
                }

                Console.WriteLine($"[GraphSimplifier] Iteration {iteration}: examined {unitigs.Count} unitig(s) " +
                    $"(tip threshold < {threshold}bp, coverage baseline {baseline:0.#}), " +
                    $"removed {tipsRemoved} tip(s), trimmed {trimmedKmerCount} low-coverage k-mer(s) from {unitigsTrimmed} unitig edge(s).");

                if (tipsRemoved == 0 && trimmedKmerCount == 0)
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

        /// <summary>
        /// unitigの両端から、カバレッジがcutoff未満のk-merが続く間だけ
        /// 1つずつ信頼できる集合から除去する(合流点に達して初めて
        /// cutoff以上のk-merが現れたら、そこで止めて先には進まない)。
        /// 先頭側と末尾側で除去範囲が重ならないよう、互いの残り長で制限する。
        /// 戻り値: 除去したk-merの総数。
        /// </summary>
        private static int TrimLowCoverageEdges(TrustedKmerIndex index, byte[] seqBytes, int kmerLength, double cutoff)
        {
            var numKmers = seqBytes.Length - kmerLength + 1;
            if (numKmers <= 0)
            {
                return 0;
            }

            var removed = 0;

            var fromStart = 0;
            while (fromStart < numKmers && index.GetCoverage(seqBytes.AsSpan(fromStart, kmerLength)) < cutoff)
            {
                fromStart++;
            }

            var fromEnd = 0;
            while (fromEnd < numKmers - fromStart && index.GetCoverage(seqBytes.AsSpan(numKmers - 1 - fromEnd, kmerLength)) < cutoff)
            {
                fromEnd++;
            }

            for (var i = 0; i < fromStart; i++)
            {
                index.RemoveTrusted(seqBytes.AsSpan(i, kmerLength));
                removed++;
            }
            for (var i = 0; i < fromEnd; i++)
            {
                index.RemoveTrusted(seqBytes.AsSpan(numKmers - 1 - i, kmerLength));
                removed++;
            }

            return removed;
        }

        /// <summary>
        /// unitigを構成する全k-merのカバレッジの単純平均。
        /// </summary>
        private static double AverageCoverage(TrustedKmerIndex index, byte[] seqBytes, int kmerLength)
        {
            ulong sum = 0;
            var count = 0;
            for (var i = 0; i + kmerLength <= seqBytes.Length; i++)
            {
                sum += index.GetCoverage(seqBytes.AsSpan(i, kmerLength));
                count++;
            }
            return count == 0 ? 0 : (double)sum / count;
        }

        /// <summary>
        /// 全unitigの平均カバレッジの長さ加重中央値。多数を占めうる短い
        /// 断片(エラー由来のtip/bubble候補そのもの)に引きずられず、
        /// ゲノムの大部分を占める正しい主経路のカバレッジ水準を推定するため、
        /// 単純平均・単純中央値ではなく塩基数で重み付けした中央値を使う。
        /// </summary>
        private static double WeightedMedianCoverage(TrustedKmerIndex index, List<string> unitigs, int kmerLength)
        {
            if (unitigs.Count == 0)
            {
                return 0;
            }

            var pairs = unitigs
                .Select(u => (Length: (long)u.Length, Coverage: AverageCoverage(index, ToNucleotideIdBytes(u), kmerLength)))
                .OrderBy(p => p.Coverage)
                .ToList();
            var totalLength = pairs.Sum(p => p.Length);
            if (totalLength == 0)
            {
                return 0;
            }

            var half = totalLength / 2.0;
            long cumulative = 0;
            foreach (var (length, coverage) in pairs)
            {
                cumulative += length;
                if (cumulative >= half)
                {
                    return coverage;
                }
            }
            return pairs[^1].Coverage;
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
