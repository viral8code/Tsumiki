using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// スキャフォールドの N で埋められたギャップを、de Bruijn グラフ上で
    /// 両側を繋ぐ経路を探して実際の塩基配列に置き換える。
    ///
    /// スキャフォールディングは「2つの contig がこの向きでこれくらい離れて
    /// 隣接している」というペアエンドの情報だけで繋いでおり、その間の配列が
    /// 何であるかは埋めていない。しかしギャップの中身は多くの場合 k-mer 集合の
    /// 中に存在している。contig が途切れたのは配列が無いからではなく、
    /// 「分岐があってどちらへ進むか決められなかった」ためであることが多いからで、
    /// その場合ギャップの両端を繋ぐ経路はグラフ上に実在する。
    ///
    /// 探索は左端の k-mer から1塩基ずつ伸ばす幅優先で行い、推定ギャップ長の
    /// 前後 margin の範囲で右端の k-mer に到達する経路を集める。
    /// 経路がちょうど1本に定まったときだけ埋める。複数見つかった場合は
    /// どれが正しいか決められないため N のまま残す(誤った配列で埋めるより、
    /// 分からないことが分かる状態のほうが下流の解析にとって安全)。
    /// </summary>
    internal static class GapFiller
    {
        /// <summary>
        /// 推定ギャップ長に対して許容する誤差(塩基)。インサートサイズ推定の
        /// ばらつきがそのままギャップ長推定のばらつきになるため、
        /// ぴったりの長さだけを探すと現実にはまず当たらない。
        /// </summary>
        private const int LengthMargin = 30;

        /// <summary>
        /// 1つのギャップあたりに展開してよい探索状態の上限。
        /// 分岐の多い領域では経路数が指数的に増えるため、上限を超えたら
        /// 「解けなかった」として諦める(時間をかけても曖昧なままのことが多い)。
        /// </summary>
        private const int MaxStatesPerGap = 200_000;

        /// <summary>
        /// これより長いギャップは探索空間が広すぎるうえ、推定長の誤差も大きく
        /// 一意に定まる見込みが薄いため対象外とする。
        /// </summary>
        private const int MaxGapLength = 500;

        internal readonly record struct Stats(int TotalGaps, int FilledGaps, int FilledBases, int AmbiguousGaps, int UnreachableGaps);

        /// <summary>
        /// scaffoldPath を読み込み、埋められるギャップを埋めて同じパスへ書き戻す。
        /// </summary>
        public static Stats Run(string scaffoldPath, TrustedKmerIndex index, int kmerLength)
        {
            List<(string Id, string Seq)> scaffolds = [];
            using (var reader = new FastaReader(scaffoldPath))
            {
                while (reader.HasNext())
                {
                    var seq = reader.NextSequence();
                    scaffolds.Add((seq.ID.TrimStart('>'), seq.Seq));
                }
            }

            var totalGaps = 0;
            var filledGaps = 0;
            var filledBases = 0;
            var ambiguous = 0;
            var unreachable = 0;

            List<(string Id, string Seq)> result = [];
            foreach (var (scaffoldId, sequence) in scaffolds)
            {
                var sb = new StringBuilder();
                var position = 0;
                while (position < sequence.Length)
                {
                    if (sequence[position] != 'N')
                    {
                        _ = sb.Append(sequence[position]);
                        position++;
                        continue;
                    }

                    // N の連続区間 = 1つのギャップ。
                    var gapStart = position;
                    while (position < sequence.Length && sequence[position] == 'N')
                    {
                        position++;
                    }
                    var gapLength = position - gapStart;
                    totalGaps++;

                    var filled = TryFill(sb, sequence, gapStart, gapLength, position, index, kmerLength, out var outcome);
                    if (filled != null)
                    {
                        _ = sb.Append(filled);
                        filledGaps++;
                        filledBases += filled.Length;
                    }
                    else
                    {
                        if (outcome == FillOutcome.Ambiguous)
                        {
                            ambiguous++;
                        }
                        else
                        {
                            unreachable++;
                        }
                        _ = sb.Append('N', gapLength);
                    }
                }
                result.Add((scaffoldId, sb.ToString()));
            }

            using (var writer = new FastaWriter(scaffoldPath))
            {
                foreach (var (scaffoldId, seq) in result)
                {
                    writer.Write(scaffoldId, seq);
                }
            }

            return new Stats(totalGaps, filledGaps, filledBases, ambiguous, unreachable);
        }

        private enum FillOutcome
        {
            Filled,
            Ambiguous,
            Unreachable,
        }

        /// <summary>
        /// ギャップの左右の足場から、その間を埋める配列を探す。
        /// 見つからない/一意に定まらない場合は null を返す。
        /// </summary>
        private static string? TryFill(
            StringBuilder left,
            string sequence,
            int gapStart,
            int gapLength,
            int gapEnd,
            TrustedKmerIndex index,
            int kmerLength,
            out FillOutcome outcome)
        {
            outcome = FillOutcome.Unreachable;

            if (gapLength > MaxGapLength || left.Length < kmerLength)
            {
                return null;
            }
            if (gapEnd + kmerLength > sequence.Length)
            {
                return null;
            }

            // 左側の足場: 既に書き出した配列の末尾 k-mer。
            var leftKmer = new byte[kmerLength];
            for (var i = 0; i < kmerLength; i++)
            {
                leftKmer[i] = Util.GetSimpleNucleotideID(left[left.Length - kmerLength + i]);
            }

            // 右側の足場: ギャップ直後の k-mer。ここへ到達できれば繋がったことになる。
            var targetKmer = new byte[kmerLength];
            for (var i = 0; i < kmerLength; i++)
            {
                targetKmer[i] = Util.GetSimpleNucleotideID(sequence[gapEnd + i]);
            }

            if (Array.IndexOf(leftKmer, Consts.InvalidBase) >= 0 || Array.IndexOf(targetKmer, Consts.InvalidBase) >= 0)
            {
                return null;
            }
            if (!index.Contains(leftKmer) || !index.Contains(targetKmer))
            {
                // 足場そのものが信頼できる k-mer 集合に無いなら探索しても意味がない。
                return null;
            }

            var minLength = Math.Max(0, gapLength - LengthMargin);
            var maxLength = gapLength + LengthMargin;

            // 幅優先で1塩基ずつ伸ばす。
            //
            // 各状態が「これまでに継ぎ足した塩基列」そのものを持つと、
            // 状態数の上限(20万)× 経路長ぶんのメモリと文字列コピーが発生する。
            // 代わりに親へのインデックスと追加した1塩基だけを持ち、
            // 解が見つかったときに親を辿って復元する。1状態あたり定数サイズで済む。
            var nodes = new List<(int Parent, byte Base)>(1024) { (-1, 0) };
            var kmers = new List<byte[]>(1024) { leftKmer };
            var depths = new List<int>(1024) { 0 };

            var found = new List<string>();
            var queue = new Queue<int>();
            queue.Enqueue(0);

            var buffer = new byte[kmerLength];

            while (queue.Count > 0)
            {
                var currentIndex = queue.Dequeue();
                var kmer = kmers[currentIndex];
                var added = depths[currentIndex];

                // added は「左の足場 k-mer の後ろに継ぎ足した塩基数」。目標 k-mer に
                // 到達した時点では、その末尾 kmerLength 塩基が目標 k-mer 自身に
                // あたる(それは元の配列に既にある)ので、ギャップを実際に埋める
                // 長さは added - kmerLength になる。
                // 打ち切りもこの「埋める長さ」で判断しないと、正解の経路を
                // 目標到達の直前で切ってしまう。
                var fillLength = added - kmerLength;
                if (fillLength > maxLength)
                {
                    continue;
                }

                if (fillLength >= minLength && kmer.AsSpan().SequenceEqual(targetKmer))
                {
                    found.Add(Reconstruct(nodes, currentIndex, fillLength));
                    if (found.Count > 1)
                    {
                        // 2本見つかった時点で一意には定まらない。
                        outcome = FillOutcome.Ambiguous;
                        return null;
                    }
                    continue;
                }

                if (nodes.Count > MaxStatesPerGap)
                {
                    outcome = FillOutcome.Ambiguous;
                    return null;
                }

                for (byte b = Consts.NucleotideID.A; b <= Consts.NucleotideID.T; b++)
                {
                    Array.Copy(kmer, 1, buffer, 0, kmerLength - 1);
                    buffer[kmerLength - 1] = b;
                    if (!index.Contains(buffer))
                    {
                        continue;
                    }
                    nodes.Add((currentIndex, b));
                    kmers.Add((byte[])buffer.Clone());
                    depths.Add(added + 1);
                    queue.Enqueue(nodes.Count - 1);
                }
            }

            if (found.Count == 1)
            {
                outcome = FillOutcome.Filled;
                return found[0];
            }

            outcome = found.Count > 1 ? FillOutcome.Ambiguous : FillOutcome.Unreachable;
            return null;
        }

        /// <summary>
        /// 親を辿って、継ぎ足した塩基列のうち先頭 fillLength 塩基を復元する。
        /// 末尾側(目標 k-mer と重なる分)は捨てる。
        /// </summary>
        private static string Reconstruct(List<(int Parent, byte Base)> nodes, int leafIndex, int fillLength)
        {
            List<byte> reversed = [];
            var index = leafIndex;
            while (index > 0)
            {
                reversed.Add(nodes[index].Base);
                index = nodes[index].Parent;
            }
            reversed.Reverse();
            return string.Concat(reversed.Take(fillLength).Select(Util.ByteToBaseString));
        }

        public static void Report(Stats stats)
        {
            if (stats.TotalGaps == 0)
            {
                Console.WriteLine("[Info] Gap filling: no gaps to fill.");
                return;
            }
            Console.WriteLine(
                $"[Info] Gap filling: {stats.FilledGaps}/{stats.TotalGaps} gap(s) closed with real sequence " +
                $"({stats.FilledBases:N0}bp of N replaced). " +
                $"{stats.AmbiguousGaps} left as N because more than one path fits, " +
                $"{stats.UnreachableGaps} because no path through the graph connects the two sides.");
        }
    }
}
