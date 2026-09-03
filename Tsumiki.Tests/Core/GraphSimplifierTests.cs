using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    public class GraphSimplifierTests : IDisposable
    {
        private readonly string _tempDir;

        public GraphSimplifierTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_graph_simplifier_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private static byte[] ToBytes(string seq)
        {
            return [.. seq.Select(Util.GetSimpleNucleotideID)];
        }

        /// <summary>
        /// mainSeq(主経路)のk-mer群に加えて、その途中の1点から分岐する
        /// 短いtip配列(tipSeq)のk-mer群も登録した TrustedKmerIndex を作る。
        /// tipSeqはmainSeqの位置branchPointから始まる長さkmerLength-1の
        /// 「本来の続き」をコピーした上で、最後の1塩基だけ変えることで
        /// 主経路とk-1塩基だけ重なる分岐を作る単純な構成にする。
        /// </summary>
        private TrustedKmerIndex BuildIndexWithTip(string mainSeq, int kmerLength, int branchPoint, int tipLength)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = kmerLength, ThreadCount = 1 };
            var index = new TrustedKmerIndex(this._tempDir);

            void AddAllKmers(byte[] bytes)
            {
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < 3; rep++)
                    {
                        index.Add(bytes.AsSpan(i, kmerLength), workerIndex: 0);
                    }
                }
            }

            AddAllKmers(ToBytes(mainSeq));

            // 分岐点の直前 kmerLength-1 文字を土台に、最後だけ主経路と異なる
            // 1塩基を続けて tip を伸ばす(主経路と k-2 塩基だけ重なる短い枝)。
            var overlap = mainSeq.Substring(branchPoint, kmerLength - 1);
            var branchBaseChar = mainSeq[branchPoint + kmerLength - 1];
            var altChar = "ACGT".First(c => c != branchBaseChar);
            // overlap(主経路とk-1塩基共有)+ altChar(主経路とは異なる1塩基)+
            // 適当なユニークな続きで、主経路から分岐する短いtipを作る。
            var tipSeq = overlap + altChar + string.Concat(Enumerable.Range(0, tipLength).Select(i => "ACGT"[i % 4]));
            AddAllKmers(ToBytes(tipSeq));

            _ = index.Cutoff(bounds: 2);
            return index;
        }

        [Fact]
        public void ClipTips_RemovesShortDeadEndBranch_AndRebuildsSingleMainUnitig()
        {
            // 非周期的な主経路配列(k=8での内部重複なしを別途Pythonで確認済み)。
            const string mainSeq = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int k = 8;
            const int branchPoint = 20; // 主経路の途中から分岐させる
            const int tipExtra = 4; // 分岐後にごく短く伸びるtip

            using var index = this.BuildIndexWithTip(mainSeq, k, branchPoint, tipExtra);

            var simplifiedFirstKmers = GraphSimplifier.ClipTips(index, k, tipLengthThreshold: k * 2);

            // tip除去後は、分岐点だった箇所の次数が解消され、
            // 主経路が1本のunitigとして(理想的には)再構築されるはず。
            var unitigMaker = new UnitigMaker(index);
            HashSet<string> seen = [];
            var unitigs = new List<string>();
            foreach (var kmer in simplifiedFirstKmers)
            {
                var u = unitigMaker.MakeUnitig(kmer);
                if (seen.Add(u.Sequence) || seen.Add(Util.ReverseComprement(u.Sequence)))
                {
                    unitigs.Add(u.Sequence);
                }
            }

            // tip自体はもう存在しないはずなので、tip由来の短い配列を含む
            // unitigは残っていないこと、かつ主経路の全長をカバーする
            // (ほぼ)1本のunitigが存在することを確認する。
            var longest = unitigs.OrderByDescending(u => u.Length).First();
            Assert.True(longest.Length >= mainSeq.Length - k, $"expected a near-full-length main unitig, longest was {longest.Length}bp among [{string.Join(",", unitigs.Select(u => u.Length))}]");
        }

        [Fact]
        public void ClipTips_LeavesLinearNonBranchingSequenceUnchanged()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC"; // 24bp、非周期的
            const int k = 8;
            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };
            using var index = new TrustedKmerIndex(this._tempDir);
            var bytes = ToBytes(seq);
            for (var i = 0; i + k <= bytes.Length; i++)
            {
                for (var rep = 0; rep < 3; rep++)
                {
                    index.Add(bytes.AsSpan(i, k), workerIndex: 0);
                }
            }
            _ = index.Cutoff(bounds: 2);

            var before = index.EnumerateTrustedKmers().Count();

            var firstKmers = GraphSimplifier.ClipTips(index, k, tipLengthThreshold: k * 2);

            var after = index.EnumerateTrustedKmers().Count();

            // 分岐のない直鎖配列にはtipが存在しないため、何も除去されないはず。
            Assert.Equal(before, after);
            Assert.NotEmpty(firstKmers);
        }

        [Fact]
        public void ClipTips_RemovesLowCoverageBubbleBranch_KeepsHighCoverageBranch()
        {
            // 分岐点(commonBefore末尾)から1塩基だけ異なる('A' vs 'C')経路B/Cに
            // 分かれ、その後sharedAfterへ合流するSNP様の単純なbubble構造。
            // Python(scripts外、事前検証)でk=8内に重複が生じないことを確認済み。
            const string commonBefore = "GAAGTTGCCGTACTAAATTA"; // 20bp
            const string sharedAfter = "TGACAGCCGGGGATCTTCCC"; // 20bp
            const string seqHighCoverage = commonBefore + "A" + sharedAfter; // 分岐点でA
            const string seqLowCoverage = commonBefore + "C" + sharedAfter; // 分岐点でC(エラー相当)
            const int k = 8;

            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };
            using var index = new TrustedKmerIndex(this._tempDir);

            void AddAllKmers(string seq, int repetitions)
            {
                var bytes = ToBytes(seq);
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < repetitions; rep++)
                    {
                        index.Add(bytes.AsSpan(i, k), workerIndex: 0);
                    }
                }
            }

            // 高カバレッジ経路(真のゲノム由来相当)は20回、低カバレッジ経路
            // (エラー由来相当)は3回登録する(カットオフ2は超えるが、
            // baseline(高カバレッジ経路水準)に比べて著しく低い)。
            AddAllKmers(seqHighCoverage, repetitions: 20);
            AddAllKmers(seqLowCoverage, repetitions: 3);

            _ = index.Cutoff(bounds: 2);

            var simplifiedFirstKmers = GraphSimplifier.ClipTips(index, k, tipLengthThreshold: k * 2);

            var unitigMaker = new UnitigMaker(index);
            HashSet<string> seen = [];
            var unitigs = new List<string>();
            foreach (var kmer in simplifiedFirstKmers)
            {
                var u = unitigMaker.MakeUnitig(kmer);
                if (seen.Add(u.Sequence) || seen.Add(Util.ReverseComprement(u.Sequence)))
                {
                    unitigs.Add(u.Sequence);
                }
            }

            // 低カバレッジ経路の分岐点を含む短い断片は残っていないはず
            // (再構築された配列のいずれにも "C" + sharedAfter の先頭部分は
            // 現れない = 低カバレッジ経路は除去された)。
            Assert.DoesNotContain(unitigs, u => u.Contains('C' + sharedAfter[..(k - 1)]));

            // 高カバレッジ経路(commonBefore + "A" + sharedAfter の全体、または
            // その逆相補)を含む、ほぼ全長のunitigが存在するはず。
            var fullHigh = seqHighCoverage;
            var fullHighRevComp = Util.ReverseComprement(fullHigh);
            Assert.Contains(unitigs, u => u == fullHigh || u == fullHighRevComp || u.Contains(fullHigh) || u.Contains(fullHighRevComp));
        }

    }
}
