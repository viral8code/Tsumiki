using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    public class GraphSimplifierTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public GraphSimplifierTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_graph_simplifier_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);

            // これらのテストは比率ベースの tip 判定を検証する
            // 他のテストが残した
            // 混合モデルの適合結果が「無条件に信頼する下限」として漏れ込み、
            // 判定を横取りしないようにする
            ConfigurationManager.A_スペクトルモデル = null;
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="seq">元の配列</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] ToBytes(string seq)
        {
            return [.. seq.Select(Util.Get_塩基ID)];
        }

        /// <summary>
        /// mainSeq(主経路) の k-mer 群に加えて、その途中の 1 点から分岐する
        /// 短いtip配列 (tipSeq) の k-mer 群も登録した TrustedKmerIndex を作る
        /// </summary>
        /// <remarks>
        /// tipSeq は mainSeq の位置 branchPoint から始まる長さ kmerLength-1 の
        /// 「本来の続き」をコピーした上で、最後の 1 塩基だけ変えることで
        /// 主経路と k-1 塩基だけ重なる分岐を作る単純な構成にする
        /// </remarks>
        private TrustedKmerIndex BuildIndexWithTip(
            string mainSeq, int kmerLength, int branchPoint, int tipLength,
            int mainRepetitions = 10, int tipRepetitions = 2)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = kmerLength, A_スレッド数 = 1 };
            var index = new TrustedKmerIndex(this._tempDir);

            void AddAllKmers(byte[] bytes, int repetitions)
            {
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < repetitions; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, kmerLength), p_ワーカー番号: 0);
                    }
                }
            }

            AddAllKmers(ToBytes(mainSeq), mainRepetitions);

            // 分岐点の直前 kmerLength-1 文字を土台に、最後だけ主経路と異なる
            // 1 塩基を続けて tip を伸ばす (主経路と k-2 塩基だけ重なる短い枝)
            var overlap = mainSeq.Substring(branchPoint, kmerLength - 1);
            var branchBaseChar = mainSeq[branchPoint + kmerLength - 1];
            var altChar = "ACGT".First(c => c != branchBaseChar);
            // overlap(主経路と k-1 塩基共有)+ altChar(主経路とは異なる 1 塩基)+
            // 適当なユニークな続きで、主経路から分岐する短い tip を作る
            var tipSeq = overlap + altChar + string.Concat(Enumerable.Range(0, tipLength).Select(i => "ACGT"[i % 4]));
            AddAllKmers(ToBytes(tipSeq), tipRepetitions);

            _ = index.V_カットオフ(p_カットオフ: 2);
            return index;
        }

        [Fact]
        public void ClipTips_RemovesShortDeadEndBranch_AndRebuildsSingleMainUnitig()
        {
            // 非周期的な主経路配列 (k=8 での内部重複なしを別途 Python で確認済み)
            const string mainSeq = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int k = 8;
            const int branchPoint = 20; // 主経路の途中から分岐させる
            const int tipExtra = 4; // 分岐後にごく短く伸びるtip

            using var index = this.BuildIndexWithTip(mainSeq, k, branchPoint, tipExtra);

            var simplifiedFirstKmers = GraphSimplifier.V_除去_tip(index, k, p_tip長閾値: k * 2);

            // tip 除去後は、分岐点だった箇所の次数が解消され、
            // 主経路が 1 本のunitigとして (理想的には) 再構築されるはず
            var unitigMaker = new UnitigMaker(index);
            HashSet<string> seen = [];
            var unitigs = new List<string>();
            foreach (var kmer in simplifiedFirstKmers)
            {
                var u = unitigMaker.Get_ユニティグ(kmer);
                if (seen.Contains(u.A_配列) || seen.Contains(Util.V_逆相補(u.A_配列)))
                {
                    continue;
                }
                _ = seen.Add(u.A_配列);
                _ = seen.Add(Util.V_逆相補(u.A_配列));
                unitigs.Add(u.A_配列);
            }

            // tip 自体はもう存在しないはずなので、tip 由来の短い配列を含む
            // unitig は残っていないこと、かつ主経路の全長をカバーする
            // (ほぼ) 1 本の unitig が存在することを確認する
            var longest = unitigs.OrderByDescending(u => u.Length).First();
            Assert.True(longest.Length >= mainSeq.Length - k, $"expected a near-full-length main unitig, longest was {longest.Length}bp among [{string.Join(",", unitigs.Select(u => u.Length))}]");
        }

        /// <summary>
        /// 行き止まりでも、カバレッジが主経路並みなら除去しないこと
        /// </summary>
        /// <remarks>
        /// カバレッジの切れ目で孤立した実配列がこの形になるため、
        /// 行き止まりというだけで消すとゲノム被覆率を落とす
        /// </remarks>
        [Fact]
        public void ClipTips_KeepsAShortDeadEndBranchWithMainPathCoverage()
        {
            const string mainSeq = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int k = 8;
            const int branchPoint = 20;
            const int tipExtra = 4;

            using var index = this.BuildIndexWithTip(
                mainSeq, k, branchPoint, tipExtra, mainRepetitions: 10, tipRepetitions: 10);

            var before = index.Get_信頼kmer一覧().Count();
            _ = GraphSimplifier.V_除去_tip(index, k, p_tip長閾値: k * 2);

            Assert.Equal(before, index.Get_信頼kmer一覧().Count());
        }

        [Fact]
        public void ClipTips_LeavesLinearNonBranchingSequenceUnchanged()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC"; // 24bp、非周期的
            const int k = 8;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            using var index = new TrustedKmerIndex(this._tempDir);
            var bytes = ToBytes(seq);
            for (var i = 0; i + k <= bytes.Length; i++)
            {
                for (var rep = 0; rep < 3; rep++)
                {
                    index.V_登録(bytes.AsSpan(i, k), p_ワーカー番号: 0);
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);

            var before = index.Get_信頼kmer一覧().Count();

            var firstKmers = GraphSimplifier.V_除去_tip(index, k, p_tip長閾値: k * 2);

            var after = index.Get_信頼kmer一覧().Count();

            // 分岐のない直鎖配列には tip が存在しないため、何も除去されないはず
            Assert.Equal(before, after);
            Assert.NotEmpty(firstKmers);
        }

        [Fact]
        public void ClipTips_RemovesLowCoverageBubbleBranch_KeepsHighCoverageBranch()
        {
            // 分岐点 (commonBefore 末尾) から 1 塩基だけ異なる ('A' vs 'C') 経路 B/C に
            // 分かれ、その後 sharedAfter へ合流する SNP 様の単純な bubble 構造
            // Python(scripts 外、事前検証) で k=8 内に重複が生じないことを確認済み
            const string commonBefore = "GAAGTTGCCGTACTAAATTA"; // 20bp
            const string sharedAfter = "TGACAGCCGGGGATCTTCCC"; // 20bp
            const string seqHighCoverage = commonBefore + "A" + sharedAfter; // 分岐点でA
            const string seqLowCoverage = commonBefore + "C" + sharedAfter; // 分岐点でC(エラー相当)
            const int k = 8;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            using var index = new TrustedKmerIndex(this._tempDir);

            void AddAllKmers(string seq, int repetitions)
            {
                var bytes = ToBytes(seq);
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < repetitions; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, k), p_ワーカー番号: 0);
                    }
                }
            }

            // 高カバレッジ経路 (真のゲノム由来相当) は 20 回、低カバレッジ経路
            // (エラー由来相当) は 3 回登録する (カットオフ 2 は超えるが、
            // baseline(高カバレッジ経路水準) に比べて著しく低い)
            AddAllKmers(seqHighCoverage, repetitions: 20);
            AddAllKmers(seqLowCoverage, repetitions: 3);

            _ = index.V_カットオフ(p_カットオフ: 2);

            var simplifiedFirstKmers = GraphSimplifier.V_除去_tip(index, k, p_tip長閾値: k * 2);

            var unitigMaker = new UnitigMaker(index);
            HashSet<string> seen = [];
            var unitigs = new List<string>();
            foreach (var kmer in simplifiedFirstKmers)
            {
                var u = unitigMaker.Get_ユニティグ(kmer);
                if (seen.Contains(u.A_配列) || seen.Contains(Util.V_逆相補(u.A_配列)))
                {
                    continue;
                }
                _ = seen.Add(u.A_配列);
                _ = seen.Add(Util.V_逆相補(u.A_配列));
                unitigs.Add(u.A_配列);
            }

            // 低カバレッジ経路の分岐点を含む短い断片は残っていないはず
            // (再構築された配列のいずれにも "C" + sharedAfter の先頭部分は
            // 現れない = 低カバレッジ経路は除去された)
            Assert.DoesNotContain(unitigs, u => u.Contains('C' + sharedAfter[..(k - 1)]));

            // 高カバレッジ経路 (commonBefore + "A" + sharedAfter の全体、または
            // その逆相補) を含む、ほぼ全長の unitig が存在するはず
            var fullHigh = seqHighCoverage;
            var fullHighRevComp = Util.V_逆相補(fullHigh);
            Assert.Contains(unitigs, u => u == fullHigh || u == fullHighRevComp || u.Contains(fullHigh) || u.Contains(fullHighRevComp));
        }

        /// <summary>
        /// 2 つの異なる経路が同じ配列へ合流する構造 (reverse bubble) で、
        /// 合流後の共有配列が複数の unitig に重複して現れないことを確認する
        /// </summary>
        /// <remarks>
        /// unitig の定義は「内部の全節点が入次数 1 かつ出次数 1 の極大パス」であり、
        /// 合流点 (入次数 2) からは別の unitig が始まらなければならない<br/>
        /// この規則が無いと、両方の経路の walk が共有配列を走り抜けてしまい、
        /// 同じ配列を 2 度出力する (実データで k-mer 延べ数が実内容の 1.43 倍に
        /// 膨らんでいた原因)
        /// </remarks>
        [Fact]
        public void MakeUnitig_StopsAtMergePoint_SoSharedSuffixIsNotDuplicated()
        {
            const string prefixA = "ATATCACACCCAACCTTCAA";
            const string prefixB = "ATGCCGTGCCCTAACGCCCT";
            const string shared = "AATCCTGCGCTAGGGGTTGCAGCGACCAGA";
            const int k = 8;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            using var index = new TrustedKmerIndex(this._tempDir);

            void AddAllKmers(string seq)
            {
                var bytes = ToBytes(seq);
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < 10; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, k), p_ワーカー番号: 0);
                    }
                }
            }

            AddAllKmers(prefixA + shared);
            AddAllKmers(prefixB + shared);

            var firstKmers = index.V_カットオフ(p_カットオフ: 2);

            var unitigMaker = new UnitigMaker(index);
            HashSet<string> seen = [];
            var unitigs = new List<string>();
            foreach (var kmer in firstKmers)
            {
                var u = unitigMaker.Get_ユニティグ(kmer);
                if (seen.Contains(u.A_配列) || seen.Contains(Util.V_逆相補(u.A_配列)))
                {
                    continue;
                }
                _ = seen.Add(u.A_配列);
                _ = seen.Add(Util.V_逆相補(u.A_配列));
                unitigs.Add(u.A_配列);
            }

            // 共有配列の先頭 k-mer(またはその逆相補) が、全 unitig を通じて
            // 延べ 1 回しか現れないこと (=重複出力されていないこと) を確認する
            var sharedKmer = shared[..k];
            var sharedKmerRc = Util.V_逆相補(sharedKmer);
            var occurrences = unitigs.Sum(u => CountOccurrences(u, sharedKmer) + CountOccurrences(u, sharedKmerRc));
            Assert.Equal(1, occurrences);
        }

        /// <summary>
        /// 部分文字列が現れる回数を数える
        /// </summary>
        /// <param name="haystack">探される側</param>
        /// <param name="needle">探す文字列</param>
        /// <returns>現れた回数</returns>
        private static int CountOccurrences(string haystack, string needle)
        {
            var count = 0;
            for (var i = 0; i + needle.Length <= haystack.Length; i++)
            {
                if (string.CompareOrdinal(haystack, i, needle, 0, needle.Length) == 0)
                {
                    count++;
                }
            }
            return count;
        }
    }
}
