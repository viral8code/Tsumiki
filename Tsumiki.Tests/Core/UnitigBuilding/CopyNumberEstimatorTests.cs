using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// カバレッジから unitig のコピー数を推定する処理の検証<br/>
    /// ゲノム中に 1 回しか現れない領域のカバレッジを基準値とすると、n 回現れる
    /// 反復配列にはリードが n 倍集まる<br/>
    /// したがってカバレッジ比を丸めれば
    /// コピー数になる<br/>
    /// これが分かると、反復配列かどうかをグラフの形ではなく
    /// 量的な根拠で判定でき、経路探索では「何回まで使ってよいか」の予算になる
    /// </summary>
    public class CopyNumberEstimatorTests : IDisposable
    {
        private readonly string _tempDir;

        public CopyNumberEstimatorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_copynumber_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);

            // これらのテストは長さ加重中央値のフォールバック経路を検証する
            // 他のテスト (KmerCutoffSelectorTests 等) が残した混合モデルの
            // 適合結果が ConfigurationManager 経由で漏れ込まないようにする
            ConfigurationManager.A_スペクトルモデル = null;
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        [Fact]
        public void Estimate_SeparatesSingleCopyFromTwoCopyAndFourCopySequences()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            // 単一コピー相当を 2 本 (長さで基準値を支配させる)、
            // 2 倍・4 倍のカバレッジで登録する配列を 1 本ずつ用意する
            var single1 = RandomSequence(400, seed: 1);
            var single2 = RandomSequence(400, seed: 2);
            var doubled = RandomSequence(120, seed: 3);
            var quadrupled = RandomSequence(120, seed: 4);

            using var index = new TrustedKmerIndex(this._tempDir);

            void Add(string seq, int depth)
            {
                var bytes = seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, k), p_ワーカー番号: 0);
                    }
                }
            }

            Add(single1, 20);
            Add(single2, 20);
            Add(doubled, 40);
            Add(quadrupled, 80);

            _ = index.V_カットオフ(p_カットオフ: 2);

            Dictionary<int, string> unitigs = new()
            {
                [1] = single1,
                [2] = single2,
                [3] = doubled,
                [4] = quadrupled,
            };
            var lengths = unitigs.ToDictionary(kv => kv.Key, kv => kv.Value.Length);

            var coverage = CopyNumberEstimator.Get_カバレッジ(index, unitigs, k);
            var result = CopyNumberEstimator.Get_推定結果(coverage, lengths);

            // 基準値は長さ加重中央値なので、長い単一コピー配列の水準になるはず
            Assert.InRange(result.A_単一コピー基準値, 15, 25);

            Assert.Equal(1, result.A_コピー数[1]);
            Assert.Equal(1, result.A_コピー数[2]);
            Assert.Equal(2, result.A_コピー数[3]);
            Assert.Equal(4, result.A_コピー数[4]);
        }

        /// <summary>
        /// カバレッジがわずかに高いだけの配列を反復と誤判定してはいけない<br/>
        /// 実データのカバレッジは領域ごとにかなりばらつくため、
        /// 1.5 倍未満は単一コピーとして扱う
        /// </summary>
        [Fact]
        public void Estimate_TreatsMildlyElevatedCoverageAsSingleCopy()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var baselineSeq = RandomSequence(400, seed: 5);
            var slightlyHigher = RandomSequence(120, seed: 6);

            using var index = new TrustedKmerIndex(this._tempDir);

            void Add(string seq, int depth)
            {
                var bytes = seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, k), p_ワーカー番号: 0);
                    }
                }
            }

            Add(baselineSeq, 20);
            Add(slightlyHigher, 26); // 1.3倍

            _ = index.V_カットオフ(p_カットオフ: 2);

            Dictionary<int, string> unitigs = new() { [1] = baselineSeq, [2] = slightlyHigher };
            var lengths = unitigs.ToDictionary(kv => kv.Key, kv => kv.Value.Length);

            var coverage = CopyNumberEstimator.Get_カバレッジ(index, unitigs, k);
            var result = CopyNumberEstimator.Get_推定結果(coverage, lengths);

            Assert.Equal(1, result.A_コピー数[2]);
        }

        /// <summary>
        /// k-mer 長より短い unitig はカバレッジを測れないが、
        /// コピー数 0 にして経路から締め出してはいけない (配列自体は存在する)
        /// </summary>
        [Fact]
        public void Estimate_UnitigShorterThanKmer_GetsCopyNumberOneRatherThanZero()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var normal = RandomSequence(300, seed: 8);
            var tooShort = RandomSequence(10, seed: 9);

            using var index = new TrustedKmerIndex(this._tempDir);
            var bytes = normal.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + k <= bytes.Length; i++)
            {
                for (var rep = 0; rep < 20; rep++)
                {
                    index.V_登録(bytes.AsSpan(i, k), p_ワーカー番号: 0);
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);

            Dictionary<int, string> unitigs = new() { [1] = normal, [2] = tooShort };
            var lengths = unitigs.ToDictionary(kv => kv.Key, kv => kv.Value.Length);

            var coverage = CopyNumberEstimator.Get_カバレッジ(index, unitigs, k);
            var result = CopyNumberEstimator.Get_推定結果(coverage, lengths);

            Assert.Equal(0, coverage[2]);
            Assert.Equal(1, result.A_コピー数[2]);
        }

        /// <summary>
        /// プラスミドのように染色体とは異なるカバレッジ水準を持つ領域は、
        /// 大域基準値との比だけで見ると多コピーの反復に見える<br/>
        /// しかし
        /// その単一コピー領域同士は分岐の無い (排他的な) 鎖で繋がっているため、
        /// 接続構造を使えば「大域とは水準が違うだけの単一コピー」だと分かる<br/>
        /// unicycler の copy depth propagation が解決する問題そのもの
        /// </summary>
        [Fact]
        public void Estimate_WithGraph_RecognisesAHighCoveragePlasmidBackboneAsSingleCopy()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            // 染色体相当 (長さで大域基準値=60 を支配する)
            var chromosome = RandomSequence(400, seed: 10);

            // プラスミド相当
            // 1 本の配列を k-1(=20) ずつ重ねて 3 本に切り出し、
            // 分岐の無い鎖 plasmid1 -> plasmid2 -> plasmid3 を作る
            var plasmidFull = RandomSequence(90, seed: 20);
            var plasmid1 = plasmidFull[..40];
            var plasmid2 = plasmidFull[20..60];
            var plasmid3 = plasmidFull[40..90];

            var fastaPath = Path.Combine(this._tempDir, "unitigs.fasta");
            using (var writer = new FastaWriter(fastaPath))
            {
                writer.V_書き込み(1, chromosome);
                writer.V_書き込み(2, plasmid1);
                writer.V_書き込み(3, plasmid2);
                writer.V_書き込み(4, plasmid3);
            }

            var contigMaker = new ContigMaker(fastaPath);
            var graph = contigMaker.Get_グラフ();

            Dictionary<int, double> coverage = new()
            {
                [1] = 60.0,
                [2] = 300.0, // 大域基準値(60)との比は5倍 -> 単独では多コピー判定
                [3] = 300.0,
                [4] = 300.0,
            };
            Dictionary<int, int> lengths = new()
            {
                [1] = chromosome.Length,
                [2] = plasmid1.Length,
                [3] = plasmid2.Length,
                [4] = plasmid3.Length,
            };

            var withoutGraph = CopyNumberEstimator.Get_推定結果(coverage, lengths);
            Assert.Equal(5, withoutGraph.A_コピー数[2]);
            Assert.Equal(5, withoutGraph.A_コピー数[3]);
            Assert.Equal(5, withoutGraph.A_コピー数[4]);

            var withGraph = CopyNumberEstimator.Get_推定結果(coverage, lengths, graph);
            Assert.Equal(1, withGraph.A_コピー数[2]);
            Assert.Equal(1, withGraph.A_コピー数[3]);
            Assert.Equal(1, withGraph.A_コピー数[4]);
            // 染色体側は元々単一コピー判定であり、接続補正の対象にもならない
            Assert.Equal(1, withGraph.A_コピー数[1]);
        }

        /// <summary>
        /// 排他的に繋がる相手がいない (孤立した) 高カバレッジ unitig は、
        /// 比較材料が無いため接続補正の対象にせず、大域基準値との比のまま残す
        /// </summary>
        [Fact]
        public void Estimate_WithGraph_LeavesAnIsolatedHighCoverageUnitigUnchanged()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            var chromosome = RandomSequence(400, seed: 11);
            var isolatedRepeat = RandomSequence(50, seed: 21); // 他のどれとも重ならない

            var fastaPath = Path.Combine(this._tempDir, "unitigs.fasta");
            using (var writer = new FastaWriter(fastaPath))
            {
                writer.V_書き込み(1, chromosome);
                writer.V_書き込み(2, isolatedRepeat);
            }

            var contigMaker = new ContigMaker(fastaPath);
            var graph = contigMaker.Get_グラフ();

            Dictionary<int, double> coverage = new() { [1] = 60.0, [2] = 300.0 };
            Dictionary<int, int> lengths = new() { [1] = chromosome.Length, [2] = isolatedRepeat.Length };

            var result = CopyNumberEstimator.Get_推定結果(coverage, lengths, graph);

            Assert.Equal(5, result.A_コピー数[2]);
        }

        /// <summary>
        /// 小さなプラスミドが分岐無しの 1 本の unitig にきれいに閉じた、
        /// もっとも典型的なケース<br/>
        /// 染色体側の成分とは一切繋がりが無い、
        /// 十分な長さを持つ「島」なので、大域基準値との比が高くても
        /// 単一コピーとみなしてよい (高コピープラスミド自身の水準で 1 コピー)
        /// </summary>
        [Fact]
        public void Estimate_WithGraph_RecognisesAnIsolatedLongUnitigAsItsOwnSingleCopyReplicon()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            var chromosome = RandomSequence(400, seed: 12);
            var plasmid = RandomSequence(600, seed: 22); // 染色体とは無関係、500bp超

            var fastaPath = Path.Combine(this._tempDir, "unitigs.fasta");
            using (var writer = new FastaWriter(fastaPath))
            {
                writer.V_書き込み(1, chromosome);
                writer.V_書き込み(2, plasmid);
            }

            var contigMaker = new ContigMaker(fastaPath);
            var graph = contigMaker.Get_グラフ();

            Dictionary<int, double> coverage = new() { [1] = 60.0, [2] = 300.0 };
            Dictionary<int, int> lengths = new() { [1] = chromosome.Length, [2] = plasmid.Length };

            var result = CopyNumberEstimator.Get_推定結果(coverage, lengths, graph);

            Assert.Equal(1, result.A_コピー数[2]);
        }
    }
}
