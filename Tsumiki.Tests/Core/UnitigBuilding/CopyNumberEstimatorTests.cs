using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// カバレッジから unitig のコピー数を推定する処理の検証
    /// </summary>
    /// <remarks>
    /// ゲノム中に 1 回しか現れない領域のカバレッジを基準値とすると、n 回現れる
    /// 反復配列にはリードが n 倍集まる<br/>
    /// したがってカバレッジ比を丸めれば
    /// コピー数になる<br/>
    /// これが分かると、反復配列かどうかをグラフの形ではなく
    /// 量的な根拠で判定でき、経路探索では「何回まで使ってよいか」の予算になる
    /// </remarks>
    public class CopyNumberEstimatorTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
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
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_length">作る長さ</param>
        /// <param name="p_seed">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_length, int p_seed)
        {
            var l_rng = new Random(p_seed);
            return string.Concat(Enumerable.Range(0, p_length).Select(_ => "ACGT"[l_rng.Next(4)]));
        }

        /// <summary>
        /// 配列の全 k-mer を、指定した深さだけ信頼できる kmer 索引へ登録する
        /// </summary>
        /// <param name="p_index">登録先の索引</param>
        /// <param name="p_seq">登録する配列</param>
        /// <param name="p_depth">各 k-mer を登録する回数</param>
        /// <param name="p_k">k-mer 長</param>
        private static void V_登録_全kmer(TrustedKmerIndex p_index, string p_seq, int p_depth, int p_k)
        {
            var l_bytes = p_seq.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + p_k <= l_bytes.Length; i++)
            {
                for (var rep = 0; rep < p_depth; rep++)
                {
                    p_index.V_登録(l_bytes.AsSpan(i, p_k), p_ワーカー番号: 0);
                }
            }
        }

        /// <summary>
        /// 単一コピーと 2 倍・4 倍コピーの配列をカバレッジ比から分離できる
        /// </summary>
        [Fact]
        public void 単一コピーと2倍_4倍コピーの配列を分離できる()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            // 単一コピー相当を 2 本 (長さで基準値を支配させる)、
            // 2 倍・4 倍のカバレッジで登録する配列を 1 本ずつ用意する
            var single1 = V_生成_乱数配列(400, p_seed: 1);
            var single2 = V_生成_乱数配列(400, p_seed: 2);
            var doubled = V_生成_乱数配列(120, p_seed: 3);
            var quadrupled = V_生成_乱数配列(120, p_seed: 4);

            using var index = new TrustedKmerIndex(this._tempDir);

            V_登録_全kmer(index, single1, 20, k);
            V_登録_全kmer(index, single2, 20, k);
            V_登録_全kmer(index, doubled, 40, k);
            V_登録_全kmer(index, quadrupled, 80, k);

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
        /// カバレッジがわずかに高いだけの配列を反復と誤判定してはいけない
        /// </summary>
        /// <remarks>
        /// 実データのカバレッジは領域ごとにかなりばらつくため、
        /// 1.5 倍未満は単一コピーとして扱う
        /// </remarks>
        [Fact]
        public void わずかに高いカバレッジの配列を単一コピーとして扱う()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var baselineSeq = V_生成_乱数配列(400, p_seed: 5);
            var slightlyHigher = V_生成_乱数配列(120, p_seed: 6);

            using var index = new TrustedKmerIndex(this._tempDir);

            V_登録_全kmer(index, baselineSeq, 20, k);
            V_登録_全kmer(index, slightlyHigher, 26, k); // 1.3倍

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
        public void kmer長より短いunitigはコピー数0でなく1になる()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var normal = V_生成_乱数配列(300, p_seed: 8);
            var tooShort = V_生成_乱数配列(10, p_seed: 9);

            using var index = new TrustedKmerIndex(this._tempDir);
            V_登録_全kmer(index, normal, 20, k);
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
        /// 大域基準値との比だけで見ると多コピーの反復に見える
        /// </summary>
        /// <remarks>
        /// しかし
        /// その単一コピー領域同士は分岐の無い (排他的な) 鎖で繋がっているため、
        /// 接続構造を使えば「大域とは水準が違うだけの単一コピー」だと分かる<br/>
        /// unicycler の copy depth propagation が解決する問題そのもの
        /// </remarks>
        [Fact]
        public void グラフを使うと高カバレッジのプラスミド骨格を単一コピーと認識できる()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            // 染色体相当 (長さで大域基準値=60 を支配する)
            var chromosome = V_生成_乱数配列(400, p_seed: 10);

            // プラスミド相当
            // 1 本の配列を k-1(=20) ずつ重ねて 3 本に切り出し、
            // 分岐の無い鎖 plasmid1 -> plasmid2 -> plasmid3 を作る
            var plasmidFull = V_生成_乱数配列(90, p_seed: 20);
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
        public void グラフを使っても孤立した高カバレッジunitigは補正されない()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            var chromosome = V_生成_乱数配列(400, p_seed: 11);
            var isolatedRepeat = V_生成_乱数配列(50, p_seed: 21); // 他のどれとも重ならない

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
        /// もっとも典型的なケース
        /// </summary>
        /// <remarks>
        /// 染色体側の成分とは一切繋がりが無い、
        /// 十分な長さを持つ「島」なので、大域基準値との比が高くても
        /// 単一コピーとみなしてよい (高コピープラスミド自身の水準で 1 コピー)
        /// </remarks>
        [Fact]
        public void グラフを使うと孤立した長いunitigを単独の単一コピーレプリコンと認識できる()
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            var chromosome = V_生成_乱数配列(400, p_seed: 12);
            var plasmid = V_生成_乱数配列(600, p_seed: 22); // 染色体とは無関係、500bp超

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
