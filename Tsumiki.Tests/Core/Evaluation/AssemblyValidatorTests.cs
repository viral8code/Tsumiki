using Tsumiki.Common;
using Tsumiki.Core.Evaluation;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかを
    /// 確かめる自己検査の検証
    /// </summary>
    /// <remarks>
    /// リファレンス配列なしで「取りこぼし」と
    /// 「出しすぎ」を検出できることを固定する<br/>
    /// 「出しすぎ」の検出は特に重要で、総延長が実際のゲノムサイズより大きく
    /// なる原因はほぼこれ (実際、修正前は同じ配列を順鎖と逆鎖の両方で出力して
    /// いて総長がちょうど 2.009 倍に膨れていた)
    /// </remarks>
    public class AssemblyValidatorTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public AssemblyValidatorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_validator_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
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
        /// この検証で使う k 長
        /// </summary>
        private const int K = 21;

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_乱数配列(int p_長さ, int p_シード)
        {
            var l_rng = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_rng.Next(4)]));
        }

        /// <summary>
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_深さ">登録する深さ</param>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_索引(int p_深さ, params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            var l_index = new TrustedKmerIndex(this._tempDir);
            foreach (var l_seq in p_配列群)
            {
                var l_bytes = l_seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + K <= l_bytes.Length; i++)
                {
                    for (var l_rep = 0; l_rep < p_深さ; l_rep++)
                    {
                        l_index.V_登録(l_bytes.AsSpan(i, K), p_ワーカー番号: 0);
                    }
                }
            }
            _ = l_index.V_カットオフ(p_カットオフ: 2);
            return l_index;
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列群">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_Fasta(string p_ファイル名, params string[] p_配列群)
        {
            var l_path = Path.Combine(this._tempDir, p_ファイル名);
            using var l_writer = new FastaWriter(l_path);
            var l_id = 1;
            foreach (var l_seq in p_配列群)
            {
                l_writer.V_書き込み($"NODE{l_id++}", l_seq);
            }
            return l_path;
        }

        /// <summary>
        /// 与えた深さで k-mer 集合へ配列を登録する
        /// </summary>
        /// <param name="p_index">登録先の k-mer 集合</param>
        /// <param name="p_seq">登録する配列</param>
        /// <param name="p_深さ">登録する深さ</param>
        private static void V_登録_複数(TrustedKmerIndex p_index, string p_seq, int p_深さ)
        {
            var l_bytes = p_seq.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + K <= l_bytes.Length; i++)
            {
                for (var l_rep = 0; l_rep < p_深さ; l_rep++)
                {
                    p_index.V_登録(l_bytes.AsSpan(i, K), p_ワーカー番号: 0);
                }
            }
        }

        /// <summary>
        /// 入力を過不足なく再現したアセンブリは取りこぼしも出しすぎも報告しないことを確かめる
        /// </summary>
        [Fact]
        public void 入力を過不足なく再現したアセンブリは取りこぼしも出しすぎも報告しない()
        {
            var l_truth = V_乱数配列(600, p_シード: 101);
            using var l_index = this.V_構築_索引(p_深さ: 20, l_truth);

            var l_path = this.V_書き出し_Fasta("perfect.fasta", l_truth);
            var l_result = AssemblyValidator.Get_検査結果(l_path, l_index, K, p_単一コピー基準値: 20)!.Value;

            Assert.Equal(0, l_result.A_取りこぼし数);
            Assert.Equal(0, l_result.A_余分な延べ数);
        }

        /// <summary>
        /// 配列の後半を欠いたアセンブリは取りこぼした kmer を報告することを確かめる
        /// </summary>
        [Fact]
        public void 配列の後半を欠いたアセンブリは取りこぼしたkmerを報告する()
        {
            var l_truth = V_乱数配列(600, p_シード: 102);
            using var l_index = this.V_構築_索引(p_深さ: 20, l_truth);

            // 後半を落としたアセンブリ
            var l_path = this.V_書き出し_Fasta("truncated.fasta", l_truth[..300]);
            var l_result = AssemblyValidator.Get_検査結果(l_path, l_index, K, p_単一コピー基準値: 20)!.Value;

            Assert.True(l_result.A_取りこぼし数 > 0, "truncated assembly should report missing k-mers");
            // 600 bp の k-mer は 580 個、そのうち前半 300 bp に含まれるのは 280 個
            Assert.Equal(580 - 280, l_result.A_取りこぼし数);
            Assert.InRange(l_result.A_取りこぼし率, 45, 55);
        }

        /// <summary>
        /// 単一コピーの配列を 2 回出力してしまった場合、カバレッジは 1 コピー分しか
        /// 無いので「出しすぎ」として検出されなければならない
        /// </summary>
        /// <remarks>
        /// これが検出できないと、
        /// 総延長が水増しされていることに気付けない
        /// </remarks>
        [Fact]
        public void 単一コピーの配列を2回出力すると出しすぎとして検出される()
        {
            var l_truth = V_乱数配列(600, p_シード: 103);
            using var l_index = this.V_構築_索引(p_深さ: 20, l_truth);

            var l_path = this.V_書き出し_Fasta("duplicated.fasta", l_truth, l_truth);
            var l_result = AssemblyValidator.Get_検査結果(l_path, l_index, K, p_単一コピー基準値: 20)!.Value;

            Assert.Equal(0, l_result.A_取りこぼし数);
            // 各 k-mer が期待の 2 倍出ているので、延べ数の半分が余分
            Assert.Equal(580, l_result.A_出しすぎkmer種類数);
            Assert.Equal(580, l_result.A_余分な延べ数);
            Assert.InRange(l_result.A_出しすぎ率, 45, 55);
        }

        /// <summary>
        /// 逆相補で出力されていても同じ配列とみなされること (正規化の確認) を確かめる
        /// </summary>
        /// <remarks>
        /// これが効いていないと、逆鎖側の contig がすべて「取りこぼし」に見えてしまう
        /// </remarks>
        [Fact]
        public void 逆相補で出力しても同じ配列とみなされる()
        {
            var l_truth = V_乱数配列(600, p_シード: 104);
            using var l_index = this.V_構築_索引(p_深さ: 20, l_truth);

            var l_path = this.V_書き出し_Fasta("revcomp.fasta", Util.V_逆相補(l_truth));
            var l_result = AssemblyValidator.Get_検査結果(l_path, l_index, K, p_単一コピー基準値: 20)!.Value;

            Assert.Equal(0, l_result.A_取りこぼし数);
            Assert.Equal(0, l_result.A_余分な延べ数);
        }

        /// <summary>
        /// 2 コピー分のカバレッジがある反復配列を 2 回出力するのは正しいことを確かめる
        /// </summary>
        /// <remarks>
        /// これを「出しすぎ」と誤判定してはいけない
        /// </remarks>
        [Fact]
        public void 二コピー分のカバレッジがある反復配列を2回出力しても出しすぎにならない()
        {
            var l_single = V_乱数配列(600, p_シード: 105);
            var l_repeat = V_乱数配列(200, p_シード: 106);

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            using var l_index = new TrustedKmerIndex(this._tempDir);

            V_登録_複数(l_index, l_single, 20);
            V_登録_複数(l_index, l_repeat, 40); // 2コピー相当のカバレッジ
            _ = l_index.V_カットオフ(p_カットオフ: 2);

            var l_path = this.V_書き出し_Fasta("repeat_twice.fasta", l_single, l_repeat, l_repeat);
            var l_result = AssemblyValidator.Get_検査結果(l_path, l_index, K, p_単一コピー基準値: 20)!.Value;

            Assert.Equal(0, l_result.A_取りこぼし数);
            Assert.Equal(0, l_result.A_余分な延べ数);
        }
    }
}
