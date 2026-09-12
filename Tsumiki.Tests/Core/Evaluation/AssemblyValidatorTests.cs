using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかを確かめる自己検査の検証
    /// </summary>
    public class AssemblyValidatorTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 21;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public AssemblyValidatorTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_validator_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 入力を過不足なく再現したアセンブリは取りこぼしも出しすぎも報告しないことを確かめる
        /// </summary>
        [Fact]
        public void V_入力を過不足なく再現したアセンブリは取りこぼしも出しすぎも報告しない()
        {
            var l_正解 = V_乱数配列(600, p_シード: 101);
            using var l_インデックス = this.V_構築_索引(p_深さ: 20, l_正解);

            var l_パス = this.V_書き出し_Fasta("perfect.fasta", l_正解);
            var l_結果 = AssemblyValidator.Get_検査結果(l_パス, l_インデックス, k長, p_単一コピー基準値: 20D)!.Value;

            Assert.Equal(0L, l_結果.A_取りこぼし数);
            Assert.Equal(0L, l_結果.A_余分な延べ数);
        }

        /// <summary>
        /// 配列の後半を欠いたアセンブリは取りこぼした kmer を報告することを確かめる
        /// </summary>
        [Fact]
        public void V_配列の後半を欠いたアセンブリは取りこぼしたkmerを報告する()
        {
            var l_正解 = V_乱数配列(600, p_シード: 102);
            using var l_インデックス = this.V_構築_索引(p_深さ: 20, l_正解);

            // 後半を落としたアセンブリ
            var l_パス = this.V_書き出し_Fasta("truncated.fasta", l_正解[..300]);
            var l_結果 = AssemblyValidator.Get_検査結果(l_パス, l_インデックス, k長, p_単一コピー基準値: 20D)!.Value;

            Assert.True(l_結果.A_取りこぼし数 > 0L, "truncated assembly should report missing k-mers");
            // 600 bp の k-mer は 580 個、そのうち前半 300 bp に含まれるのは 280 個
            Assert.Equal(580 - 280, l_結果.A_取りこぼし数);
            Assert.InRange(l_結果.A_取りこぼし率, 45D, 55D);
        }

        /// <summary>
        /// 単一コピーの配列を 2 回出力してしまった場合、カバレッジは 1 コピー分しか無いので「出しすぎ」として検出されなければならない
        /// </summary>
        /// <remarks>
        /// これが検出できないと、総延長が水増しされていることに気付けない
        /// </remarks>
        [Fact]
        public void V_単一コピーの配列を2回出力すると出しすぎとして検出される()
        {
            var l_正解 = V_乱数配列(600, p_シード: 103);
            using var l_インデックス = this.V_構築_索引(p_深さ: 20, l_正解);

            var l_パス = this.V_書き出し_Fasta("duplicated.fasta", l_正解, l_正解);
            var l_結果 = AssemblyValidator.Get_検査結果(l_パス, l_インデックス, k長, p_単一コピー基準値: 20D)!.Value;

            Assert.Equal(0L, l_結果.A_取りこぼし数);
            // 各 k-mer が期待の 2 倍出ているので、延べ数の半分が余分
            Assert.Equal(580L, l_結果.A_出しすぎkmer種類数);
            Assert.Equal(580L, l_結果.A_余分な延べ数);
            Assert.InRange(l_結果.A_出しすぎ率, 45D, 55D);
        }

        /// <summary>
        /// 逆相補で出力されていても同じ配列とみなされること (正規化の確認) を確かめる
        /// </summary>
        /// <remarks>
        /// これが効いていないと、逆鎖側の contig がすべて「取りこぼし」に見えてしまう
        /// </remarks>
        [Fact]
        public void V_逆相補で出力しても同じ配列とみなされる()
        {
            var l_正解 = V_乱数配列(600, p_シード: 104);
            using var l_インデックス = this.V_構築_索引(p_深さ: 20, l_正解);

            var l_パス = this.V_書き出し_Fasta("revcomp.fasta", Util.V_逆相補(l_正解));
            var l_結果 = AssemblyValidator.Get_検査結果(l_パス, l_インデックス, k長, p_単一コピー基準値: 20D)!.Value;

            Assert.Equal(0L, l_結果.A_取りこぼし数);
            Assert.Equal(0L, l_結果.A_余分な延べ数);
        }

        /// <summary>
        /// 2 コピー分のカバレッジがある反復配列を 2 回出力するのは正しいことを確かめる
        /// </summary>
        /// <remarks>
        /// これを「出しすぎ」と誤判定してはいけない
        /// </remarks>
        [Fact]
        public void V_二コピー分のカバレッジがある反復配列を2回出力しても出しすぎにならない()
        {
            var l_single = V_乱数配列(600, p_シード: 105);
            var l_反復配列 = V_乱数配列(200, p_シード: 106);

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            V_登録_複数(l_インデックス, l_single, 20);
            V_登録_複数(l_インデックス, l_反復配列, 40); // 2 コピー相当のカバレッジ
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            var l_パス = this.V_書き出し_Fasta("repeat_twice.fasta", l_single, l_反復配列, l_反復配列);
            var l_結果 = AssemblyValidator.Get_検査結果(l_パス, l_インデックス, k長, p_単一コピー基準値: 20D)!.Value;

            Assert.Equal(0L, l_結果.A_取りこぼし数);
            Assert.Equal(0L, l_結果.A_余分な延べ数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_乱数配列(int p_長さ, int p_シード)
        {
            var l_乱数生成器 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数生成器.Next(4)]));
        }

        /// <summary>
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_深さ">登録する深さ</param>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_索引(int p_深さ, params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            foreach (var l_配列 in p_配列群)
            {
                var l_バイト列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + k長 <= l_バイト列.Length; i++)
                {
                    for (var l_反復回数 = 0; l_反復回数 < p_深さ; l_反復回数++)
                    {
                        l_インデックス.V_登録(l_バイト列.AsSpan(i, k長));
                    }
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列群">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_Fasta(string p_ファイル名, params string[] p_配列群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new FastaWriter(l_パス);
            var l_ID = 1;
            foreach (var l_配列 in p_配列群)
            {
                l_書き込み.V_書き込み($"NODE{l_ID++}", l_配列);
            }
            return l_パス;
        }

        /// <summary>
        /// 与えた深さで k-mer 集合へ配列を登録する
        /// </summary>
        /// <param name="p_インデックス">登録先の k-mer 集合</param>
        /// <param name="p_配列">登録する配列</param>
        /// <param name="p_深さ">登録する深さ</param>
        private static void V_登録_複数(TrustedKmerIndex p_インデックス, string p_配列, int p_深さ)
        {
            var l_バイト列 = p_配列.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + k長 <= l_バイト列.Length; i++)
            {
                for (var l_反復回数 = 0; l_反復回数 < p_深さ; l_反復回数++)
                {
                    p_インデックス.V_登録(l_バイト列.AsSpan(i, k長));
                }
            }
        }

        #endregion

    }
}
