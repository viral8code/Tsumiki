using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Preprocessing;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 3' 末端の品質トリムの検証
    /// </summary>
    /// <remarks>
    /// 末尾の崩れた区間だけを切り、品質の良い区間の孤立した低品質塩基では切らないことを固定する
    /// </remarks>
    public class QualityTrimTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// Phred+33
        /// </summary>
        private const int オフセット = 33;

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
        public QualityTrimTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_quality_trim_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 末尾に続く低品質区間だけを切る
        /// </summary>
        /// <param name="p_品質">各塩基の Phred スコア (空白区切り)</param>
        /// <param name="p_残す長さ">残るはずの長さ</param>
        [Theory]
        [InlineData("40 40 40 40 40 40", 6)]
        [InlineData("40 40 40 40 2 2", 4)]
        [InlineData("40 40 40 5 40 40", 6)]
        [InlineData("40 40 40 40 5 30 2 2", 4)]
        [InlineData("40 40 40 40 10 30 2 2", 6)]
        [InlineData("2 2 2 2", 1)]
        public void 末尾の低品質区間だけを切る(string p_品質, int p_残す長さ)
        {
            var l_クオリティ = new string([.. p_品質.Split(' ').Select(x => (char)(int.Parse(x) + オフセット))]);
            Assert.Equal(p_残す長さ, Preprocessor.Get_品質トリム後の長さ(l_クオリティ, オフセット, 20));
        }

        /// <summary>
        /// 閾値 0 では何も切らない
        /// </summary>
        [Fact]
        public void 閾値0では切らない()
        {
            var l_結果 = new ペア前処理結果("ACGT", "!!!!", "TTTT", "!!!!", false, 0);
            Assert.Equal(l_結果, Preprocessor.Get_品質トリム済み(l_結果, オフセット, 0));
        }

        /// <summary>
        /// 前処理を通すと両側の末尾が切れ、切った塩基数を数える
        /// </summary>
        /// <remarks>
        /// 重なりの無い (断片がリード 2 本より長い) ペアでも切れることを確かめる
        /// </remarks>
        [Fact]
        public void 前処理を通すと両側の末尾が切れる()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_スレッド数 = 2, A_品質トリム閾値 = 20 };
            var l_乱数 = new Random(2501);
            var l_配列1 = new string([.. Enumerable.Range(0, 100).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_配列2 = new string([.. Enumerable.Range(0, 100).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_品質1 = new string('I', 80) + new string('#', 20);
            var l_品質2 = new string('I', 90) + new string('#', 10);
            var l_入力1 = Path.Combine(this._作業ディレクトリ, "in_1.fq");
            var l_入力2 = Path.Combine(this._作業ディレクトリ, "in_2.fq");
            File.WriteAllText(l_入力1, $"@r1/1\n{l_配列1}\n+\n{l_品質1}\n");
            File.WriteAllText(l_入力2, $"@r1/2\n{l_配列2}\n+\n{l_品質2}\n");
            var l_出力1 = Path.Combine(this._作業ディレクトリ, "out_1.fq");
            var l_出力2 = Path.Combine(this._作業ディレクトリ, "out_2.fq");

            var l_統計 = Preprocessor.V_前処理_リードファイル(l_入力1, l_入力2, l_出力1, l_出力2, オフセット);

            using var l_読み込み1 = new FastqReader(l_出力1);
            using var l_読み込み2 = new FastqReader(l_出力2);
            var (_, l_出た配列1, _) = l_読み込み1.Get_次のレコード();
            var (_, l_出た配列2, _) = l_読み込み2.Get_次のレコード();
            Assert.Equal(l_配列1[..80], l_出た配列1);
            Assert.Equal(l_配列2[..90], l_出た配列2);
            Assert.Equal(30, l_統計.A_品質トリム塩基数);
            Assert.Equal(200, l_統計.A_総塩基数);
        }

        #endregion
    }
}
