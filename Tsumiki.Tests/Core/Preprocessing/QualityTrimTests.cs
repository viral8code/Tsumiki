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
    public class QualityTrimTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki quality trim tests
        /// </summary>
        private const string C_項目_tsumiki_quality_trim_tests = "tsumiki_quality_trim_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目
        /// </summary>
        private const string C_項目 = "!!!!";

        /// <summary>
        /// 塩基配列 TTTT
        /// </summary>
        private const string C_塩基配列_TTTT = "TTTT";

        /// <summary>
        /// ファイル名 in 1 fq
        /// </summary>
        private const string C_ファイル名_in_1_fq = "in_1.fq";

        /// <summary>
        /// ファイル名 in 2 fq
        /// </summary>
        private const string C_ファイル名_in_2_fq = "in_2.fq";

        /// <summary>
        /// ファイル名 out 1 fq
        /// </summary>
        private const string C_ファイル名_out_1_fq = "out_1.fq";

        /// <summary>
        /// ファイル名 out 2 fq
        /// </summary>
        private const string C_ファイル名_out_2_fq = "out_2.fq";

        /// <summary>
        /// Phred+33
        /// </summary>
        private const int C_オフセット = 33;

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
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_quality_trim_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        public void V_末尾の低品質区間だけを切る(string p_品質, int p_残す長さ)
        {
            var l_クオリティ = new string([..p_品質.Split(' ').Select(x => (char)(int.Parse(x) + C_オフセット))]);
            Assert.Equal(p_残す長さ, Preprocessor.Get_品質トリム後の長さ(l_クオリティ, C_オフセット, 20));
        }

        /// <summary>
        /// 閾値 0 では何も切らない
        /// </summary>
        [Fact]
        public void V_閾値0では切らない()
        {
            var l_結果 = new ペア前処理結果(Consts.塩基文字, C_項目, C_塩基配列_TTTT, C_項目, false, 0);
            Assert.Equal(l_結果, Preprocessor.Get_品質トリム済み(l_結果, C_オフセット, 0));
        }

        /// <summary>
        /// 前処理を通すと両側の末尾が切れ、切った塩基数を数える
        /// </summary>
        [Fact]
        public void V_前処理を通すと両側の末尾が切れる()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_スレッド数 = 2,
                A_品質トリム閾値 = 20
            };
            var l_乱数 = new Random(2501);
            var l_基準配列 = new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_比較配列 = new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_品質1 = new string('I', 80) + new string('#', 20);
            var l_品質2 = new string('I', 90) + new string('#', 10);
            var l_入力1 = Path.Combine(this._作業ディレクトリ, C_ファイル名_in_1_fq);
            var l_入力2 = Path.Combine(this._作業ディレクトリ, C_ファイル名_in_2_fq);
            File.WriteAllText(l_入力1, $"@r1/1\n{l_基準配列}\n+\n{l_品質1}\n");
            File.WriteAllText(l_入力2, $"@r1/2\n{l_比較配列}\n+\n{l_品質2}\n");
            var l_出力1 = Path.Combine(this._作業ディレクトリ, C_ファイル名_out_1_fq);
            var l_出力2 = Path.Combine(this._作業ディレクトリ, C_ファイル名_out_2_fq);
            var l_統計 = Preprocessor.V_前処理_リードファイル(l_入力1, l_入力2, l_出力1, l_出力2, C_オフセット);
            using var l_読み込み1 = new FastqReader(l_出力1);
            using var l_読み込み2 = new FastqReader(l_出力2);
            var (_, l_出た基準配列, _) = l_読み込み1.Get_次のレコード();
            var (_, l_出た比較配列, _) = l_読み込み2.Get_次のレコード();
            Assert.Equal(l_基準配列[..80], l_出た基準配列);
            Assert.Equal(l_比較配列[..90], l_出た比較配列);
            Assert.Equal(30L, l_統計.A_品質トリム塩基数);
            Assert.Equal(200L, l_統計.A_総塩基数);
        }

        #endregion
    }
}
