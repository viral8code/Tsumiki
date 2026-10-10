using Tsumiki.Models.Correction;

namespace Tsumiki.Tests.Models.Correction
{
    /// <summary>
    /// 誤りの無い区間の数え方と、そこからのカバレッジの見積もりを検証する
    /// </summary>
    public class 無誤り区間の度数Tests : IDisposable
    {
        #region 定数

        /// <summary>
        /// ファイル名 none
        /// </summary>
        private const string C_ファイル名_none = ".none";

        #endregion

        #region 内部変数

        /// <summary>
        /// 書き出し先
        /// </summary>
        private readonly string _パス = Path.GetTempFileName();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 信頼できる窓の連続を、窓の数 + k - 1 の長さの区間として数える
        /// </summary>
        [Fact]
        public void V_信頼できる窓の連続を区間として数える()
        {
            var l_度数 = new 無誤り区間の度数(3)
            {
                A_単一コピー平均 = 10D
            };
            l_度数.V_追加([true, true, false, true, true, true ], false);
            l_度数.V_追加([true, true, false, true, true, true ], true);
            Assert.Equal(10D, l_度数.Get_予測カバレッジ(3), 6);
            Assert.Equal(2D, l_度数.Get_予測カバレッジ(5), 6);
        }

        /// <summary>
        /// 訂正で取れる k-mer が増えた分だけ、訂正前のカバレッジより高く見積もる
        /// </summary>
        [Fact]
        public void V_訂正で増えたkmerの分だけ高く見積もる()
        {
            var l_度数 = new 無誤り区間の度数(3)
            {
                A_単一コピー平均 = 10D
            };
            l_度数.V_追加([true, true, false, true, true, true ], false);
            l_度数.V_追加([true, true, true, true, true, true ], true);
            Assert.Equal(12D, l_度数.Get_予測カバレッジ(3), 6);
        }

        /// <summary>
        /// 書き出したものを読み込むと同じ見積もりになる
        /// </summary>
        [Fact]
        public void V_書き出して読み込むと同じ見積もりになる()
        {
            var l_度数 = new 無誤り区間の度数(3)
            {
                A_単一コピー平均 = 26.5D
            };
            l_度数.V_追加([true, false, true, true ], false);
            l_度数.V_追加([true, true, true, true ], true);
            l_度数.V_書き出し(this._パス);
            var l_読込 = 無誤り区間の度数.Get_読込(this._パス);
            Assert.NotNull(l_読込);
            Assert.Equal(3, l_読込.A_k長);
            Assert.Equal(l_度数.Get_予測カバレッジ(4), l_読込.Get_予測カバレッジ(4), 9);
        }

        /// <summary>
        /// ファイルが無ければ null を返す
        /// </summary>
        [Fact]
        public void V_ファイルが無ければnullを返す()
        {
            Assert.Null(無誤り区間の度数.Get_読込(this._パス + C_ファイル名_none));
        }

        /// <summary>
        /// 書き出し先を片付ける
        /// </summary>
        public void Dispose()
        {
            File.Delete(this._パス);
            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
