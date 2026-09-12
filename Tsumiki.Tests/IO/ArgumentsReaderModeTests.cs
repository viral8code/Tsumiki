using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// -mode (積極性のプリセット) が -pu/-pc を正しく束ねること、および後ろに書いた個別指定で上書きできることの検証
    /// </summary>
    public class ArgumentsReaderModeTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 実在するだけのダミーのリードのパス
        /// </summary>
        private readonly string _ダミーリードパス;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public ArgumentsReaderModeTests()
        {
            this._ダミーリードパス = Path.GetTempFileName();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (File.Exists(this._ダミーリードパス))
            {
                File.Delete(this._ダミーリードパス);
            }
        }

        /// <summary>
        /// -mode に保守的を指定すると、両方の閾値が保守的モードの値へ引き上がることを検証する
        /// </summary>
        [Fact]
        public void V_modeに保守的を指定すると両方の閾値が上がる()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-mode", Consts.積極性モード名.保守的]);

            Assert.Equal(Consts.保守的モードのペア結合閾値, l_引数.A_ペア結合閾値);
            Assert.Equal(Consts.保守的モードのペア支持数閾値, l_引数.A_ペア支持数閾値);
        }

        /// <summary>
        /// -mode に積極的を指定すると、両方の閾値が積極的モードの値へ引き下がることを検証する
        /// </summary>
        [Fact]
        public void V_modeに積極的を指定すると両方の閾値が下がる()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-mode", Consts.積極性モード名.積極的]);

            Assert.Equal(Consts.積極的モードのペア結合閾値, l_引数.A_ペア結合閾値);
            Assert.Equal(Consts.積極的モードのペア支持数閾値, l_引数.A_ペア支持数閾値);
        }

        /// <summary>
        /// -mode に標準を指定すると、素の既定値と一致することを検証する
        /// </summary>
        [Fact]
        public void V_modeに標準を指定すると素の既定値と一致する()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-mode", Consts.積極性モード名.標準]);

            Assert.Equal(Consts.ペア結合閾値の既定値, l_引数.A_ペア結合閾値);
            Assert.Equal(Consts.ペア支持数閾値の既定値, l_引数.A_ペア支持数閾値);
        }

        /// <summary>
        /// -mode の後ろに個別の -pu/-pc を書けば、そちらで上書きできる (通常の CLI 引数と同じく、後に書いたものが勝つ)
        /// </summary>
        [Fact]
        public void V_modeの後ろに書いた個別指定が優先される()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-mode", Consts.積極性モード名.保守的, "-pu", "0.5"]);

            Assert.Equal(0.5M, l_引数.A_ペア結合閾値);
            // -pc は指定していないので保守的モードの値のまま
            Assert.Equal(Consts.保守的モードのペア支持数閾値, l_引数.A_ペア支持数閾値);
        }

        #endregion

    }
}
