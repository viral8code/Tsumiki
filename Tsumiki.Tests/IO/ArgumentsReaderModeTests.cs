using Tsumiki.Common;
using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// -mode(積極性のプリセット) が -pu/-pc を正しく束ねること、
    /// および後ろに書いた個別指定で上書きできることの検証
    /// </summary>
    public class ArgumentsReaderModeTests : IDisposable
    {
        private readonly string _dummyReadPath;

        public ArgumentsReaderModeTests()
        {
            this._dummyReadPath = Path.GetTempFileName();
        }

        public void Dispose()
        {
            if (File.Exists(this._dummyReadPath))
            {
                File.Delete(this._dummyReadPath);
            }
        }

        [Fact]
        public void Mode_Conservative_RaisesBothThresholds()
        {
            var args = ArgumentsReader.Get_実行時引数(
                ["-1", this._dummyReadPath, "-mode", Consts.積極性モード名.保守的]);

            Assert.Equal(Consts.保守的モードのペア結合閾値, args.A_ペア結合閾値);
            Assert.Equal(Consts.保守的モードのペア支持数閾値, args.A_ペア支持数閾値);
        }

        [Fact]
        public void Mode_Bold_LowersBothThresholds()
        {
            var args = ArgumentsReader.Get_実行時引数(
                ["-1", this._dummyReadPath, "-mode", Consts.積極性モード名.積極的]);

            Assert.Equal(Consts.積極的モードのペア結合閾値, args.A_ペア結合閾値);
            Assert.Equal(Consts.積極的モードのペア支持数閾値, args.A_ペア支持数閾値);
        }

        [Fact]
        public void Mode_Normal_MatchesThePlainDefaults()
        {
            var args = ArgumentsReader.Get_実行時引数(
                ["-1", this._dummyReadPath, "-mode", Consts.積極性モード名.標準]);

            Assert.Equal(Consts.ペア結合閾値の既定値, args.A_ペア結合閾値);
            Assert.Equal(Consts.ペア支持数閾値の既定値, args.A_ペア支持数閾値);
        }

        /// <summary>
        /// -mode の後ろに個別の -pu/-pc を書けば、そちらで上書きできる
        /// (通常の CLI 引数と同じく、後に書いたものが勝つ)
        /// </summary>
        [Fact]
        public void Mode_FollowedByExplicitThreshold_TheExplicitOneWins()
        {
            var args = ArgumentsReader.Get_実行時引数(
                ["-1", this._dummyReadPath, "-mode", Consts.積極性モード名.保守的, "-pu", "0.5"]);

            Assert.Equal(0.5M, args.A_ペア結合閾値);
            // -pc は指定していないので保守的モードの値のまま
            Assert.Equal(Consts.保守的モードのペア支持数閾値, args.A_ペア支持数閾値);
        }
    }
}
