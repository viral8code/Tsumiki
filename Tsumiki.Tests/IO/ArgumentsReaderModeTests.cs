using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.UnitigBuilding;

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
        /// 入力だけの CLI が標準機能を有効にすることを検証する
        /// </summary>
        [Fact]
        public void V_入力だけで標準プロファイルを使用する()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス]);
            Assert.True(l_引数.A_Is前処理);
            Assert.True(l_引数.A_Isエラー訂正);
            Assert.True(l_引数.A_Isマルチk);
            Assert.True(l_引数.A_Is引き継ぎ);
            Assert.True(l_引数.A_IsSuperRead作成);
            Assert.True(l_引数.A_Is反復rMer検証);
            Assert.True(l_引数.A_Is局所アセンブリ);
            Assert.True(l_引数.A_IsGFA出力);
            Assert.True(l_引数.A_Isポリッシュ);
            Assert.True(l_引数.A_Is環状閉鎖検証);
            Assert.True(l_引数.A_Is救済kmer使用);
            Assert.False(l_引数.A_Isマージ);
            Assert.True(l_引数.A_Is低カバレッジ端トリミング);
            Assert.Equal(コピー数基準の出所.Weighted, l_引数.A_コピー数基準の出所);
        }

        /// <summary>
        /// 従来プロファイルと後続の個別指定を検証する
        /// </summary>
        [Fact]
        public void V_従来プロファイルに個別設定を追加できる()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-profile", "legacy", "-cnb", "weighted", "-nt"]);
            Assert.False(l_引数.A_Is前処理);
            Assert.False(l_引数.A_Isエラー訂正);
            Assert.False(l_引数.A_Isマルチk);
            Assert.False(l_引数.A_IsSuperRead作成);
            Assert.False(l_引数.A_Isポリッシュ);
            Assert.Equal(コピー数基準の出所.Weighted, l_引数.A_コピー数基準の出所);
            Assert.False(l_引数.A_Is低カバレッジ端トリミング);
            Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-profile", "invalid"]));
        }

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

        [Fact]
        public void V_cnbとntは解析され複製と設定表示に残る()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-cnb", "weighted", "-nt"]);
            var l_複製 = l_引数.Get_複製();

            Assert.Equal(コピー数基準の出所.Weighted, l_複製.A_コピー数基準の出所);
            Assert.False(l_複製.A_Is低カバレッジ端トリミング);
            Assert.Contains("copy-number baseline requested : Weighted", l_複製.ToString());
            Assert.Contains("trim low-coverage graph ends : False", l_複製.ToString());
            Assert.Contains("-cnb", HelpText.Get_ヘルプ());
            Assert.Contains("-nt", HelpText.Get_ヘルプ());
        }

        #endregion

    }
}
