using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
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
        /// 既定で有効な処理は、それぞれの -n で始まるフラグで 1 つずつ切れる
        /// </summary>
        /// <param name="p_キー">切るフラグ</param>
        /// <param name="p_項目">切れるはずの設定の名前</param>
        [Theory]
        [InlineData("-npp", nameof(Parameters.A_Is前処理))]
        [InlineData("-nec", nameof(Parameters.A_Isエラー訂正))]
        [InlineData("-nmk", nameof(Parameters.A_Isマルチk))]
        [InlineData("-nsr", nameof(Parameters.A_IsSuperRead作成))]
        [InlineData("-nrv", nameof(Parameters.A_Is反復rMer検証))]
        [InlineData("-nla", nameof(Parameters.A_Is局所アセンブリ))]
        [InlineData("-ngfa", nameof(Parameters.A_IsGFA出力))]
        [InlineData("-npo", nameof(Parameters.A_Isポリッシュ))]
        [InlineData("-ncc", nameof(Parameters.A_Is環状閉鎖検証))]
        [InlineData("-nmy", nameof(Parameters.A_Is救済kmer使用))]
        [InlineData("-nc", nameof(Parameters.A_Is引き継ぎ))]
        [InlineData("-nt", nameof(Parameters.A_Is低カバレッジ端トリミング))]
        public void 既定で有効な処理を1つずつ切れる(string p_キー, string p_項目)
        {
            string[] l_全項目 =
            [
                nameof(Parameters.A_Is前処理), nameof(Parameters.A_Isエラー訂正), nameof(Parameters.A_Isマルチk), nameof(Parameters.A_IsSuperRead作成),
                nameof(Parameters.A_Is反復rMer検証), nameof(Parameters.A_Is局所アセンブリ), nameof(Parameters.A_IsGFA出力), nameof(Parameters.A_Isポリッシュ),
                nameof(Parameters.A_Is環状閉鎖検証), nameof(Parameters.A_Is救済kmer使用), nameof(Parameters.A_Is引き継ぎ), nameof(Parameters.A_Is低カバレッジ端トリミング),
            ];

            var l_引数 = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, p_キー]);

            foreach (var l_項目 in l_全項目)
            {
                var l_値 = (bool)typeof(Parameters).GetProperty(l_項目)!.GetValue(l_引数)!;
                Assert.Equal(l_項目 != p_項目, l_値);
            }
        }

        /// <summary>
        /// 廃止した -profile と、既定で有効になった処理を有効にするだけのフラグは受け付けない
        /// </summary>
        /// <param name="p_引数">渡す引数 (空白区切り)</param>
        [Theory]
        [InlineData("-profile legacy")]
        [InlineData("-po")]
        [InlineData("-mk")]
        public void 廃止したフラグは受け付けない(string p_引数)
        {
            Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, .. p_引数.Split(' ')]));
        }

        /// <summary>
        /// 複数の -k は試して選ぶ指定なので、マルチ k を切る指定とは併用できない
        /// </summary>
        [Fact]
        public void 複数のkとマルチkなしは併用できない()
        {
            var l_例外 = Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-nmk", "-k", "31,63"]));
            Assert.Contains("-nmk", l_例外.Message, StringComparison.Ordinal);

            var l_単一k = ArgumentsReader.Get_実行時引数(["-1", this._ダミーリードパス, "-nmk", "-k", "31"]);
            Assert.False(l_単一k.A_Isマルチk);
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
