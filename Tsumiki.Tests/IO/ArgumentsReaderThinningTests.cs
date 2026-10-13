using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// -cov (目標の深度) の読み取り、既定値、0 の意味、負の値の拒否、設定表示への反映の検証
    /// </summary>
    public class ArgumentsReaderThinningTests : IDisposable
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
        public ArgumentsReaderThinningTests()
        {
            this._ダミーリードパス = Path.GetTempFileName();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// -cov 150 は目標の深度 150 として読まれる
        /// </summary>
        [Fact]
        public void V_covは目標の深度として読まれる()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, this._ダミーリードパス, Consts.引数キー.目標深度, "150"]);
            Assert.Equal(150, l_引数.A_目標深度);
        }

        /// <summary>
        /// -cov を書かなければ既定値の 100 になる
        /// </summary>
        [Fact]
        public void V_covを省略すると既定値100になる()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, this._ダミーリードパス]);
            Assert.Equal(Consts.目標深度の既定値, l_引数.A_目標深度);
            Assert.Equal(100, l_引数.A_目標深度);
        }

        /// <summary>
        /// -cov 0 は間引かない指定として 0 になる
        /// </summary>
        [Fact]
        public void V_cov0は間引かない指定になる()
        {
            var l_引数 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, this._ダミーリードパス, Consts.引数キー.目標深度, "0"]);
            Assert.Equal(0, l_引数.A_目標深度);
        }

        /// <summary>
        /// 負の値は受け付けない
        /// </summary>
        [Fact]
        public void V_cov負の値は例外になる()
        {
            Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, this._ダミーリードパス, Consts.引数キー.目標深度, "-1"]));
        }

        /// <summary>
        /// 設定表示には目標の深度が現れ、値が変わると表示も変わる
        /// </summary>
        [Fact]
        public void V_設定表示に目標の深度が現れる()
        {
            var l_既定 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, this._ダミーリードパス]);
            var l_変更 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, this._ダミーリードパス, Consts.引数キー.目標深度, "150"]);
            Assert.Contains("target depth", l_既定.ToString(), StringComparison.Ordinal);
            Assert.Contains("target depth (0 = keep all) : 100", l_既定.ToString(), StringComparison.Ordinal);
            Assert.Contains("target depth (0 = keep all) : 150", l_変更.ToString(), StringComparison.Ordinal);
        }

        /// <summary>
        /// 一時ファイルを片付ける
        /// </summary>
        public void Dispose()
        {
            if (File.Exists(this._ダミーリードパス))
            {
                File.Delete(this._ダミーリードパス);
            }
        }

        #endregion
    }
}
