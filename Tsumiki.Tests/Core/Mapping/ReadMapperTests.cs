using Tsumiki.Commons;
using Tsumiki.Cores.Mapping;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core.Mapping
{
    /// <summary>
    /// ReadMapper の種チェインと整列を固定する
    /// </summary>
    public class ReadMapperTests
    {
        #region コンストラクタ

        /// <summary>
        /// テストに必要な共有状態を初期化する
        /// </summary>
        public ReadMapperTests()
        {
            ConfigurationManager.A_実行時引数 = new Parameters();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 順鎖のリードを元の位置へ配置できることを確かめる
        /// </summary>
        [Fact]
        public void Get_配置_順鎖のリードを配置する()
        {
            var l_参照 = Get_乱数配列(500, 1);
            var l_マッパー = new ReadMapper([l_参照]);
            var l_配置 = l_マッパー.Get_配置(l_参照.Substring(120, 100));

            Assert.Equal(0, l_配置.A_配列番号);
            Assert.False(l_配置.A_逆鎖か);
            Assert.True(l_配置.A_信頼度 > 0);
            Assert.Equal(100, l_配置.A_整列位置群.Count);
            Assert.Equal(120, l_配置.A_整列位置群[0].A_参照位置);
        }

        /// <summary>
        /// 逆相補リードを元の位置へ配置できることを確かめる
        /// </summary>
        [Fact]
        public void Get_配置_逆鎖のリードを配置する()
        {
            var l_参照 = Get_乱数配列(500, 2);
            var l_マッパー = new ReadMapper([l_参照]);
            var l_リード = Util.V_逆相補(l_参照.Substring(140, 100));
            var l_配置 = l_マッパー.Get_配置(l_リード);

            Assert.Equal(0, l_配置.A_配列番号);
            Assert.True(l_配置.A_逆鎖か);
            Assert.True(l_配置.A_信頼度 > 0);
            Assert.Equal(140, l_配置.A_整列位置群[0].A_参照位置);
        }

        /// <summary>
        /// 小さな挿入を含むリードでも対応する塩基を配置できることを確かめる
        /// </summary>
        [Fact]
        public void Get_配置_挿入を含むリードを配置する()
        {
            var l_参照 = Get_乱数配列(500, 3);
            var l_元のリード = l_参照.Substring(160, 100);
            var l_リード = l_元のリード.Insert(50, "A");
            var l_マッパー = new ReadMapper([l_参照]);
            var l_配置 = l_マッパー.Get_配置(l_リード);

            Assert.Equal(0, l_配置.A_配列番号);
            Assert.True(l_配置.A_信頼度 > 0);
            Assert.Equal(100, l_配置.A_整列位置群.Count);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ"></param>
        /// <param name="p_種"></param>
        /// <returns>塩基配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string l_塩基 = "ACGT";
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => l_塩基[l_乱数.Next(4)]));
        }

        #endregion
    }
}
