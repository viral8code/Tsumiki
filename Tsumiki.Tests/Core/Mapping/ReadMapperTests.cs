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
        #region 定数

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        #endregion

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
            Assert.False(l_配置.A_Is逆鎖);
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
            Assert.True(l_配置.A_Is逆鎖);
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
            var l_リード = l_元のリード.Insert(50, C_項目_A);
            var l_マッパー = new ReadMapper([l_参照]);
            var l_配置 = l_マッパー.Get_配置(l_リード);
            Assert.Equal(0, l_配置.A_配列番号);
            Assert.True(l_配置.A_信頼度 > 0);
            Assert.Equal(100, l_配置.A_整列位置群.Count);
        }

        /// <summary>
        /// 参照の末端からはみ出すリードでも例外を出さずに配置できることを確かめる
        /// </summary>
        [Fact]
        public void Get_配置_参照末端からはみ出すリードを配置する()
        {
            var l_参照 = Get_乱数配列(500, 4);
            var l_リード = l_参照.Substring(450, 50) + Get_乱数配列(100, 5);
            var l_マッパー = new ReadMapper([l_参照]);
            var l_配置 = l_マッパー.Get_配置(l_リード);
            Assert.True(l_配置.A_整列位置群.Count == 0 || l_配置.A_整列位置群[0].A_参照位置 >= 440);
        }

        /// <summary>
        /// 断片長の上限より離れた 2 コピーの反復に入った片方も、相方の近くのコピーに置いてペアとして一意に決まることを確かめる
        /// </summary>
        [Fact]
        public void Get_組んだ配置_反復に入った相方を近くのコピーに置く()
        {
            var l_反復 = Get_乱数配列(150, 11);
            var l_参照 = Get_乱数配列(1_000, 12) + l_反復 + Get_乱数配列(6_000, 13) + l_反復 + Get_乱数配列(1_000, 14);
            var l_マッパー = new ReadMapper([l_参照]);
            var l_順リード = l_参照.Substring(800, 100);
            var l_逆リード = Util.V_逆相補(l_参照.Substring(1_020, 100));
            Assert.Equal(0, ReadMapper.Get_最良の配置(l_マッパー.Get_配置候補群(l_逆リード)).A_信頼度);
            var l_組 = ReadMapper.Get_組んだ配置(l_マッパー.Get_配置候補群(l_順リード), l_マッパー.Get_配置候補群(l_逆リード), 5_000);
            Assert.NotNull(l_組);
            Assert.Equal(1_020, l_組.Value.A_配置2.A_整列位置群[0].A_参照位置);
            Assert.True(l_組.Value.A_配置1.A_信頼度 >= 20);
            Assert.True(l_組.Value.A_配置2.A_信頼度 >= 20);
        }

        /// <summary>
        /// 相方が単独ではずっと良く当たる別の場所を持つなら、得点を大きく譲ってまでペアに組まないことを確かめる
        /// </summary>
        [Fact]
        public void Get_組んだ配置_相方の本来の場所より大きく劣る組は採らない()
        {
            var l_左 = Get_乱数配列(1_000, 21);
            var l_右 = Get_乱数配列(1_000, 22);
            var l_別の場所 = Get_乱数配列(1_000, 23);
            var l_マッパー = new ReadMapper([l_左 + l_右, l_別の場所]);
            var l_順リード = l_左.Substring(700, 100);
            var l_逆リード = Util.V_逆相補(l_右.Substring(0, 70) + l_別の場所.Substring(500, 80));
            Assert.Contains(l_マッパー.Get_配置候補群(l_逆リード), x => x.A_配列番号 == 0 && x.A_スコア >= 30);
            Assert.Null(ReadMapper.Get_組んだ配置(l_マッパー.Get_配置候補群(l_順リード), l_マッパー.Get_配置候補群(l_逆リード), 5_000));
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
            const string l_塩基 = Consts.塩基文字;
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => l_塩基[l_乱数.Next(4)]));
        }

        #endregion
    }
}
