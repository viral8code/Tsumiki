using System.Reflection;
using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Evidence;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 局所組み立ての長さの余裕幅の引数を省いたときの既定値と、既定値が従来どおり 30 として働くことの検証
    /// </summary>
    public class LocalAssemblerMarginTests
    {
        #region コンストラクタ

        /// <summary>
        /// 実行時引数を既定に戻す
        /// </summary>
        public LocalAssemblerMarginTests()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = 21,
                A_スレッド数 = 1
            };
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 余裕幅の引数は省略可能で、既定値は従来の 30 である
        /// </summary>
        [Fact]
        public void 余裕幅の既定値は従来どおり30()
        {
            Assert.Equal(30, Consts.ギャップ充填の長さの余裕幅);
            foreach (var l_名前 in new[] { "Get_固定kの局所結果", "Get_適応kの局所結果" })
            {
                var l_メソッド = typeof(LocalAssembler).GetMethod(l_名前, BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic);
                Assert.NotNull(l_メソッド);
                var l_最後 = l_メソッド!.GetParameters().Last();
                Assert.Equal("p_長さの余裕幅", l_最後.Name);
                Assert.True(l_最後.HasDefaultValue);
                Assert.Equal(30, l_最後.DefaultValue);
            }
        }

        /// <summary>
        /// 引数を省いた結果と、30 を明示した結果が同じになる
        /// </summary>
        [Fact]
        public void 余裕幅を省いた結果と30を明示した結果が同じになる()
        {
            const int l_基準k = 37;
            var l_左 = V_生成_ランダム配列(300, 501);
            var l_充填 = V_生成_ランダム配列(90, 502);
            var l_右 = V_生成_ランダム配列(300, 503);
            var l_正解 = l_左 + l_充填 + l_右;
            var l_局所リード = new List<読取証拠>();
            for (var i = l_左.Length - 35; i <= l_左.Length + l_充填.Length - 20; i += 5)
            {
                l_局所リード.Add(new 読取証拠(l_正解.Substring(i, 50), string.Empty));
            }

            var l_ギャップ = new 局所ギャップ(0, l_左.Length, l_充填.Length, l_左, l_右);
            var l_省略 = LocalAssembler.Get_適応kの局所結果(l_ギャップ, l_局所リード, l_基準k, out var l_省略判定);
            var l_明示 = LocalAssembler.Get_適応kの局所結果(l_ギャップ, l_局所リード, l_基準k, out var l_明示判定, 30);

            Assert.Equal(l_明示, l_省略);
            Assert.Equal(l_明示判定, l_省略判定);
            Assert.Equal(l_充填, l_省略);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 乱数の種を固定した ACGT の配列を作る
        /// </summary>
        /// <param name="p_長さ">長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>配列</returns>
        private static string V_生成_ランダム配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
        }

        #endregion
    }
}
