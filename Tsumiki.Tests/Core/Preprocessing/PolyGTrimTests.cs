using Tsumiki.Cores.Preprocessing;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 3' 末端の poly-G を削る処理 (Preprocessor.Get_polyG末尾の長さ) の検証
    /// </summary>
    public class PolyGTrimTests
    {
        #region 定数

        /// <summary>
        /// poly-G の付かない配列
        /// </summary>
        private const string C_配列 = "ACGTTGCATGCAACGTAGCTAGCATCGATCGATTACA";

        /// <summary>
        /// Phred オフセット
        /// </summary>
        private const int C_Phredオフセット = 33;

        /// <summary>
        /// 作業ディレクトリの名前の頭
        /// </summary>
        private const string C_作業の頭 = "tsumiki_polyg_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 in_1.fq
        /// </summary>
        private const string C_ファイル名_in_1_fq = "in_1.fq";

        /// <summary>
        /// ファイル名 in_2.fq
        /// </summary>
        private const string C_ファイル名_in_2_fq = "in_2.fq";

        /// <summary>
        /// ファイル名 out_1.fq
        /// </summary>
        private const string C_ファイル名_out_1_fq = "out_1.fq";

        /// <summary>
        /// ファイル名 out_2.fq
        /// </summary>
        private const string C_ファイル名_out_2_fq = "out_2.fq";

        /// <summary>
        /// 読み 2 の配列
        /// </summary>
        private const string C_読み2 = "TTGACCATGGCATTAGGCATCAGGATCCATGCAATGCC";

        /// <summary>
        /// poly-G の付いた読み 1 の末尾
        /// </summary>
        private const string C_polyG25 = "GGGGGGGGGGGGGGGGGGGGGGGGG";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 3' 末端の G の連続の長さに応じて、削る長さが決まる (10 以上だけ削る)
        /// </summary>
        /// <param name="p_末尾">配列の末尾に足す配列</param>
        /// <param name="p_期待値">削る長さ</param>
        [Theory]
        [InlineData("GGGGGGGGGGGG", 12)]
        [InlineData("GGGGGGGGG", 0)]
        [InlineData("GGGGGAGGGGGG", 12)]
        [InlineData("", 0)]
        [InlineData("GGGGGGGGGGGGT", 13)]
        public void V_末尾のpolyGの長さを求める(string p_末尾, int p_期待値)
        {
            Assert.Equal(p_期待値, Preprocessor.Get_polyG末尾の長さ(C_配列 + p_末尾));
        }

        /// <summary>
        /// 全体が G のリードは 1 塩基を残して削る
        /// </summary>
        [Fact]
        public void V_全体がGなら1塩基残す()
        {
            Assert.Equal(150, Preprocessor.Get_polyG末尾の長さ(new string('G', 151)));
        }

        /// <summary>
        /// 前処理のファイル処理で、poly-G を削った配列とクオリティが同じ長さになり、統計に数えられる
        /// </summary>
        [Fact]
        public void V_前処理で配列とクオリティを同じ長さで削る()
        {
            var l_作業 = Path.Combine(Path.GetTempPath(), C_作業の頭 + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(l_作業);
            try
            {
                var l_入力1 = Path.Combine(l_作業, C_ファイル名_in_1_fq);
                var l_入力2 = Path.Combine(l_作業, C_ファイル名_in_2_fq);
                var l_出力1 = Path.Combine(l_作業, C_ファイル名_out_1_fq);
                var l_出力2 = Path.Combine(l_作業, C_ファイル名_out_2_fq);
                var l_読み1 = C_配列 + C_polyG25;
                var l_読み2 = C_読み2;
                File.WriteAllText(l_入力1, $"@r/1\n{l_読み1}\n+\n{new string('I', l_読み1.Length)}\n");
                File.WriteAllText(l_入力2, $"@r/2\n{l_読み2}\n+\n{new string('I', l_読み2.Length)}\n");

                var l_統計 = Preprocessor.V_前処理_リードファイル(l_入力1, l_入力2, l_出力1, l_出力2, C_Phredオフセット);

                var l_行 = File.ReadAllLines(l_出力1);
                Assert.Equal(C_配列, l_行[1]);
                Assert.Equal(C_配列.Length, l_行[3].Length);
                Assert.Equal(1, l_統計.A_polyGを削ったリード数);
                Assert.Equal(25L, l_統計.A_polyGで削った塩基数);
            }
            finally
            {
                Directory.Delete(l_作業, recursive: true);
            }
        }

        #endregion
    }
}
