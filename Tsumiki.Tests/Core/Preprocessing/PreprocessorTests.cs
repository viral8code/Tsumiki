using Tsumiki.Cores.Preprocessing;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// R1 と RC (R2) の重なり解析による、アダプタリードスルーのトリムとペア相互訂正 (fastp 型の前処理) を固定する
    /// </summary>
    public class PreprocessorTests
    {
        #region 定数

        /// <summary>
        /// 塩基配列 ACGTGCATTGCAGTCAGGCTAACG CATGCAGTACGGCTTA
        /// </summary>
        private const string C_塩基配列_ACGTGCATTGCAGTCAGGCTAACG_CATGCAGTACGGCTTA = "ACGTGCATTGCAGTCAGGCTAACGGTTCCAAGGTCATGCAGTACGGCTTA";

        /// <summary>
        /// 塩基配列 TTGGCACCGGATCCAAGTTGGCCA ATCGTAGGCCTTAACG
        /// </summary>
        private const string C_塩基配列_TTGGCACCGGATCCAAGTTGGCCA_ATCGTAGGCCTTAACG = "TTGGCACCGGATCCAAGTTGGCCAATGGCTTACGATCGTAGGCCTTAACG";

        /// <summary>
        /// 塩基配列 TGACCTGAAGCTTAGGCATCGGTAACCTTGGACGTCAGTA
        /// </summary>
        private const string C_塩基配列_TGACCTGAAGCTTAGGCATCGGTAACCTTGGACGTCAGTA = "TGACCTGAAGCTTAGGCATCGGTAACCTTGGACGTCAGTA";

        /// <summary>
        /// 塩基配列 AGATCGGAAG
        /// </summary>
        private const string C_塩基配列_AGATCGGAAG = "AGATCGGAAG";

        /// <summary>
        /// 塩基配列 TTTTTTTTTT
        /// </summary>
        private const string C_塩基配列_TTTTTTTTTT = "TTTTTTTTTT";

        /// <summary>
        /// 塩基配列 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        /// <summary>
        /// Phred オフセット
        /// </summary>
        private const int C_Phredオフセット = 33;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// フラグメント長がリード長を超える通常のペアは変更しないことを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_フラグメント長がリード長を超える通常のペアは変更しない()
        {
            const string l_基準配列 = C_塩基配列_ACGTGCATTGCAGTCAGGCTAACG_CATGCAGTACGGCTTA;
            const string l_比較配列 = C_塩基配列_TTGGCACCGGATCCAAGTTGGCCA_ATCGTAGGCCTTAACG;
            var l_クオリティ1 = V_高品質クオリティ(l_基準配列.Length);
            var l_クオリティ2 = V_高品質クオリティ(l_比較配列.Length);
            var l_結果 = Preprocessor.Get_前処理結果(l_基準配列, l_クオリティ1, l_比較配列, l_クオリティ2, C_Phredオフセット);
            Assert.Equal(l_基準配列, l_結果.A_順リード配列);
            Assert.Equal(l_比較配列, l_結果.A_逆リード配列);
            Assert.False(l_結果.A_Hasアダプタ検出);
            Assert.Equal(0, l_結果.A_訂正塩基数);
        }

        /// <summary>
        /// アダプタリードスルーを検出し、フラグメント長までトリムすることを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_アダプタリードスルーを検出してフラグメント長までトリムする()
        {
            const string l_真の断片 = C_塩基配列_TGACCTGAAGCTTAGGCATCGGTAACCTTGGACGTCAGTA;
            const string l_アダプタ1 = C_塩基配列_AGATCGGAAG;
            const string l_アダプタ2 = C_塩基配列_TTTTTTTTTT;
            var l_基準配列 = l_真の断片 + l_アダプタ1;
            var l_比較配列 = Tsumiki.Commons.Util.V_逆相補(l_真の断片) + l_アダプタ2;
            var l_クオリティ1 = V_高品質クオリティ(l_基準配列.Length);
            var l_クオリティ2 = V_高品質クオリティ(l_比較配列.Length);
            var l_結果 = Preprocessor.Get_前処理結果(l_基準配列, l_クオリティ1, l_比較配列, l_クオリティ2, C_Phredオフセット);
            Assert.True(l_結果.A_Hasアダプタ検出);
            Assert.Equal(l_真の断片, l_結果.A_順リード配列);
            Assert.Equal(Tsumiki.Commons.Util.V_逆相補(l_真の断片), l_結果.A_逆リード配列);
            Assert.Equal(l_真の断片.Length, l_結果.A_順リード品質.Length);
            Assert.Equal(l_真の断片.Length, l_結果.A_逆リード品質.Length);
        }

        /// <summary>
        /// 一方が高信頼で他方が低信頼のときだけ、低信頼側を上書きすることを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_一方が高信頼で他方が低信頼のときだけ低信頼側を上書きする()
        {
            const string l_真の断片 = C_塩基配列_ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT;
            const int l_エラー位置 = 5;
            var l_誤った塩基 = l_真の断片[l_エラー位置] == 'A' ? 'C' : 'A';
            var l_変更塩基列 = l_真の断片.ToCharArray();
            l_変更塩基列[l_エラー位置] = l_誤った塩基;
            var l_基準配列 = new string(l_変更塩基列);
            var l_比較配列 = l_真の断片;
            var l_クオリティ1 = V_位置だけ変更したクオリティ(l_基準配列.Length, p_基本スコア: 35, l_エラー位置, p_その位置のスコア: 5);
            var l_対応する読み2の位置 = l_比較配列.Length - 1 - l_エラー位置;
            var l_クオリティ2 = V_位置だけ変更したクオリティ(l_比較配列.Length, p_基本スコア: 35, l_対応する読み2の位置, p_その位置のスコア: 35);
            var l_結果 = Preprocessor.Get_前処理結果(l_基準配列, l_クオリティ1, l_比較配列, l_クオリティ2, C_Phredオフセット);
            Assert.False(l_結果.A_Hasアダプタ検出);
            Assert.Equal(1, l_結果.A_訂正塩基数);
            Assert.Equal(l_真の断片, l_結果.A_順リード配列);
            Assert.Equal(l_比較配列, l_結果.A_逆リード配列);
        }

        /// <summary>
        /// 曖昧塩基を含む位置は書き換えないことを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_曖昧塩基を含む位置は書き換えない()
        {
            const string l_真の断片 = C_塩基配列_ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT;
            const int l_N位置 = 5;
            var l_変更塩基列 = l_真の断片.ToCharArray();
            l_変更塩基列[l_N位置] = 'N';
            var l_基準配列 = new string(l_変更塩基列);
            var l_比較配列 = l_真の断片;
            var l_クオリティ1 = V_位置だけ変更したクオリティ(l_基準配列.Length, p_基本スコア: 35, l_N位置, p_その位置のスコア: 2);
            var l_対応する読み2の位置 = l_比較配列.Length - 1 - l_N位置;
            var l_クオリティ2 = V_位置だけ変更したクオリティ(l_比較配列.Length, p_基本スコア: 35, l_対応する読み2の位置, p_その位置のスコア: 35);
            var l_結果 = Preprocessor.Get_前処理結果(l_基準配列, l_クオリティ1, l_比較配列, l_クオリティ2, C_Phredオフセット);
            Assert.Equal(0, l_結果.A_訂正塩基数);
            Assert.Equal('N', l_結果.A_順リード配列[l_N位置]);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 全体が高品質なクオリティ文字列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_スコア">与える Phred スコア</param>
        /// <returns>クオリティ文字列</returns>
        private static string V_高品質クオリティ(int p_長さ, int p_スコア = 35)
        {
            return new string((char)(p_スコア + C_Phredオフセット), p_長さ);
        }

        /// <summary>
        /// 全体が同じスコアのクオリティ文字列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_スコア">与える Phred スコア</param>
        /// <returns>クオリティ文字列</returns>
        private static string V_一様クオリティ(int p_長さ, int p_スコア)
        {
            return new string((char)(p_スコア + C_Phredオフセット), p_長さ);
        }

        /// <summary>
        /// 1 箇所だけスコアを変えたクオリティ文字列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_基本スコア">全体に与える Phred スコア</param>
        /// <param name="p_位置">変える位置</param>
        /// <param name="p_その位置のスコア">その位置に与える Phred スコア</param>
        /// <returns>クオリティ文字列</returns>
        private static string V_位置だけ変更したクオリティ(int p_長さ, int p_基本スコア, int p_位置, int p_その位置のスコア)
        {
            var l_文字 = V_一様クオリティ(p_長さ, p_基本スコア).ToCharArray();
            l_文字[p_位置] = (char)(p_その位置のスコア + C_Phredオフセット);
            return new string(l_文字);
        }

        #endregion
    }
}
