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
        /// Phred オフセット
        /// </summary>
        private const int Phredオフセット = 33;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// フラグメント長がリード長を超える通常のペアは変更しないことを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_フラグメント長がリード長を超える通常のペアは変更しない()
        {
            // 40 bp を大きく超える、互いに無関係な 60 bp の配列 (=重なりが存在しない、典型的なケース)
            const string l_配列1 = "ACGTGCATTGCAGTCAGGCTAACGGTTCCAAGGTCATGCAGTACGGCTTA";
            const string l_配列2 = "TTGGCACCGGATCCAAGTTGGCCAATGGCTTACGATCGTAGGCCTTAACG";
            var l_クオリティ1 = V_高品質クオリティ(l_配列1.Length);
            var l_クオリティ2 = V_高品質クオリティ(l_配列2.Length);

            var l_結果 = Preprocessor.Get_前処理結果(l_配列1, l_クオリティ1, l_配列2, l_クオリティ2, Phredオフセット);

            Assert.Equal(l_配列1, l_結果.A_配列1);
            Assert.Equal(l_配列2, l_結果.A_配列2);
            Assert.False(l_結果.A_アダプタを検出したか);
            Assert.Equal(0, l_結果.A_訂正塩基数);
        }

        /// <summary>
        /// アダプタリードスルーを検出し、フラグメント長までトリムすることを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_アダプタリードスルーを検出してフラグメント長までトリムする()
        {
            // read1 = 真の断片 + アダプタ、read2 = RC (真の断片) + 別のアダプタ、という
            // 構成で「フラグメント長がリード長より短く、両端から重なって読んだ」
            // 状況を作る
            // 繰り返し配列だと周期性でオフセットが一意に定まらないため、
            // 非周期な配列を使う
            const string l_真の断片 = "TGACCTGAAGCTTAGGCATCGGTAACCTTGGACGTCAGTA"; // 40 bp, non-repetitive
            const string l_アダプタ1 = "AGATCGGAAG";
            const string l_アダプタ2 = "TTTTTTTTTT";
            var l_配列1 = l_真の断片 + l_アダプタ1; // 50 bp
            var l_配列2 = Tsumiki.Commons.Util.V_逆相補(l_真の断片) + l_アダプタ2; // 50 bp

            var l_クオリティ1 = V_高品質クオリティ(l_配列1.Length);
            var l_クオリティ2 = V_高品質クオリティ(l_配列2.Length);

            var l_結果 = Preprocessor.Get_前処理結果(l_配列1, l_クオリティ1, l_配列2, l_クオリティ2, Phredオフセット);

            Assert.True(l_結果.A_アダプタを検出したか);
            Assert.Equal(l_真の断片, l_結果.A_配列1);
            Assert.Equal(Tsumiki.Commons.Util.V_逆相補(l_真の断片), l_結果.A_配列2);
            Assert.Equal(l_真の断片.Length, l_結果.A_クオリティ1.Length);
            Assert.Equal(l_真の断片.Length, l_結果.A_クオリティ2.Length);
        }

        /// <summary>
        /// 一方が高信頼で他方が低信頼のときだけ、低信頼側を上書きすることを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_一方が高信頼で他方が低信頼のときだけ低信頼側を上書きする()
        {
            const string l_真の断片 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 48 bp、自己 RC 対称
            const int l_エラー位置 = 5;
            var l_誤った塩基 = l_真の断片[l_エラー位置] == 'A' ? 'C' : 'A';

            var l_変更塩基列 = l_真の断片.ToCharArray();
            l_変更塩基列[l_エラー位置] = l_誤った塩基;
            var l_配列1 = new string(l_変更塩基列);
            var l_配列2 = l_真の断片; // read2 は誤りなし

            // read1 のエラー位置だけ低信頼、read2 は対応する位置 (逆向きなので末尾側) だけ高信頼にする
            // 他は両方とも高信頼にしておき、他の位置で誤訂正が起きないことも確認する
            var l_クオリティ1 = V_位置だけ変更したクオリティ(l_配列1.Length, p_基本スコア: 35, l_エラー位置, p_その位置のスコア: 5);
            var l_対応する読み2の位置 = l_配列2.Length - 1 - l_エラー位置;
            var l_クオリティ2 = V_位置だけ変更したクオリティ(l_配列2.Length, p_基本スコア: 35, l_対応する読み2の位置, p_その位置のスコア: 35);

            var l_結果 = Preprocessor.Get_前処理結果(l_配列1, l_クオリティ1, l_配列2, l_クオリティ2, Phredオフセット);

            Assert.False(l_結果.A_アダプタを検出したか);
            Assert.Equal(1, l_結果.A_訂正塩基数);
            Assert.Equal(l_真の断片, l_結果.A_配列1);
            Assert.Equal(l_配列2, l_結果.A_配列2);
        }

        /// <summary>
        /// 曖昧塩基を含む位置は書き換えないことを検証する
        /// </summary>
        [Fact]
        public void Get_前処理結果_曖昧塩基を含む位置は書き換えない()
        {
            const string l_真の断片 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 48 bp
            const int l_N位置 = 5;

            var l_変更塩基列 = l_真の断片.ToCharArray();
            l_変更塩基列[l_N位置] = 'N';
            var l_配列1 = new string(l_変更塩基列);
            var l_配列2 = l_真の断片;

            var l_クオリティ1 = V_位置だけ変更したクオリティ(l_配列1.Length, p_基本スコア: 35, l_N位置, p_その位置のスコア: 2);
            var l_対応する読み2の位置 = l_配列2.Length - 1 - l_N位置;
            var l_クオリティ2 = V_位置だけ変更したクオリティ(l_配列2.Length, p_基本スコア: 35, l_対応する読み2の位置, p_その位置のスコア: 35);

            var l_結果 = Preprocessor.Get_前処理結果(l_配列1, l_クオリティ1, l_配列2, l_クオリティ2, Phredオフセット);

            Assert.Equal(0, l_結果.A_訂正塩基数);
            Assert.Equal('N', l_結果.A_配列1[l_N位置]);
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
            return new string((char)(p_スコア + Phredオフセット), p_長さ);
        }

        /// <summary>
        /// 全体が同じスコアのクオリティ文字列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_スコア">与える Phred スコア</param>
        /// <returns>クオリティ文字列</returns>
        private static string V_一様クオリティ(int p_長さ, int p_スコア)
        {
            return new string((char)(p_スコア + Phredオフセット), p_長さ);
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
            l_文字[p_位置] = (char)(p_その位置のスコア + Phredオフセット);
            return new string(l_文字);
        }

        #endregion

    }
}
