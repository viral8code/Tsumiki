using Tsumiki.Core.Preprocessing;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// R1 と RC(R2) の重なり解析による、アダプタリードスルーのトリムと
    /// ペア相互訂正 (fastp 型の前処理) を固定する
    /// </summary>
    public class PreprocessorTests
    {
        private const int Phredオフセット = 33;

        private static string 高品質クオリティ(int p_長さ, int p_スコア = 35)
        {
            return new string((char)(p_スコア + Phredオフセット), p_長さ);
        }

        private static string 一様クオリティ(int p_長さ, int p_スコア)
        {
            return new string((char)(p_スコア + Phredオフセット), p_長さ);
        }

        private static string 位置だけ変更したクオリティ(int p_長さ, int p_基本スコア, int p_位置, int p_その位置のスコア)
        {
            var l_文字 = 一様クオリティ(p_長さ, p_基本スコア).ToCharArray();
            l_文字[p_位置] = (char)(p_その位置のスコア + Phredオフセット);
            return new string(l_文字);
        }

        [Fact]
        public void Get_前処理結果_フラグメント長がリード長を超える通常のペアは変更しない()
        {
            // 40 bp を大きく超える、互いに無関係な60bpの配列 (=重なりが存在しない、典型的なケース)
            const string 配列1 = "ACGTGCATTGCAGTCAGGCTAACGGTTCCAAGGTCATGCAGTACGGCTTA";
            const string 配列2 = "TTGGCACCGGATCCAAGTTGGCCAATGGCTTACGATCGTAGGCCTTAACG";
            var l_クオリティ1 = 高品質クオリティ(配列1.Length);
            var l_クオリティ2 = 高品質クオリティ(配列2.Length);

            var l_結果 = Preprocessor.Get_前処理結果(配列1, l_クオリティ1, 配列2, l_クオリティ2, Phredオフセット);

            Assert.Equal(配列1, l_結果.A_配列1);
            Assert.Equal(配列2, l_結果.A_配列2);
            Assert.False(l_結果.A_アダプタを検出したか);
            Assert.Equal(0, l_結果.A_訂正塩基数);
        }

        [Fact]
        public void Get_前処理結果_アダプタリードスルーを検出してフラグメント長までトリムする()
        {
            // read1 = 真の断片 + アダプタ、read2 = RC(真の断片) + 別のアダプタ、という
            // 構成で「フラグメント長がリード長より短く、両端から重なって読んだ」
            // 状況を作る
            // 繰り返し配列だと周期性でオフセットが一意に定まらないため、
            // 非周期な配列を使う
            const string 真の断片 = "TGACCTGAAGCTTAGGCATCGGTAACCTTGGACGTCAGTA"; // 40bp, non-repetitive
            const string アダプタ1 = "AGATCGGAAG";
            const string アダプタ2 = "TTTTTTTTTT";
            var 配列1 = 真の断片 + アダプタ1; // 50bp
            var 配列2 = Tsumiki.Common.Util.V_逆相補(真の断片) + アダプタ2; // 50bp

            var l_クオリティ1 = 高品質クオリティ(配列1.Length);
            var l_クオリティ2 = 高品質クオリティ(配列2.Length);

            var l_結果 = Preprocessor.Get_前処理結果(配列1, l_クオリティ1, 配列2, l_クオリティ2, Phredオフセット);

            Assert.True(l_結果.A_アダプタを検出したか);
            Assert.Equal(真の断片, l_結果.A_配列1);
            Assert.Equal(Tsumiki.Common.Util.V_逆相補(真の断片), l_結果.A_配列2);
            Assert.Equal(真の断片.Length, l_結果.A_クオリティ1.Length);
            Assert.Equal(真の断片.Length, l_結果.A_クオリティ2.Length);
        }

        [Fact]
        public void Get_前処理結果_一方が高信頼で他方が低信頼のときだけ低信頼側を上書きする()
        {
            const string 真の断片 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 48bp、自己RC対称
            const int エラー位置 = 5;
            var 誤った塩基 = 真の断片[エラー位置] == 'A' ? 'C' : 'A';

            var l_配列1 = 真の断片.ToCharArray();
            l_配列1[エラー位置] = 誤った塩基;
            var 配列1 = new string(l_配列1);
            var 配列2 = 真の断片; // read2 は誤りなし

            // read1 のエラー位置だけ低信頼、read2 は対応する位置 (逆向きなので末尾側) だけ高信頼にする
            // 他は両方とも高信頼にしておき、他の位置で誤訂正が起きないことも確認する
            var l_クオリティ1 = 位置だけ変更したクオリティ(配列1.Length, p_基本スコア: 35, エラー位置, p_その位置のスコア: 5);
            var l_対応する読み2の位置 = 配列2.Length - 1 - エラー位置;
            var l_クオリティ2 = 位置だけ変更したクオリティ(配列2.Length, p_基本スコア: 35, l_対応する読み2の位置, p_その位置のスコア: 35);

            var l_結果 = Preprocessor.Get_前処理結果(配列1, l_クオリティ1, 配列2, l_クオリティ2, Phredオフセット);

            Assert.False(l_結果.A_アダプタを検出したか);
            Assert.Equal(1, l_結果.A_訂正塩基数);
            Assert.Equal(真の断片, l_結果.A_配列1);
            Assert.Equal(配列2, l_結果.A_配列2);
        }

        [Fact]
        public void Get_前処理結果_曖昧塩基を含む位置は書き換えない()
        {
            const string 真の断片 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 48bp
            const int N位置 = 5;

            var l_配列1 = 真の断片.ToCharArray();
            l_配列1[N位置] = 'N';
            var 配列1 = new string(l_配列1);
            var 配列2 = 真の断片;

            var l_クオリティ1 = 位置だけ変更したクオリティ(配列1.Length, p_基本スコア: 35, N位置, p_その位置のスコア: 2);
            var l_対応する読み2の位置 = 配列2.Length - 1 - N位置;
            var l_クオリティ2 = 位置だけ変更したクオリティ(配列2.Length, p_基本スコア: 35, l_対応する読み2の位置, p_その位置のスコア: 35);

            var l_結果 = Preprocessor.Get_前処理結果(配列1, l_クオリティ1, 配列2, l_クオリティ2, Phredオフセット);

            Assert.Equal(0, l_結果.A_訂正塩基数);
            Assert.Equal('N', l_結果.A_配列1[N位置]);
        }
    }
}
