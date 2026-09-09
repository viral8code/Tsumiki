using Tsumiki.Common;
using Tsumiki.Core.Preprocessing;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// パック経路 (k &lt;= 64) と逐次経路が、同じリードに対して同じ訂正を返すことの確認<br/>
    /// パック経路は窓の評価を 2 bit 演算に置き換えた最適化なので、
    /// 結果が 1 文字でも違えば最適化が壊れている
    /// </summary>
    public class ErrorCorrectorPackedPathTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        public ErrorCorrectorPackedPathTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_ec_packed_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] Get_塩基列(string p_配列)
        {
            return [.. p_配列.Select(Util.Get_塩基ID)];
        }

        /// <summary>
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_真の配列">元になる配列</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex Get_インデックス(string p_真の配列, int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._一時ディレクトリ);
            var l_塩基列 = Get_塩基列(p_真の配列);
            for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
            {
                for (var l_回 = 0; l_回 < 3; l_回++)
                {
                    l_インデックス.V_登録(l_塩基列.AsSpan(i, p_k長), p_ワーカー番号: 0);
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2);
            return l_インデックス;
        }

        [Theory]
        [InlineData(21)]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(63)]
        public void パック経路は逐次経路と同じ訂正を返す(int p_k長)
        {
            var l_乱数 = new Random(p_k長);
            const string 塩基 = "ACGT";
            var l_真の配列 = string.Concat(
                Enumerable.Range(0, 4000).Select(_ => 塩基[l_乱数.Next(4)]));

            using var l_インデックス = this.Get_インデックス(l_真の配列, p_k長);

            for (var l_回 = 0; l_回 < 200; l_回++)
            {
                var l_開始 = l_乱数.Next(l_真の配列.Length - 150);
                var l_リード = l_真の配列.Substring(l_開始, 150).ToCharArray();

                // 置換エラーと曖昧塩基を混ぜる
                for (var l_誤り = 0; l_誤り < l_乱数.Next(4); l_誤り++)
                {
                    l_リード[l_乱数.Next(l_リード.Length)] = 塩基[l_乱数.Next(4)];
                }
                if (l_回 % 10 == 0)
                {
                    l_リード[l_乱数.Next(l_リード.Length)] = 'N';
                }

                var l_塩基列 = Get_塩基列(new string(l_リード));
                var l_パック = ErrorCorrector.Get_訂正結果(l_塩基列, l_インデックス, p_k長);
                var l_逐次 = ErrorCorrector.Get_訂正結果_逐次(
                    [.. l_塩基列], l_インデックス, p_k長, p_最大反復数: 10);

                Assert.Equal(l_逐次.A_訂正数, l_パック.A_訂正数);
                Assert.Equal(l_逐次.A_塩基列, l_パック.A_塩基列);
            }
        }
    }
}
