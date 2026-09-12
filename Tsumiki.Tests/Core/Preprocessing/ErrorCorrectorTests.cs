using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ErrorCorrector によるリードの誤り訂正の検証
    /// </summary>
    public class ErrorCorrectorTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public ErrorCorrectorTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_error_corrector_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 単一の置換エラーを真の配列へ訂正することを検証する
        /// </summary>
        [Fact]
        public void V_単一の置換エラーを真の配列へ訂正する()
        {
            const string l_正解配列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 49 bp
            const int l_k長 = 15;
            using var l_インデックス = this.V_構築_信頼できるインデックス(l_正解配列, l_k長);

            var l_変異配列 = l_正解配列.ToCharArray();
            l_変異配列[20] = l_変異配列[20] == 'A' ? 'C' : 'A'; // 真の配列と異なる塩基に置換
            var l_変異塩基列 = Get_塩基列(new string(l_変異配列));

            var l_結果 = ErrorCorrector.Get_訂正結果(l_変異塩基列, l_インデックス, l_k長);

            Assert.Equal(l_正解配列, Get_配列(l_結果.A_塩基列));
            Assert.Equal(1, l_結果.A_訂正数);
        }

        /// <summary>
        /// エラーが無ければリードを変更せず、訂正数もゼロになることを検証する
        /// </summary>
        [Fact]
        public void V_エラーが無ければ変更せず訂正数もゼロになる()
        {
            const string l_正解配列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int l_k長 = 15;
            using var l_インデックス = this.V_構築_信頼できるインデックス(l_正解配列, l_k長);

            var l_結果 = ErrorCorrector.Get_訂正結果(Get_塩基列(l_正解配列), l_インデックス, l_k長);

            Assert.Equal(l_正解配列, Get_配列(l_結果.A_塩基列));
            Assert.Equal(0, l_結果.A_訂正数);
        }

        /// <summary>
        /// 十分に離れた 2 箇所のエラーを両方とも訂正することを検証する
        /// </summary>
        [Fact]
        public void V_十分離れた2箇所のエラーを両方とも訂正する()
        {
            const string l_正解配列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 70 bp
            const int l_k長 = 15;
            using var l_インデックス = this.V_構築_信頼できるインデックス(l_正解配列, l_k長);

            var l_変異配列 = l_正解配列.ToCharArray();
            l_変異配列[10] = l_変異配列[10] == 'A' ? 'G' : 'A';
            l_変異配列[55] = l_変異配列[55] == 'A' ? 'G' : 'A';
            var l_変異塩基列 = Get_塩基列(new string(l_変異配列));

            var l_結果 = ErrorCorrector.Get_訂正結果(l_変異塩基列, l_インデックス, l_k長);

            Assert.Equal(l_正解配列, Get_配列(l_結果.A_塩基列));
            Assert.Equal(2, l_結果.A_訂正数);
        }

        /// <summary>
        /// k 長より短いリードはそのまま返すことを検証する
        /// </summary>
        [Fact]
        public void V_k長より短いリードはそのまま返す()
        {
            const string l_正解配列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int l_k長 = 15;
            using var l_インデックス = this.V_構築_信頼できるインデックス(l_正解配列, l_k長);

            var l_短いリード = Get_塩基列("ACGT");
            var l_結果 = ErrorCorrector.Get_訂正結果(l_短いリード, l_インデックス, l_k長);

            Assert.Equal("ACGT", Get_配列(l_結果.A_塩基列));
            Assert.Equal(0, l_結果.A_訂正数);
        }

        /// <summary>
        /// 無効な塩基の位置は書き換えないことを検証する
        /// </summary>
        [Fact]
        public void V_無効な塩基の位置は書き換えない()
        {
            const string l_正解配列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int l_k長 = 15;
            using var l_インデックス = this.V_構築_信頼できるインデックス(l_正解配列, l_k長);

            var l_曖昧塩基付き配列 = l_正解配列.ToCharArray();
            l_曖昧塩基付き配列[20] = 'N';
            var l_曖昧塩基付きバイト列 = Get_塩基列(new string(l_曖昧塩基付き配列));

            var l_結果 = ErrorCorrector.Get_訂正結果(l_曖昧塩基付きバイト列, l_インデックス, l_k長);

            Assert.Equal(Consts.無効な塩基, l_結果.A_塩基列[20]);
        }

        /// <summary>
        /// 訂正処理が入力配列そのものを変更しないことを検証する
        /// </summary>
        [Fact]
        public void V_入力配列を変更しない()
        {
            const string l_正解配列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int l_k長 = 15;
            using var l_インデックス = this.V_構築_信頼できるインデックス(l_正解配列, l_k長);

            var l_変異配列 = l_正解配列.ToCharArray();
            l_変異配列[20] = l_変異配列[20] == 'A' ? 'C' : 'A';
            var l_変異塩基列 = Get_塩基列(new string(l_変異配列));
            var l_元配列 = (byte[])l_変異塩基列.Clone();

            _ = ErrorCorrector.Get_訂正結果(l_変異塩基列, l_インデックス, l_k長);

            Assert.Equal(l_元配列, l_変異塩基列);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] Get_塩基列(string p_配列)
        {
            return [.. p_配列.Select(c => c switch { 'A' => Consts.塩基ID.A, 'C' => Consts.塩基ID.C, 'G' => Consts.塩基ID.G, 'T' => Consts.塩基ID.T, 'N' => Consts.無効な塩基, _ => throw new InvalidOperationException(), })];
        }

        /// <summary>
        /// 塩基 ID 列を配列へ変換する
        /// </summary>
        /// <param name="p_塩基列">元の塩基 ID 列</param>
        /// <returns>配列</returns>
        private static string Get_配列(byte[] p_塩基列)
        {
            return string.Join(string.Empty, p_塩基列.Select(Util.V_変換_塩基文字));
        }

        /// <summary>
        /// "true" 配列の全 k-mer (順鎖・逆鎖) をカットオフ以上登録した TrustedKmerIndex を構築する
        /// </summary>
        /// <param name="p_真の配列">元になる配列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_スレッド数">使うスレッド数</param>
        /// <returns>構築した信頼できる k-mer インデックス</returns>
        private TrustedKmerIndex V_構築_信頼できるインデックス(string p_真の配列, int p_k長, int p_スレッド数 = 1)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = p_スレッド数 };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_塩基列 = Get_塩基列(p_真の配列);
            for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
            {
                // カットオフ (2) を超えるよう複数回登録する
                for (var l_回 = 0; l_回 < 3; l_回++)
                {
                    l_インデックス.V_登録(l_塩基列.AsSpan(i, p_k長));
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        #endregion

    }
}
