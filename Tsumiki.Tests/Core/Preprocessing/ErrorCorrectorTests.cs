using Tsumiki.Common;
using Tsumiki.Core.Preprocessing;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ErrorCorrector によるリードの誤り訂正の検証
    /// </summary>
    public class ErrorCorrectorTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public ErrorCorrectorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_error_corrector_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] Get_塩基列(string p_配列)
        {
            return [.. p_配列.Select(c => c switch
            {
                'A' => Consts.塩基ID.A,
                'C' => Consts.塩基ID.C,
                'G' => Consts.塩基ID.G,
                'T' => Consts.塩基ID.T,
                'N' => Consts.無効な塩基,
                _ => throw new InvalidOperationException(),
            })];
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
        /// "true" 配列の全 k-mer (順鎖・逆鎖) をカットオフ以上登録した
        /// TrustedKmerIndex を構築する
        /// </summary>
        /// <param name="p_真の配列">元になる配列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_スレッド数">使うスレッド数</param>
        /// <returns>構築した信頼できる k-mer インデックス</returns>
        private TrustedKmerIndex V_構築_信頼できるインデックス(string p_真の配列, int p_k長, int p_スレッド数 = 1)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = p_スレッド数 };
            var l_インデックス = new TrustedKmerIndex(this._tempDir);
            var l_塩基列 = Get_塩基列(p_真の配列);
            for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
            {
                // カットオフ (2) を超えるよう複数回登録する
                for (var l_回 = 0; l_回 < 3; l_回++)
                {
                    l_インデックス.V_登録(l_塩基列.AsSpan(i, p_k長), p_ワーカー番号: 0);
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2);
            return l_インデックス;
        }

        /// <summary>
        /// 単一の置換エラーを真の配列へ訂正することを検証する
        /// </summary>
        [Fact]
        public void 単一の置換エラーを真の配列へ訂正する()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 49 bp
            const int k = 15;
            using var index = this.V_構築_信頼できるインデックス(trueSeq, k);

            var mutated = trueSeq.ToCharArray();
            mutated[20] = mutated[20] == 'A' ? 'C' : 'A'; // 真の配列と異なる塩基に置換
            var mutatedBytes = Get_塩基列(new string(mutated));

            var result = ErrorCorrector.Get_訂正結果(mutatedBytes, index, k);

            Assert.Equal(trueSeq, Get_配列(result.A_塩基列));
            Assert.Equal(1, result.A_訂正数);
        }

        /// <summary>
        /// エラーが無ければリードを変更せず、訂正数もゼロになることを検証する
        /// </summary>
        [Fact]
        public void エラーが無ければ変更せず訂正数もゼロになる()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.V_構築_信頼できるインデックス(trueSeq, k);

            var result = ErrorCorrector.Get_訂正結果(Get_塩基列(trueSeq), index, k);

            Assert.Equal(trueSeq, Get_配列(result.A_塩基列));
            Assert.Equal(0, result.A_訂正数);
        }

        /// <summary>
        /// 十分に離れた 2 箇所のエラーを両方とも訂正することを検証する
        /// </summary>
        [Fact]
        public void 十分離れた2箇所のエラーを両方とも訂正する()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 70 bp
            const int k = 15;
            using var index = this.V_構築_信頼できるインデックス(trueSeq, k);

            var mutated = trueSeq.ToCharArray();
            mutated[10] = mutated[10] == 'A' ? 'G' : 'A';
            mutated[55] = mutated[55] == 'A' ? 'G' : 'A';
            var mutatedBytes = Get_塩基列(new string(mutated));

            var result = ErrorCorrector.Get_訂正結果(mutatedBytes, index, k);

            Assert.Equal(trueSeq, Get_配列(result.A_塩基列));
            Assert.Equal(2, result.A_訂正数);
        }

        /// <summary>
        /// k 長より短いリードはそのまま返すことを検証する
        /// </summary>
        [Fact]
        public void k長より短いリードはそのまま返す()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.V_構築_信頼できるインデックス(trueSeq, k);

            var shortRead = Get_塩基列("ACGT");
            var result = ErrorCorrector.Get_訂正結果(shortRead, index, k);

            Assert.Equal("ACGT", Get_配列(result.A_塩基列));
            Assert.Equal(0, result.A_訂正数);
        }

        /// <summary>
        /// 無効な塩基の位置は書き換えないことを検証する
        /// </summary>
        [Fact]
        public void 無効な塩基の位置は書き換えない()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.V_構築_信頼できるインデックス(trueSeq, k);

            var withN = trueSeq.ToCharArray();
            withN[20] = 'N';
            var bytesWithN = Get_塩基列(new string(withN));

            var result = ErrorCorrector.Get_訂正結果(bytesWithN, index, k);

            Assert.Equal(Consts.無効な塩基, result.A_塩基列[20]);
        }

        /// <summary>
        /// 訂正処理が入力配列そのものを変更しないことを検証する
        /// </summary>
        [Fact]
        public void 入力配列を変更しない()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.V_構築_信頼できるインデックス(trueSeq, k);

            var mutated = trueSeq.ToCharArray();
            mutated[20] = mutated[20] == 'A' ? 'C' : 'A';
            var mutatedBytes = Get_塩基列(new string(mutated));
            var original = (byte[])mutatedBytes.Clone();

            _ = ErrorCorrector.Get_訂正結果(mutatedBytes, index, k);

            Assert.Equal(original, mutatedBytes);
        }
    }
}
