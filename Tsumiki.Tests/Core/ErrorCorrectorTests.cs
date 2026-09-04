using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    public class ErrorCorrectorTests : IDisposable
    {
        private readonly string _tempDir;

        public ErrorCorrectorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_error_corrector_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private static byte[] ToBytes(string seq)
        {
            return [.. seq.Select(c => c switch
            {
                'A' => Consts.塩基ID.A,
                'C' => Consts.塩基ID.C,
                'G' => Consts.塩基ID.G,
                'T' => Consts.塩基ID.T,
                'N' => Consts.無効な塩基,
                _ => throw new InvalidOperationException(),
            })];
        }

        private static string ToSeq(byte[] bytes)
        {
            return string.Join(string.Empty, bytes.Select(Util.V_変換_塩基文字));
        }

        /// <summary>
        /// "true"配列の全k-mer(順鎖・逆鎖)をカットオフ以上登録した
        /// TrustedKmerIndex を構築する。
        /// </summary>
        private TrustedKmerIndex BuildTrustedIndex(string trueSeq, int kmerLength, int threadCount = 1)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = kmerLength, A_スレッド数 = threadCount };
            var index = new TrustedKmerIndex(this._tempDir);
            var bytes = ToBytes(trueSeq);
            for (var i = 0; i + kmerLength <= bytes.Length; i++)
            {
                // カットオフ(2)を超えるよう複数回登録する。
                for (var rep = 0; rep < 3; rep++)
                {
                    index.V_登録(bytes.AsSpan(i, kmerLength), p_ワーカー番号: 0);
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);
            return index;
        }

        [Fact]
        public void CorrectRead_FixesSingleSubstitutionError_BackToTrueSequence()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 49bp
            const int k = 15;
            using var index = this.BuildTrustedIndex(trueSeq, k);

            var mutated = trueSeq.ToCharArray();
            mutated[20] = mutated[20] == 'A' ? 'C' : 'A'; // 真の配列と異なる塩基に置換
            var mutatedBytes = ToBytes(new string(mutated));

            var result = ErrorCorrector.Get_訂正結果(mutatedBytes, index, k);

            Assert.Equal(trueSeq, ToSeq(result.A_塩基列));
            Assert.Equal(1, result.A_訂正数);
        }

        [Fact]
        public void CorrectRead_NoErrors_LeavesReadUnchangedWithZeroCorrections()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.BuildTrustedIndex(trueSeq, k);

            var result = ErrorCorrector.Get_訂正結果(ToBytes(trueSeq), index, k);

            Assert.Equal(trueSeq, ToSeq(result.A_塩基列));
            Assert.Equal(0, result.A_訂正数);
        }

        [Fact]
        public void CorrectRead_TwoWellSeparatedErrors_FixesBoth()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 70bp
            const int k = 15;
            using var index = this.BuildTrustedIndex(trueSeq, k);

            var mutated = trueSeq.ToCharArray();
            mutated[10] = mutated[10] == 'A' ? 'G' : 'A';
            mutated[55] = mutated[55] == 'A' ? 'G' : 'A';
            var mutatedBytes = ToBytes(new string(mutated));

            var result = ErrorCorrector.Get_訂正結果(mutatedBytes, index, k);

            Assert.Equal(trueSeq, ToSeq(result.A_塩基列));
            Assert.Equal(2, result.A_訂正数);
        }

        [Fact]
        public void CorrectRead_ShorterThanKmerLength_ReturnsUnchanged()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.BuildTrustedIndex(trueSeq, k);

            var shortRead = ToBytes("ACGT");
            var result = ErrorCorrector.Get_訂正結果(shortRead, index, k);

            Assert.Equal("ACGT", ToSeq(result.A_塩基列));
            Assert.Equal(0, result.A_訂正数);
        }

        [Fact]
        public void CorrectRead_NeverModifiesInvalidBasePositions()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.BuildTrustedIndex(trueSeq, k);

            var withN = trueSeq.ToCharArray();
            withN[20] = 'N';
            var bytesWithN = ToBytes(new string(withN));

            var result = ErrorCorrector.Get_訂正結果(bytesWithN, index, k);

            Assert.Equal(Consts.無効な塩基, result.A_塩基列[20]);
        }

        [Fact]
        public void CorrectRead_DoesNotMutateInputArray()
        {
            const string trueSeq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            const int k = 15;
            using var index = this.BuildTrustedIndex(trueSeq, k);

            var mutated = trueSeq.ToCharArray();
            mutated[20] = mutated[20] == 'A' ? 'C' : 'A';
            var mutatedBytes = ToBytes(new string(mutated));
            var original = (byte[])mutatedBytes.Clone();

            _ = ErrorCorrector.Get_訂正結果(mutatedBytes, index, k);

            Assert.Equal(original, mutatedBytes);
        }
    }
}
