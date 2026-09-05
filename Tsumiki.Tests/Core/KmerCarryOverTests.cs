using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 前段の k の配列を次の k へ引き継ぐ処理の検証。
    ///
    /// 引き継ぎで最も壊れやすいのはカバレッジである。名目値で埋めると
    /// コピー数推定・低カバレッジ端のトリミング・自己検査がまとめて狂う。
    /// 連結が保たれることと同じくらい、カバレッジが保たれることを固定する。
    /// </summary>
    public class KmerCarryOverTests : IDisposable
    {
        private readonly string _tempDir;

        public KmerCarryOverTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_carryover_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        private TrustedKmerIndex BuildIndex(int kmerLength, int depth, params string[] sequences)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = kmerLength, A_スレッド数 = 1 };
            // インデックスごとに作業ディレクトリを分ける(同じ場所を使うと
            // 一時ファイルの後始末が互いに干渉する)。
            var l_作業 = Path.Combine(this._tempDir, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業);
            var index = new TrustedKmerIndex(l_作業);
            foreach (var seq in sequences)
            {
                var bytes = seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, kmerLength), p_ワーカー番号: 0);
                    }
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);
            return index;
        }

        private string WriteFasta(string name, params string[] sequences)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new FastaWriter(path);
            var id = 1;
            foreach (var seq in sequences)
            {
                writer.V_書き込み($"NODE{id++}", seq);
            }
            return path;
        }

        [Fact]
        public void Prepare_RecordsTheCoverageOfEveryKmerInTheSequence()
        {
            const int k = 21;
            var truth = RandomSequence(2_000, seed: 701);
            using var index = this.BuildIndex(k, depth: 17, truth);

            var path = this.WriteFasta("prep.fasta", truth);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(path, index, k);

            var l_項目 = Assert.Single(l_引き継ぎ);
            Assert.Equal(truth, l_項目.A_配列);
            Assert.Equal(k, l_項目.A_k長);
            Assert.Equal(truth.Length - k + 1, l_項目.A_カバレッジ.Length);
            Assert.All(l_項目.A_カバレッジ, x => Assert.Equal(17, x));
        }

        /// <summary>
        /// 引き継いだ k-mer のカバレッジが、k の差ぶんスケールされること。
        /// スケールしないと、引き継いだ領域だけカバレッジが高く見えて
        /// コピー数を過大に推定する。
        /// </summary>
        [Fact]
        public void CarriedCoverage_IsScaledByTheNumberOfKmersPerRead()
        {
            // 前段 k=21 でカバレッジ 100、リード長 150。
            // 次段 k=101 では 1リードあたり 130 本から 50 本へ減るので、
            // 100 * 50 / 130 = 38。
            var l_引き継ぎ = new 引き継ぎ配列(
                RandomSequence(300, seed: 702), [.. Enumerable.Repeat(100, 280)], 21);

            var l_カバレッジ = KmerCarryOver.Get_引き継ぐカバレッジ(
                l_引き継ぎ, p_位置: 0, p_k長: 101, p_リード長: 150);

            Assert.Equal((ulong)Math.Round(100.0 * 50 / 130), l_カバレッジ);
        }

        /// <summary>
        /// 引き継ぐ k-mer は、それを構成する前段の k-mer の最小値を超えないこと。
        /// 長い k-mer は短い k-mer をすべて含むので、最も弱い部分より強くはなれない。
        /// </summary>
        [Fact]
        public void CarriedCoverage_TakesTheWeakestConstituentKmer()
        {
            var l_カバレッジ = new int[280];
            Array.Fill(l_カバレッジ, 100);
            l_カバレッジ[10] = 7;

            var l_引き継ぎ = new 引き継ぎ配列(RandomSequence(300, seed: 703), l_カバレッジ, 21);

            // 位置 0 から k=41 の窓は前段の位置 0..20 を含むので、7 が効く。
            var l_弱い部分を含む = KmerCarryOver.Get_引き継ぐカバレッジ(
                l_引き継ぎ, p_位置: 0, p_k長: 41, p_リード長: null);
            // 位置 30 の窓は 30..50 なので 7 を含まない。
            var l_含まない = KmerCarryOver.Get_引き継ぐカバレッジ(
                l_引き継ぎ, p_位置: 30, p_k長: 41, p_リード長: null);

            Assert.Equal(7UL, l_弱い部分を含む);
            Assert.Equal(100UL, l_含まない);
        }

        /// <summary>
        /// 引き継ぎの本題。カバレッジが薄くて次の k では観測されなかった領域が、
        /// 前段の配列から復元されること。
        /// </summary>
        [Fact]
        public void CarryOver_RestoresKmersThatTheLargerKDidNotObserve()
        {
            const int 前段のk = 21;
            const int 次のk = 41;
            var truth = RandomSequence(3_000, seed: 711);

            // 次の k では中央の領域が観測されていない状況を作る。
            var l_左 = truth[..1_200];
            var l_右 = truth[1_800..];
            using var l_次段 = this.BuildIndex(次のk, depth: 20, l_左, l_右);

            // 前段は全体を観測している。
            using var l_前段 = this.BuildIndex(前段のk, depth: 20, truth);
            var l_パス = this.WriteFasta("carry.fasta", truth);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, 前段のk);

            var l_塩基列 = truth.Select(Util.Get_塩基ID).ToArray();
            var l_中央 = l_塩基列.AsSpan(1_400, 次のk);
            Assert.False(l_次段.Get_含まれるか(l_中央));

            var l_追加数 = KmerCarryOver.V_引き継ぎ(l_引き継ぎ, l_次段, 次のk, p_リード長: null);

            Assert.True(l_追加数 > 0);
            Assert.True(l_次段.Get_含まれるか(l_中央));
            Assert.True(l_次段.Get_カバレッジ(l_中央) > 0);
        }

        /// <summary>
        /// 既に観測されている k-mer のカバレッジは書き換えないこと。
        /// 実際のリードによる観測のほうが、前段からの推定より確かである。
        /// </summary>
        [Fact]
        public void CarryOver_DoesNotOverwriteCoverageThatWasActuallyObserved()
        {
            const int 前段のk = 21;
            const int 次のk = 31;
            var truth = RandomSequence(2_000, seed: 721);

            using var l_次段 = this.BuildIndex(次のk, depth: 40, truth);
            using var l_前段 = this.BuildIndex(前段のk, depth: 5, truth);

            var l_塩基列 = truth.Select(Util.Get_塩基ID).ToArray();
            var l_観測済み = l_塩基列.AsSpan(100, 次のk);
            var l_元のカバレッジ = l_次段.Get_カバレッジ(l_観測済み);

            var l_パス = this.WriteFasta("nooverwrite.fasta", truth);
            _ = KmerCarryOver.V_引き継ぎ(
                KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, 前段のk), l_次段, 次のk, p_リード長: null);

            Assert.Equal(l_元のカバレッジ, l_次段.Get_カバレッジ(l_観測済み));
        }

        /// <summary>
        /// 2コピーの反復配列は、引き継いでも2コピー相当のカバレッジを保つこと。
        /// ここが崩れるとコピー数推定が壊れ、反復配列の扱いが総崩れになる。
        /// </summary>
        [Fact]
        public void CarryOver_PreservesTheRelativeCoverageOfRepeats()
        {
            const int 前段のk = 21;
            const int 次のk = 41;
            var 単一 = RandomSequence(1_500, seed: 731);
            var 反復 = RandomSequence(1_500, seed: 732);

            // 反復側だけ倍の深さで観測されている状況。
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 前段のk, A_スレッド数 = 1 };
            var l_前段の作業 = Path.Combine(this._tempDir, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_前段の作業);
            using var l_前段 = new TrustedKmerIndex(l_前段の作業);
            foreach (var (l_配列, l_深さ) in new[] { (単一, 20), (反復, 40) })
            {
                var l_bytes = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + 前段のk <= l_bytes.Length; i++)
                {
                    for (var rep = 0; rep < l_深さ; rep++)
                    {
                        l_前段.V_登録(l_bytes.AsSpan(i, 前段のk), p_ワーカー番号: 0);
                    }
                }
            }
            _ = l_前段.V_カットオフ(p_カットオフ: 2);

            // 次の k は何も観測していない。
            using var l_次段 = this.BuildIndex(次のk, depth: 20, RandomSequence(2_000, seed: 733));

            var l_パス = this.WriteFasta("repeat_carry.fasta", 単一, 反復);
            _ = KmerCarryOver.V_引き継ぎ(
                KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, 前段のk), l_次段, 次のk, p_リード長: null);

            var l_単一の位置 = 単一.Select(Util.Get_塩基ID).ToArray().AsSpan(50, 次のk);
            var l_反復の位置 = 反復.Select(Util.Get_塩基ID).ToArray().AsSpan(50, 次のk);

            // 2倍の関係が保たれていること。
            Assert.Equal(2 * l_次段.Get_カバレッジ(l_単一の位置), l_次段.Get_カバレッジ(l_反復の位置));
        }

        /// <summary>
        /// 短い断片は引き継がないこと。連結の役に立たないうえ、
        /// エラー由来の残骸である可能性が相対的に高い。
        /// </summary>
        [Fact]
        public void Prepare_SkipsShortSequences()
        {
            const int k = 21;
            var 長い = RandomSequence(2_000, seed: 741);
            var 短い = RandomSequence(120, seed: 742);
            using var index = this.BuildIndex(k, depth: 20, 長い, 短い);

            var path = this.WriteFasta("short.fasta", 長い, 短い);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(path, index, k);

            var l_項目 = Assert.Single(l_引き継ぎ);
            Assert.Equal(長い, l_項目.A_配列);
        }
    }
}
