using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ContigMaker.FindDominantUnitig が返す LastMatchEndOffset が、
    /// 「read内での位置」ではなく「unitig内での正しい位置」を指すことを
    /// 直接検証する。以前は kmerDict が (unitigId のみ) しか保持しておらず、
    /// FindDominantUnitig は read 内の最終ヒットk-merの終端位置(read基準)を
    /// そのまま unitig 内終端位置として誤用していた。unitig が read より
    /// 十分短い間は両者がたまたま近い値になり問題が表面化しなかったが、
    /// tip clipping導入後にunitigが大幅に長くなり、インサートサイズ自動推定
    /// (この値を使う)が明後日の値(中央値30bp等)を返すようになったことで発覚した。
    /// </summary>
    public class ContigMakerFindDominantUnitigTests : IDisposable
    {
        private readonly string _tempDir;

        public ContigMakerFindDominantUnitigTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_contigmaker_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        // k=8で内部にk-mer重複のないことを確認済みの100bp配列。
        private const string UnitigSeq = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAG";

        private ContigMaker CreateContigMakerWithSingleUnitig(int kmerLength)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = kmerLength, ThreadCount = 1 };
            var unitigsPath = Path.Combine(this._tempDir, "unitigs.fasta");
            File.WriteAllText(unitigsPath, $">1\n{UnitigSeq}\n");
            return new ContigMaker(unitigsPath);
        }

        [Fact]
        public void FindDominantUnitig_ForwardMatch_ReturnsUnitigCoordinateEndOffset_NotReadCoordinate()
        {
            var contigMaker = this.CreateContigMakerWithSingleUnitig(kmerLength: 8);

            // read = unitig の [40,70) 部分(30bp)。read自身の長さ(30)ではなく、
            // unitig内での終端位置(70)が返るはず。
            var read = UnitigSeq.Substring(40, 30);

            var hit = contigMaker.FindDominantUnitig(read);

            Assert.Equal(1, hit.UnitigId); // 正の値 = 順鎖でのヒット
            Assert.Equal(70, hit.LastMatchEndOffset);
            Assert.Equal(UnitigSeq.Length, hit.UnitigLength);
        }

        [Fact]
        public void FindDominantUnitig_ReverseComplementMatch_ReturnsReverseOrientedUnitigCoordinate()
        {
            var contigMaker = this.CreateContigMakerWithSingleUnitig(kmerLength: 8);

            // 元の [40,70) を逆相補した read。unitig全体を逆相補した向きで見ると、
            // 元の区間 [40,70) は [100-70, 100-40) = [30,60) に写る。
            var read = Util.ReverseComprement(UnitigSeq.Substring(40, 30));

            var hit = contigMaker.FindDominantUnitig(read);

            Assert.Equal(-1, hit.UnitigId); // 負の値 = 逆鎖でのヒット
            Assert.Equal(60, hit.LastMatchEndOffset);
        }

        [Fact]
        public void FindDominantUnitig_MatchAtVeryEndOfUnitig_ReturnsFullUnitigLength()
        {
            var contigMaker = this.CreateContigMakerWithSingleUnitig(kmerLength: 8);

            var read = UnitigSeq[^20..]; // unitigの末尾20bp

            var hit = contigMaker.FindDominantUnitig(read);

            Assert.Equal(1, hit.UnitigId);
            Assert.Equal(UnitigSeq.Length, hit.LastMatchEndOffset);
            Assert.Equal(0, hit.RemainingLength);
        }
    }
}
