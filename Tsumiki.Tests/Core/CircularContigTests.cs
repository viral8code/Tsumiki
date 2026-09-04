using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状に閉じた複製単位(細菌の染色体・プラスミドはいずれも環状)を
    /// 組み上げられた場合に、それを検出して名前で示し、かつ円周の長さが
    /// 正しくなることを検証する。
    ///
    /// 環状経路では末尾 unitig が「始点 unitig と重なる k-1 塩基」を自分の
    /// 末尾に持っている。線状の連結では次の unitig 側から重なりを取り除くが、
    /// 環状では「次」が既に出力済みの始点なので取り除く相手がおらず、
    /// そのままだと円周が k-1 塩基ぶん長く出てしまう。
    /// </summary>
    public class CircularContigTests : IDisposable
    {
        private readonly string _tempDir;
        private readonly string _originalCurrentDirectory;

        public CircularContigTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_circular_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
            this._originalCurrentDirectory = Environment.CurrentDirectory;
        }

        public void Dispose()
        {
            Environment.CurrentDirectory = this._originalCurrentDirectory;
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        // 90bp の環状配列。k=8 で3本の unitig 全体を通じて重複する
        // 正規化 k-mer が無いことを確認済み。
        private const string Circle = "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATA";

        // 隣り合う unitig が k-1 = 7 塩基ずつ重なり、末尾 unitig の末尾 7 塩基が
        // 先頭 unitig の先頭 7 塩基と一致する(= 環が閉じる)ように切り分けたもの。
        private const string UnitigA = "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGC";
        private const string UnitigB = "GGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTA";
        private const string UnitigC = "TCCGCTAGACATAATCGTTCTGCCTATATCATTGG";

        [Fact]
        public void UniteContigs_ClosedCircle_IsMarkedCircular_AndHasExactCircumferenceLength()
        {
            const int k = 8;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var unitigsPath = Path.Combine(this._tempDir, "unitigs.fasta");
            File.WriteAllText(unitigsPath, $">1\n{UnitigA}\n>2\n{UnitigB}\n>3\n{UnitigC}\n");

            var contigPath = Path.Combine(this._tempDir, "contigs.fasta");
            var contigMaker = new ContigMaker(unitigsPath);

            // リードを与えなくても、環状の3本は各頂点の出次数がちょうど1なので
            // 相互一意性を満たし、そのまま1周に結合されるはず。
            contigMaker.V_結合_コンティグ(contigPath, p_優勢閾値: 0.8m, p_最小証拠数: 1);

            List<(string A_ID, string A_配列)> contigs = [];
            using (var reader = new FastaReader(contigPath))
            {
                while (reader.Get_続きがあるか())
                {
                    var seq = reader.Get_次の配列();
                    contigs.Add((seq.A_ID.TrimStart('>'), seq.A_配列));
                }
            }

            var contig = Assert.Single(contigs);
            Assert.Contains("circular", contig.A_ID);

            // 重なりを二重に数えず、円周ちょうどの長さになっていること。
            Assert.Equal(Circle.Length, contig.A_配列.Length);

            // 配列としても、環状配列のいずれかの回転(またはその逆相補)に
            // 一致していなければならない。
            var doubled = Circle + Circle;
            var doubledRevComp = Util.V_逆相補(Circle) + Util.V_逆相補(Circle);
            Assert.True(
                doubled.Contains(contig.A_配列) || doubledRevComp.Contains(contig.A_配列),
                $"assembled circle did not match any rotation of the true circle: {contig.A_配列}");
        }

        [Fact]
        public void UniteContigs_LinearPath_IsNotMarkedCircular()
        {
            const int k = 8;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            // 環を閉じる最後の unitig を外し、A -> B の線状経路だけにする。
            var unitigsPath = Path.Combine(this._tempDir, "unitigs_linear.fasta");
            File.WriteAllText(unitigsPath, $">1\n{UnitigA}\n>2\n{UnitigB}\n");

            var contigPath = Path.Combine(this._tempDir, "contigs_linear.fasta");
            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.V_結合_コンティグ(contigPath, p_優勢閾値: 0.8m, p_最小証拠数: 1);

            List<(string A_ID, string A_配列)> contigs = [];
            using (var reader = new FastaReader(contigPath))
            {
                while (reader.Get_続きがあるか())
                {
                    var seq = reader.Get_次の配列();
                    contigs.Add((seq.A_ID.TrimStart('>'), seq.A_配列));
                }
            }

            var contig = Assert.Single(contigs);
            Assert.DoesNotContain("circular", contig.A_ID);
            // A(38bp) + B の重なりを除いた分(38 - 7 = 31bp)= 69bp。
            Assert.Equal(UnitigA.Length + UnitigB.Length - (k - 1), contig.A_配列.Length);
        }
    }
}
