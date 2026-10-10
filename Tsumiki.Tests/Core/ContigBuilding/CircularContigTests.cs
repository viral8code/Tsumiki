using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状に閉じた複製単位 (細菌の染色体・プラスミドはいずれも環状) を組み上げられた場合に、それを検出して名前で示し、かつ円周の長さが正しくなることを検証する
    /// </summary>
    public class CircularContigTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki circular tests
        /// </summary>
        private const string C_項目_tsumiki_circular_tests = "tsumiki_circular_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 unitigs fasta
        /// </summary>
        private const string C_ファイル名_unitigs_fasta = "unitigs.fasta";

        /// <summary>
        /// ファイル名 contigs fasta
        /// </summary>
        private const string C_ファイル名_contigs_fasta = "contigs.fasta";

        /// <summary>
        /// ファイル名 unitigs linear fasta
        /// </summary>
        private const string C_ファイル名_unitigs_linear_fasta = "unitigs_linear.fasta";

        /// <summary>
        /// ファイル名 contigs linear fasta
        /// </summary>
        private const string C_ファイル名_contigs_linear_fasta = "contigs_linear.fasta";

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k = 21;

        /// <summary>
        /// 円周
        /// </summary>
        private const int C_円周 = 1_200;

        /// <summary>
        /// 環状の複製単位そのもの
        /// </summary>
        private static readonly string C_Circle = Get_乱数配列(C_円周, p_種: 20_250_908);

        /// <summary>
        /// 環を 3 分割したうちの 1 本目
        /// </summary>
        private static readonly string C_入口unitig = C_Circle[..(400 + C_k - 1)];

        /// <summary>
        /// 環を 3 分割したうちの 2 本目
        /// </summary>
        private static readonly string C_代替入口unitig = C_Circle[400..(800 + C_k - 1)];

        /// <summary>
        /// 環を 3 分割したうちの 3 本目、先頭へ戻る重なりを含む
        /// </summary>
        private static readonly string C_出口unitig = C_Circle[800..] + C_Circle[..(C_k - 1)];

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        /// <summary>
        /// 実行前のカレントディレクトリ
        /// </summary>
        private readonly string _元のカレントディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public CircularContigTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_circular_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            this._元のカレントディレクトリ = Environment.CurrentDirectory;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            Environment.CurrentDirectory = this._元のカレントディレクトリ;
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 閉じた環状経路を結合すると circular の印が付き円周ちょうどの長さになる
        /// </summary>
        [Fact]
        public void V_閉じた環は環状と判定され円周ちょうどの長さになる()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k,
                A_スレッド数 = 1
            };
            var l_unitigパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_unitigs_fasta);
            File.WriteAllText(l_unitigパス, $">1\n{C_入口unitig}\n>2\n{C_代替入口unitig}\n>3\n{C_出口unitig}\n");
            var l_contigパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_contigs_fasta);
            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_結合_Contig(l_contigパス, p_優勢閾値: 0.8M, p_最小証拠数: 1UL);
            List<(string A_ID, string A_配列)> l_contig群 = [];
            using (var l_読み込み = new FastaReader(l_contigパス))
            {
                while (l_読み込み.Has続き())
                {
                    var l_配列 = l_読み込み.Get_次の配列();
                    l_contig群.Add((l_配列.A_ID.TrimStart('>'), l_配列.A_配列));
                }
            }

            var l_contig = Assert.Single(l_contig群);
            Assert.Contains(Consts.環状の目印, l_contig.A_ID);
            Assert.Equal(C_円周, l_contig.A_配列.Length);
            var l_二重配列 = C_Circle + C_Circle;
            var l_二重逆相補配列 = Util.V_逆相補(C_Circle) + Util.V_逆相補(C_Circle);
            Assert.True(l_二重配列.Contains(l_contig.A_配列) || l_二重逆相補配列.Contains(l_contig.A_配列), $"assembled circle did not match any rotation of the true circle: {l_contig.A_配列}");
        }

        /// <summary>
        /// 環を閉じない線状経路には circular の印が付かない
        /// </summary>
        [Fact]
        public void V_環を閉じない経路は環状と判定されない()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k,
                A_スレッド数 = 1
            };
            var l_unitigパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_unitigs_linear_fasta);
            File.WriteAllText(l_unitigパス, $">1\n{C_入口unitig}\n>2\n{C_代替入口unitig}\n");
            var l_contigパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_contigs_linear_fasta);
            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_結合_Contig(l_contigパス, p_優勢閾値: 0.8M, p_最小証拠数: 1UL);
            List<(string A_ID, string A_配列)> l_contig群 = [];
            using (var l_読み込み = new FastaReader(l_contigパス))
            {
                while (l_読み込み.Has続き())
                {
                    var l_配列 = l_読み込み.Get_次の配列();
                    l_contig群.Add((l_配列.A_ID.TrimStart('>'), l_配列.A_配列));
                }
            }

            var l_contig = Assert.Single(l_contig群);
            Assert.DoesNotContain(Consts.環状の目印, l_contig.A_ID);
            Assert.Equal(C_入口unitig.Length + C_代替入口unitig.Length - (C_k - 1), l_contig.A_配列.Length);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string l_塩基 = Consts.塩基文字;
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => l_塩基[l_乱数.Next(4)]));
        }

        #endregion
    }
}
