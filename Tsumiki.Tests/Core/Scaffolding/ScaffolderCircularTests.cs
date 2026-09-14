using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状であることの目印が scaffold まで残ることを固定する
    /// </summary>
    /// <remarks>
    /// 環状かどうかは ContigMaker が名前に書き込み、AssemblyScorer と閉じ目の検証がその名前を根拠に数える<br/>
    /// 間のスキャフォールディングで名前を付け替えて目印を落とすと、下流は黙って「環状は 0 本」と答える
    /// </remarks>
    public class ScaffolderCircularTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// k 長
        /// </summary>
        private const int k長 = 21;

        /// <summary>
        /// 円周
        /// </summary>
        private const int 円周 = 1_200;

        // 環をちょうど 1 周する 3 本
        // 隣り合う unitig が k-1 塩基ずつ重なり、
        // 末尾の k-1 塩基が先頭の k-1 塩基と一致する
        // 複製単位として数えてもらえる長さ (Consts.環状として数える最小長) を
        // 超えるようにしないと、環状の目印が付かない

        /// <summary>
        /// 環
        /// </summary>
        private static readonly string 環 = Get_乱数配列(円周, p_種: 20_250_908);

        /// <summary>
        /// unitig A
        /// </summary>
        private static readonly string unitigA = 環[..(400 + k長 - 1)];

        /// <summary>
        /// unitig B
        /// </summary>
        private static readonly string unitigB = 環[400..(800 + k長 - 1)];

        /// <summary>
        /// unitig C
        /// </summary>
        private static readonly string unitigC = 環[800..] + 環[..(k長 - 1)];

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public ScaffolderCircularTests()
        {
            this._一時ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_scaffold_circ_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

        #endregion

        #region 公開メソッド

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
        /// 単独で出た環状 contig が、スキャフォールディング後も環状の目印を名前に保つことを検証する
        /// </summary>
        [Fact]
        public void V_実行_単独で出た環状contigは環状のまま名前に残る()
        {
            var l_scaffoldパス = this.Get_scaffold出力();

            Assert.True(File.Exists(l_scaffoldパス), "scaffoldが出力されていない");

            var l_エントリ = Assert.Single(FastaReader.Get_全エントリ(l_scaffoldパス));
            Assert.Contains(Consts.環状の目印, l_エントリ.A_ID, StringComparison.OrdinalIgnoreCase);

            // 環状の目印を落とすと、ここが 0 になって完全長の判定が通らなくなる
            Assert.Equal(1, CompletenessValidator.Get_環状本数(l_scaffoldパス));
        }

        /// <summary>
        /// contig 側の環状判定が、スキャフォールディング後もそのまま引き継がれることを検証する
        /// </summary>
        [Fact]
        public void V_実行_contig側の環状判定がそのまま引き継がれる()
        {
            var l_scaffoldパス = this.Get_scaffold出力();
            var l_contigパス = Path.Combine(this._一時ディレクトリ, "contigs.fasta");

            Assert.Equal(CompletenessValidator.Get_環状本数(l_contigパス), CompletenessValidator.Get_環状本数(l_scaffoldパス));
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
            const string l_塩基 = "ACGT";
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => l_塩基[l_乱数.Next(4)]));
        }

        /// <summary>
        /// スキャフォールディングまで通して、その出力を返す
        /// </summary>
        /// <returns>scaffold の配列</returns>
        private string Get_scaffold出力()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = k長,
                A_スレッド数 = 1,
                A_インサートサイズ = 40,
            };

            var l_unitigパス = Path.Combine(this._一時ディレクトリ, "unitigs.fasta");
            File.WriteAllText(l_unitigパス, $">1\n{unitigA}\n>2\n{unitigB}\n>3\n{unitigC}\n");

            var l_contigパス = Path.Combine(this._一時ディレクトリ, "contigs.fasta");
            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_結合_Contig(l_contigパス, p_優勢閾値: 0.8M, p_最小証拠数: 1UL);

            var l_scaffoldパス = Path.Combine(this._一時ディレクトリ, "scaffolds.fasta");
            new Scaffolder(l_contig構築, l_contigパス, p_リード長: 30)
                .V_実行(l_scaffoldパス);
            return l_scaffoldパス;
        }

        #endregion

    }
}
