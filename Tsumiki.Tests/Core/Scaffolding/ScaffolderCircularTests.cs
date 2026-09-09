using Tsumiki.Common;
using Tsumiki.Core.Evaluation;
using Tsumiki.Core.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状であることの目印がスキャフォールドまで残ることを固定する<br/>
    /// 環状かどうかは ContigMaker が名前に書き込み、AssemblyScorer と
    /// 閉じ目の検証がその名前を根拠に数える<br/>
    /// 間のスキャフォールディングで
    /// 名前を付け替えて目印を落とすと、下流は黙って「環状は 0 本」と答える
    /// </summary>
    public class ScaffolderCircularTests : IDisposable
    {
        /// <summary>
        /// k 長
        /// </summary>
        private const int k長 = 21;

        /// <summary>
        /// 円周
        /// </summary>
        private const int 円周 = 1200;

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        public ScaffolderCircularTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_scaffold_circ_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
            GC.SuppressFinalize(this);
        }

        // 環をちょうど 1 周する 3 本
        // 隣り合う unitig が k-1 塩基ずつ重なり、
        // 末尾の k-1 塩基が先頭の k-1 塩基と一致する
        // 複製単位として数えてもらえる長さ (Consts.環状として数える最小長) を
        // 超えるようにしないと、環状の目印が付かない
        /// <summary>
        /// 環
        /// </summary>
        private static readonly string 環 = Get_乱数配列(円周, p_種: 20250908);

        /// <summary>
        /// ユニティグ A
        /// </summary>
        private static readonly string ユニティグA = 環[..(400 + k長 - 1)];

        /// <summary>
        /// ユニティグ B
        /// </summary>
        private static readonly string ユニティグB = 環[400..(800 + k長 - 1)];

        /// <summary>
        /// ユニティグ C
        /// </summary>
        private static readonly string ユニティグC = 環[800..] + 環[..(k長 - 1)];

        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string 塩基 = "ACGT";
            return string.Concat(
                Enumerable.Range(0, p_長さ).Select(_ => 塩基[l_乱数.Next(4)]));
        }

        private string Get_スキャフォールド出力()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = k長,
                A_スレッド数 = 1,
                A_インサートサイズ = 40,
            };

            var l_ユニティグパス = Path.Combine(this._一時ディレクトリ, "unitigs.fasta");
            File.WriteAllText(
                l_ユニティグパス, $">1\n{ユニティグA}\n>2\n{ユニティグB}\n>3\n{ユニティグC}\n");

            var l_コンティグパス = Path.Combine(this._一時ディレクトリ, "contigs.fasta");
            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);
            l_コンティグ構築.V_結合_コンティグ(l_コンティグパス, p_優勢閾値: 0.8M, p_最小証拠数: 1);

            var l_スキャフォールドパス = Path.Combine(this._一時ディレクトリ, "scaffolds.fasta");
            new Scaffolder(l_コンティグ構築, l_コンティグパス, p_リード長: 30)
                .V_実行(l_スキャフォールドパス);
            return l_スキャフォールドパス;
        }

        [Fact]
        public void V_実行_単独で出た環状コンティグは環状のまま名前に残る()
        {
            var l_スキャフォールドパス = this.Get_スキャフォールド出力();

            Assert.True(File.Exists(l_スキャフォールドパス), "スキャフォールドが出力されていない");

            var l_エントリ = Assert.Single(FastaReader.Get_全エントリ(l_スキャフォールドパス));
            Assert.Contains(Consts.環状の目印, l_エントリ.A_ID, StringComparison.OrdinalIgnoreCase);

            // 環状の目印を落とすと、ここが 0 になって完全長の判定が通らなくなる
            Assert.Equal(1, CompletenessValidator.Get_環状本数(l_スキャフォールドパス));
        }

        [Fact]
        public void V_実行_コンティグ側の環状判定がそのまま引き継がれる()
        {
            var l_スキャフォールドパス = this.Get_スキャフォールド出力();
            var l_コンティグパス = Path.Combine(this._一時ディレクトリ, "contigs.fasta");

            Assert.Equal(
                CompletenessValidator.Get_環状本数(l_コンティグパス),
                CompletenessValidator.Get_環状本数(l_スキャフォールドパス));
        }
    }
}
