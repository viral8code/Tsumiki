using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状であることの目印がスキャフォールドまで残ることを固定する。
    ///
    /// 環状かどうかは ContigMaker が名前に書き込み、AssemblyScorer と
    /// 閉じ目の検証がその名前を根拠に数える。間のスキャフォールディングで
    /// 名前を付け替えて目印を落とすと、下流は黙って「環状は0本」と答える。
    /// </summary>
    public class ScaffolderCircularTests : IDisposable
    {
        private const int k長 = 8;

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

        // 90bp の環をちょうど1周する3本。隣り合う unitig が k-1 = 7 塩基ずつ
        // 重なり、末尾の 7 塩基が先頭の 7 塩基と一致する。
        private const string ユニティグA = "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGC";
        private const string ユニティグB = "GGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTA";
        private const string ユニティグC = "TCCGCTAGACATAATCGTTCTGCCTATATCATTGG";

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
            l_コンティグ構築.V_結合_コンティグ(l_コンティグパス, p_優勢閾値: 0.8m, p_最小証拠数: 1);

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

            // 環状の目印を落とすと、ここが 0 になって完全長の判定が通らなくなる。
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
