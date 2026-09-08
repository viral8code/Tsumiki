using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// k ごとのチェックポイントの保存と再開を固定する。
    ///
    /// 再利用してよい条件を緩めると、違う設定で作った成果物が混ざる。
    /// 署名が一致したときだけ返すことを確かめる。
    /// </summary>
    public class CheckpointStoreTests : IDisposable
    {
        private readonly string _一時ディレクトリ;

        private readonly string _リード1のパス;

        private readonly string _リード2のパス;

        public CheckpointStoreTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_ckpt_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);

            this._リード1のパス = Path.Combine(this._一時ディレクトリ, "r1.fq");
            this._リード2のパス = Path.Combine(this._一時ディレクトリ, "r2.fq");
            File.WriteAllText(this._リード1のパス, "@a\nACGT\n+\nIIII\n");
            File.WriteAllText(this._リード2のパス, "@a\nACGT\n+\nIIII\n");
        }

        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
            GC.SuppressFinalize(this);
        }

        private Parameters Get_引数()
        {
            return new Parameters
            {
                A_リード1のパス = this._リード1のパス,
                A_リード2のパス = this._リード2のパス,
            };
        }

        private アセンブリ実行結果 Get_結果()
        {
            var l_ユニティグ = Path.Combine(this._一時ディレクトリ, "unitigs.fasta");
            var l_コンティグ = Path.Combine(this._一時ディレクトリ, "contigs.fasta");
            var l_スキャフォールド = Path.Combine(this._一時ディレクトリ, "scaffolds.fasta");
            File.WriteAllText(l_ユニティグ, ">1\nACGT\n");
            File.WriteAllText(l_コンティグ, ">1\nACGT\n");
            File.WriteAllText(l_スキャフォールド, ">1\nACGT\n");
            return new アセンブリ実行結果(
                31, l_ユニティグ, l_コンティグ, l_スキャフォールド, 3, 42.5, null,
                new 整合性検査結果(1000, 900, 800, 10, 2, 5));
        }

        [Fact]
        public void Get_再開結果_同じ署名なら保存した内容をそのまま返す()
        {
            var l_引数 = this.Get_引数();
            var l_署名 = CheckpointStore.Get_署名(l_引数, 31);
            var l_結果 = this.Get_結果();
            List<引き継ぎ配列> l_引き継ぎ = [new 引き継ぎ配列("ACGTACGT", [3, 4, 5], 31)];

            CheckpointStore.V_保存(this._一時ディレクトリ, l_署名, l_結果, l_引き継ぎ);

            List<引き継ぎ配列> l_読み直し = [];
            var l_再開 = CheckpointStore.Get_再開結果(this._一時ディレクトリ, l_署名, l_読み直し);

            Assert.NotNull(l_再開);
            Assert.Equal(31, l_再開!.A_k長);
            Assert.Equal(3ul, l_再開.A_kmerカットオフ);
            Assert.Equal(42.5, l_再開.A_単一コピー基準値);
            Assert.Equal(l_結果.A_スキャフォールドパス, l_再開.A_スキャフォールドパス);
            Assert.Equal(10, l_再開.A_整合性検査!.Value.A_取りこぼし数);

            var l_1本 = Assert.Single(l_読み直し);
            Assert.Equal("ACGTACGT", l_1本.A_配列);
            Assert.Equal([3, 4, 5], l_1本.A_カバレッジ);
            Assert.Equal(31, l_1本.A_k長);
        }

        [Fact]
        public void Get_再開結果_署名が違えば再利用しない()
        {
            var l_引数 = this.Get_引数();
            CheckpointStore.V_保存(
                this._一時ディレクトリ, CheckpointStore.Get_署名(l_引数, 31), this.Get_結果(), null);

            Assert.Null(CheckpointStore.Get_再開結果(
                this._一時ディレクトリ, CheckpointStore.Get_署名(l_引数, 41), null));
        }

        [Fact]
        public void Get_署名_結果に効くオプションが変われば署名も変わる()
        {
            var l_基準 = CheckpointStore.Get_署名(this.Get_引数(), 31);

            var l_SuperRead付き = this.Get_引数();
            l_SuperRead付き.A_SuperReadを作るか = true;
            Assert.NotEqual(l_基準, CheckpointStore.Get_署名(l_SuperRead付き, 31));

            var l_閾値違い = this.Get_引数();
            l_閾値違い.A_ペア支持数閾値 = 99;
            Assert.NotEqual(l_基準, CheckpointStore.Get_署名(l_閾値違い, 31));
        }

        [Fact]
        public void Get_署名_結果に影響しないオプションでは署名を変えない()
        {
            var l_基準 = CheckpointStore.Get_署名(this.Get_引数(), 31);

            // スレッド数と言語は決定的な結果を変えない。
            var l_別条件 = this.Get_引数();
            l_別条件.A_スレッド数 = 1;
            l_別条件.A_言語 = 言語.英語;

            Assert.Equal(l_基準, CheckpointStore.Get_署名(l_別条件, 31));
        }

        [Fact]
        public void Get_再開結果_成果物が消えていれば再利用しない()
        {
            var l_引数 = this.Get_引数();
            var l_署名 = CheckpointStore.Get_署名(l_引数, 31);
            var l_結果 = this.Get_結果();
            CheckpointStore.V_保存(this._一時ディレクトリ, l_署名, l_結果, null);

            File.Delete(l_結果.A_コンティグパス);

            Assert.Null(CheckpointStore.Get_再開結果(this._一時ディレクトリ, l_署名, null));
        }

        [Fact]
        public void Get_再開結果_チェックポイントが無ければnullを返す()
        {
            Assert.Null(CheckpointStore.Get_再開結果(this._一時ディレクトリ, "signature", null));
        }
    }
}
