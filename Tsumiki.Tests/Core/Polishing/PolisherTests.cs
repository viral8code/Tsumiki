using Tsumiki.Commons;
using Tsumiki.Cores.Polishing;
using Tsumiki.Cores.Mapping;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列へのリード再マッピングと多数決による置換訂正を固定する
    /// </summary>
    /// <remarks>
    /// 直すべきものを直すことと同じくらい、根拠が無い位置を動かさないことが
    /// 重要なので、両方を確かめる
    /// </remarks>
    public class PolisherTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        public PolisherTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_polish_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters { A_スレッド数 = 2 };
        }

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
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string 塩基 = "ACGT";
            return string.Concat(
                Enumerable.Range(0, p_長さ).Select(_ => 塩基[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ID">配列 ID</param>
        /// <param name="p_配列">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_FASTA(string p_ID, string p_配列)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, "assembly.fasta");
            using var l_書き込み = new FastaWriter(l_パス);
            l_書き込み.V_書き込み(p_ID, p_配列);
            return l_パス;
        }

        /// <summary>
        /// 真の配列から等間隔にリードを切り出して FASTQ にする
        /// </summary>
        /// <param name="p_真の配列">元になる配列</param>
        /// <param name="p_リード長">切り出すリードの長さ</param>
        /// <param name="p_刻み">開始位置の刻み幅</param>
        /// <param name="p_開始">切り出しを始める位置</param>
        /// <param name="p_終了">切り出しを終える位置</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_FASTQ(
            string p_真の配列, int p_リード長, int p_刻み, int p_開始 = 0, int? p_終了 = null)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, $"reads{p_開始}.fq");
            var l_終了 = p_終了 ?? p_真の配列.Length;
            var l_品質 = new string('I', p_リード長);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_番号 = 0;
            for (var i = p_開始; i + p_リード長 <= l_終了; i += p_刻み)
            {
                var l_配列 = p_真の配列.Substring(i, p_リード長);

                // 逆鎖側の経路も通るよう、1 本おきに逆相補で出す
                if (l_番号 % 2 == 1)
                {
                    l_配列 = Util.V_逆相補(l_配列);
                }
                l_書き込み.WriteLine($"@read{l_番号++}");
                l_書き込み.WriteLine(l_配列);
                l_書き込み.WriteLine("+");
                l_書き込み.WriteLine(l_品質);
            }
            return l_パス;
        }

        /// <summary>
        /// リードが支持する塩基へ置換を直すことを確かめる
        /// </summary>
        [Fact]
        public void Get_磨いた結果_リードが支持する塩基へ置換を直す()
        {
            var l_真の配列 = Get_乱数配列(2000, 1);
            var l_誤り位置 = 1000;
            var l_文字 = l_真の配列.ToCharArray();
            l_文字[l_誤り位置] = l_文字[l_誤り位置] == 'A' ? 'C' : 'A';
            var l_誤りを含む配列 = new string(l_文字);

            var l_FASTA = this.V_書き出し_FASTA("contig1", l_誤りを含む配列);
            var l_FASTQ = this.V_書き出し_FASTQ(l_真の配列, p_リード長: 150, p_刻み: 10);
            var l_出力 = Path.Combine(this._一時ディレクトリ, "polished.fasta");
            var l_配置 = new ReadMapper([l_誤りを含む配列]).Get_配置(l_真の配列.Substring(950, 150));
            Console.WriteLine(l_配置.A_整列位置群.First(x => x.A_参照位置 == l_誤り位置));

            var l_統計 = Polisher.Get_磨いた結果(l_FASTA, l_FASTQ, null, l_出力);

            Console.WriteLine(l_統計);

            Assert.NotNull(l_統計);
            Console.WriteLine(l_真の配列.Substring(l_誤り位置 - 3, 7));
            Console.WriteLine(FastaReader.Get_全エントリ(l_出力)[0].A_配列.Substring(l_誤り位置 - 3, 7));
            Assert.Equal(1, l_統計!.Value.A_訂正した塩基数);

            var l_結果 = FastaReader.Get_全エントリ(l_出力);
            _ = Assert.Single(l_結果);
            Assert.Equal(l_真の配列, l_結果[0].A_配列);
        }

        /// <summary>
        /// リードが届かない位置は動かさないことを確かめる
        /// </summary>
        [Fact]
        public void Get_磨いた結果_リードが届かない位置は動かさない()
        {
            var l_真の配列 = Get_乱数配列(2000, 2);
            // 後半に誤りを置き、前半にしかリードを与えない
            var l_誤り位置 = 1600;
            var l_文字 = l_真の配列.ToCharArray();
            l_文字[l_誤り位置] = l_文字[l_誤り位置] == 'G' ? 'T' : 'G';
            var l_誤りを含む配列 = new string(l_文字);

            var l_FASTA = this.V_書き出し_FASTA("contig1", l_誤りを含む配列);
            var l_FASTQ = this.V_書き出し_FASTQ(
                l_真の配列, p_リード長: 150, p_刻み: 10, p_開始: 0, p_終了: 800);
            var l_出力 = Path.Combine(this._一時ディレクトリ, "polished.fasta");

            var l_統計 = Polisher.Get_磨いた結果(l_FASTA, l_FASTQ, null, l_出力);

            Assert.NotNull(l_統計);
            Assert.Equal(0, l_統計!.Value.A_訂正した塩基数);

            var l_結果 = FastaReader.Get_全エントリ(l_出力);
            Assert.Equal(l_誤りを含む配列, l_結果[0].A_配列);
        }

        /// <summary>
        /// リードが届かない範囲は深度不足として数えることを確かめる
        /// </summary>
        [Fact]
        public void Get_磨いた結果_リードが届かない範囲は深度不足として数える()
        {
            var l_真の配列 = Get_乱数配列(2000, 3);
            var l_FASTA = this.V_書き出し_FASTA("contig1", l_真の配列);
            var l_FASTQ = this.V_書き出し_FASTQ(
                l_真の配列, p_リード長: 150, p_刻み: 10, p_開始: 0, p_終了: 1000);
            var l_出力 = Path.Combine(this._一時ディレクトリ, "polished.fasta");

            var l_統計 = Polisher.Get_磨いた結果(l_FASTA, l_FASTQ, null, l_出力);

            Assert.NotNull(l_統計);
            Assert.Equal(2000, l_統計!.Value.A_評価できた位置数);

            // 覆われていない後半のぶんが深度不足として出る
            Assert.True(l_統計.Value.A_深度不足率 > 0.4, $"深度不足率={l_統計.Value.A_深度不足率}");
        }
    }
}
