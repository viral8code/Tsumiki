using Tsumiki.Common;
using Tsumiki.Core.Evaluation;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状の閉じ目を元リードで裏付けられるかの判定を固定する。
    ///
    /// グラフ上で閉じたという主張と、閉じ目を実際に読んだという証拠は別物で、
    /// ここが甘いと線状の断片が完全長として通ってしまう。
    /// </summary>
    public class CircularClosureVerifierTests : IDisposable
    {
        private readonly string _一時ディレクトリ;

        public CircularClosureVerifierTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_closure_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters { A_スレッド数 = 2 };
        }

        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
            GC.SuppressFinalize(this);
        }

        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string 塩基 = "ACGT";
            return string.Concat(
                Enumerable.Range(0, p_長さ).Select(_ => 塩基[l_乱数.Next(4)]));
        }

        private string V_書き出し_FASTA(string p_ID, string p_配列)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, "assembly.fasta");
            using var l_書き込み = new FastaWriter(l_パス);
            l_書き込み.V_書き込み(p_ID, p_配列);
            return l_パス;
        }

        private string V_書き出し_FASTQ(IEnumerable<string> p_リード群)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, "reads.fq");
            using var l_書き込み = new StreamWriter(l_パス);
            var l_番号 = 0;
            foreach (var l_リード in p_リード群)
            {
                l_書き込み.WriteLine($"@read{l_番号++}");
                l_書き込み.WriteLine(l_リード);
                l_書き込み.WriteLine("+");
                l_書き込み.WriteLine(new string('I', l_リード.Length));
            }
            return l_パス;
        }

        /// <summary>
        /// 末尾と先頭を跨ぐリード。閉じ目の左右へ 50bp ずつ踏み込む。
        /// </summary>
        private static string Get_閉じ目を跨ぐリード(string p_配列)
        {
            return string.Concat(p_配列.AsSpan(p_配列.Length - 50), p_配列.AsSpan(0, 50));
        }

        [Fact]
        public void Get_検証結果_閉じ目を跨ぐリードが足りていれば支持される()
        {
            var l_配列 = Get_乱数配列(500, 11);
            var l_FASTA = this.V_書き出し_FASTA("scaffold1_circular", l_配列);

            // 必要本数(5本)を超える本数を与える。
            var l_跨ぐリード = Get_閉じ目を跨ぐリード(l_配列);
            var l_FASTQ = this.V_書き出し_FASTQ(Enumerable.Repeat(l_跨ぐリード, 6));

            var l_結果 = CircularClosureVerifier.Get_検証結果(l_FASTA, l_FASTQ, null);

            var l_1件 = Assert.Single(l_結果);
            Assert.Equal(6, l_1件.A_跨いだリード数);
            Assert.True(l_1件.A_支持されたか);
        }

        [Fact]
        public void Get_検証結果_逆相補で読まれたリードも跨いだものとして数える()
        {
            var l_配列 = Get_乱数配列(500, 12);
            var l_FASTA = this.V_書き出し_FASTA("scaffold1_circular", l_配列);

            var l_跨ぐリード = Util.V_逆相補(Get_閉じ目を跨ぐリード(l_配列));
            var l_FASTQ = this.V_書き出し_FASTQ(Enumerable.Repeat(l_跨ぐリード, 6));

            var l_結果 = CircularClosureVerifier.Get_検証結果(l_FASTA, l_FASTQ, null);

            Assert.True(Assert.Single(l_結果).A_支持されたか);
        }

        [Fact]
        public void Get_検証結果_閉じ目を跨がないリードだけなら支持されない()
        {
            var l_配列 = Get_乱数配列(500, 13);
            var l_FASTA = this.V_書き出し_FASTA("scaffold1_circular", l_配列);

            // 内部だけを読んだリード。グラフ上は閉じていても閉じ目の証拠にはならない。
            var l_FASTQ = this.V_書き出し_FASTQ(
                Enumerable.Range(0, 20).Select(i => l_配列.Substring(i * 10, 100)));

            var l_結果 = CircularClosureVerifier.Get_検証結果(l_FASTA, l_FASTQ, null);

            var l_1件 = Assert.Single(l_結果);
            Assert.Equal(0, l_1件.A_跨いだリード数);
            Assert.False(l_1件.A_支持されたか);
        }

        [Fact]
        public void Get_検証結果_環状でない配列は検証の対象にしない()
        {
            var l_配列 = Get_乱数配列(500, 14);
            var l_FASTA = this.V_書き出し_FASTA("scaffold1", l_配列);
            var l_FASTQ = this.V_書き出し_FASTQ([Get_閉じ目を跨ぐリード(l_配列)]);

            Assert.Empty(CircularClosureVerifier.Get_検証結果(l_FASTA, l_FASTQ, null));
        }

        [Fact]
        public void Get_検証結果_同じリードが何度跨いで見えても1本と数える()
        {
            // 閉じ目の窓を2回含むリードを作る(短い環状を2周ぶん読んだ形)。
            var l_配列 = Get_乱数配列(120, 15);
            var l_FASTA = this.V_書き出し_FASTA("plasmid_circular", l_配列);
            var l_2周 = l_配列 + l_配列 + l_配列;
            var l_FASTQ = this.V_書き出し_FASTQ([l_2周]);

            Assert.Equal(1, Assert.Single(CircularClosureVerifier.Get_検証結果(
                l_FASTA, l_FASTQ, null)).A_跨いだリード数);
        }
    }
}
