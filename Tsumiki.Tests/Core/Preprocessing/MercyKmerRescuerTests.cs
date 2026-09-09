using Tsumiki.Common;
using Tsumiki.Core.Preprocessing;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// カットオフで落ちた k-mer の救済を固定する<br/>
    /// 救うべきもの(信頼できる k-mer に挟まれた低頻度)と、救ってはいけないもの
    /// (端に生えているだけの低頻度)の線引きが要点
    /// </summary>
    public class MercyKmerRescuerTests : IDisposable
    {
        private const int k長 = 21;

        private readonly string _一時ディレクトリ;

        public MercyKmerRescuerTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_mercy_" + Guid.NewGuid().ToString("N"));
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

        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string 塩基 = "ACGT";
            return string.Concat(
                Enumerable.Range(0, p_長さ).Select(_ => 塩基[l_乱数.Next(4)]));
        }

        private string V_書き出し_FASTQ(string p_名前, IEnumerable<string> p_リード群)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, p_名前);
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
        /// 指定した窓だけ観測回数を 1 にし、残りを 5 にした k-mer インデックスを作る<br/>
        /// カットオフ 2 で、その窓だけが落ちた状態になる
        /// </summary>
        private TrustedKmerIndex Get_穴のあるインデックス(string p_配列, int p_穴の開始, int p_穴の長さ)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._一時ディレクトリ);
            var l_塩基列 = p_配列.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + k長 <= l_塩基列.Length; i++)
            {
                var l_穴か = i >= p_穴の開始 && i < p_穴の開始 + p_穴の長さ;
                for (var l_回 = 0; l_回 < (l_穴か ? 1 : 5); l_回++)
                {
                    l_インデックス.V_登録(l_塩基列.AsSpan(i, k長), p_ワーカー番号: 0);
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2);
            return l_インデックス;
        }

        private static bool Get_含まれるか(TrustedKmerIndex p_インデックス, string p_配列, int p_位置)
        {
            var l_kmer = p_配列.Substring(p_位置, k長).Select(Util.Get_塩基ID).ToArray();
            return p_インデックス.Get_含まれるか(l_kmer);
        }

        [Fact]
        public void Get_救済数_信頼できるkmerに挟まれた低頻度kmerを救う()
        {
            var l_配列 = Get_乱数配列(300, 21);
            const int 穴の開始 = 100;
            const int 穴の長さ = 3;

            using var l_インデックス = this.Get_穴のあるインデックス(l_配列, 穴の開始, 穴の長さ);
            Assert.False(Get_含まれるか(l_インデックス, l_配列, 穴の開始));

            // 穴を跨いで両側の信頼できる窓まで届くリードを 2 本与える
            var l_リード = l_配列.Substring(穴の開始 - 30, 100);
            var l_FASTQ = this.V_書き出し_FASTQ("reads.fq", [l_リード, l_リード]);
            var l_引数 = new Parameters { A_リード1のパス = l_FASTQ, A_スレッド数 = 2, A_k長 = k長 };
            ConfigurationManager.A_実行時引数 = l_引数;

            var l_救済数 = MercyKmerRescuer.Get_救済数(l_引数, l_インデックス, k長);

            Assert.Equal(穴の長さ, l_救済数);
            for (var i = 0; i < 穴の長さ; i++)
            {
                Assert.True(
                    Get_含まれるか(l_インデックス, l_配列, 穴の開始 + i),
                    $"穴の窓 {i} が救済されていない");
            }
        }

        [Fact]
        public void Get_救済数_観測が1回だけなら救わない()
        {
            var l_配列 = Get_乱数配列(300, 22);
            using var l_インデックス = this.Get_穴のあるインデックス(l_配列, 100, 3);

            var l_FASTQ = this.V_書き出し_FASTQ("reads.fq", [l_配列.Substring(70, 100)]);
            var l_引数 = new Parameters { A_リード1のパス = l_FASTQ, A_スレッド数 = 2, A_k長 = k長 };
            ConfigurationManager.A_実行時引数 = l_引数;

            Assert.Equal(0, MercyKmerRescuer.Get_救済数(l_引数, l_インデックス, k長));
        }

        [Fact]
        public void Get_救済数_片側しか信頼できない端の低頻度kmerは救わない()
        {
            var l_配列 = Get_乱数配列(300, 23);
            using var l_インデックス = this.Get_穴のあるインデックス(l_配列, 100, 3);

            // 穴の左側だけを含み、右側の信頼できる窓まで届かないリード
            var l_リード = l_配列.Substring(80, k長 + 22);
            var l_FASTQ = this.V_書き出し_FASTQ("reads.fq", [l_リード, l_リード]);
            var l_引数 = new Parameters { A_リード1のパス = l_FASTQ, A_スレッド数 = 2, A_k長 = k長 };
            ConfigurationManager.A_実行時引数 = l_引数;

            Assert.Equal(0, MercyKmerRescuer.Get_救済数(l_引数, l_インデックス, k長));
        }

        [Fact]
        public void Get_救済数_kが64を超える場合は何もしない()
        {
            var l_配列 = Get_乱数配列(300, 24);
            using var l_インデックス = this.Get_穴のあるインデックス(l_配列, 100, 3);

            var l_FASTQ = this.V_書き出し_FASTQ("reads.fq", [l_配列]);
            var l_引数 = new Parameters { A_リード1のパス = l_FASTQ, A_スレッド数 = 1, A_k長 = k長 };
            ConfigurationManager.A_実行時引数 = l_引数;

            Assert.Equal(0, MercyKmerRescuer.Get_救済数(l_引数, l_インデックス, p_k長: 65));
        }
    }
}
