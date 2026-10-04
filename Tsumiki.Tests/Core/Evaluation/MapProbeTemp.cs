using System.Diagnostics;
using System.Security.Cryptography;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Polishing;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core.Evaluation
{
    /// <summary>
    /// 一時的な確認用 (既存のアセンブリにポリッシュと継ぎ目の評価をかけ、時間と出力のハッシュを残す)
    /// </summary>
    public class MapProbeTemp
    {
        /// <summary>
        /// 環境変数 MP_ASM・MP_R1・MP_R2・MP_DIR があるときだけ流す
        /// </summary>
        [Fact]
        public void 実データで流す()
        {
            var l_アセンブリ = Environment.GetEnvironmentVariable("MP_ASM");
            if (string.IsNullOrEmpty(l_アセンブリ))
            {
                return;
            }

            var l_出力先 = Environment.GetEnvironmentVariable("MP_DIR")!;
            ConfigurationManager.A_実行時引数 = new Parameters { A_スレッド数 = 16 };
            (string, string)[] l_ライブラリ群 = [(Environment.GetEnvironmentVariable("MP_R1")!, Environment.GetEnvironmentVariable("MP_R2")!)];

            var l_時計 = Stopwatch.StartNew();
            _ = Polisher.Get_磨いた結果(l_アセンブリ, l_ライブラリ群, Path.Combine(l_出力先, "polished.fasta"));
            var l_ポリッシュ秒 = l_時計.Elapsed.TotalSeconds;

            l_時計.Restart();
            _ = JunctionRiskEvaluator.Get_評価結果(l_アセンブリ, l_ライブラリ群, Path.Combine(l_出力先, "junctions.tsv"));
            var l_継ぎ目秒 = l_時計.Elapsed.TotalSeconds;

            File.WriteAllLines(Path.Combine(l_出力先, "probe.txt"), [
                $"polish_s={l_ポリッシュ秒:F1} sha={Get_ハッシュ(Path.Combine(l_出力先, "polished.fasta"))}",
                $"junction_s={l_継ぎ目秒:F1} sha={Get_ハッシュ(Path.Combine(l_出力先, "junctions.tsv"))}",
            ]);
        }

        /// <summary>
        /// ファイルの SHA-256 の先頭 16 文字
        /// </summary>
        /// <param name="p_パス"></param>
        /// <returns></returns>
        private static string Get_ハッシュ(string p_パス)
        {
            return Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(p_パス)))[..16];
        }
    }
}
