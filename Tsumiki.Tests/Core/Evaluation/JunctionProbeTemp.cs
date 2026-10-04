using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core.Evaluation
{
    /// <summary>
    /// 一時的な確認用 (既存のアセンブリに C# の継ぎ目の評価だけをかける)
    /// </summary>
    public class JunctionProbeTemp
    {
        /// <summary>
        /// 環境変数 JP_ASM・JP_R1・JP_R2・JP_OUT があるときだけ流す
        /// </summary>
        [Fact]
        public void 実データで流す()
        {
            var l_アセンブリ = Environment.GetEnvironmentVariable("JP_ASM");
            if (string.IsNullOrEmpty(l_アセンブリ))
            {
                return;
            }

            ConfigurationManager.A_実行時引数 = new Parameters { A_スレッド数 = 16 };
            _ = JunctionRiskEvaluator.Get_評価結果(l_アセンブリ, [(Environment.GetEnvironmentVariable("JP_R1")!, Environment.GetEnvironmentVariable("JP_R2")!)], Environment.GetEnvironmentVariable("JP_OUT")!);
        }
    }
}
