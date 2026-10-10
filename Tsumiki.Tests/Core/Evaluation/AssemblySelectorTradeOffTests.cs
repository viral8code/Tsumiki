using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 完全性と正確性が k に対して逆向きに動く場面での選択を固定する
    /// </summary>
    public class AssemblySelectorTradeOffTests
    {
        #region 公開メソッド

        /// <summary>
        /// 完全性・正確性の差が揺らぎの範囲に収まるときは連続性で決めることを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_揺らぎの範囲の差では連続性で決める()
        {
            var l_選択 = AssemblySelector.Get_最良(Get_実データの候補());
            Assert.NotNull(l_選択);
            Assert.Equal(93, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 正確性に実質的な差があれば連続性より正確性を優先することを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_正確性に実質的な差があればそちらを優先する()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(21, 38_261L, 0.9843D, 0.9994D), Get_候補(93, 277_063L, 0.9881D, 0.9700D), ]);
            Assert.Equal(21, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性で足切りされた候補は連続性がどれほど高くても選ばれないことを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_完全性で足切りされた候補は連続性に関わらず選ばれない()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(21, 38_261L, 0.9843D, 0.9994D), Get_候補(93, 900_000L, 0.9500D, 0.9994D), ]);
            Assert.Equal(21, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 統合で正確性が大きく落ちた候補は退けることを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_統合で正確性が大きく落ちた候補は退ける()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(21, 38_261L, 0.9843D, 0.9994D), Get_候補(21, 44_916L, 0.9851D, 0.8697D), ]);
            Assert.Equal(38_261L, l_選択!.Value.A_評価.A_NG50);
        }

        /// <summary>
        /// 同点幅に収まる差は同じ段になることを確かめる
        /// </summary>
        [Fact]
        public void Get_段表_同点幅に収まる差は同じ段になる()
        {
            var l_段 = AssemblySelector.Get_段表([0.9994D, 0.9976D, 0.9970D]);
            Assert.Equal(0, l_段[0.9994D]);
            Assert.Equal(0, l_段[0.9976D]);
            Assert.Equal(0, l_段[0.9970D]);
        }

        /// <summary>
        /// 同点幅を超える差があれば段が下がることを確かめる
        /// </summary>
        [Fact]
        public void Get_段表_同点幅を超える差で段が下がる()
        {
            var l_段 = AssemblySelector.Get_段表([0.9994D, 0.9900D]);
            Assert.Equal(0, l_段[0.9994D]);
            Assert.True(l_段[0.9900D] > 0);
        }

        /// <summary>
        /// 隣との差で切るので、最良値からの距離が同点幅を超えていても間が詰まっていれば同じ段になることを確かめる
        /// </summary>
        [Fact]
        public void Get_段表_隣が詰まっていれば最良から離れていても同じ段になる()
        {
            var l_段 = AssemblySelector.Get_段表([0.9650D, 0.9612D, 0.9582D]);
            Assert.Equal(0, l_段[0.9650D]);
            Assert.Equal(0, l_段[0.9612D]);
            Assert.Equal(0, l_段[0.9582D]);
        }

        /// <summary>
        /// 完全性が同点の範囲なら、連続性で選ばれることを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_完全性が同点なら連続性で選ぶ()
        {
            var l_選択 = AssemblySelector.Get_最良([
                Get_候補(53, 134_913L, 0.9582D, 0.9994D),
                Get_候補(87, 55_706L, 0.9612D, 0.9978D),
                Get_候補(139, 55_618L, 0.9650D, 0.9958D),
            ]);
            Assert.Equal(53, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// NG50 の差が小さく、完全性がはっきり高い候補へ乗り換えることを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_NG50の差が小さく完全性がはっきり高ければ乗り換える()
        {
            var l_選択 = AssemblySelector.Get_最良([
                Get_候補(21, 77_083L, 0.9666D, 0.9999D),
                Get_候補(41, 112_758L, 0.9709D, 0.9998D),
                Get_候補(59, 104_681L, 0.9746D, 0.9995D),
                Get_候補(83, 104_705L, 0.9756D, 0.9992D),
                Get_候補(119, 104_741L, 0.9762D, 0.9989D),
            ]);
            Assert.Equal(119, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// NG50 が大きく劣る候補へは、完全性が高くても乗り換えないことを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_NG50が大きく劣れば完全性が高くても乗り換えない()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(63, 150_000L, 0.9700D, 0.9990D), Get_候補(127, 110_000L, 0.9740D, 0.9990D), ]);
            Assert.Equal(63, l_選択!.Value.A_実行結果.A_k長);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 検証用のアセンブリ候補を組み立てる
        /// </summary>
        /// <param name="p_k長">候補の k 長</param>
        /// <param name="p_NG50">候補の NG50</param>
        /// <param name="p_完全性">候補の完全性</param>
        /// <param name="p_正確性">候補の正確性</param>
        /// <returns>実行結果と評価の組</returns>
        private static (アセンブリ実行結果, アセンブリ評価) Get_候補(int p_k長, long p_NG50, double p_完全性, double p_正確性)
        {
            const long l_期待延べ数 = 1_000_000L;
            var l_実行結果 = new アセンブリ実行結果(p_k長, $"k{p_k長}_unitigs.fasta", $"k{p_k長}_contigs.fasta", $"k{p_k長}_scaffolds.fasta", 2UL, 20.0D);
            var l_評価 = new アセンブリ評価(A_期待延べ数: l_期待延べ数, A_欠損延べ数: (long)Math.Round(l_期待延べ数 * (1D - p_完全性)), A_過剰延べ数: (long)Math.Round(l_期待延べ数 * (1D - p_正確性)), A_総延長: 7_200_000L, A_本数: 100, A_NG50: p_NG50);
            return (l_実行結果, l_評価);
        }

        /// <summary>
        /// 完全性と正確性が逆方向へ動く 6 候補
        /// </summary>
        /// <returns></returns>
        private static List<(アセンブリ実行結果, アセンブリ評価)> Get_実データの候補()
        {
            return[
                Get_候補(21, 38_261L, 0.9843D, 0.9994D),
                Get_候補(29, 108_311L, 0.9859D, 0.9990D),
                Get_候補(43, 173_868L, 0.9868D, 0.9987D),
                Get_候補(63, 229_134L, 0.9875D, 0.9983D),
                Get_候補(93, 277_063L, 0.9881D, 0.9976D),
                Get_候補(135, 260_804L, 0.9886D, 0.9970D),
            ];
        }

        #endregion
    }
}
