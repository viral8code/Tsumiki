using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 完全性と正確性が k に対して逆向きに動く場面での選択を固定する
    /// </summary>
    /// <remarks>
    /// k を上げるとグラフが解けて完全性が上がり、反復を正しく複製したぶん「出しすぎ」が増えて正確性が下がる<br/>
    /// この 2 つを生の値で順に見ると、先に見たほうの端の k が機械的に選ばれるだけで比較にならない<br/>
    /// 揺らぎの範囲を同点として扱い、実質的な差があるときだけ効くようにする
    /// </remarks>
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

            // 完全性・正確性の差はいずれも同点幅に収まるので NG50 が決める
            // 正確性を生の値で先に見ていた頃は、3 指標のうち 2 つで最下位の
            // k=21 が 0.18 ポイントの差だけで選ばれていた
            Assert.NotNull(l_選択);
            Assert.Equal(93, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 正確性に実質的な差があれば連続性より正確性を優先することを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_正確性に実質的な差があればそちらを優先する()
        {
            // 連続性で勝る候補が、正確性で同点幅をはっきり超えて劣る場合
            var l_選択 = AssemblySelector.Get_最良([ Get_候補(21, 38_261L, 0.9843D, 0.9994D), Get_候補(93, 277_063L, 0.9881D, 0.9700D), ]);

            Assert.Equal(21, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性で足切りされた候補は連続性がどれほど高くても選ばれないことを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_完全性で足切りされた候補は連続性に関わらず選ばれない()
        {
            // 反復を飛ばして繋いだ結果、連続性は跳ね上がるが配列を大きく落とした候補
            var l_選択 = AssemblySelector.Get_最良([ Get_候補(21, 38_261L, 0.9843D, 0.9994D), Get_候補(93, 900_000L, 0.9500D, 0.9994D), ]);

            Assert.Equal(21, l_選択!.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 統合で正確性が大きく落ちた候補は退けることを確かめる
        /// </summary>
        [Fact]
        public void Get_最良_統合で正確性が大きく落ちた候補は退ける()
        {
            // 実データで観測された統合の失敗 (99.94% -> 86.97%)
            var l_選択 = AssemblySelector.Get_最良([ Get_候補(21, 38_261L, 0.9843D, 0.9994D), Get_候補(21, 44_916L, 0.9851D, 0.8697D), ]);

            Assert.Equal(38_261L, l_選択!.Value.A_評価.A_NG50);
        }

        /// <summary>
        /// 同点幅に収まる差は同じ段になることを確かめる
        /// </summary>
        [Fact]
        public void Get_段_同点幅に収まる差は同じ段になる()
        {
            Assert.Equal(0, AssemblySelector.Get_段(0.9994D, 0.9994D));
            Assert.Equal(0, AssemblySelector.Get_段(0.9976D, 0.9994D));
            Assert.Equal(0, AssemblySelector.Get_段(0.9970D, 0.9994D));

            // 同点幅 (0.5 ポイント) を超えれば段が下がる
            Assert.True(AssemblySelector.Get_段(0.9900D, 0.9994D) > 0);
        }

        /// <summary>
        /// 基準より良い値は 0 段のままにすることを確かめる
        /// </summary>
        [Fact]
        public void Get_段_基準より良い値は0段のままにする()
        {
            // 丸めの都合で基準をわずかに上回っても、負の段にはしない
            Assert.Equal(0, AssemblySelector.Get_段(0.9999D, 0.9994D));
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
        /// 7.4 Mbp の実データで観測された 6 候補
        /// </summary>
        /// <remarks>
        /// 完全性は単調増加、正確性は単調減少し、NG50 は k=93 で最大になる
        /// </remarks>
        /// <returns></returns>
        private static List<(アセンブリ実行結果, アセンブリ評価)> Get_実データの候補()
        {
            return
            [
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
