using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 複数のアセンブリ候補から 1 つを選ぶ規則の検証
    /// </summary>
    public class AssemblySelectorTests
    {
        #region 公開メソッド

        /// <summary>
        /// 候補が空なら null を返すことを確かめる
        /// </summary>
        [Fact]
        public void V_候補が空ならnullを返す()
        {
            Assert.Null(AssemblySelector.Get_最良([]));
        }

        /// <summary>
        /// 候補が 1 つならそれを返すことを確かめる
        /// </summary>
        [Fact]
        public void V_候補が1つならそれを返す()
        {
            var l_候補 = Get_候補(31, 50_000L, 0.97D);
            var l_選択 = AssemblySelector.Get_最良([l_候補]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性が同程度なら、連続性が高いほうを採る
        /// </summary>
        [Fact]
        public void V_完全性が同程度なら最も連続性が高い候補を選ぶ()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 63_058L, 0.980D), Get_候補(45, 151_085L, 0.981D), Get_候補(63, 175_674L, 0.979D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(63, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// R. sphaeroides の実測値そのもの
        /// </summary>
        [Fact]
        public void V_実測値のばらつきでも既知の最良kを選ぶ()
        {
            var l_選択 = AssemblySelector.Get_最良([
                Get_候補(31, 53_893L, 0.9708D),
                Get_候補(45, 36_387L, 0.9684D),
                Get_候補(55, 19_347L, 0.9528D),
                Get_候補(63, 15_750L, 0.9316D),
            ]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// これが二段階にした理由
        /// </summary>
        [Fact]
        public void V_連続性が高くても完全性が許容差を超えて落ちる候補は選ばれない()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 8_300L, 0.933D), Get_候補(63, 16_300L, 0.675D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 同じ領域を重複して出している候補も、連続性が高くても採らないことを確かめる
        /// </summary>
        [Fact]
        public void V_連続性が高くても重複した候補は選ばれない()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 50_000L, 0.97D, p_正確性: 0.99D), Get_候補(63, 90_000L, 0.97D, p_正確性: 0.60D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性には 2 つの段階がある
        /// </summary>
        [Fact]
        public void V_完全性の差が同点とみなす幅を超えていれば連続性より先に負ける()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 10_000L, 0.99D), Get_候補(63, 90_000L, 0.99D - AssemblySelector.C_同点とみなす差 - 0.001D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 差が同点とみなす幅に収まっていれば、完全性では決めずに連続性で決める
        /// </summary>
        [Fact]
        public void V_完全性の差が同点とみなす幅に収まれば連続性で決める()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 10_000L, 0.99D), Get_候補(63, 90_000L, 0.99D - AssemblySelector.C_同点とみなす差 + 0.001D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(63, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性の差が許容差をわずかに超えると候補が選ばれないことを確かめる
        /// </summary>
        [Fact]
        public void V_完全性の差が許容差をわずかに超えると候補が選ばれない()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 10_000L, 0.99D), Get_候補(63, 90_000L, 0.99D - AssemblySelector.C_完全性の許容差 - 0.001D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 提案 H: 完全性・正確性が同程度でも、より多くの複製単位を環状に閉じられた候補を、NG50 より優先して選ぶこと
        /// </summary>
        [Fact]
        public void V_より多くの複製単位を環状に閉じた候補はNG50より優先される()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(63, 200_000L, 0.97D, p_環状本数: 0, p_環状化率: 0.0D), Get_候補(31, 50_000L, 0.97D, p_環状本数: 2, p_環状化率: 0.98D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(31, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 環状本数が同じなら、閉じた総塩基がゲノム推定サイズに占める割合 (環状化率) で比べる
        /// </summary>
        [Fact]
        public void V_環状本数が同じなら環状化率が高い候補を選ぶ()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 200_000L, 0.97D, p_環状本数: 1, p_環状化率: 0.30D), Get_候補(63, 50_000L, 0.97D, p_環状本数: 1, p_環状化率: 0.95D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(63, l_選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 環状化の状況が全く同じ (両方 0) なら、これまでどおり NG50 で決める (既存の挙動を壊していないことの確認)
        /// </summary>
        [Fact]
        public void V_どの候補も環状でなければNG50で決める()
        {
            var l_選択 = AssemblySelector.Get_最良([Get_候補(31, 50_000L, 0.97D), Get_候補(63, 90_000L, 0.97D), ]);
            Assert.NotNull(l_選択);
            Assert.Equal(63, l_選択.Value.A_実行結果.A_k長);
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
        /// <param name="p_環状本数">候補の環状本数</param>
        /// <param name="p_環状化率">候補の環状化率</param>
        /// <returns>実行結果と評価の組</returns>
        private static (アセンブリ実行結果, アセンブリ評価) Get_候補(int p_k長, long p_NG50, double p_完全性, double p_正確性 = 1.0D, int p_環状本数 = 0, double p_環状化率 = 0D)
        {
            const long l_期待延べ数 = 1_000_000L;
            var l_実行結果 = new アセンブリ実行結果(p_k長, $"k{p_k長}_unitigs.fasta", $"k{p_k長}_contigs.fasta", $"k{p_k長}_scaffolds.fasta", 2UL, 20.0D);
            var l_評価 = new アセンブリ評価(A_期待延べ数: l_期待延べ数, A_欠損延べ数: (long)(l_期待延べ数 * (1D - p_完全性)), A_過剰延べ数: (long)(l_期待延べ数 * (1D - p_正確性)), A_総延長: 5_000_000L, A_本数: 100, A_NG50: p_NG50, A_環状本数: p_環状本数, A_環状化率: p_環状化率);
            return (l_実行結果, l_評価);
        }

        #endregion
    }
}
