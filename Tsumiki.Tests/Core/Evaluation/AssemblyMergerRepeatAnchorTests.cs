using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 橋渡しの両側の錨の塊に、推定コピー数が 1 以下の位置があるかの判定の検証
    /// </summary>
    public class AssemblyMergerRepeatAnchorTests
    {
        #region 公開メソッド

        /// <summary>
        /// 塊の中に推定コピー数 1 の位置があれば true になる
        /// </summary>
        [Fact]
        public void V_塊に1コピーの位置があれば一意な錨とみなす()
        {
            var l_塊 = new List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)>
            {
                (0, 0, 100, true),
                (1, 0, 101, true),
                (2, 0, 102, true),
            };
            Func<string, int, int> l_期待コピー数 = (p_配列, p_位置) => p_位置 switch
            {
                0 => 2,
                1 => 1,
                _ => 3,
            };
            Assert.True(AssemblyMerger.Has一意な錨(l_塊, Consts.塩基文字, l_期待コピー数));
        }

        /// <summary>
        /// 塊のすべての位置が 2 コピー以上なら、一意な錨はないので false になる
        /// </summary>
        [Fact]
        public void V_塊のすべてが複数コピーなら一意な錨ではない()
        {
            var l_塊 = new List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)>
            {
                (0, 0, 100, true),
                (1, 0, 101, true),
                (2, 0, 102, true),
            };
            Func<string, int, int> l_期待コピー数 = (p_配列, p_位置) => 2;
            Assert.False(AssemblyMerger.Has一意な錨(l_塊, Consts.塩基文字, l_期待コピー数));
        }

        /// <summary>
        /// 塊が空なら、一意な錨はないので false になる
        /// </summary>
        [Fact]
        public void V_空の塊は一意な錨ではない()
        {
            var l_塊 = new List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)>();
            Func<string, int, int> l_期待コピー数 = (p_配列, p_位置) => 1;
            Assert.False(AssemblyMerger.Has一意な錨(l_塊, Consts.塩基文字, l_期待コピー数));
        }

        /// <summary>
        /// 推定コピー数 0 も「1 以下」に含むので、一意な錨として true になる
        /// </summary>
        [Fact]
        public void V_推定コピー数0の位置も一意な錨に含む()
        {
            var l_塊 = new List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)>
            {
                (0, 0, 100, true),
                (1, 0, 101, true),
            };
            Func<string, int, int> l_期待コピー数 = (p_配列, p_位置) => p_位置 == 0 ? 0 : 2;
            Assert.True(AssemblyMerger.Has一意な錨(l_塊, Consts.塩基文字, l_期待コピー数));
        }

        /// <summary>
        /// 独立支持 2 本の通常の候補だけなら、相互一意として 2 → 4 が確定する
        /// </summary>
        [Fact]
        public void V_通常の候補だけなら相互一意に確定する()
        {
            var l_候補 = new List<橋渡し候補>
            {
                new(2, 4, Consts.塩基文字, 21),
                new(2, 4, Consts.塩基文字, 29),
            };
            var l_確定 = AssemblyMerger.Get_相互一意な橋渡し(l_候補, new List<橋渡し候補>(), 4, 2, new HashSet<(int, int)>(), null, out _, out _);
            Assert.True(l_確定.TryGetValue(2, out var l_橋渡し));
            Assert.Equal(4, l_橋渡し.A_終点);
        }

        /// <summary>
        /// 通常の候補に、同じ始点から別の行き先へ向かう競合だけの候補が 1 本あると、始点 2 は行き先が 2 つになり確定しない
        /// </summary>
        [Fact]
        public void V_競合だけの候補があると行き先が2つになり確定しない()
        {
            var l_候補 = new List<橋渡し候補>
            {
                new(2, 4, Consts.塩基文字, 21),
                new(2, 4, Consts.塩基文字, 29),
            };
            var l_競合だけの候補 = new List<橋渡し候補>
            {
                new(2, 6, Consts.塩基文字, 21, 0, true),
            };
            var l_確定 = AssemblyMerger.Get_相互一意な橋渡し(l_候補, l_競合だけの候補, 4, 2, new HashSet<(int, int)>(), null, out _, out _);
            Assert.False(l_確定.ContainsKey(2));
        }

        /// <summary>
        /// 競合だけの候補しかなければ、支持が無いので何も確定しない
        /// </summary>
        [Fact]
        public void V_競合だけの候補だけでは何も確定しない()
        {
            var l_競合だけの候補 = new List<橋渡し候補>
            {
                new(2, 6, Consts.塩基文字, 21, 0, true),
            };
            var l_確定 = AssemblyMerger.Get_相互一意な橋渡し(new List<橋渡し候補>(), l_競合だけの候補, 4, 2, new HashSet<(int, int)>(), null, out _, out _);
            Assert.Empty(l_確定);
        }

        #endregion
    }
}
