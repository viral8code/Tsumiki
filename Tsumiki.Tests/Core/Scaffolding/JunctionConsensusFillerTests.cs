using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列の継ぎ目を、跨ぐリードの品質を見た多数決で埋める処理 (JunctionConsensusFiller) の検証
    /// </summary>
    public class JunctionConsensusFillerTests
    {
        #region 定数

        /// <summary>
        /// 項目 S
        /// </summary>
        private const string C_項目_S = "S";

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 1000;

        /// <summary>
        /// リードの長さ
        /// </summary>
        private const int C_リード長 = 150;

        /// <summary>
        /// 左の片の終わり (G[..500))
        /// </summary>
        private const int C_左の終わり = 500;

        /// <summary>
        /// 右の片の始まり (G[530..))
        /// </summary>
        private const int C_右の始まり = 530;

        /// <summary>
        /// 左右が重なる継ぎ目の、左の片の終わり (G[..520))
        /// </summary>
        private const int C_重なる左の終わり = 520;

        /// <summary>
        /// 左右が重なる継ぎ目の、右の片の始まり (G[480..))
        /// </summary>
        private const int C_重なる右の始まり = 480;

        /// <summary>
        /// 食い違うリードの続きを作る乱数の種
        /// </summary>
        private const int C_食い違う続きの種 = 7;

        /// <summary>
        /// Phred+33 のオフセット
        /// </summary>
        private const int C_オフセット = 33;

        /// <summary>
        /// 品質 40 の文字
        /// </summary>
        private const char C_高い品質 = 'I';

        /// <summary>
        /// 品質 2 の文字
        /// </summary>
        private const char C_低い品質 = '#';

        /// <summary>
        /// 品質 15 の文字
        /// </summary>
        private const char C_品質15 = '0';

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 未確認の繋ぎ目の印を跨ぐリードが 5 本 (うち 2 本は逆向き) あれば、G と同じ 1 本に埋まり、埋めた数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_印を跨ぐリードが5本_Gに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            var (l_結果, l_統計) = Get_実行(l_配列, Get_跨ぐリード群(l_ゲノム, 5, 2));
            Assert.Single(l_結果);
            Assert.Equal(C_項目_S, l_結果[0].A_ID);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_対象数);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 左右の片が 40 塩基重なる (左 = G[..520)、右 = G[480..)) ときも、G と同じ 1 本に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_左右の片が重なる_Gに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..520], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[480..]);
            var (l_結果, l_統計) = Get_実行(l_配列, Get_跨ぐリード群(l_ゲノム, 5, 2));
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 間に読み違いを持つリードが 6 本あっても、その塩基の品質が低ければ数えず、正しいリード 3 本の合意で G に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_品質の低い読み違い_数えずGに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            var l_リード群 = Get_跨ぐリード群(l_ゲノム, 3, 0);
            for (var i = 0; i < 6; i++)
            {
                var l_始め = 410 + (5 * i);
                var l_塩基 = l_ゲノム.Substring(l_始め, C_リード長).ToCharArray();
                var l_品質 = new string(C_高い品質, C_リード長).ToCharArray();
                foreach (var l_位置 in new[] { 505, 512, 518, 524, })
                {
                    l_塩基[l_位置 - l_始め] = Get_別の塩基(l_塩基[l_位置 - l_始め]);
                    l_品質[l_位置 - l_始め] = C_低い品質;
                }

                l_リード群.Add((new string(l_塩基), new string(l_品質)));
            }

            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 跨ぐリードが 2 本だけなら埋めず、配列はそのまま、票が足りない数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_跨ぐリードが2本_埋めない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            var (l_結果, l_統計) = Get_実行(l_配列, Get_跨ぐリード群(l_ゲノム, 2, 0));
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_票が足りない数);
        }

        /// <summary>
        /// 間の長さが違う跨ぐリードが 3 本ずつあれば埋めず、票が割れた数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_間の長さが割れる_埋めない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            var l_挿入のあるゲノム = string.Concat(l_ゲノム[..515], "ACG", l_ゲノム[515..]);
            var l_リード群 = Get_跨ぐリード群(l_ゲノム, 3, 0);
            l_リード群.AddRange(Get_跨ぐリード群(l_挿入のあるゲノム, 3, 0));
            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_票が割れた数);
        }

        /// <summary>
        /// 間が 7 塩基の単位 10 回の縦列反復なら、跨ぐリードがあっても埋めず、縦列反復の数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_縦列反復の中_埋めない()
        {
            var l_元 = Get_ゲノム();
            var l_ゲノム = string.Concat(l_元[..500], string.Concat(Enumerable.Repeat("ACCGTTA", 10)), l_元[500..]);
            var l_配列 = string.Concat(l_ゲノム[..510], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[560..]);
            List<(string A_配列, string A_品質)> l_リード群 = [];
            for (var i = 0; i < 5; i++)
            {
                l_リード群.Add((l_ゲノム.Substring(440 + (5 * i), C_リード長), new string(C_高い品質, C_リード長)));
            }

            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_縦列反復数);
        }

        /// <summary>
        /// N の連続 (50 個) の継ぎ目も、跨ぐリードの合意で G に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_Nの連続_Gに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], new string(Consts.塩基コード.N, 50), l_ゲノム[C_右の始まり..]);
            var (l_結果, l_統計) = Get_実行(l_配列, Get_跨ぐリード群(l_ゲノム, 5, 2));
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 左の片と品質 40 の塩基が 3 か所違うリードは片に置かず、跨ぐリード無しの数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_確かな不一致が3か所_置かない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            List<(string A_配列, string A_品質)> l_リード群 = [];
            for (var i = 0; i < 5; i++)
            {
                var l_始め = 410 + (5 * i);
                var l_塩基 = l_ゲノム.Substring(l_始め, C_リード長).ToCharArray();
                foreach (var l_位置 in new[] { 460, 470, 480, })
                {
                    l_塩基[l_位置 - l_始め] = Get_別の塩基(l_塩基[l_位置 - l_始め]);
                }

                l_リード群.Add((new string(l_塩基), new string(C_高い品質, C_リード長)));
            }

            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_跨ぐリード無し数);
        }

        /// <summary>
        /// 1 本の配列に継ぎ目が 2 つあり、どちらにも跨ぐリードがあれば、両方とも埋まって G に戻る
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_継ぎ目が2つ_両方埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..300], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[310..C_左の終わり], new string(Consts.塩基コード.N, 100), l_ゲノム[C_右の始まり..]);
            List<(string A_配列, string A_品質)> l_リード群 = [];
            for (var i = 0; i < 5; i++)
            {
                l_リード群.Add((l_ゲノム.Substring(230 + (5 * i), C_リード長), new string(C_高い品質, C_リード長)));
            }

            l_リード群.AddRange(Get_跨ぐリード群(l_ゲノム, 5, 0));
            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(2, l_統計.A_対象数);
            Assert.Equal(2, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 間の塩基がどのリードも品質 15 で確かな票が無くても、5 本が揃えば品質の合計で決まり G に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_間がどのリードも品質15_品質の合計でGに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            var l_リード群 = Get_間の品質を下げたリード群(l_ゲノム, 5, -1);
            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 品質 15 の列で、ある塩基が 3 本・別の塩基が 1 本なら (品質の合計 45 対 15 で 8 割に届かない) 埋めず、合意が弱い数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_品質15の列が割れる_埋めない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            var l_リード群 = Get_間の品質を下げたリード群(l_ゲノム, 4, 3);
            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_合意が弱い数);
        }

        /// <summary>
        /// 逆向きのリードが N (品質 2) を含んでいても、逆相補にして置き、G に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_逆向きのリードがNを含む_Gに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_右の始まり..]);
            List<(string A_配列, string A_品質)> l_リード群 = [];
            foreach (var (l_塩基列, l_品質列) in Get_跨ぐリード群(l_ゲノム, 5, 5))
            {
                var l_塩基 = l_塩基列.ToCharArray();
                var l_品質 = l_品質列.ToCharArray();
                l_塩基[10] = Consts.塩基コード.N;
                l_品質[10] = C_低い品質;
                l_リード群.Add((new string(l_塩基), new string(l_品質)));
            }

            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 左 = G[..520)、右 = G[480..) で跨ぐリードが 5 本あっても、左にだけ置けて右と食い違うリードが 12 本あれば埋めず、片側のリードが食い違う数は 1
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_片側のリードが12本食い違う_埋めない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_重なる左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_重なる右の始まり..]);
            var l_リード群 = Get_跨ぐリード群(l_ゲノム, 5, 0);
            l_リード群.AddRange(Get_食い違うリード群(l_ゲノム, 12));
            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_食い違い数);
        }

        /// <summary>
        /// 食い違うリードが 9 本なら (下限 10 に届かない)、G と同じ 1 本に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_片側のリードが9本食い違う_Gに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_重なる左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_重なる右の始まり..]);
            var l_リード群 = Get_跨ぐリード群(l_ゲノム, 5, 0);
            l_リード群.AddRange(Get_食い違うリード群(l_ゲノム, 9));
            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
            Assert.Equal(0, l_統計.A_食い違い数);
        }

        /// <summary>
        /// 食い違うリード 12 本に、左にだけ置けて右と確かな塩基で合うリードが 4 本あれば (12 &lt; 3 × 5)、G と同じ 1 本に埋まる
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_食い違う12本と合う4本_Gに埋まる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..C_重なる左の終わり], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_重なる右の始まり..]);
            var l_リード群 = Get_跨ぐリード群(l_ゲノム, 5, 0);
            l_リード群.AddRange(Get_食い違うリード群(l_ゲノム, 12));
            for (var i = 0; i < 4; i++)
            {
                var l_塩基 = l_ゲノム[440..560].ToCharArray();
                var l_品質 = new string(C_高い品質, l_塩基.Length).ToCharArray();
                for (var l_位置 = C_重なる右の始まり; l_位置 < 560; l_位置 += 10)
                {
                    l_塩基[l_位置 - 440] = Get_別の塩基(l_塩基[l_位置 - 440]);
                    l_品質[l_位置 - 440] = C_低い品質;
                }

                l_リード群.Add((new string(l_塩基), new string(l_品質)));
            }

            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_埋めた数);
        }

        /// <summary>
        /// 150 塩基の片 G[530..680) が右の片 G[500..) の範囲に収まる (D = −180 で、重なりが片より長い) ときは埋めず、片が重なりに収まる数は 1 (J1c までは埋まり、G[500..530) が落ちていた)
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_片が重なりに収まる_埋めない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列 = string.Concat(l_ゲノム[..400], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[530..680], Consts.未確認の繋ぎ目.ToString(), l_ゲノム[C_左の終わり..]);
            List<(string A_配列, string A_品質)> l_リード群 = [];
            for (var i = 0; i < 5; i++)
            {
                l_リード群.Add((l_ゲノム.Substring(450 + (5 * i), C_リード長), new string(C_高い品質, C_リード長)));
            }

            var (l_結果, l_統計) = Get_実行(l_配列, l_リード群);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(2, l_統計.A_対象数);
            Assert.Equal(0, l_統計.A_埋めた数);
            Assert.Equal(1, l_統計.A_跨ぐリード無し数);
            Assert.Equal(1, l_統計.A_片が収まる数);
        }

        /// <summary>
        /// 継ぎ目の無い配列はそのまま返し、リードを読まない
        /// </summary>
        [Fact]
        public void Get_埋めた配列群_継ぎ目なし_そのままでリードを読まない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_Is読んだ = false;
            var l_結果 = JunctionConsensusFiller.Get_埋めた配列群([(C_項目_S, l_ゲノム)], () =>
            {
                l_Is読んだ = true;
                return [];
            }, out var l_統計);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
            Assert.False(l_Is読んだ);
            Assert.Equal(0, l_統計.A_対象数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 本の配列とリード (品質は Phred+33) で埋める処理を実行する
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_リード群">リードの配列と品質</param>
        /// <returns>埋めた後の配列群と統計</returns>
        private static (List<(string A_ID, string A_配列)> A_結果, 継ぎ目合意統計 A_統計) Get_実行(string p_配列, List<(string A_配列, string A_品質)> p_リード群)
        {
            var l_結果 = JunctionConsensusFiller.Get_埋めた配列群([(C_項目_S, p_配列)], () => p_リード群.Select(static x => (x.A_配列, x.A_品質, C_オフセット)), out var l_統計);
            return (l_結果, l_統計);
        }

        /// <summary>
        /// G[410..) から 10 塩基ずつずらした、G[480..550) を含む長さ 150 のリードを作る (後ろの p_逆向きの数 本は逆相補にする)
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_本数">リードの本数</param>
        /// <param name="p_逆向きの数">逆相補にする本数</param>
        /// <returns>リードの配列と品質 (すべて品質 40)</returns>
        private static List<(string A_配列, string A_品質)> Get_跨ぐリード群(string p_ゲノム, int p_本数, int p_逆向きの数)
        {
            List<(string A_配列, string A_品質)> l_結果 = [];
            for (var i = 0; i < p_本数; i++)
            {
                var l_リード = p_ゲノム.Substring(410 + (10 * (i % 5)) + (i / 5), C_リード長);
                l_結果.Add((i >= p_本数 - p_逆向きの数 ? Util.V_逆相補(l_リード) : l_リード, new string(C_高い品質, C_リード長)));
            }

            return l_結果;
        }

        /// <summary>
        /// G[410..) から 10 塩基ずつずらした長さ 150 のリードで、間 G[500..530) の塩基の品質を 15 にしたものを作る (p_塩基を変える番号 番目のリードだけ、G[515] を別の塩基にする)
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_本数">リードの本数</param>
        /// <param name="p_塩基を変える番号">G[515] を別の塩基にするリードの番号、変えなければ -1</param>
        /// <returns>リードの配列と品質</returns>
        private static List<(string A_配列, string A_品質)> Get_間の品質を下げたリード群(string p_ゲノム, int p_本数, int p_塩基を変える番号)
        {
            List<(string A_配列, string A_品質)> l_結果 = [];
            for (var i = 0; i < p_本数; i++)
            {
                var l_始め = 410 + (10 * i);
                var l_塩基 = p_ゲノム.Substring(l_始め, C_リード長).ToCharArray();
                var l_品質 = new string(C_高い品質, C_リード長).ToCharArray();
                for (var l_位置 = C_左の終わり; l_位置 < C_右の始まり; l_位置++)
                {
                    l_品質[l_位置 - l_始め] = C_品質15;
                }

                if (i == p_塩基を変える番号)
                {
                    l_塩基[515 - l_始め] = Get_別の塩基(l_塩基[515 - l_始め]);
                }

                l_結果.Add((new string(l_塩基), new string(l_品質)));
            }

            return l_結果;
        }

        /// <summary>
        /// G[400..520) の後に、G[520..600) とは違う乱数の 80 塩基を付けたリード (すべて品質 40) を作る
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_本数">リードの本数</param>
        /// <returns>リードの配列と品質</returns>
        private static List<(string A_配列, string A_品質)> Get_食い違うリード群(string p_ゲノム, int p_本数)
        {
            var l_乱数 = new Random(C_食い違う続きの種);
            var l_続き = new char[80];
            for (var i = 0; i < l_続き.Length; i++)
            {
                l_続き[i] = Consts.塩基文字[l_乱数.Next(4)];
            }

            var l_リード = string.Concat(p_ゲノム[400..C_重なる左の終わり], new string(l_続き));
            return [.. Enumerable.Repeat((l_リード, new string(C_高い品質, l_リード.Length)), p_本数)];
        }

        /// <summary>
        /// 決まった乱数の種で G を作る
        /// </summary>
        /// <returns>ゲノム</returns>
        private static string Get_ゲノム()
        {
            var l_乱数 = new Random(20261109);
            var l_配列 = new char[C_ゲノム長];
            for (var i = 0; i < C_ゲノム長; i++)
            {
                l_配列[i] = Consts.塩基文字[l_乱数.Next(4)];
            }

            return new string(l_配列);
        }

        /// <summary>
        /// 別の塩基を返す
        /// </summary>
        /// <param name="p_塩基">元の塩基</param>
        /// <returns>別の塩基</returns>
        private static char Get_別の塩基(char p_塩基)
        {
            return p_塩基 == 'A' ? 'C' : 'A';
        }

        #endregion
    }
}
