using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 配列の端どうしを、跨ぐリードで繋ぐ処理 (ReadEndJoiner) の検証
    /// </summary>
    public class ReadEndJoinerTests
    {
        #region 定数

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        /// <summary>
        /// 項目 B
        /// </summary>
        private const string C_項目_B = "B";

        /// <summary>
        /// 項目 C
        /// </summary>
        private const string C_項目_C = "C";

        /// <summary>
        /// 項目 B2
        /// </summary>
        private const string C_項目_B2 = "B2";

        /// <summary>
        /// 項目 D
        /// </summary>
        private const string C_項目_D = "D";

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 150;

        /// <summary>
        /// リードの開始位置の間隔
        /// </summary>
        private const int C_リードの間隔 = 10;

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 6000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 間が 60 塩基空いた 2 本は、跨ぐリードで 1 本に繋がり、ゲノムと一致する
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_隙間あり_1本になりゲノムと一致する()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_ゲノム[..3000]),
                (C_項目_B, l_ゲノム[3060..])
            };
            var (l_結果, l_候補数, l_繋いだ数) = Get_実行(l_配列群, Get_リード群(l_ゲノム));
            Assert.Equal(1, l_候補数);
            Assert.Equal(1, l_繋いだ数);
            Assert.Single(l_結果);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
        }

        /// <summary>
        /// 100 塩基重なる 2 本は、1 本に繋がり、ゲノムと一致する
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_重なりあり_1本になりゲノムと一致する()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_ゲノム[..3000]),
                (C_項目_B, l_ゲノム[2900..])
            };
            var (l_結果, _, _) = Get_実行(l_配列群, Get_リード群(l_ゲノム));
            Assert.Single(l_結果);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
        }

        /// <summary>
        /// 2 本目を逆相補にして渡しても、1 本に繋がり、ゲノムかその逆相補と一致する
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_逆向き_1本になりゲノムかその逆相補と一致する()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_ゲノム[..3000]),
                (C_項目_B, Util.V_逆相補_曖昧塩基あり(l_ゲノム[3060..]))
            };
            var (l_結果, _, _) = Get_実行(l_配列群, Get_リード群(l_ゲノム));
            Assert.Single(l_結果);
            var l_配列 = l_結果[0].A_配列;
            Assert.True(l_配列 == l_ゲノム || l_配列 == Util.V_逆相補_曖昧塩基あり(l_ゲノム));
        }

        /// <summary>
        /// 3 本の鎖は、入力の順によらず 1 本に繋がり、ゲノムと一致する<br/>
        /// ID は鎖の最初の配列 (A) の ID になる
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_3本の鎖_入力の順によらず1本になりIDは先頭のAになる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_C, l_ゲノム[4030..]),
                (C_項目_A, l_ゲノム[..2000]),
                (C_項目_B, l_ゲノム[2040..4000]),
            };
            var (l_結果, _, l_繋いだ数) = Get_実行(l_配列群, Get_リード群(l_ゲノム));
            Assert.Equal(2, l_繋いだ数);
            Assert.Single(l_結果);
            Assert.Equal(C_項目_A, l_結果[0].A_ID);
            Assert.Equal(l_ゲノム, l_結果[0].A_配列);
        }

        /// <summary>
        /// 隙間を跨ぐリードが 1 本だけなら、票が足りず繋がない<br/>
        /// 配列は元のまま 2 本残る
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_跨ぐリードが1本だけ_繋がず元のまま()
        {
            var l_ゲノム = Get_ゲノム();
            var l_A = l_ゲノム[..3000];
            var l_B = l_ゲノム[3060..];
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_A),
                (C_項目_B, l_B)
            };
            var l_リード群 = Get_リード群(l_ゲノム, s => s < 2855 || s >= 2975 || s == 2940);
            var (l_結果, _, l_繋いだ数) = Get_実行(l_配列群, l_リード群);
            Assert.Equal(0, l_繋いだ数);
            Assert.Equal(2, l_結果.Count);
            Assert.Equal(l_配列群, l_結果);
        }

        /// <summary>
        /// B と同じ先頭 300 塩基を持つ別の配列 B' があると、B の先頭の 25-mer が一意でなくなり、繋がない
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_相手が2つ_繋がずAは元のまま()
        {
            var l_ゲノム = Get_ゲノム();
            var l_乱数 = new Random(777);
            var l_B2 = l_ゲノム[3060..3360] + Get_ランダムな塩基(l_乱数, 2700);
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_ゲノム[..3000]),
                (C_項目_B, l_ゲノム[3060..]),
                (C_項目_B2, l_B2),
            };
            var (l_結果, _, l_繋いだ数) = Get_実行(l_配列群, Get_リード群(l_ゲノム));
            Assert.Equal(0, l_繋いだ数);
            Assert.Contains(l_結果, x => x.A_ID == C_項目_A && x.A_配列 == l_ゲノム[..3000]);
        }

        /// <summary>
        /// 長さ 100 の配列は繋がない<br/>
        /// 入力と同じ 2 本が、同じ順で出る
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_短すぎる配列_繋がない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_ゲノム[..3000]),
                (C_項目_B, l_ゲノム[3060..3160])
            };
            var (l_結果, _, l_繋いだ数) = Get_実行(l_配列群, Get_リード群(l_ゲノム));
            Assert.Equal(0, l_繋いだ数);
            Assert.Equal(l_配列群, l_結果);
        }

        /// <summary>
        /// 同じ入力を 2 回渡すと、出力の並び・ID・配列が同一になる
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_同じ入力なら出力が同一()
        {
            var l_ゲノム = Get_ゲノム();
            var l_リード群 = Get_リード群(l_ゲノム);
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_C, l_ゲノム[4030..]),
                (C_項目_A, l_ゲノム[..2000]),
                (C_項目_B, l_ゲノム[2040..4000]),
                (C_項目_D, l_ゲノム[..3000]),
            };
            var (l_結果1, l_候補1, l_繋いだ1) = Get_実行(l_配列群, l_リード群);
            var (l_結果2, l_候補2, l_繋いだ2) = Get_実行(l_配列群, l_リード群);
            Assert.Equal(l_結果1, l_結果2);
            Assert.Equal(l_候補1, l_候補2);
            Assert.Equal(l_繋いだ1, l_繋いだ2);
        }

        /// <summary>
        /// 縦に 2 回並んだ反復 (600 塩基) の 1 コピー目で終わる配列と、2 コピー目の後から始まる配列は、端が反復なので繋がない (繋ぐと 2 コピー目を飛ばす)
        /// </summary>
        [Fact]
        public void Get_繋いだ配列群_縦に並んだ反復を挟む端_繋がない()
        {
            var l_乱数 = new Random(20240602);
            var l_反復 = Get_ランダムな塩基(l_乱数, 600);
            var l_ゲノム = Get_ランダムな塩基(l_乱数, 2000) + l_反復 + l_反復 + Get_ランダムな塩基(l_乱数, 2000);
            var l_配列群 = new List<(string A_ID, string A_配列)>
            {
                (C_項目_A, l_ゲノム[..2600]),
                (C_項目_B, l_ゲノム[3200..])
            };
            var l_リード群 = Get_リード群(l_ゲノム);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            var l_一意の出現数 = Scaffolder.Get_一意の出現数(l_索引, l_配列群.Select(x => x.A_配列));
            var l_結果 = ReadEndJoiner.Get_繋いだ配列群(l_配列群, l_索引, C_リード長, l_一意の出現数, out _, out var l_繋いだ数);
            Assert.True(l_一意の出現数 > 0);
            Assert.Equal(0, l_繋いだ数);
            Assert.Equal(l_配列群, l_結果);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を固定した乱数で、一意な 6,000 塩基のゲノムを作る
        /// </summary>
        /// <returns>ゲノム</returns>
        private static string Get_ゲノム()
        {
            return Get_ランダムな塩基(new Random(20240601), C_ゲノム長);
        }

        /// <summary>
        /// 指定した乱数で、A/C/G/T からなる塩基列を作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <param name="p_長さ">長さ</param>
        /// <returns>塩基列</returns>
        private static string Get_ランダムな塩基(Random p_乱数, int p_長さ)
        {
            const string l_塩基 = Consts.塩基文字;
            var l_配列 = new char[p_長さ];
            for (var i = 0; i < p_長さ; i++)
            {
                l_配列[i] = l_塩基[p_乱数.Next(4)];
            }

            return new string(l_配列);
        }

        /// <summary>
        /// ゲノムから、長さ 150 のリードを 10 塩基おきに取る<br/>
        /// 開始位置の 10 の位が奇数のものは逆相補にする
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_採る">開始位置を採るかどうか、省略時は全部採る</param>
        /// <returns>リード群</returns>
        private static List<string> Get_リード群(string p_ゲノム, Func<int, bool>? p_採る = null)
        {
            var l_リード群 = new List<string>();
            for (var s = 0; s + C_リード長 <= p_ゲノム.Length; s += C_リードの間隔)
            {
                if (p_採る is not null && !p_採る(s))
                {
                    continue;
                }

                var l_リード = p_ゲノム.Substring(s, C_リード長);
                l_リード群.Add((s / C_リードの間隔) % 2 == 1 ? Util.V_逆相補_曖昧塩基あり(l_リード) : l_リード);
            }

            return l_リード群;
        }

        /// <summary>
        /// リード群から索引を作り、配列群を繋ぐ
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_リード群">リード群</param>
        /// <returns>繋いだ後の配列、候補数、繋いだ数</returns>
        private static (List<(string A_ID, string A_配列)> A_結果, int A_候補数, int A_繋いだ数) Get_実行(List<(string A_ID, string A_配列)> p_配列群, List<string> p_リード群)
        {
            var l_索引 = ReadMinimizerIndex.V_構築(() => p_リード群);
            var l_結果 = ReadEndJoiner.Get_繋いだ配列群(p_配列群, l_索引, C_リード長, 0, out var l_候補数, out var l_繋いだ数);
            return (l_結果, l_候補数, l_繋いだ数);
        }

        #endregion
    }
}
