using Tsumiki.Commons;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 複数のアセンブリ候補から 1 つを選ぶ
    /// </summary>
    internal static class AssemblySelector
    {
        #region 定数

        /// <summary>
        /// 足切りに使う完全性の許容差
        /// </summary>
        /// <remarks>
        /// 取りこぼしはそのままゲノムの欠落なので、正確性より厳しく見る (7 Mbp 級なら 1% で 70 kbp に相当する)
        /// </remarks>
        public const double 完全性の許容差 = 0.01D;

        /// <summary>
        /// 足切りに使う正確性の許容差
        /// </summary>
        public const double 正確性の許容差 = 0.05D;

        /// <summary>
        /// 完全性・正確性で同点とみなす差
        /// </summary>
        /// <remarks>
        /// どちらもカバレッジからの期待コピー数の丸めに依存するため、この程度の差は候補の優劣ではなく推定の揺らぎとみなす<br/>
        /// 特に正確性は、反復を正しく複製したときにも (丸めが 1 つ下に落ちれば) 下がる向きに動く
        /// </remarks>
        public const double 同点とみなす差 = 0.005D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 候補から最良のものを選ぶ
        /// </summary>
        /// <param name="p_候補"></param>
        /// <remarks>
        /// 候補が空なら null
        /// </remarks>
        /// <returns></returns>
        public static (アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)? Get_最良(IReadOnlyList<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補)
        {
            if (p_候補.Count == 0)
            {
                return null;
            }

            if (p_候補.Count == 1)
            {
                return p_候補[0];
            }

            var l_最良の完全性 = p_候補.Max(x => x.A_評価.A_完全性);
            var l_最良の正確性 = p_候補.Max(x => x.A_評価.A_正確性);

            // 配列を落として連続性を買う取引を、許容差を超えては認めない
            var l_残った候補 = p_候補
                .Where(x => x.A_評価.A_完全性 >= l_最良の完全性 - 完全性の許容差)
                .Where(x => x.A_評価.A_正確性 >= l_最良の正確性 - 正確性の許容差)
                .ToList();

            if (l_残った候補.Count == 0)
            {
                l_残った候補 = [.. p_候補];
            }

            // 足切り後の最良値を基準に段を切る
            // 段が同じなら差は無かったものとして
            // 次の観点へ進む
            // 基準を全候補ではなく残った候補から取り直すのは、
            // 足切りで落ちたものに段の刻み位置を左右されないようにするため
            var l_基準の完全性 = l_残った候補.Max(x => x.A_評価.A_完全性);
            var l_基準の正確性 = l_残った候補.Max(x => x.A_評価.A_正確性);

            return l_残った候補
                .OrderByDescending(x => x.A_評価.A_環状本数)
                .ThenByDescending(x => x.A_評価.A_環状化率)
                .ThenBy(x => Get_段(x.A_評価.A_完全性, l_基準の完全性))
                .ThenBy(x => Get_段(x.A_評価.A_正確性, l_基準の正確性))
                .ThenByDescending(x => x.A_評価.A_NG50)
                .ThenByDescending(x => x.A_評価.A_完全性)
                .ThenByDescending(x => x.A_評価.A_正確性)
                .First();
        }

        /// <summary>
        /// 基準値からどれだけ離れているかを、同点とみなす幅で刻んだ段
        /// </summary>
        /// <param name="p_値"></param>
        /// <param name="p_基準値"></param>
        /// <remarks>
        /// 0 が基準と同等で、大きいほど劣る
        /// </remarks>
        /// <returns></returns>
        public static int Get_段(double p_値, double p_基準値)
        {
            return (int)Math.Floor(Math.Max(0D, p_基準値 - p_値) / 同点とみなす差);
        }

        /// <summary>
        /// 候補の一覧を出力する
        /// </summary>
        /// <param name="p_候補"></param>
        /// <param name="p_採用したもの"></param>
        /// <remarks>
        /// 自動選択の妥当性を利用者が確かめられるようにする
        /// </remarks>
        public static void V_出力_候補一覧(IReadOnlyList<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補, アセンブリ実行結果 p_採用したもの)
        {
            Logger.V_出力(メッセージID.候補一覧の見出し);
            foreach (var (l_実行結果, l_評価) in p_候補.OrderBy(x => x.A_実行結果.A_k長))
            {
                var l_印 = l_実行結果.A_k長 == p_採用したもの.A_k長 ? " <- selected" : string.Empty;
                Logger.V_出力(メッセージID.候補一覧の明細, l_実行結果.A_k長, l_評価, l_印);
            }
        }

        #endregion
    }
}
