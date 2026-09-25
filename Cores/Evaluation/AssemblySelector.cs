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
        public const double C_完全性の許容差 = 0.01D;

        /// <summary>
        /// 足切りに使う正確性の許容差
        /// </summary>
        public const double C_正確性の許容差 = 0.05D;

        /// <summary>
        /// 完全性・正確性で同点とみなす差
        /// </summary>
        public const double C_同点とみなす差 = 0.005D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 候補から最良のものを選ぶ
        /// </summary>
        /// <param name="p_候補"></param>
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

            var l_残った候補 = p_候補
                .Where(x => x.A_評価.A_完全性 >= l_最良の完全性 - C_完全性の許容差)
                .Where(x => x.A_評価.A_正確性 >= l_最良の正確性 - C_正確性の許容差)
                .ToList();

            if (l_残った候補.Count == 0)
            {
                l_残った候補 = [.. p_候補];
            }

            var l_完全性の段 = Get_段表(l_残った候補.Select(x => x.A_評価.A_完全性));
            var l_正確性の段 = Get_段表(l_残った候補.Select(x => x.A_評価.A_正確性));

            return l_残った候補
                .OrderByDescending(x => x.A_評価.A_環状本数)
                .ThenByDescending(x => x.A_評価.A_環状化率)
                .ThenBy(x => l_完全性の段[x.A_評価.A_完全性])
                .ThenBy(x => l_正確性の段[x.A_評価.A_正確性])
                .ThenByDescending(x => x.A_評価.A_NG50)
                .ThenByDescending(x => x.A_評価.A_完全性)
                .ThenByDescending(x => x.A_評価.A_正確性)
                .First();
        }

        /// <summary>
        /// 値ごとの段を、良い順に並べて隣との差が同点幅を超えたところで切って作る
        /// </summary>
        /// <param name="p_値群">候補の値</param>
        /// <returns>値から段を引く表</returns>
        public static Dictionary<double, int> Get_段表(IEnumerable<double> p_値群)
        {
            var l_降順 = p_値群.Distinct().OrderByDescending(x => x).ToList();
            Dictionary<double, int> l_表 = [];
            var l_段 = 0;
            for (var i = 0; i < l_降順.Count; i++)
            {
                if (i > 0 && l_降順[i - 1] - l_降順[i] > C_同点とみなす差)
                {
                    l_段++;
                }

                l_表[l_降順[i]] = l_段;
            }

            return l_表;
        }

        /// <summary>
        /// 候補の一覧を出力する
        /// </summary>
        /// <param name="p_候補"></param>
        /// <param name="p_採用したもの"></param>
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
