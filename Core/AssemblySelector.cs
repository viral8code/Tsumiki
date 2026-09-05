using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// 複数のアセンブリ候補から1つを選ぶ。
    ///
    /// 連続性・完全性・正確性を掛け合わせた単一のスコアでは選べない。
    /// 連続性の利得が完全性の損失を上回りうるため、配列を大きく落とした
    /// 誤アセンブリのほうが高い点になる。順序を明示した二段階にしてある。
    /// </summary>
    internal static class AssemblySelector
    {
        /// <summary>
        /// 足切りに使う完全性の許容差。k を変えたときに生じる差はこの範囲に収まり、
        /// 配列を飛ばした誤アセンブリはこれよりはるかに大きく落ちる。
        /// </summary>
        public const double 完全性の許容差 = 0.05;

        public const double 正確性の許容差 = 0.05;

        /// <summary>
        /// 候補から最良のものを選ぶ。候補が空なら null。
        /// </summary>
        public static (アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)? Get_最良(
            IReadOnlyList<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補)
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

            // 配列を落として連続性を買う取引を、許容差を超えては認めない。
            var l_残った候補 = p_候補
                .Where(x => x.A_評価.A_完全性 >= l_最良の完全性 - 完全性の許容差)
                .Where(x => x.A_評価.A_正確性 >= l_最良の正確性 - 正確性の許容差)
                .ToList();

            if (l_残った候補.Count == 0)
            {
                l_残った候補 = [.. p_候補];
            }

            return l_残った候補
                .OrderByDescending(x => x.A_評価.A_NG50)
                .ThenByDescending(x => x.A_評価.A_完全性)
                .First();
        }

        /// <summary>
        /// 候補の一覧を出力する。自動選択の妥当性を利用者が確かめられるようにする。
        /// </summary>
        public static void V_出力_候補一覧(
            IReadOnlyList<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補,
            アセンブリ実行結果 p_採用したもの)
        {
            Console.WriteLine("[Multi-k] Candidate assemblies (reference-free evaluation against a common anchor k-mer set):");
            foreach (var (l_実行結果, l_評価) in p_候補.OrderBy(x => x.A_実行結果.A_k長))
            {
                var l_印 = l_実行結果.A_k長 == p_採用したもの.A_k長 ? " <- selected" : string.Empty;
                Console.WriteLine($"[Multi-k]   k={l_実行結果.A_k長,3}: {l_評価}{l_印}");
            }
        }
    }
}
