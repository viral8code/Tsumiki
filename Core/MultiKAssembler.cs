using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 複数の k 長でアセンブリし、リファレンス無しの評価で最良のものを選ぶ。
    ///
    /// 最適な k はゲノムの反復配列の量で決まる。反復が短ければ k を上げるほど
    /// 跨げて有利になり、長ければ跨げないまま k-mer カバレッジとゲノム被覆を
    /// 失うだけになる。反復の量はリードからは事前に分からないため、試すしかない。
    /// </summary>
    internal static class MultiKAssembler
    {
        private const int 試すk長の下限 = 21;

        /// <summary>
        /// 複数の k で実行し、最良の結果を返す。
        /// どの k でもアセンブリできなかった場合は null。
        /// </summary>
        public static アセンブリ実行結果? Get_実行結果(
            Parameters p_引数, string p_一時ディレクトリ, int? p_リード長)
        {
            var l_k候補 = Get_k候補一覧(p_引数, p_リード長);
            Console.WriteLine($"[Multi-k] Trying k = {string.Join(", ", l_k候補)}");

            var l_実行結果一覧 = new List<アセンブリ実行結果>();
            アセンブリ実行結果? l_直前 = null;

            // 前段の配列を次の k へ渡す。k を上げるとカバレッジが痩せて
            // グラフが千切れるが、前段は既にその領域を通っている。
            List<引き継ぎ配列> l_引き継ぎ = [];
            List<引き継ぎ配列> l_次への引き継ぎ = [];

            foreach (var l_k長 in l_k候補)
            {
                if (Get_薄すぎるか(l_直前, l_k長, p_リード長, p_引数, out var l_予測))
                {
                    Console.WriteLine(
                        $"[Multi-k] Skipping k={l_k長}: predicted k-mer coverage {l_予測:F1}x is below " +
                        $"{Consts.マルチkの最小kmerカバレッジ:F1}x, so this k cannot produce a usable graph.");
                    continue;
                }

                Console.WriteLine();
                Console.WriteLine($"[Multi-k] ===== Assembling with k={l_k長} =====");
                var l_結果 = AssemblyPipeline.Get_実行結果(
                    p_引数, l_k長, p_一時ディレクトリ, $"k{l_k長}_", p_リード長,
                    p_引数.A_引き継ぐか ? l_引き継ぎ : null,
                    p_引数.A_引き継ぐか ? l_次への引き継ぎ : null);
                if (l_結果 is null)
                {
                    Console.WriteLine($"[Multi-k] k={l_k長} produced no assembly; skipping it.");
                    continue;
                }
                l_実行結果一覧.Add(l_結果);
                l_直前 = l_結果;
                l_引き継ぎ = [.. l_次への引き継ぎ];
            }

            if (l_実行結果一覧.Count == 0)
            {
                return null;
            }
            if (l_実行結果一覧.Count == 1)
            {
                Console.WriteLine("[Multi-k] Only one k produced an assembly; using it without comparison.");
                V_複製_採用した結果(l_実行結果一覧[0]);
                return l_実行結果一覧[0];
            }

            var l_候補 = Get_評価済み候補(p_引数, l_実行結果一覧, l_k候補[0], p_一時ディレクトリ, p_リード長);
            if (l_候補.Count == 0)
            {
                // 評価できない以上、根拠のある選択はできない。
                Console.WriteLine("[Multi-k] Could not evaluate the candidates; falling back to the largest k.");
                var l_代替 = l_実行結果一覧[^1];
                V_複製_採用した結果(l_代替);
                return l_代替;
            }

            var l_最良 = AssemblySelector.Get_最良(l_候補)!.Value;
            AssemblySelector.V_出力_候補一覧(l_候補, l_最良.A_実行結果);
            Console.WriteLine($"[Multi-k] Selected k={l_最良.A_実行結果.A_k長}.");

            var l_採用 = (p_引数.A_マージするか
                    ? Get_統合結果(p_引数, l_最良, l_候補, l_k候補[0], p_一時ディレクトリ, p_リード長)
                    : null)
                ?? l_最良.A_実行結果;
            V_複製_採用した結果(l_採用);
            return l_採用;
        }

        /// <summary>
        /// 骨格に他の k の配列を統合し、良くなっていれば統合結果を返す。
        /// 良くならなければ null を返して骨格をそのまま使う。
        ///
        /// 統合は誤った連結を持ち込みうるので、必ず同じ物差しで測り直して
        /// 骨格に勝ったときだけ採る。勝敗の判定は候補選びと同じ規則に任せる。
        /// </summary>
        private static アセンブリ実行結果? Get_統合結果(
            Parameters p_引数,
            (アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価) p_最良,
            List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補,
            int p_アンカーk長,
            string p_一時ディレクトリ,
            int? p_リード長)
        {
            Console.WriteLine();
            Console.WriteLine("[Merge] Merging the other k values into the selected assembly");

            var l_統合パス = "merged_" + Consts.スキャフォールドファイル名;
            var l_全候補 = p_候補.Select(x => x.A_実行結果).ToList();
            if (!AssemblyMerger.V_統合(p_最良.A_実行結果, l_全候補, p_アンカーk長, l_統合パス))
            {
                return null;
            }

            var l_統合結果 = p_最良.A_実行結果 with
            {
                A_コンティグパス = l_統合パス,
                A_スキャフォールドパス = l_統合パス,
            };

            var l_統合の評価 = Get_統合の評価(
                p_引数, l_統合結果, p_アンカーk長, p_一時ディレクトリ, p_リード長);
            if (l_統合の評価 is null)
            {
                Console.WriteLine("[Merge] Could not evaluate the merged assembly; keeping the selected one.");
                return null;
            }

            Console.WriteLine($"[Merge]   before: {p_最良.A_評価}");
            Console.WriteLine($"[Merge]   merged: {l_統合の評価}");

            var l_勝者 = AssemblySelector.Get_最良([p_最良, (l_統合結果, l_統合の評価)])!.Value;
            if (l_勝者.A_実行結果.A_最終パス != l_統合パス)
            {
                Console.WriteLine("[Merge] The merged assembly did not beat the selected one; keeping the selected one.");
                return null;
            }

            Console.WriteLine("[Merge] Using the merged assembly.");
            return l_統合結果;
        }

        /// <summary>
        /// 統合結果を、候補と同じアンカー k-mer 集合で測り直す。
        /// アンカーは選択時に作ったものと同じ条件で作り直す。
        /// </summary>
        private static アセンブリ評価? Get_統合の評価(
            Parameters p_引数,
            アセンブリ実行結果 p_統合結果,
            int p_アンカーk長,
            string p_一時ディレクトリ,
            int? p_リード長)
        {
            var l_作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"anchor{p_アンカーk長}_merge");
            _ = Directory.CreateDirectory(l_作業ディレクトリ);

            p_引数.Set_推定k長(p_アンカーk長);
            using var l_アンカー = new TrustedKmerIndex(l_作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_アンカー;

            V_読込_リード(p_引数, l_アンカー);
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_アンカー);
            _ = l_アンカー.V_カットオフ(p_引数.A_kmerカットオフ);

            var l_解析 = KmerHistogram.Get_解析結果(l_アンカー.A_出現回数ヒストグラム);
            return l_解析 is null
                ? null
                : AssemblyScorer.Get_評価(
                    p_統合結果.A_最終パス, l_アンカー, p_アンカーk長,
                    l_解析.A_ピーク出現回数, l_解析.A_推定ゲノムサイズ);
        }

        /// <summary>
        /// この k ではカバレッジが薄すぎて試すだけ無駄か。
        /// 判断できる材料が無い(まだ1つも走っていない、リード長が不明、
        /// -k で明示指定された)場合は捨てない。
        /// </summary>
        private static bool Get_薄すぎるか(
            アセンブリ実行結果? p_直前, int p_k長, int? p_リード長, Parameters p_引数, out double p_予測)
        {
            p_予測 = 0;
            if (p_直前 is null || p_リード長 is not { } l_リード長 || p_引数.A_k長一覧.Count > 0)
            {
                return false;
            }

            p_予測 = Get_予測kmerカバレッジ(
                p_直前.A_単一コピー基準値, p_直前.A_k長, p_k長, l_リード長);
            return p_予測 < Consts.マルチkの最小kmerカバレッジ;
        }

        /// <summary>
        /// 全候補を、共通のアンカー k-mer 集合に対して評価する。
        ///
        /// k が違えば k-mer 集合の大きさも意味も変わるため、各アセンブリを
        /// 自身の k で測ったのでは比較にならない。カウントのパスが1回増えるが、
        /// 共通の物差しが無ければ比較そのものが成立しない。
        /// </summary>
        private static List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> Get_評価済み候補(
            Parameters p_引数,
            List<アセンブリ実行結果> p_実行結果一覧,
            int p_アンカーk長,
            string p_一時ディレクトリ,
            int? p_リード長)
        {
            Console.WriteLine();
            Console.WriteLine($"[Multi-k] Building the common anchor k-mer set (k={p_アンカーk長}) for comparison");

            var l_作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"anchor{p_アンカーk長}");
            _ = Directory.CreateDirectory(l_作業ディレクトリ);

            p_引数.Set_推定k長(p_アンカーk長);
            using var l_アンカー = new TrustedKmerIndex(l_作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_アンカー;

            V_読込_リード(p_引数, l_アンカー);

            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_アンカー);
            _ = l_アンカー.V_カットオフ(p_引数.A_kmerカットオフ);
            KmerHistogram.V_出力_スペクトル(l_アンカー.A_出現回数ヒストグラム, p_アンカーk長, p_リード長);

            // 山の位置が単一コピーのカバレッジ、面積÷山がゲノムサイズになる。
            // 前者はコピー数の換算に、後者は NG50 の分母に使う。
            var l_解析 = KmerHistogram.Get_解析結果(l_アンカー.A_出現回数ヒストグラム);
            if (l_解析 is null)
            {
                Console.WriteLine("[Multi-k] The anchor k-mer spectrum is not bimodal; candidates cannot be compared.");
                return [];
            }

            var l_候補 = new List<(アセンブリ実行結果, アセンブリ評価)>();
            foreach (var l_実行結果 in p_実行結果一覧)
            {
                var l_評価 = AssemblyScorer.Get_評価(
                    l_実行結果.A_最終パス, l_アンカー, p_アンカーk長,
                    l_解析.A_ピーク出現回数, l_解析.A_推定ゲノムサイズ);
                if (l_評価 is not null)
                {
                    l_候補.Add((l_実行結果, l_評価));
                }
            }
            return l_候補;
        }

        private static void V_読込_リード(Parameters p_引数, TrustedKmerIndex p_kmerインデックス)
        {
            if (p_引数.A_曖昧塩基を許容するか)
            {
                KmerCounting.V_読込_リードファイル_曖昧塩基あり(p_引数.A_リード1のパス, p_kmerインデックス);
            }
            else
            {
                KmerCounting.V_読込_リードファイル(p_引数.A_リード1のパス, p_kmerインデックス);
            }

            if (string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                return;
            }

            if (p_引数.A_曖昧塩基を許容するか)
            {
                KmerCounting.V_読込_リードファイル_曖昧塩基あり(p_引数.A_リード2のパス, p_kmerインデックス);
            }
            else
            {
                KmerCounting.V_読込_リードファイル(p_引数.A_リード2のパス, p_kmerインデックス);
            }
        }

        /// <summary>
        /// 採用した k の生成物を接頭辞の無い名前へ複製する。
        /// 各 k の生成物は、選択の妥当性を後から確かめられるよう残す。
        /// </summary>
        private static void V_複製_採用した結果(アセンブリ実行結果 p_結果)
        {
            V_複製(p_結果.A_ユニティグパス, Consts.ユニティグファイル名);
            V_複製(p_結果.A_コンティグパス, Consts.コンティグファイル名);
            if (p_結果.A_スキャフォールドパス is { } l_スキャフォールドパス)
            {
                V_複製(l_スキャフォールドパス, Consts.スキャフォールドファイル名);
            }
        }

        private static void V_複製(string p_元, string p_先)
        {
            if (p_元 != p_先 && File.Exists(p_元))
            {
                File.Copy(p_元, p_先, overwrite: true);
            }
        }

        /// <summary>
        /// 試す k の一覧。-k にカンマ区切りで指定されていればそれをそのまま使う。
        ///
        /// 自動の場合は 21 からリード長の <see cref="Consts.マルチk上限のリード長比"/> 倍までを
        /// 等比で刻む。等比にするのは、k を 21 から 42 に上げたときと 110 から 131 に
        /// 上げたときでは、跨げるようになる反復配列の範囲が桁で違うため。
        /// 上限をリード長近くまで取るのは、カバレッジが十分あれば
        /// リード長に近い k のほうが良い場合があるため。
        /// </summary>
        public static List<int> Get_k候補一覧(Parameters p_引数, int? p_リード長)
        {
            if (p_引数.A_k長一覧.Count > 0)
            {
                return [.. p_引数.A_k長一覧];
            }

            if (p_リード長 is not { } l_リード長)
            {
                return [p_引数.A_k長];
            }

            var l_上限 = Get_奇数((int)(l_リード長 * Consts.マルチk上限のリード長比));
            var l_下限 = Consts.マルチkの下限;
            if (l_下限 >= l_上限)
            {
                return [Get_奇数(Math.Min(l_上限, l_リード長 - 1))];
            }

            var l_比 = Math.Pow((double)l_上限 / l_下限, 1.0 / (Consts.マルチkで試す個数 - 1));
            var l_候補 = new SortedSet<int>();
            for (var i = 0; i < Consts.マルチkで試す個数; i++)
            {
                _ = l_候補.Add(Get_奇数((int)Math.Round(l_下限 * Math.Pow(l_比, i))));
            }
            return [.. l_候補];
        }

        /// <summary>
        /// 直前の k での単一コピーカバレッジから、次の k でのカバレッジを予測する。
        /// 1リードから取れる k-mer は リード長 - k + 1 本なので、その比で縮む。
        /// </summary>
        public static double Get_予測kmerカバレッジ(
            double p_直前の基準値, int p_直前のk長, int p_次のk長, int p_リード長)
        {
            var l_直前の本数 = p_リード長 - p_直前のk長 + 1;
            var l_次の本数 = p_リード長 - p_次のk長 + 1;
            return l_直前の本数 <= 0 || l_次の本数 <= 0
                ? 0
                : p_直前の基準値 * l_次の本数 / l_直前の本数;
        }

        /// <summary>
        /// 奇数へ切り下げる。偶数の k は k-mer 自身がその逆相補と一致しうるため、
        /// 正規形が縮退して隣接判定が壊れる。
        /// </summary>
        private static int Get_奇数(int p_値)
        {
            return p_値 % 2 == 0 ? p_値 - 1 : p_値;
        }
    }
}
