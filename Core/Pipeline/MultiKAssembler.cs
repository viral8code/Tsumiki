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
            Logger.V_出力(メッセージID.試すk一覧, string.Join(", ", l_k候補));

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
                    Logger.V_出力(メッセージID.kが薄すぎて省略, l_k長, l_予測, Consts.マルチkの最小kmerカバレッジ);
                    continue;
                }

                Logger.V_出力_空行();
                Logger.V_出力(メッセージID.kの開始見出し, l_k長);
                var l_結果 = AssemblyPipeline.Get_実行結果(
                    p_引数, l_k長, p_一時ディレクトリ, p_リード長,
                    p_引数.A_引き継ぐか ? l_引き継ぎ : null,
                    p_引数.A_引き継ぐか ? l_次への引き継ぎ : null);
                if (l_結果 is null)
                {
                    Logger.V_出力(メッセージID.kでアセンブリできず, l_k長);
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
                Logger.V_出力(メッセージID.単一のkのみ成功);
                return l_実行結果一覧[0];
            }

            // アンカーの k-mer カウントは、候補の評価と(-mg指定時の)統合評価の
            // 両方が同じ条件(同じ k・同じカットオフ)で必要とする。生リードの
            // 走査とカウントはコストが高いため、1つのインデックスを両方で
            // 使い回す(以前は統合評価のたびに同じ内容をもう一度数え直していた)。
            var l_アンカーk長 = l_k候補[0];
            var l_アンカー作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"anchor{l_アンカーk長}");
            _ = Directory.CreateDirectory(l_アンカー作業ディレクトリ);

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.アンカーkmer集合の構築, l_アンカーk長);
            p_引数.Set_推定k長(l_アンカーk長);
            using var l_アンカー = new TrustedKmerIndex(l_アンカー作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_アンカー;

            KmerCounting.V_読込_リードペア(p_引数, l_アンカー);
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_アンカー);
            _ = l_アンカー.V_カットオフ(p_引数.A_kmerカットオフ);
            KmerHistogram.V_出力_スペクトル(l_アンカー.A_出現回数ヒストグラム, l_アンカーk長, p_リード長);

            // 山の位置が単一コピーのカバレッジ、面積÷山がゲノムサイズになる。
            // 前者はコピー数の換算に、後者は NG50 の分母に使う。
            var l_解析 = KmerHistogram.Get_解析結果(l_アンカー.A_出現回数ヒストグラム);
            if (l_解析 is null)
            {
                Logger.V_出力(メッセージID.アンカースペクトルが二峰でない);
                Logger.V_出力(メッセージID.候補を評価できない);
                var l_代替 = l_実行結果一覧[^1];
                return l_代替;
            }

            var l_候補 = Get_評価済み候補(l_実行結果一覧, l_アンカー, l_アンカーk長, l_解析);
            if (l_候補.Count == 0)
            {
                // 評価できない以上、根拠のある選択はできない。
                Logger.V_出力(メッセージID.候補を評価できない);
                var l_代替 = l_実行結果一覧[^1];
                return l_代替;
            }

            var l_最良 = AssemblySelector.Get_最良(l_候補)!.Value;
            AssemblySelector.V_出力_候補一覧(l_候補, l_最良.A_実行結果);
            Logger.V_出力(メッセージID.採用したk, l_最良.A_実行結果.A_k長);

            var l_採用 = (p_引数.A_マージするか
                    ? Get_統合結果(l_最良, l_候補, l_アンカー, l_アンカーk長, l_解析, p_一時ディレクトリ)
                    : null)
                ?? l_最良.A_実行結果;
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
            (アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価) p_最良,
            List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補,
            TrustedKmerIndex p_アンカー,
            int p_アンカーk長,
            スペクトル解析結果 p_解析,
            string p_一時ディレクトリ)
        {
            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.統合開始);

            var l_統合パス = Path.Combine(p_一時ディレクトリ, "merged_" + Consts.スキャフォールドファイル名);
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

            // 統合評価は候補評価と全く同じアンカー k-mer 集合(同じ k・同じ
            // カットオフで数え終えた既存のインデックス)を使い回す。以前は
            // ここで生リードの走査・カウント・カットオフをもう一度
            // やり直しており、同一の結果を得るためだけに重複したコストを
            // 払っていた。
            var l_統合の評価 = AssemblyScorer.Get_評価(
                l_統合結果.A_最終パス, p_アンカー, p_アンカーk長,
                p_解析.A_ピーク出現回数, p_解析.A_推定ゲノムサイズ);
            if (l_統合の評価 is null)
            {
                Logger.V_出力(メッセージID.統合結果を評価できない);
                return null;
            }

            Logger.V_出力(メッセージID.統合前の評価, p_最良.A_評価);
            Logger.V_出力(メッセージID.統合後の評価, l_統合の評価);

            var l_勝者 = AssemblySelector.Get_最良([p_最良, (l_統合結果, l_統合の評価)])!.Value;
            if (l_勝者.A_実行結果.A_最終パス != l_統合パス)
            {
                Logger.V_出力(メッセージID.統合が骨格に勝てず);
                return null;
            }

            Logger.V_出力(メッセージID.統合結果を採用);
            return l_統合結果;
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
        /// 自身の k で測ったのでは比較にならない。アンカーは呼び出し側が
        /// 既に構築済みのものを渡す(-mg指定時の統合評価とも共有するため)。
        /// </summary>
        private static List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> Get_評価済み候補(
            List<アセンブリ実行結果> p_実行結果一覧,
            TrustedKmerIndex p_アンカー,
            int p_アンカーk長,
            スペクトル解析結果 p_解析)
        {
            var l_候補 = new List<(アセンブリ実行結果, アセンブリ評価)>();
            foreach (var l_実行結果 in p_実行結果一覧)
            {
                var l_評価 = AssemblyScorer.Get_評価(
                    l_実行結果.A_最終パス, p_アンカー, p_アンカーk長,
                    p_解析.A_ピーク出現回数, p_解析.A_推定ゲノムサイズ);
                if (l_評価 is not null)
                {
                    l_候補.Add((l_実行結果, l_評価));
                }
            }
            return l_候補;
        }

        /// <summary>
        /// 採用した k の生成物を接頭辞の無い名前へ複製する。
        /// 各 k の生成物は、選択の妥当性を後から確かめられるよう残す。
        /// </summary>
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
