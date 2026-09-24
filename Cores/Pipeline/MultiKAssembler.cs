using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 複数の k 長でアセンブリし、リファレンス無しの評価で最良のものを選ぶ
    /// </summary>
    internal static class MultiKAssembler
    {
        #region 定数

        /// <summary>
        /// 自動で試す k のリード長に対する上限比
        /// </summary>
        private const double マルチk上限のリード長比 = 0.9D;

        /// <summary>
        /// 自動で試す k の下限
        /// </summary>
        private const int マルチkの下限 = 21;

        /// <summary>
        /// アンカー k を候補の最小 k から下げる量
        /// </summary>
        private const int アンカーk長の候補からの差 = 2;

        /// <summary>
        /// アンカー k の下限
        /// </summary>
        private const int アンカーk長の下限 = 11;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 複数の k で実行し、最良の結果を返す
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_原入力">反復検査に使う加工前の入力</param>
        /// <returns></returns>
        public static アセンブリ実行結果? Get_実行結果(Parameters p_引数, string p_一時ディレクトリ, int? p_リード長, Parameters? p_原入力 = null)
        {
            var l_k候補 = Get_k候補一覧(p_引数, p_リード長);
            Logger.V_出力(メッセージID.試すk一覧, string.Join(", ", l_k候補));

            var l_実行結果一覧 = new List<アセンブリ実行結果>();
            アセンブリ実行結果? l_直前 = null;

            List<引き継ぎ配列> l_引き継ぎ = [];
            List<引き継ぎ配列> l_次への引き継ぎ = [];

            List<引き継ぎ配列> l_合成リードの控え = [];

            foreach (var l_k長 in l_k候補)
            {
                if (Is薄すぎる(l_直前, l_k長, p_リード長, p_引数, out var l_予測))
                {
                    Logger.V_出力(メッセージID.kが薄すぎて省略, l_k長, l_予測, Consts.マルチkの最小kmerカバレッジ);
                    continue;
                }

                Logger.V_出力_空行();
                Logger.V_出力(メッセージID.kの開始見出し, l_k長);
                var l_結果 = AssemblyPipeline.Get_実行結果(p_引数, l_k長, p_一時ディレクトリ, p_リード長, p_引数.A_Is引き継ぎ ? l_引き継ぎ : null, Is次段引き継ぎ必要(p_引数, l_k候補, l_k長) ? l_次への引き継ぎ : null, l_合成リードの控え, p_原入力);
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

            var l_アンカーk長 = Get_アンカーk長(l_k候補);
            var l_アンカー作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"anchor{l_アンカーk長}");
            _ = Directory.CreateDirectory(l_アンカー作業ディレクトリ);

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.アンカーkmer集合の構築, l_アンカーk長);
            p_引数.Set_推定k長(l_アンカーk長);
            using var l_アンカー = new TrustedKmerIndex(l_アンカー作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_アンカー;

            KmerCounting.V_読込_リードペア(p_引数, l_アンカー);
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_アンカー);
            l_アンカー.V_適用_カットオフ(p_引数.A_kmerカットオフ);
            KmerHistogram.V_出力_スペクトル(l_アンカー.A_出現回数ヒストグラム, l_アンカーk長, p_リード長);

            var l_解析 = KmerHistogram.Get_解析結果(l_アンカー.A_出現回数ヒストグラム);
            if (l_解析 is null)
            {
                Logger.V_出力(メッセージID.アンカースペクトルが二峰でない);
                Logger.V_出力(メッセージID.候補を評価できない);
                return l_実行結果一覧[^1];
            }

            var l_候補 = Get_評価済み候補(l_実行結果一覧, l_アンカー, l_アンカーk長, l_解析);
            if (l_候補.Count == 0)
            {
                Logger.V_出力(メッセージID.候補を評価できない);
                return l_実行結果一覧[^1];
            }

            if (l_実行結果一覧.Count == 1)
            {
                Logger.V_出力(メッセージID.単一のkのみ成功);
                return l_候補[0].A_実行結果 with { A_固定アンカー評価 = l_候補[0].A_評価 };
            }

            var l_最良 = AssemblySelector.Get_最良(l_候補)!.Value;
            AssemblySelector.V_出力_候補一覧(l_候補, l_最良.A_実行結果);
            Logger.V_出力(メッセージID.採用したk, l_最良.A_実行結果.A_k長);

            var l_採用 = (p_引数.A_Isマージ
                    ? Get_統合結果(l_最良, l_候補, l_アンカー, l_アンカーk長, l_解析, p_一時ディレクトリ)
                    : null)
                ?? l_最良.A_実行結果 with { A_固定アンカー評価 = l_最良.A_評価 };
            return l_採用;
        }

        /// <summary>
        /// 単一 k 指定 (マルチk不使用) で得た結果に、固定アンカー k-mer 集合による独立評価を付ける
        /// </summary>
        /// <param name="p_結果">単一 k で得たアセンブリ実行結果</param>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_一時ディレクトリ">一時ディレクトリ</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <returns>評価を付けた結果、アンカースペクトルが二峰でない等で測れなければ元の結果をそのまま返す</returns>
        public static アセンブリ実行結果 Get_固定アンカー評価を付与(アセンブリ実行結果 p_結果, Parameters p_引数, string p_一時ディレクトリ, int? p_リード長)
        {
            var l_アンカーk長 = Get_アンカーk長([p_結果.A_k長]);
            var l_アンカー作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"anchor{l_アンカーk長}");
            _ = Directory.CreateDirectory(l_アンカー作業ディレクトリ);

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.アンカーkmer集合の構築, l_アンカーk長);
            p_引数.Set_推定k長(l_アンカーk長);
            using var l_アンカー = new TrustedKmerIndex(l_アンカー作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_アンカー;

            KmerCounting.V_読込_リードペア(p_引数, l_アンカー);
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_アンカー);
            l_アンカー.V_適用_カットオフ(p_引数.A_kmerカットオフ);
            KmerHistogram.V_出力_スペクトル(l_アンカー.A_出現回数ヒストグラム, l_アンカーk長, p_リード長);

            var l_解析 = KmerHistogram.Get_解析結果(l_アンカー.A_出現回数ヒストグラム);
            if (l_解析 is null)
            {
                Logger.V_出力(メッセージID.アンカースペクトルが二峰でない);
                Logger.V_出力(メッセージID.候補を評価できない);
                return p_結果;
            }

            var l_評価 = AssemblyScorer.Get_評価(p_結果.A_最終パス, l_アンカー, l_アンカーk長, l_解析.A_単一コピー基準値, l_解析.A_推定ゲノムサイズ, l_解析.A_単一コピー上限);
            if (l_評価 is null)
            {
                Logger.V_出力(メッセージID.候補を評価できない);
                return p_結果;
            }
            return p_結果 with { A_固定アンカー評価 = l_評価 };
        }

        /// <summary>
        /// 試す k の一覧
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_リード長">リード長、不明なら null</param>
        /// <returns>試す k の一覧</returns>
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

            var l_上限 = Get_奇数((int)(l_リード長 * マルチk上限のリード長比));
            var l_下限 = マルチkの下限;
            if (l_下限 >= l_上限)
            {
                return [Get_奇数(Math.Min(l_上限, l_リード長 - 1))];
            }

            var l_比 = Math.Pow((double)l_上限 / l_下限, 1D / (Consts.マルチkで試す個数 - 1));
            var l_候補 = new SortedSet<int>();
            for (var i = 0; i < Consts.マルチkで試す個数; i++)
            {
                _ = l_候補.Add(Get_奇数((int)Math.Round(l_下限 * Math.Pow(l_比, i))));
            }
            return [.. l_候補];
        }

        /// <summary>
        /// 次段で使用する引き継ぎ配列を準備する必要があるか判定する
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_k候補"></param>
        /// <param name="p_現在k"></param>
        /// <returns></returns>
        internal static bool Is次段引き継ぎ必要(Parameters p_引数, IReadOnlyList<int> p_k候補, int p_現在k)
        {
            return p_引数.A_Is引き継ぎ && p_k候補.Any(x => x > p_現在k);
        }

        /// <summary>
        /// 候補を評価する物差しの k
        /// </summary>
        /// <param name="p_k候補"></param>
        /// <returns></returns>
        public static int Get_アンカーk長(IReadOnlyList<int> p_k候補)
        {
            var l_k長 = p_k候補[0] - アンカーk長の候補からの差;

            if (l_k長 % 2 == 0)
            {
                l_k長--;
            }
            return Math.Max(アンカーk長の下限, l_k長);
        }

        /// <summary>
        /// 直前の k での単一コピーカバレッジから、次の k でのカバレッジを予測する
        /// </summary>
        /// <param name="p_直前の基準値"></param>
        /// <param name="p_直前のk長"></param>
        /// <param name="p_次のk長"></param>
        /// <param name="p_リード長"></param>
        /// <returns></returns>
        public static double Get_予測kmerカバレッジ(double p_直前の基準値, int p_直前のk長, int p_次のk長, int p_リード長)
        {
            var l_直前の本数 = p_リード長 - p_直前のk長 + 1;
            var l_次の本数 = p_リード長 - p_次のk長 + 1;
            return l_直前の本数 <= 0 || l_次の本数 <= 0
                ? 0D
                : p_直前の基準値 * l_次の本数 / l_直前の本数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 骨格に他の k の配列を統合し、良くなっていれば統合結果を返す
        /// </summary>
        /// <param name="p_最良"></param>
        /// <param name="p_候補"></param>
        /// <param name="p_アンカー"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_解析"></param>
        /// <param name="p_一時ディレクトリ"></param>
        /// <returns></returns>
        private static アセンブリ実行結果? Get_統合結果((アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価) p_最良, List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補, TrustedKmerIndex p_アンカー, int p_アンカーk長, スペクトル解析結果 p_解析, string p_一時ディレクトリ)
        {
            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.統合開始);

            var l_統合パス = Path.Combine(p_一時ディレクトリ, "merged_" + Consts.Scaffoldファイル名);
            var l_全候補 = p_候補.Select(x => x.A_実行結果).ToList();
            if (!AssemblyMerger.Try統合(p_最良.A_実行結果, l_全候補, p_アンカーk長, l_統合パス))
            {
                return null;
            }

            var l_統合contigパス = Path.Combine(p_一時ディレクトリ, "merged_" + AssemblyPipeline.Contigファイル名);
            V_書き出し_N分割(l_統合パス, l_統合contigパス);

            var l_統合結果 = p_最良.A_実行結果 with
            {
                A_contigパス = l_統合contigパス,
                A_scaffoldパス = l_統合パス,
            };

            var l_統合の評価 = AssemblyScorer.Get_評価(l_統合結果.A_最終パス, p_アンカー, p_アンカーk長, p_解析.A_単一コピー基準値, p_解析.A_推定ゲノムサイズ, p_解析.A_単一コピー上限);
            if (l_統合の評価 is null)
            {
                Logger.V_出力(メッセージID.統合結果を評価できない);
                return null;
            }
            l_統合結果 = l_統合結果 with { A_固定アンカー評価 = l_統合の評価 };

            Logger.V_出力(メッセージID.統合前の評価, p_最良.A_評価);
            Logger.V_出力(メッセージID.統合後の評価, l_統合の評価);

            var (l_実行結果, l_評価) = AssemblySelector.Get_最良([p_最良, (l_統合結果, l_統合の評価)])!.Value;
            if (l_実行結果.A_最終パス != l_統合パス)
            {
                Logger.V_出力(メッセージID.統合が骨格に勝てず);
                return null;
            }

            Logger.V_出力(メッセージID.統合結果を採用);
            return l_統合結果;
        }

        /// <summary>
        /// 配列を N の連続で分断して書き出す
        /// </summary>
        /// <param name="p_入力パス">分断する FASTA</param>
        /// <param name="p_出力パス">書き出す FASTA</param>
        internal static void V_書き出し_N分割(string p_入力パス, string p_出力パス)
        {
            using var l_書き込み = new FastaWriter(p_出力パス);
            var l_連番 = 1;
            foreach (var (A_ID, A_配列) in FastaReader.Get_全エントリ(p_入力パス))
            {
                foreach (var l_片 in A_配列.Split('N', StringSplitOptions.RemoveEmptyEntries))
                {
                    l_書き込み.V_書き込み($"NODE{l_連番}", l_片);
                    l_連番++;
                }
            }
        }

        /// <summary>
        /// この k ではカバレッジが薄すぎて試すだけ無駄か
        /// </summary>
        /// <param name="p_直前"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_引数"></param>
        /// <param name="p_予測"></param>
        /// <returns></returns>
        private static bool Is薄すぎる(アセンブリ実行結果? p_直前, int p_k長, int? p_リード長, Parameters p_引数, out double p_予測)
        {
            p_予測 = 0D;
            if (p_直前 is null || p_リード長 is not { } l_リード長 || p_引数.A_k長一覧.Count > 0)
            {
                return false;
            }

            p_予測 = Get_予測kmerカバレッジ(p_直前.A_単一コピー基準値, p_直前.A_k長, p_k長, l_リード長);
            return p_予測 < Consts.マルチkの最小kmerカバレッジ;
        }

        /// <summary>
        /// 全候補を、共通のアンカー k-mer 集合に対して評価する
        /// </summary>
        /// <param name="p_実行結果一覧"></param>
        /// <param name="p_アンカー"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_解析"></param>
        /// <returns></returns>
        private static List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> Get_評価済み候補(List<アセンブリ実行結果> p_実行結果一覧, TrustedKmerIndex p_アンカー, int p_アンカーk長, スペクトル解析結果 p_解析)
        {
            var l_候補 = new List<(アセンブリ実行結果, アセンブリ評価)>();
            foreach (var l_実行結果 in p_実行結果一覧)
            {
                var l_評価 = AssemblyScorer.Get_評価(l_実行結果.A_最終パス, p_アンカー, p_アンカーk長, p_解析.A_単一コピー基準値, p_解析.A_推定ゲノムサイズ, p_解析.A_単一コピー上限);
                if (l_評価 is not null)
                {
                    l_候補.Add((l_実行結果, l_評価));
                }
            }
            return l_候補;
        }

        /// <summary>
        /// 奇数へ切り下げる
        /// </summary>
        /// <param name="p_値"></param>
        /// <returns></returns>
        private static int Get_奇数(int p_値)
        {
            return p_値 % 2 == 0 ? p_値 - 1 : p_値;
        }

        #endregion
    }
}
