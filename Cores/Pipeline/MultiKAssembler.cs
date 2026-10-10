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
        private const double C_マルチk上限のリード長比 = 0.9D;

        /// <summary>
        /// 統合の橋渡しを見送る、掛かる継ぎ目の誤りの確率の下限 (統合で生まれる継ぎ目は骨格の継ぎ目より誤りが 20 倍ほど多く、0.1 以上では防げる誤りが失う正しい繋がりを上回る)
        /// </summary>
        private const double C_橋渡しを見送る確率 = 0.1D;

        /// <summary>
        /// 継ぎ目が橋渡しに掛かるとみなす、橋渡しの両端からの距離
        /// </summary>
        private const int C_橋渡しに掛かる余白 = 50;

        /// <summary>
        /// 自動で試す k の下限
        /// </summary>
        private const int C_マルチkの下限 = 21;

        /// <summary>
        /// アンカー k を候補の最小 k から下げる量
        /// </summary>
        private const int C_アンカーk長の候補からの差 = 2;

        /// <summary>
        /// アンカー k の下限
        /// </summary>
        private const int C_アンカーk長の下限 = 11;

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

            var l_再開署名 = 中間データ置き場.A_Is有効 ? null : Get_再開署名(p_引数, l_k候補);
            var l_Is再利用中 = p_引数.A_Is再開 && l_再開署名 is not null;
            var l_Is合成リード保存済み = false;

            foreach (var l_k長 in l_k候補)
            {
                if (Is薄すぎる(l_直前, l_k長, p_リード長, p_引数, out var l_予測))
                {
                    Logger.V_出力(メッセージID.kが薄すぎて省略, l_k長, l_予測, Consts.マルチkの最小kmerカバレッジ);
                    continue;
                }

                var l_k作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"k{l_k長}");
                if (l_Is再利用中 && Try再利用(p_引数, l_k作業ディレクトリ, p_一時ディレクトリ, l_再開署名!, l_合成リードの控え, out var l_再利用結果, out var l_再利用引き継ぎ))
                {
                    Logger.V_出力(メッセージID.再開_kの結果を再利用, l_k長);
                    l_Is合成リード保存済み |= l_合成リードの控え.Count > 0;
                    l_実行結果一覧.Add(l_再利用結果!);
                    l_直前 = l_再利用結果;
                    l_引き継ぎ = l_再利用引き継ぎ;
                    continue;
                }

                l_Is再利用中 = false;
                if (l_再開署名 is not null && Directory.Exists(l_k作業ディレクトリ))
                {
                    MultiKCheckpoint.V_削除_k(l_k作業ディレクトリ);
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

                if (l_再開署名 is not null)
                {
                    if (!l_Is合成リード保存済み && l_合成リードの控え.Count > 0)
                    {
                        MultiKCheckpoint.V_保存_合成リード(p_一時ディレクトリ, l_再開署名, l_合成リードの控え);
                        l_Is合成リード保存済み = true;
                    }

                    var l_合成リード数 = l_合成リードの控え.Count > 0 && l_引き継ぎ.Count >= l_合成リードの控え.Count && ReferenceEquals(l_引き継ぎ[^l_合成リードの控え.Count], l_合成リードの控え[0]) ? l_合成リードの控え.Count : 0;
                    MultiKCheckpoint.V_保存_k(l_k作業ディレクトリ, l_再開署名, l_結果, l_引き継ぎ[..^l_合成リード数], l_合成リード数);
                }
            }

            RepeatRMerVerifier.V_解放_共有索引();

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
                    ? Get_統合結果(l_最良, l_候補, l_アンカー, l_アンカーk長, l_解析, p_一時ディレクトリ, (p_原入力 ?? p_引数).A_ライブラリ群)
                    : null)
                ?? l_最良.A_実行結果 with { A_固定アンカー評価 = l_最良.A_評価 };
            return Get_補完結果(l_採用, l_候補, l_アンカー, l_アンカーk長, l_解析, p_一時ディレクトリ) ?? l_採用;
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

            var l_上限 = Get_カバレッジで決めた上限(p_引数, Get_奇数((int)(l_リード長 * C_マルチk上限のリード長比)));
            var l_下限 = C_マルチkの下限;
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
        /// エラー訂正で数えた誤りの無い区間から、予測カバレッジが試す価値のある下限を保てる最大の k を返す
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_上限">リード長から決めた上限</param>
        /// <returns>決めた上限、見積もれなければ p_上限 のまま</returns>
        private static int Get_カバレッジで決めた上限(Parameters p_引数, int p_上限)
        {
            var l_度数群 = p_引数.A_無誤り区間の度数群;
            if (l_度数群.Count == 0)
            {
                return p_上限;
            }

            var l_最小k長 = l_度数群.Max(x => x.A_k長);
            for (var l_k長 = p_上限; l_k長 >= l_最小k長; l_k長 -= 2)
            {
                var l_予測 = l_度数群.Sum(x => x.Get_予測カバレッジ(l_k長));
                if (l_予測 >= Consts.マルチkの最小kmerカバレッジ)
                {
                    if (l_k長 < p_上限)
                    {
                        Logger.V_出力(メッセージID.k上限をカバレッジで決定, l_k長, l_予測, Consts.マルチkの最小kmerカバレッジ);
                    }

                    return l_k長;
                }
            }

            return p_上限;
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
            var l_k長 = p_k候補[0] - C_アンカーk長の候補からの差;

            if (l_k長 % 2 == 0)
            {
                l_k長--;
            }

            return Math.Max(C_アンカーk長の下限, l_k長);
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

        /// <summary>
        /// 配列を N の連続で分断して書き出す
        /// </summary>
        /// <param name="p_入力パス">分断する FASTA</param>
        /// <param name="p_出力パス">書き出す FASTA</param>
        public static void V_書き出し_N分割(string p_入力パス, string p_出力パス)
        {
            using var l_書き込み = new FastaWriter(p_出力パス);
            var l_連番 = 1;
            foreach (var (A_ID, A_配列) in FastaReader.Get_全エントリ(p_入力パス))
            {
                foreach (var l_片 in A_配列.Split(['N', Consts.未確認の繋ぎ目], StringSplitOptions.RemoveEmptyEntries))
                {
                    l_書き込み.V_書き込み($"NODE{l_連番}", l_片);
                    l_連番++;
                }
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k ごとの記録が前回と同じ入力・設定・k の候補一覧で作られたかを見分ける値を作る
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_k候補">試す k の一覧</param>
        /// <returns></returns>
        private static string Get_再開署名(Parameters p_引数, List<int> p_k候補)
        {
            var l_設定 = p_引数.Get_複製();
            if (!l_設定.A_Iskmerカットオフ明示指定)
            {
                l_設定.Set_推定kmerカットオフ(1UL);
            }

            return StageCheckpoint.Get_入力署名(l_設定) + "\n" + string.Join(",", p_k候補);
        }

        /// <summary>
        /// 前回の実行で保存したその k の結果と引き継ぎを読み戻す
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_k作業ディレクトリ"></param>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_署名"></param>
        /// <param name="p_合成リードの控え">まだ空なら保存した合成リードで埋める</param>
        /// <param name="p_結果"></param>
        /// <param name="p_引き継ぎ">次の k への引き継ぎ (合成リードを含む)</param>
        /// <returns>読み戻せたら true</returns>
        private static bool Try再利用(Parameters p_引数, string p_k作業ディレクトリ, string p_一時ディレクトリ, string p_署名, List<引き継ぎ配列> p_合成リードの控え, out アセンブリ実行結果? p_結果, out List<引き継ぎ配列> p_引き継ぎ)
        {
            p_引き継ぎ = [];
            if (!MultiKCheckpoint.Try読込_k(p_k作業ディレクトリ, p_署名, out p_結果, out var l_引き継ぎ, out var l_合成リード数) || p_結果 is null)
            {
                return false;
            }

            if (l_合成リード数 > 0 && p_合成リードの控え.Count == 0)
            {
                if (!MultiKCheckpoint.Try読込_合成リード(p_一時ディレクトリ, p_署名, l_合成リード数, out var l_合成リード))
                {
                    p_結果 = null;
                    return false;
                }

                p_合成リードの控え.AddRange(l_合成リード);
            }

            if (l_合成リード数 > 0 && p_合成リードの控え.Count != l_合成リード数)
            {
                p_結果 = null;
                return false;
            }

            p_引き継ぎ = l_合成リード数 > 0 ? [.. l_引き継ぎ, .. p_合成リードの控え] : l_引き継ぎ;
            p_引数.Set_推定k長(p_結果.A_k長);
            p_引数.Set_推定kmerカットオフ(p_結果.A_kmerカットオフ);
            AmbiguityRecorder.V_読込(p_k作業ディレクトリ, p_結果.A_k長);
            return true;
        }

        /// <summary>
        /// 配列の位置から始まるアンカー k-mer の、リードのカバレッジで見た期待コピー数を返す関数を作る
        /// </summary>
        /// <param name="p_アンカー"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_解析"></param>
        /// <returns></returns>
        private static Func<string, int, int> Get_期待コピー数(TrustedKmerIndex p_アンカー, int p_アンカーk長, スペクトル解析結果 p_解析)
        {
            return (p_配列, p_位置) =>
            {
                Span<byte> l_kmer = stackalloc byte[p_アンカーk長];
                for (var i = 0; i < p_アンカーk長; i++)
                {
                    l_kmer[i] = Util.Get_塩基ID(p_配列[p_位置 + i]);
                }

                return KmerHistogram.Get_期待コピー数(p_アンカー.Get_カバレッジ(l_kmer), p_解析.A_単一コピー基準値, p_解析.A_単一コピー上限);
            };
        }

        /// <summary>
        /// 統合した配列の継ぎ目を元リードで評価し、誤りの確率が 橋渡しを見送る確率 以上の継ぎ目に掛かる橋渡しの始点の頂点を返す
        /// </summary>
        /// <param name="p_統合パス">統合した配列</param>
        /// <param name="p_橋渡しの場所">統合した配列の中の橋渡しした配列の場所</param>
        /// <param name="p_ライブラリ群">元リード</param>
        /// <returns>見送る始点の頂点 (評価できなければ空)</returns>
        private static HashSet<int> Get_危ない橋渡し(string p_統合パス, List<(string A_配列名, int A_開始, int A_終了, int A_始点)> p_橋渡しの場所, IReadOnlyList<(string A_リード1, string A_リード2)> p_ライブラリ群)
        {
            HashSet<int> l_見送る = [];
            if (p_橋渡しの場所.Count == 0 || JunctionRiskEvaluator.Get_評価(p_統合パス, p_ライブラリ群) is not { } l_評価群)
            {
                return l_見送る;
            }

            var l_場所 = p_橋渡しの場所.ToLookup(x => x.A_配列名.TrimStart('>'));
            var l_掛かる数 = 0;
            var l_最高 = 0D;
            foreach (var l_評価 in l_評価群)
            {
                foreach (var (A_配列名, A_開始, A_終了, A_始点) in l_場所[l_評価.A_候補.A_配列名])
                {
                    if (l_評価.A_候補.A_開始 >= A_終了 + C_橋渡しに掛かる余白 || l_評価.A_候補.A_終了 <= A_開始 - C_橋渡しに掛かる余白)
                    {
                        continue;
                    }

                    l_掛かる数++;
                    l_最高 = Math.Max(l_最高, l_評価.A_誤りの確率);
                    if (l_評価.A_誤りの確率 >= C_橋渡しを見送る確率)
                    {
                        _ = l_見送る.Add(A_始点);
                    }
                }
            }

            Logger.V_出力(メッセージID.橋渡しの継ぎ目を評価, p_橋渡しの場所.Count, l_評価群.Count, l_掛かる数, l_最高);
            return l_見送る;
        }

        /// <summary>
        /// 骨格に他の k の配列を統合し、良くなっていれば統合結果を返す
        /// </summary>
        /// <param name="p_最良"></param>
        /// <param name="p_候補"></param>
        /// <param name="p_アンカー"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_解析"></param>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_ライブラリ群">継ぎ目の誤り確率を見る元リード</param>
        /// <returns></returns>
        private static アセンブリ実行結果? Get_統合結果((アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価) p_最良, List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補, TrustedKmerIndex p_アンカー, int p_アンカーk長, スペクトル解析結果 p_解析, string p_一時ディレクトリ, IReadOnlyList<(string A_リード1, string A_リード2)> p_ライブラリ群)
        {
            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.統合開始);

            var l_統合パス = Path.Combine(p_一時ディレクトリ, AssemblyWorkspace.C_統合接頭辞 + Consts.Scaffoldファイル名);
            var l_全候補 = p_候補.Select(x => x.A_実行結果).ToList();
            var l_期待コピー数 = Get_期待コピー数(p_アンカー, p_アンカーk長, p_解析);
            List<(string A_配列名, int A_開始, int A_終了, int A_始点)> l_橋渡しの場所 = [];
            if (!AssemblyMerger.Try統合(p_最良.A_実行結果, l_全候補, p_アンカーk長, l_統合パス, l_期待コピー数, p_橋渡しの場所: l_橋渡しの場所))
            {
                return null;
            }

            var l_見送る始点 = Get_危ない橋渡し(l_統合パス, l_橋渡しの場所, p_ライブラリ群);
            if (l_見送る始点.Count > 0)
            {
                Logger.V_出力(メッセージID.危ない橋渡しを見送る, l_見送る始点.Count, C_橋渡しを見送る確率);
                if (!AssemblyMerger.Try統合(p_最良.A_実行結果, l_全候補, p_アンカーk長, l_統合パス, l_期待コピー数, p_見送る始点: l_見送る始点))
                {
                    return null;
                }
            }

            var l_統合contigパス = Path.Combine(p_一時ディレクトリ, AssemblyWorkspace.C_統合接頭辞 + AssemblyPipeline.C_Contigファイル名);
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
        /// 採用したアセンブリに無い配列を他の k から補い、良くなっていれば補った結果を返す
        /// </summary>
        /// <param name="p_採用">統合まで済んだ採用結果</param>
        /// <param name="p_候補">評価済みの全 k の候補</param>
        /// <param name="p_アンカー">アンカー k-mer 集合</param>
        /// <param name="p_アンカーk長">アンカー k 長</param>
        /// <param name="p_解析">アンカーのスペクトル解析結果</param>
        /// <param name="p_一時ディレクトリ">一時ディレクトリ</param>
        /// <returns>補った結果、補う配列が無いか良くならなければ null</returns>
        private static アセンブリ実行結果? Get_補完結果(アセンブリ実行結果 p_採用, List<(アセンブリ実行結果 A_実行結果, アセンブリ評価 A_評価)> p_候補, TrustedKmerIndex p_アンカー, int p_アンカーk長, スペクトル解析結果 p_解析, string p_一時ディレクトリ)
        {
            if (p_採用.A_固定アンカー評価 is not { } l_補完前の評価)
            {
                return null;
            }

            var l_補完パス = Path.Combine(p_一時ディレクトリ, AssemblyWorkspace.C_補完接頭辞 + Consts.Scaffoldファイル名);
            var (l_本数, l_延長) = AssemblyComplementer.Get_補完(p_採用, [.. p_候補.Select(x => x.A_実行結果)], x => p_アンカー.Haskmer_正規形(UInt128.Zero, x), p_アンカーk長, l_補完パス);
            if (l_本数 == 0)
            {
                return null;
            }

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.補完した配列, l_本数, l_延長);

            var l_補完contigパス = Path.Combine(p_一時ディレクトリ, AssemblyWorkspace.C_補完接頭辞 + AssemblyPipeline.C_Contigファイル名);
            V_書き出し_N分割(l_補完パス, l_補完contigパス);
            var l_補完結果 = p_採用 with
            {
                A_contigパス = l_補完contigパス,
                A_scaffoldパス = l_補完パス,
            };

            var l_補完後の評価 = AssemblyScorer.Get_評価(l_補完パス, p_アンカー, p_アンカーk長, p_解析.A_単一コピー基準値, p_解析.A_推定ゲノムサイズ, p_解析.A_単一コピー上限);
            if (l_補完後の評価 is null)
            {
                return null;
            }

            Logger.V_出力(メッセージID.補完後の評価, l_補完後の評価);
            var (l_実行結果, _) = AssemblySelector.Get_最良([(p_採用, l_補完前の評価), (l_補完結果, l_補完後の評価)])!.Value;
            if (l_実行結果.A_最終パス != l_補完パス)
            {
                Logger.V_出力(メッセージID.補完が勝てず);
                return null;
            }

            Logger.V_出力(メッセージID.補完結果を採用);
            return l_補完結果 with { A_固定アンカー評価 = l_補完後の評価 };
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
