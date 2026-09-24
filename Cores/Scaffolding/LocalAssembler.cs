using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;
using Tsumiki.Models.Evidence;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Reporting;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// GapFiller が埋められなかった scaffold のギャップを、局所アセンブリ (MEGAHIT/IDBA の localasm 型) で埋める
    /// </summary>
    internal static class LocalAssembler
    {
        #region 定数

        /// <summary>
        /// ギャップの両端からアンカーとして使う長さ
        /// </summary>
        private const int アンカー長 = 300;

        /// <summary>
        /// 1 ギャップに集める局所リードの上限
        /// </summary>
        /// <remarks>
        /// アンカーが反復配列と重なると際限なくリードが集まりうるため、暴走を防ぐ
        /// </remarks>
        private const int 局所リード数の上限 = 4_000;

        /// <summary>
        /// 並列に照合するペアの 1 まとまりの数
        /// </summary>
        private const int 照合のバッチサイズ = 16_384;

        /// <summary>
        /// 局所リードの reservoir sampling に使う乱数の種
        /// </summary>
        /// <remarks>
        /// 固定することで、同じ入力なら常に同じ選択結果になる (再現性のため)
        /// </remarks>
        private const int 局所リード選択の乱数種 = 20_260_914;

        /// <summary>
        /// 近傍thread拡張回収の reservoir sampling に使う乱数のシード値
        /// </summary>
        private const int 近傍拡張の乱数シード値 = 20_260_915;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 残ったギャップを、その両端に付いたリードだけで組み直して埋める
        /// </summary>
        /// <param name="p_scaffoldパス">対象の scaffold のパス</param>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>局所アセンブリの集計</returns>
        public static 局所アセンブリ統計 V_充填_ギャップ(string p_scaffoldパス, IReadOnlyList<(string A_リード1, string A_リード2)> p_ライブラリ群, int p_k長)
        {
            using var l_計測 = new StageTimer($"local-assembly k={p_k長}");
            var l_scaffold群 = FastaReader.Get_全エントリ(p_scaffoldパス);

            var l_ギャップ一覧 = Get_対象ギャップ一覧(l_scaffold群, p_k長);
            if (l_ギャップ一覧.Count == 0)
            {
                return new 局所アセンブリ統計(0, 0, 0, 0, 0, 0);
            }

            var l_アンカー索引 = Get_アンカー索引(l_ギャップ一覧, p_k長);
            var l_種長 = Math.Min(31, p_k長);
            HashSet<ulong> l_種集合 = [];
            foreach (var l_ギャップ in l_ギャップ一覧)
            {
                foreach (var l_アンカー in new[] { l_ギャップ.A_左アンカー, l_ギャップ.A_右アンカー })
                {
                    var l_窓 = new RollingKmer(l_種長);
                    foreach (var l_塩基 in l_アンカー)
                    {
                        if (l_窓.Try追加(l_塩基, out var l_キー))
                        {
                            _ = l_種集合.Add((ulong)l_キー.A_下位);
                        }
                    }
                }
            }
            var l_局所リード = Get_局所リード(l_アンカー索引, p_ライブラリ群, p_k長, l_ギャップ一覧.Count, l_種集合, l_種長);

            var l_結果 = new string?[l_ギャップ一覧.Count];
            var l_判定群 = new ギャップ充填判定?[l_ギャップ一覧.Count];

            // ギャップごとのミニアセンブリは互いに独立で、結果を番号の位置へ書くので並列にしても出力は変わらない
            var l_並列設定 = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数) };
            _ = Parallel.For(0, l_ギャップ一覧.Count, l_並列設定, g =>
            {
                if (l_局所リード[g].Count == 0)
                {
                    return;
                }
                l_結果[g] = Get_適応kの局所結果(l_ギャップ一覧[g], l_局所リード[g], p_k長, out var l_判定);
                l_判定群[g] = l_判定;
            });

            // 段階的な近傍thread回収: 1回目でリードは集まったが解決しなかったgapだけを対象に、
            // 既に回収した局所リード自体をシードにしたもう1パスで近傍のリードを拡張回収し、再試行する
            // (アンカー直接ヒットの先まで回収範囲を広げる。難所だけを対象にするため全体をやり直さない)
            var l_未解決gap = Enumerable.Range(0, l_ギャップ一覧.Count).Where(g => l_結果[g] is null && l_局所リード[g].Count > 0).ToList();
            if (l_未解決gap.Count > 0)
            {
                V_拡張回収_近傍thread(l_未解決gap, l_局所リード, p_ライブラリ群, p_k長);
                _ = Parallel.ForEach(l_未解決gap, l_並列設定, g =>
                {
                    l_結果[g] = Get_適応kの局所結果(l_ギャップ一覧[g], l_局所リード[g], p_k長, out var l_判定);
                    l_判定群[g] = l_判定;
                });
            }

            var l_埋めた数 = 0;
            var l_埋めた塩基数 = 0;
            var l_リード無し数 = 0;
            var l_一意でない数 = 0;
            var l_到達不能数 = 0;

            for (var g = 0; g < l_ギャップ一覧.Count; g++)
            {
                if (l_結果[g] is { } l_埋め)
                {
                    l_埋めた数++;
                    l_埋めた塩基数 += l_埋め.Length;
                    continue;
                }

                var l_ギャップ = l_ギャップ一覧[g];
                var l_場所 = $"scaffold{l_ギャップ.A_足場番号}:{l_ギャップ.A_開始}-{l_ギャップ.A_開始 + l_ギャップ.A_長さ}";
                var l_安定ID = AmbiguityRecorder.Get_安定ID(l_ギャップ.A_左アンカー, l_ギャップ.A_右アンカー);

                if (l_判定群[g] is not { } l_判定)
                {
                    l_リード無し数++;
                    AmbiguityRecorder.V_記録(曖昧箇所の種別.リード無し, l_場所, p_安定ID: l_安定ID);
                }
                else if (l_判定 == ギャップ充填判定.一意でない)
                {
                    l_一意でない数++;
                    AmbiguityRecorder.V_記録(曖昧箇所の種別.経路が一意でない, l_場所, p_安定ID: l_安定ID);
                }
                else
                {
                    l_到達不能数++;
                    AmbiguityRecorder.V_記録(l_判定 == ギャップ充填判定.探索打切り ? 曖昧箇所の種別.探索打切り : 曖昧箇所の種別.到達不能, l_場所, p_安定ID: l_安定ID);
                }
            }

            V_書き戻し(p_scaffoldパス, l_scaffold群, l_ギャップ一覧, l_結果);

            return new 局所アセンブリ統計(l_ギャップ一覧.Count, l_埋めた数, l_埋めた塩基数, l_リード無し数, l_一意でない数, l_到達不能数);
        }

        /// <summary>
        /// 局所アセンブリの結果をログへ出力する
        /// </summary>
        /// <param name="p_統計">局所アセンブリの結果</param>
        public static void V_出力_統計(局所アセンブリ統計 p_統計)
        {
            if (p_統計.A_対象ギャップ数 == 0)
            {
                Logger.V_出力(メッセージID.局所アセンブリ_対象なし);
                return;
            }
            Logger.V_出力(メッセージID.局所アセンブリ統計, p_統計.A_埋めたギャップ数, p_統計.A_対象ギャップ数, p_統計.A_埋めた塩基数, p_統計.A_局所リードが集まらなかった数, p_統計.A_一意に定まらなかった数, p_統計.A_到達できなかった数);
        }

        /// <summary>
        /// 指定した文脈長で局所グラフの一意な充填配列を求める
        /// </summary>
        /// <param name="p_ギャップ">充填対象</param>
        /// <param name="p_局所リード">局所に集めたリード</param>
        /// <param name="p_k長">文脈長</param>
        /// <param name="p_判定">探索の結果</param>
        /// <returns>一意な充填配列、未確定なら null</returns>
        internal static string? Get_固定kの局所結果(局所ギャップ p_ギャップ, List<読取証拠> p_局所リード, int p_k長, out ギャップ充填判定 p_判定)
        {
            if (p_ギャップ.A_左アンカー.Length < p_k長 || p_ギャップ.A_右アンカー.Length < p_k長)
            {
                p_判定 = ギャップ充填判定.到達不能;
                return null;
            }

            var l_集合 = new LocalKmerSet(p_k長);
            V_登録_全kmer(l_集合, p_ギャップ.A_左アンカー, p_k長);
            V_登録_全kmer(l_集合, p_ギャップ.A_右アンカー, p_k長);
            foreach (var l_証拠 in p_局所リード)
            {
                V_登録_全kmer(l_集合, l_証拠.A_配列, p_k長);
            }

            var l_左のkmer = Get_kmerバイト列(p_ギャップ.A_左アンカー, p_ギャップ.A_左アンカー.Length - p_k長, p_k長);
            var l_目標kmer = Get_kmerバイト列(p_ギャップ.A_右アンカー, 0, p_k長);
            if (l_左のkmer is null || l_目標kmer is null)
            {
                p_判定 = ギャップ充填判定.到達不能;
                return null;
            }

            var l_最小長 = Math.Max(0, p_ギャップ.A_長さ - Consts.ギャップ充填の長さの余裕幅);
            var l_最大長 = p_ギャップ.A_長さ + Consts.ギャップ充填の長さの余裕幅;
            (var l_経路, p_判定) = ConstrainedPathFinder.Get_経路(l_左のkmer, l_目標kmer, l_最小長, l_最大長, l_集合, p_k長);
            return l_経路;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 埋める対象になるギャップを集めて返す
        /// </summary>
        /// <param name="p_scaffold群">対象の scaffold</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>対象のギャップ</returns>
        private static List<局所ギャップ> Get_対象ギャップ一覧(List<(string A_ID, string A_配列)> p_scaffold群, int p_k長)
        {
            List<局所ギャップ> l_結果 = [];
            for (var s = 0; s < p_scaffold群.Count; s++)
            {
                var l_配列 = p_scaffold群[s].A_配列;
                var l_i = 0;
                while (l_i < l_配列.Length)
                {
                    if (l_配列[l_i] != 'N')
                    {
                        l_i++;
                        continue;
                    }
                    var l_開始 = l_i;
                    while (l_i < l_配列.Length && l_配列[l_i] == 'N')
                    {
                        l_i++;
                    }
                    var l_長さ = l_i - l_開始;

                    // 両端に最低 k 長ぶんの足場が要る (左右のアンカー k-mer を取るため)
                    if (l_長さ <= Consts.ギャップ充填のギャップ長上限 && l_開始 >= p_k長 && l_i + p_k長 <= l_配列.Length)
                    {
                        var l_左長 = Math.Min(アンカー長, l_開始);
                        var l_右長 = Math.Min(アンカー長, l_配列.Length - l_i);
                        l_結果.Add(new 局所ギャップ(s, l_開始, l_長さ, l_配列.Substring(l_開始 - l_左長, l_左長), l_配列.Substring(l_i, l_右長)));
                    }
                }
            }
            return l_結果;
        }

        /// <summary>
        /// 全ギャップの左右アンカーから k-mer 索引を作る
        /// </summary>
        /// <remarks>
        /// キーは正規形 (順鎖・逆鎖どちらでも同じキーに寄る)
        /// </remarks>
        /// <param name="p_ギャップ一覧"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static Dictionary<KmerKey, List<int>> Get_アンカー索引(List<局所ギャップ> p_ギャップ一覧, int p_k長)
        {
            Dictionary<KmerKey, List<int>> l_索引 = [];
            for (var g = 0; g < p_ギャップ一覧.Count; g++)
            {
                V_登録_kmer列(l_索引, p_ギャップ一覧[g].A_左アンカー, p_k長, g);
                V_登録_kmer列(l_索引, p_ギャップ一覧[g].A_右アンカー, p_k長, g);
            }
            return l_索引;
        }

        /// <summary>
        /// 配列の k-mer を、どのギャップの端かと一緒に索引へ登録する
        /// </summary>
        /// <param name="p_索引">登録先の索引</param>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_ギャップ番号">ギャップの番号</param>
        private static void V_登録_kmer列(Dictionary<KmerKey, List<int>> p_索引, string p_配列, int p_k長, int p_ギャップ番号)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (Has曖昧塩基(p_配列, i, p_k長))
                {
                    continue;
                }

                var l_鍵 = new KmerKey(p_配列.AsSpan(i, p_k長)).Get_正規形();
                if (!p_索引.TryGetValue(l_鍵, out var l_一覧))
                {
                    l_一覧 = [];
                    p_索引[l_鍵] = l_一覧;
                }

                if (!l_一覧.Contains(p_ギャップ番号))
                {
                    l_一覧.Add(p_ギャップ番号);
                }
            }
        }

        /// <summary>
        /// 指定した範囲に曖昧塩基が含まれるか
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始">調べ始める位置</param>
        /// <param name="p_長さ">調べる長さ</param>
        /// <returns>含まれれば true</returns>
        private static bool Has曖昧塩基(string p_配列, int p_開始, int p_長さ)
        {
            for (var j = 0; j < p_長さ; j++)
            {
                if (Util.Is曖昧塩基(p_配列[p_開始 + j]))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// 生リードを 1 回走査し、アンカーに触れたリードをギャップごとに集める
        /// </summary>
        /// <remarks>
        /// アンカーに触れたペアは両方のリードを局所リード集合に入れる (相方が、まだ組み込まれていない領域を読んでいる可能性があるため)
        /// </remarks>
        /// <param name="p_アンカー索引"></param>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <param name="p_k長"></param>
        /// <param name="p_ギャップ数"></param>
        /// <param name="p_種集合">アンカー内の短い正準キー</param>
        /// <param name="p_種長">種の長さ</param>
        /// <returns></returns>
        private static List<読取証拠>[] Get_局所リード(Dictionary<KmerKey, List<int>> p_アンカー索引, IReadOnlyList<(string A_リード1, string A_リード2)> p_ライブラリ群, int p_k長, int p_ギャップ数, HashSet<ulong> p_種集合, int p_種長)
        {
            var l_局所リード = new List<読取証拠>[p_ギャップ数];
            for (var g = 0; g < p_ギャップ数; g++)
            {
                l_局所リード[g] = [];
            }

            // 上限到達後にどの候補を残すかを、先着順ではなく一様な確率で決める (reservoir sampling)
            // ための、ギャップごとの遭遇数
            var l_遭遇数 = new int[p_ギャップ数];
            var l_乱数 = new Random(局所リード選択の乱数種);

            foreach (var (l_リード1のパス, l_リード2のパス) in p_ライブラリ群)
            {

                if (!string.IsNullOrWhiteSpace(l_リード1のパス) && !string.IsNullOrWhiteSpace(l_リード2のパス)
                    && 中間データ置き場.Is存在(l_リード1のパス) && 中間データ置き場.Is存在(l_リード2のパス))
                {
                    foreach (var (A_ID1, A_配列1, l_一致1, A_ID2, A_配列2, l_一致2) in Get_照合済みペア列(l_リード1のパス, l_リード2のパス, p_アンカー索引, p_k長, p_種集合, p_種長))
                    {
                        var l_pairID1 = Util.Get_ペア共通ID(A_ID1);
                        var l_pairID2 = Util.Get_ペア共通ID(A_ID2);

                        // 対応 ID が崩れた FASTQ を mate として混ぜると、無関係な配列を局所グラフへ持ち込む。
                        // その場合は各リード自身が当たったギャップだけへ、pair 情報を持たない単独読み取りとして入れる。
                        if (l_pairID1 != l_pairID2)
                        {
                            V_追加_局所リード(l_局所リード, l_遭遇数, l_乱数, l_一致1, new 読取証拠(A_配列1, ""));
                            V_追加_局所リード(l_局所リード, l_遭遇数, l_乱数, l_一致2, new 読取証拠(A_配列2, ""));
                            continue;
                        }

                        var l_証拠1 = new 読取証拠(A_配列1, l_pairID1);
                        var l_証拠2 = new 読取証拠(A_配列2, l_pairID2);
                        var l_ペアの一致 = l_一致1.Concat(l_一致2).ToHashSet();
                        foreach (var l_g in l_ペアの一致)
                        {
                            l_遭遇数[l_g]++;

                            // ペアを途中で切らない。残り 1 枠なら、直接アンカーに当たった側だけを優先する。
                            var l_残り枠 = 局所リード数の上限 - l_局所リード[l_g].Count;
                            if (l_残り枠 >= 2)
                            {
                                l_局所リード[l_g].Add(l_証拠1);
                                l_局所リード[l_g].Add(l_証拠2);
                            }
                            else if (l_残り枠 == 1)
                            {
                                l_局所リード[l_g].Add(l_一致1.Contains(l_g) ? l_証拠1 : l_証拠2);
                            }
                            else
                            {
                                // 上限到達後は、これまで遭遇したペアの中から一様な確率で選ばれるよう
                                // 代表 1 本を reservoir sampling で入れ替える (先着順のバイアスを避ける)
                                var l_置き換え位置 = l_乱数.Next(l_遭遇数[l_g]);
                                if (l_置き換え位置 < 局所リード数の上限)
                                {
                                    l_局所リード[l_g][l_置き換え位置] = l_一致1.Contains(l_g) ? l_証拠1 : l_証拠2;
                                }
                            }
                        }
                    }
                }
                else
                {
                    foreach (var l_パス in new[] { l_リード1のパス, l_リード2のパス })
                    {
                        if (!中間データ置き場.Is存在(l_パス))
                        {
                            continue;
                        }
                        foreach (var l_リード in FastqReader.Get_生リード列(l_パス))
                        {
                            V_追加_局所リード(l_局所リード, l_遭遇数, l_乱数, Get_一致するギャップ_候補選別付き(p_アンカー索引, l_リード, p_k長, p_種集合, p_種長), new 読取証拠(l_リード, ""));
                        }
                    }
                }
            }

            return l_局所リード;
        }

        /// <summary>
        /// 1回目でリードは集まったが解決しなかったgapについて、既に回収した局所リード自体をシードに、
        /// もう1パスだけ近傍のリードを追加回収する (anchor直接ヒットの先までの read thread 拡張)
        /// </summary>
        /// <param name="p_未解決gap">対象の gap 番号一覧 (p_局所リード のインデックス)</param>
        /// <param name="p_局所リード">gap ごとの局所リード (該当 gap 分を直接更新する)</param>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <param name="p_k長"></param>
        /// <remarks>
        /// 対象を「1回目で失敗した gap」だけに絞ることで、解決済みの gap に無駄な追加コストをかけない<br/>
        /// シードがアンカーではなく回収済みリードそのものになる点だけが元の回収と異なり、
        /// 一致判定・上限・reservoir sampling の仕組みは共通のヘルパーをそのまま再利用する
        /// </remarks>
        private static void V_拡張回収_近傍thread(List<int> p_未解決gap, List<読取証拠>[] p_局所リード, IReadOnlyList<(string A_リード1, string A_リード2)> p_ライブラリ群, int p_k長)
        {
            Dictionary<KmerKey, List<int>> l_拡張索引 = [];
            var l_種長 = Math.Min(31, p_k長);
            HashSet<ulong> l_種集合 = [];
            for (var i = 0; i < p_未解決gap.Count; i++)
            {
                foreach (var l_証拠 in p_局所リード[p_未解決gap[i]])
                {
                    V_登録_kmer列(l_拡張索引, l_証拠.A_配列, p_k長, i);

                    var l_窓 = new RollingKmer(l_種長);
                    foreach (var l_塩基 in l_証拠.A_配列)
                    {
                        if (l_窓.Try追加(l_塩基, out var l_キー))
                        {
                            _ = l_種集合.Add((ulong)l_キー.A_下位);
                        }
                    }
                }
            }
            if (l_拡張索引.Count == 0)
            {
                return;
            }

            var l_拡張プール = new List<読取証拠>[p_未解決gap.Count];
            for (var i = 0; i < p_未解決gap.Count; i++)
            {
                l_拡張プール[i] = [];
            }
            var l_遭遇数 = new int[p_未解決gap.Count];
            var l_乱数 = new Random(近傍拡張の乱数シード値);

            foreach (var (l_リード1のパス, l_リード2のパス) in p_ライブラリ群)
            {

                if (!string.IsNullOrWhiteSpace(l_リード1のパス) && !string.IsNullOrWhiteSpace(l_リード2のパス)
                    && 中間データ置き場.Is存在(l_リード1のパス) && 中間データ置き場.Is存在(l_リード2のパス))
                {
                    foreach (var (A_ID1, A_配列1, l_一致1, A_ID2, A_配列2, l_一致2) in Get_照合済みペア列(l_リード1のパス, l_リード2のパス, l_拡張索引, p_k長, l_種集合, l_種長))
                    {
                        var l_pairID1 = Util.Get_ペア共通ID(A_ID1);
                        var l_pairID2 = Util.Get_ペア共通ID(A_ID2);
                        var l_pairID = l_pairID1 == l_pairID2 ? l_pairID1 : "";

                        V_追加_局所リード(l_拡張プール, l_遭遇数, l_乱数, l_一致1, new 読取証拠(A_配列1, l_pairID));
                        V_追加_局所リード(l_拡張プール, l_遭遇数, l_乱数, l_一致2, new 読取証拠(A_配列2, l_pairID));
                    }
                }
                else
                {
                    foreach (var l_パス in new[] { l_リード1のパス, l_リード2のパス })
                    {
                        if (!中間データ置き場.Is存在(l_パス))
                        {
                            continue;
                        }
                        foreach (var l_リード in FastqReader.Get_生リード列(l_パス))
                        {
                            V_追加_局所リード(l_拡張プール, l_遭遇数, l_乱数, Get_一致するギャップ_候補選別付き(l_拡張索引, l_リード, p_k長, l_種集合, l_種長), new 読取証拠(l_リード, ""));
                        }
                    }
                }
            }

            for (var i = 0; i < p_未解決gap.Count; i++)
            {
                if (l_拡張プール[i].Count == 0)
                {
                    continue;
                }
                var l_g = p_未解決gap[i];
                var l_既存配列 = p_局所リード[l_g].Select(x => x.A_配列).ToHashSet(StringComparer.Ordinal);
                foreach (var l_証拠 in l_拡張プール[i])
                {
                    if (p_局所リード[l_g].Count >= 局所リード数の上限)
                    {
                        break;
                    }
                    if (l_既存配列.Add(l_証拠.A_配列))
                    {
                        p_局所リード[l_g].Add(l_証拠);
                    }
                }
            }
        }

        /// <summary>
        /// ペアのリードを読み進め、それぞれが当たるギャップを並列に照合して読み込み順に返す
        /// </summary>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <param name="p_索引"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_種集合"></param>
        /// <param name="p_種長"></param>
        /// <remarks>
        /// 照合は読み取りだけで互いに独立だが、局所リードへの追加は乱数による入れ替えを含み順序に依存するため、照合だけをまとめて並列にし、追加は呼び出し側が順に行う<br/>
        /// 読み込みは 1 本の流れでしか進められないので、あるバッチを照合している間に次のバッチを読み込んでおく
        /// </remarks>
        /// <returns></returns>
        private static IEnumerable<(string A_ID1, string A_配列1, HashSet<int> A_一致1, string A_ID2, string A_配列2, HashSet<int> A_一致2)> Get_照合済みペア列(string p_リード1のパス, string p_リード2のパス, Dictionary<KmerKey, List<int>> p_索引, int p_k長, HashSet<ulong> p_種集合, int p_種長)
        {
            var l_並列設定 = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数) };

            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);
            var l_次の読み込み = Task.Run(() => Get_ペアのバッチ(l_読み込み1, l_読み込み2));
            try
            {
                while (true)
                {
                    var l_バッチ = l_次の読み込み.Result;
                    if (l_バッチ.Count == 0)
                    {
                        yield break;
                    }
                    l_次の読み込み = Task.Run(() => Get_ペアのバッチ(l_読み込み1, l_読み込み2));

                    _ = Parallel.For(0, l_バッチ.Count, l_並列設定, i =>
                    {
                        var l_項目 = l_バッチ[i];
                        l_バッチ[i] = l_項目 with
                        {
                            A_一致1 = Get_一致するギャップ_候補選別付き(p_索引, l_項目.A_配列1, p_k長, p_種集合, p_種長),
                            A_一致2 = Get_一致するギャップ_候補選別付き(p_索引, l_項目.A_配列2, p_k長, p_種集合, p_種長),
                        };
                    });

                    foreach (var l_項目 in l_バッチ)
                    {
                        yield return l_項目;
                    }
                }
            }
            finally
            {
                // 途中で列挙をやめられても、読み込み中のリーダーを閉じない
                try
                {
                    l_次の読み込み.Wait();
                }
                catch (AggregateException)
                {
                }
            }
        }

        /// <summary>
        /// read1/read2 を照合のバッチサイズぶん読み込む
        /// </summary>
        /// <param name="p_読み込み1"></param>
        /// <param name="p_読み込み2"></param>
        /// <returns>読み込んだペア、読み終えていれば空</returns>
        private static List<(string A_ID1, string A_配列1, HashSet<int> A_一致1, string A_ID2, string A_配列2, HashSet<int> A_一致2)> Get_ペアのバッチ(FastqReader p_読み込み1, FastqReader p_読み込み2)
        {
            List<(string A_ID1, string A_配列1, HashSet<int> A_一致1, string A_ID2, string A_配列2, HashSet<int> A_一致2)> l_バッチ = new(照合のバッチサイズ);
            while (l_バッチ.Count < 照合のバッチサイズ && p_読み込み1.Has続き() && p_読み込み2.Has続き())
            {
                var (A_ID1, A_配列1, _) = p_読み込み1.Get_次のレコード();
                var (A_ID2, A_配列2, _) = p_読み込み2.Get_次のレコード();
                l_バッチ.Add((A_ID1, A_配列1, [], A_ID2, A_配列2, []));
            }
            return l_バッチ;
        }

        private static HashSet<int> Get_一致するギャップ_候補選別付き(Dictionary<KmerKey, List<int>> p_アンカー索引, string p_リード, int p_k長, HashSet<ulong> p_種集合, int p_種長)
        {
            return p_リード.Length >= p_k長 && Has種一致(p_リード, p_種集合, p_種長)
                ? Get_一致するギャップ(p_アンカー索引, p_リード, p_k長)
                : [];
        }

        /// <summary>
        /// 1 本ぶんの証拠をギャップの局所リードへ追加する
        /// </summary>
        /// <remarks>
        /// 上限に達した後も遭遇数を数え続け、reservoir sampling で一様な確率での入れ替えに使う<br/>
        /// これにより、反復配列でリードが際限なく集まる状況でも、常に先着した本数だけが残る先着順バイアスを避ける
        /// </remarks>
        /// <param name="p_局所リード"></param>
        /// <param name="p_遭遇数">ギャップごとの遭遇数 (reservoir sampling 用)</param>
        /// <param name="p_乱数">選択に使う乱数 (固定シードで決定的にする)</param>
        /// <param name="p_ギャップ群"></param>
        /// <param name="p_証拠"></param>
        private static void V_追加_局所リード(List<読取証拠>[] p_局所リード, int[] p_遭遇数, Random p_乱数, IEnumerable<int> p_ギャップ群, 読取証拠 p_証拠)
        {
            foreach (var l_g in p_ギャップ群)
            {
                p_遭遇数[l_g]++;
                if (p_局所リード[l_g].Count < 局所リード数の上限)
                {
                    p_局所リード[l_g].Add(p_証拠);
                    continue;
                }

                var l_置き換え位置 = p_乱数.Next(p_遭遇数[l_g]);
                if (l_置き換え位置 < 局所リード数の上限)
                {
                    p_局所リード[l_g][l_置き換え位置] = p_証拠;
                }
            }
        }

        /// <summary>リードがアンカーと短い完全一致キーを共有するか調べる</summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_種集合">アンカーの正準キー集合</param>
        /// <param name="p_種長">キーの長さ</param>
        /// <returns>候補になれば true</returns>
        internal static bool Has種一致(string p_リード, HashSet<ulong> p_種集合, int p_種長)
        {
            var l_窓 = new RollingKmer(p_種長);
            foreach (var l_塩基 in p_リード)
            {
                if (l_窓.Try追加(l_塩基, out var l_キー) && p_種集合.Contains((ulong)l_キー.A_下位))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// そのリードが端に当たるギャップの番号を返す
        /// </summary>
        /// <param name="p_索引">ギャップの端の k-mer 索引</param>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>当たったギャップの番号</returns>
        private static HashSet<int> Get_一致するギャップ(Dictionary<KmerKey, List<int>> p_索引, string p_リード, int p_k長)
        {
            HashSet<int>? l_見つかった = null;
            for (var i = 0; i + p_k長 <= p_リード.Length; i++)
            {
                if (Has曖昧塩基(p_リード, i, p_k長))
                {
                    continue;
                }

                var l_鍵 = new KmerKey(p_リード.AsSpan(i, p_k長)).Get_正規形();
                if (p_索引.TryGetValue(l_鍵, out var l_一覧))
                {
                    l_見つかった ??= [];
                    foreach (var l_g in l_一覧)
                    {
                        _ = l_見つかった.Add(l_g);
                    }
                }
            }
            return l_見つかった ?? [];
        }

        /// <summary>
        /// 1 ギャップぶんのミニアセンブリ
        /// </summary>
        /// <remarks>
        /// 左右アンカー配列+局所リードだけから使い捨ての LocalKmerSet を作り、GapFiller と同じ制約付き探索で左アンカー末尾から右アンカー先頭までの経路を探す<br/>
        /// 集めた k-mer は局所リード数の上限ぶんしかなくインメモリで完結するため、TrustedKmerIndex のようなディスク経由のシャード集計は使わない (ギャップの数だけ繰り返すには重すぎる)
        /// </remarks>
        /// <param name="p_ギャップ"></param>
        /// <param name="p_局所リード"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_判定"></param>
        /// <returns></returns>
        internal static string? Get_適応kの局所結果(局所ギャップ p_ギャップ, List<読取証拠> p_局所リード, int p_k長, out ギャップ充填判定 p_判定)
        {
            var l_異なるリード = p_局所リード.DistinctBy(x => x.A_配列, StringComparer.Ordinal).ToList();

            // base-k の結果も、fallback 候補と全く同じ「複数リードの支持」ゲートを通す
            // これを素通りさせると、橋渡し全体をまたぐリードが 1 本しかなくても
            // (あるいはアンカー由来の k-mer だけでも) 採用されてしまう
            var l_経路 = Get_固定kの局所結果(p_ギャップ, p_局所リード, p_k長, out p_判定);
            if (l_経路 is not null && Has経路支持(p_ギャップ, l_経路, l_異なるリード, p_k長))
            {
                return l_経路;
            }

            string? l_一致した経路 = null;
            var l_候補k = p_判定 == ギャップ充填判定.到達不能
                ? Get_低k候補(p_k長)
                : [p_k長 + 10, p_k長 + 20];
            var l_競合あり = false;
            foreach (var l_局所k in l_候補k)
            {
                if (l_異なるリード.Count(x => x.A_配列.Length >= l_局所k + 1) < 2)
                {
                    continue;
                }

                var l_候補 = Get_固定kの局所結果(p_ギャップ, l_異なるリード, l_局所k, out var l_局所判定);
                if (l_局所判定 == ギャップ充填判定.一意でない)
                {
                    l_競合あり = true;
                    continue;
                }

                if (l_候補 is null)
                {
                    continue;
                }

                // 支持の無い候補は無視するだけにとどめ、既に支持された結論を
                // 曖昧扱いにしない (support の有無を先に確かめてから不一致を見る)
                if (!Has経路支持(p_ギャップ, l_候補, l_異なるリード, l_局所k))
                {
                    continue;
                }

                if (l_一致した経路 is not null && l_一致した経路 != l_候補)
                {
                    p_判定 = ギャップ充填判定.一意でない;
                    return null;
                }
                l_一致した経路 = l_候補;
            }
            if (l_一致した経路 is not null)
            {
                p_判定 = ギャップ充填判定.充填済み;
                return l_一致した経路;
            }
            p_判定 = l_競合あり ? ギャップ充填判定.一意でない : ギャップ充填判定.到達不能;
            return null;
        }

        /// <summary>
        /// 高 k で非連結だった局所グラフを救済する低 k 候補
        /// </summary>
        internal static IReadOnlyList<int> Get_低k候補(int p_k長)
        {
            return [.. new[] { p_k長 - 10, p_k長 - 16, 27, 21 }
                .Where(x => x >= 15 && x < p_k長)
                .Distinct()
                .OrderByDescending(x => x)];
        }

        /// <summary>
        /// 候補経路の各辺が、アンカー合成ではなく複数の独立した分子で観測されていることを確かめる
        /// </summary>
        /// <remarks>
        /// 支持は読み取り本数ではなく独立性キー (pair ID があればそれ、無ければ配列自身) の種類数で数える<br/>
        /// mate1/mate2 は同一分子由来なので、2 本読めても 1 件としてしか数えない
        /// </remarks>
        private static bool Has経路支持(局所ギャップ p_ギャップ, string p_経路, IReadOnlyList<読取証拠> p_リード群, int p_k長)
        {
            var l_窓長 = p_k長 + 1;
            var l_接続 = p_ギャップ.A_左アンカー[^p_k長..] + p_経路 + p_ギャップ.A_右アンカー[..p_k長];

            // 窓ごとに全リードを部分文字列検索すると、経路長 × リード数 × リード長になる
            // 窓とその逆相補を先に集め、リードの窓を 1 回ずつ引いて支持した独立性キーを貯める
            Dictionary<string, HashSet<string>> l_窓別支持 = new(StringComparer.Ordinal);
            List<(string A_窓, string A_逆窓)> l_窓一覧 = [];
            for (var i = 0; i + l_窓長 <= l_接続.Length; i++)
            {
                var l_窓 = l_接続.Substring(i, l_窓長);
                var l_逆窓 = Util.V_逆相補(l_窓);
                l_窓一覧.Add((l_窓, l_逆窓));
                _ = l_窓別支持.TryAdd(l_窓, new HashSet<string>(StringComparer.Ordinal));
                _ = l_窓別支持.TryAdd(l_逆窓, new HashSet<string>(StringComparer.Ordinal));
            }

            var l_参照 = l_窓別支持.GetAlternateLookup<ReadOnlySpan<char>>();
            foreach (var l_リード in p_リード群)
            {
                for (var i = 0; i + l_窓長 <= l_リード.A_配列.Length; i++)
                {
                    if (l_参照.TryGetValue(l_リード.A_配列.AsSpan(i, l_窓長), out var l_支持))
                    {
                        _ = l_支持.Add(l_リード.A_独立性キー);
                    }
                }
            }

            foreach (var (l_窓, l_逆窓) in l_窓一覧)
            {
                var l_支持 = l_窓別支持[l_窓];
                var l_独立支持数 = ReferenceEquals(l_支持, l_窓別支持[l_逆窓]) ? l_支持.Count : l_支持.Union(l_窓別支持[l_逆窓], StringComparer.Ordinal).Count();
                if (l_独立支持数 < 2)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// 配列のすべての k-mer を索引へ登録する
        /// </summary>
        /// <param name="p_索引">登録先の索引</param>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_k長">k 長</param>
        private static void V_登録_全kmer(LocalKmerSet p_索引, string p_配列, int p_k長)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (Has曖昧塩基(p_配列, i, p_k長))
                {
                    continue;
                }

                var l_kmer = Get_kmerバイト列(p_配列, i, p_k長);
                if (l_kmer is not null)
                {
                    p_索引.V_登録(l_kmer);
                }
            }
        }

        /// <summary>
        /// 指定した位置の k-mer を塩基 ID 列にして返す
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始">k-mer の開始位置</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>塩基 ID 列、曖昧塩基を含む場合は null</returns>
        private static byte[]? Get_kmerバイト列(string p_配列, int p_開始, int p_k長)
        {
            if (p_開始 < 0 || p_開始 + p_k長 > p_配列.Length)
            {
                return null;
            }
            var l_kmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[p_開始 + i]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    return null;
                }
                l_kmer[i] = l_塩基ID;
            }
            return l_kmer;
        }

        /// <summary>
        /// 埋まったギャップを反映して scaffold を書き直す
        /// </summary>
        /// <param name="p_scaffoldパス">書き出し先</param>
        /// <param name="p_scaffold群">対象の scaffold</param>
        /// <param name="p_ギャップ一覧">埋めたギャップ</param>
        /// <param name="p_結果"></param>
        private static void V_書き戻し(string p_scaffoldパス, List<(string A_ID, string A_配列)> p_scaffold群, List<局所ギャップ> p_ギャップ一覧, string?[] p_結果)
        {
            var l_scaffold別ギャップ = p_ギャップ一覧
                .Select((l_ギャップ, l_番号) => (l_ギャップ, l_番号))
                .GroupBy(x => x.l_ギャップ.A_足場番号)
                .ToDictionary(x => x.Key, x => x.ToList());

            using var l_書き込み = new FastaWriter(p_scaffoldパス);
            for (var s = 0; s < p_scaffold群.Count; s++)
            {
                var (l_ID, l_配列) = p_scaffold群[s];
                if (!l_scaffold別ギャップ.TryGetValue(s, out var l_該当))
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                    continue;
                }

                var l_出力 = new StringBuilder();
                var l_直前終端 = 0;
                foreach (var (l_ギャップ, l_番号) in l_該当)
                {
                    _ = l_出力.Append(l_配列, l_直前終端, l_ギャップ.A_開始 - l_直前終端);
                    var l_埋め = p_結果[l_番号];
                    _ = l_埋め != null ? l_出力.Append(l_埋め) : l_出力.Append('N', l_ギャップ.A_長さ);
                    l_直前終端 = l_ギャップ.A_開始 + l_ギャップ.A_長さ;
                }
                _ = l_出力.Append(l_配列, l_直前終端, l_配列.Length - l_直前終端);
                l_書き込み.V_書き込み(l_ID, l_出力.ToString());
            }
        }

        #endregion
    }
}
