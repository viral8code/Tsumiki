using System.Text;
using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Evidence;
using Tsumiki.IO;
using Tsumiki.Models.ContigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Reporting;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 確定した contig を読み直し、ペアエンド由来の隣接で N 埋め連結する
    /// </summary>
    /// <param name="p_contig構築">確定辺・配置情報を持つ contig 構築器</param>
    /// <param name="p_contigファイルパス">読み直す contig ファイルのパス</param>
    /// <param name="p_リード長">リード長、不明なら null</param>
    internal class Scaffolder(ContigMaker p_contig構築, string p_contigファイルパス, int? p_リード長)
    {
        #region 定数

        /// <summary>
        /// 同一 unitig 内標本を信頼してよい「unitig 長 / 推定フラグメント長」の下限比
        /// </summary>
        private const int 偏りが無いとみなす長さ比 = 10;

        /// <summary>
        /// scaffold 辺に必要なペア数
        /// </summary>
        private const ulong Scaffold支持数の下限 = 3UL;

        /// <summary>
        /// インサートサイズ推定に必要な標本数
        /// </summary>
        private const int インサートサイズ標本数の下限 = 30;

        /// <summary>
        /// 採用した辺をライブラリ別の支持と共に書き出すファイル名 (scaffold と同じ場所に置く)
        /// </summary>
        private const string 採用辺の書き出し名 = "scaffold_edges.tsv";

        /// <summary>
        /// 支持の下限を満たす候補が 2 本以上ある頂点の候補を書き出すファイル名 (scaffold と同じ場所に置く)
        /// </summary>
        private const string 競合候補の書き出し名 = "scaffold_candidates.tsv";

        /// <summary>
        /// scaffold に挿入する最小ギャップ長
        /// </summary>
        private const int ギャップ長の下限 = 1;

        /// <summary>
        /// k-1 に満たない重なりを畳むときに求める最短の一致長
        /// </summary>
        private const int 短い重なりの下限 = 15;

        #endregion

        #region 内部変数

        /// <summary>
        /// contig 配列
        /// </summary>
        private readonly Dictionary<int, string> _contig配列 = [];

        /// <summary>
        /// contig 名
        /// </summary>
        private readonly Dictionary<int, string> _contig名 = [];

        #endregion

        #region プロパティ

        /// <summary>
        /// 自動推定された (あるいは CLI で明示指定された) インサートサイズ
        /// </summary>
        public int? A_有効インサートサイズ { get; private set; }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// scaffolding を実行し、指定パスに結果を書き出す
        /// </summary>
        /// <param name="p_scaffoldパス"></param>
        public void V_実行(string p_scaffoldパス)
        {
            var l_ライブラリ数 = Math.Max(1, p_contig構築.A_ペアのライブラリ数);
            var l_インサートサイズ群 = new int[l_ライブラリ数];
            var l_使えるライブラリ = 0;
            for (var l_ライブラリ = 0; l_ライブラリ < l_ライブラリ数; l_ライブラリ++)
            {
                if (this.TryGet_インサートサイズ(l_ライブラリ, out var l_値))
                {
                    l_インサートサイズ群[l_ライブラリ] = l_値;
                    l_使えるライブラリ++;
                    this.A_有効インサートサイズ ??= l_値;
                }
            }

            if (l_使えるライブラリ == 0)
            {
                Logger.V_出力(メッセージID.Scaffolding省略_インサートサイズ不明);
                return;
            }
            Logger.V_出力(メッセージID.Scaffolding開始_インサートサイズ, this.A_有効インサートサイズ!.Value);

            this.V_読込_Contig();

            if (this._contig配列.Count == 0)
            {
                Logger.V_出力(メッセージID.Scaffolding省略_Contigなし);
                return;
            }

            var l_配置 = p_contig構築.A_unitig配置;

            var l_contig数 = this._contig配列.Keys.Count == 0 ? 0 : this._contig配列.Keys.Max();
            var l_頂点数 = (l_contig数 + 1) << 1;

            var l_隣接 = new List<Scaffold候補>[l_頂点数];
            for (var i = 0; i < l_頂点数; i++)
            {
                l_隣接[i] = [];
            }

            var l_対称化群 = new Dictionary<(int, int), (ulong A_支持数, List<int> A_既知長標本)>[l_ライブラリ数];
            var l_モデル群 = new PairedDistanceModel[l_ライブラリ数];
            var l_較正器群 = new 証拠較正器[l_ライブラリ数];
            var l_内部を指した数 = 0;
            var l_未配置を指した数 = 0;
            var l_unitig長 = p_contig構築.A_unitig長.Values.Select(x => (long)x).ToList();

            for (var l_ライブラリ = 0; l_ライブラリ < l_ライブラリ数; l_ライブラリ++)
            {
                var l_ペア経路 = l_ライブラリ < p_contig構築.A_ペアのライブラリ数
                    ? p_contig構築.A_ペア経路群[l_ライブラリ]
                    : new Dictionary<(int, int), List<int>>();
                l_対称化群[l_ライブラリ] = this.Get_対称化した辺(l_配置, l_ペア経路, ref l_内部を指した数, ref l_未配置を指した数);

                var l_標本 = l_ライブラリ < p_contig構築.A_同一unitig標本群.Count
                    ? p_contig構築.A_同一unitig標本群[l_ライブラリ]
                    : [];
                var l_リード長群 = ConfigurationManager.A_実行時引数.A_ライブラリのリード長;
                var l_この長さ = l_ライブラリ < l_リード長群.Count && l_リード長群[l_ライブラリ] > 0
                    ? l_リード長群[l_ライブラリ]
                    : p_リード長 ?? (l_インサートサイズ群[l_ライブラリ] > 0 ? l_インサートサイズ群[l_ライブラリ] : this.A_有効インサートサイズ!.Value);
                l_モデル群[l_ライブラリ] = new PairedDistanceModel(l_標本, l_この長さ);
                l_較正器群[l_ライブラリ] = 証拠較正器.Get_較正器(l_標本, l_この長さ, l_unitig長);

                if (l_ライブラリ数 > 1)
                {
                    Logger.V_出力_そのまま(FormattableString.Invariant(
                        $"[Info] ライブラリ {l_ライブラリ + 1}: インサートサイズ {l_インサートサイズ群[l_ライブラリ]:N0}、リード長 {l_この長さ:N0}、断片長標本 {l_標本.Count:N0} 件、ペア辺 {l_対称化群[l_ライブラリ].Count:N0} 本"));
                }
            }

            if (l_内部を指した数 > 0)
            {
                Logger.V_出力(メッセージID.内部を指したペア候補, l_内部を指した数);
            }

            if (l_未配置を指した数 > 0)
            {
                Logger.V_出力(メッセージID.未配置を指したペア候補, l_未配置を指した数);
            }

            var l_候補キー = l_対称化群.SelectMany(x => x.Keys).ToHashSet();
            var l_ライブラリ別本数 = l_ライブラリ数 > 1 ? new Dictionary<(int, int), (int[] A_本数, double[] A_期待)>() : null;

            foreach (var (l_始点, l_終点) in l_候補キー)
            {
                var l_本数群 = l_ライブラリ別本数 is null ? null : new int[l_ライブラリ数];
                var l_最良本数 = 0;
                var l_最良ギャップ = 0;
                var l_最良比 = 0D;
                for (var l_ライブラリ = 0; l_ライブラリ < l_ライブラリ数; l_ライブラリ++)
                {
                    if (!l_対称化群[l_ライブラリ].TryGetValue((l_始点, l_終点), out var l_項目))
                    {
                        continue;
                    }

                    var (l_一貫した本数, l_ギャップ長) = l_モデル群[l_ライブラリ].Get_一貫した支持(l_項目.A_既知長標本);
                    l_本数群?[l_ライブラリ] = l_一貫した本数;
                    if (l_一貫した本数 <= l_最良本数)
                    {
                        continue;
                    }

                    l_最良本数 = l_一貫した本数;
                    l_最良ギャップ = l_ギャップ長;
                    l_最良比 = l_較正器群[l_ライブラリ].Get_正規化済み支持((ulong)l_一貫した本数, this.Get_Contig長(l_始点), this.Get_Contig長(l_終点), Math.Max(0, l_ギャップ長));
                }

                l_隣接[l_始点].Add(new Scaffold候補(l_終点, (ulong)l_最良本数, l_最良ギャップ, l_最良比));
                if (l_本数群 is not null)
                {
                    var l_期待群 = new double[l_ライブラリ数];
                    for (var l_ライブラリ = 0; l_ライブラリ < l_ライブラリ数; l_ライブラリ++)
                    {
                        l_期待群[l_ライブラリ] = l_較正器群[l_ライブラリ].Get_期待本数(this.Get_Contig長(l_始点), this.Get_Contig長(l_終点), Math.Max(0, l_最良ギャップ));
                    }
                    l_ライブラリ別本数![(l_始点, l_終点)] = (l_本数群, l_期待群);
                }
            }

            var l_優勢閾値 = ConfigurationManager.A_実行時引数.A_ペア結合閾値;
            var l_最小証拠数 = Scaffold支持数の下限;

            Logger.V_出力(メッセージID.Scaffold候補辺数, l_候補キー.Count, Messages.Get_文言(l_較正器群.Any(x => x.A_Is使用可能) ? メッセージID.理想本数モデルあり : メッセージID.理想本数モデルなし));

            var l_確定辺 = new (int A_行き先, int A_ギャップ長)?[l_頂点数];
            for (var v = 2; v < l_頂点数; v++)
            {
                V_確定_Scaffold辺(l_隣接, v, l_優勢閾値, l_最小証拠数, l_確定辺);
            }

            if (l_ライブラリ別本数 is not null)
            {
                this.V_書き出し_競合候補(Path.Combine(Path.GetDirectoryName(Path.GetFullPath(p_scaffoldパス))!, 競合候補の書き出し名), l_隣接, l_確定辺, l_ライブラリ別本数, l_最小証拠数);
            }

            var l_確定数 = 0;
            for (var v = 2; v < l_頂点数; v++)
            {
                if (l_確定辺[v] != null)
                {
                    l_確定数++;
                }
            }

            var l_候補辺 = ((int A_行き先, int A_ギャップ長)?[])l_確定辺.Clone();
            var l_相互一意で棄却した数 = 0;
            for (var v = 2; v < l_頂点数; v++)
            {
                if (l_候補辺[v] is not { } l_辺)
                {
                    continue;
                }

                var l_双子 = l_辺.A_行き先 ^ 1;
                if (l_双子 >= l_頂点数 || l_候補辺[l_双子] is not { } l_戻りの辺 || l_戻りの辺.A_行き先 != (v ^ 1))
                {
                    l_確定辺[v] = null;
                    l_相互一意で棄却した数++;
                    AmbiguityRecorder.V_記録(曖昧箇所の種別.経路が一意でない, AmbiguityRecorder.Get_場所名(v, "contig"));
                }
            }
            Logger.V_出力(メッセージID.閾値後のscaffold辺, l_確定数, l_相互一意で棄却した数, l_確定数 - l_相互一意で棄却した数);

            if (l_ライブラリ別本数 is not null)
            {
                this.V_書き出し_採用辺(Path.Combine(Path.GetDirectoryName(Path.GetFullPath(p_scaffoldパス))!, 採用辺の書き出し名), l_確定辺, l_ライブラリ別本数);
            }

            var l_始点群 = new List<int>();
            for (var v = 2; v < l_頂点数; v++)
            {
                if (this._contig配列.ContainsKey(v >> 1) && (v ^ 1) < l_頂点数 && l_確定辺[v ^ 1] == null)
                {
                    l_始点群.Add(v);
                }
            }

            List<(string A_配列, bool A_Is環状)> l_scaffold群 = [];
            var l_訪問済み = new bool[l_頂点数];
            foreach (var l_始点 in l_始点群)
            {
                if (l_訪問済み[l_始点])
                {
                    continue;
                }

                var l_scaffold = this.Get_Scaffold配列(l_確定辺, l_始点, l_訪問済み, out var l_連結数);
                if (l_scaffold != null)
                {
                    l_scaffold群.Add((l_scaffold, l_連結数 == 1 && this.Is環状(l_始点 >> 1)));
                }
            }

            for (var l_contigID = 1; l_contigID <= l_contig数; l_contigID++)
            {
                var l_順鎖 = l_contigID << 1;
                var l_逆鎖 = (l_contigID << 1) | 1;
                if (l_順鎖 < l_頂点数 && !l_訪問済み[l_順鎖] && !l_訪問済み[l_逆鎖] && this._contig配列.TryGetValue(l_contigID, out var l_配列))
                {
                    l_scaffold群.Add((l_配列, this.Is環状(l_contigID)));
                    l_訪問済み[l_順鎖] = true;
                    l_訪問済み[l_逆鎖] = true;
                }
            }

            using var l_書き込み = new FastaWriter(p_scaffoldパス);
            var l_scaffoldID = 1;
            var l_総延長 = 0L;
            foreach (var (l_配列, l_Is環状) in l_scaffold群)
            {
                var l_名前 = l_Is環状
                    ? $"SCAFFOLD{l_scaffoldID}_{Consts.環状の目印}"
                    : $"SCAFFOLD{l_scaffoldID}";
                l_書き込み.V_書き込み(l_名前, l_配列);
                l_scaffoldID++;
                l_総延長 += l_配列.Length;
            }

            Logger.V_出力(メッセージID.Scaffold出力完了, l_scaffold群.Count, l_総延長, p_scaffoldパス);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 ライブラリのペア経路を contig 末端の辺へ畳み、双子側も含めて対称化する
        /// </summary>
        /// <param name="p_配置">unitig の contig 上の配置</param>
        /// <param name="p_ペア経路">そのライブラリのペア経路</param>
        /// <param name="p_内部を指した数">contig 内部を指した観測の数</param>
        /// <param name="p_未配置を指した数">contig に載っていない unitig を指した観測の数</param>
        /// <returns>対称化した辺の集計</returns>
        private Dictionary<(int, int), (ulong A_支持数, List<int> A_既知長標本)> Get_対称化した辺(
            IReadOnlyDictionary<int, Unitig配置> p_配置,
            IReadOnlyDictionary<(int, int), List<int>> p_ペア経路,
            ref int p_内部を指した数,
            ref int p_未配置を指した数)
        {
            Dictionary<(int, int), (ulong A_支持数, List<int> A_既知長標本)> l_辺の集計 = [];

            foreach (var (l_キー, l_標本) in p_ペア経路)
            {
                var (l_始点unitig, l_終点unitig) = l_キー;

                if (!TryGet_Contig末端頂点(p_配置, l_始点unitig, p_Is出口側: true, out var l_始点頂点))
                {
                    if (!p_配置.ContainsKey(Math.Abs(l_始点unitig)))
                    {
                        p_未配置を指した数++;
                    }
                    else
                    {
                        p_内部を指した数++;
                    }
                    continue;
                }

                if (!TryGet_Contig末端頂点(p_配置, l_終点unitig, p_Is出口側: false, out var l_終点頂点))
                {
                    if (!p_配置.ContainsKey(Math.Abs(l_終点unitig)))
                    {
                        p_未配置を指した数++;
                    }
                    else
                    {
                        p_内部を指した数++;
                    }
                    continue;
                }

                if (l_始点頂点 >> 1 == l_終点頂点 >> 1)
                {
                    continue;
                }

                var l_辺キー = (l_始点頂点, l_終点頂点);
                if (l_辺の集計.TryGetValue(l_辺キー, out var l_既存))
                {
                    l_既存.A_既知長標本.AddRange(l_標本);
                    l_辺の集計[l_辺キー] = (l_既存.A_支持数 + (ulong)l_標本.Count, l_既存.A_既知長標本);
                }
                else
                {
                    l_辺の集計[l_辺キー] = ((ulong)l_標本.Count, [.. l_標本]);
                }
            }

            Dictionary<(int, int), (ulong A_支持数, List<int> A_既知長標本)> l_対称化 = [];
            foreach (var ((l_始点, l_終点), (l_支持数, l_標本)) in l_辺の集計)
            {
                foreach (var l_キー in new[] { (l_始点, l_終点), (l_終点 ^ 1, l_始点 ^ 1) })
                {
                    if (l_対称化.TryGetValue(l_キー, out var l_累積))
                    {
                        l_累積.A_既知長標本.AddRange(l_標本);
                        l_対称化[l_キー] = (l_累積.A_支持数 + l_支持数, l_累積.A_既知長標本);
                    }
                    else
                    {
                        l_対称化[l_キー] = (l_支持数, [.. l_標本]);
                    }
                }
            }
            return l_対称化;
        }

        /// <summary>
        /// インサートサイズを確定する
        /// </summary>
        /// <param name="p_インサートサイズ"></param>
        /// <returns></returns>
        private bool TryGet_インサートサイズ(int p_ライブラリ, out int p_インサートサイズ)
        {
            if (ConfigurationManager.A_実行時引数.A_インサートサイズ is { } l_指定値)
            {
                p_インサートサイズ = l_指定値;
                return true;
            }

            var l_ラベル = p_contig構築.A_ペアのライブラリ数 > 1 ? FormattableString.Invariant($" (ライブラリ {p_ライブラリ + 1})") : string.Empty;

            var l_同一unitig標本 = p_ライブラリ < p_contig構築.A_同一unitig標本群.Count
                ? p_contig構築.A_同一unitig標本群[p_ライブラリ]
                : [];
            if (l_同一unitig標本.Count >= インサートサイズ標本数の下限)
            {
                var l_推定値 = StatsUtil.Get_中央値(l_同一unitig標本);
                var l_unitigN50 = Get_UnitigN50(p_contig構築.A_unitig長);
                if (l_推定値 > 0 && l_unitigN50 >= (long)l_推定値 * 偏りが無いとみなす長さ比)
                {
                    p_インサートサイズ = l_推定値;
                    Logger.V_出力(メッセージID.インサートサイズ推定_同一unitig, p_インサートサイズ, l_同一unitig標本.Count, l_unitigN50, 偏りが無いとみなす長さ比);
                    return true;
                }
            }

            var l_確定辺標本 = p_ライブラリ < p_contig構築.A_確定辺標本群.Count
                ? p_contig構築.A_確定辺標本群[p_ライブラリ]
                : [];
            if (l_確定辺標本.Count >= インサートサイズ標本数の下限)
            {
                p_インサートサイズ = StatsUtil.Get_中央値(l_確定辺標本);
                Logger.V_出力(メッセージID.インサートサイズ推定_確定辺, p_インサートサイズ, l_確定辺標本.Count);
                return true;
            }

            Logger.V_出力_そのまま(FormattableString.Invariant(
                $"[Info] インサートサイズを推定できる標本が足りない{l_ラベル}: 同一 unitig {l_同一unitig標本.Count:N0} 件、確定辺 {l_確定辺標本.Count:N0} 件 (下限 {インサートサイズ標本数の下限})"));
            p_インサートサイズ = 0;
            return false;
        }

        /// <summary>
        /// 採用した辺と、それを支えた一貫したペアの本数をライブラリ別に書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        /// <param name="p_確定辺">頂点ごとの採用辺</param>
        /// <param name="p_ライブラリ別本数">辺ごとの、ライブラリ別の一貫した支持の本数と期待本数</param>
        private void V_書き出し_採用辺(string p_パス, (int A_行き先, int A_ギャップ長)?[] p_確定辺, Dictionary<(int, int), (int[] A_本数, double[] A_期待)> p_ライブラリ別本数)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            l_書き込み.WriteLine("from\tfrom_end\tto\tto_end\tgap\tsupport_by_library\texpected_by_library");
            for (var v = 2; v < p_確定辺.Length; v++)
            {
                if (p_確定辺[v] is not { } l_辺 || !p_ライブラリ別本数.TryGetValue((v, l_辺.A_行き先), out var l_支持))
                {
                    continue;
                }
                var l_始点名 = this._contig名.GetValueOrDefault(v >> 1, string.Empty);
                var l_終点名 = this._contig名.GetValueOrDefault(l_辺.A_行き先 >> 1, string.Empty);
                l_書き込み.WriteLine(FormattableString.Invariant($"{l_始点名}\t{v & 1}\t{l_終点名}\t{l_辺.A_行き先 & 1}\t{l_辺.A_ギャップ長}\t{string.Join(",", l_支持.A_本数)}\t{string.Join(",", l_支持.A_期待.Select(x => x.ToString("F1", System.Globalization.CultureInfo.InvariantCulture)))}"));
            }
        }

        /// <summary>
        /// 支持の下限を満たす候補が 2 本以上ある頂点について、全候補を書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        /// <param name="p_隣接">頂点ごとの候補</param>
        /// <param name="p_確定辺">頂点ごとの採用辺 (相互一意の前)</param>
        /// <param name="p_ライブラリ別本数">辺ごとの、ライブラリ別の一貫した支持の本数と期待本数</param>
        /// <param name="p_最小証拠数">候補として数える支持の下限</param>
        private void V_書き出し_競合候補(string p_パス, List<Scaffold候補>[] p_隣接, (int A_行き先, int A_ギャップ長)?[] p_確定辺, Dictionary<(int, int), (int[] A_本数, double[] A_期待)> p_ライブラリ別本数, ulong p_最小証拠数)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            l_書き込み.WriteLine("from\tfrom_end\tto\tto_end\tsupport\tratio\tsupport_by_library\texpected_by_library\tchosen");
            for (var v = 2; v < p_隣接.Length; v++)
            {
                var l_候補群 = p_隣接[v].Where(x => x.A_支持数 >= p_最小証拠数).ToList();
                if (l_候補群.Count < 2)
                {
                    continue;
                }

                var l_始点名 = this._contig名.GetValueOrDefault(v >> 1, string.Empty);
                foreach (var l_候補 in l_候補群)
                {
                    var l_終点名 = this._contig名.GetValueOrDefault(l_候補.A_行き先 >> 1, string.Empty);
                    var (l_本数, l_期待) = p_ライブラリ別本数.TryGetValue((v, l_候補.A_行き先), out var l_支持) ? l_支持 : ([], []);
                    var l_Is採用 = p_確定辺[v] is { } l_辺 && l_辺.A_行き先 == l_候補.A_行き先;
                    l_書き込み.WriteLine(FormattableString.Invariant($"{l_始点名}\t{v & 1}\t{l_終点名}\t{l_候補.A_行き先 & 1}\t{l_候補.A_支持数}\t{l_候補.A_期待に対する比:F2}\t{string.Join(",", l_本数)}\t{string.Join(",", l_期待.Select(x => x.ToString("F1", System.Globalization.CultureInfo.InvariantCulture)))}\t{(l_Is採用 ? 1 : 0)}"));
                }
            }
        }

        /// <summary>
        /// unitig の N50
        /// </summary>
        /// <param name="p_unitig長"></param>
        /// <returns></returns>
        private static long Get_UnitigN50(IReadOnlyDictionary<int, int> p_unitig長)
        {
            return StatsUtil.Get_N50([.. p_unitig長.Values.Select(x => (long)x)]).A_N50;
        }

        /// <summary>
        /// contig を読み込んで、長さと配列を手元に持つ
        /// </summary>
        private void V_読込_Contig()
        {
            using var l_読み込み = new FastaReader(p_contigファイルパス);
            var l_ID = 1;
            while (l_読み込み.Has続き())
            {
                var l_配列エントリ = l_読み込み.Get_次の配列();
                this._contig名[l_ID] = l_配列エントリ.A_ID.TrimStart('>');
                this._contig配列[l_ID] = l_配列エントリ.A_配列;
                l_ID++;
            }
        }

        /// <summary>
        /// 符号付き unitig ID が contig の末端に配置されているかを判定し、配置されていれば対応する contig 頂点を返す
        /// </summary>
        /// <param name="p_配置"></param>
        /// <param name="p_符号付きunitigID"></param>
        /// <param name="p_Is出口側"></param>
        /// <param name="p_頂点番号"></param>
        /// <returns></returns>
        private static bool TryGet_Contig末端頂点(IReadOnlyDictionary<int, Unitig配置> p_配置, int p_符号付きunitigID, bool p_Is出口側, out int p_頂点番号)
        {
            p_頂点番号 = 0;
            var l_unitigID = Math.Abs(p_符号付きunitigID);
            var l_Is順鎖 = p_符号付きunitigID > 0;

            if (!p_配置.TryGetValue(l_unitigID, out var l_配置情報))
            {
                return false;
            }

            var l_Is実効順鎖 = l_Is順鎖 != l_配置情報.A_Iswalk中逆鎖;

            var l_Is該当端 = p_Is出口側
                ? l_Is実効順鎖 ? l_配置情報.A_IsContig末尾 : l_配置情報.A_IsContig先頭
                : l_Is実効順鎖 ? l_配置情報.A_IsContig先頭 : l_配置情報.A_IsContig末尾;
            if (!l_Is該当端)
            {
                return false;
            }

            var l_Is最終配列順鎖 = l_配置情報.A_IsContig逆相補 ? !l_Is実効順鎖 : l_Is実効順鎖;

            p_頂点番号 = (l_配置情報.A_ContigID << 1) | (l_Is最終配列順鎖 ? 0 : 1);
            return true;
        }

        /// <summary>
        /// 候補のなかで優勢なものを、その頂点の scaffold 辺として確定する
        /// </summary>
        /// <param name="p_隣接">頂点ごとの候補</param>
        /// <param name="p_頂点">確定させる頂点</param>
        /// <param name="p_優勢閾値">優勢とみなす比</param>
        /// <param name="p_最小証拠数">確定に要求する支持数</param>
        /// <param name="p_確定辺">確定した辺の書き留め先</param>
        private static void V_確定_Scaffold辺(List<Scaffold候補>[] p_隣接, int p_頂点, decimal p_優勢閾値, ulong p_最小証拠数, (int A_行き先, int A_ギャップ長)?[] p_確定辺)
        {
            var l_最良 = Get_優勢な候補(p_隣接[p_頂点], p_優勢閾値, p_最小証拠数);
            if (l_最良 is not { } l_辺)
            {
                p_確定辺[p_頂点] = null;
                return;
            }

            p_確定辺[p_頂点] = (l_辺.A_行き先, l_辺.A_ギャップ長);
        }

        /// <summary>
        /// 支持数と期待本数比の下限を満たし、その中で優勢比を超える辺を返す
        /// </summary>
        /// <param name="p_候補"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <returns></returns>
        internal static Scaffold候補? Get_優勢な候補(IReadOnlyList<Scaffold候補> p_候補, decimal p_優勢閾値, ulong p_最小証拠数)
        {
            var l_候補 = p_候補.Where(x => x.A_支持数 >= p_最小証拠数).ToList();
            if (l_候補.Count == 0)
            {
                return null;
            }

            var l_合計 = l_候補.Sum(x => x.A_期待に対する比);
            var l_最良 = l_候補.OrderByDescending(x => x.A_期待に対する比).First();
            return l_合計 <= 0D || (decimal)(l_最良.A_期待に対する比 / l_合計) < p_優勢閾値 ? null : l_最良;
        }

        /// <summary>
        /// 頂点に対応する contig の長さを返す
        /// </summary>
        /// <param name="p_頂点">向き付きの頂点番号</param>
        /// <returns>contig の長さ</returns>
        private long Get_Contig長(int p_頂点)
        {
            return this._contig配列.TryGetValue(p_頂点 >> 1, out var l_配列) ? l_配列.Length : 0;
        }

        /// <summary>
        /// contig を 1 本の scaffold へ連ねる
        /// </summary>
        /// <param name="p_確定辺">確定した辺の書き留め先</param>
        /// <param name="p_始点"></param>
        /// <param name="p_訪問済み"></param>
        /// <param name="p_連結したcontig数"></param>
        /// <returns></returns>
        private string? Get_Scaffold配列((int A_行き先, int A_ギャップ長)?[] p_確定辺, int p_始点, bool[] p_訪問済み, out int p_連結したcontig数)
        {
            p_連結したcontig数 = 0;
            var l_contigID = p_始点 >> 1;
            var l_Is逆鎖 = (p_始点 & 1) == 1;
            if (!this._contig配列.TryGetValue(l_contigID, out var l_配列))
            {
                return null;
            }
            p_連結したcontig数 = 1;

            var l_出力 = new StringBuilder(l_Is逆鎖 ? Util.V_逆相補(l_配列) : l_配列);
            var l_現在 = p_始点;
            V_記録_訪問済み(p_訪問済み, l_現在);
            while (p_確定辺[l_現在] is { } l_辺 && !p_訪問済み[l_辺.A_行き先])
            {
                var l_次のcontigID = l_辺.A_行き先 >> 1;
                var l_Is次が逆鎖 = (l_辺.A_行き先 & 1) == 1;
                if (!this._contig配列.TryGetValue(l_次のcontigID, out var l_次の配列))
                {
                    break;
                }

                var l_次の向き付き配列 = l_Is次が逆鎖 ? Util.V_逆相補(l_次の配列) : l_次の配列;
                var l_重なり長 = l_辺.A_ギャップ長 <= 0 ? Get_畳める重なり長(l_出力, l_次の向き付き配列) : 0;
                if (l_重なり長 > 0)
                {
                    _ = l_出力.Append(l_次の向き付き配列, l_重なり長, l_次の向き付き配列.Length - l_重なり長);
                }
                else
                {
                    _ = l_出力.Append('N', Math.Max(ギャップ長の下限, l_辺.A_ギャップ長));
                    _ = l_出力.Append(l_次の向き付き配列);
                }

                l_現在 = l_辺.A_行き先;
                p_連結したcontig数++;
                V_記録_訪問済み(p_訪問済み, l_現在);
            }

            return l_出力.ToString();
        }

        /// <summary>
        /// 連結の際に畳んでよい重なりの長さを返す
        /// </summary>
        /// <param name="p_出力">ここまでの scaffold 配列</param>
        /// <param name="p_次の配列">繋ぐ向きに直した次の contig 配列</param>
        /// <returns>畳んでよい重なりの長さ (k-1 が一致しなければ、それより短く 15 塩基以上で完全に一致する最長のもの)、畳めないなら 0</returns>
        internal static int Get_畳める重なり長(StringBuilder p_出力, string p_次の配列)
        {
            var l_最長 = Math.Min(ConfigurationManager.A_実行時引数.A_k長 - 1, Math.Min(p_出力.Length, p_次の配列.Length));
            if (l_最長 <= 0)
            {
                return 0;
            }

            var l_k引く1 = ConfigurationManager.A_実行時引数.A_k長 - 1;
            var l_末尾 = p_出力.ToString(p_出力.Length - l_最長, l_最長);
            for (var l_長さ = l_最長; l_長さ == l_k引く1 || l_長さ >= 短い重なりの下限; l_長さ--)
            {
                if (l_末尾.AsSpan(l_最長 - l_長さ).SequenceEqual(p_次の配列.AsSpan(0, l_長さ)))
                {
                    return l_長さ;
                }
            }
            return 0;
        }

        /// <summary>
        /// その contig が環状に閉じたものとして作られたか
        /// </summary>
        /// <param name="p_contigID"></param>
        /// <returns></returns>
        private bool Is環状(int p_contigID)
        {
            return this._contig名.TryGetValue(p_contigID, out var l_名前)
                && l_名前.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// 頂点番号が指す contig の両方の向きの頂点を訪問済みにする
        /// </summary>
        /// <param name="p_訪問済み"></param>
        /// <param name="p_頂点番号"></param>
        private static void V_記録_訪問済み(bool[] p_訪問済み, int p_頂点番号)
        {
            var l_contigID = p_頂点番号 >> 1;
            var l_順鎖 = l_contigID << 1;
            var l_逆鎖 = l_順鎖 | 1;
            if (l_逆鎖 < p_訪問済み.Length)
            {
                p_訪問済み[l_順鎖] = true;
                p_訪問済み[l_逆鎖] = true;
            }
        }

        #endregion
    }
}
