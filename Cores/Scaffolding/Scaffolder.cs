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
using Tsumiki.Utilities;

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
        /// 複数ライブラリのインサートサイズをログに並べるときの区切り
        /// </summary>
        private const string C_インサートサイズの区切り = " / ";

        /// <summary>
        /// 同一 unitig 内標本を信頼してよい「unitig 長 / 推定フラグメント長」の下限比
        /// </summary>
        private const int C_偏りが無いとみなす長さ比 = 10;

        /// <summary>
        /// scaffold 辺に必要なペア数
        /// </summary>
        private const ulong C_Scaffold支持数の下限 = 3UL;

        /// <summary>
        /// 期待本数がこれ未満の候補は、期待に対する比で優劣を決められないため生の支持数で比べる
        /// </summary>
        /// <remarks>
        /// リード長に近い短い contig などは期待本数がほぼ 0 になり、比が 0 か極端に大きくなる
        /// </remarks>
        private const double C_信頼できる期待本数の下限 = 1.0D;

        /// <summary>
        /// インサートサイズ推定に必要な標本数
        /// </summary>
        private const int C_インサートサイズ標本数の下限 = 30;

        /// <summary>
        /// 採用した辺をライブラリ別の支持と共に書き出すファイル名 (scaffold と同じ場所に置く)
        /// </summary>
        private const string C_採用辺の書き出し名 = "scaffold_edges.tsv";

        /// <summary>
        /// 支持の下限を満たす候補が 2 本以上ある頂点の候補を書き出すファイル名 (scaffold と同じ場所に置く)
        /// </summary>
        private const string C_競合候補の書き出し名 = "scaffold_candidates.tsv";

        /// <summary>
        /// 繋ぎ目を確かめるのに要る、繋いだ配列を含むリードの箇所数
        /// </summary>
        private const int C_繋ぎ目の支持数の下限 = 2;

        /// <summary>
        /// 繋ぎ目を確かめる配列に、重なりの両側から含める塩基数の目安
        /// </summary>
        private const int C_繋ぎ目の余白 = 16;

        /// <summary>
        /// 繋ぎ目を確かめる配列に、重なりの両側から含める塩基数の下限
        /// </summary>
        private const int C_繋ぎ目の余白の下限 = 4;

        /// <summary>
        /// 繋ぎ目をリードで埋めるとき、左の末尾・右の先頭から削ってみる長さの上限 (末尾の読み違いを落とす)
        /// </summary>
        private const int C_埋めるときに削る上限 = 20;

        /// <summary>
        /// リードの続きに右の片が現れたとみなす、右の先頭の長さ
        /// </summary>
        internal const int C_右の錨長 = 25;

        /// <summary>
        /// 繋ぎ目を埋めるのに要る、同じ埋め方をするリードの数
        /// </summary>
        private const int C_埋めるのに要るリード数 = 3;

        /// <summary>
        /// 繋ぎ目を埋めるとき、首位の埋め方が次点の何倍以上要るか
        /// </summary>
        private const int C_首位の優勢比 = 3;

        /// <summary>
        /// 錨の続きを見るリードの数の上限
        /// </summary>
        private const int C_続きを見るリード数の上限 = 64;

        /// <summary>
        /// 繋ぎ目を埋めるときに見る、左右の片の長さ
        /// </summary>
        private const int C_埋めるときに見る長さ = 500;

        /// <summary>
        /// 縦に並んだ反復の単位とみなす、右の錨との不一致数の上限
        /// </summary>
        private const int C_似たとみなす不一致数 = 3;

        /// <summary>
        /// 錨より手前を左の片と比べる長さ
        /// </summary>
        private const int C_手前を比べる長さ = 20;

        /// <summary>
        /// 左の片の末尾がゲノムに 2 回あるかを確かめるのに使う、リードの錨より手前の長さ
        /// </summary>
        private const int C_場所を確かめる手前の長さ = 40;

        /// <summary>
        /// 繋ぎ目の近くで縦に並んだ反復を調べる配列をずらす間隔
        /// </summary>
        private const int C_縦の反復を調べる間隔 = 4;

        /// <summary>
        /// リードに出てくる数が一意な配列の何倍以上なら反復とみなすか
        /// </summary>
        private const double C_反復とみなす出現数の比 = 1.8;

        /// <summary>
        /// 一意な配列がリードに出てくる数を測るのに、配列から問い合わせを取る間隔
        /// </summary>
        private const int C_一意の出現数を測る間隔 = 1_000;

        /// <summary>
        /// 一意な配列がリードに出てくる数を測るときの数え上げの上限
        /// </summary>
        private const int C_一意の出現数の上限 = 1_000;

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

        /// <summary>
        /// k-1 の重なりで畳んだ繋ぎ目の数
        /// </summary>
        private int _k引く1で畳んだ数;

        /// <summary>
        /// 重なりを確かめられず、未確認の繋ぎ目の印で繋いだ繋ぎ目の数
        /// </summary>
        private int _確かめられなかった数;

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

            Logger.V_出力(メッセージID.Scaffolding開始_インサートサイズ, string.Join(C_インサートサイズの区切り, l_インサートサイズ群.Where(x => x > 0)));

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
                l_対称化群[l_ライブラリ] = Get_対称化した辺(l_配置, l_ペア経路, ref l_内部を指した数, ref l_未配置を指した数);

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
                    Logger.V_出力(メッセージID.Scaffoldingライブラリ別の概要, l_ライブラリ + 1, l_インサートサイズ群[l_ライブラリ], l_この長さ, l_標本.Count, l_対称化群[l_ライブラリ].Count);
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
                var l_最良ライブラリ = -1;
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
                    l_最良ライブラリ = l_ライブラリ;
                    l_最良比 = l_較正器群[l_ライブラリ].Get_正規化済み支持((ulong)l_一貫した本数, this.Get_Contig長(l_始点), this.Get_Contig長(l_終点), Math.Max(0, l_ギャップ長));
                }

                var l_最良期待本数 = l_最良ライブラリ >= 0 && l_較正器群[l_最良ライブラリ].A_Is使用可能
                    ? l_較正器群[l_最良ライブラリ].Get_期待本数(this.Get_Contig長(l_始点), this.Get_Contig長(l_終点), Math.Max(0, l_最良ギャップ))
                    : double.NaN;
                l_隣接[l_始点].Add(new Scaffold候補(l_終点, (ulong)l_最良本数, l_最良ギャップ, l_最良比, l_最良期待本数));
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
            var l_最小証拠数 = C_Scaffold支持数の下限;

            Logger.V_出力(メッセージID.Scaffold候補辺数, l_候補キー.Count, Messages.Get_文言(l_較正器群.Any(x => x.A_Is使用可能) ? メッセージID.理想本数モデルあり : メッセージID.理想本数モデルなし));

            var l_確定辺 = new (int A_行き先, int A_ギャップ長)?[l_頂点数];
            var l_生の支持数で判定数 = 0;
            var l_生の支持数で採用数 = 0;
            for (var v = 2; v < l_頂点数; v++)
            {
                var l_生の支持数で判定 = Is生の支持数で判定(l_隣接[v], l_最小証拠数);
                V_確定_Scaffold辺(l_隣接, v, l_優勢閾値, l_最小証拠数, l_確定辺);
                if (l_生の支持数で判定)
                {
                    l_生の支持数で判定数++;
                    if (l_確定辺[v] != null)
                    {
                        l_生の支持数で採用数++;
                    }
                }
            }

            Logger.V_出力(メッセージID.Scaffold生の支持数で判定, l_生の支持数で判定数, l_生の支持数で採用数);

            if (l_ライブラリ別本数 is not null)
            {
                this.V_書き出し_競合候補(Path.Combine(Path.GetDirectoryName(Path.GetFullPath(p_scaffoldパス))!, C_競合候補の書き出し名), l_隣接, l_確定辺, l_ライブラリ別本数, l_最小証拠数);
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
                this.V_書き出し_採用辺(Path.Combine(Path.GetDirectoryName(Path.GetFullPath(p_scaffoldパス))!, C_採用辺の書き出し名), l_確定辺, l_ライブラリ別本数);
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
            Logger.V_出力(メッセージID.Scaffold繋ぎ目の判定, this._k引く1で畳んだ数, this._確かめられなかった数);
        }

        /// <summary>
        /// 支持数と期待本数比の下限を満たし、その中で優勢比を超える辺を返す
        /// </summary>
        /// <param name="p_候補"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <returns></returns>
        public static Scaffold候補? Get_優勢な候補(IReadOnlyList<Scaffold候補> p_候補, decimal p_優勢閾値, ulong p_最小証拠数)
        {
            var l_候補 = p_候補.Where(x => x.A_支持数 >= p_最小証拠数).ToList();
            if (l_候補.Count == 0)
            {
                return null;
            }

            if (Is生の支持数で判定(l_候補, p_最小証拠数))
            {
                var l_支持合計 = l_候補.Sum(x => (long)x.A_支持数);
                var l_支持最良 = l_候補.OrderByDescending(x => x.A_支持数).First();
                return l_支持合計 == 0L || (decimal)l_支持最良.A_支持数 / l_支持合計 < p_優勢閾値 ? null : l_支持最良;
            }

            var l_合計 = l_候補.Sum(x => x.A_期待に対する比);
            var l_最良 = l_候補.OrderByDescending(x => x.A_期待に対する比).First();
            return l_合計 <= 0D || (decimal)(l_最良.A_期待に対する比 / l_合計) < p_優勢閾値 ? null : l_最良;
        }

        /// <summary>
        /// 生の支持数で判定するかを返す。支持数の下限を満たす候補のうち、期待に対する比で最良の候補の期待本数が C_信頼できる期待本数の下限 未満か
        /// (比そのものがあてにならない)、期待本数が下限未満の別の候補が最良の候補以上の支持数を持つ (比では評価できないが、同じだけの証拠がある競合) ときに真
        /// </summary>
        /// <param name="p_候補">頂点の候補</param>
        /// <param name="p_最小証拠数">確定に要求する支持数</param>
        /// <returns>生の支持数で判定するなら真</returns>
        internal static bool Is生の支持数で判定(IReadOnlyList<Scaffold候補> p_候補, ulong p_最小証拠数)
        {
            var l_候補 = p_候補.Where(x => x.A_支持数 >= p_最小証拠数).ToList();
            if (l_候補.Count == 0)
            {
                return false;
            }

            var l_最良の番号 = 0;
            for (var i = 1; i < l_候補.Count; i++)
            {
                if (l_候補[i].A_期待に対する比 > l_候補[l_最良の番号].A_期待に対する比)
                {
                    l_最良の番号 = i;
                }
            }

            var l_最良 = l_候補[l_最良の番号];
            return Is期待本数が小さい(l_最良) || l_候補.Where((x, i) => i != l_最良の番号).Any(x => Is期待本数が小さい(x) && x.A_支持数 >= l_最良.A_支持数);
        }

        /// <summary>
        /// 候補の期待本数が分かっていて、C_信頼できる期待本数の下限 未満か
        /// </summary>
        /// <param name="p_候補">候補</param>
        /// <returns>下限未満なら真</returns>
        private static bool Is期待本数が小さい(Scaffold候補 p_候補)
        {
            return !double.IsNaN(p_候補.A_期待本数) && p_候補.A_期待本数 < C_信頼できる期待本数の下限;
        }

        /// <summary>
        /// k-1 より短い重なり (0 を含む) のうち、繋いだ配列がリードに出てくるものがただ 1 つなら、その長さを返す
        /// </summary>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_出力">ここまでの scaffold 配列</param>
        /// <param name="p_次の配列">繋ぐ向きに直した次の contig 配列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_リード長">リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安 (Get_一意の出現数)、分からなければ 0</param>
        /// <returns>確かめた重なりの長さ、無いか 2 つ以上なら null</returns>
        public static int? Get_リードで確かめた重なり長(ReadMinimizerIndex p_リード索引, StringBuilder p_出力, string p_次の配列, int p_k長, int p_リード長, int p_一意の出現数)
        {
            var l_最長 = Math.Min(p_k長 - 2, Math.Min(p_出力.Length, p_次の配列.Length) - C_繋ぎ目の余白の下限);
            if (l_最長 < 0)
            {
                return null;
            }

            var l_末尾 = p_出力.ToString(p_出力.Length - Math.Min(p_出力.Length, l_最長 + ReadMinimizerIndex.C_最短の問い合わせ長), Math.Min(p_出力.Length, l_最長 + ReadMinimizerIndex.C_最短の問い合わせ長));
            int? l_見つけた長さ = null;
            for (var l_長さ = l_最長; l_長さ >= 0; l_長さ--)
            {
                if (!l_末尾.AsSpan(l_末尾.Length - l_長さ).SequenceEqual(p_次の配列.AsSpan(0, l_長さ)))
                {
                    continue;
                }

                var l_余白 = Math.Max(C_繋ぎ目の余白, (ReadMinimizerIndex.C_最短の問い合わせ長 - l_長さ + 1) / 2);
                l_余白 = Math.Min(l_余白, Math.Min((p_リード長 - l_長さ) / 2, Math.Min(l_末尾.Length - l_長さ, p_次の配列.Length - l_長さ)));
                if (l_余白 < C_繋ぎ目の余白の下限 || l_長さ + (2 * l_余白) < ReadMinimizerIndex.C_最短の問い合わせ長)
                {
                    continue;
                }

                var l_繋いだ配列 = string.Concat(l_末尾.AsSpan(l_末尾.Length - l_長さ - l_余白), p_次の配列.AsSpan(l_長さ, l_余白));
                if (p_リード索引.Get_出現数(l_繋いだ配列, C_繋ぎ目の支持数の下限) < C_繋ぎ目の支持数の下限)
                {
                    continue;
                }

                if (l_見つけた長さ is not null)
                {
                    return null;
                }

                l_見つけた長さ = l_長さ;
            }

            if (l_見つけた長さ is not { } l_重なり長)
            {
                return null;
            }

            var l_左 = p_出力.ToString(p_出力.Length - Math.Min(p_出力.Length, C_埋めるときに見る長さ), Math.Min(p_出力.Length, C_埋めるときに見る長さ));
            l_左 = l_左[(l_左.AsSpan().LastIndexOfAny('N', 'n', Consts.未確認の繋ぎ目) + 1)..];
            var l_畳んだ配列 = string.Concat(l_左, p_次の配列.AsSpan(l_重なり長, Math.Min(p_次の配列.Length, C_埋めるときに見る長さ) - l_重なり長));
            return Is両脇まで跨ぐリードが無い(p_リード索引, l_畳んだ配列, l_左.Length - l_重なり長, l_左.Length, p_リード長, p_一意の出現数) ? null : l_重なり長;
        }

        /// <summary>
        /// 繋いだ配列の、繋ぎ目の区間 (左右の片が重なる所) とその両脇の 最短の問い合わせ長 ずつをそのまま含むリードが 繋ぎ目の支持数の下限 に満たないかを返す。
        /// その長さがリード長を超えるときは、両脇と区間の両端 (どれも 最短の問い合わせ長) のどれかが反復 (Is反復の中) かを返す
        /// (区間が反復だと、別のコピーの境目で繋いだ繋ぎ目の短い配列もリードに出てくるが、両脇まで含めた配列は出てこないため)
        /// </summary>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_繋いだ配列">繋ぎ目の前後を繋いだ配列</param>
        /// <param name="p_区間の始まり">繋ぎ目の区間の始まり</param>
        /// <param name="p_区間の終わり">繋ぎ目の区間の終わり (含まない)</param>
        /// <param name="p_リード長">リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安、分からなければ 0</param>
        /// <returns>繋ぐのに要る跨ぐリードが無ければ true</returns>
        private static bool Is両脇まで跨ぐリードが無い(ReadMinimizerIndex p_リード索引, string p_繋いだ配列, int p_区間の始まり, int p_区間の終わり, int p_リード長, int p_一意の出現数)
        {
            const int l_錨長 = ReadMinimizerIndex.C_最短の問い合わせ長;
            var l_左の始まり = p_区間の始まり - l_錨長;
            var l_跨ぐ長さ = p_区間の終わり + l_錨長 - l_左の始まり;
            return l_左の始まり >= 0 && p_区間の終わり + l_錨長 <= p_繋いだ配列.Length && (l_跨ぐ長さ <= p_リード長
                ? p_リード索引.Get_出現数(p_繋いだ配列.AsSpan(l_左の始まり, l_跨ぐ長さ), C_繋ぎ目の支持数の下限) < C_繋ぎ目の支持数の下限
                : Is反復の中(p_リード索引, p_繋いだ配列, l_左の始まり, p_一意の出現数) || Is反復の中(p_リード索引, p_繋いだ配列, p_区間の終わり, p_一意の出現数)
                    || Is反復の中(p_リード索引, p_繋いだ配列, p_区間の始まり, p_一意の出現数) || Is反復の中(p_リード索引, p_繋いだ配列, p_区間の終わり - l_錨長, p_一意の出現数));
        }

        /// <summary>
        /// 配列の p_位置 からの 最短の問い合わせ長 の配列が反復か (近くに似た配列があるか、リードに出てくる数が一意な配列の 反復とみなす出現数の比 倍以上か)
        /// </summary>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_配列"></param>
        /// <param name="p_位置"></param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安、分からなければ 0</param>
        /// <returns></returns>
        internal static bool Is反復の中(ReadMinimizerIndex p_リード索引, string p_配列, int p_位置, int p_一意の出現数)
        {
            if (Has似た配列(p_配列, p_配列.AsSpan(p_位置, C_右の錨長), p_位置))
            {
                return true;
            }

            if (p_一意の出現数 <= 0)
            {
                return false;
            }

            var l_下限 = (int)Math.Ceiling(p_一意の出現数 * C_反復とみなす出現数の比);
            return p_リード索引.Get_出現数(p_配列.AsSpan(p_位置, ReadMinimizerIndex.C_最短の問い合わせ長), l_下限) >= l_下限;
        }

        /// <summary>
        /// 配列群から 一意の出現数を測る間隔 おきに取った 最短の問い合わせ長 の配列 (N を含まないもの) がリードに出てくる数の中央値を返す
        /// </summary>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_配列群"></param>
        /// <returns>中央値、測れなければ 0</returns>
        public static int Get_一意の出現数(ReadMinimizerIndex p_リード索引, IEnumerable<string> p_配列群)
        {
            List<int> l_数群 = [];
            foreach (var l_配列 in p_配列群)
            {
                for (var l_位置 = 0; l_位置 + ReadMinimizerIndex.C_最短の問い合わせ長 <= l_配列.Length; l_位置 += C_一意の出現数を測る間隔)
                {
                    var l_断片 = l_配列.AsSpan(l_位置, ReadMinimizerIndex.C_最短の問い合わせ長);
                    if (l_断片.IndexOfAnyExcept("ACGT") < 0)
                    {
                        l_数群.Add(p_リード索引.Get_出現数(l_断片, C_一意の出現数の上限));
                    }
                }
            }

            l_数群.Sort();
            return l_数群.Count == 0 ? 0 : l_数群[l_数群.Count / 2];
        }

        /// <summary>
        /// 繋ぎ目を、両側を跨ぐリードの続きで埋める。左の末尾を錨にして右へ読む向きで埋められなければ、右の先頭を錨にして左へ読む向きでも試す
        /// (GC に富む所の読み違いは読む向きで出方が違い、片方の向きのリードは崩れていても逆の向きのリードは読めていることがある)
        /// </summary>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_出力">ここまでの配列</param>
        /// <param name="p_次の配列">繋ぐ向きに直した次の片</param>
        /// <param name="p_リード長">リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安 (Get_一意の出現数)、分からなければ 0</param>
        /// <param name="p_要るリード数">埋めるのに要る、同じ埋め方をするリードの数 (既定は C_埋めるのに要るリード数)</param>
        /// <param name="p_削る上限">左の末尾・右の先頭から削ってみる長さの上限 (既定は C_埋めるときに削る上限)</param>
        /// <returns>左の末尾から削る長さ・間に入れる配列・右の先頭から削る長さ、埋められなければ null</returns>
        public static (int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)? Get_リードで埋めた繋ぎ目(ReadMinimizerIndex p_リード索引, StringBuilder p_出力, string p_次の配列, int p_リード長, int p_一意の出現数, int p_要るリード数 = C_埋めるのに要るリード数, int p_削る上限 = C_埋めるときに削る上限)
        {
            var l_左 = p_出力.ToString(p_出力.Length - Math.Min(p_出力.Length, C_埋めるときに見る長さ), Math.Min(p_出力.Length, C_埋めるときに見る長さ));
            l_左 = l_左[(l_左.AsSpan().LastIndexOfAny('N', 'n', Consts.未確認の繋ぎ目) + 1)..];
            var l_右 = p_次の配列[..Math.Min(p_次の配列.Length, C_埋めるときに見る長さ)];
            return Get_片側から埋める方法(p_リード索引, l_左, l_右, p_リード長, p_一意の出現数, p_要るリード数, p_削る上限) is { } l_右へ
                ? l_右へ
                : Get_片側から埋める方法(p_リード索引, Util.V_逆相補_曖昧塩基あり(l_右), Util.V_逆相補_曖昧塩基あり(l_左), p_リード長, p_一意の出現数, p_要るリード数, p_削る上限) is { } l_左へ
                ? (l_左へ.A_右から削る長さ, Util.V_逆相補_曖昧塩基あり(l_左へ.A_埋める配列), l_左へ.A_左から削る長さ)
                : null;
        }

        /// <summary>
        /// 左の末尾を少しずつ削って錨にし、錨を含むリードの続きに右の先頭 (これも少しずつ削る) が現れるまでの配列で繋ぎ目を埋める。
        /// 錨より手前も左の片と (読み違いを除いて) 一致するリードだけを使う (錨がゲノムの別の場所にもあると、そこから来たリードが別の続きを持ち込むため)。
        /// さらに、手前が 場所を確かめる手前の長さ まで左の片と一致するのに、採る埋め方とは揃って別の続きを持つリードが 埋めるのに要るリード数 以上あれば、左の片の末尾がゲノムに 2 回ある所なので埋めない。
        /// 右の先頭が続きに 1 回だけ現れ、そこから先の続きも右の片とそのまま一致し、その一致が錨の先へ 繋ぎ目の余白 以上及ぶリードだけを数える
        /// (隙間の中にある右の先頭と同じ配列で間を飛ばさないため。錨の直後で終わるリードは、錨と右の先頭が同じ配列なら証拠なしに一致してしまう)。
        /// 同じ埋め方をするリードが 埋めるのに要るリード数 以上あり、次点の 首位の優勢比 倍以上のときだけ採る (末尾の読み違い・重なり・短い隙間をまとめて扱う)。
        /// 埋めた所が縦に並んだ反復の中 (Is縦の反復の中) か、左右の片が重なる埋め方で両脇まで跨ぐリードが無い (Is両脇まで跨ぐリードが無い) ときは埋めない
        /// </summary>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_左">左の片の末尾</param>
        /// <param name="p_右">右の片の先頭</param>
        /// <param name="p_リード長">リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安、分からなければ 0</param>
        /// <param name="p_要るリード数">埋めるのに要る、同じ埋め方をするリードの数</param>
        /// <param name="p_削る上限">左の末尾・右の先頭から削ってみる長さの上限</param>
        /// <returns>左の末尾から削る長さ・間に入れる配列・右の先頭から削る長さ、埋められなければ null</returns>
        private static (int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)? Get_片側から埋める方法(ReadMinimizerIndex p_リード索引, string p_左, string p_右, int p_リード長, int p_一意の出現数, int p_要るリード数, int p_削る上限)
        {
            const int l_錨長 = ReadMinimizerIndex.C_最短の問い合わせ長;
            for (var l_削る = 0; l_削る <= p_削る上限 && p_左.Length >= l_削る + l_錨長 + C_繋ぎ目の余白; l_削る++)
            {
                var l_錨 = p_左.Substring(p_左.Length - l_削る - l_錨長, l_錨長);
                var l_手前 = p_左[..(p_左.Length - l_削る - l_錨長)];
                var l_前後群 = p_リード索引.Get_前後群(l_錨, p_リード長, C_続きを見るリード数の上限).Where(x => Is手前が同じ場所(l_手前, x.A_前, C_手前を比べる長さ)).ToList();
                if (l_前後群.Count < p_要るリード数)
                {
                    continue;
                }

                Dictionary<(int A_右から削る長さ, string A_埋める配列), int> l_票 = [];
                var l_投票 = new (int A_右から削る長さ, string A_埋める配列)?[l_前後群.Count];
                for (var l_番号 = 0; l_番号 < l_前後群.Count; l_番号++)
                {
                    var l_続き = l_前後群[l_番号].A_続き;
                    var l_錨から = l_錨 + l_続き;
                    for (var l_右を削る = 0; l_右を削る <= p_削る上限 && l_右を削る + C_右の錨長 <= p_右.Length; l_右を削る++)
                    {
                        var l_右の錨 = p_右.AsSpan(l_右を削る, C_右の錨長);
                        var l_位置 = l_錨から.AsSpan().IndexOf(l_右の錨);
                        if (l_位置 < 0)
                        {
                            continue;
                        }

                        var l_残り = Math.Min(l_錨から.Length - l_位置, p_右.Length - l_右を削る);
                        if (l_位置 + l_残り >= l_錨長 + C_繋ぎ目の余白 && l_錨から.AsSpan(l_位置 + 1).IndexOf(l_右の錨) < 0 && l_錨から.AsSpan(l_位置, l_残り).SequenceEqual(p_右.AsSpan(l_右を削る, l_残り)))
                        {
                            var l_キー = (l_右を削る, l_錨から[..l_位置]);
                            l_票[l_キー] = l_票.GetValueOrDefault(l_キー) + 1;
                            l_投票[l_番号] = l_キー;
                        }

                        break;
                    }
                }

                if (l_票.Count == 0)
                {
                    continue;
                }

                var l_並び = l_票.OrderByDescending(x => x.Value).ToList();
                var l_次点 = l_並び.Count > 1 ? l_並び[1].Value : 0;
                var (l_右を削る長さ, l_埋める配列) = l_並び[0].Key;
                if (l_並び[0].Value < p_要るリード数 || l_並び[0].Value < C_首位の優勢比 * l_次点)
                {
                    return null;
                }

                var l_見えたはずの長さ = Math.Max(C_右の錨長, l_埋める配列.Length - l_錨長 + C_右の錨長);
                var l_別の続きの最多 = Enumerable.Range(0, l_前後群.Count)
                    .Where(x => l_投票[x] != l_並び[0].Key && l_前後群[x].A_続き.Length >= l_見えたはずの長さ && l_前後群[x].A_前.Length >= C_場所を確かめる手前の長さ && Is手前が同じ場所(l_手前, l_前後群[x].A_前, C_場所を確かめる手前の長さ))
                    .GroupBy(x => l_前後群[x].A_続き[..C_右の錨長]).Select(x => x.Count()).DefaultIfEmpty(0).Max();
                if (l_別の続きの最多 >= p_要るリード数)
                {
                    return null;
                }

                var l_繋いだ配列 = string.Concat(p_左.AsSpan(0, p_左.Length - l_削る - l_錨長), l_埋める配列, p_右.AsSpan(l_右を削る長さ));
                var l_右の錨の位置 = p_左.Length - l_削る - l_錨長 + l_埋める配列.Length;
                var l_左の終わり = p_左.Length - l_削る;
                return Is縦の反復の中(l_繋いだ配列, p_左.Length - l_削る - l_錨長, l_右の錨の位置)
                    || (l_右の錨の位置 < l_左の終わり && Is両脇まで跨ぐリードが無い(p_リード索引, l_繋いだ配列, l_右の錨の位置, l_左の終わり, p_リード長, p_一意の出現数))
                    ? null
                    : (l_削る + l_錨長, l_埋める配列, l_右を削る長さ);
            }

            return null;
        }

        /// <summary>
        /// リードの錨の直前 (p_比べる長さ まで) が、左の片の錨の直前と 似たとみなす不一致数 以内で一致するか (読み違いは許し、錨のゲノムの別のコピーから来たリードを除く)
        /// </summary>
        /// <param name="p_左の手前">左の片の、錨より手前</param>
        /// <param name="p_リードの手前">リードの、錨より手前</param>
        /// <param name="p_比べる長さ"></param>
        /// <returns></returns>
        private static bool Is手前が同じ場所(string p_左の手前, string p_リードの手前, int p_比べる長さ)
        {
            var l_長さ = Math.Min(p_比べる長さ, Math.Min(p_左の手前.Length, p_リードの手前.Length));
            var l_不一致 = 0;
            for (var i = 1; i <= l_長さ && l_不一致 <= C_似たとみなす不一致数; i++)
            {
                if (p_左の手前[^i] != p_リードの手前[^i])
                {
                    l_不一致++;
                }
            }

            return l_不一致 <= C_似たとみなす不一致数;
        }

        /// <summary>
        /// 繋いだ配列の、錨の始まりから右の錨の先 (錨の長さ分) までの区間を 縦の反復を調べる間隔 ずつずらした 右の錨長 の配列のどれかに、似た配列が別の場所にあるか。
        /// あれば縦に並んだ反復の中の繋ぎ目で、リードで跨げず単位の数を取り違えるので埋めない
        /// </summary>
        /// <param name="p_繋いだ配列"></param>
        /// <param name="p_錨の位置">繋いだ配列での左の錨の始まり</param>
        /// <param name="p_右の錨の位置">繋いだ配列での右の錨の始まり</param>
        /// <returns></returns>
        private static bool Is縦の反復の中(string p_繋いだ配列, int p_錨の位置, int p_右の錨の位置)
        {
            var l_終わり = Math.Min(p_繋いだ配列.Length - C_右の錨長, p_右の錨の位置 + ReadMinimizerIndex.C_最短の問い合わせ長);
            for (var l_位置 = Math.Min(p_錨の位置, p_右の錨の位置); l_位置 <= l_終わり; l_位置 += C_縦の反復を調べる間隔)
            {
                if (Has似た配列(p_繋いだ配列, p_繋いだ配列.AsSpan(l_位置, C_右の錨長), l_位置))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// 配列の中に、錨と 似たとみなす不一致数 以内で一致する箇所が (錨そのものの位置を除いて) あるか
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_錨"></param>
        /// <param name="p_除く位置">錨そのものの位置、無ければ -1</param>
        /// <returns></returns>
        private static bool Has似た配列(string p_配列, ReadOnlySpan<char> p_錨, int p_除く位置)
        {
            for (var i = 0; i + p_錨.Length <= p_配列.Length; i++)
            {
                if (i == p_除く位置)
                {
                    continue;
                }

                var l_不一致 = 0;
                for (var j = 0; j < p_錨.Length && l_不一致 <= C_似たとみなす不一致数; j++)
                {
                    if (p_配列[i + j] != p_錨[j])
                    {
                        l_不一致++;
                    }
                }

                if (l_不一致 <= C_似たとみなす不一致数)
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// 連結の際に畳んでよい重なりの長さを返す
        /// </summary>
        /// <param name="p_出力">ここまでの scaffold 配列</param>
        /// <param name="p_次の配列">繋ぐ向きに直した次の contig 配列</param>
        /// <returns>畳んでよい重なりの長さ、畳めないなら 0</returns>
        public static int Get_畳める重なり長(StringBuilder p_出力, string p_次の配列)
        {
            var l_重なり長 = ConfigurationManager.A_実行時引数.A_k長 - 1;
            return l_重なり長 <= 0 || p_出力.Length < l_重なり長 || p_次の配列.Length < l_重なり長
                ? 0
                : p_出力.ToString(p_出力.Length - l_重なり長, l_重なり長) == p_次の配列[..l_重なり長] ? l_重なり長 : 0;
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
        private static Dictionary<(int, int), (ulong A_支持数, List<int> A_既知長標本)> Get_対称化した辺(
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
        /// <param name="p_ライブラリ"></param>
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
            if (l_同一unitig標本.Count >= C_インサートサイズ標本数の下限)
            {
                var l_推定値 = StatsUtil.Get_中央値(l_同一unitig標本);
                var l_unitigN50 = Get_UnitigN50(p_contig構築.A_unitig長);
                if (l_推定値 > 0 && l_unitigN50 >= (long)l_推定値 * C_偏りが無いとみなす長さ比)
                {
                    p_インサートサイズ = l_推定値;
                    Logger.V_出力(メッセージID.インサートサイズ推定_同一unitig, p_インサートサイズ, l_同一unitig標本.Count, l_unitigN50, C_偏りが無いとみなす長さ比);
                    return true;
                }
            }

            var l_確定辺標本 = p_ライブラリ < p_contig構築.A_確定辺標本群.Count
                ? p_contig構築.A_確定辺標本群[p_ライブラリ]
                : [];
            if (l_確定辺標本.Count >= C_インサートサイズ標本数の下限)
            {
                p_インサートサイズ = StatsUtil.Get_中央値(l_確定辺標本);
                Logger.V_出力(メッセージID.インサートサイズ推定_確定辺, p_インサートサイズ, l_確定辺標本.Count);
                return true;
            }

            Logger.V_出力_そのまま(FormattableString.Invariant(
                $"[Info] インサートサイズを推定できる標本が足りない{l_ラベル}: 同一 unitig {l_同一unitig標本.Count:N0} 件、確定辺 {l_確定辺標本.Count:N0} 件 (下限 {C_インサートサイズ標本数の下限})"));
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
                var l_重なり長 = l_辺.A_ギャップ長 <= 0 ? this.Get_繋ぎ目の重なり長(l_出力, l_次の向き付き配列) : null;
                if (l_重なり長 is { } l_長さ)
                {
                    _ = l_出力.Append(l_次の向き付き配列, l_長さ, l_次の向き付き配列.Length - l_長さ);
                }
                else
                {
                    _ = l_辺.A_ギャップ長 <= 0 ? l_出力.Append(Consts.未確認の繋ぎ目) : l_出力.Append('N', l_辺.A_ギャップ長);
                    _ = l_出力.Append(l_次の向き付き配列);
                }

                l_現在 = l_辺.A_行き先;
                p_連結したcontig数++;
                V_記録_訪問済み(p_訪問済み, l_現在);
            }

            return l_出力.ToString();
        }

        /// <summary>
        /// 重なると見込まれる繋ぎ目で、畳む重なりの長さを決める
        /// </summary>
        /// <param name="p_出力">ここまでの scaffold 配列</param>
        /// <param name="p_次の配列">繋ぐ向きに直した次の contig 配列</param>
        /// <returns>畳む重なりの長さ、k-1 の重なりが一致しなければ null (最終出力でリードで確かめる)</returns>
        private int? Get_繋ぎ目の重なり長(StringBuilder p_出力, string p_次の配列)
        {
            var l_k引く1 = Get_畳める重なり長(p_出力, p_次の配列);
            if (l_k引く1 > 0)
            {
                this._k引く1で畳んだ数++;
                return l_k引く1;
            }

            this._確かめられなかった数++;
            return null;
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
