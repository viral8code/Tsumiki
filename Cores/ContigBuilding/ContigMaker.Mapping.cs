using System.Runtime.InteropServices;
using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.IO;
using Tsumiki.Models.ContigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、unitig への k-mer 索引構築とリードマッピングを担う部分
    /// </summary>
    internal partial class ContigMaker
    {
        #region 定数

        /// <summary>
        /// 曖昧 kmer の番兵
        /// </summary>
        private const int 曖昧kmerの番兵 = int.MinValue;

        /// <summary>
        /// read1 と read2 の並びを 1 本に繋ぐのに要る、共有する unitig の数
        /// </summary>
        private const int ペア経路を繋ぐ最小の重なり = 2;

        #endregion

        #region 内部変数

        /// <summary>
        /// k-mer から、それが載る unitig と開始位置を引く辞書
        /// </summary>
        private readonly Dictionary<KmerKey, (int A_unitigID, int A_開始位置)> _kmer辞書;

        /// <summary>
        /// unitig 長
        /// </summary>
        private readonly Dictionary<int, int> _unitig長;

        /// <summary>
        /// 正逆両鎖の unitig 配列
        /// </summary>
        private readonly List<string> _unitig配列;

        /// <summary>
        /// リードが跨いだ unitig の組と、その本数
        /// </summary>
        private readonly Dictionary<(int, int), ulong> _リード隣接;

        /// <summary>
        /// ライブラリごとの、ペアが跨いだ unitig の組と、その間に通った頂点
        /// </summary>
        private readonly List<Dictionary<(int, int), List<int>>> _ペア経路群 = [];

        /// <summary>
        /// 前段 k の確定済み経路 (scaffold/contig 全体) が跨いだ unitig の組と、その本数
        /// </summary>
        private readonly Dictionary<(int, int), ulong> _経路引き継ぎ隣接;

        /// <summary>
        /// リードとペアが通った 3 unitig 以上の並び (符号付き ID、正準の向き) と、その独立な観測数
        /// </summary>
        private readonly Dictionary<経路キー, ulong> _経路集計;

        /// <summary>
        /// 前段 k の確定済み経路が通った 3 unitig 以上の並び
        /// </summary>
        private readonly Dictionary<経路キー, ulong> _引き継ぎ経路集計;

        /// <summary>
        /// unitig 配置
        /// </summary>
        private readonly Dictionary<int, Unitig配置> _unitig配置 = [];

        #endregion

        #region コンストラクタ

        /// <summary>
        /// コンストラクタ
        /// </summary>
        /// <param name="p_unitigファイルパス"></param>
        public ContigMaker(string p_unitigファイルパス)
        {
            this._kmer辞書 = [];
            this._unitig長 = [];
            this._unitig配列 = [string.Empty, string.Empty];
            this._リード隣接 = [];

            this._経路引き継ぎ隣接 = [];
            this._経路集計 = [];
            this._引き継ぎ経路集計 = [];
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            using FastaReader l_読み込み = new(p_unitigファイルパス);
            var l_ID = 1;
            var l_短すぎるunitig数 = 0;
            var l_曖昧数 = 0;
            while (l_読み込み.Has続き())
            {
                var l_unitig = l_読み込み.Get_次の配列();
                this._unitig長[l_ID] = l_unitig.A_配列.Length;
                this._unitig配列.Add(l_unitig.A_配列);
                this._unitig配列.Add(Util.V_逆相補(l_unitig.A_配列));

                if (l_unitig.A_配列.Length < l_k長)
                {
                    l_短すぎるunitig数++;
                    l_ID++;
                    continue;
                }
                for (var i = l_k長; i <= l_unitig.A_配列.Length; i++)
                {
                    var l_開始位置 = i - l_k長;
                    var l_キー = new KmerKey(l_unitig.A_配列.AsSpan(l_開始位置, l_k長));
                    var l_逆鎖キー = l_キー.Get_逆相補();

                    var l_逆鎖開始位置 = l_unitig.A_配列.Length - i;
                    l_曖昧数 += V_登録_kmer(this._kmer辞書, l_キー, l_ID, l_開始位置);
                    l_曖昧数 += V_登録_kmer(this._kmer辞書, l_逆鎖キー, -l_ID, l_逆鎖開始位置);
                }
                l_ID++;
            }

            if (l_短すぎるunitig数 > 0)
            {
                Logger.V_出力(メッセージID.短すぎるunitigの除外, l_短すぎるunitig数);
            }

            if (l_曖昧数 > 0)
            {
                Logger.V_出力(メッセージID.曖昧なkmer登録, l_曖昧数);
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// unitig 間の隣接を de Bruijn グラフから厳密に構築する
        /// </summary>
        /// <returns></returns>
        public UnitigGraph Get_グラフ()
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            return UnitigGraph.Get_グラフ(this._unitig配列, this._kmer辞書, l_k長, 曖昧kmerの番兵);
        }

        /// <summary>
        /// ここまでに貼り付けたリードとペアが通った並びの索引を作る
        /// </summary>
        /// <returns></returns>
        internal ReadPathIndex Get_経路索引()
        {
            return ReadPathIndex.Get_索引(this._経路集計);
        }

        /// <summary>
        /// リードを unitig へ貼り付け、隣接とペアの支持を集める
        /// </summary>
        /// <param name="p_リードパス">貼り付けるリードのパス</param>
        public void V_マッピング_リード(string p_リードパス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];
            var l_ローカル経路 = new Dictionary<経路キー, ulong>[l_スレッド数];
            var l_作業域 = new リード走査作業域[l_スレッド数];
            for (var i = 0; i < l_スレッド数; i++)
            {
                l_ローカル隣接[i] = [];
                l_ローカル経路[i] = [];
                l_作業域[i] = new リード走査作業域();
            }

            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, 塩基列控え.Get_塩基列(p_リードパス), (l_リード, l_ワーカー番号) => this.V_処理_1リード(l_リード, l_ローカル隣接[l_ワーカー番号], l_ローカル経路[l_ワーカー番号], l_作業域[l_ワーカー番号]));

            V_統合_隣接(this._リード隣接, l_ローカル隣接);
            V_統合_経路(this._経路集計, l_ローカル経路);
        }

        /// <summary>
        /// ペアエンドから unitig 間の隣接を検出する
        /// </summary>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        public void V_マッピング_ペアリード(string p_リード1のパス, string p_リード2のパス, int p_ライブラリ番号 = 0)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];
            var l_ローカル経路 = new Dictionary<経路キー, ulong>[l_スレッド数];
            var l_作業域1 = new リード走査作業域[l_スレッド数];
            var l_作業域2 = new リード走査作業域[l_スレッド数];

            var l_ローカルペア経路 = new Dictionary<(int, int), List<int>>[l_スレッド数];

            var l_ローカル同一向き標本 = new List<int>[l_スレッド数];
            var l_ローカル逆向き標本 = new List<int>[l_スレッド数];
            for (var i = 0; i < l_スレッド数; i++)
            {
                l_ローカル隣接[i] = [];
                l_ローカル経路[i] = [];
                l_作業域1[i] = new リード走査作業域();
                l_作業域2[i] = new リード走査作業域();
                l_ローカルペア経路[i] = [];
                l_ローカル同一向き標本[i] = [];
                l_ローカル逆向き標本[i] = [];
            }

            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, Get_ペアリード列(p_リード1のパス, p_リード2のパス), (l_ペア, l_ワーカー番号) => this.V_処理_1ペア(l_ペア.A_リード1, l_ペア.A_リード2, l_ローカル隣接[l_ワーカー番号], l_ローカル経路[l_ワーカー番号], l_作業域1[l_ワーカー番号], l_作業域2[l_ワーカー番号], l_ローカルペア経路[l_ワーカー番号], l_ローカル同一向き標本[l_ワーカー番号], l_ローカル逆向き標本[l_ワーカー番号]));

            V_統合_隣接(this._リード隣接, l_ローカル隣接);
            V_統合_経路(this._経路集計, l_ローカル経路);

            var l_ペア経路 = this.Get_ペア経路(p_ライブラリ番号);
            foreach (var l_ローカルペア in l_ローカルペア経路)
            {
                foreach (var (l_キー, l_値群) in l_ローカルペア)
                {
                    if (l_ペア経路.TryGetValue(l_キー, out var l_一覧))
                    {
                        l_一覧.AddRange(l_値群);
                    }
                    else
                    {
                        l_ペア経路[l_キー] = [.. l_値群];
                    }
                }
            }

            var l_同一向き合計 = l_ローカル同一向き標本.Sum(x => x.Count);
            var l_逆向き合計 = l_ローカル逆向き標本.Sum(x => x.Count);

            IEnumerable<List<int>> l_採用する標本群;
            string l_採用ラベル;

            if (l_同一向き合計 == 0 && l_逆向き合計 == 0)
            {
                l_採用する標本群 = [];
                l_採用ラベル = "none";
            }
            else if (l_同一向き合計 >= l_逆向き合計)
            {
                l_採用する標本群 = l_ローカル同一向き標本;
                l_採用ラベル = "same-orientation";
            }
            else
            {
                l_採用する標本群 = l_ローカル逆向き標本;
                l_採用ラベル = "opposite-orientation";
            }

            var l_同一unitig標本 = new List<int>();
            foreach (var l_標本 in l_採用する標本群)
            {
                l_同一unitig標本.AddRange(l_標本);
            }
            this.A_インサートサイズ標本.AddRange(l_同一unitig標本);
            this.A_同一unitig標本.AddRange(l_同一unitig標本);
            this.Get_同一unitig標本(p_ライブラリ番号).AddRange(l_同一unitig標本);

            var l_ペア支持数 = l_ペア経路.Values.Sum(x => x.Count);
            Logger.V_出力(メッセージID.ペア隣接候補数, l_ペア経路.Count, l_ペア支持数);
            Logger.V_出力(メッセージID.同一unitigのペア向き集計, l_同一向き合計, l_逆向き合計, l_採用ラベル, l_同一unitig標本.Count);
            if (l_同一unitig標本.Count > 0)
            {
                Logger.V_出力(メッセージID.同一unitigの断片長分布, Get_分布要約(l_同一unitig標本));
                Logger.V_出力(メッセージID.同一unitigの断片長中央値, StatsUtil.Get_中央値(l_同一unitig標本), l_同一unitig標本.Count);
            }
        }

        /// <summary>
        /// 前段 k で確定した経路 (scaffold/contig 全体の配列) を、この k の unitig グラフへ再マッピングして隣接の由来にする
        /// </summary>
        /// <param name="p_引き継ぎ経路群">前段 k の確定済み配列</param>
        public void V_マッピング_引き継ぎ経路(IEnumerable<string> p_引き継ぎ経路群)
        {
            var l_件数 = 0;
            var l_作業域 = new リード走査作業域();
            foreach (var l_配列 in p_引き継ぎ経路群)
            {
                this.V_処理_1リード(l_配列, this._経路引き継ぎ隣接, this._引き継ぎ経路集計, l_作業域);
                l_件数++;
            }
            Logger.V_出力(メッセージID.経路引き継ぎのマッピング数, l_件数, this._経路引き継ぎ隣接.Count);
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// 1 本のリードが代表としてどの unitig にマップされるかを判定する
        /// </summary>
        /// <param name="p_リード"></param>
        /// <returns></returns>
        public 代表Unitigヒット Get_代表Unitig(string p_リード)
        {
            return this.Get_走査結果(p_リード, null, new リード走査作業域());
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer 辞書へ 1 件登録する
        /// </summary>
        /// <param name="p_辞書"></param>
        /// <param name="p_キー"></param>
        /// <param name="p_ID"></param>
        /// <param name="p_位置"></param>
        /// <returns></returns>
        private static int V_登録_kmer(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.Item1 == 曖昧kmerの番兵 || l_既存.Item1 == p_ID)
                {
                    return 0;
                }
                p_辞書[p_キー] = (曖昧kmerの番兵, 0);
                return 1;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
            return 0;
        }

        /// <summary>
        /// ワーカーごとの隣接を 1 つにまとめる
        /// </summary>
        /// <param name="p_統合先"></param>
        /// <param name="p_ローカル群"></param>
        private static void V_統合_隣接(Dictionary<(int, int), ulong> p_統合先, Dictionary<(int, int), ulong>[] p_ローカル群)
        {
            foreach (var l_ローカル in p_ローカル群)
            {
                foreach (var (l_キー, l_値) in l_ローカル)
                {
                    p_統合先[l_キー] = p_統合先.GetValueOrDefault(l_キー) + l_値;
                }
            }
        }

        /// <summary>
        /// ワーカーごとの並びの集計を 1 つにまとめる
        /// </summary>
        /// <param name="p_統合先"></param>
        /// <param name="p_ローカル群"></param>
        private static void V_統合_経路(Dictionary<経路キー, ulong> p_統合先, Dictionary<経路キー, ulong>[] p_ローカル群)
        {
            foreach (var l_ローカル in p_ローカル群)
            {
                foreach (var (l_キー, l_値) in l_ローカル)
                {
                    p_統合先[l_キー] = p_統合先.GetValueOrDefault(l_キー) + l_値;
                }
            }
        }

        /// <summary>
        /// read1/read2 を同時に読み進めて対応するペアを返す
        /// </summary>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <returns></returns>
        private static IEnumerable<(string A_リード1, string A_リード2)> Get_ペアリード列(string p_リード1のパス, string p_リード2のパス)
        {
            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);

            var l_Is不一致警告済み = false;
            while (l_読み込み1.Has続き() && l_読み込み2.Has続き())
            {
                var (A_ID1, A_配列1, _) = l_読み込み1.Get_次のレコード();
                var (A_ID2, A_配列2, _) = l_読み込み2.Get_次のレコード();

                if (Util.Get_ペア共通ID(A_ID1) != Util.Get_ペア共通ID(A_ID2))
                {
                    if (!l_Is不一致警告済み)
                    {
                        Logger.V_出力(メッセージID.ペアリードIDの不一致, A_ID1, A_ID2);
                        l_Is不一致警告済み = true;
                    }

                    yield return (A_配列1, string.Empty);
                    yield return (A_配列2, string.Empty);
                    continue;
                }

                yield return (A_配列1, A_配列2);
            }

            while (l_読み込み1.Has続き())
            {
                yield return (l_読み込み1.Get_次のレコード().A_配列, string.Empty);
            }
            while (l_読み込み2.Has続き())
            {
                yield return (l_読み込み2.Get_次のレコード().A_配列, string.Empty);
            }
        }

        /// <summary>
        /// ペア 1 組を処理する
        /// </summary>
        /// <param name="p_リード1"></param>
        /// <param name="p_リード2"></param>
        /// <param name="p_ローカル隣接"></param>
        /// <param name="p_ローカル経路"></param>
        /// <param name="p_作業域1"></param>
        /// <param name="p_作業域2"></param>
        /// <param name="p_ローカルペア経路"></param>
        /// <param name="p_同一向き標本"></param>
        /// <param name="p_逆向き標本"></param>
        private void V_処理_1ペア(string p_リード1, string p_リード2, Dictionary<(int, int), ulong> p_ローカル隣接, Dictionary<経路キー, ulong> p_ローカル経路, リード走査作業域 p_作業域1, リード走査作業域 p_作業域2, Dictionary<(int, int), List<int>> p_ローカルペア経路, List<int> p_同一向き標本, List<int> p_逆向き標本)
        {
            var l_ヒット1 = this.Get_走査結果(p_リード1, p_ローカル隣接, p_作業域1);
            var l_ヒット2 = this.Get_走査結果(p_リード2, p_ローカル隣接, p_作業域2);
            V_集計_ペアの並び(p_作業域1.A_経路, p_作業域2.A_経路, p_ローカル経路);

            if (l_ヒット1.A_unitigID == 0 || l_ヒット2.A_unitigID == 0)
            {
                return;
            }

            if (Math.Abs(l_ヒット1.A_unitigID) == Math.Abs(l_ヒット2.A_unitigID))
            {
                V_収集_同一unitig標本(l_ヒット1, l_ヒット2, p_リード1, p_リード2, p_同一向き標本, p_逆向き標本);
            }
            else
            {
                V_収集_ペア経路(l_ヒット1, l_ヒット2, p_リード1, p_リード2, p_ローカルペア経路);
            }
        }

        /// <summary>
        /// 単一リード 1 本を処理する
        /// </summary>
        /// <param name="p_リード"></param>
        /// <param name="p_ローカル隣接"></param>
        /// <param name="p_ローカル経路"></param>
        /// <param name="p_作業域"></param>
        private void V_処理_1リード(string p_リード, Dictionary<(int, int), ulong> p_ローカル隣接, Dictionary<経路キー, ulong> p_ローカル経路, リード走査作業域 p_作業域)
        {
            _ = this.Get_走査結果(p_リード, p_ローカル隣接, p_作業域);
            V_集計_並び(CollectionsMarshal.AsSpan(p_作業域.A_経路), p_ローカル経路);
        }

        /// <summary>
        /// 並びが 3 unitig 以上なら正準の向きで 1 件数える
        /// </summary>
        /// <param name="p_符号付きID列"></param>
        /// <param name="p_ローカル経路"></param>
        private static void V_集計_並び(ReadOnlySpan<int> p_符号付きID列, Dictionary<経路キー, ulong> p_ローカル経路)
        {
            if (p_符号付きID列.Length < ReadPathIndex.最短の頂点数)
            {
                return;
            }
            var l_キー = new 経路キー(ReadPathIndex.Get_正準経路(p_符号付きID列));
            p_ローカル経路[l_キー] = p_ローカル経路.GetValueOrDefault(l_キー) + 1UL;
        }

        /// <summary>
        /// 1 組のペアの並びを、同じ部分を二重に数えないようにして数える
        /// </summary>
        /// <param name="p_経路1">read1 の並び</param>
        /// <param name="p_経路2">read2 の並び (read2 自身の向き)</param>
        /// <param name="p_ローカル経路"></param>
        private static void V_集計_ペアの並び(List<int> p_経路1, List<int> p_経路2, Dictionary<経路キー, ulong> p_ローカル経路)
        {
            var l_長さ1 = p_経路1.Count;
            var l_長さ2 = p_経路2.Count;
            if (l_長さ1 < ReadPathIndex.最短の頂点数 && l_長さ2 < ReadPathIndex.最短の頂点数 && (l_長さ1 < ペア経路を繋ぐ最小の重なり || l_長さ2 < ペア経路を繋ぐ最小の重なり))
            {
                return;
            }

            var l_並び1 = CollectionsMarshal.AsSpan(p_経路1);
            var l_並び2 = l_長さ2 <= 64 ? stackalloc int[l_長さ2] : new int[l_長さ2];
            for (var i = 0; i < l_長さ2; i++)
            {
                l_並び2[i] = -p_経路2[l_長さ2 - 1 - i];
            }

            var l_重なり = Get_並びの重なり(l_並び1, l_並び2);
            if (l_重なり >= ペア経路を繋ぐ最小の重なり)
            {
                var l_繋いだ並び = new int[l_長さ1 + l_長さ2 - l_重なり];
                l_並び1.CopyTo(l_繋いだ並び);
                l_並び2[l_重なり..].CopyTo(l_繋いだ並び.AsSpan(l_長さ1));
                if (!Has重複unitig(l_繋いだ並び))
                {
                    V_集計_並び(l_繋いだ並び, p_ローカル経路);
                    return;
                }
            }

            if (l_長さ1 >= l_長さ2 && l_並び1.IndexOf(l_並び2) >= 0)
            {
                V_集計_並び(l_並び1, p_ローカル経路);
                return;
            }
            if (l_長さ2 > l_長さ1 && l_並び2.IndexOf(l_並び1) >= 0)
            {
                V_集計_並び(l_並び2, p_ローカル経路);
                return;
            }
            V_集計_並び(l_並び1, p_ローカル経路);
            V_集計_並び(l_並び2, p_ローカル経路);
        }

        /// <summary>
        /// 前の並びの末尾と後の並びの先頭が一致する最長の要素数
        /// </summary>
        /// <param name="p_前"></param>
        /// <param name="p_後"></param>
        /// <returns></returns>
        private static int Get_並びの重なり(ReadOnlySpan<int> p_前, ReadOnlySpan<int> p_後)
        {
            for (var l_重なり = Math.Min(p_前.Length, p_後.Length); l_重なり > 0; l_重なり--)
            {
                if (p_前[^l_重なり..].SequenceEqual(p_後[..l_重なり]))
                {
                    return l_重なり;
                }
            }
            return 0;
        }

        /// <summary>
        /// 並びに同じ unitig がどちらかの向きで 2 度以上現れるか
        /// </summary>
        /// <param name="p_並び"></param>
        /// <returns></returns>
        private static bool Has重複unitig(ReadOnlySpan<int> p_並び)
        {
            for (var i = 0; i < p_並び.Length; i++)
            {
                for (var j = i + 1; j < p_並び.Length; j++)
                {
                    if (Math.Abs(p_並び[i]) == Math.Abs(p_並び[j]))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        /// <summary>
        /// 1 本のリードを k-mer 索引で走査し、隣接を数え、通った unitig の並びと代表 unitig を求める
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_ローカル隣接">このワーカーが集めた隣接、null なら数えない</param>
        /// <param name="p_作業域">通った unitig の並びを受け取る</param>
        /// <returns>最多得票の unitig</returns>
        private 代表Unitigヒット Get_走査結果(string p_リード, Dictionary<(int, int), ulong>? p_ローカル隣接, リード走査作業域 p_作業域)
        {
            p_作業域.V_初期化();
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            if (p_リード.Length < l_k長)
            {
                return 代表Unitigヒット.A_ヒットなし;
            }

            var l_経路 = p_作業域.A_経路;
            var l_票数 = p_作業域.A_票数;
            var l_終端位置 = p_作業域.A_終端位置;
            var l_曖昧塩基数 = 0;

            for (var i = 0; i < l_k長; i++)
            {
                if (Util.Is曖昧塩基(p_リード[i]))
                {
                    l_曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Is曖昧塩基(p_リード[i - l_k長]))
                {
                    l_曖昧塩基数--;
                }

                if (l_曖昧塩基数 != 0)
                {
                    continue;
                }

                var l_キー = new KmerKey(p_リード.AsSpan(i - l_k長, l_k長));
                if (!this._kmer辞書.TryGetValue(l_キー, out var l_項目) || l_項目.A_unitigID == 曖昧kmerの番兵)
                {
                    continue;
                }

                var l_ID = l_項目.A_unitigID;

                var l_終端 = l_項目.A_開始位置 + l_k長;
                if (l_経路.Count > 0 && l_経路[^1] == l_ID)
                {
                    l_票数[^1]++;
                    l_終端位置[^1] = l_終端;
                    continue;
                }

                if (l_経路.Count > 0 && p_ローカル隣接 is not null)
                {
                    var l_経路キー = (l_経路[^1], l_ID);
                    p_ローカル隣接[l_経路キー] = p_ローカル隣接.GetValueOrDefault(l_経路キー) + 1UL;
                }
                l_経路.Add(l_ID);
                l_票数.Add(1);
                l_終端位置.Add(l_終端);
            }

            if (l_経路.Count == 0)
            {
                return 代表Unitigヒット.A_ヒットなし;
            }

            var l_最良 = 0;
            var l_最多得票 = 0;
            var l_最良の終端 = 0;
            for (var s = 0; s < l_経路.Count; s++)
            {
                var l_合計 = 0;
                var l_最後の終端 = 0;
                for (var t = 0; t < l_経路.Count; t++)
                {
                    if (l_経路[t] == l_経路[s])
                    {
                        l_合計 += l_票数[t];
                        l_最後の終端 = l_終端位置[t];
                    }
                }
                if (l_合計 > l_最多得票)
                {
                    l_最良 = l_経路[s];
                    l_最多得票 = l_合計;
                    l_最良の終端 = l_最後の終端;
                }
            }

            var l_unitig長 = this._unitig長.GetValueOrDefault(Math.Abs(l_最良), 0);
            return new 代表Unitigヒット(l_最良, l_最良の終端, l_unitig長);
        }

        #endregion

    }
}
