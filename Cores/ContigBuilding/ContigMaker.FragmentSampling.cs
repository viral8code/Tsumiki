using Tsumiki.Commons;
using Tsumiki.Models.ContigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、フラグメント長・インサートサイズの標本収集を担う部分
    /// </summary>
    /// <remarks>
    /// (unitig へのマッピングは ContigMaker.Mapping.cs、contig 結合そのものは ContigMaker.cs を参照)
    /// </remarks>
    internal partial class ContigMaker
    {
        #region プロパティ

        /// <summary>
        /// 標本抽出された距離の一覧
        /// </summary>
        /// <remarks>
        /// Scaffolder は出所によるバイアスの違いを見るため、この結合ではなく個別の一覧を優先する
        /// </remarks>
        public List<int> A_インサートサイズ標本 { get; } = [];

        /// <summary>
        /// 単一 unitig 内で両リードがヒットしたペアからの標本
        /// </summary>
        /// <remarks>
        /// unitig がフラグメント長より短いと短いフラグメントに偏る
        /// </remarks>
        public List<int> A_同一unitig標本 { get; } = [];

        /// <summary>
        /// unitig 同士が k-1 オーバーラップで直接結合されたペアからの標本
        /// </summary>
        /// <remarks>
        /// 同一 unitig 標本のような長さバイアスを受けない
        /// </remarks>
        public List<int> A_確定辺標本 { get; } = [];

        /// <summary>
        /// ペアエンド由来の隣接候補
        /// </summary>
        /// <remarks>
        /// キーは (始点, 終点) の unitig ID (符号は向き)、値は各観測ペアの既知長の一覧<br/>
        /// Scaffolder から参照される
        /// </remarks>
        public IReadOnlyDictionary<(int, int), List<int>> A_ペア経路 => this._ペア経路;

        /// <summary>
        /// unitig ID (1 始まり、符号なし) からその塩基長を引く
        /// </summary>
        /// <remarks>
        /// Scaffolder が contig 側の末端 unitig の長さを参照する際に使う
        /// </remarks>
        public IReadOnlyDictionary<int, int> A_unitig長 => this._unitig長;

        /// <summary>
        /// unitig 配置
        /// </summary>
        public IReadOnlyDictionary<int, Unitig配置> A_unitig配置 => this._unitig配置;

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 両リードが同一 unitig にマップされたペアから、フラグメント長の標本を取る
        /// </summary>
        /// <param name="p_ヒット1"></param>
        /// <param name="p_ヒット2"></param>
        /// <param name="p_リード1"></param>
        /// <param name="p_リード2"></param>
        /// <param name="p_同一向き標本"></param>
        /// <param name="p_逆向き標本"></param>
        /// <remarks>
        /// 2 ヒットの順鎖座標の差はフラグメント長ではなく、2 リードに挟まれた内側の未読区間 (inner distance) である<br/>
        /// FR 配置では「フラグメント長 = 内側距離 + 両リード長」なので、ここで足し戻して以降の推定値の単位をフラグメント長に揃える
        /// </remarks>
        private static void V_収集_同一unitig標本(代表Unitigヒット p_ヒット1, 代表Unitigヒット p_ヒット2, string p_リード1, string p_リード2, List<int> p_同一向き標本, List<int> p_逆向き標本)
        {
            if ((p_ヒット1.A_unitigID > 0) == (p_ヒット2.A_unitigID > 0))
            {
                // 同じ向き同士 (FF/RR 相当)
                // 両リードの内側の端はどちらも
                // 同じ側を向いているため、差は「開始位置の差」に相当する
                // 下流側リード 1 本分を足すとフラグメント長になる
                var l_内側距離 = Math.Abs(Get_順鎖座標(p_ヒット1) - Get_順鎖座標(p_ヒット2));
                var l_フラグメント長 = l_内側距離 + Math.Max(p_リード1.Length, p_リード2.Length);

                if (l_フラグメント長 > 0)
                {
                    p_同一向き標本.Add(l_フラグメント長);
                }
            }
            else
            {
                // 互いに逆向き (FR 相当、Illumina ペアエンドの通常配置)
                // 順鎖側ヒットのリードがフラグメントの左端、
                // 逆鎖側ヒットのリードが右端を占める
                var l_Isヒット1順鎖 = p_ヒット1.A_unitigID > 0;
                var l_順鎖側の端 = Get_順鎖座標(l_Isヒット1順鎖 ? p_ヒット1 : p_ヒット2);
                var l_逆鎖側の端 = Get_順鎖座標(l_Isヒット1順鎖 ? p_ヒット2 : p_ヒット1);
                var l_順鎖側リード長 = l_Isヒット1順鎖 ? p_リード1.Length : p_リード2.Length;
                var l_逆鎖側リード長 = l_Isヒット1順鎖 ? p_リード2.Length : p_リード1.Length;

                // フラグメントの左端 = 順鎖リードの開始位置、
                // 右端 = 逆鎖リードの終了位置
                var l_フラグメント長 = (l_逆鎖側の端 + l_逆鎖側リード長) - (l_順鎖側の端 - l_順鎖側リード長);

                if (l_フラグメント長 > 0)
                {
                    p_逆向き標本.Add(l_フラグメント長);
                }
            }
        }

        /// <summary>
        /// 両リードが別々の unitig にマップされたペアを隣接候補として記録する
        /// </summary>
        /// <param name="p_ヒット1"></param>
        /// <param name="p_ヒット2"></param>
        /// <param name="p_リード1"></param>
        /// <param name="p_リード2"></param>
        /// <param name="p_ローカルペア経路"></param>
        /// <remarks>
        /// read2 は逆鎖側から読まれるため、read1 の向きへ揃えるには read2 側の unitig ID の符号を反転させる<br/>
        /// 記録するのは「フラグメントのうち既に見えている分の長さ」 (read1 長 + unitig1 末端までの残り + unitig2 先頭からの残り + read2 長) で、ギャップ長 G との間に フラグメント長 = 既知長 + G が常に成り立つ (直接 k-1 で結合された場合は G = - (k-1))
        /// </remarks>
        private static void V_収集_ペア経路(代表Unitigヒット p_ヒット1, 代表Unitigヒット p_ヒット2, string p_リード1, string p_リード2, Dictionary<(int, int), List<int>> p_ローカルペア経路)
        {
            var l_キー = (p_ヒット1.A_unitigID, -p_ヒット2.A_unitigID);

            var l_残り1 = p_ヒット1.A_末尾までの残り長;
            var l_残り2 = Get_反転後残長(p_ヒット2);
            var l_既知長 = l_残り1 + l_残り2 + p_リード1.Length + p_リード2.Length;

            if (p_ローカルペア経路.TryGetValue(l_キー, out var l_一覧))
            {
                l_一覧.Add(l_既知長);
            }
            else
            {
                p_ローカルペア経路[l_キー] = [l_既知長];
            }
        }

        /// <summary>
        /// ヒットを、その unitig を逆向きに見た座標系での残り長に変換する
        /// </summary>
        /// <param name="p_ヒット"></param>
        /// <remarks>
        /// 元の向きでの先頭からの既知長が、逆向きでの残り長にそのまま相当する
        /// </remarks>
        /// <returns></returns>
        private static int Get_反転後残長(代表Unitigヒット p_ヒット)
        {
            return Math.Max(0, p_ヒット.A_最終一致終端位置);
        }

        /// <summary>
        /// ヒットの終端位置を常に順鎖座標系へ揃える
        /// </summary>
        /// <param name="p_ヒット"></param>
        /// <remarks>
        /// 同一 unitig 上の 2 ヒット間の距離を求めるのに使う
        /// </remarks>
        /// <returns></returns>
        private static int Get_順鎖座標(代表Unitigヒット p_ヒット)
        {
            return p_ヒット.A_unitigID > 0 ? p_ヒット.A_最終一致終端位置 : p_ヒット.A_unitig長 - p_ヒット.A_最終一致終端位置;
        }

        /// <summary>
        /// 結合が確定した辺について、ペア経路の既知長から「フラグメント長 = 既知長 - (k-1) 」を計算して標本に積む
        /// </summary>
        /// <param name="p_結合"></param>
        private void V_収集_確定辺標本(int[] p_結合)
        {
            var l_重なり長 = ConfigurationManager.A_実行時引数.A_k長 - 1;
            List<int> l_確定辺標本 = [];

            for (var v = 2; v < p_結合.Length; v++)
            {
                var l_次 = p_結合[v];

                if (l_次 < 0)
                {
                    continue;
                }

                // 頂点番号 -> 符号付き unitig ID
                var l_始点unitig = (v >> 1) * ((v & 1) == 0 ? 1 : -1);
                var l_終点unitig = (l_次 >> 1) * ((l_次 & 1) == 0 ? 1 : -1);

                if (!this._ペア経路.TryGetValue((l_始点unitig, l_終点unitig), out var l_既知長標本))
                {
                    continue;
                }

                foreach (var l_既知長 in l_既知長標本)
                {
                    // 直接結合された辺では 2 つの unitig が k-1 塩基重なるので、
                    // 未知区間の長さは G = - (k-1)
                    // よって
                    // フラグメント長 = 既知長 - (k-1)
                    var l_フラグメント長 = l_既知長 - l_重なり長;
                    if (l_フラグメント長 > 0)
                    {
                        l_確定辺標本.Add(l_フラグメント長);
                    }
                }
            }

            this.A_インサートサイズ標本.AddRange(l_確定辺標本);
            this.A_確定辺標本.AddRange(l_確定辺標本);

            Logger.V_出力(メッセージID.確定辺標本数, l_確定辺標本.Count);

            if (l_確定辺標本.Count > 0)
            {
                // このプールは「unitig 同士が k-1 オーバーラップで直接結合された」
                // ペアのみを対象とするため、同一 unitig 標本のような
                // 「フラグメントが 1 つの unitig に収まる必要がある」制約が
                // なく、短い unitig による短フラグメントへの偏りを受けにくい
                Logger.V_出力(メッセージID.確定辺標本の中央値, StatsUtil.Get_中央値(l_確定辺標本), l_確定辺標本.Count);
            }
        }

        /// <summary>
        /// フラグメント長分布の分位点
        /// </summary>
        /// <remarks>
        /// 橋渡しできる未知区間の長さを決めるのは中央値ではなく分布の上側の裾なので、そこまで出す
        /// </remarks>
        /// <param name="p_値一覧"></param>
        /// <returns></returns>
        private static string Get_分布要約(List<int> p_値一覧)
        {
            var l_整列済み = p_値一覧.OrderBy(x => x).ToList();
            var l_p01 = StatsUtil.Get_分位点(l_整列済み, 0.01D);
            var l_p10 = StatsUtil.Get_分位点(l_整列済み, 0.10D);
            var l_p25 = StatsUtil.Get_分位点(l_整列済み, 0.25D);
            var l_p50 = StatsUtil.Get_分位点(l_整列済み, 0.50D);
            var l_p75 = StatsUtil.Get_分位点(l_整列済み, 0.75D);
            var l_p90 = StatsUtil.Get_分位点(l_整列済み, 0.90D);
            var l_p99 = StatsUtil.Get_分位点(l_整列済み, 0.99D);
            return $"p1={l_p01}, p10={l_p10}, p25={l_p25}, p50={l_p50}, p75={l_p75}, p90={l_p90}, p99={l_p99}, max={l_整列済み[^1]}";
        }

        #endregion
    }
}
