using Tsumiki.Common;
using Tsumiki.Model.ContigBuilding;
using Tsumiki.Model.Foundation;
using Tsumiki.Model.UnitigBuilding;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、フラグメント長・インサートサイズの標本収集を担う部分<br/>
    /// (unitig へのマッピングは ContigMaker.Mapping.cs、contig 結合そのものは
    /// ContigMaker.cs を参照)
    /// </summary>
    internal partial class ContigMaker
    {
        /// <summary>
        /// 標本抽出された距離の一覧<br/>
        /// Scaffolder は出所によるバイアスの違いを
        /// 見るため、この結合ではなく個別の一覧を優先する
        /// </summary>
        public List<int> A_インサートサイズ標本 { get; } = [];

        /// <summary>
        /// 単一 unitig 内で両リードがヒットしたペアからの標本<br/>
        /// unitig がフラグメント長より短いと短いフラグメントに偏る
        /// </summary>
        public List<int> A_同一ユニティグ標本 { get; } = [];

        /// <summary>
        /// unitig 同士が k-1 オーバーラップで直接結合されたペアからの標本<br/>
        /// 同一ユニティグ標本のような長さバイアスを受けない
        /// </summary>
        public List<int> A_確定辺標本 { get; } = [];

        /// <summary>
        /// ペアエンド由来の隣接候補<br/>
        /// キーは (始点, 終点) の unitig ID(符号は向き)、
        /// 値は各観測ペアの既知長の一覧<br/>
        /// Scaffolder から参照される
        /// </summary>
        public IReadOnlyDictionary<(int, int), List<int>> A_ペア経路 => this._ペア経路;

        /// <summary>
        /// unitig ID(1 始まり、符号なし) からその塩基長を引く<br/>
        /// Scaffolder が
        /// contig 側の末端 unitig の長さを参照する際に使う
        /// </summary>
        public IReadOnlyDictionary<int, int> A_ユニティグ長 => this._ユニティグ長;

        public IReadOnlyDictionary<int, ユニティグ配置> A_ユニティグ配置 => this._ユニティグ配置;

        /// <summary>
        /// 両リードが同一 unitig にマップされたペアから、フラグメント長の標本を取る<br/>
        /// 2 ヒットの順鎖座標の差はフラグメント長ではなく、2 リードに挟まれた
        /// 内側の未読区間 (inner distance) である<br/>
        /// FR 配置では
        /// フラグメント長 = 内側距離 + 両リード長 なので、ここで足し戻して
        /// 以降の推定値の単位をフラグメント長に揃える
        /// </summary>
        private static void V_収集_同一ユニティグ標本(
            代表ユニティグヒット p_ヒット1, 代表ユニティグヒット p_ヒット2,
            string p_リード1, string p_リード2,
            List<int> p_同一向き標本, List<int> p_逆向き標本)
        {
            if ((p_ヒット1.A_ユニティグID > 0) == (p_ヒット2.A_ユニティグID > 0))
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
                var l_ヒット1が順鎖か = p_ヒット1.A_ユニティグID > 0;
                var l_順鎖側の端 = Get_順鎖座標(l_ヒット1が順鎖か ? p_ヒット1 : p_ヒット2);
                var l_逆鎖側の端 = Get_順鎖座標(l_ヒット1が順鎖か ? p_ヒット2 : p_ヒット1);
                var l_順鎖側リード長 = l_ヒット1が順鎖か ? p_リード1.Length : p_リード2.Length;
                var l_逆鎖側リード長 = l_ヒット1が順鎖か ? p_リード2.Length : p_リード1.Length;

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
        /// 両リードが別々の unitig にマップされたペアを隣接候補として記録する<br/>
        /// read2 は逆鎖側から読まれるため、read1 の向きへ揃えるには
        /// read2 側の unitig ID の符号を反転させる<br/>
        /// 記録するのは「フラグメントのうち既に見えている分の長さ」
        /// (read1長 + unitig1末端までの残り + unitig2先頭からの残り + read2長) で、
        /// ギャップ長 G との間に フラグメント長 = 既知長 + G が常に成り立つ
        /// (直接 k-1 で結合された場合は G = -(k-1))
        /// </summary>
        private static void V_収集_ペア経路(
            代表ユニティグヒット p_ヒット1, 代表ユニティグヒット p_ヒット2,
            string p_リード1, string p_リード2,
            Dictionary<(int, int), List<int>> p_ローカルペア経路)
        {
            var l_キー = (p_ヒット1.A_ユニティグID, -p_ヒット2.A_ユニティグID);

            var l_残り1 = p_ヒット1.A_末尾までの残り長;
            var l_残り2 = Get_反転後の残り長(p_ヒット2);
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
        /// ヒットを、その unitig を逆向きに見た座標系での残り長に変換する<br/>
        /// 元の向きでの先頭からの既知長が、逆向きでの残り長にそのまま相当する
        /// </summary>
        private static int Get_反転後の残り長(代表ユニティグヒット p_ヒット)
        {
            return Math.Max(0, p_ヒット.A_最終一致終端位置);
        }

        /// <summary>
        /// ヒットの終端位置を常に順鎖座標系へ揃える<br/>
        /// 同一 unitig 上の 2 ヒット間の距離を求めるのに使う
        /// </summary>
        private static int Get_順鎖座標(代表ユニティグヒット p_ヒット)
        {
            return p_ヒット.A_ユニティグID > 0
                ? p_ヒット.A_最終一致終端位置
                : p_ヒット.A_ユニティグ長 - p_ヒット.A_最終一致終端位置;
        }

        /// <summary>
        /// 結合が確定した辺について、ペア経路の既知長から
        /// フラグメント長 = 既知長 - (k-1) を計算して標本に積む
        /// </summary>
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
                var l_始点ユニティグ = (v >> 1) * ((v & 1) == 0 ? 1 : -1);
                var l_終点ユニティグ = (l_次 >> 1) * ((l_次 & 1) == 0 ? 1 : -1);

                if (!this._ペア経路.TryGetValue((l_始点ユニティグ, l_終点ユニティグ), out var l_既知長標本))
                {
                    continue;
                }

                foreach (var l_既知長 in l_既知長標本)
                {
                    // 直接結合された辺では 2 つの unitig が k-1 塩基重なるので、
                    // 未知区間の長さは G = -(k-1)
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
        /// フラグメント長分布の分位点<br/>
        /// 橋渡しできる未知区間の長さを決めるのは
        /// 中央値ではなく分布の上側の裾なので、そこまで出す
        /// </summary>
        private static string Get_分布要約(List<int> p_値一覧)
        {
            var l_整列済み = p_値一覧.OrderBy(x => x).ToList();
            return $"p1={StatsUtil.Get_分位点(l_整列済み, 0.01)}, p10={StatsUtil.Get_分位点(l_整列済み, 0.10)}, " +
                $"p25={StatsUtil.Get_分位点(l_整列済み, 0.25)}, p50={StatsUtil.Get_分位点(l_整列済み, 0.50)}, " +
                $"p75={StatsUtil.Get_分位点(l_整列済み, 0.75)}, p90={StatsUtil.Get_分位点(l_整列済み, 0.90)}, " +
                $"p99={StatsUtil.Get_分位点(l_整列済み, 0.99)}, max={l_整列済み[^1]}";
        }
    }
}
