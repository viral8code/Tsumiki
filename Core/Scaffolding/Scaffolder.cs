using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// 確定した contig を読み直し、ペアエンド由来の隣接で N 埋め連結する。
    /// 出力は新規ファイルで、contigs.fasta 自体は変更しない。
    /// </summary>
    internal class Scaffolder(ContigMaker p_コンティグ構築, string p_コンティグファイルパス, int? p_リード長)
    {
        /// <summary>
        /// 同一 unitig 内標本を信頼してよい「unitig長 / 推定フラグメント長」の下限比。
        /// unitig がフラグメントより短いと両端が収まるペアしか観測できず
        /// 短い側へ偏るが、この倍率以上に長ければ打ち切りは事実上起きない。
        /// </summary>
        private const int 偏りが無いとみなす長さ比 = 10;

        // contig ID(FastaWriter が振った 1 始まりの ID) -> 配列本体。
        private readonly Dictionary<int, string> _コンティグ配列 = [];

        // contig ID -> ID 文字列(先頭 ">" の次に書かれていた文字列。"NODE1" 等)。
        // 出力時に元の命名をある程度踏襲するために保持する。
        private readonly Dictionary<int, string> _コンティグ名 = [];

        /// <summary>
        /// 自動推定された(あるいは CLI で明示指定された)インサートサイズ。
        /// 推定に失敗した場合は null のままとなり、その場合スキャフォールディングは
        /// 行われない。
        /// </summary>
        public int? A_有効インサートサイズ { get; private set; }

        /// <summary>
        /// スキャフォールディングを実行し、指定パスに結果を書き出す。
        /// インサートサイズが(指定・推定いずれの方法でも)確定できなかった場合は、
        /// その旨をログに出力して何もせずに戻る(ファイルは作成されない)。
        /// </summary>
        public void V_実行(string p_スキャフォールドパス)
        {
            if (!this.Get_インサートサイズ(out var l_インサートサイズ))
            {
                Logger.V_出力(メッセージID.スキャフォールディング省略_インサートサイズ不明);
                return;
            }
            this.A_有効インサートサイズ = l_インサートサイズ;
            Logger.V_出力(メッセージID.スキャフォールディング開始_インサートサイズ, l_インサートサイズ);

            this.V_読込_コンティグ();

            if (this._コンティグ配列.Count == 0)
            {
                Logger.V_出力(メッセージID.スキャフォールディング省略_コンティグなし);
                return;
            }

            var l_配置 = p_コンティグ構築.A_ユニティグ配置;
            var l_ペア経路 = p_コンティグ構築.A_ペア経路;

            // contig 単位の頂点空間を作る。unitig 同様、各 contig を
            // 「順方向」「逆方向」の2頂点として扱う。
            // 頂点番号 = コンティグID << 1 (順方向) / コンティグID << 1 | 1 (逆方向)
            var l_コンティグ数 = this._コンティグ配列.Keys.Count == 0 ? 0 : this._コンティグ配列.Keys.Max();
            var l_頂点数 = (l_コンティグ数 + 1) << 1;

            var l_隣接 = new List<スキャフォールド候補>[l_頂点数];
            for (var i = 0; i < l_頂点数; i++)
            {
                l_隣接[i] = [];
            }

            var l_辺の集計 = new Dictionary<(int, int), (ulong A_支持数, List<int> A_既知長標本)>();

            var l_内部を指した数 = 0;
            var l_未配置を指した数 = 0;

            foreach (var (l_キー, l_標本) in l_ペア経路)
            {
                var (l_始点ユニティグ, l_終点ユニティグ) = l_キー;

                if (!Get_コンティグ末端頂点(l_配置, l_始点ユニティグ, p_出口側か: true, out var l_始点頂点))
                {
                    if (!l_配置.ContainsKey(Math.Abs(l_始点ユニティグ)))
                    {
                        l_未配置を指した数++;
                    }
                    else
                    {
                        l_内部を指した数++;
                    }
                    continue;
                }

                if (!Get_コンティグ末端頂点(l_配置, l_終点ユニティグ, p_出口側か: false, out var l_終点頂点))
                {
                    if (!l_配置.ContainsKey(Math.Abs(l_終点ユニティグ)))
                    {
                        l_未配置を指した数++;
                    }
                    else
                    {
                        l_内部を指した数++;
                    }
                    continue;
                }

                // 自己ループ(同一 contig の同一末端同士)は無視する。
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

            if (l_内部を指した数 > 0)
            {
                Logger.V_出力(メッセージID.内部を指したペア候補, l_内部を指した数);
            }
            if (l_未配置を指した数 > 0)
            {
                Logger.V_出力(メッセージID.未配置を指したペア候補, l_未配置を指した数);
            }

            // v→w と双子 w^1→v^1 は同一の隣接だが、ペアエンドの観測は
            // 片方の向きにしか記録されない。対称化しないと逆鎖側の支持がゼロになり、
            // 相互一意性の検査が常に落ちる。各観測は一方のキーにしか入っていないので
            // 和を取っても二重計上にはならない。
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

            var l_モデル = new PairedDistanceModel(p_コンティグ構築.A_同一ユニティグ標本, p_リード長 ?? l_インサートサイズ);
            var l_較正器 = 証拠較正器.Get_較正器(
                p_コンティグ構築.A_同一ユニティグ標本, p_リード長 ?? l_インサートサイズ,
                p_コンティグ構築.A_ユニティグ長.Values.Select(x => (long)x));

            foreach (var ((l_始点, l_終点), (_, l_標本)) in l_対称化)
            {
                var (l_一貫した本数, l_ギャップ長) = l_モデル.Get_一貫した支持(l_標本);

                // 期待は接合点から1フラグメント長ぶんの窓しか効かないので、
                // 重なっている(ギャップが負)場合は接している場合と同じとみなす。
                var l_期待に対する比 = l_較正器.Get_正規化済み支持(
                    (ulong)l_一貫した本数, this.Get_コンティグ長(l_始点), this.Get_コンティグ長(l_終点), Math.Max(0, l_ギャップ長));

                l_隣接[l_始点].Add(new スキャフォールド候補(l_終点, (ulong)l_一貫した本数, l_ギャップ長, l_期待に対する比));
            }

            var l_優勢閾値 = ConfigurationManager.A_実行時引数.A_ペア結合閾値;
            var l_最小証拠数 = Consts.スキャフォールド支持数の下限;

            Logger.V_出力(
                メッセージID.スキャフォールド候補辺数,
                l_辺の集計.Count,
                Messages.Get_文言(l_較正器.A_使えるか
                    ? メッセージID.理想本数モデルあり
                    : メッセージID.理想本数モデルなし));

            // 各頂点について、最多支持の辺1本だけを残す。
            var l_確定辺 = new (int A_行き先, int A_ギャップ長)?[l_頂点数];
            for (var v = 2; v < l_頂点数; v++)
            {
                V_確定_スキャフォールド辺(l_隣接, v, l_優勢閾値, l_最小証拠数, l_確定辺);
            }

            var l_確定数 = 0;
            for (var v = 2; v < l_頂点数; v++)
            {
                if (l_確定辺[v] != null)
                {
                    l_確定数++;
                }
            }

            // 相互一意な辺だけを採用する。v→w を繋いでよいのは
            // 「v の唯一の行き先が w」であり、かつ「w の唯一の来訪元が v」で
            // あるときに限る。後者は逆鎖対称性より 確定辺[w^1] が v^1 を
            // 指すことと同値。これを課さないと、複数の contig が同じ次の contig を
            // 指した場合に先着1本だけが繋がれ、残りは黙って千切れる
            // (どれが正しいかの根拠がないまま1本を選ぶことになる)。
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
                    AmbiguityRecorder.V_記録(
                        曖昧箇所の種別.経路が一意でない,
                        AmbiguityRecorder.Get_場所名(v, "contig"));
                }
            }
            Logger.V_出力(メッセージID.閾値後のスキャフォールド辺, l_確定数, l_相互一意で棄却した数, l_確定数 - l_相互一意で棄却した数);

            // 「入ってくる結合を持たない」頂点が経路の始点。v への結合が
            // 存在することは、逆鎖対称性より 確定辺[v^1] != null と同値。
            var l_始点群 = new List<int>();
            for (var v = 2; v < l_頂点数; v++)
            {
                if (this._コンティグ配列.ContainsKey(v >> 1) && (v ^ 1) < l_頂点数 && l_確定辺[v ^ 1] == null)
                {
                    l_始点群.Add(v);
                }
            }

            // 配列と、それが環状に閉じた複製単位そのものかどうか。
            // 単独の contig がそのまま1本のスキャフォールドになった場合だけ
            // 環状を引き継ぐ。他の contig を継ぎ足した時点で、それはもう
            // 閉じた環ではない。
            List<(string A_配列, bool A_環状か)> l_スキャフォールド群 = [];
            var l_訪問済み = new bool[l_頂点数];
            foreach (var l_始点 in l_始点群)
            {
                // 始点群には同一 contig の順鎖/逆鎖の両方の頂点が独立に
                // 含まれうる。先に処理された方の walk が両方向を訪問済みに
                // するため、後から来た方はここでスキップしないと、同じ contig を
                // 起点とするスキャフォールドが二重に生成されてしまう
                // (contig 数の水増し・配列の重複の原因)。
                if (l_訪問済み[l_始点])
                {
                    continue;
                }
                var l_スキャフォールド = this.Get_スキャフォールド配列(
                    l_確定辺, l_始点, l_訪問済み, out var l_連結数);
                if (l_スキャフォールド != null)
                {
                    l_スキャフォールド群.Add((
                        l_スキャフォールド,
                        l_連結数 == 1 && this.Get_環状か(l_始点 >> 1)));
                }
            }

            // まだ訪問されていない(=孤立した、あるいは循環に巻き込まれた)contig を
            // 単独スキャフォールドとして出力する。
            for (var l_コンティグID = 1; l_コンティグID <= l_コンティグ数; l_コンティグID++)
            {
                var l_順鎖 = l_コンティグID << 1;
                var l_逆鎖 = (l_コンティグID << 1) | 1;
                if (l_順鎖 < l_頂点数 && !l_訪問済み[l_順鎖] && !l_訪問済み[l_逆鎖]
                    && this._コンティグ配列.TryGetValue(l_コンティグID, out var l_配列))
                {
                    l_スキャフォールド群.Add((l_配列, this.Get_環状か(l_コンティグID)));
                    l_訪問済み[l_順鎖] = true;
                    l_訪問済み[l_逆鎖] = true;
                }
            }

            using var l_書き込み = new FastaWriter(p_スキャフォールドパス);
            var l_スキャフォールドID = 1;
            long l_総延長 = 0;
            foreach (var (l_配列, l_環状か) in l_スキャフォールド群)
            {
                var l_名前 = l_環状か
                    ? $"SCAFFOLD{l_スキャフォールドID}_{Consts.環状の目印}"
                    : $"SCAFFOLD{l_スキャフォールドID}";
                l_書き込み.V_書き込み(l_名前, l_配列);
                l_スキャフォールドID++;
                l_総延長 += l_配列.Length;
            }

            Logger.V_出力(メッセージID.スキャフォールド出力完了, l_スキャフォールド群.Count, l_総延長, p_スキャフォールドパス);
        }

        /// <summary>
        /// インサートサイズを確定する。明示指定があればそれを使う。
        ///
        /// 未指定なら2種類の標本群から選ぶ。同一 unitig 内標本は打ち切りバイアスを
        /// 持つが unitig が十分長ければ起きず、標本数が桁違いに多い。
        /// 確定辺由来は unitig 長に縛られない代わりに標本数が極端に少なく、
        /// 誤結合や誤マッピングの影響を受けやすい。
        /// </summary>
        private bool Get_インサートサイズ(out int p_インサートサイズ)
        {
            if (ConfigurationManager.A_実行時引数.A_インサートサイズ is { } l_指定値)
            {
                p_インサートサイズ = l_指定値;
                return true;
            }

            var l_同一ユニティグ標本 = p_コンティグ構築.A_同一ユニティグ標本;
            if (l_同一ユニティグ標本.Count >= Consts.インサートサイズ標本数の下限)
            {
                var l_推定値 = StatsUtil.Get_中央値(l_同一ユニティグ標本);
                var l_ユニティグN50 = Get_ユニティグN50(p_コンティグ構築.A_ユニティグ長);
                if (l_推定値 > 0 && l_ユニティグN50 >= (long)l_推定値 * 偏りが無いとみなす長さ比)
                {
                    p_インサートサイズ = l_推定値;
                    Logger.V_出力(メッセージID.インサートサイズ推定_同一ユニティグ, p_インサートサイズ, l_同一ユニティグ標本.Count, l_ユニティグN50, 偏りが無いとみなす長さ比);
                    return true;
                }
            }

            var l_確定辺標本 = p_コンティグ構築.A_確定辺標本;
            if (l_確定辺標本.Count >= Consts.インサートサイズ標本数の下限)
            {
                p_インサートサイズ = StatsUtil.Get_中央値(l_確定辺標本);
                Logger.V_出力(メッセージID.インサートサイズ推定_確定辺, p_インサートサイズ, l_確定辺標本.Count);
                return true;
            }

            var l_全標本 = p_コンティグ構築.A_インサートサイズ標本;
            if (l_全標本.Count < Consts.インサートサイズ標本数の下限)
            {
                Logger.V_出力(メッセージID.インサートサイズ推定_標本不足, Consts.インサートサイズ標本数の下限, l_確定辺標本.Count, l_全標本.Count);
                p_インサートサイズ = 0;
                return false;
            }

            p_インサートサイズ = StatsUtil.Get_中央値(l_全標本);
            Logger.V_出力(メッセージID.インサートサイズ推定_全標本, p_インサートサイズ, l_全標本.Count, l_確定辺標本.Count);
            return true;
        }

        /// <summary>
        /// unitig の N50。打ち切りバイアスの有無の判断に使う。
        /// 平均ではなく N50 を使うのは、本数では短い断片が多くても
        /// ペアが実際に観測される場所は長い unitig に偏るため。
        /// </summary>
        private static long Get_ユニティグN50(IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            return StatsUtil.Get_N50([.. p_ユニティグ長.Values.Select(x => (long)x)]).A_N50;
        }

        private void V_読込_コンティグ()
        {
            using var l_読み込み = new FastaReader(p_コンティグファイルパス);
            var l_ID = 1;
            while (l_読み込み.Get_続きがあるか())
            {
                var l_配列エントリ = l_読み込み.Get_次の配列();
                this._コンティグ名[l_ID] = l_配列エントリ.A_ID.TrimStart('>');
                this._コンティグ配列[l_ID] = l_配列エントリ.A_配列;
                l_ID++;
            }
        }

        /// <summary>
        /// 符号付き unitig ID が contig の末端に配置されているかを判定し、
        /// 配置されていれば対応する contig 頂点を返す。
        ///
        /// 出口側(読み進める起点)として有効なのは「順鎖かつ contig 内で末尾」
        /// または「逆鎖かつ先頭」、入口側はその逆。
        /// contig が正規化で逆相補化されていると walk 順の先頭/末尾の意味が
        /// 反転するため、その分も考慮して向きを決める。
        /// </summary>
        private static bool Get_コンティグ末端頂点(
            IReadOnlyDictionary<int, ユニティグ配置> p_配置,
            int p_符号付きユニティグID,
            bool p_出口側か,
            out int p_頂点番号)
        {
            p_頂点番号 = 0;
            var l_ユニティグID = Math.Abs(p_符号付きユニティグID);
            var l_順鎖か = p_符号付きユニティグID > 0;

            if (!p_配置.TryGetValue(l_ユニティグID, out var l_配置情報))
            {
                return false;
            }

            // unitig 自身が walk 中に逆鎖として使われていた場合、ペア経路上の
            // 向きは「unitig 単体の元の向き」を基準にしているため、
            // walk 内での実効的な向きに変換する。
            var l_実効的に順鎖か = l_順鎖か != l_配置情報.A_walk中で逆鎖か;

            var l_該当する端にあるか = p_出口側か
                ? l_実効的に順鎖か ? l_配置情報.A_コンティグ末尾か : l_配置情報.A_コンティグ先頭か
                : l_実効的に順鎖か ? l_配置情報.A_コンティグ先頭か : l_配置情報.A_コンティグ末尾か;
            if (!l_該当する端にあるか)
            {
                return false;
            }

            // contig 全体が正規化のために逆相補化されている場合、
            // 「walk 順で見た先頭/末尾」と「実際の contigs.fasta 上の先頭/末尾」が
            // 入れ替わる。スキャフォールディングは contigs.fasta 上の配列
            // (=実際に出力された向き)を基準に扱うため、ここで反転させる。
            var l_最終配列で順鎖か = l_配置情報.A_コンティグが逆相補か ? !l_実効的に順鎖か : l_実効的に順鎖か;

            p_頂点番号 = (l_配置情報.A_コンティグID << 1) | (l_最終配列で順鎖か ? 0 : 1);
            return true;
        }

        private static void V_確定_スキャフォールド辺(
            List<スキャフォールド候補>[] p_隣接,
            int p_頂点,
            decimal p_優勢閾値,
            ulong p_最小証拠数,
            (int A_行き先, int A_ギャップ長)?[] p_確定辺)
        {
            var l_最良 = Get_優勢な候補(p_隣接[p_頂点], p_優勢閾値, p_最小証拠数);
            if (l_最良 is not { } l_辺)
            {
                p_確定辺[p_頂点] = null;
                return;
            }

            p_確定辺[p_頂点] = (l_辺.A_行き先, Math.Max(Consts.ギャップ長の下限, l_辺.A_ギャップ長));
        }

        /// <summary>
        /// 支持数と期待本数比の下限を満たし、その中で優勢比を超える辺を返す。
        ///
        /// 期待本数と比べるのは、辺が長く距離が近いほど多く観測されるという
        /// 幾何的な偏りを外すため。観測本数だけを固定の下限と比べると、
        /// 期待が数本の場所と数百本の場所を同じ物差しで測ることになる。
        /// </summary>
        internal static スキャフォールド候補? Get_優勢な候補(
            IReadOnlyList<スキャフォールド候補> p_候補,
            decimal p_優勢閾値,
            ulong p_最小証拠数)
        {
            var l_候補 = p_候補.Where(x => x.A_支持数 >= p_最小証拠数).ToList();
            if (l_候補.Count == 0)
            {
                return null;
            }

            var l_合計 = l_候補.Sum(x => x.A_期待に対する比);
            var l_最良 = l_候補.OrderByDescending(x => x.A_期待に対する比).First();
            if (l_合計 <= 0 || (decimal)(l_最良.A_期待に対する比 / l_合計) < p_優勢閾値)
            {
                return null;
            }

            return l_最良;
        }

        private long Get_コンティグ長(int p_頂点)
        {
            return this._コンティグ配列.TryGetValue(p_頂点 >> 1, out var l_配列) ? l_配列.Length : 0;
        }

        /// <summary>
        /// 標本群から挿入する N の数を決める。各標本は既知長で、
        /// フラグメント長 = 標本 + ギャップ長 が成り立つため、
        /// ギャップ長 = インサートサイズ - 標本 の中央値を採る。
        /// 推定が負や 0 でも隣接の事実自体には証拠があるので、下限で丸めて
        /// 少なくとも1つの N を残す。
        /// </summary>
        private int Get_推定ギャップ長(List<int> p_既知長標本)
        {
            var l_インサートサイズ = this.A_有効インサートサイズ ?? 0;
            if (p_既知長標本.Count == 0)
            {
                return Consts.ギャップ長の下限;
            }

            var l_ギャップ候補 = p_既知長標本.Select(x => l_インサートサイズ - x).ToList();
            return Math.Max(Consts.ギャップ長の下限, StatsUtil.Get_中央値(l_ギャップ候補));
        }

        /// <summary>
        /// contig を1本のスキャフォールドへ連ねる。p_連結したコンティグ数 は
        /// 実際に繋いだ本数で、1 なら元の contig がそのまま出ていることを意味する。
        /// </summary>
        private string? Get_スキャフォールド配列(
            (int A_行き先, int A_ギャップ長)?[] p_確定辺, int p_始点, bool[] p_訪問済み,
            out int p_連結したコンティグ数)
        {
            p_連結したコンティグ数 = 0;
            var l_コンティグID = p_始点 >> 1;
            var l_逆鎖か = (p_始点 & 1) == 1;
            if (!this._コンティグ配列.TryGetValue(l_コンティグID, out var l_配列))
            {
                return null;
            }
            p_連結したコンティグ数 = 1;

            var l_出力 = new StringBuilder(l_逆鎖か ? Util.V_逆相補(l_配列) : l_配列);
            var l_現在 = p_始点;
            // 頂点を「消費」した(=いずれかの向きでスキャフォールドに組み込んだ)際は、
            // その contig の両方の向きの頂点を訪問済みにする。
            // 片方の頂点だけを訪問済みにすると、同じ contig の反対向きの頂点が
            // 別の開始点や「未訪問の孤立 contig」判定で再度使われてしまう
            // (同じ contig が2回出力される)おそれがあるため。
            V_記録_訪問済み(p_訪問済み, l_現在);
            while (p_確定辺[l_現在] is { } l_辺 && !p_訪問済み[l_辺.A_行き先])
            {
                var l_次のコンティグID = l_辺.A_行き先 >> 1;
                var l_次が逆鎖か = (l_辺.A_行き先 & 1) == 1;
                if (!this._コンティグ配列.TryGetValue(l_次のコンティグID, out var l_次の配列))
                {
                    break;
                }

                _ = l_出力.Append('N', l_辺.A_ギャップ長);
                _ = l_出力.Append(l_次が逆鎖か ? Util.V_逆相補(l_次の配列) : l_次の配列);

                l_現在 = l_辺.A_行き先;
                p_連結したコンティグ数++;
                V_記録_訪問済み(p_訪問済み, l_現在);
            }

            return l_出力.ToString();
        }

        /// <summary>その contig が環状に閉じたものとして作られたか。</summary>
        private bool Get_環状か(int p_コンティグID)
        {
            return this._コンティグ名.TryGetValue(p_コンティグID, out var l_名前)
                && l_名前.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// 頂点番号が指す contig の両方の向きの頂点を訪問済みにする。
        /// </summary>
        private static void V_記録_訪問済み(bool[] p_訪問済み, int p_頂点番号)
        {
            var l_コンティグID = p_頂点番号 >> 1;
            var l_順鎖 = l_コンティグID << 1;
            var l_逆鎖 = l_順鎖 | 1;
            if (l_逆鎖 < p_訪問済み.Length)
            {
                p_訪問済み[l_順鎖] = true;
                p_訪問済み[l_逆鎖] = true;
            }
        }
    }
}
