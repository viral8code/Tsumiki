using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Evidence;
using Tsumiki.Cores.Output;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.IO;
using Tsumiki.Models.ContigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Reporting;
using Tsumiki.Utilities;

namespace Tsumiki.Core
{
    /// <summary>
    /// unitig を辺で結合して contig を組み立てる
    /// </summary>
    internal partial class ContigMaker
    {
        #region 定数

        /// <summary>
        /// 簡略化ラウンド上限
        /// </summary>
        private const int ラウンド数上限 = 5;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// unitig グラフから辺を選び結合を確定して、contig を FASTA へ書き出す
        /// </summary>
        /// <param name="p_コンティグパス">出力先の FASTA パス</param>
        /// <param name="p_優勢閾値">分岐選択で優勢とみなす正規化支持の割合</param>
        /// <param name="p_最小証拠数">分岐選択に必要な最小の証拠数</param>
        /// <param name="p_コピー数">
        /// unitig ID -> 推定コピー数<br/>
        /// 先読み探索で「この unitig を何回まで通ってよいか」の予算に使う<br/>
        /// 渡さない場合はすべて 1 コピーとして扱い、先読み探索も控えめになる
        /// </param>
        /// <param name="p_バブル敗者への引き継ぎ先">
        /// 渡すと、バブル除去で外れた側の経路の配列 (careful_bubble) をここへ集める<br/>
        /// 呼び出し側がマルチ k の次の k への引き継ぎに足すことを想定している
        /// </param>
        /// <param name="p_リード長">
        /// 分岐選択・先読みスコアを生カウントではなく期待本数との比で測るための較正器の構築に使う<br/>
        /// 渡さない場合は生カウント方式になる
        /// </param>
        /// <param name="p_r_mer検証器">
        /// 渡すと、短い反復解決の対応付けを r-mer で検証する拒否権 (ABySS RResolver 型) を課す<br/>
        /// 詳細は UnitigGraph.V_解決_短い反復 を参照
        /// </param>
        /// <param name="p_GFAパス">
        /// 渡すと、バブル除去・反復解決を終えたあとの unitig グラフを GFA1 形式でこのパスへ書き出す (Bandage 等のビューア向け)
        /// </param>
        public void V_結合_コンティグ(string p_コンティグパス, decimal p_優勢閾値, ulong p_最小証拠数, IReadOnlyDictionary<int, int>? p_コピー数 = null, List<string>? p_バブル敗者への引き継ぎ先 = null, int? p_リード長 = null, RepeatRMerVerifier? p_r_mer検証器 = null, string? p_GFAパス = null)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_重なり長 = l_k長 - 1;

            var l_ユニティグ配列 = this._ユニティグ配列;

            // 隣接は de Bruijn グラフから厳密に導く (UnitigGraph の説明を参照)
            // リードマッピング由来の隣接情報は「辺を作る」ためではなく、
            // 分岐点でどの辺を選ぶかの「重み」としてのみ使う
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ配列, this._kmer辞書, l_k長, 曖昧kmerの番兵);

            var l_辺数 = 0;
            var l_分岐頂点数 = 0;
            for (var v = 2; v < l_グラフ.A_出辺.Count; v++)
            {
                l_辺数 += l_グラフ.A_出辺[v].Count;
                if (l_グラフ.A_出辺[v].Count > 1)
                {
                    l_分岐頂点数++;
                }
            }
            Logger.V_出力(メッセージID.デブルーイングラフの要約, l_辺数, l_分岐頂点数, l_グラフ.A_出辺.Count - 2);

            var (l_支持, l_ペア連結) = this.Get_辺重み(l_グラフ);

            // 跨げる見込みのある長さの上限
            // フラグメント長の実測中央値を使う
            // (これより長い反復は、そもそも両端を別々の unitig に載せたペアが存在しえない)
            // 標本が無い場合は控えめな既定値
            var l_反復長の上限 = this.A_同一ユニティグ標本.Count > 0 ? StatsUtil.Get_中央値(this.A_同一ユニティグ標本) : l_k長 * 4;

            V_簡略化ラウンド(l_グラフ, l_ユニティグ配列, l_支持, l_ペア連結, l_反復長の上限, p_優勢閾値, p_最小証拠数, p_r_mer検証器, p_バブル敗者への引き継ぎ先);

            // 支持を生カウントではなく期待本数との比で測るための較正器
            // 短い辺には厳しすぎ、長い辺には緩すぎる固定閾値のバイアスを外す
            // (較正器が使えない場合は生カウントへフォールバックし、挙動は従来と完全に一致する)
            var l_較正器 = 証拠較正器.Get_較正器(this.A_同一ユニティグ標本, p_リード長, this._ユニティグ長.Values.Select(x => (long)x));

            var l_選択 = this.Get_辺選択(l_グラフ, l_支持, l_較正器, p_コピー数, p_優勢閾値, p_最小証拠数);
            var l_結合 = Get_結合確定(l_グラフ, l_選択, p_コピー数);

            // 1 歩だけを見る相互一意性の判定では決めきれなかった分岐を、
            // 数 kb 先まで複数経路を並行して伸ばして (ビームサーチ) 解けるだけ解く
            // 分岐の直後だけを見ると五分五分でも、少し先まで進めると片方だけが
            // ペアエンドの証拠と整合する、という状況を拾える
            var l_先読みで解決した数 = BeamSearchExtender.V_延長_先読み(l_グラフ, l_ユニティグ配列, l_結合, l_ペア連結, p_コピー数 ?? new Dictionary<int, int>(), l_反復長の上限, p_優勢閾値, p_最小証拠数, l_較正器);

            if (l_先読みで解決した数 > 0)
            {
                Logger.V_出力(メッセージID.先読みで解決した分岐数, l_先読みで解決した数 >> 1);
            }

            if (p_GFAパス is not null)
            {
                GfaWriter.V_出力(p_GFAパス, l_ユニティグ配列, l_グラフ, l_k長, p_コピー数);
                Logger.V_出力(メッセージID.GFA出力完了, p_GFAパス);
            }

            this.V_収集_確定辺標本(l_結合);

            this.V_walk実行してFASTA書き出し(l_グラフ, l_ユニティグ配列, l_結合, l_重なり長, p_コンティグパス);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// リード隣接・ペア経路から、辺選択に使う支持数 (逆鎖対称に集計) と反復解決に使うペア連結を組み立てる
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <returns></returns>
        private (Dictionary<(int, int), ulong> A_支持, Dictionary<(int, int), ulong> A_ペア連結) Get_辺重み(UnitigGraph p_グラフ)
        {
            // リード由来の支持数を逆鎖対称に集計する
            // 辺 v→w と w^1→v^1 は
            // 同一の物理的な隣接を表すため、重みも同一でなければ順鎖側と
            // 逆鎖側で異なる経路が選ばれ、同じ領域が 2 通りに組み立てられてしまう
            Dictionary<(int, int), ulong> l_支持 = [];
            foreach (var ((l_始点, l_終点), l_件数) in this._リード隣接)
            {
                if (l_始点 == l_終点)
                {
                    continue;
                }
                var l_始点番号 = Get_頂点番号(l_始点);
                var l_終点番号 = Get_頂点番号(l_終点);
                l_支持[(l_始点番号, l_終点番号)] = l_支持.GetValueOrDefault((l_始点番号, l_終点番号)) + l_件数;
                l_支持[(l_終点番号 ^ 1, l_始点番号 ^ 1)] = l_支持.GetValueOrDefault((l_終点番号 ^ 1, l_始点番号 ^ 1)) + l_件数;
            }

            // 1 本のリードでは跨げない長さの反復も、フラグメント長なら跨げる
            // 参照するのはグラフ上に実在する辺だけなので、隣接していない
            // unitig 対がペア経路に入っていても選択には影響しない
            Dictionary<(int, int), ulong> l_ペア連結 = [];
            var l_ペア支持を足した数 = 0;
            foreach (var ((l_始点, l_終点), l_標本) in this._ペア経路)
            {
                if (l_始点 == l_終点)
                {
                    continue;
                }
                var l_始点番号 = Get_頂点番号(l_始点);
                var l_終点番号 = Get_頂点番号(l_終点);
                var l_件数 = (ulong)l_標本.Count;
                l_支持[(l_始点番号, l_終点番号)] = l_支持.GetValueOrDefault((l_始点番号, l_終点番号)) + l_件数;
                l_支持[(l_終点番号 ^ 1, l_始点番号 ^ 1)] = l_支持.GetValueOrDefault((l_終点番号 ^ 1, l_始点番号 ^ 1)) + l_件数;
                l_ペア連結[(l_始点番号, l_終点番号)] = l_ペア連結.GetValueOrDefault((l_始点番号, l_終点番号)) + l_件数;
                l_ペア連結[(l_終点番号 ^ 1, l_始点番号 ^ 1)] = l_ペア連結.GetValueOrDefault((l_終点番号 ^ 1, l_始点番号 ^ 1)) + l_件数;
                l_ペア支持を足した数++;
            }
            Logger.V_出力(メッセージID.分岐選択の重み内訳, this._リード隣接.Count, l_ペア支持を足した数);

            return (l_支持, l_ペア連結);
        }

        /// <summary>
        /// 単純バブルを潰してから辺を選ぶ
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_ユニティグ配列"></param>
        /// <param name="p_支持"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_反復長の上限"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <param name="p_r_mer検証器"></param>
        /// <param name="p_バブル敗者への引き継ぎ先"></param>
        private static void V_簡略化ラウンド(UnitigGraph p_グラフ, List<string> p_ユニティグ配列, Dictionary<(int, int), ulong> p_支持, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, int p_反復長の上限, decimal p_優勢閾値, ulong p_最小証拠数, RepeatRMerVerifier? p_r_mer検証器, List<string>? p_バブル敗者への引き継ぎ先)
        {
            var l_除去バブル数 = 0;
            var l_解決した反復数 = 0;
            for (var l_ラウンド = 1; l_ラウンド <= ラウンド数上限; l_ラウンド++)
            {
                var l_今回のバブル数 = p_グラフ.V_除去_単純バブル(p_ユニティグ配列, p_支持, ConfigurationManager.A_実行時引数.A_k長, p_バブル敗者への引き継ぎ先);
                var l_今回の反復数 = p_グラフ.V_解決_短い反復(p_ユニティグ配列, p_支持, p_ペア連結, p_反復長の上限, p_優勢閾値, p_最小証拠数, p_r_mer検証器);
                l_除去バブル数 += l_今回のバブル数;
                l_解決した反復数 += l_今回の反復数;

                if (l_今回のバブル数 == 0 && l_今回の反復数 == 0)
                {
                    Logger.V_出力(メッセージID.単純化の収束, l_ラウンド);
                    break;
                }

                if (l_ラウンド == ラウンド数上限)
                {
                    Logger.V_出力(メッセージID.単純化の打ち切り, ラウンド数上限);
                }
            }
            if (l_除去バブル数 > 0)
            {
                Logger.V_出力(メッセージID.バブル除去数, l_除去バブル数);
            }
            Logger.V_出力(メッセージID.反復解決数, l_解決した反復数, p_反復長の上限);
        }

        /// <summary>
        /// 各頂点について「出て行く先」を高々 1 つに絞る
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_支持"></param>
        /// <param name="p_較正器"></param>
        /// <param name="p_コピー数"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <returns></returns>
        private int[] Get_辺選択(UnitigGraph p_グラフ, Dictionary<(int, int), ulong> p_支持, 証拠較正器 p_較正器, IReadOnlyDictionary<int, int>? p_コピー数, decimal p_優勢閾値, ulong p_最小証拠数)
        {
            var l_選択 = new int[p_グラフ.A_出辺.Count];
            Array.Fill(l_選択, -1);
            var l_一意な頂点数 = 0;
            var l_支持で解決した数 = 0;
            var l_反復由来で未解決の数 = 0;
            for (var v = 2; v < p_グラフ.A_出辺.Count; v++)
            {
                var l_出辺 = p_グラフ.A_出辺[v];

                if (l_出辺.Count == 0)
                {
                    continue;
                }

                if (l_出辺.Count == 1)
                {
                    l_選択[v] = l_出辺[0];
                    l_一意な頂点数++;
                    continue;
                }

                // 分岐元が多コピーだと、リード支持はどのコピー由来か区別できない
                // 各コピーの続きが 1 頂点に集まるため支持が全ての行き先に付き、
                // 行き先を選ぶ根拠にならない
                // 正しい続きは反復の外側の
                // 単一コピー領域からのペアエンドでしか決められない
                if (p_コピー数 is not null && p_コピー数.GetValueOrDefault(v >> 1, 1) > 1)
                {
                    l_反復由来で未解決の数++;
                    AmbiguityRecorder.V_記録(曖昧箇所の種別.反復の内側, AmbiguityRecorder.Get_場所名(v));
                    continue;
                }

                var l_始点長 = this._ユニティグ長.GetValueOrDefault(v >> 1, 0);
                var l_最良 = -1;
                var l_最良の生本数 = 0UL;
                var l_最良の正規化 = double.NegativeInfinity;
                var l_正規化合計 = 0D;
                foreach (var w in l_出辺)
                {
                    var l_件数 = p_支持.GetValueOrDefault((v, w));
                    var l_終点長 = this._ユニティグ長.GetValueOrDefault(w >> 1, 0);

                    // 較正器が使えない場合は生カウントをそのまま正規化値として扱う
                    // これにより以下の判定式は較正器が無かった従来のロジックと
                    // 完全に同じ結果になる
                    var l_正規化 = p_較正器.A_Is使用可能 ? p_較正器.Get_正規化済み支持(l_件数, l_始点長, l_終点長, p_ギャップ長: 0) : l_件数;
                    l_正規化合計 += l_正規化;
                    if (l_正規化 > l_最良の正規化)
                    {
                        l_最良の正規化 = l_正規化;
                        l_最良の生本数 = l_件数;
                        l_最良 = w;
                    }
                }

                if (l_最良 >= 0 && l_最良の生本数 >= p_最小証拠数 && l_正規化合計 > 0D && (decimal)(l_最良の正規化 / l_正規化合計) >= p_優勢閾値)
                {
                    l_選択[v] = l_最良;
                    l_支持で解決した数++;
                }
                else
                {
                    // 支持が足りないのか、上位が割れているのかで意味が違う
                    var l_種別 = l_最良の生本数 < p_最小証拠数 ? 曖昧箇所の種別.支持なし : 曖昧箇所の種別.僅差;
                    var l_場所 = AmbiguityRecorder.Get_場所名(v);
                    var l_首位の支持 = l_最良の正規化 is double.NegativeInfinity ? 0D : l_最良の正規化;
                    var l_次点の支持 = l_正規化合計 - (l_最良の正規化 is double.NegativeInfinity ? 0D : l_最良の正規化);
                    AmbiguityRecorder.V_記録(l_種別, l_場所, l_首位の支持, l_次点の支持, (long)l_最良の生本数);
                }
            }
            Logger.V_出力(メッセージID.辺選択の内訳, l_一意な頂点数, l_支持で解決した数, l_反復由来で未解決の数);

            return l_選択;
        }

        /// <summary>
        /// v→w を結合してよいのは v の唯一の行き先が w で、かつ w の唯一の来訪元が v のときだけ (後者は逆鎖対称性より 選択[w^1] == v^1)
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_選択"></param>
        /// <param name="p_コピー数"></param>
        /// <remarks>
        /// これを欠くと、同じ行き先を指す複数の unitig のうち先着だけが 結合され、残りが根拠なく千切れる
        /// </remarks>
        /// <returns></returns>
        private static int[] Get_結合確定(UnitigGraph p_グラフ, int[] p_選択, IReadOnlyDictionary<int, int>? p_コピー数)
        {
            var l_結合 = new int[p_グラフ.A_出辺.Count];
            Array.Fill(l_結合, -1);
            var l_結合数 = 0;
            var l_反復通り抜けで棄却した数 = 0;
            for (var v = 2; v < p_グラフ.A_出辺.Count; v++)
            {
                var l_終点 = p_選択[v];
                if (l_終点 < 0 || p_選択[l_終点 ^ 1] != (v ^ 1))
                {
                    continue;
                }

                // 結合は逆鎖側と対で成立する
                // 片側だけ許すと結合の対称性が
                // 崩れ、walk の始点判定が壊れるため、
                // どちらかが通り抜け不可なら対ごと採用しない
                if (!p_グラフ.Is通過可能(p_コピー数, v) || !p_グラフ.Is通過可能(p_コピー数, l_終点 ^ 1))
                {
                    l_反復通り抜けで棄却した数++;
                    continue;
                }
                l_結合[v] = l_終点;
                l_結合数++;
            }

            // 分岐を 1 つも持たない環状の複製単位は、unitig 1 本の末尾が自分の
            // 先頭へ戻るだけの形になる
            // この自己ループはグラフの辺として持てないため、
            // 他に出入りが無い孤立した頂点に限って結合として明示し、
            // walk に環として閉じさせる (閉じないと重なりの k-1 塩基が余分に残り、環状であることも分からないまま線状の断片になる)
            var l_孤立した環の数 = 0;
            foreach (var l_始点 in p_グラフ.A_自己ループ)
            {
                if (l_結合[l_始点] != -1 || p_グラフ.A_出辺[l_始点].Count != 0 || p_グラフ.Get_入次数(l_始点) != 0)
                {
                    continue;
                }
                l_結合[l_始点] = l_始点;
                l_孤立した環の数++;
            }

            if (l_孤立した環の数 > 0)
            {
                Logger.V_出力(メッセージID.孤立した環状の複製単位, l_孤立した環の数 >> 1);
            }

            if (l_反復通り抜けで棄却した数 > 0)
            {
                Logger.V_出力(メッセージID.反復通り抜けで棄却した結合数, l_反復通り抜けで棄却した数 >> 1);
            }
            Logger.V_出力(メッセージID.相互一意で残った結合数, l_結合数, l_結合数 >> 1);

            return l_結合;
        }

        /// <summary>
        /// 確定した結合を辿って contig を組み立て、FASTA へ書き出す
        /// </summary>
        /// <param name="p_グラフ">unitig グラフ</param>
        /// <param name="p_ユニティグ配列">unitig ID 順の配列</param>
        /// <param name="p_結合">頂点ごとの結合先</param>
        /// <param name="p_重なり長">隣り合う unitig が共有する長さ</param>
        /// <param name="p_コンティグパス">書き出し先</param>
        private void V_walk実行してFASTA書き出し(UnitigGraph p_グラフ, List<string> p_ユニティグ配列, int[] p_結合, int p_重なり長, string p_コンティグパス)
        {
            // 双子 (v と v^1) は同一 unitig の裏表なので、unitig 単位で訪問済みを管理する
            // これを頂点単位でやっていたため、順鎖側の walk と逆鎖側の
            // walk が同じ unitig を別々に出力し、contig 総長が unitig 総長の
            // ちょうど 2 倍に膨れていた
            var l_ユニティグ数 = (p_ユニティグ配列.Count - 2) >> 1;
            var l_訪問済み = new bool[l_ユニティグ数 + 1];

            List<string> l_コンティグ群 = [];
            List<List<int>> l_walk順群 = [];
            List<bool> l_環状フラグ群 = [];

            // 結合グラフ上で「入ってくる結合を持たない」頂点が経路の始点
            // v への結合が存在することは、逆鎖対称性より 結合[v^1] != -1 と同値
            for (var v = 2; v < p_グラフ.A_出辺.Count; v++)
            {
                if (p_結合[v ^ 1] != -1 || l_訪問済み[v >> 1])
                {
                    continue;
                }
                V_実行_walk(p_ユニティグ配列, p_結合, l_訪問済み, p_重なり長, v, l_コンティグ群, l_walk順群, l_環状フラグ群);
            }

            // 始点を持たない=循環している経路を拾う (環状ゲノム/プラスミド等)
            for (var v = 2; v < p_グラフ.A_出辺.Count; v += 2)
            {
                if (l_訪問済み[v >> 1])
                {
                    continue;
                }
                V_実行_walk(p_ユニティグ配列, p_結合, l_訪問済み, p_重なり長, v, l_コンティグ群, l_walk順群, l_環状フラグ群);
            }

            using var l_書き込み = new FastaWriter(p_コンティグパス);
            var l_ID = 1;
            var l_総延長 = 0L;
            for (var c = 0; c < l_コンティグ群.Count; c++)
            {
                var l_コンティグ = l_コンティグ群[c];
                var l_walk順 = l_walk順群[c];
                var l_逆相補 = Util.V_逆相補(l_コンティグ);
                var l_Is逆相補採用 = string.CompareOrdinal(l_コンティグ, l_逆相補) > 0;

                // 環状に閉じた contig は、その複製単位 (染色体・プラスミド) を
                // 完全に組み上げられたことを意味するため、名前に明示する
                // 閉じていても、複製単位と呼べる長さが無ければ目印は付けない
                // ホモポリマー由来の 1 bp の閉路まで環状のレプリコンとして数えると、
                // 候補選択も完全性の判定もその雑音に従ってしまう
                var l_Is複製単位 = l_環状フラグ群[c] && l_コンティグ群[c].Length >= Consts.環状として数える最小長;
                var l_名前 = l_Is複製単位 ? $"NODE{l_ID}_{Consts.環状の目印}" : $"NODE{l_ID}";
                var l_出力配列 = l_Is逆相補採用 ? l_逆相補 : l_コンティグ;

                if (l_環状フラグ群[c])
                {
                    // 環状配列は開始位置が任意 (walk がどこから始まったかの
                    // 産物でしかない)
                    // 下流の比較・再現性のため、辞書式順序で
                    // 最小になる回転へ正規化する (鎖の向きは上の比較で既に
                    // 決めているため、ここでは回転のみ
                    // 鎖の選択まで変えると
                    // 下の unitig 配置 (逆相補か) の記録と食い違う)
                    l_出力配列 = Util.Get_最小回転(l_出力配列);
                }
                l_書き込み.V_書き込み(l_名前, l_出力配列);

                // 逆相補を採用した場合、出力配列は walk 順と逆向きになる
                // 位置 (先頭/末尾) の解釈は Scaffolder 側で反転させる
                for (var w = 0; w < l_walk順.Count; w++)
                {
                    var l_頂点番号 = l_walk順[w];
                    this._ユニティグ配置[l_頂点番号 >> 1] = new ユニティグ配置(l_ID, l_Is逆相補採用, w, l_walk順.Count, (l_頂点番号 & 1) == 1);
                }

                l_ID++;
                l_総延長 += l_コンティグ.Length;
            }
            Logger.V_出力(メッセージID.コンティグ総延長, l_総延長);

            var l_環状コンティグ = Enumerable.Range(0, l_コンティグ群.Count)
                .Where(x => l_環状フラグ群[x] && l_コンティグ群[x].Length >= Consts.環状として数える最小長)
                .ToList();
            var l_短すぎる閉路 = Enumerable.Range(0, l_コンティグ群.Count)
                .Count(x => l_環状フラグ群[x] && l_コンティグ群[x].Length < Consts.環状として数える最小長);

            if (l_短すぎる閉路 > 0)
            {
                Logger.V_出力(メッセージID.短すぎる閉路, l_短すぎる閉路, Consts.環状として数える最小長);
            }

            if (l_環状コンティグ.Count > 0)
            {
                var l_長さ一覧 = string.Join(", ", l_環状コンティグ.Select(x => $"{l_コンティグ群[x].Length}bp"));
                Logger.V_出力(メッセージID.環状コンティグあり, l_環状コンティグ.Count, l_長さ一覧);
            }
            else
            {
                Logger.V_出力(メッセージID.環状コンティグなし);
            }
        }

        #endregion
    }
}
