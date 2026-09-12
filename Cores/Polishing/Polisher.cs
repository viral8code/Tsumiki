using Tsumiki.Commons;
using Tsumiki.Cores.Mapping;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Polishing;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Polishing
{
    /// <summary>
    /// 最終配列に元リードを貼り直し、各位置の塩基の多数決で置換を直す
    /// </summary>
    internal static class Polisher
    {
        #region 定数

        /// <summary>
        /// リードの置き場所を探す種にする長さ
        /// </summary>
        /// <remarks>
        /// 短いほど反復配列で曖昧になり、長いほどエラーを 1 つ含んだだけで種が潰れる
        /// </remarks>
        private const int シード長 = 21;

        /// <summary>
        /// 投票に使う整列として意味を持つ最小の対応塩基数
        /// </summary>
        private const int 最小整列長 = シード長;

        /// <summary>
        /// 置換を認めるのに必要な、その位置の深度
        /// </summary>
        private const int 訂正に必要な深度 = 5;

        /// <summary>
        /// 置換を認めるのに必要な、対立塩基の占有率
        /// </summary>
        /// <remarks>
        /// 元の塩基が少数派というだけでは足りず、対立塩基が明確に多数派でなければ動かさない
        /// </remarks>
        private const double 訂正に必要な占有率 = 0.7D;

        /// <summary>
        /// 深度が不足しているとみなす、中央値に対する比
        /// </summary>
        private const double 深度不足とみなす比 = 0.2D;

        /// <summary>
        /// 深度のヒストグラムを取る上限
        /// </summary>
        /// <remarks>
        /// これ以上は同じ枠に入れる
        /// </remarks>
        private const int 深度ヒストグラムの上限 = 65_535;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_FASTAパス を磨いて p_出力パス へ書き出す
        /// </summary>
        /// <param name="p_FASTAパス"></param>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <param name="p_出力パス"></param>
        /// <param name="p_Is訂正">false なら配列を変更せず深度を再測定する</param>
        /// <remarks>
        /// 磨く対象が無い (配列が空、種が 1 つも取れない) 場合は null を返す
        /// </remarks>
        /// <returns></returns>
        public static ポリッシュ統計? Get_磨いた結果(string p_FASTAパス, string p_リード1のパス, string? p_リード2のパス, string p_出力パス, bool p_Is訂正 = true)
        {
            var l_エントリ群 = FastaReader.Get_全エントリ(p_FASTAパス);
            if (l_エントリ群.Count == 0)
            {
                return null;
            }

            var l_配列群 = l_エントリ群.Select(x => x.A_配列.ToCharArray()).ToList();
            var l_総延長 = l_配列群.Sum(x => (long)x.Length);

            Logger.V_出力(メッセージID.ポリッシュの索引構築, l_エントリ群.Count, l_総延長);
            var l_マッパー = new ReadMapper([.. l_配列群.Select(x => new string(x))]);
            if (l_配列群.All(x => x.Length < シード長))
            {
                Logger.V_出力(メッセージID.ポリッシュの種が無い, シード長);
                return null;
            }

            // 位置ごとの塩基の得票
            // 塩基 ID は 1..4 なので (位置 * 4 + 塩基 ID - 1)
            var l_得票 = l_配列群.Select(x => new int[x.Length * 4]).ToArray();

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_マップ数 = new long[l_スレッド数];
            var l_棄却数 = new long[l_スレッド数];

            Logger.V_出力(メッセージID.ポリッシュのマッピング開始);
            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, FastqReader.Get_生リード列(p_リード1のパス, p_リード2のパス), (l_リード, l_ワーカー番号) =>
            {
                if (Try集計_塩基票(l_リード, l_マッパー, l_得票))
                {
                    l_マップ数[l_ワーカー番号]++;
                }
                else
                {
                    l_棄却数[l_ワーカー番号]++;
                }
            });

            var l_中央値 = Get_深度中央値(l_配列群, l_得票);
            Logger.V_出力(メッセージID.ポリッシュの深度中央値, l_中央値);

            var l_訂正数 = V_訂正_多数決(l_配列群, l_得票, l_中央値, out var l_深度不足数, out var l_評価位置数, p_Is訂正);

            using (var l_書き込み = new FastaWriter(p_出力パス))
            {
                for (var i = 0; i < l_エントリ群.Count; i++)
                {
                    l_書き込み.V_書き込み(l_エントリ群[i].A_ID, new string(l_配列群[i]));
                }
            }

            return new ポリッシュ統計(A_配列数: l_エントリ群.Count, A_総延長: l_総延長, A_マップされたリード数: l_マップ数.Sum(), A_棄却されたリード数: l_棄却数.Sum(), A_訂正した塩基数: l_訂正数, A_深度不足の位置数: l_深度不足数, A_評価できた位置数: l_評価位置数, A_深度の中央値: l_中央値);
        }

        /// <summary>
        /// ポリッシュの結果をログへ出力する
        /// </summary>
        /// <param name="p_統計">ポリッシュの結果</param>
        public static void V_出力_統計(ポリッシュ統計? p_統計)
        {
            if (p_統計 is not { } l_統計)
            {
                Logger.V_出力(メッセージID.ポリッシュを行えず);
                return;
            }
            Logger.V_出力(メッセージID.ポリッシュのマッピング結果, l_統計.A_マップされたリード数, l_統計.A_棄却されたリード数);
            Logger.V_出力(メッセージID.ポリッシュの訂正結果, l_統計.A_訂正した塩基数, l_統計.A_総延長);
            Logger.V_出力(メッセージID.ポリッシュの深度不足, l_統計.A_深度不足の位置数, l_統計.A_評価できた位置数, l_統計.A_深度不足率 * 100D);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 本のリードを置ける場所へ置き、各位置の得票を加算する
        /// </summary>
        /// <param name="p_リード"></param>
        /// <param name="p_マッパー"></param>
        /// <param name="p_得票"></param>
        /// <remarks>
        /// 置けたら true<br/>
        /// 得票は複数のワーカーが同じ配列を触るため Interlocked で足す<br/>
        /// 加算は順序に依らないので、並列でも結果は毎回同じになる
        /// </remarks>
        /// <returns></returns>
        private static bool Try集計_塩基票(string p_リード, ReadMapper p_マッパー, int[][] p_得票)
        {
            var l_配置 = p_マッパー.Get_配置(p_リード);
            if (l_配置.A_配列番号 < 0 || l_配置.A_信頼度 == 0 || l_配置.A_整列位置群.Count < 最小整列長)
            {
                return false;
            }

            foreach (var l_整列位置 in l_配置.A_整列位置群)
            {
                var l_塩基ID = Util.Get_塩基ID(p_リード[l_整列位置.A_リード位置]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    continue;
                }

                // A_リード位置 は元のリードの向きでの添字なので、逆鎖に載ったリードは
                // 参照と同じ向きにするため相補を取ってから投票する
                if (l_配置.A_Is逆鎖)
                {
                    l_塩基ID = Util.Get_相補塩基ID(l_塩基ID);
                }
                _ = Interlocked.Increment(ref p_得票[l_配置.A_配列番号][(l_整列位置.A_参照位置 * 4) + l_塩基ID - 1]);
            }
            return l_配置.A_整列位置群.Count >= 最小整列長;
        }

        /// <summary>
        /// リードが載っている位置での深度の中央値
        /// </summary>
        /// <param name="p_配列群"></param>
        /// <param name="p_得票"></param>
        /// <remarks>
        /// ヒストグラムから求めるのは、位置ごとの値をすべて並べるとゲノムサイズぶんの配列をもう 1 本持つことになるため<br/>
        /// 深度 0 の位置を混ぜないのが要点<br/>
        /// 混ぜると、覆われていない範囲が広いアセンブリほど中央値が 0 へ引き寄せられ、本来そこを咎めるはずの深度不足の判定が何も引っ掛けなくなる
        /// </remarks>
        /// <returns></returns>
        private static double Get_深度中央値(List<char[]> p_配列群, int[][] p_得票)
        {
            var l_ヒストグラム = new long[深度ヒストグラムの上限 + 1];
            var l_総数 = 0L;
            for (var i = 0; i < p_配列群.Count; i++)
            {
                for (var l_位置 = 0; l_位置 < p_配列群[i].Length; l_位置++)
                {
                    if (p_配列群[i][l_位置] == 'N')
                    {
                        continue;
                    }
                    var l_深度 = Get_深度(p_得票[i], l_位置);
                    if (l_深度 == 0)
                    {
                        continue;
                    }
                    l_ヒストグラム[Math.Min(深度ヒストグラムの上限, l_深度)]++;
                    l_総数++;
                }
            }
            if (l_総数 == 0L)
            {
                return 0D;
            }

            var l_累積 = 0L;
            for (var l_深度 = 0; l_深度 <= 深度ヒストグラムの上限; l_深度++)
            {
                l_累積 += l_ヒストグラム[l_深度];
                if (l_累積 * 2L >= l_総数)
                {
                    return l_深度;
                }
            }
            return 0D;
        }

        /// <summary>
        /// その位置に載ったリードの本数を返す
        /// </summary>
        /// <param name="p_得票">位置と塩基ごとの得票</param>
        /// <param name="p_位置">調べる位置</param>
        /// <returns>載ったリードの本数</returns>
        private static int Get_深度(int[] p_得票, int p_位置)
        {
            var l_起点 = p_位置 * 4;
            return p_得票[l_起点] + p_得票[l_起点 + 1] + p_得票[l_起点 + 2] + p_得票[l_起点 + 3];
        }

        /// <summary>
        /// 得票の多数決で置換を適用する
        /// </summary>
        /// <param name="p_配列群"></param>
        /// <param name="p_得票"></param>
        /// <param name="p_深度の中央値"></param>
        /// <param name="p_深度不足数"></param>
        /// <param name="p_評価位置数"></param>
        /// <param name="p_Is訂正">置換を適用するか</param>
        /// <remarks>
        /// 併せて深度不足の位置を数える<br/>
        /// N の位置は触らない<br/>
        /// ギャップの長さは推定値であり、そこを塩基で埋めるのは多数決の仕事ではない
        /// </remarks>
        /// <returns></returns>
        private static long V_訂正_多数決(List<char[]> p_配列群, int[][] p_得票, double p_深度の中央値, out long p_深度不足数, out long p_評価位置数, bool p_Is訂正)
        {
            var l_深度不足の閾値 = p_深度の中央値 * 深度不足とみなす比;

            var l_訂正数 = 0L;
            p_深度不足数 = 0L;
            p_評価位置数 = 0L;

            for (var i = 0; i < p_配列群.Count; i++)
            {
                var l_配列 = p_配列群[i];
                var l_票 = p_得票[i];
                for (var l_位置 = 0; l_位置 < l_配列.Length; l_位置++)
                {
                    if (l_配列[l_位置] == 'N')
                    {
                        continue;
                    }
                    p_評価位置数++;

                    var l_深度 = Get_深度(l_票, l_位置);
                    if (l_深度 == 0 || l_深度 < l_深度不足の閾値)
                    {
                        p_深度不足数++;
                    }
                    if (!p_Is訂正 || l_深度 < 訂正に必要な深度)
                    {
                        continue;
                    }

                    var l_起点 = l_位置 * 4;
                    var l_最多の塩基ID = Consts.塩基ID.A;
                    var l_最多得票 = l_票[l_起点];
                    for (var l_塩基ID = Consts.塩基ID.C; l_塩基ID <= Consts.塩基ID.T; l_塩基ID++)
                    {
                        if (l_票[l_起点 + l_塩基ID - 1] > l_最多得票)
                        {
                            l_最多得票 = l_票[l_起点 + l_塩基ID - 1];
                            l_最多の塩基ID = l_塩基ID;
                        }
                    }

                    var l_最多の塩基 = Util.Get_塩基文字(l_最多の塩基ID);
                    if (l_最多の塩基 == l_配列[l_位置] || l_最多得票 < l_深度 * 訂正に必要な占有率)
                    {
                        continue;
                    }
                    l_配列[l_位置] = l_最多の塩基;
                    l_訂正数++;
                }
            }
            return l_訂正数;
        }

        #endregion
    }
}
