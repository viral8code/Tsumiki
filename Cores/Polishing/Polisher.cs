using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Polishing;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Polishing
{
    /// <summary>
    /// 最終配列に元リードを貼り直し、各位置の塩基の多数決で置換を直す
    /// </summary>
    /// <remarks>
    /// グラフから組み立てた配列の塩基は k-mer 集合が根拠であり、カットオフを
    /// 通り抜けたエラー k-mer がそのまま残ることがある<br/>
    /// リードそのものの
    /// 多数決は k-mer とは独立した証拠で、そこを最後に均す<br/>
    /// 直すのは置換だけとする<br/>
    /// ここの照合は ungapped なので挿入・欠失は
    /// そもそも検出できず、検出できないものを直そうとすれば配列を壊す<br/>
    /// 併せて位置ごとの深度が得られる<br/>
    /// 連結の裏付けが無い接合点はその前後で
    /// 深度が不連続になるため、完全長の判定にも使う
    /// </remarks>
    /// <summary>
    /// 種索引の 1 件
    /// </summary>
    /// <remarks>
    /// A_配列番号 が負の値なら複数箇所に当たる曖昧な種
    /// </remarks>
    internal readonly record struct 種の位置(int A_配列番号, int A_位置, bool A_逆鎖);

    internal static class Polisher
    {
        /// <summary>
        /// リードの置き場所を探す種にする長さ
        /// </summary>
        /// <remarks>
        /// 短いほど反復配列で曖昧になり、
        /// 長いほどエラーを 1 つ含んだだけで種が潰れる
        /// </remarks>
        private const int シード長 = 31;

        /// <summary>
        /// 参照側で種を登録する間隔
        /// </summary>
        /// <remarks>
        /// リード長はこれよりはるかに長いため、
        /// 間引いてもリードのどこかは必ず登録済みの位置に重なる<br/>
        /// 全位置を持つと索引がゲノムサイズそのものの規模になる
        /// </remarks>
        private const int シード間隔 = 8;

        /// <summary>
        /// 1 本のリードにつき試す種ヒットの数
        /// </summary>
        /// <remarks>
        /// 反復配列では種が当たっても
        /// 照合に落ちることが続くため、諦める上限を決める
        /// </remarks>
        private const int 試すヒット数 = 4;

        /// <summary>
        /// 照合を認める不一致の割合
        /// </summary>
        /// <remarks>
        /// これを超えたら別の場所とみなす
        /// </remarks>
        private const double 許容不一致率 = 0.1D;

        /// <summary>
        /// ungapped 照合として意味を持つ最小の重なり長
        /// </summary>
        private const int 最小重なり長 = シード長;

        /// <summary>
        /// 置換を認めるのに必要な、その位置の深度
        /// </summary>
        private const int 訂正に必要な深度 = 5;

        /// <summary>
        /// 置換を認めるのに必要な、対立塩基の占有率
        /// </summary>
        /// <remarks>
        /// 元の塩基が少数派というだけでは足りず、対立塩基が明確に
        /// 多数派でなければ動かさない
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
        private const int 深度ヒストグラムの上限 = 65535;

        /// <summary>
        /// 複数の位置に当たった種であることを示す番兵
        /// </summary>
        private const int 曖昧な種の番兵 = -1;

        /// <summary>
        /// p_FASTAパス を磨いて p_出力パス へ書き出す
        /// </summary>
        /// <remarks>
        /// 磨く対象が無い (配列が空、種が 1 つも取れない) 場合は null を返す
        /// </remarks>
        public static ポリッシュ統計? Get_磨いた結果(
            string p_FASTAパス, string p_リード1のパス, string? p_リード2のパス, string p_出力パス)
        {
            var l_エントリ群 = FastaReader.Get_全エントリ(p_FASTAパス);
            if (l_エントリ群.Count == 0)
            {
                return null;
            }

            var l_配列群 = l_エントリ群.Select(x => x.A_配列.ToCharArray()).ToList();
            var l_総延長 = l_配列群.Sum(x => (long)x.Length);

            Logger.V_出力(メッセージID.ポリッシュの索引構築, l_エントリ群.Count, l_総延長);
            var l_種索引 = Get_種索引(l_配列群);
            if (l_種索引.Count == 0)
            {
                Logger.V_出力(メッセージID.ポリッシュの種が無い, シード長);
                return null;
            }

            // 位置ごとの塩基の得票
            // 塩基 ID は 1..4 なので (位置 * 4 + 塩基ID - 1)
            var l_得票 = l_配列群.Select(x => new int[x.Length * 4]).ToArray();

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_マップ数 = new long[l_スレッド数];
            var l_棄却数 = new long[l_スレッド数];

            Logger.V_出力(メッセージID.ポリッシュのマッピング開始);
            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 256,
                FastqReader.Get_生リード列(p_リード1のパス, p_リード2のパス),
                (l_リード, l_ワーカー番号) =>
                {
                    if (V_貼り付け_1リード(l_リード, l_種索引, l_配列群, l_得票))
                    {
                        l_マップ数[l_ワーカー番号]++;
                    }
                    else
                    {
                        l_棄却数[l_ワーカー番号]++;
                    }
                });

            var l_中央値 = Get_深度の中央値(l_配列群, l_得票);
            Logger.V_出力(メッセージID.ポリッシュの深度中央値, l_中央値);

            var l_訂正数 = V_訂正_多数決(
                l_配列群, l_得票, l_中央値, out var l_深度不足数, out var l_評価位置数);

            using (var l_書き込み = new FastaWriter(p_出力パス))
            {
                for (var i = 0; i < l_エントリ群.Count; i++)
                {
                    l_書き込み.V_書き込み(l_エントリ群[i].A_ID, new string(l_配列群[i]));
                }
            }

            return new ポリッシュ統計(
                A_配列数: l_エントリ群.Count,
                A_総延長: l_総延長,
                A_マップされたリード数: l_マップ数.Sum(),
                A_棄却されたリード数: l_棄却数.Sum(),
                A_訂正した塩基数: l_訂正数,
                A_深度不足の位置数: l_深度不足数,
                A_評価できた位置数: l_評価位置数,
                A_深度の中央値: l_中央値);
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
            Logger.V_出力(
                メッセージID.ポリッシュのマッピング結果,
                l_統計.A_マップされたリード数, l_統計.A_棄却されたリード数);
            Logger.V_出力(メッセージID.ポリッシュの訂正結果, l_統計.A_訂正した塩基数, l_統計.A_総延長);
            Logger.V_出力(
                メッセージID.ポリッシュの深度不足,
                l_統計.A_深度不足の位置数, l_統計.A_評価できた位置数, l_統計.A_深度不足率 * 100);
        }

        /// <summary>
        /// 参照配列から種索引を作る
        /// </summary>
        /// <remarks>
        /// 順鎖と逆相補を別のキーで登録し、
        /// リードがどちらの向きで載ったかを引けるようにする<br/>
        /// 複数の位置に当たる種は反復配列由来なので、曖昧として使わない
        /// </remarks>
        private static Dictionary<UInt128, 種の位置> Get_種索引(List<char[]> p_配列群)
        {
            Dictionary<UInt128, 種の位置> l_索引 = [];
            for (var i = 0; i < p_配列群.Count; i++)
            {
                var l_配列 = new string(p_配列群[i]);
                for (var l_位置 = 0; l_位置 + シード長 <= l_配列.Length; l_位置 += シード間隔)
                {
                    if (!KmerPacking.Get_パック(l_配列, l_位置, シード長, out var l_順鎖))
                    {
                        continue;
                    }
                    V_登録_種(l_索引, l_順鎖, new 種の位置(i, l_位置, false));
                    V_登録_種(
                        l_索引, KmerPacking.Get_逆相補(l_順鎖, シード長), new 種の位置(i, l_位置, true));
                }
            }
            return l_索引;
        }

        /// <summary>
        /// リードの貼り付け位置を探す種を索引へ登録する
        /// </summary>
        /// <param name="p_索引">登録先の索引</param>
        /// <param name="p_キー">種の k-mer</param>
        /// <param name="p_値">その種が指す位置</param>
        private static void V_登録_種(
            Dictionary<UInt128, 種の位置> p_索引, UInt128 p_キー, 種の位置 p_値)
        {
            if (p_索引.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.A_配列番号 != 曖昧な種の番兵)
                {
                    p_索引[p_キー] = l_既存 with { A_配列番号 = 曖昧な種の番兵 };
                }
                return;
            }
            p_索引[p_キー] = p_値;
        }

        /// <summary>
        /// 1 本のリードを置ける場所へ置き、各位置の得票を加算する
        /// </summary>
        /// <remarks>
        /// 置けたら true<br/>
        /// 得票は複数のワーカーが同じ配列を触るため Interlocked で足す<br/>
        /// 加算は順序に依らないので、並列でも結果は毎回同じになる
        /// </remarks>
        private static bool V_貼り付け_1リード(
            string p_リード,
            Dictionary<UInt128, 種の位置> p_種索引,
            List<char[]> p_配列群,
            int[][] p_得票)
        {
            if (p_リード.Length < シード長)
            {
                return false;
            }

            string? l_逆相補 = null;
            var l_マスク = (UInt128.One << (2 * シード長)) - 1;
            UInt128 l_順鎖 = 0;
            var l_直近の曖昧位置 = -1;
            var l_試した回数 = 0;

            for (var i = 0; i < p_リード.Length; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_リード[i]);
                var l_有効か = l_塩基ID is >= Consts.塩基ID.A and <= Consts.塩基ID.T;
                l_順鎖 = ((l_順鎖 << 2) | (UInt128)(l_有効か ? l_塩基ID - 1 : 0)) & l_マスク;
                if (!l_有効か)
                {
                    l_直近の曖昧位置 = i;
                }

                var l_開始 = i - シード長 + 1;
                if (l_開始 < 0 || l_直近の曖昧位置 >= l_開始)
                {
                    continue;
                }
                if (!p_種索引.TryGetValue(l_順鎖, out var l_位置) || l_位置.A_配列番号 == 曖昧な種の番兵)
                {
                    continue;
                }

                // 逆鎖に当たった場合、リードを逆相補にすると参照と同じ向きになる
                // そのとき種はリードの後ろから数えた位置に移る
                string l_照合するリード;
                int l_参照開始;
                if (l_位置.A_逆鎖)
                {
                    l_逆相補 ??= Util.V_逆相補_曖昧塩基あり(p_リード);
                    l_照合するリード = l_逆相補;
                    l_参照開始 = l_位置.A_位置 - (p_リード.Length - l_開始 - シード長);
                }
                else
                {
                    l_照合するリード = p_リード;
                    l_参照開始 = l_位置.A_位置 - l_開始;
                }

                if (V_照合(l_照合するリード, p_配列群[l_位置.A_配列番号], l_参照開始, p_得票[l_位置.A_配列番号]))
                {
                    return true;
                }
                if (++l_試した回数 >= 試すヒット数)
                {
                    return false;
                }
            }
            return false;
        }

        /// <summary>
        /// リードを参照の指定位置へ ungapped に重ね、不一致が許容内なら得票を加算する
        /// </summary>
        /// <remarks>
        /// 参照からはみ出す部分は切り詰める
        /// </remarks>
        private static bool V_照合(
            string p_リード, char[] p_参照, int p_参照開始, int[] p_得票)
        {
            var l_左 = Math.Max(0, p_参照開始);
            var l_右 = Math.Min(p_参照.Length, p_参照開始 + p_リード.Length);
            var l_重なり = l_右 - l_左;
            if (l_重なり < 最小重なり長)
            {
                return false;
            }

            var l_許容不一致数 = (int)(l_重なり * 許容不一致率);
            var l_不一致数 = 0;
            for (var i = l_左; i < l_右; i++)
            {
                var l_参照塩基 = p_参照[i];
                var l_リード塩基 = p_リード[i - p_参照開始];
                if (l_参照塩基 == 'N' || l_リード塩基 == 'N' || l_参照塩基 == l_リード塩基)
                {
                    continue;
                }
                if (++l_不一致数 > l_許容不一致数)
                {
                    return false;
                }
            }

            for (var i = l_左; i < l_右; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_リード[i - p_参照開始]);
                if (l_塩基ID is >= Consts.塩基ID.A and <= Consts.塩基ID.T)
                {
                    _ = Interlocked.Increment(ref p_得票[(i * 4) + l_塩基ID - 1]);
                }
            }
            return true;
        }

        /// <summary>
        /// リードが載っている位置での深度の中央値
        /// </summary>
        /// <remarks>
        /// ヒストグラムから求めるのは、
        /// 位置ごとの値をすべて並べるとゲノムサイズぶんの配列をもう 1 本
        /// 持つことになるため<br/>
        /// 深度 0 の位置を混ぜないのが要点<br/>
        /// 混ぜると、覆われていない範囲が
        /// 広いアセンブリほど中央値が 0 へ引き寄せられ、本来そこを咎めるはずの
        /// 深度不足の判定が何も引っ掛けなくなる
        /// </remarks>
        private static double Get_深度の中央値(List<char[]> p_配列群, int[][] p_得票)
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
            if (l_総数 == 0)
            {
                return 0;
            }

            var l_累積 = 0L;
            for (var l_深度 = 0; l_深度 <= 深度ヒストグラムの上限; l_深度++)
            {
                l_累積 += l_ヒストグラム[l_深度];
                if (l_累積 * 2 >= l_総数)
                {
                    return l_深度;
                }
            }
            return 0;
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
        /// <remarks>
        /// 併せて深度不足の位置を数える<br/>
        /// N の位置は触らない<br/>
        /// ギャップの長さは推定値であり、
        /// そこを塩基で埋めるのは多数決の仕事ではない
        /// </remarks>
        private static long V_訂正_多数決(
            List<char[]> p_配列群, int[][] p_得票, double p_深度の中央値,
            out long p_深度不足数, out long p_評価位置数)
        {
            var l_深度不足の閾値 = p_深度の中央値 * 深度不足とみなす比;

            var l_訂正数 = 0L;
            p_深度不足数 = 0;
            p_評価位置数 = 0;

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
                    if (l_深度 < l_深度不足の閾値)
                    {
                        p_深度不足数++;
                    }
                    if (l_深度 < 訂正に必要な深度)
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
    }
}
