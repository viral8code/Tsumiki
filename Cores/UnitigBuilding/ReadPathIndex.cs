using Tsumiki.Core;
using Tsumiki.Models.ContigBuilding;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// リードが通り抜けた向き付き unitig の並びを数え、どの入口からどの出口へ抜けたかを引けるようにする
    /// </summary>
    /// <remarks>
    /// 隣り合う 2 頂点の組では、反復 R を挟んだ入口と出口の対応が消える<br/>
    /// 並びは向き付き頂点番号 (2*ID+鎖) で持ち、逆鎖側から読んだ出現も同じ並びとして数える
    /// </remarks>
    internal sealed class ReadPathIndex
    {
        #region 定数

        /// <summary>
        /// 索引に入れる並びの最短の頂点数
        /// </summary>
        /// <remarks>
        /// 2 頂点の並びは隣接の支持と同じ情報しか持たない
        /// </remarks>
        public const int 最短の頂点数 = 3;

        #endregion

        #region 内部変数

        /// <summary>
        /// 並び
        /// </summary>
        private readonly List<int[]> _経路 = [];

        /// <summary>
        /// 並びごとの独立なリード・ペアの数
        /// </summary>
        private readonly List<ulong> _件数 = [];

        /// <summary>
        /// unitig ID から、それを含む並びの番号を引く
        /// </summary>
        private readonly Dictionary<int, List<int>> _unitig別の経路 = [];

        #endregion

        #region プロパティ

        /// <summary>
        /// 並びの種類数
        /// </summary>
        public int A_経路数 => this._経路.Count;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 符号付き unitig ID で集計した並びから索引を作る
        /// </summary>
        /// <param name="p_集計"></param>
        /// <returns></returns>
        public static ReadPathIndex Get_索引(IEnumerable<KeyValuePair<経路キー, ulong>> p_集計)
        {
            var l_索引 = new ReadPathIndex();
            foreach (var (l_キー, l_件数) in p_集計)
            {
                var l_頂点列 = new int[l_キー.A_頂点列.Length];
                for (var i = 0; i < l_頂点列.Length; i++)
                {
                    l_頂点列[i] = ContigMaker.Get_頂点番号(l_キー.A_頂点列[i]);
                }
                l_索引.V_追加(l_頂点列, l_件数);
            }
            return l_索引;
        }

        /// <summary>
        /// 符号付き unitig ID の並びを、逆鎖から読んでも同じになる向きに揃える
        /// </summary>
        /// <param name="p_符号付きID列"></param>
        /// <returns></returns>
        public static int[] Get_正準経路(ReadOnlySpan<int> p_符号付きID列)
        {
            var l_長さ = p_符号付きID列.Length;
            var l_Is逆鎖採用 = false;
            for (var i = 0; i < l_長さ; i++)
            {
                var l_順 = p_符号付きID列[i];
                var l_逆 = -p_符号付きID列[l_長さ - 1 - i];
                if (l_順 != l_逆)
                {
                    l_Is逆鎖採用 = l_逆 < l_順;
                    break;
                }
            }

            var l_結果 = new int[l_長さ];
            for (var i = 0; i < l_長さ; i++)
            {
                l_結果[i] = l_Is逆鎖採用 ? -p_符号付きID列[l_長さ - 1 - i] : p_符号付きID列[i];
            }
            return l_結果;
        }

        /// <summary>
        /// 向き付き頂点番号の並びを 1 件足す
        /// </summary>
        /// <param name="p_頂点列"></param>
        /// <param name="p_件数"></param>
        public void V_追加(int[] p_頂点列, ulong p_件数)
        {
            var l_番号 = this._経路.Count;
            this._経路.Add(p_頂点列);
            this._件数.Add(p_件数);
            foreach (var v in p_頂点列)
            {
                this.V_登録_unitig別(v >> 1, l_番号);
            }
        }

        /// <summary>
        /// 部分列をどちらかの向きで連続して含む並びの件数の合計
        /// </summary>
        /// <param name="p_部分列">向き付き頂点番号の並び</param>
        /// <remarks>
        /// 1 つの並びに部分列が何度現れても 1 回と数える (1 本のリードは 1 回の観測)
        /// </remarks>
        /// <returns></returns>
        public ulong Get_出現数(ReadOnlySpan<int> p_部分列)
        {
            var l_長さ = p_部分列.Length;
            if (l_長さ == 0 || !this._unitig別の経路.TryGetValue(p_部分列[l_長さ / 2] >> 1, out var l_番号群))
            {
                return 0UL;
            }

            Span<int> l_逆部分列 = l_長さ <= 64 ? stackalloc int[l_長さ] : new int[l_長さ];
            for (var i = 0; i < l_長さ; i++)
            {
                l_逆部分列[i] = p_部分列[l_長さ - 1 - i] ^ 1;
            }

            var l_合計 = 0UL;
            foreach (var l_番号 in l_番号群)
            {
                var l_経路 = this._経路[l_番号].AsSpan();
                if (l_経路.IndexOf(p_部分列) >= 0 || l_経路.IndexOf(l_逆部分列) >= 0)
                {
                    l_合計 += this._件数[l_番号];
                }
            }
            return l_合計;
        }

        /// <summary>
        /// 上流の並びに続けて読まれた行き先の件数が、候補の 1 つに偏っていればそれを返す
        /// </summary>
        /// <param name="p_上流">起点から分岐元までの頂点</param>
        /// <param name="p_候補">分岐元の行き先</param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <returns>決まらなければ null</returns>
        public int? Get_優勢な行き先(IReadOnlyList<int> p_上流, IReadOnlyList<int> p_候補, decimal p_優勢閾値, ulong p_最小証拠数)
        {
            var l_部分列 = new int[p_上流.Count + 1];
            for (var i = 0; i < p_上流.Count; i++)
            {
                l_部分列[i] = p_上流[i];
            }

            var l_最良 = -1;
            var l_最良件数 = 0UL;
            var l_合計 = 0UL;
            var l_Is同点 = false;
            foreach (var w in p_候補)
            {
                l_部分列[^1] = w;
                var l_件数 = this.Get_出現数(l_部分列);
                l_合計 += l_件数;
                if (l_件数 > l_最良件数)
                {
                    (l_最良, l_最良件数, l_Is同点) = (w, l_件数, false);
                }
                else if (l_件数 == l_最良件数 && l_件数 > 0UL)
                {
                    l_Is同点 = true;
                }
            }

            return l_最良 < 0 || l_Is同点 || l_最良件数 < p_最小証拠数 || (decimal)l_最良件数 / l_合計 < p_優勢閾値 ? null : l_最良;
        }

        /// <summary>
        /// 反復の鎖を複製して入口 1 本と出口 1 本を移したあと、移した側を通った並びを複製の頂点へ書き換える
        /// </summary>
        /// <param name="p_鎖">元の鎖の頂点 (入口側から順に)</param>
        /// <param name="p_複製鎖">p_鎖 と同じ順の複製の頂点</param>
        /// <param name="p_移す入口"></param>
        /// <param name="p_移す出口"></param>
        /// <param name="p_入口群">付け替え前に元の鎖へ入っていた頂点</param>
        /// <param name="p_出口群">付け替え前に元の鎖から出ていた頂点</param>
        /// <remarks>
        /// 入口側と出口側で移した側・残した側の判定が食い違う並びは、どちらのコピーを通ったか決められないので件数を 0 にする<br/>
        /// 鎖の外に手がかりの無い並びは元の頂点のまま残す (入口と出口の対応の証拠にはならない)
        /// </remarks>
        public void V_付け替え_複製(IReadOnlyList<int> p_鎖, IReadOnlyList<int> p_複製鎖, int p_移す入口, int p_移す出口, IReadOnlyCollection<int> p_入口群, IReadOnlyCollection<int> p_出口群)
        {
            SortedSet<int> l_候補 = [];
            foreach (var l_頂点 in p_鎖)
            {
                if (this._unitig別の経路.TryGetValue(l_頂点 >> 1, out var l_番号群))
                {
                    l_候補.UnionWith(l_番号群);
                }
            }

            foreach (var l_番号 in l_候補)
            {
                var l_結果 = V_書き換え_経路(this._経路[l_番号], p_鎖, p_複製鎖, p_移す入口, p_移す出口, p_入口群, p_出口群);
                if (l_結果 < 0)
                {
                    this._件数[l_番号] = 0UL;
                }
                else if (l_結果 > 0)
                {
                    foreach (var l_頂点 in p_複製鎖)
                    {
                        this.V_登録_unitig別(l_頂点 >> 1, l_番号);
                    }
                }
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// unitig ID から並びの番号を引けるようにする
        /// </summary>
        /// <param name="p_unitigID"></param>
        /// <param name="p_番号"></param>
        private void V_登録_unitig別(int p_unitigID, int p_番号)
        {
            if (!this._unitig別の経路.TryGetValue(p_unitigID, out var l_番号群))
            {
                l_番号群 = [];
                this._unitig別の経路[p_unitigID] = l_番号群;
            }
            if (l_番号群.Count == 0 || l_番号群[^1] != p_番号)
            {
                l_番号群.Add(p_番号);
            }
        }

        /// <summary>
        /// 1 本の並びの中で鎖を通った区間を、移した側なら複製の頂点へ書き換える
        /// </summary>
        /// <param name="p_経路"></param>
        /// <param name="p_鎖"></param>
        /// <param name="p_複製鎖"></param>
        /// <param name="p_移す入口"></param>
        /// <param name="p_移す出口"></param>
        /// <param name="p_入口群"></param>
        /// <param name="p_出口群"></param>
        /// <returns>書き換えたら 1、変更なしなら 0、判定が食い違ったら -1</returns>
        private static int V_書き換え_経路(int[] p_経路, IReadOnlyList<int> p_鎖, IReadOnlyList<int> p_複製鎖, int p_移す入口, int p_移す出口, IReadOnlyCollection<int> p_入口群, IReadOnlyCollection<int> p_出口群)
        {
            var l_鎖長 = p_鎖.Count;
            var l_Is変更 = false;
            for (var j = 0; j < p_経路.Length; j++)
            {
                var l_順位置 = Get_鎖内位置(p_鎖, p_経路[j]);
                var l_逆位置 = Get_鎖内位置(p_鎖, p_経路[j] ^ 1);
                if (l_順位置 < 0 && l_逆位置 < 0)
                {
                    continue;
                }

                // 逆鎖から読んだ区間では、並びの左が鎖の出口側になる
                var l_Is順 = l_順位置 >= 0;
                var l_歩幅 = l_Is順 ? 1 : -1;
                var l_開始 = j;
                var l_開始位置 = l_Is順 ? l_順位置 : l_逆位置;
                var l_現在位置 = l_開始位置;
                while (j + 1 < p_経路.Length)
                {
                    var l_次位置 = l_現在位置 + l_歩幅;
                    if (l_次位置 < 0 || l_次位置 >= l_鎖長 || p_経路[j + 1] != (l_Is順 ? p_鎖[l_次位置] : p_鎖[l_次位置] ^ 1))
                    {
                        break;
                    }
                    j++;
                    l_現在位置 = l_次位置;
                }

                int? l_入口 = null;
                int? l_出口 = null;
                if (l_Is順)
                {
                    l_入口 = l_開始位置 == 0 && l_開始 > 0 ? p_経路[l_開始 - 1] : null;
                    l_出口 = l_現在位置 == l_鎖長 - 1 && j + 1 < p_経路.Length ? p_経路[j + 1] : null;
                }
                else
                {
                    l_出口 = l_開始位置 == l_鎖長 - 1 && l_開始 > 0 ? p_経路[l_開始 - 1] ^ 1 : null;
                    l_入口 = l_現在位置 == 0 && j + 1 < p_経路.Length ? p_経路[j + 1] ^ 1 : null;
                }

                var l_入口判定 = Get_側判定(l_入口, p_移す入口, p_入口群);
                var l_出口判定 = Get_側判定(l_出口, p_移す出口, p_出口群);
                if (l_入口判定 * l_出口判定 < 0)
                {
                    return -1;
                }
                if (l_入口判定 + l_出口判定 <= 0)
                {
                    continue;
                }

                for (var x = l_開始; x <= j; x++)
                {
                    var l_位置 = l_開始位置 + (l_歩幅 * (x - l_開始));
                    p_経路[x] = l_Is順 ? p_複製鎖[l_位置] : p_複製鎖[l_位置] ^ 1;
                }
                l_Is変更 = true;
            }
            return l_Is変更 ? 1 : 0;
        }

        /// <summary>
        /// 鎖の中での位置
        /// </summary>
        /// <param name="p_鎖"></param>
        /// <param name="p_頂点"></param>
        /// <returns>含まれなければ -1</returns>
        private static int Get_鎖内位置(IReadOnlyList<int> p_鎖, int p_頂点)
        {
            for (var i = 0; i < p_鎖.Count; i++)
            {
                if (p_鎖[i] == p_頂点)
                {
                    return i;
                }
            }
            return -1;
        }

        /// <summary>
        /// 鎖の外側の頂点が、移した側か残した側か
        /// </summary>
        /// <param name="p_頂点"></param>
        /// <param name="p_移す頂点"></param>
        /// <param name="p_全頂点"></param>
        /// <returns>移した側なら 1、残した側なら -1、手がかりにならなければ 0</returns>
        private static int Get_側判定(int? p_頂点, int p_移す頂点, IReadOnlyCollection<int> p_全頂点)
        {
            return p_頂点 is not { } l_頂点 ? 0 : l_頂点 == p_移す頂点 ? 1 : p_全頂点.Contains(l_頂点) ? -1 : 0;
        }

        #endregion
    }
}
