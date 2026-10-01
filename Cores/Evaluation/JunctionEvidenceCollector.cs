using Tsumiki.Commons;
using Tsumiki.Cores.Mapping;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Mapping;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 最終配列へ元リードをペアのまま当て直し、継ぎ目の評価に使う証拠を集める
    /// </summary>
    internal static class JunctionEvidenceCollector
    {
        #region 定数

        /// <summary>
        /// 一意に当たったとみなす信頼度
        /// </summary>
        internal const int C_一意とみなす信頼度 = 20;

        /// <summary>
        /// 切れ端として数える、整列から外れたリード端の長さ
        /// </summary>
        private const int C_切れ端の最小長 = 20;

        /// <summary>
        /// ペアの断片として認める長さの上限
        /// </summary>
        private const int C_断片長の上限 = 5_000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_配列群 へ p_ライブラリ群 のリードを当て、証拠を集める
        /// </summary>
        /// <param name="p_配列群">当てる先の配列</param>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組 (リード 2 が空なら片側だけのライブラリ)</param>
        /// <returns>集めた証拠</returns>
        public static 継ぎ目の証拠 Get_証拠(IReadOnlyList<string> p_配列群, IReadOnlyList<(string A_リード1, string A_リード2)> p_ライブラリ群)
        {
            var l_マッパー = new ReadMapper(p_配列群);
            var l_深さの差分 = p_配列群.Select(x => new int[x.Length + 1]).ToArray();
            var l_一意な深さの差分 = p_配列群.Select(x => new int[x.Length + 1]).ToArray();
            var l_左の切れ端 = p_配列群.Select(x => new int[x.Length + 1]).ToArray();
            var l_右の切れ端 = p_配列群.Select(x => new int[x.Length + 1]).ToArray();

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_読み = Enumerable.Range(0, l_スレッド数).Select(_ => new List<(int, int, int)>()).ToArray();
            var l_断片 = Enumerable.Range(0, l_スレッド数).Select(_ => new List<(int, int, int)>()).ToArray();
            var l_外れ錨 = Enumerable.Range(0, l_スレッド数).Select(_ => new List<(int, int, int)>()).ToArray();

            foreach (var (l_リード1, l_リード2) in p_ライブラリ群)
            {
                var l_Isペア = !string.IsNullOrWhiteSpace(l_リード2);
                var l_供給元 = l_Isペア ? FastqReader.Get_ペア塩基列(l_リード1, l_リード2) : FastqReader.Get_生リード列(l_リード1).Select(x => (x, string.Empty));
                ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, l_供給元, (l_組, l_ワーカー番号) =>
                {
                    var (l_当たり1, l_当たり2) = Get_ペアの当たり(l_マッパー, l_組.Item1, l_組.Item2);
                    foreach (var l_当たり in new[] { l_当たり1, l_当たり2 })
                    {
                        if (l_当たり is { } l_有効)
                        {
                            V_加算_深さと切れ端(l_有効, l_深さの差分, l_一意な深さの差分, l_左の切れ端, l_右の切れ端);
                            l_読み[l_ワーカー番号].Add((l_有効.A_配列番号, l_有効.A_開始, l_有効.A_終了));
                        }
                    }

                    if (l_Isペア)
                    {
                        V_分類_ペア(l_当たり1, l_当たり2, l_断片[l_ワーカー番号], l_外れ錨[l_ワーカー番号]);
                    }
                });
            }

            return new 継ぎ目の証拠
            {
                A_配列群 = p_配列群,
                A_深さ = [.. l_深さの差分.Select(Get_累積)],
                A_一意な深さ = [.. l_一意な深さの差分.Select(Get_累積)],
                A_左の切れ端 = l_左の切れ端,
                A_右の切れ端 = l_右の切れ端,
                A_読み = Get_配列ごとに整列(l_読み, p_配列群.Count),
                A_断片 = Get_配列ごとに整列(l_断片, p_配列群.Count),
                A_外れ錨 = Get_配列ごとに整列(l_外れ錨, p_配列群.Count),
            };
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ペアの 2 本を当てる。向かい合わせに組める置き方があればそれを選び (反復に入った片方も相方の近くのコピーに置く)、無ければ 1 本ずつ最良の場所に置く
        /// </summary>
        /// <param name="p_マッパー"></param>
        /// <param name="p_リード1"></param>
        /// <param name="p_リード2">片側だけのライブラリなら空文字</param>
        /// <returns>2 本の当たり、当たらなかった側は null</returns>
        private static (リードの当たり? A_当たり1, リードの当たり? A_当たり2) Get_ペアの当たり(ReadMapper p_マッパー, string p_リード1, string p_リード2)
        {
            var l_候補1 = string.IsNullOrEmpty(p_リード1) ? [] : p_マッパー.Get_配置候補群(p_リード1);
            var l_候補2 = string.IsNullOrEmpty(p_リード2) ? [] : p_マッパー.Get_配置候補群(p_リード2);
            if (ReadMapper.Get_組んだ配置(l_候補1, l_候補2, C_断片長の上限) is { } l_組)
            {
                return (Get_当たり(l_組.A_配置1, p_リード1.Length), Get_当たり(l_組.A_配置2, p_リード2.Length));
            }

            return (Get_当たり(ReadMapper.Get_最良の配置(l_候補1), p_リード1.Length), Get_当たり(ReadMapper.Get_最良の配置(l_候補2), p_リード2.Length));
        }

        /// <summary>
        /// 配置を、配列上の整列範囲・両端の切れ端・向き・信頼度にまとめる
        /// </summary>
        /// <param name="p_配置"></param>
        /// <param name="p_リード長"></param>
        /// <returns>配置が無ければ null</returns>
        private static リードの当たり? Get_当たり(リード配置 p_配置, int p_リード長)
        {
            if (p_配置.A_配列番号 < 0 || p_配置.A_整列位置群.Count == 0)
            {
                return null;
            }

            var l_最初 = p_配置.A_整列位置群[0];
            var l_最後 = p_配置.A_整列位置群[^1];
            var l_最初の照合位置 = p_配置.A_Is逆鎖 ? p_リード長 - 1 - l_最初.A_リード位置 : l_最初.A_リード位置;
            var l_最後の照合位置 = p_配置.A_Is逆鎖 ? p_リード長 - 1 - l_最後.A_リード位置 : l_最後.A_リード位置;
            return new リードの当たり(p_配置.A_配列番号, l_最初.A_参照位置, l_最後.A_参照位置 + 1, l_最初の照合位置, p_リード長 - 1 - l_最後の照合位置, p_配置.A_Is逆鎖, p_配置.A_信頼度);
        }

        /// <summary>
        /// 当たったリードの深さと切れ端を加える
        /// </summary>
        /// <param name="p_当たり"></param>
        /// <param name="p_深さの差分"></param>
        /// <param name="p_一意な深さの差分"></param>
        /// <param name="p_左の切れ端"></param>
        /// <param name="p_右の切れ端"></param>
        private static void V_加算_深さと切れ端(リードの当たり p_当たり, int[][] p_深さの差分, int[][] p_一意な深さの差分, int[][] p_左の切れ端, int[][] p_右の切れ端)
        {
            var i = p_当たり.A_配列番号;
            _ = Interlocked.Increment(ref p_深さの差分[i][p_当たり.A_開始]);
            _ = Interlocked.Decrement(ref p_深さの差分[i][p_当たり.A_終了]);
            if (p_当たり.A_信頼度 >= C_一意とみなす信頼度)
            {
                _ = Interlocked.Increment(ref p_一意な深さの差分[i][p_当たり.A_開始]);
                _ = Interlocked.Decrement(ref p_一意な深さの差分[i][p_当たり.A_終了]);
            }

            if (p_当たり.A_左の切れ端 >= C_切れ端の最小長)
            {
                _ = Interlocked.Increment(ref p_左の切れ端[i][p_当たり.A_開始]);
            }

            if (p_当たり.A_右の切れ端 >= C_切れ端の最小長)
            {
                _ = Interlocked.Increment(ref p_右の切れ端[i][p_当たり.A_終了]);
            }
        }

        /// <summary>
        /// ペアの 2 本を、向き正しく組めた断片か、相方が正しく組めていない錨かに分ける
        /// </summary>
        /// <param name="p_当たり1"></param>
        /// <param name="p_当たり2"></param>
        /// <param name="p_断片">(配列番号, 断片の始まり, 断片の終わり) を足す先</param>
        /// <param name="p_外れ錨">(配列番号, 位置, 向き) を足す先</param>
        private static void V_分類_ペア(リードの当たり? p_当たり1, リードの当たり? p_当たり2, List<(int, int, int)> p_断片, List<(int, int, int)> p_外れ錨)
        {
            if (p_当たり1 is { } l_1 && p_当たり2 is { } l_2 && l_1.A_配列番号 == l_2.A_配列番号 && l_1.A_Is逆鎖 != l_2.A_Is逆鎖
                && Math.Min(l_1.A_信頼度, l_2.A_信頼度) >= C_一意とみなす信頼度)
            {
                var (l_左, l_右) = l_1.A_Is逆鎖 ? (l_2, l_1) : (l_1, l_2);
                if (l_左.A_開始 <= l_右.A_終了 && l_右.A_終了 - l_左.A_開始 < C_断片長の上限)
                {
                    p_断片.Add((l_左.A_配列番号, l_左.A_開始, l_右.A_終了));
                    return;
                }
            }

            foreach (var l_当たり in new[] { p_当たり1, p_当たり2 })
            {
                if (l_当たり is { } l_有効 && l_有効.A_信頼度 >= C_一意とみなす信頼度)
                {
                    p_外れ錨.Add(l_有効.A_Is逆鎖 ? (l_有効.A_配列番号, l_有効.A_終了, -1) : (l_有効.A_配列番号, l_有効.A_開始, 1));
                }
            }
        }

        /// <summary>
        /// 差分の累積で位置ごとの数にする
        /// </summary>
        /// <param name="p_差分">長さは配列長 + 1</param>
        /// <returns>長さは配列長</returns>
        private static int[] Get_累積(int[] p_差分)
        {
            var l_数 = new int[p_差分.Length - 1];
            var l_累積 = 0;
            for (var i = 0; i < l_数.Length; i++)
            {
                l_累積 += p_差分[i];
                l_数[i] = l_累積;
            }

            return l_数;
        }

        /// <summary>
        /// ワーカーごとの (配列番号, 値 1, 値 2) を配列ごとに分け、値の順に並べる
        /// </summary>
        /// <param name="p_ワーカーごと"></param>
        /// <param name="p_配列数"></param>
        /// <returns>配列ごとの (値 1, 値 2)</returns>
        private static (int, int)[][] Get_配列ごとに整列(List<(int, int, int)>[] p_ワーカーごと, int p_配列数)
        {
            var l_件数 = new int[p_配列数];
            foreach (var l_一覧 in p_ワーカーごと)
            {
                foreach (var (l_配列番号, _, _) in l_一覧)
                {
                    l_件数[l_配列番号]++;
                }
            }

            var l_結果 = l_件数.Select(x => new (int, int)[x]).ToArray();
            Array.Clear(l_件数);
            foreach (var l_一覧 in p_ワーカーごと)
            {
                foreach (var (l_配列番号, l_値1, l_値2) in l_一覧)
                {
                    l_結果[l_配列番号][l_件数[l_配列番号]++] = (l_値1, l_値2);
                }

                l_一覧.Clear();
            }

            foreach (var l_配列 in l_結果)
            {
                Array.Sort(l_配列);
            }

            return l_結果;
        }

        #endregion
    }
}
