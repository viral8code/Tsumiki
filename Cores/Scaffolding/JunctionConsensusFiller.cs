using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 最終配列の継ぎ目 (未確認の繋ぎ目の印と N の連続) を、両側の片を跨ぐリードの品質を見た多数決で埋める<br/>
    /// 両側の片に置けるリードを間の長さで票にし、最多の長さのリードの確かな塩基で列ごとに合意を取る
    /// </summary>
    internal static class JunctionConsensusFiller
    {
        #region 定数

        /// <summary>
        /// 対象にする N の連続の長さの上限
        /// </summary>
        public const int C_対象のNの長さ上限 = 500;

        /// <summary>
        /// 継ぎ目の両側から片として見る長さの上限
        /// </summary>
        private const int C_片の長さ = 300;

        /// <summary>
        /// 両側の片がこれより短い継ぎ目は対象にしない
        /// </summary>
        private const int C_片の最短長 = 60;

        /// <summary>
        /// 種を取る、片の継ぎ目側の長さ
        /// </summary>
        private const int C_種を取る長さ = 120;

        /// <summary>
        /// 種の長さ
        /// </summary>
        private const int C_種長 = 21;

        /// <summary>
        /// 1 つの継ぎ目に集めるリードの上限 (読み込み順)
        /// </summary>
        private const int C_集めるリードの上限 = 1_000;

        /// <summary>
        /// 確かな塩基とみなす品質の下限 (Phred)
        /// </summary>
        private const int C_確かな品質 = 20;

        /// <summary>
        /// 片に置くときに許す、確かな塩基の不一致の数
        /// </summary>
        private const int C_許す不一致の数 = 2;

        /// <summary>
        /// 最多の間の長さに要るリードの数
        /// </summary>
        private const int C_要るリード数 = 3;

        /// <summary>
        /// 最多の間の長さのリード数が、次点の何倍以上なら採るか
        /// </summary>
        private const int C_首位の優勢比 = 3;

        /// <summary>
        /// 合意で置き直す、片の継ぎ目側の長さの下限
        /// </summary>
        private const int C_余白の最短 = 20;

        /// <summary>
        /// 両側の片が重なるとき、重なりの半分に足す余白
        /// </summary>
        private const int C_重なりに足す余白 = 10;

        /// <summary>
        /// 品質の合計で列を決めるときに数える、塩基の品質の下限 (Phred)
        /// </summary>
        private const int C_低い品質の下限 = 5;

        /// <summary>
        /// 品質の合計で列を決めるときに、最多の塩基に要る品質の合計
        /// </summary>
        private const int C_品質の合計の下限 = 40;

        /// <summary>
        /// 合意の列に要る、最多の確かな塩基の数
        /// </summary>
        private const int C_列に要る本数 = 2;

        /// <summary>
        /// 合意の列で、最多の確かな塩基が占める割合の下限
        /// </summary>
        private const double C_合意の割合 = 0.8D;

        /// <summary>
        /// 縦に並んだ反復を調べる、合意の区間の外側の長さ
        /// </summary>
        private const int C_反復を調べる外側の長さ = 25;

        /// <summary>
        /// 片側だけに置けるリードを、繋いだ形のもう一方の片と比べるのに要る重なりの長さ
        /// </summary>
        private const int C_食い違いを数える重なり = 30;

        /// <summary>
        /// 繋いだ形と食い違う片側だけのリードがこの数以上なら、埋めるのを見送る候補にする
        /// </summary>
        private const int C_食い違うリードの下限 = 10;

        /// <summary>
        /// 食い違う片側だけのリードが、合うもの + 1 の何倍以上なら埋めないか
        /// </summary>
        private const int C_食い違いの優勢比 = 3;

        /// <summary>
        /// 種の照合を並列にする 1 まとまりのリード数
        /// </summary>
        private const int C_照合のバッチサイズ = 16_384;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群の継ぎ目を、跨ぐリードの合意で埋める
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_リード群">リードの配列・品質・Phred オフセットを流す関数 (継ぎ目があるときだけ 1 回呼ぶ)</param>
        /// <param name="p_統計">結果の統計</param>
        /// <returns>埋めた後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_埋めた配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, Func<IEnumerable<(string A_配列, string A_品質, int A_Phredオフセット)>> p_リード群, out 継ぎ目合意統計 p_統計)
        {
            var l_継ぎ目群 = Get_継ぎ目群(p_配列群);
            if (l_継ぎ目群.Count == 0)
            {
                p_統計 = default;
                return [.. p_配列群];
            }

            var l_候補群 = Get_候補リード群(l_継ぎ目群, p_リード群());
            var l_埋め方群 = new 継ぎ目の埋め方?[l_継ぎ目群.Count];
            var l_判定群 = new 継ぎ目合意判定[l_継ぎ目群.Count];
            var l_内訳群 = new 継ぎ目合意の内訳[l_継ぎ目群.Count];
            var l_並列設定 = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数)
            };
            _ = Parallel.For(0, l_継ぎ目群.Count, l_並列設定, j =>
            {
                l_埋め方群[j] = Get_埋め方(l_継ぎ目群[j], l_候補群[j], out var l_判定, out var l_内訳);
                l_判定群[j] = l_判定;
                l_内訳群[j] = l_内訳;
            });
            for (var j = 0; j < l_継ぎ目群.Count; j++)
            {
                var l_継ぎ目 = l_継ぎ目群[j];
                var l_内訳 = l_内訳群[j];
                Logger.V_出力(メッセージID.継ぎ目の合意の内訳, p_配列群[l_継ぎ目.A_配列番号].A_ID, l_継ぎ目.A_開始, l_継ぎ目.A_左の片.Length, l_継ぎ目.A_右の片.Length, l_候補群[j].Count, l_内訳.A_跨ぐ数, l_内訳.A_間, l_内訳.A_食い違う数, l_内訳.A_合う数, l_判定群[j]);
            }

            var l_結果 = Get_書き換えた配列群(p_配列群, l_継ぎ目群, l_埋め方群, out var l_埋めた数);
            p_統計 = new 継ぎ目合意統計(l_継ぎ目群.Count, l_埋めた数, l_判定群.Count(static x => x == 継ぎ目合意判定.跨ぐリード無し), l_判定群.Count(static x => x == 継ぎ目合意判定.票が足りない), l_判定群.Count(static x => x == 継ぎ目合意判定.票が割れた), l_判定群.Count(static x => x == 継ぎ目合意判定.合意が弱い), l_判定群.Count(static x => x == 継ぎ目合意判定.縦列反復), l_判定群.Count(static x => x == 継ぎ目合意判定.片側のリードが食い違う), l_判定群.Count(static x => x == 継ぎ目合意判定.片が重なりに収まる));
            return l_結果;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 両側の片が十分に長い継ぎ目を、配列の順・位置の順に集める
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <returns>継ぎ目の一覧</returns>
        private static List<合意の継ぎ目> Get_継ぎ目群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群)
        {
            List<合意の継ぎ目> l_結果 = [];
            for (var s = 0; s < p_配列群.Count; s++)
            {
                var l_配列 = p_配列群[s].A_配列;
                var l_前の終わり = 0;
                var i = 0;
                while (i < l_配列.Length)
                {
                    if (!Util.Isギャップ文字(l_配列[i]))
                    {
                        i++;
                        continue;
                    }

                    var l_終わり = Util.Get_ギャップの終わり(l_配列, i);
                    var l_次の始まり = l_終わり;
                    while (l_次の始まり < l_配列.Length && !Util.Isギャップ文字(l_配列[l_次の始まり]))
                    {
                        l_次の始まり++;
                    }

                    var l_左の長さ = Math.Min(C_片の長さ, i - l_前の終わり);
                    var l_右の長さ = Math.Min(C_片の長さ, l_次の始まり - l_終わり);
                    var l_Is対象の長さ = l_配列[i] == Consts.未確認の繋ぎ目 || l_終わり - i <= C_対象のNの長さ上限;
                    if (l_Is対象の長さ && l_左の長さ >= C_片の最短長 && l_右の長さ >= C_片の最短長)
                    {
                        l_結果.Add(new 合意の継ぎ目(s, i, l_終わり, l_配列.Substring(i - l_左の長さ, l_左の長さ), l_配列.Substring(l_終わり, l_右の長さ)));
                    }

                    l_前の終わり = l_終わり;
                    i = l_終わり;
                }
            }

            return l_結果;
        }

        /// <summary>
        /// リードを 1 回読み、種に当たるリードを継ぎ目ごとに集める
        /// </summary>
        /// <param name="p_継ぎ目群">継ぎ目の一覧</param>
        /// <param name="p_リード群">リードの配列・品質・Phred オフセット</param>
        /// <returns>継ぎ目ごとの候補のリード (配列と Phred 値)</returns>
        private static List<(string A_配列, byte[] A_品質)>[] Get_候補リード群(IReadOnlyList<合意の継ぎ目> p_継ぎ目群, IEnumerable<(string A_配列, string A_品質, int A_Phredオフセット)> p_リード群)
        {
            var l_種表 = Get_集める種表(p_継ぎ目群);
            var l_候補群 = Enumerable.Range(0, p_継ぎ目群.Count).Select(static _ => new List<(string A_配列, byte[] A_品質)>()).ToArray();
            var l_並列設定 = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数)
            };
            var l_バッチ = new List<(string A_配列, string A_品質, int A_Phredオフセット)>(C_照合のバッチサイズ);
            foreach (var l_リード in p_リード群)
            {
                l_バッチ.Add(l_リード);
                if (l_バッチ.Count == C_照合のバッチサイズ)
                {
                    V_振り分け_候補(l_バッチ, l_種表, l_候補群, l_並列設定);
                    l_バッチ.Clear();
                }
            }

            V_振り分け_候補(l_バッチ, l_種表, l_候補群, l_並列設定);
            return l_候補群;
        }

        /// <summary>
        /// 全継ぎ目の両側の片の継ぎ目側から、種の正準キーと、それを持つ継ぎ目の番号の表を作る
        /// </summary>
        /// <param name="p_継ぎ目群">継ぎ目の一覧</param>
        /// <returns>種の正準キーから継ぎ目の番号の一覧を引く表</returns>
        private static Dictionary<UInt128, List<int>> Get_集める種表(IReadOnlyList<合意の継ぎ目> p_継ぎ目群)
        {
            Dictionary<UInt128, List<int>> l_表 = [];
            for (var j = 0; j < p_継ぎ目群.Count; j++)
            {
                var l_左 = p_継ぎ目群[j].A_左の片;
                var l_右 = p_継ぎ目群[j].A_右の片;
                foreach (var l_窓 in new[] { l_左[Math.Max(0, l_左.Length - C_種を取る長さ)..], l_右[..Math.Min(l_右.Length, C_種を取る長さ)], })
                {
                    var l_転がし = new RollingKmer(C_種長);
                    foreach (var l_塩基 in l_窓)
                    {
                        if (!l_転がし.Is成功_追加(l_塩基, out var l_キー))
                        {
                            continue;
                        }

                        if (!l_表.TryGetValue(l_キー.A_下位, out var l_番号群))
                        {
                            l_番号群 = [];
                            l_表[l_キー.A_下位] = l_番号群;
                        }

                        if (l_番号群.Count == 0 || l_番号群[^1] != j)
                        {
                            l_番号群.Add(j);
                        }
                    }
                }
            }

            return l_表;
        }

        /// <summary>
        /// 1 まとまりのリードの種の照合を並列に行い、当たった継ぎ目の候補へ読み込み順に足す
        /// </summary>
        /// <param name="p_バッチ">リードの配列・品質・Phred オフセット</param>
        /// <param name="p_種表">種の正準キーから継ぎ目の番号の一覧を引く表</param>
        /// <param name="p_候補群">継ぎ目ごとの候補のリード (足す先)</param>
        /// <param name="p_並列設定">並列の設定</param>
        private static void V_振り分け_候補(List<(string A_配列, string A_品質, int A_Phredオフセット)> p_バッチ, Dictionary<UInt128, List<int>> p_種表, List<(string A_配列, byte[] A_品質)>[] p_候補群, ParallelOptions p_並列設定)
        {
            var l_当たり群 = new List<int>?[p_バッチ.Count];
            _ = Parallel.For(0, p_バッチ.Count, p_並列設定, i => l_当たり群[i] = Get_当たる継ぎ目群(p_バッチ[i].A_配列, p_種表));
            for (var i = 0; i < p_バッチ.Count; i++)
            {
                if (l_当たり群[i] is not { } l_番号群)
                {
                    continue;
                }

                byte[]? l_品質 = null;
                foreach (var j in l_番号群)
                {
                    if (p_候補群[j].Count >= C_集めるリードの上限)
                    {
                        continue;
                    }

                    l_品質 ??= Get_Phred値(p_バッチ[i].A_品質, p_バッチ[i].A_Phredオフセット, p_バッチ[i].A_配列.Length);
                    p_候補群[j].Add((p_バッチ[i].A_配列, l_品質));
                }
            }
        }

        /// <summary>
        /// リードの種が当たる継ぎ目の番号を返す
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_種表">種の正準キーから継ぎ目の番号の一覧を引く表</param>
        /// <returns>当たった継ぎ目の番号、無ければ null</returns>
        private static List<int>? Get_当たる継ぎ目群(string p_リード, Dictionary<UInt128, List<int>> p_種表)
        {
            List<int>? l_結果 = null;
            var l_転がし = new RollingKmer(C_種長);
            foreach (var l_塩基 in p_リード)
            {
                if (!l_転がし.Is成功_追加(l_塩基, out var l_キー) || !p_種表.TryGetValue(l_キー.A_下位, out var l_番号群))
                {
                    continue;
                }

                l_結果 ??= [];
                foreach (var j in l_番号群)
                {
                    if (!l_結果.Contains(j))
                    {
                        l_結果.Add(j);
                    }
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 品質の文字列を Phred 値の列にする (足りない分は 0)
        /// </summary>
        /// <param name="p_品質">品質の文字列</param>
        /// <param name="p_オフセット">Phred オフセット</param>
        /// <param name="p_長さ">リードの長さ</param>
        /// <returns>Phred 値の列</returns>
        private static byte[] Get_Phred値(string p_品質, int p_オフセット, int p_長さ)
        {
            var l_値 = new byte[p_長さ];
            for (var i = 0; i < p_長さ && i < p_品質.Length; i++)
            {
                l_値[i] = (byte)Math.Clamp(p_品質[i] - p_オフセット, 0, byte.MaxValue);
            }

            return l_値;
        }

        /// <summary>
        /// 1 つの継ぎ目を、跨ぐリードの間の長さの票と、列ごとの確かな塩基の多数決で埋める方法を決める
        /// </summary>
        /// <param name="p_継ぎ目">継ぎ目</param>
        /// <param name="p_候補">候補のリード (配列と Phred 値)</param>
        /// <param name="p_判定">判定</param>
        /// <param name="p_内訳">判定に使ったリードの数</param>
        /// <returns>埋め方、埋めなければ null</returns>
        private static 継ぎ目の埋め方? Get_埋め方(合意の継ぎ目 p_継ぎ目, List<(string A_配列, byte[] A_品質)> p_候補, out 継ぎ目合意判定 p_判定, out 継ぎ目合意の内訳 p_内訳)
        {
            var l_左 = p_継ぎ目.A_左の片;
            var l_右 = p_継ぎ目.A_右の片;
            var l_跨ぐリード群 = Get_跨ぐリード群(l_左, l_右, p_候補, out var l_片側だけのリード群);
            p_内訳 = new 継ぎ目合意の内訳(l_跨ぐリード群.Count, null, 0, 0);
            if (l_跨ぐリード群.Count == 0)
            {
                p_判定 = 継ぎ目合意判定.跨ぐリード無し;
                return null;
            }

            var l_票 = l_跨ぐリード群.GroupBy(static x => x.A_間).Select(static x => (A_間: x.Key, A_数: x.Count())).OrderByDescending(static x => x.A_数).ThenBy(static x => x.A_間).ToList();
            var l_次点 = l_票.Count > 1 ? l_票[1].A_数 : 0;
            if (l_票[0].A_数 < C_要るリード数)
            {
                p_判定 = 継ぎ目合意判定.票が足りない;
                return null;
            }

            if (l_票[0].A_数 < C_首位の優勢比 * l_次点)
            {
                p_判定 = 継ぎ目合意判定.票が割れた;
                return null;
            }

            var l_間 = l_票[0].A_間;
            if (-l_間 >= Math.Min(l_左.Length, l_右.Length))
            {
                p_内訳 = new 継ぎ目合意の内訳(l_跨ぐリード群.Count, l_間, 0, 0);
                p_判定 = 継ぎ目合意判定.片が重なりに収まる;
                return null;
            }

            var l_食い違う数 = Get_食い違うリード数(l_左, l_右, l_間, l_片側だけのリード群, out var l_合う数);
            p_内訳 = new 継ぎ目合意の内訳(l_跨ぐリード群.Count, l_間, l_食い違う数, l_合う数);
            if (l_食い違う数 >= C_食い違うリードの下限 && l_食い違う数 >= C_食い違いの優勢比 * (l_合う数 + 1))
            {
                p_判定 = 継ぎ目合意判定.片側のリードが食い違う;
                return null;
            }

            var l_余白 = Math.Max(C_余白の最短, ((-l_間 + 1) / 2) + C_重なりに足す余白);
            if (l_余白 > l_左.Length - C_種長 || l_余白 > l_右.Length - C_種長 || Get_合意(l_跨ぐリード群.Where(x => x.A_間 == l_間), l_余白, (2 * l_余白) + l_間) is not { } l_合意)
            {
                p_判定 = 継ぎ目合意判定.合意が弱い;
                return null;
            }

            var l_合意の始まり = l_左.Length - l_余白;
            var l_繋いだ配列 = string.Concat(l_左.AsSpan(0, l_合意の始まり), l_合意, l_右.AsSpan(l_余白));
            if (Scaffolder.Is縦の反復の中(l_繋いだ配列, Math.Max(0, l_合意の始まり - C_反復を調べる外側の長さ), l_合意の始まり + l_合意.Length))
            {
                p_判定 = 継ぎ目合意判定.縦列反復;
                return null;
            }

            p_判定 = 継ぎ目合意判定.埋めた;
            return new 継ぎ目の埋め方(l_余白, l_合意, l_余白);
        }

        /// <summary>
        /// 候補のリード (両方の向き) のうち、左右の両方の片に置けるものを、間の長さ (負なら重なり) と共に返す<br/>
        /// どちらの向きでも両方には置けず、片側にだけ置けたリードは、最初に置けた向きと片で別に返す
        /// </summary>
        /// <param name="p_左">左の片</param>
        /// <param name="p_右">右の片</param>
        /// <param name="p_候補">候補のリード (配列と Phred 値)</param>
        /// <param name="p_片側だけのリード群">片側にだけ置けたリードの、向きを揃えた配列と Phred 値・置けたのが左の片か・置く差</param>
        /// <returns>跨ぐリードの、間の長さ・向きを揃えた配列と Phred 値・リード上の左の片の終わりの位置</returns>
        private static List<(int A_間, string A_配列, byte[] A_品質, int A_左の終わり)> Get_跨ぐリード群(string p_左, string p_右, List<(string A_配列, byte[] A_品質)> p_候補, out List<(string A_配列, byte[] A_品質, bool A_Is左, int A_差)> p_片側だけのリード群)
        {
            var l_左の種表 = Get_片の種表(p_左, p_左.Length - Math.Min(p_左.Length, C_種を取る長さ), p_左.Length).GetAlternateLookup<ReadOnlySpan<char>>();
            var l_右の種表 = Get_片の種表(p_右, 0, Math.Min(p_右.Length, C_種を取る長さ)).GetAlternateLookup<ReadOnlySpan<char>>();
            List<(int A_間, string A_配列, byte[] A_品質, int A_左の終わり)> l_結果 = [];
            p_片側だけのリード群 = [];
            foreach (var (l_配列, l_品質) in p_候補)
            {
                (string A_配列, byte[] A_品質, bool A_Is左, int A_差)? l_片側 = null;
                var l_Is跨ぐ = false;
                foreach (var l_Is逆 in new[] { false, true, })
                {
                    var l_向きの配列 = l_Is逆 ? Util.V_逆相補_曖昧塩基あり(l_配列) : l_配列;
                    var l_向きの品質 = l_Is逆 ? [.. l_品質.AsEnumerable().Reverse()] : l_品質;
                    var l_左の差 = Get_置き方(l_向きの配列, l_向きの品質, p_左, l_左の種表);
                    var l_右の差 = Get_置き方(l_向きの配列, l_向きの品質, p_右, l_右の種表);
                    if (l_左の差 is { } l_左 && l_右の差 is { } l_右)
                    {
                        var l_左の終わり = p_左.Length - l_左;
                        l_結果.Add((-l_右 - l_左の終わり, l_向きの配列, l_向きの品質, l_左の終わり));
                        l_Is跨ぐ = true;
                        break;
                    }

                    l_片側 ??= l_左の差 is { } l_左だけ ? (l_向きの配列, l_向きの品質, true, l_左だけ) : l_右の差 is { } l_右だけ ? (l_向きの配列, l_向きの品質, false, l_右だけ) : null;
                }

                if (!l_Is跨ぐ && l_片側 is { } l_片側のリード)
                {
                    p_片側だけのリード群.Add(l_片側のリード);
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 片側だけに置けたリードを、間の長さ p_間 で繋いだ形のもう一方の片に当て、重なりが <see cref="C_食い違いを数える重なり"/> 以上のものの食い違いを数える<br/>
        /// 確かな塩基の不一致が <see cref="C_許す不一致の数"/> を超えれば食い違う、それ以下なら合うとする
        /// </summary>
        /// <param name="p_左">左の片</param>
        /// <param name="p_右">右の片</param>
        /// <param name="p_間">間の長さ (負なら重なり)</param>
        /// <param name="p_片側だけのリード群">片側にだけ置けたリードの、向きを揃えた配列と Phred 値・置けたのが左の片か・置く差</param>
        /// <param name="p_合う数">繋いだ形と合うリードの数</param>
        /// <returns>繋いだ形と食い違うリードの数</returns>
        private static int Get_食い違うリード数(string p_左, string p_右, int p_間, List<(string A_配列, byte[] A_品質, bool A_Is左, int A_差)> p_片側だけのリード群, out int p_合う数)
        {
            var l_食い違う数 = 0;
            p_合う数 = 0;
            foreach (var (l_配列, l_品質, l_Is左, l_差) in p_片側だけのリード群)
            {
                var l_不一致 = l_Is左 ? Get_確かな不一致数(l_配列, l_品質, p_右, l_差 - p_左.Length - p_間, out var l_重なり) : Get_確かな不一致数(l_配列, l_品質, p_左, l_差 + p_左.Length + p_間, out l_重なり);
                if (l_重なり < C_食い違いを数える重なり)
                {
                    continue;
                }

                if (l_不一致 > C_許す不一致の数)
                {
                    l_食い違う数++;
                }
                else
                {
                    p_合う数++;
                }
            }

            return l_食い違う数;
        }

        /// <summary>
        /// 片の [p_始め, p_終わり) にある種のうち、その範囲で 1 回だけ出て ACGT だけからなるものの、片の上の位置の表を作る
        /// </summary>
        /// <param name="p_片">片の配列</param>
        /// <param name="p_始め">範囲の始め</param>
        /// <param name="p_終わり">範囲の終わり (含まない)</param>
        /// <returns>種から片の上の位置を引く表</returns>
        private static Dictionary<string, int> Get_片の種表(string p_片, int p_始め, int p_終わり)
        {
            Dictionary<string, int> l_表 = new(StringComparer.Ordinal);
            HashSet<string> l_重複 = new(StringComparer.Ordinal);
            for (var i = p_始め; i + C_種長 <= p_終わり; i++)
            {
                var l_種 = p_片.Substring(i, C_種長);
                if (l_種.AsSpan().IndexOfAnyExcept(Consts.塩基文字) >= 0)
                {
                    continue;
                }

                if (!l_表.TryAdd(l_種, i))
                {
                    _ = l_重複.Add(l_種);
                }
            }

            foreach (var l_種 in l_重複)
            {
                _ = l_表.Remove(l_種);
            }

            return l_表;
        }

        /// <summary>
        /// リードを片に置く差 (片の位置 − リードの位置) を、種の当たりの票で決め、確かな塩基の不一致が少なければ返す
        /// </summary>
        /// <param name="p_リード">向きを揃えたリードの配列</param>
        /// <param name="p_品質">向きを揃えたリードの Phred 値</param>
        /// <param name="p_片">片の配列</param>
        /// <param name="p_種表">種から片の上の位置を引く表</param>
        /// <returns>置く差、置けなければ null</returns>
        private static int? Get_置き方(string p_リード, byte[] p_品質, string p_片, Dictionary<string, int>.AlternateLookup<ReadOnlySpan<char>> p_種表)
        {
            Dictionary<int, int> l_差の票 = [];
            for (var i = 0; i + C_種長 <= p_リード.Length; i++)
            {
                if (p_種表.TryGetValue(p_リード.AsSpan(i, C_種長), out var l_位置))
                {
                    l_差の票[l_位置 - i] = l_差の票.GetValueOrDefault(l_位置 - i) + 1;
                }
            }

            if (l_差の票.Count == 0)
            {
                return null;
            }

            var l_並び = l_差の票.OrderByDescending(static x => x.Value).ThenBy(static x => x.Key).ToList();
            if (l_並び.Count > 1 && l_並び[1].Value * 2 >= l_並び[0].Value)
            {
                return null;
            }

            var l_差 = l_並び[0].Key;
            var l_不一致 = Get_確かな不一致数(p_リード, p_品質, p_片, l_差, out var l_重なり);
            return l_重なり >= C_種長 && l_不一致 <= C_許す不一致の数 ? l_差 : null;
        }

        /// <summary>
        /// リードを差 p_差 (片の位置 − リードの位置) で片に重ね、重なりの長さと、確かな塩基の不一致の数を返す
        /// </summary>
        /// <param name="p_リード">向きを揃えたリードの配列</param>
        /// <param name="p_品質">向きを揃えたリードの Phred 値</param>
        /// <param name="p_片">片の配列</param>
        /// <param name="p_差">置く差</param>
        /// <param name="p_重なり">リードと片が重なる長さ</param>
        /// <returns>重なりの中の、品質が <see cref="C_確かな品質"/> 以上の塩基の不一致の数</returns>
        private static int Get_確かな不一致数(string p_リード, byte[] p_品質, string p_片, int p_差, out int p_重なり)
        {
            p_重なり = 0;
            var l_不一致 = 0;
            for (var i = Math.Max(0, -p_差); i < p_リード.Length && i + p_差 < p_片.Length; i++)
            {
                p_重なり++;
                if (p_リード[i] != p_片[i + p_差] && p_品質[i] >= C_確かな品質)
                {
                    l_不一致++;
                }
            }

            return l_不一致;
        }

        /// <summary>
        /// 跨ぐリードの [左の片の終わり − p_余白, 左の片の終わり − p_余白 + p_長さ) を列ごとに数え、多数決で配列を決める<br/>
        /// 列ごとに、確かな塩基の票で決まればその塩基、決まらなければ品質の合計で決める (<see cref="Get_列の塩基"/>)
        /// </summary>
        /// <param name="p_リード群">間の長さが同じ跨ぐリード</param>
        /// <param name="p_余白">左の片の終わりから手前に置き直す長さ</param>
        /// <param name="p_長さ">決める配列の長さ</param>
        /// <returns>決めた配列、決まらない列があれば null</returns>
        private static string? Get_合意(IEnumerable<(int A_間, string A_配列, byte[] A_品質, int A_左の終わり)> p_リード群, int p_余白, int p_長さ)
        {
            if (p_長さ <= 0)
            {
                return null;
            }

            var l_票 = new int[p_長さ, Consts.塩基文字.Length];
            var l_品質の合計 = new int[p_長さ, Consts.塩基文字.Length];
            foreach (var (_, l_配列, l_品質, l_左の終わり) in p_リード群)
            {
                for (var t = 0; t < p_長さ; t++)
                {
                    var i = l_左の終わり - p_余白 + t;
                    if (i < 0 || i >= l_配列.Length || l_品質[i] < C_低い品質の下限)
                    {
                        continue;
                    }

                    var l_番号 = Consts.塩基文字.IndexOf(l_配列[i]);
                    if (l_番号 < 0)
                    {
                        continue;
                    }

                    l_品質の合計[t, l_番号] += l_品質[i];
                    if (l_品質[i] >= C_確かな品質)
                    {
                        l_票[t, l_番号]++;
                    }
                }
            }

            var l_合意 = new StringBuilder(p_長さ);
            for (var t = 0; t < p_長さ; t++)
            {
                if (Get_列の塩基(l_票, l_品質の合計, t) is not { } l_番号)
                {
                    return null;
                }

                _ = l_合意.Append(Consts.塩基文字[l_番号]);
            }

            return l_合意.ToString();
        }

        /// <summary>
        /// 合意の 1 列の塩基を決める<br/>
        /// 確かな塩基の最多が <see cref="C_列に要る本数"/> 以上かつ <see cref="C_合意の割合"/> 以上ならその塩基、そうでなければ、品質の合計の最多が <see cref="C_品質の合計の下限"/> 以上かつ <see cref="C_合意の割合"/> 以上で、ほかの塩基に確かな票が無ければその塩基
        /// </summary>
        /// <param name="p_票">列・塩基ごとの確かな塩基の数</param>
        /// <param name="p_品質の合計">列・塩基ごとの、品質が <see cref="C_低い品質の下限"/> 以上の塩基の品質の合計</param>
        /// <param name="p_列">列</param>
        /// <returns>塩基の番号 (<see cref="Consts.塩基文字"/> の位置)、決まらなければ null</returns>
        private static int? Get_列の塩基(int[,] p_票, int[,] p_品質の合計, int p_列)
        {
            var l_票の最多 = Get_最多の番号(p_票, p_列, out var l_票の合計);
            if (p_票[p_列, l_票の最多] >= C_列に要る本数 && p_票[p_列, l_票の最多] >= C_合意の割合 * l_票の合計)
            {
                return l_票の最多;
            }

            var l_合計の最多 = Get_最多の番号(p_品質の合計, p_列, out var l_品質の総和);
            var l_最多の合計 = p_品質の合計[p_列, l_合計の最多];
            return l_最多の合計 >= C_品質の合計の下限 && l_最多の合計 >= C_合意の割合 * l_品質の総和 && l_票の合計 == p_票[p_列, l_合計の最多] ? l_合計の最多 : null;
        }

        /// <summary>
        /// 列の中で値が最多の塩基の番号を返す (同じなら番号の小さいもの)
        /// </summary>
        /// <param name="p_表">列・塩基ごとの値</param>
        /// <param name="p_列">列</param>
        /// <param name="p_合計">その列の値の合計</param>
        /// <returns>最多の塩基の番号</returns>
        private static int Get_最多の番号(int[,] p_表, int p_列, out int p_合計)
        {
            var l_最多 = 0;
            p_合計 = 0;
            for (var b = 0; b < Consts.塩基文字.Length; b++)
            {
                p_合計 += p_表[p_列, b];
                if (p_表[p_列, b] > p_表[p_列, l_最多])
                {
                    l_最多 = b;
                }
            }

            return l_最多;
        }

        /// <summary>
        /// 埋め方の決まった継ぎ目を配列に当てる (前の継ぎ目で削った所と重なる継ぎ目は当てない)
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_継ぎ目群">継ぎ目の一覧 (配列の順・位置の順)</param>
        /// <param name="p_埋め方群">継ぎ目ごとの埋め方</param>
        /// <param name="p_埋めた数">当てた継ぎ目の数</param>
        /// <returns>書き換えた配列の ID と配列</returns>
        private static List<(string A_ID, string A_配列)> Get_書き換えた配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<合意の継ぎ目> p_継ぎ目群, 継ぎ目の埋め方?[] p_埋め方群, out int p_埋めた数)
        {
            p_埋めた数 = 0;
            var l_配列ごと = Enumerable.Range(0, p_継ぎ目群.Count).Where(x => p_埋め方群[x] is not null).GroupBy(x => p_継ぎ目群[x].A_配列番号).ToDictionary(static x => x.Key, static x => x.ToList());
            List<(string A_ID, string A_配列)> l_結果 = [];
            for (var s = 0; s < p_配列群.Count; s++)
            {
                var (l_ID, l_配列) = p_配列群[s];
                if (!l_配列ごと.TryGetValue(s, out var l_番号群))
                {
                    l_結果.Add((l_ID, l_配列));
                    continue;
                }

                var l_出力 = new StringBuilder(l_配列.Length);
                var l_読んだ位置 = 0;
                foreach (var j in l_番号群)
                {
                    var l_継ぎ目 = p_継ぎ目群[j];
                    var l_埋め方 = p_埋め方群[j]!.Value;
                    var l_左の端 = l_継ぎ目.A_開始 - l_埋め方.A_左から削る長さ;
                    if (l_左の端 < l_読んだ位置)
                    {
                        continue;
                    }

                    _ = l_出力.Append(l_配列, l_読んだ位置, l_左の端 - l_読んだ位置).Append(l_埋め方.A_埋める配列);
                    l_読んだ位置 = l_継ぎ目.A_終わり + l_埋め方.A_右から削る長さ;
                    p_埋めた数++;
                }

                _ = l_出力.Append(l_配列, l_読んだ位置, l_配列.Length - l_読んだ位置);
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            return l_結果;
        }

        #endregion
    }
}
