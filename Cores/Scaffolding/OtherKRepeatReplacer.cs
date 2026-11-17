using System.Runtime.InteropServices;
using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 最終配列のうち、ほかの k の contig と錨の間の長さが食い違う所を、リードに続けて出てくるほうの形に差し替える<br/>
    /// 錨は、最終配列の中で両鎖を通して 1 回だけ出て、ほかの k の contig の中でも 1 回だけ出る C_錨長 の k-mer
    /// </summary>
    internal static class OtherKRepeatReplacer
    {
        #region 定数

        /// <summary>
        /// 錨の長さ<br/>
        /// 2 bit で 64 bit に収まり、逆相補と同じにならない奇数
        /// </summary>
        internal const int C_錨長 = 31;

        /// <summary>
        /// 差し替える区間の長さの上限<br/>
        /// 最終配列の側とほかの k の側のどちらも、両端の錨の始まりの間でこの長さまで
        /// </summary>
        private const int C_区間の長さの上限 = 2000;

        /// <summary>
        /// 差し替える区間の、最終配列とほかの k の長さの差の下限<br/>
        /// これより小さい食い違いは、ポリッシュに任せる
        /// </summary>
        private const int C_長さの差の下限 = 11;

        /// <summary>
        /// 隣り合う食い違いを 1 つの区間にまとめる、間の錨の並びの長さの上限<br/>
        /// 縦の反復の写しどうしの違いで、反復の中にも錨が残る所を 1 つにまとめる
        /// </summary>
        private const int C_まとめる間の上限 = 100;

        /// <summary>
        /// リードで確かめる窓の長さの、リード長に対する比の分子<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分子 = 7;

        /// <summary>
        /// リードで確かめる窓の長さの、リード長に対する比の分母<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分母 = 10;

        /// <summary>
        /// 差し替える形を確かめるときに要る、リードの出現数の下限<br/>
        /// PairBridgeFiller.Is写しがリードに続けて出る の窓の支持数と同じ
        /// </summary>
        private const int C_支持数の下限 = 2;

        /// <summary>
        /// 差し替えた所のログに出す、両脇の配列の長さの上限<br/>
        /// 評価のあとで、最終配列の中の場所を探すのに使う
        /// </summary>
        private const int C_ログの脇の長さ = 60;

        /// <summary>
        /// 錨の表で、2 回以上出る k-mer の印
        /// </summary>
        private const long C_一意でない = -1L;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群のうち、ほかの k の contig と錨の間の長さが食い違う区間を、リードに続けて出てくるほかの k の形に差し替える<br/>
        /// 最終配列の形がリードに続けて出てくる区間、ほかの k の形がただ 1 つに決まらない区間、ほかの k の形の両端の錨の間の長さをリードで確かめられない区間は、そのまま残す
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_ほかのkの配列群">ほかの k の contig と scaffold の配列 (大文字)</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_統計">差し替えた結果の統計</param>
        /// <returns>差し替えた後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_差し替えた配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<string> p_ほかのkの配列群, ReadMinimizerIndex p_リード索引, int p_リード長, out ほかのkで差し替えた統計 p_統計)
        {
            var l_錨表 = Get_一意な錨の表(p_配列群);
            List<(int A_配列番号, int A_始め, int A_終わり, string A_差し替える配列)> l_食い違い群 = [];
            foreach (var l_ほかの配列 in p_ほかのkの配列群)
            {
                V_集める_食い違い(l_ほかの配列, l_錨表, p_配列群, l_食い違い群);
            }

            var l_窓長 = Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_確かめる窓の分子 / C_確かめる窓の分母);
            int l_区間数 = 0, l_片が短い数 = 0, l_骨格がリードにある数 = 0, l_決まらない数 = 0, l_錨を確かめられない数 = 0;
            List<(string A_ID, string A_配列)> l_結果 = [];
            var l_差し替えた数 = 0;
            for (var l_番号 = 0; l_番号 < p_配列群.Count; l_番号++)
            {
                var (l_ID, l_配列) = p_配列群[l_番号];
                var l_区間群 = Get_区間群(l_食い違い群.Where(x => x.A_配列番号 == l_番号));
                var l_出力 = new StringBuilder(l_配列.Length);
                var l_位置 = 0;
                foreach (var (l_始め, l_終わり, l_形群) in l_区間群)
                {
                    l_区間数++;
                    var l_差し替える形 = Get_差し替える形(l_配列, l_始め, l_終わり, l_形群, p_リード索引, p_リード長, l_窓長, ref l_片が短い数, ref l_骨格がリードにある数, ref l_決まらない数, ref l_錨を確かめられない数);
                    if (l_差し替える形 is null)
                    {
                        continue;
                    }

                    l_差し替えた数++;
                    Logger.V_出力(メッセージID.ほかのkで差し替えた所の内訳, l_ID, l_始め, l_終わり - l_始め, l_差し替える形.Length, l_配列[Math.Max(0, l_始め - C_ログの脇の長さ)..l_始め], l_配列[l_終わり..Math.Min(l_配列.Length, l_終わり + C_ログの脇の長さ)]);
                    _ = l_出力.Append(l_配列.AsSpan(l_位置, l_始め - l_位置)).Append(l_差し替える形);
                    l_位置 = l_終わり;
                }

                _ = l_出力.Append(l_配列.AsSpan(l_位置));
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            p_統計 = new ほかのkで差し替えた統計(l_区間数, l_差し替えた数, l_片が短い数, l_骨格がリードにある数, l_決まらない数, l_錨を確かめられない数);
            return l_結果;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列群の C_錨長 の k-mer の正準形ごとに、両鎖を通して 1 回だけ出るものの場所を表にする<br/>
        /// 場所は (配列番号 &lt;&lt; 32) | (位置 &lt;&lt; 1) | (順向きが正準形なら 1)。2 回以上出るものは C_一意でない
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <returns>正準形から場所への表</returns>
        internal static Dictionary<ulong, long> Get_一意な錨の表(IReadOnlyList<(string A_ID, string A_配列)> p_配列群)
        {
            var l_表 = new Dictionary<ulong, long>(p_配列群.Sum(static x => x.A_配列.Length));
            for (var l_番号 = 0; l_番号 < p_配列群.Count; l_番号++)
            {
                foreach (var (l_位置, l_正準, l_順向き) in Get_錨群(p_配列群[l_番号].A_配列))
                {
                    ref var l_場所 = ref CollectionsMarshal.GetValueRefOrAddDefault(l_表, l_正準, out var l_Isある);
                    l_場所 = l_Isある ? C_一意でない : ((long)l_番号 << 32) | ((long)l_位置 << 1) | (l_順向き ? 1L : 0L);
                }
            }

            return l_表;
        }

        /// <summary>
        /// 配列の A・C・G・T だけからなる C_錨長 の k-mer ごとに、位置、正準形、順向きが正準形かを列挙する
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>k-mer の位置、正準形、順向きが正準形か</returns>
        internal static IEnumerable<(int A_位置, ulong A_正準, bool A_順向き)> Get_錨群(string p_配列)
        {
            const ulong l_マスク = (1UL << (2 * C_錨長)) - 1UL;
            var l_順 = 0UL;
            var l_逆 = 0UL;
            var l_連続 = 0;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_値 = ReadEndJoiner.Get_塩基値(p_配列[i]);
                if (l_値 < 0)
                {
                    l_連続 = 0;
                    continue;
                }

                l_順 = ((l_順 << 2) | (uint)l_値) & l_マスク;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_錨長 - 1)));
                l_連続++;
                if (l_連続 >= C_錨長)
                {
                    yield return (i - C_錨長 + 1, Math.Min(l_順, l_逆), l_順 < l_逆);
                }
            }
        }

        /// <summary>
        /// ほかの k の配列 1 本について、両方の向きで錨の並びを作り、最終配列と錨の間の長さが食い違う所を集める<br/>
        /// 錨は、最終配列で 1 回だけ出て、この配列の中でも 1 回だけ出る k-mer。並びは最終配列の配列ごとに分けて、別の配列に当たる錨を飛ばす
        /// </summary>
        /// <param name="p_ほかの配列">ほかの k の contig か scaffold の配列</param>
        /// <param name="p_錨表">最終配列の一意な錨の表</param>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_食い違い群">食い違いを足す先 (配列番号、最終配列での始めと終わり、差し替える配列)</param>
        private static void V_集める_食い違い(string p_ほかの配列, Dictionary<ulong, long> p_錨表, IReadOnlyList<(string A_ID, string A_配列)> p_配列群, List<(int A_配列番号, int A_始め, int A_終わり, string A_差し替える配列)> p_食い違い群)
        {
            List<(int A_位置, ulong A_正準, bool A_順向き, long A_場所)> l_当たり = [];
            foreach (var (l_位置, l_正準, l_順向き) in Get_錨群(p_ほかの配列))
            {
                if (p_錨表.TryGetValue(l_正準, out var l_場所) && l_場所 != C_一意でない)
                {
                    l_当たり.Add((l_位置, l_正準, l_順向き, l_場所));
                }
            }

            var l_数 = l_当たり.CountBy(static x => x.A_正準).ToDictionary();
            var l_一意 = l_当たり.Where(x => l_数[x.A_正準] == 1).ToList();
            var l_同じ向き = l_一意.Where(static x => x.A_順向き == ((x.A_場所 & 1L) == 1L)).Select(static x => (A_q: x.A_位置, A_配列番号: (int)(x.A_場所 >> 32), A_t: (int)((x.A_場所 >> 1) & int.MaxValue))).ToList();
            foreach (var l_並び in l_同じ向き.GroupBy(static x => x.A_配列番号))
            {
                V_集める_並びの食い違い([.. l_並び], p_ほかの配列, p_配列群, p_食い違い群);
            }

            var l_逆の向き = l_一意.Where(static x => x.A_順向き != ((x.A_場所 & 1L) == 1L)).Select(x => (A_q: p_ほかの配列.Length - C_錨長 - x.A_位置, A_配列番号: (int)(x.A_場所 >> 32), A_t: (int)((x.A_場所 >> 1) & int.MaxValue))).Reverse().ToList();
            if (l_逆の向き.Count < 2)
            {
                return;
            }

            var l_逆相補 = Util.V_逆相補_曖昧塩基あり(p_ほかの配列);
            foreach (var l_並び in l_逆の向き.GroupBy(static x => x.A_配列番号))
            {
                V_集める_並びの食い違い([.. l_並び], l_逆相補, p_配列群, p_食い違い群);
            }
        }

        /// <summary>
        /// ほかの k の配列 (最終配列と同じ向き) の上で前から並べた錨で、隣り合う 2 つが最終配列の同じ配列に同じ順で並ぶのに、間の長さが違う所を集める<br/>
        /// 食い違いどうしの間の錨の並びが C_まとめる間の上限 以下なら、1 つの区間にまとめる<br/>
        /// 区間の長さがどちらも C_区間の長さの上限 以下で、長さの差が C_長さの差の下限 以上、最終配列の側に N が無く、ほかの k の側が A・C・G・T だけのものを足す
        /// </summary>
        /// <param name="p_並び">錨の並び (ほかの k の配列での位置、最終配列の配列番号と位置)</param>
        /// <param name="p_向きの配列">最終配列と同じ向きにしたほかの k の配列</param>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_食い違い群">食い違いを足す先</param>
        private static void V_集める_並びの食い違い(List<(int A_q, int A_配列番号, int A_t)> p_並び, string p_向きの配列, IReadOnlyList<(string A_ID, string A_配列)> p_配列群, List<(int A_配列番号, int A_始め, int A_終わり, string A_差し替える配列)> p_食い違い群)
        {
            var l_組の始め = -1;
            var l_組の終わり = -1;
            for (var i = 0; i + 1 < p_並び.Count; i++)
            {
                var (l_前, l_後) = (p_並び[i], p_並び[i + 1]);
                if (l_前.A_配列番号 != l_後.A_配列番号 || l_後.A_t <= l_前.A_t || l_後.A_q - l_前.A_q > C_区間の長さの上限 || l_後.A_t - l_前.A_t > C_区間の長さの上限)
                {
                    V_足す_食い違い(p_並び, l_組の始め, l_組の終わり, p_向きの配列, p_配列群, p_食い違い群);
                    l_組の始め = -1;
                    continue;
                }

                if (l_後.A_q - l_前.A_q == l_後.A_t - l_前.A_t)
                {
                    continue;
                }

                if (l_組の始め >= 0 && l_前.A_q - p_並び[l_組の終わり + 1].A_q <= C_まとめる間の上限)
                {
                    l_組の終わり = i;
                    continue;
                }

                V_足す_食い違い(p_並び, l_組の始め, l_組の終わり, p_向きの配列, p_配列群, p_食い違い群);
                (l_組の始め, l_組の終わり) = (i, i);
            }

            V_足す_食い違い(p_並び, l_組の始め, l_組の終わり, p_向きの配列, p_配列群, p_食い違い群);
        }

        /// <summary>
        /// 錨の並びの l_組の始め 番目から l_組の終わり + 1 番目の錨までを 1 つの食い違いの区間にして、条件を満たせば足す<br/>
        /// 区間は、始めの錨の始まりから終わりの錨の始まりまで (終わりの錨は含まない)
        /// </summary>
        /// <param name="p_並び">錨の並び</param>
        /// <param name="p_組の始め">区間の始めの錨の番号 (負なら何もしない)</param>
        /// <param name="p_組の終わり">区間の終わりの錨の 1 つ前の番号</param>
        /// <param name="p_向きの配列">最終配列と同じ向きにしたほかの k の配列</param>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_食い違い群">食い違いを足す先</param>
        private static void V_足す_食い違い(List<(int A_q, int A_配列番号, int A_t)> p_並び, int p_組の始め, int p_組の終わり, string p_向きの配列, IReadOnlyList<(string A_ID, string A_配列)> p_配列群, List<(int A_配列番号, int A_始め, int A_終わり, string A_差し替える配列)> p_食い違い群)
        {
            if (p_組の始め < 0)
            {
                return;
            }

            var (l_前, l_後) = (p_並び[p_組の始め], p_並び[p_組の終わり + 1]);
            var (l_q長, l_t長) = (l_後.A_q - l_前.A_q, l_後.A_t - l_前.A_t);
            if (l_q長 > C_区間の長さの上限 || l_t長 > C_区間の長さの上限 || Math.Abs(l_q長 - l_t長) < C_長さの差の下限)
            {
                return;
            }

            var l_配列 = p_配列群[l_前.A_配列番号].A_配列;
            var l_差し替える配列 = p_向きの配列.Substring(l_前.A_q, l_q長);
            if (l_配列.AsSpan(l_前.A_t, l_t長).IndexOfAnyExcept(Consts.塩基文字) >= 0 || l_差し替える配列.AsSpan().IndexOfAnyExcept(Consts.塩基文字) >= 0)
            {
                return;
            }

            p_食い違い群.Add((l_前.A_配列番号, l_前.A_t, l_後.A_t, l_差し替える配列));
        }

        /// <summary>
        /// 1 本の最終配列の食い違いを、最終配列の上で重なるものどうし 1 つの区間にまとめる<br/>
        /// 区間ごとに、区間全体をほかの k の形にした配列 (区間の中の食い違いの外は最終配列のまま) を重ねずに集める
        /// </summary>
        /// <param name="p_食い違い群">1 本の最終配列の食い違い</param>
        /// <returns>区間の始め、終わり、ほかの k の形の群 (始めの順)</returns>
        internal static List<(int A_始め, int A_終わり, List<(int A_始め, int A_終わり, string A_差し替える配列)> A_形群)> Get_区間群(IEnumerable<(int A_配列番号, int A_始め, int A_終わり, string A_差し替える配列)> p_食い違い群)
        {
            List<(int A_始め, int A_終わり, List<(int A_始め, int A_終わり, string A_差し替える配列)> A_形群)> l_区間群 = [];
            foreach (var l_食い違い in p_食い違い群.OrderBy(static x => x.A_始め).ThenBy(static x => x.A_終わり))
            {
                if (l_区間群.Count > 0 && l_食い違い.A_始め < l_区間群[^1].A_終わり)
                {
                    var l_最後 = l_区間群[^1];
                    l_最後.A_形群.Add((l_食い違い.A_始め, l_食い違い.A_終わり, l_食い違い.A_差し替える配列));
                    l_区間群[^1] = (l_最後.A_始め, Math.Max(l_最後.A_終わり, l_食い違い.A_終わり), l_最後.A_形群);
                }
                else
                {
                    l_区間群.Add((l_食い違い.A_始め, l_食い違い.A_終わり, [(l_食い違い.A_始め, l_食い違い.A_終わり, l_食い違い.A_差し替える配列)]));
                }
            }

            return l_区間群;
        }

        /// <summary>
        /// 1 つの区間を差し替える形を求める<br/>
        /// 両脇の片 (隣の N の列か配列の端まで) がどちらも窓長 - 1 以上あり、最終配列の形が PairBridgeFiller.Is写しがリードに続けて出る を満たさず、<br/>
        /// ほかの k の形 (区間全体) のうち満たすものがただ 1 つで、Is錨の間がリードで確かめられる なら、それを返す
        /// </summary>
        /// <param name="p_配列">最終配列</param>
        /// <param name="p_始め">区間の始め</param>
        /// <param name="p_終わり">区間の終わり (含まない)</param>
        /// <param name="p_形群">区間の中の食い違い (始め、終わり、差し替える配列)</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_窓長">リードで確かめる窓の長さ</param>
        /// <param name="p_片が短い数">両脇の片が窓長 - 1 より短い数 (加算する)</param>
        /// <param name="p_骨格がリードにある数">最終配列の形がリードに続けて出てくる数 (加算する)</param>
        /// <param name="p_決まらない数">リードに続けて出てくるほかの k の形がただ 1 つに決まらない数 (加算する)</param>
        /// <param name="p_錨を確かめられない数">ただ 1 つに決まったが、両端の錨の間の長さをリードで確かめられない数 (加算する)</param>
        /// <returns>区間を差し替える配列、差し替えなければ null</returns>
        private static string? Get_差し替える形(string p_配列, int p_始め, int p_終わり, List<(int A_始め, int A_終わり, string A_差し替える配列)> p_形群, ReadMinimizerIndex p_リード索引, int p_リード長, int p_窓長, ref int p_片が短い数, ref int p_骨格がリードにある数, ref int p_決まらない数, ref int p_錨を確かめられない数)
        {
            var l_片の長さ = p_窓長 - 1;
            var l_左の片の始め = p_配列.LastIndexOf(Consts.塩基コード.N, Math.Max(0, p_始め - 1)) + 1;
            var l_右の片の終わり = p_配列.IndexOf(Consts.塩基コード.N, p_終わり);
            l_右の片の終わり = l_右の片の終わり < 0 ? p_配列.Length : l_右の片の終わり;
            if (p_始め - l_左の片の始め < l_片の長さ || l_右の片の終わり - p_終わり < l_片の長さ)
            {
                p_片が短い数++;
                return null;
            }

            var l_左 = p_配列.Substring(p_始め - l_片の長さ, l_片の長さ);
            var l_右 = p_配列.Substring(p_終わり, l_片の長さ);
            if (PairBridgeFiller.Is写しがリードに続けて出る(p_リード索引, l_左, p_配列[p_始め..p_終わり], l_右, p_窓長))
            {
                p_骨格がリードにある数++;
                return null;
            }

            var l_形群 = p_形群.Select(x => (A_食い違い: x, A_形: string.Concat(p_配列.AsSpan(p_始め, x.A_始め - p_始め), x.A_差し替える配列, p_配列.AsSpan(x.A_終わり, p_終わり - x.A_終わり)))).DistinctBy(static x => x.A_形, StringComparer.Ordinal);
            var l_通った = l_形群.Where(x => PairBridgeFiller.Is写しがリードに続けて出る(p_リード索引, l_左, x.A_形, l_右, p_窓長)).Take(2).ToList();
            if (l_通った.Count != 1)
            {
                p_決まらない数++;
                return null;
            }

            var (l_食い違い, l_形) = l_通った[0];
            if (!Is錨の間がリードで確かめられる(l_左 + l_形 + l_右, l_左.Length + l_食い違い.A_始め - p_始め, Math.Min(l_食い違い.A_終わり - l_食い違い.A_始め, l_食い違い.A_差し替える配列.Length) + C_錨長, p_リード索引, p_リード長, p_窓長))
            {
                p_錨を確かめられない数++;
                return null;
            }

            return l_形;
        }

        /// <summary>
        /// ほかの k の形で、食い違いの左の錨から、最終配列とほかの k の区間の短いほうに C_錨長 を足した長さ (確かめる長さ) の配列が、リードに C_支持数の下限 か所以上あるか<br/>
        /// この配列は、短いほうの形なら両端の錨を、長いほうの形なら左の錨と、短いほうの形で右の錨があるはずの所を含むので、間の長さを取り違えた形はリードに無い<br/>
        /// 確かめる長さが窓長以下なら、窓の確かめで済んでいるので真。リード長より長ければ確かめられないので偽
        /// </summary>
        /// <param name="p_繋いだ配列">左の片、ほかの k の形 (区間全体)、右の片を繋いだ配列</param>
        /// <param name="p_錨の位置">繋いだ配列での、食い違いの左の錨の始まり</param>
        /// <param name="p_確かめる長さ">確かめる配列の長さ</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_窓長">リードで確かめる窓の長さ</param>
        /// <returns>確かめられれば真</returns>
        internal static bool Is錨の間がリードで確かめられる(string p_繋いだ配列, int p_錨の位置, int p_確かめる長さ, ReadMinimizerIndex p_リード索引, int p_リード長, int p_窓長)
        {
            return p_確かめる長さ <= p_窓長
                || (p_確かめる長さ <= p_リード長 && p_リード索引.Get_出現数(p_繋いだ配列.AsSpan(p_錨の位置, p_確かめる長さ), C_支持数の下限) >= C_支持数の下限);
        }

        #endregion
    }
}
