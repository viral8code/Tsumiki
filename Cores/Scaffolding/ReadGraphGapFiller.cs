using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 最終配列に残った N の列を、リードに 2 本以上ある k-mer だけを辿る 1 本の道で埋める<br/>
    /// k はリード長の 70%、50%、35% の順に、道が見つからないときだけ下げる。道が 2 本以上なら、短い順に並べた道をリード長の窓で比べて 1 本を選ぶ。<br/>
    /// 1 本の道は、埋める形なら左の錨から右の錨までを、畳む形 (間が負) なら畳んだ形の深い錨から深い錨までを、リードで確かめられるときだけ埋める
    /// </summary>
    internal static class ReadGraphGapFiller
    {
        #region 定数

        /// <summary>
        /// 道を辿る k-mer の長さを試す順の、リード長に対する百分率<br/>
        /// 最初は PairBridgeFiller.Is写しがリードに続けて出る の窓と同じ 70%。端の k-mer がリードに無いか道が無いときだけ、次へ下げる
        /// </summary>
        private static readonly int[] C_k長の百分率群 = [70, 50, 35];

        /// <summary>
        /// 百分率の分母
        /// </summary>
        private const int C_百分率の分母 = 100;

        /// <summary>
        /// 道に使う k-mer に要る、リードの出現数の下限
        /// </summary>
        private const int C_支持数の下限 = 2;

        /// <summary>
        /// 端の k-mer がリードに無いとき、片の端から内へずらす長さの上限<br/>
        /// 最終の配列の片の端には、継ぎ目の誤りが残ることがある
        /// </summary>
        private const int C_端から削る上限 = 200;

        /// <summary>
        /// 端の k-mer を内へずらす幅
        /// </summary>
        private const int C_端から削る幅 = 10;

        /// <summary>
        /// 埋める間の長さの上限
        /// </summary>
        private const int C_間の長さの上限 = 2000;

        /// <summary>
        /// 道を 1 塩基進める回数の上限<br/>
        /// 縦の反復の周回などで道が増え続けるときに打ち切る
        /// </summary>
        private const int C_進める回数の上限 = 20_000;

        /// <summary>
        /// 埋めた所のログに出す、両脇の配列の長さ<br/>
        /// 評価のあとで、最終配列の中の場所を探すのに使う
        /// </summary>
        private const int C_ログの脇の長さ = 60;

        /// <summary>
        /// 道が 2 本以上のとき、短い順に並べる道の数の上限
        /// </summary>
        private const int C_並べる道の上限 = 12;

        /// <summary>
        /// 道を短い順に並べるとき、道を 1 塩基進める回数の上限
        /// </summary>
        private const int C_並べる回数の上限 = 100_000;

        /// <summary>
        /// 道どうしを比べる窓の長さの、リード長に対する百分率<br/>
        /// リード長と同じだと窓を丸ごと含むリードが少なく、窓にリードが無いことが当てにならない
        /// </summary>
        private const int C_窓の百分率 = 85;

        /// <summary>
        /// ほかの道にだけある窓にリードが無いことで、その道を外すのに要る窓の数の下限
        /// </summary>
        private const int C_窓の下限数 = 10;

        /// <summary>
        /// 窓がほかの道の形にあるかを調べるとき、道の両側に付ける片の長さ
        /// </summary>
        private const int C_広い形の脇の長さ = 2000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群の N の列を、リードに 2 本以上ある k-mer の 1 本の道で埋める<br/>
        /// N の列ごとに、左の片 (手前の N の列か配列の端から。手前の N の列を埋めたときは、そこで削った所から) と右の片 (次の N の列か配列の端まで) を Get_kを下げて埋める形 に渡し、<br/>
        /// 形が決まれば、左の片を道の始まりの前まで、右の片を道の終わりのあとから残して、道で繋ぐ。埋めなかった N の列は理由を、窓で比べた N の列はその結果を、ログに出す
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_統計">埋めた結果の統計</param>
        /// <param name="p_埋めた番号">配列 ID ごとの、埋めた N の列の番号 (埋める前の配列での 0 始まり)</param>
        /// <returns>埋めた後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_埋めた配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, ReadMinimizerIndex p_リード索引, int p_リード長, out リードの道で埋めた統計 p_統計, out Dictionary<string, HashSet<int>> p_埋めた番号)
        {
            int[] l_k長群 = [.. C_k長の百分率群.Select(x => Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * x / C_百分率の分母)).Distinct()];
            var l_錨表 = OtherKRepeatReplacer.Get_一意な錨の表(p_配列群);
            p_埋めた番号 = [];
            var l_理由ごとの数 = new int[Enum.GetValues<リードの道で埋めない理由>().Length];
            int l_対象数 = 0, l_埋めた数 = 0, l_畳んだ数 = 0, l_小さいkで埋めた数 = 0, l_窓で選んだ数 = 0;
            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in p_配列群)
            {
                var l_列群 = OverlapGapFolder.Get_Nの列群(l_配列);
                var l_出力 = new StringBuilder(l_配列.Length);
                var l_次の始め = 0;
                for (var l_番号 = 0; l_番号 < l_列群.Count; l_番号++)
                {
                    var (l_始め, l_終わり) = l_列群[l_番号];
                    var l_左 = l_配列[l_次の始め..l_始め];
                    var l_右 = l_配列[l_終わり..(l_番号 + 1 < l_列群.Count ? l_列群[l_番号 + 1].A_始め : l_配列.Length)];
                    l_対象数++;
                    var (l_得た形, l_k長, l_理由, l_窓) = Get_kを下げて埋める形(l_左, l_右, l_錨表, p_リード索引, p_リード長, l_k長群);
                    if (l_窓 is { } l_記録)
                    {
                        if (l_記録.A_選んだ順 > 0)
                        {
                            l_窓で選んだ数++;
                            Logger.V_出力(メッセージID.リードの道を窓で選んだ, l_ID, l_始め, l_終わり - l_始め, l_記録.A_道の数, l_記録.A_窓長, l_記録.A_選んだ順, l_記録.A_相手の窓の最少);
                        }
                        else
                        {
                            Logger.V_出力(メッセージID.リードの道を窓で選ばなかった, l_ID, l_始め, l_終わり - l_始め, l_記録.A_道の数, l_記録.A_Is打ち切り, l_記録.A_畳む形の数, l_記録.A_窓長, l_記録.A_条件を満たす数);
                        }
                    }

                    if (l_得た形 is not { } l_形)
                    {
                        l_理由ごとの数[(int)l_理由]++;
                        Logger.V_出力(メッセージID.リードの道で埋めなかったNの列の内訳, l_ID, l_始め, l_終わり - l_始め, l_理由, l_k長);
                        _ = l_出力.Append(l_配列.AsSpan(l_次の始め, l_終わり - l_次の始め));
                        l_次の始め = l_終わり;
                        continue;
                    }

                    var l_残す左 = l_左[..^l_形.A_左の削り];
                    var l_間 = l_形.A_道.Length - (2 * l_k長);
                    l_埋めた数++;
                    l_畳んだ数 += l_間 - l_形.A_左の削り - l_形.A_右の削り < 0 ? 1 : 0;
                    l_小さいkで埋めた数 += l_k長 != l_k長群[0] ? 1 : 0;
                    Logger.V_出力(メッセージID.リードの道で埋めたNの列の内訳, l_ID, l_始め, l_終わり - l_始め, l_間, l_形.A_左の削り, l_形.A_右の削り, l_残す左[^Math.Min(l_残す左.Length, C_ログの脇の長さ)..], l_右.Substring(l_形.A_右の削り, Math.Min(l_右.Length - l_形.A_右の削り, C_ログの脇の長さ)), l_k長);
                    _ = l_出力.Append(l_残す左.AsSpan(0, l_残す左.Length - l_k長)).Append(l_形.A_道);
                    if (!p_埋めた番号.TryGetValue(l_ID, out var l_番号群))
                    {
                        l_番号群 = [];
                        p_埋めた番号[l_ID] = l_番号群;
                    }

                    _ = l_番号群.Add(l_番号);
                    l_次の始め = l_終わり + l_形.A_右の削り + l_k長;
                }

                _ = l_出力.Append(l_配列.AsSpan(l_次の始め));
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            p_統計 = new リードの道で埋めた統計(l_k長群, l_対象数, l_埋めた数, l_畳んだ数, l_小さいkで埋めた数, l_窓で選んだ数, l_理由ごとの数[(int)リードの道で埋めない理由.片が短い], l_理由ごとの数[(int)リードの道で埋めない理由.端がリードに無い], l_理由ごとの数[(int)リードの道で埋めない理由.道が無い], l_理由ごとの数[(int)リードの道で埋めない理由.道が決まらない], l_理由ごとの数[(int)リードの道で埋めない理由.畳む形を確かめられない], l_理由ごとの数[(int)リードの道で埋めない理由.錨から錨までを確かめられない]);
            return l_結果;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 つの N の列を、k 長群の順に Get_埋める形 で試す<br/>
        /// 端の k-mer がリードに無いか道が無ければ次の k を試し、ほかの理由で埋めないときと、最後の k でも埋まらないときは止める
        /// </summary>
        /// <param name="p_左">左の片</param>
        /// <param name="p_右">右の片</param>
        /// <param name="p_錨表">アセンブリの一意な錨の表 (OtherKRepeatReplacer.Get_一意な錨の表)</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_k長群">試す k-mer の長さ (試す順)</param>
        /// <returns>埋める形 (埋めなければ null)、使ったか最後に試した k、埋めなかった理由 (埋めるときは既定値)、その k で道を窓で比べた結果 (比べなければ null)</returns>
        private static ((int A_左の削り, string A_道, int A_右の削り)? A_形, int A_k長, リードの道で埋めない理由 A_理由, 窓で道を選んだ記録? A_窓) Get_kを下げて埋める形(string p_左, string p_右, Dictionary<ulong, long> p_錨表, ReadMinimizerIndex p_索引, int p_リード長, int[] p_k長群)
        {
            var l_理由 = default(リードの道で埋めない理由);
            foreach (var l_k長 in p_k長群)
            {
                var l_形 = Get_埋める形(p_左, p_右, p_錨表, p_索引, p_リード長, l_k長, out l_理由, out var l_窓);
                if (l_形 is not null || l_理由 is not (リードの道で埋めない理由.端がリードに無い or リードの道で埋めない理由.道が無い))
                {
                    return (l_形, l_k長, l_理由, l_窓);
                }
            }

            return (null, p_k長群[^1], l_理由, null);
        }

        /// <summary>
        /// 1 つの N の列を、1 つの k で埋める形を求める<br/>
        /// 両脇の片がどちらも PairBridgeFiller.C_反復を調べる長さ (と k 長 + C_端から削る上限) 以上あるとき、片の端の k-mer (リードに無ければ内へずらす) を始まりと終わりにして道を探す。<br/>
        /// 道が 1 本なら、埋める形は左の錨から右の錨までが OtherKDeepGapBridger.Is錨から錨までがリードにある を、畳む形は Is畳んだ形がリードにある を満たすときだけ返す。<br/>
        /// 道が 2 本以上か打ち切りなら、Get_短い順の道群 で並べ直し、Get_窓で選んだ道 で選べた道を返す
        /// </summary>
        /// <param name="p_左">左の片</param>
        /// <param name="p_右">右の片</param>
        /// <param name="p_錨表">アセンブリの一意な錨の表 (OtherKRepeatReplacer.Get_一意な錨の表)</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_k長">道を辿る k-mer の長さ</param>
        /// <param name="p_理由">埋めなかった理由 (形を返すときは既定値)</param>
        /// <param name="p_窓">道を窓で比べた結果 (比べなければ null)</param>
        /// <returns>左の片の末尾から削る長さ、道 (始まりの k-mer から終わりの k-mer まで)、右の片の先頭から削る長さ。埋めなければ null</returns>
        internal static (int A_左の削り, string A_道, int A_右の削り)? Get_埋める形(string p_左, string p_右, Dictionary<ulong, long> p_錨表, ReadMinimizerIndex p_索引, int p_リード長, int p_k長, out リードの道で埋めない理由 p_理由, out 窓で道を選んだ記録? p_窓)
        {
            p_理由 = default;
            p_窓 = null;
            var l_下限 = Math.Max(PairBridgeFiller.C_反復を調べる長さ, p_k長 + C_端から削る上限);
            if (p_左.Length < l_下限 || p_右.Length < l_下限)
            {
                p_理由 = リードの道で埋めない理由.片が短い;
                return null;
            }

            if (Get_端の削り(p_左, true, p_索引, p_k長) is not { } l_左の削り || Get_端の削り(p_右, false, p_索引, p_k長) is not { } l_右の削り)
            {
                p_理由 = リードの道で埋めない理由.端がリードに無い;
                return null;
            }

            var l_始まり = p_左.Substring(p_左.Length - l_左の削り - p_k長, p_k長);
            var l_終わり = p_右.Substring(l_右の削り, p_k長);
            var l_長さの上限 = (2 * p_k長) + l_左の削り + l_右の削り + C_間の長さの上限;
            var l_道群 = Get_道群(l_始まり, l_終わり, l_長さの上限, p_索引, out var l_Is打ち切り);
            if (l_道群.Count == 0 && !l_Is打ち切り)
            {
                p_理由 = リードの道で埋めない理由.道が無い;
                return null;
            }

            if (l_道群.Count != 1 || l_Is打ち切り)
            {
                var l_並べた道群 = Get_短い順の道群(l_始まり, l_終わり, l_長さの上限, p_索引, out var l_Is並べて打ち切り);
                var (l_選んだ道, l_記録) = Get_窓で選んだ道(l_並べた道群, l_Is並べて打ち切り, p_左[..^l_左の削り], p_右[l_右の削り..], p_k長, l_左の削り + l_右の削り, p_索引, p_リード長);
                p_窓 = l_記録;
                if (l_選んだ道 is null)
                {
                    p_理由 = リードの道で埋めない理由.道が決まらない;
                    return null;
                }

                return (l_左の削り, l_選んだ道, l_右の削り);
            }

            var l_道 = l_道群[0];
            var l_錨長 = OtherKRepeatReplacer.C_錨長;
            if (l_道.Length - (2 * p_k長) - l_左の削り - l_右の削り < 0)
            {
                if (!Is畳んだ形がリードにある(p_左[..^l_左の削り], l_道, p_右[l_右の削り..], p_k長, p_錨表, p_索引, p_リード長))
                {
                    p_理由 = リードの道で埋めない理由.畳む形を確かめられない;
                    return null;
                }
            }
            else if (!OtherKDeepGapBridger.Is錨から錨までがリードにある(l_道.Substring(p_k長 - l_錨長, l_道.Length - (2 * p_k長) + (2 * l_錨長)), p_索引, p_リード長, p_k長))
            {
                p_理由 = リードの道で埋めない理由.錨から錨までを確かめられない;
                return null;
            }

            return (l_左の削り, l_道, l_右の削り);
        }

        /// <summary>
        /// 畳む形で、畳んだ形 (削った後の左の片 + 道[k..] + 削った後の右の片[k..]) の、左の片の深い錨から右の片の深い錨までがリードにあるか<br/>
        /// 錨はどちらも片の端に一番近いものを使い、右の錨は左の片と重なる所 (2k − 道の長さ) より後ろに限る。<br/>
        /// 錨から錨までがリード長以下で、リードに C_支持数の下限 か所以上あれば真。どちらかの錨が無ければ偽
        /// </summary>
        /// <param name="p_左">削った後の左の片 (末尾は道の始まりの k-mer)</param>
        /// <param name="p_道">道 (始まりの k-mer から終わりの k-mer まで)</param>
        /// <param name="p_右">削った後の右の片 (先頭は道の終わりの k-mer)</param>
        /// <param name="p_k長">道を辿った k-mer の長さ</param>
        /// <param name="p_錨表">アセンブリの一意な錨の表</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <returns>確かめられれば真</returns>
        private static bool Is畳んだ形がリードにある(string p_左, string p_道, string p_右, int p_k長, Dictionary<ulong, long> p_錨表, ReadMinimizerIndex p_索引, int p_リード長)
        {
            var l_重なり = (2 * p_k長) - p_道.Length;
            var l_錨長 = OtherKRepeatReplacer.C_錨長;
            var l_左の錨群 = OtherKEdgeBridger.Get_錨の位置群(p_左, p_錨表, true);
            var l_右の錨 = OtherKEdgeBridger.Get_錨の位置群(p_右, p_錨表, false).Where(x => x >= l_重なり).DefaultIfEmpty(-1).First();
            if (l_左の錨群.Count == 0 || l_右の錨 < 0)
            {
                return false;
            }

            var l_長さ = p_左.Length - l_左の錨群[0] - l_重なり + l_右の錨 + l_錨長;
            if (l_長さ > p_リード長)
            {
                return false;
            }

            var l_右から足す長さ = Math.Clamp(l_右の錨 + l_錨長 - p_k長, 0, p_右.Length - p_k長);
            var l_畳んだ形 = string.Concat(p_左.AsSpan(l_左の錨群[0]), p_道.AsSpan(p_k長), p_右.AsSpan(p_k長, l_右から足す長さ));
            return p_索引.Get_出現数(l_畳んだ形.AsSpan(0, l_長さ), C_支持数の下限) >= C_支持数の下限;
        }

        /// <summary>
        /// 片の端の k-mer がリードに C_支持数の下限 本以上あるまで、端から C_端から削る幅 ずつ内へずらした長さを求める
        /// </summary>
        /// <param name="p_片">片</param>
        /// <param name="p_Is末尾">末尾の k-mer なら真、先頭の k-mer なら偽</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_k長">k-mer の長さ</param>
        /// <returns>端から削る長さ、C_端から削る上限 までずらしても無ければ null</returns>
        private static int? Get_端の削り(string p_片, bool p_Is末尾, ReadMinimizerIndex p_索引, int p_k長)
        {
            for (var l_削り = 0; l_削り <= C_端から削る上限; l_削り += C_端から削る幅)
            {
                var l_端のkmer = p_片.AsSpan(p_Is末尾 ? p_片.Length - l_削り - p_k長 : l_削り, p_k長);
                if (l_端のkmer.IndexOfAnyExcept(Consts.塩基文字) < 0 && p_索引.Get_出現数(l_端のkmer, C_支持数の下限) >= C_支持数の下限)
                {
                    return l_削り;
                }
            }

            return null;
        }

        /// <summary>
        /// 始まりの k-mer から、次の k-mer がリードに C_支持数の下限 本以上ある塩基へ 1 つずつ進み (深さ優先)、終わりの k-mer で終わる道を集める<br/>
        /// 道の長さが上限に達した枝は止め、道が 2 本見つかるか、進めた回数が C_進める回数の上限 を超えたら打ち切る
        /// </summary>
        /// <param name="p_始まり">始まりの k-mer</param>
        /// <param name="p_終わり">終わりの k-mer</param>
        /// <param name="p_長さの上限">道の長さの上限</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_Is打ち切り">進めた回数が上限を超えたら真</param>
        /// <returns>道 (始まりの k-mer から終わりの k-mer までの配列)。2 本見つかったら 2 本</returns>
        internal static List<string> Get_道群(string p_始まり, string p_終わり, int p_長さの上限, ReadMinimizerIndex p_索引, out bool p_Is打ち切り)
        {
            var l_k長 = p_始まり.Length;
            var l_道 = new StringBuilder(p_始まり, p_長さの上限);
            var l_窓 = new char[l_k長];
            Stack<(int A_長さ, char A_塩基)> l_枝 = new();
            List<string> l_道群 = [];
            V_積む_次の塩基(l_道, l_窓, p_索引, l_枝);
            var l_回数 = 0;
            p_Is打ち切り = false;
            while (l_枝.Count > 0)
            {
                if (++l_回数 > C_進める回数の上限)
                {
                    p_Is打ち切り = true;
                    break;
                }

                var (l_長さ, l_塩基) = l_枝.Pop();
                l_道.Length = l_長さ;
                _ = l_道.Append(l_塩基);
                if (Is末尾が一致(l_道, p_終わり))
                {
                    l_道群.Add(l_道.ToString());
                    if (l_道群.Count > 1)
                    {
                        break;
                    }

                    continue;
                }

                if (l_道.Length < p_長さの上限)
                {
                    V_積む_次の塩基(l_道, l_窓, p_索引, l_枝);
                }
            }

            return l_道群;
        }

        /// <summary>
        /// 始まりの k-mer から、次の k-mer がリードに C_支持数の下限 本以上ある塩基へ 1 つずつ進み (幅優先)、終わりの k-mer で終わる道を短い順に集める<br/>
        /// 道の長さが上限に達した枝は止め、道が C_並べる道の上限 本見つかるか、進めた回数が C_並べる回数の上限 を超えたら打ち切る。<br/>
        /// 縦の反復では周回のたびに同じ k-mer を問い合わせるので、k-mer がリードにあるかは表に覚える
        /// </summary>
        /// <param name="p_始まり">始まりの k-mer</param>
        /// <param name="p_終わり">終わりの k-mer</param>
        /// <param name="p_長さの上限">道の長さの上限</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_Is打ち切り">進めた回数が上限を超えたら真</param>
        /// <returns>道 (始まりの k-mer から終わりの k-mer までの配列) を短い順に</returns>
        internal static List<string> Get_短い順の道群(string p_始まり, string p_終わり, int p_長さの上限, ReadMinimizerIndex p_索引, out bool p_Is打ち切り)
        {
            List<(int A_親, char A_塩基, int A_長さ, string A_kmer)> l_節群 = [(-1, default, p_始まり.Length, p_始まり)];
            Queue<int> l_待ち = new();
            l_待ち.Enqueue(0);
            List<string> l_道群 = [];
            Dictionary<string, bool> l_ある表 = [];
            var l_回数 = 0;
            p_Is打ち切り = false;
            while (l_待ち.Count > 0 && l_道群.Count < C_並べる道の上限)
            {
                if (++l_回数 > C_並べる回数の上限)
                {
                    p_Is打ち切り = true;
                    break;
                }

                var l_番号 = l_待ち.Dequeue();
                var (_, _, l_長さ, l_kmer) = l_節群[l_番号];
                if (l_長さ >= p_長さの上限)
                {
                    continue;
                }

                foreach (var l_塩基 in Consts.塩基文字)
                {
                    var l_次 = l_kmer[1..] + l_塩基;
                    if (!l_ある表.TryGetValue(l_次, out var l_Isある))
                    {
                        l_Isある = p_索引.Get_出現数(l_次, C_支持数の下限) >= C_支持数の下限;
                        l_ある表[l_次] = l_Isある;
                    }

                    if (!l_Isある)
                    {
                        continue;
                    }

                    l_節群.Add((l_番号, l_塩基, l_長さ + 1, l_次));
                    if (l_次 != p_終わり)
                    {
                        l_待ち.Enqueue(l_節群.Count - 1);
                        continue;
                    }

                    l_道群.Add(Get_道の配列(l_節群, l_節群.Count - 1));
                    if (l_道群.Count >= C_並べる道の上限)
                    {
                        break;
                    }
                }
            }

            return l_道群;
        }

        /// <summary>
        /// 幅優先の節を親へ辿り、始まりの k-mer からその節までの道の配列を組む
        /// </summary>
        /// <param name="p_節群">幅優先の節 (親の番号、足した塩基、道の長さ、末尾の k-mer。0 番は始まりの k-mer)</param>
        /// <param name="p_番号">道の終わりの節の番号</param>
        /// <returns>道の配列</returns>
        private static string Get_道の配列(List<(int A_親, char A_塩基, int A_長さ, string A_kmer)> p_節群, int p_番号)
        {
            var l_塩基群 = new char[p_節群[p_番号].A_長さ];
            for (var l_番号 = p_番号; p_節群[l_番号].A_親 >= 0; l_番号 = p_節群[l_番号].A_親)
            {
                l_塩基群[p_節群[l_番号].A_長さ - 1] = p_節群[l_番号].A_塩基;
            }

            p_節群[0].A_kmer.CopyTo(l_塩基群);
            return new string(l_塩基群);
        }

        /// <summary>
        /// 短い順に並べた道を、長さ リード長 × C_窓の百分率 / 100 の窓で比べ、条件を満たす道が 1 本だけならそれを選ぶ。並べた道に畳む形があれば選ばない<br/>
        /// 条件は、その道の形の窓がどれもリードに C_支持数の下限 か所以上あり、ほかのどの道についても、ほかの道の形にだけある窓がリードに 1 か所も無く、<br/>
        /// その窓が C_窓の下限数 個以上あるか、その道の形にだけある窓があること。形は、片の端 リード長 bp + 道の間 + 片の端 リード長 bp
        /// </summary>
        /// <param name="p_道群">短い順に並べた道 (始まりの k-mer から終わりの k-mer まで)</param>
        /// <param name="p_Is打ち切り">並べる途中で打ち切ったか</param>
        /// <param name="p_左">削った後の左の片 (末尾は道の始まりの k-mer)</param>
        /// <param name="p_右">削った後の右の片 (先頭は道の終わりの k-mer)</param>
        /// <param name="p_k長">道を辿った k-mer の長さ</param>
        /// <param name="p_削りの和">両側で削った長さの和</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <returns>選んだ道 (選ばなければ null) と、比べた結果</returns>
        private static (string? A_道, 窓で道を選んだ記録 A_記録) Get_窓で選んだ道(List<string> p_道群, bool p_Is打ち切り, string p_左, string p_右, int p_k長, int p_削りの和, ReadMinimizerIndex p_索引, int p_リード長)
        {
            var l_窓長 = Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_窓の百分率 / C_百分率の分母);
            var l_畳む形の数 = p_道群.Count(x => x.Length - (2 * p_k長) - p_削りの和 < 0);
            if (p_道群.Count < 2 || l_畳む形の数 > 0)
            {
                return (null, new 窓で道を選んだ記録(p_道群.Count, p_Is打ち切り, l_畳む形の数, 0, 0, l_窓長, 0));
            }

            var l_左の脇 = p_左[^Math.Min(p_左.Length, p_リード長)..];
            var l_右の脇 = p_右[..Math.Min(p_右.Length, p_リード長)];
            var l_左の広い脇 = p_左[^Math.Min(p_左.Length, C_広い形の脇の長さ)..];
            var l_右の広い脇 = p_右[..Math.Min(p_右.Length, C_広い形の脇の長さ)];
            string[] l_形群 = [.. p_道群.Select(x => string.Concat(l_左の脇, x.AsSpan(p_k長, x.Length - (2 * p_k長)), l_右の脇))];
            string[] l_広い形群 = [.. p_道群.Select(x => string.Concat(l_左の広い脇, x.AsSpan(p_k長, x.Length - (2 * p_k長)), l_右の広い脇))];
            Dictionary<string, int> l_数の表 = [];
            List<(int A_番号, int A_最少)> l_満たす道群 = [];
            for (var i = 0; i < l_形群.Length; i++)
            {
                if (Get_窓群(l_形群[i], l_窓長).Any(x => Get_窓の出現数(x, l_数の表, p_索引) < C_支持数の下限))
                {
                    continue;
                }

                var l_Is満たす = true;
                var l_最少 = int.MaxValue;
                for (var j = 0; j < l_形群.Length && l_Is満たす; j++)
                {
                    if (j == i)
                    {
                        continue;
                    }

                    string[] l_相手だけの窓群 = [.. Get_窓群(l_形群[j], l_窓長).Where(x => !l_広い形群[i].Contains(x, StringComparison.Ordinal))];
                    l_Is満たす = !l_相手だけの窓群.Any(x => Get_窓の出現数(x, l_数の表, p_索引) > 0)
                        && (l_相手だけの窓群.Length >= C_窓の下限数 || Get_窓群(l_形群[i], l_窓長).Any(x => !l_広い形群[j].Contains(x, StringComparison.Ordinal)));
                    l_最少 = Math.Min(l_最少, l_相手だけの窓群.Length);
                }

                if (l_Is満たす)
                {
                    l_満たす道群.Add((i, l_最少));
                }
            }

            return l_満たす道群.Count == 1
                ? (p_道群[l_満たす道群[0].A_番号], new 窓で道を選んだ記録(p_道群.Count, p_Is打ち切り, 0, 1, l_満たす道群[0].A_番号 + 1, l_窓長, l_満たす道群[0].A_最少))
                : (null, new 窓で道を選んだ記録(p_道群.Count, p_Is打ち切り, 0, l_満たす道群.Count, 0, l_窓長, 0));
        }

        /// <summary>
        /// 配列の長さ p_窓長 の窓を、先頭から 1 塩基ずつずらして並べる
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_窓長">窓の長さ</param>
        /// <returns>窓の配列</returns>
        private static IEnumerable<string> Get_窓群(string p_配列, int p_窓長)
        {
            for (var l_始め = 0; l_始め + p_窓長 <= p_配列.Length; l_始め++)
            {
                yield return p_配列.Substring(l_始め, p_窓長);
            }
        }

        /// <summary>
        /// 窓がリードかその逆相補に出てくる箇所の数を、C_支持数の下限 まで数える (表に覚えて、同じ窓は 1 回だけ問い合わせる)。A・C・G・T 以外を含む窓は 0
        /// </summary>
        /// <param name="p_窓">窓</param>
        /// <param name="p_数の表">数えた窓の表</param>
        /// <param name="p_索引">リードの索引</param>
        /// <returns>出てくる箇所の数 (C_支持数の下限 で頭打ち)</returns>
        private static int Get_窓の出現数(string p_窓, Dictionary<string, int> p_数の表, ReadMinimizerIndex p_索引)
        {
            if (!p_数の表.TryGetValue(p_窓, out var l_数))
            {
                l_数 = p_窓.AsSpan().IndexOfAnyExcept(Consts.塩基文字) < 0 ? p_索引.Get_出現数(p_窓, C_支持数の下限) : 0;
                p_数の表[p_窓] = l_数;
            }

            return l_数;
        }

        /// <summary>
        /// 道の末尾 k - 1 塩基のあとに各塩基を付けた k-mer のうち、リードに C_支持数の下限 本以上あるものの塩基を、枝として積む
        /// </summary>
        /// <param name="p_道">今の道</param>
        /// <param name="p_窓">k-mer を組む作業域 (長さ k)</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_枝">枝を積む先 (積む前の道の長さと塩基)</param>
        private static void V_積む_次の塩基(StringBuilder p_道, char[] p_窓, ReadMinimizerIndex p_索引, Stack<(int A_長さ, char A_塩基)> p_枝)
        {
            var l_k長 = p_窓.Length;
            p_道.CopyTo(p_道.Length - l_k長 + 1, p_窓.AsSpan(0, l_k長 - 1), l_k長 - 1);
            foreach (var l_塩基 in Consts.塩基文字)
            {
                p_窓[^1] = l_塩基;
                if (p_索引.Get_出現数(p_窓, C_支持数の下限) >= C_支持数の下限)
                {
                    p_枝.Push((p_道.Length, l_塩基));
                }
            }
        }

        /// <summary>
        /// 道が k より長く、末尾が終わりの k-mer と一致するか
        /// </summary>
        /// <param name="p_道">今の道</param>
        /// <param name="p_終わり">終わりの k-mer</param>
        /// <returns>一致すれば真</returns>
        private static bool Is末尾が一致(StringBuilder p_道, string p_終わり)
        {
            if (p_道.Length <= p_終わり.Length)
            {
                return false;
            }

            var l_始め = p_道.Length - p_終わり.Length;
            for (var i = 0; i < p_終わり.Length; i++)
            {
                if (p_道[l_始め + i] != p_終わり[i])
                {
                    return false;
                }
            }

            return true;
        }

        #endregion
    }
}
