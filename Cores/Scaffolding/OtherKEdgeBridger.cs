using Tsumiki.Commons;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// ペアの橋渡しで反復の端として除いた辺を、ほかの k の contig の配列で繋ぐ形を決める<br/>
    /// u の末尾と x の先頭の近くで、アセンブリの中で 1 回だけ出る錨を取り、ほかの k の contig が両方の錨を順に含めば、錨の間の配列で繋ぐ
    /// </summary>
    internal static class OtherKEdgeBridger
    {
        #region 定数

        /// <summary>
        /// 錨を探す、u の末尾と x の先頭の長さ<br/>
        /// 反復の端では端の近くの k-mer が一意でないので、端から離れた所まで探す
        /// </summary>
        private const int C_錨を探す長さ = 2000;

        /// <summary>
        /// 試す錨どうしの間隔の下限<br/>
        /// 端の近くの継ぎ目の k-mer がほかの k に無いときに、少し離れた錨を試すため
        /// </summary>
        private const int C_錨の間隔 = 25;

        /// <summary>
        /// 片側で試す錨の数の上限
        /// </summary>
        private const int C_試す錨の数 = 20;

        /// <summary>
        /// 繋ぐ配列をリードで確かめる窓の長さの、リード長に対する比の分子<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分子 = 7;

        /// <summary>
        /// 繋ぐ配列をリードで確かめる窓の長さの、リード長に対する比の分母<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分母 = 10;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// u の末尾と x の先頭を、ほかの k の contig の配列で繋ぐ形を求める<br/>
        /// u の末尾から近い順と x の先頭から近い順に、アセンブリで 1 回だけ出る錨 (OtherKRepeatReplacer の錨の表) を C_錨の間隔 おきに取り、ほかの k にどちらの向きでも無い錨を除いて、<br/>
        /// ほかの k の連結配列で左の錨のあとに右の錨が続く最初の組で、錨の間の配列を集める。<br/>
        /// そのうち、u の左の錨までと x の右の錨からの片で挟んで PairBridgeFiller.Is写しがリードに続けて出る を満たすものが、ただ 1 つなら返す
        /// </summary>
        /// <param name="p_u文字列">u の文字列 (この向きの末尾を繋ぐ)</param>
        /// <param name="p_x文字列">x の文字列 (この向きの先頭を繋ぐ)</param>
        /// <param name="p_錨表">アセンブリの一意な錨の表 (OtherKRepeatReplacer.Get_一意な錨の表)</param>
        /// <param name="p_ほかのk">ほかの k の contig を区切りで連結した配列</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <returns>u の末尾から削る長さ、間に入れる配列、x の先頭から削る長さ、左右の錨の位置。繋がなければ null</returns>
        public static (int A_u側, string A_埋める配列, int A_x側, int A_左の錨, int A_右の錨)? Get_繋ぐ形(string p_u文字列, string p_x文字列, Dictionary<ulong, long> p_錨表, string p_ほかのk, ReadMinimizerIndex p_索引, int p_リード長)
        {
            var l_窓長 = Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_確かめる窓の分子 / C_確かめる窓の分母);
            var l_左群 = Get_錨の位置群(p_u文字列, p_錨表, true).Where(x => Isほかのkにある(p_ほかのk, p_u文字列.Substring(x, OtherKRepeatReplacer.C_錨長)));
            var l_右群 = Get_錨の位置群(p_x文字列, p_錨表, false).Where(x => Isほかのkにある(p_ほかのk, p_x文字列.Substring(x, OtherKRepeatReplacer.C_錨長))).ToList();
            foreach (var l_左 in l_左群)
            {
                var l_左の錨 = p_u文字列.Substring(l_左, OtherKRepeatReplacer.C_錨長);
                foreach (var l_右 in l_右群)
                {
                    var l_右の錨 = p_x文字列.Substring(l_右, OtherKRepeatReplacer.C_錨長);
                    HashSet<string> l_間群 = new(StringComparer.Ordinal);
                    OtherKGapFiller.V_集める_錨の間(p_ほかのk, l_左の錨, l_右の錨, false, l_間群);
                    OtherKGapFiller.V_集める_錨の間(p_ほかのk, Util.V_逆相補(l_右の錨), Util.V_逆相補(l_左の錨), true, l_間群);
                    if (l_間群.Count == 0)
                    {
                        continue;
                    }

                    var l_左の片 = p_u文字列[..(l_左 + OtherKRepeatReplacer.C_錨長)];
                    var l_右の片 = p_x文字列[l_右..];
                    var l_通った = l_間群.Where(x => PairBridgeFiller.Is写しがリードに続けて出る(p_索引, l_左の片, x, l_右の片, l_窓長)).Take(2).ToList();
                    return l_通った.Count == 1 ? (p_u文字列.Length - l_左の片.Length, l_通った[0], l_右, l_左, l_右) : null;
                }
            }

            return null;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 錨か、その逆相補が、ほかの k の連結配列にあるか<br/>
        /// 無い錨は、どの錨と組んでも錨の間を集められないので、組を試す前に除く
        /// </summary>
        /// <param name="p_ほかのk">ほかの k の contig を区切りで連結した配列</param>
        /// <param name="p_錨">錨</param>
        /// <returns>あれば真</returns>
        private static bool Isほかのkにある(string p_ほかのk, string p_錨)
        {
            return p_ほかのk.Contains(p_錨, StringComparison.Ordinal) || p_ほかのk.Contains(Util.V_逆相補(p_錨), StringComparison.Ordinal);
        }

        /// <summary>
        /// 配列の末尾 (または先頭) の C_錨を探す長さ の中で、アセンブリで 1 回だけ出る錨の位置を、端に近い順に C_錨の間隔 おきに C_試す錨の数 まで並べる<br/>
        /// 錨は、前後 C_錨長 - 1 bp の中で始まる 31-mer もすべて 1 回だけのものに限る。反復の境目に掛かる錨は、ほとんどが反復の配列なので、ほかの k の別の所に当たりうる
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_錨表">アセンブリの一意な錨の表</param>
        /// <param name="p_Is末尾">末尾から探すなら真、先頭から探すなら偽</param>
        /// <returns>錨の始まりの位置 (配列の中の位置)</returns>
        internal static List<int> Get_錨の位置群(string p_配列, Dictionary<ulong, long> p_錨表, bool p_Is末尾)
        {
            var l_始め = p_Is末尾 ? Math.Max(0, p_配列.Length - C_錨を探す長さ) : 0;
            var l_範囲 = p_配列.Substring(l_始め, Math.Min(C_錨を探す長さ, p_配列.Length));
            var l_一意 = new bool[l_範囲.Length];
            foreach (var (l_位置, l_正準, _) in OtherKRepeatReplacer.Get_錨群(l_範囲))
            {
                l_一意[l_位置] = p_錨表.TryGetValue(l_正準, out var l_場所) && l_場所 >= 0;
            }

            var l_余白 = OtherKRepeatReplacer.C_錨長 - 1;
            var l_深い = Enumerable.Range(l_余白, Math.Max(0, l_範囲.Length - OtherKRepeatReplacer.C_錨長 + 1 - (2 * l_余白))).Where(p => Enumerable.Range(p - l_余白, (2 * l_余白) + 1).All(i => l_一意[i]));
            List<int> l_位置群 = [];
            foreach (var l_位置 in p_Is末尾 ? l_深い.OrderByDescending(static x => x) : l_深い.OrderBy(static x => x))
            {
                if (l_位置群.Count > 0 && Math.Abs(l_始め + l_位置 - l_位置群[^1]) < C_錨の間隔)
                {
                    continue;
                }

                l_位置群.Add(l_始め + l_位置);
                if (l_位置群.Count == C_試す錨の数)
                {
                    break;
                }
            }

            return l_位置群;
        }

        #endregion
    }
}
