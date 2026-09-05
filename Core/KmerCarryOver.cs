using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 前段の k で組み上がった配列を、次の k の k-mer 集合へ引き継ぐ。
    ///
    /// k を上げるとカバレッジが痩せてグラフが千切れるが、前段の配列は
    /// その領域を既に通っている。配列を渡せば連結が保たれる。
    ///
    /// 渡すのは配列であって、繋ぐという決定ではない。決定を渡すと前段の
    /// 誤アセンブリをそのまま継承するが、配列を渡すだけなら次の k が
    /// 自分の証拠で経路を決め直せる。
    /// </summary>
    internal static class KmerCarryOver
    {
        /// <summary>
        /// 引き継ぐ配列の最小長。前段で短く切れた断片は連結の役に立たないうえ、
        /// エラー由来の残骸である可能性が相対的に高い。
        /// </summary>
        private const int 引き継ぐ配列の最小長 = 500;

        /// <summary>
        /// 引き継ぎ元の配列とカバレッジを、その k の成果物から作る。
        /// k-mer インデックスが生きているうちにしか作れない。
        /// </summary>
        public static List<引き継ぎ配列> Get_引き継ぎ配列(
            string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            List<引き継ぎ配列> l_結果 = [];
            using var l_読み込み = new FastaReader(p_FASTAパス);

            while (l_読み込み.Get_続きがあるか())
            {
                var l_配列 = l_読み込み.Get_次の配列().A_配列;
                if (l_配列.Length < Math.Max(引き継ぐ配列の最小長, p_k長))
                {
                    continue;
                }

                var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                var l_カバレッジ = new int[l_配列.Length - p_k長 + 1];
                for (var i = 0; i < l_カバレッジ.Length; i++)
                {
                    l_カバレッジ[i] = (int)Math.Min(
                        int.MaxValue, p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長)));
                }
                l_結果.Add(new 引き継ぎ配列(l_配列, l_カバレッジ, p_k長));
            }
            return l_結果;
        }

        /// <summary>
        /// 引き継ぎ配列のうち、この k の集合に無い k-mer を足す。
        /// 既にある k-mer は触らない(実際のリード由来の観測を優先する)。
        /// 戻り値は足した k-mer の数。
        /// </summary>
        public static int V_引き継ぎ(
            IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列,
            TrustedKmerIndex p_kmerインデックス,
            int p_k長,
            int? p_リード長)
        {
            var l_追加数 = 0;
            foreach (var l_引き継ぎ in p_引き継ぎ配列)
            {
                if (l_引き継ぎ.A_配列.Length < p_k長)
                {
                    continue;
                }

                var l_塩基列 = l_引き継ぎ.A_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    if (Array.IndexOf(l_塩基列, Consts.無効な塩基, i, p_k長) >= 0)
                    {
                        continue;
                    }
                    var l_カバレッジ = Get_引き継ぐカバレッジ(l_引き継ぎ, i, p_k長, p_リード長);
                    if (l_カバレッジ > 0
                        && p_kmerインデックス.V_追加_信頼kmer(l_塩基列.AsSpan(i, p_k長), l_カバレッジ))
                    {
                        l_追加数++;
                    }
                }
            }
            return l_追加数;
        }

        /// <summary>
        /// この k-mer に与えるカバレッジ。
        ///
        /// 前段の k-mer のうちこの窓に重なるものの最小値を取る。長い k-mer は
        /// 構成する短い k-mer すべてを含むので、最も弱い部分より強くはなれない。
        ///
        /// そのうえで k の差ぶんスケールする。1リードから取れる k-mer は
        /// リード長 - k + 1 本なので、k を上げれば同じ座位のカバレッジは
        /// その比で下がる。スケールしないと、引き継いだ領域だけカバレッジが
        /// 高く見えてコピー数を過大に推定する。
        /// </summary>
        public static ulong Get_引き継ぐカバレッジ(
            引き継ぎ配列 p_引き継ぎ, int p_位置, int p_k長, int? p_リード長)
        {
            var l_終端 = Math.Min(p_引き継ぎ.A_カバレッジ.Length - 1, p_位置 + p_k長 - p_引き継ぎ.A_k長);
            if (p_位置 > l_終端)
            {
                return 0;
            }

            var l_最小 = int.MaxValue;
            for (var i = p_位置; i <= l_終端; i++)
            {
                l_最小 = Math.Min(l_最小, p_引き継ぎ.A_カバレッジ[i]);
            }
            if (l_最小 <= 0)
            {
                return 0;
            }

            if (p_リード長 is not { } l_リード長 || l_リード長 <= p_k長)
            {
                return (ulong)l_最小;
            }

            var l_前段の本数 = l_リード長 - p_引き継ぎ.A_k長 + 1;
            var l_今の本数 = l_リード長 - p_k長 + 1;
            if (l_前段の本数 <= 0 || l_今の本数 <= 0)
            {
                return (ulong)l_最小;
            }
            return (ulong)Math.Max(1, (long)Math.Round((double)l_最小 * l_今の本数 / l_前段の本数));
        }
    }
}
