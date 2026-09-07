using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// 信頼できる k-mer 集合の中で、2つの k-mer を繋ぐ経路を幅優先探索で探す。
    /// 経路がちょうど1本、かつ追加した塩基数が指定範囲に収まるときだけ結果を返す。
    /// 複数見つかった、あるいは1本も見つからない場合はどれが正しいか決められない
    /// ため null を返す。誤った配列で埋めるより、分からないことが分かる状態のほうが
    /// 下流の解析にとって安全、という方針そのものは呼び出し側の目的(ギャップ充填/
    /// リードペアの橋渡し)によらず共通なので、探索エンジンをここへ切り出している。
    /// </summary>
    internal static class ConstrainedPathFinder
    {
        /// <summary>
        /// 展開してよい探索状態の上限。分岐の多い領域では経路数が指数的に
        /// 増えるため、上限を超えたら「解けなかった」として諦める
        /// (時間をかけても曖昧なままのことが多い)。
        /// </summary>
        public const int 既定状態数上限 = 200_000;

        /// <summary>
        /// p_左のkmer から1塩基ずつ伸ばし、p_目標kmer に一致する状態のうち、
        /// 追加した塩基数(=末尾 k-mer 同士の重なりを除いた、新たに埋まる長さ)が
        /// [p_最小長, p_最大長] に収まるものを探す。
        /// </summary>
        public static (string? A_経路, ギャップ充填判定 A_判定) Get_経路(
            byte[] p_左のkmer, byte[] p_目標kmer, int p_最小長, int p_最大長,
            TrustedKmerIndex p_kmerインデックス, int p_k長, int p_状態数上限 = 既定状態数上限)
        {
            // 各状態が「これまでに継ぎ足した塩基列」そのものを持つと、
            // 状態数の上限 × 経路長ぶんのメモリと文字列コピーが発生する。
            // 代わりに親へのインデックスと追加した1塩基だけを持ち、
            // 解が見つかったときに親を辿って復元する。1状態あたり定数サイズで済む。
            var l_節点 = new List<(int A_親, byte A_塩基)>(1024) { (-1, 0) };
            var l_kmer群 = new List<byte[]>(1024) { p_左のkmer };
            var l_深さ群 = new List<int>(1024) { 0 };

            var l_見つかった経路 = new List<string>();
            var l_キュー = new Queue<int>();
            l_キュー.Enqueue(0);

            var l_作業バッファ = new byte[p_k長];

            while (l_キュー.Count > 0)
            {
                var l_現在 = l_キュー.Dequeue();
                var l_現在のkmer = l_kmer群[l_現在];
                var l_継ぎ足した数 = l_深さ群[l_現在];

                // 継ぎ足した数は「左の k-mer の後ろに継ぎ足した塩基数」。目標 k-mer に
                // 到達した時点では、その末尾 k長 塩基が目標 k-mer 自身に
                // あたる(呼び出し側が既に知っている)ので、実際に新しく埋まる
                // 長さは 継ぎ足した数 - k長 になる。打ち切りもこの「埋める長さ」で
                // 判断しないと、正解の経路を目標到達の直前で切ってしまう。
                var l_埋める長さ = l_継ぎ足した数 - p_k長;
                if (l_埋める長さ > p_最大長)
                {
                    continue;
                }

                if (l_埋める長さ >= p_最小長 && l_現在のkmer.AsSpan().SequenceEqual(p_目標kmer))
                {
                    l_見つかった経路.Add(Get_復元経路(l_節点, l_現在, l_埋める長さ));
                    if (l_見つかった経路.Count > 1)
                    {
                        // 2本見つかった時点で一意には定まらない。
                        return (null, ギャップ充填判定.一意でない);
                    }
                    continue;
                }

                if (l_節点.Count > p_状態数上限)
                {
                    return (null, ギャップ充填判定.一意でない);
                }

                for (byte l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
                {
                    Array.Copy(l_現在のkmer, 1, l_作業バッファ, 0, p_k長 - 1);
                    l_作業バッファ[p_k長 - 1] = l_塩基;
                    if (!p_kmerインデックス.Get_含まれるか(l_作業バッファ))
                    {
                        continue;
                    }
                    l_節点.Add((l_現在, l_塩基));
                    l_kmer群.Add((byte[])l_作業バッファ.Clone());
                    l_深さ群.Add(l_継ぎ足した数 + 1);
                    l_キュー.Enqueue(l_節点.Count - 1);
                }
            }

            return l_見つかった経路.Count == 1
                ? (l_見つかった経路[0], ギャップ充填判定.充填済み)
                : (null, l_見つかった経路.Count > 1 ? ギャップ充填判定.一意でない : ギャップ充填判定.到達不能);
        }

        /// <summary>
        /// 親を辿って、継ぎ足した塩基列のうち先頭 p_埋める長さ 塩基を復元する。
        /// 末尾側(目標 k-mer と重なる分)は捨てる。
        /// </summary>
        private static string Get_復元経路(List<(int A_親, byte A_塩基)> p_節点, int p_末端, int p_埋める長さ)
        {
            List<byte> l_逆順 = [];
            var l_位置 = p_末端;
            while (l_位置 > 0)
            {
                l_逆順.Add(p_節点[l_位置].A_塩基);
                l_位置 = p_節点[l_位置].A_親;
            }
            l_逆順.Reverse();
            return string.Concat(l_逆順.Take(p_埋める長さ).Select(Util.V_変換_塩基文字));
        }
    }
}
