using Tsumiki.Commons;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 信頼できる k-mer 集合の中で、2 つの k-mer を繋ぐ経路を幅優先探索で探す
    /// </summary>
    /// <remarks>
    /// 経路がちょうど 1 本、かつ追加した塩基数が指定範囲に収まるときだけ結果を返す<br/>
    /// 複数見つかった、あるいは 1 本も見つからない場合はどれが正しいか決められないため null を返す<br/>
    /// 誤った配列で埋めるより、分からないことが分かる状態のほうが下流の解析にとって安全、という方針そのものは呼び出し側の目的 (ギャップ充填/リードペアの橋渡し) によらず共通なので、探索エンジンをここへ切り出している
    /// </remarks>
    internal static class ConstrainedPathFinder
    {
        #region 定数

        /// <summary>
        /// 展開してよい探索状態の上限
        /// </summary>
        /// <remarks>
        /// 分岐の多い領域では経路数が指数的に増えるため、上限を超えたら「解けなかった」として諦める (時間をかけても曖昧なままのことが多い)
        /// </remarks>
        public const int 既定状態数上限 = 200_000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_左のkmer から 1 塩基ずつ伸ばし、p_目標kmer に一致する状態のうち、追加した塩基数 (=末尾 k-mer 同士の重なりを除いた、新たに埋まる長さ) が[p_最小長, p_最大長] に収まるものを探す
        /// </summary>
        /// <param name="p_左のkmer"></param>
        /// <param name="p_目標kmer"></param>
        /// <param name="p_最小長"></param>
        /// <param name="p_最大長"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_状態数上限"></param>
        /// <returns></returns>
        public static (string? A_経路, ギャップ充填判定 A_判定) Get_経路(byte[] p_左のkmer, byte[] p_目標kmer, int p_最小長, int p_最大長, IKmerLookup p_kmerインデックス, int p_k長, int p_状態数上限 = 既定状態数上限)
        {
            // 各状態が「これまでに継ぎ足した塩基列」そのものを持つと、
            // 状態数の上限 × 経路長ぶんのメモリと文字列コピーが発生する
            // 代わりに親へのインデックスと追加した 1 塩基だけを持ち、
            // 解が見つかったときに親を辿って復元する
            // 1 状態あたり定数サイズで済む
            var l_節点 = new List<(int A_親, byte A_塩基)>(1_024) { (-1, 0) };
            var l_kmer群 = new List<byte[]>(1_024) { p_左のkmer };
            var l_深さ群 = new List<int>(1_024) { 0 };

            // 同じ k-mer に同じ深さで別経路から着いた状態を作らないための記録
            // これをしないと分岐点の数だけ経路が掛け算で増え、同じ部分木を
            // 何度も展開する (反復配列では容易に状態数上限に達する)
            // 一意性の判定は捨てられないので、重ねて到達された状態には印を付け、
            // 目標へ届いた経路がその印を通っていたら一意でないとみなす
            // k <= 64 なら k-mer は UInt128 に詰められるので、鍵の生成に
            // 割り当てが要らない
            // それを超える k だけ文字列に落とす
            var l_到達済み = p_k長 <= 64 ? new Dictionary<(UInt128, int), int>(1_024) : null;
            var l_到達済み_長いk = p_k長 <= 64 ? null : new Dictionary<(string, int), int>(1_024);
            var l_多重到達 = new List<bool>(1_024) { false };

            var l_見つかった経路 = new List<string>();
            var l_キュー = new Queue<int>();
            l_キュー.Enqueue(0);

            var l_作業バッファ = new byte[p_k長];

            while (l_キュー.Count > 0)
            {
                var l_現在 = l_キュー.Dequeue();
                var l_現在のkmer = l_kmer群[l_現在];
                var l_継ぎ足した数 = l_深さ群[l_現在];

                // 継ぎ足した数は「左の k-mer の後ろに継ぎ足した塩基数」
                // 目標 k-mer に
                // 到達した時点では、その末尾 k 長 塩基が目標 k-mer 自身に
                // あたる (呼び出し側が既に知っている) ので、実際に新しく埋まる
                // 長さは 継ぎ足した数 - k 長 になる
                // 打ち切りもこの「埋める長さ」で
                // 判断しないと、正解の経路を目標到達の直前で切ってしまう
                var l_埋める長さ = l_継ぎ足した数 - p_k長;
                if (l_埋める長さ > p_最大長)
                {
                    continue;
                }

                if (l_埋める長さ >= p_最小長 && l_現在のkmer.AsSpan().SequenceEqual(p_目標kmer))
                {
                    if (Has多重到達(l_節点, l_多重到達, l_現在))
                    {
                        return (null, ギャップ充填判定.一意でない);
                    }

                    l_見つかった経路.Add(Get_復元経路(l_節点, l_現在, l_埋める長さ));
                    if (l_見つかった経路.Count > 1)
                    {
                        // 2 本見つかった時点で一意には定まらない
                        return (null, ギャップ充填判定.一意でない);
                    }
                    continue;
                }

                if (l_節点.Count > p_状態数上限)
                {
                    return (null, ギャップ充填判定.一意でない);
                }

                for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
                {
                    Array.Copy(l_現在のkmer, 1, l_作業バッファ, 0, p_k長 - 1);
                    l_作業バッファ[p_k長 - 1] = l_塩基;

                    if (!p_kmerインデックス.Haskmer(l_作業バッファ))
                    {
                        continue;
                    }

                    var l_深さ = l_継ぎ足した数 + 1;
                    if (l_到達済み is not null)
                    {
                        var l_鍵 = (TryGet_パック(l_作業バッファ), l_深さ);
                        if (l_到達済み.TryGetValue(l_鍵, out var l_既存))
                        {
                            l_多重到達[l_既存] = true;
                            continue;
                        }
                        l_到達済み[l_鍵] = l_節点.Count;
                    }
                    else
                    {
                        var l_鍵 = (Get_状態の鍵(l_作業バッファ), l_深さ);
                        if (l_到達済み_長いk!.TryGetValue(l_鍵, out var l_既存))
                        {
                            l_多重到達[l_既存] = true;
                            continue;
                        }
                        l_到達済み_長いk[l_鍵] = l_節点.Count;
                    }
                    l_節点.Add((l_現在, l_塩基));
                    l_kmer群.Add((byte[])l_作業バッファ.Clone());
                    l_深さ群.Add(l_深さ);
                    l_多重到達.Add(false);
                    l_キュー.Enqueue(l_節点.Count - 1);
                }
            }

            return l_見つかった経路.Count == 1
                ? (l_見つかった経路[0], ギャップ充填判定.充填済み)
                : (null, l_見つかった経路.Count > 1 ? ギャップ充填判定.一意でない : ギャップ充填判定.到達不能);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer を 1 塩基 2 ビットで詰める (k &lt;= 64 でのみ使える)
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static UInt128 TryGet_パック(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_パック = 0;
            foreach (var l_塩基 in p_kmer)
            {
                l_パック = (l_パック << 2) | (uint)(l_塩基 - 1);
            }
            return l_パック;
        }

        /// <summary>
        /// k &gt; 64 で k-mer を鍵にするための文字列表現
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static string Get_状態の鍵(ReadOnlySpan<byte> p_kmer)
        {
            var l_文字 = new char[p_kmer.Length];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                l_文字[i] = (char)p_kmer[i];
            }
            return new string(l_文字);
        }

        /// <summary>
        /// この状態までの経路上に、別経路からも到達された状態があるか
        /// </summary>
        /// <param name="p_節点"></param>
        /// <param name="p_多重到達"></param>
        /// <param name="p_末端"></param>
        /// <returns></returns>
        private static bool Has多重到達(List<(int A_親, byte A_塩基)> p_節点, List<bool> p_多重到達, int p_末端)
        {
            for (var l_位置 = p_末端; l_位置 >= 0; l_位置 = p_節点[l_位置].A_親)
            {
                if (p_多重到達[l_位置])
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// 親を辿って、継ぎ足した塩基列のうち先頭 p_埋める長さ 塩基を復元する
        /// </summary>
        /// <param name="p_節点"></param>
        /// <param name="p_末端"></param>
        /// <param name="p_埋める長さ"></param>
        /// <remarks>
        /// 末尾側 (目標 k-mer と重なる分) は捨てる
        /// </remarks>
        /// <returns></returns>
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

        #endregion
    }
}
