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

        #region 内部変数

        /// <summary>
        /// スレッドごとに使い回すパック値探索の作業域
        /// </summary>
        [ThreadStatic]
        private static 経路探索作業域? _作業域;

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
            return p_k長 <= TrustedKmerIndex.パック値のk上限
                ? Get_経路_パック値(p_左のkmer, p_目標kmer, p_最小長, p_最大長, p_kmerインデックス, p_k長, p_状態数上限)
                : Get_経路_参照(p_左のkmer, p_目標kmer, p_最小長, p_最大長, p_kmerインデックス, p_k長, p_状態数上限);
        }

        /// <summary>
        /// 塩基列のまま状態を持つ探索 (パック値に収まらない k と、パック値版の検証用)
        /// </summary>
        /// <param name="p_左のkmer"></param>
        /// <param name="p_目標kmer"></param>
        /// <param name="p_最小長"></param>
        /// <param name="p_最大長"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_状態数上限"></param>
        /// <returns></returns>
        internal static (string? A_経路, ギャップ充填判定 A_判定) Get_経路_参照(byte[] p_左のkmer, byte[] p_目標kmer, int p_最小長, int p_最大長, IKmerLookup p_kmerインデックス, int p_k長, int p_状態数上限 = 既定状態数上限)
        {
            // 状態は親へのインデックスと追加した 1 塩基だけを持ち、解が見つかったときに親を辿って復元する
            var l_節点 = new List<(int A_親, byte A_塩基)>(1_024) { (-1, 0) };
            var l_kmer群 = new List<byte[]>(1_024) { p_左のkmer };
            var l_深さ群 = new List<int>(1_024) { 0 };

            // 同じ k-mer に同じ深さで別経路から着いた状態は作らず、重ねて到達された印を付ける
            // 目標へ届いた経路がその印を通っていたら一意でないとみなす
            var l_到達済み = new Dictionary<(string, int), int>(1_024);
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

                // 目標に着いた時点の末尾 k 塩基は目標 k-mer 自身なので、新しく埋まる長さは 継ぎ足した数 - k になり、打ち切りもこの長さで判断する
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
                        return (null, ギャップ充填判定.一意でない);
                    }
                    continue;
                }

                // 伸ばした先はどれも長さの上限を超えて捨てられる
                if (l_埋める長さ >= p_最大長)
                {
                    continue;
                }

                if (l_節点.Count > p_状態数上限)
                {
                    return (null, ギャップ充填判定.探索打切り);
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
                    var l_鍵 = (Get_状態の鍵(l_作業バッファ), l_深さ);
                    if (l_到達済み.TryGetValue(l_鍵, out var l_既存))
                    {
                        l_多重到達[l_既存] = true;
                        continue;
                    }
                    l_到達済み[l_鍵] = l_節点.Count;
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
        /// 状態の k-mer を順鎖・逆鎖のパック値で持ち、1 塩基ずつ転がして探す (k &lt;= 128)
        /// </summary>
        /// <param name="p_左のkmer"></param>
        /// <param name="p_目標kmer"></param>
        /// <param name="p_最小長"></param>
        /// <param name="p_最大長"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_状態数上限"></param>
        /// <remarks>
        /// 探索の順序と判定は Get_経路_参照 と同じ<br/>
        /// 塩基列のまま持つと、展開のたびに候補ごとの詰め直しと状態ごとの配列確保が走り、リードペアの数だけ繰り返す橋渡しでそれが支配的になる
        /// </remarks>
        /// <returns></returns>
        private static (string? A_経路, ギャップ充填判定 A_判定) Get_経路_パック値(byte[] p_左のkmer, byte[] p_目標kmer, int p_最小長, int p_最大長, IKmerLookup p_kmerインデックス, int p_k長, int p_状態数上限)
        {
            var l_下位マスク = p_k長 >= 64 ? UInt128.MaxValue : ((UInt128)1 << (2 * p_k長)) - 1;
            var l_上位マスク = p_k長 <= 64 ? 0 : p_k長 >= 128 ? UInt128.MaxValue : ((UInt128)1 << ((2 * p_k長) - 128)) - 1;
            var l_先頭シフト = 2 * (p_k長 - 1);

            var l_目標 = TrustedKmerIndex.TryGet_パック_長(p_目標kmer);
            var l_左順 = TrustedKmerIndex.TryGet_パック_長(p_左のkmer);
            var l_左逆 = Get_逆相補パック(p_左のkmer);

            var l_作業域 = Get_作業域();
            var l_節点 = l_作業域.A_節点;
            var l_状態群 = l_作業域.A_状態群;
            var l_深さ群 = l_作業域.A_深さ群;
            var l_到達済み = l_作業域.A_到達済み;
            var l_多重到達 = l_作業域.A_多重到達;
            var l_キュー = l_作業域.A_キュー;
            l_節点.Add((-1, 0));
            l_状態群.Add((l_左順.A_上位, l_左順.A_下位, l_左逆.A_上位, l_左逆.A_下位));
            l_深さ群.Add(0);
            l_多重到達.Add(false);
            l_キュー.Enqueue(0);

            var l_見つかった経路 = new List<string>();

            while (l_キュー.Count > 0)
            {
                var l_現在 = l_キュー.Dequeue();
                var (l_順上, l_順下, l_逆上, l_逆下) = l_状態群[l_現在];
                var l_継ぎ足した数 = l_深さ群[l_現在];

                var l_埋める長さ = l_継ぎ足した数 - p_k長;
                if (l_埋める長さ > p_最大長)
                {
                    continue;
                }

                if (l_埋める長さ >= p_最小長 && l_順上 == l_目標.A_上位 && l_順下 == l_目標.A_下位)
                {
                    if (Has多重到達(l_節点, l_多重到達, l_現在))
                    {
                        return (null, ギャップ充填判定.一意でない);
                    }

                    l_見つかった経路.Add(Get_復元経路(l_節点, l_現在, l_埋める長さ));
                    if (l_見つかった経路.Count > 1)
                    {
                        return (null, ギャップ充填判定.一意でない);
                    }
                    continue;
                }

                // 伸ばした先はどれも長さの上限を超えて捨てられる
                if (l_埋める長さ >= p_最大長)
                {
                    continue;
                }

                if (l_節点.Count > p_状態数上限)
                {
                    return (null, ギャップ充填判定.探索打切り);
                }

                // 左へ 1 塩基ずらした後続の共通部分は塩基によらないので先に作る
                var l_ずらし順上 = ((l_順上 << 2) | (l_順下 >> 126)) & l_上位マスク;
                var l_ずらし順下 = (l_順下 << 2) & l_下位マスク;
                var l_ずらし逆下 = (l_逆下 >> 2) | (l_逆上 << 126);
                var l_ずらし逆上 = l_逆上 >> 2;

                for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
                {
                    var l_コドン = (UInt128)(l_塩基 - 1);
                    var l_新順下 = l_ずらし順下 | l_コドン;
                    var l_新逆上 = l_ずらし逆上;
                    var l_新逆下 = l_ずらし逆下;
                    if (l_先頭シフト >= 128)
                    {
                        l_新逆上 |= (3 - l_コドン) << (l_先頭シフト - 128);
                    }
                    else
                    {
                        l_新逆下 |= (3 - l_コドン) << l_先頭シフト;
                    }

                    var l_Is順鎖 = l_ずらし順上 < l_新逆上 || (l_ずらし順上 == l_新逆上 && l_新順下 <= l_新逆下);
                    if (!(l_Is順鎖 ? p_kmerインデックス.Haskmer_正規形(l_ずらし順上, l_新順下) : p_kmerインデックス.Haskmer_正規形(l_新逆上, l_新逆下)))
                    {
                        continue;
                    }

                    var l_深さ = l_継ぎ足した数 + 1;
                    var l_鍵 = (l_ずらし順上, l_新順下, l_深さ);
                    if (l_到達済み.TryGetValue(l_鍵, out var l_既存))
                    {
                        l_多重到達[l_既存] = true;
                        continue;
                    }
                    l_到達済み[l_鍵] = l_節点.Count;
                    l_節点.Add((l_現在, l_塩基));
                    l_状態群.Add((l_ずらし順上, l_新順下, l_新逆上, l_新逆下));
                    l_深さ群.Add(l_深さ);
                    l_多重到達.Add(false);
                    l_キュー.Enqueue(l_節点.Count - 1);
                }
            }

            return l_見つかった経路.Count == 1
                ? (l_見つかった経路[0], ギャップ充填判定.充填済み)
                : (null, l_見つかった経路.Count > 1 ? ギャップ充填判定.一意でない : ギャップ充填判定.到達不能);
        }

        /// <summary>
        /// このスレッドの作業域を空にして返す
        /// </summary>
        /// <returns></returns>
        private static 経路探索作業域 Get_作業域()
        {
            if (_作業域 is not { } l_作業域 || l_作業域.A_節点.Count > 経路探索作業域.作り直す状態数)
            {
                _作業域 = new 経路探索作業域();
                return _作業域;
            }
            l_作業域.V_初期化();
            return l_作業域;
        }

        /// <summary>
        /// k-mer の逆相補を右詰めでパックする (k &lt;= 128)
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static (UInt128 A_上位, UInt128 A_下位) Get_逆相補パック(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_上位 = 0;
            UInt128 l_下位 = 0;
            for (var i = p_kmer.Length - 1; i >= 0; i--)
            {
                l_上位 = (l_上位 << 2) | (l_下位 >> 126);
                l_下位 = (l_下位 << 2) | (UInt128)(4 - p_kmer[i]);
            }
            return (l_上位, l_下位);
        }

        /// <summary>
        /// k-mer を鍵にするための文字列表現
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
