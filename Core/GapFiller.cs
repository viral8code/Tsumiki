using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// スキャフォールドの N で埋められたギャップを、de Bruijn グラフ上で
    /// 両側を繋ぐ経路を探して実際の塩基配列に置き換える。
    ///
    /// スキャフォールディングは「2つの contig がこの向きでこれくらい離れて
    /// 隣接している」というペアエンドの情報だけで繋いでおり、その間の配列が
    /// 何であるかは埋めていない。しかしギャップの中身は多くの場合 k-mer 集合の
    /// 中に存在している。contig が途切れたのは配列が無いからではなく、
    /// 「分岐があってどちらへ進むか決められなかった」ためであることが多いからで、
    /// その場合ギャップの両端を繋ぐ経路はグラフ上に実在する。
    ///
    /// 探索は左端の k-mer から1塩基ずつ伸ばす幅優先で行い、推定ギャップ長の
    /// 前後 余裕幅 の範囲で右端の k-mer に到達する経路を集める。
    /// 経路がちょうど1本に定まったときだけ埋める。複数見つかった場合は
    /// どれが正しいか決められないため N のまま残す(誤った配列で埋めるより、
    /// 分からないことが分かる状態のほうが下流の解析にとって安全)。
    /// </summary>
    internal static class GapFiller
    {
        /// <summary>
        /// 推定ギャップ長に対して許容する誤差(塩基)。インサートサイズ推定の
        /// ばらつきがそのままギャップ長推定のばらつきになるため、
        /// ぴったりの長さだけを探すと現実にはまず当たらない。
        /// </summary>
        private const int 長さの余裕幅 = 30;

        /// <summary>
        /// 1つのギャップあたりに展開してよい探索状態の上限。
        /// 分岐の多い領域では経路数が指数的に増えるため、上限を超えたら
        /// 「解けなかった」として諦める(時間をかけても曖昧なままのことが多い)。
        /// </summary>
        private const int ギャップあたりの状態数上限 = 200_000;

        /// <summary>
        /// これより長いギャップは探索空間が広すぎるうえ、推定長の誤差も大きく
        /// 一意に定まる見込みが薄いため対象外とする。
        /// </summary>
        private const int ギャップ長の上限 = 500;

        /// <summary>
        /// スキャフォールドを読み込み、埋められるギャップを埋めて同じパスへ書き戻す。
        /// </summary>
        public static ギャップ充填統計 V_充填_ギャップ(
            string p_スキャフォールドパス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            List<(string A_ID, string A_配列)> l_スキャフォールド群 = [];
            using (var l_読み込み = new FastaReader(p_スキャフォールドパス))
            {
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_配列エントリ = l_読み込み.Get_次の配列();
                    l_スキャフォールド群.Add((l_配列エントリ.A_ID.TrimStart('>'), l_配列エントリ.A_配列));
                }
            }

            var l_総ギャップ数 = 0;
            var l_埋めたギャップ数 = 0;
            var l_埋めた塩基数 = 0;
            var l_一意でない数 = 0;
            var l_到達不能数 = 0;

            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in l_スキャフォールド群)
            {
                var l_出力 = new StringBuilder();
                var l_位置 = 0;
                while (l_位置 < l_配列.Length)
                {
                    if (l_配列[l_位置] != 'N')
                    {
                        _ = l_出力.Append(l_配列[l_位置]);
                        l_位置++;
                        continue;
                    }

                    // N の連続区間 = 1つのギャップ。
                    var l_ギャップ開始 = l_位置;
                    while (l_位置 < l_配列.Length && l_配列[l_位置] == 'N')
                    {
                        l_位置++;
                    }
                    var l_ギャップ長 = l_位置 - l_ギャップ開始;
                    l_総ギャップ数++;

                    var l_埋めた配列 = Get_ギャップを埋める配列(
                        l_出力, l_配列, l_ギャップ長, l_位置, p_kmerインデックス, p_k長, out var l_判定);
                    if (l_埋めた配列 != null)
                    {
                        _ = l_出力.Append(l_埋めた配列);
                        l_埋めたギャップ数++;
                        l_埋めた塩基数 += l_埋めた配列.Length;
                    }
                    else
                    {
                        if (l_判定 == ギャップ充填判定.一意でない)
                        {
                            l_一意でない数++;
                        }
                        else
                        {
                            l_到達不能数++;
                        }
                        _ = l_出力.Append('N', l_ギャップ長);
                    }
                }
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            using (var l_書き込み = new FastaWriter(p_スキャフォールドパス))
            {
                foreach (var (l_ID, l_配列) in l_結果)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            return new ギャップ充填統計(l_総ギャップ数, l_埋めたギャップ数, l_埋めた塩基数, l_一意でない数, l_到達不能数);
        }

        /// <summary>
        /// ギャップの左右の足場から、その間を埋める配列を探す。
        /// 見つからない/一意に定まらない場合は null を返す。
        /// </summary>
        private static string? Get_ギャップを埋める配列(
            StringBuilder p_左側の出力,
            string p_配列,
            int p_ギャップ長,
            int p_ギャップ終端,
            TrustedKmerIndex p_kmerインデックス,
            int p_k長,
            out ギャップ充填判定 p_判定)
        {
            p_判定 = ギャップ充填判定.到達不能;

            if (p_ギャップ長 > ギャップ長の上限 || p_左側の出力.Length < p_k長)
            {
                return null;
            }
            if (p_ギャップ終端 + p_k長 > p_配列.Length)
            {
                return null;
            }

            // 左側の足場: 既に書き出した配列の末尾 k-mer。
            var l_左のkmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_左のkmer[i] = Util.Get_塩基ID(p_左側の出力[p_左側の出力.Length - p_k長 + i]);
            }

            // 右側の足場: ギャップ直後の k-mer。ここへ到達できれば繋がったことになる。
            var l_目標kmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_目標kmer[i] = Util.Get_塩基ID(p_配列[p_ギャップ終端 + i]);
            }

            if (Array.IndexOf(l_左のkmer, Consts.無効な塩基) >= 0 || Array.IndexOf(l_目標kmer, Consts.無効な塩基) >= 0)
            {
                return null;
            }
            if (!p_kmerインデックス.Get_含まれるか(l_左のkmer) || !p_kmerインデックス.Get_含まれるか(l_目標kmer))
            {
                // 足場そのものが信頼できる k-mer 集合に無いなら探索しても意味がない。
                return null;
            }

            var l_最小長 = Math.Max(0, p_ギャップ長 - 長さの余裕幅);
            var l_最大長 = p_ギャップ長 + 長さの余裕幅;

            // 幅優先で1塩基ずつ伸ばす。
            //
            // 各状態が「これまでに継ぎ足した塩基列」そのものを持つと、
            // 状態数の上限 × 経路長ぶんのメモリと文字列コピーが発生する。
            // 代わりに親へのインデックスと追加した1塩基だけを持ち、
            // 解が見つかったときに親を辿って復元する。1状態あたり定数サイズで済む。
            var l_節点 = new List<(int A_親, byte A_塩基)>(1024) { (-1, 0) };
            var l_kmer群 = new List<byte[]>(1024) { l_左のkmer };
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

                // 継ぎ足した数は「左の足場 k-mer の後ろに継ぎ足した塩基数」。目標 k-mer に
                // 到達した時点では、その末尾 k長 塩基が目標 k-mer 自身に
                // あたる(それは元の配列に既にある)ので、ギャップを実際に埋める
                // 長さは 継ぎ足した数 - k長 になる。
                // 打ち切りもこの「埋める長さ」で判断しないと、正解の経路を
                // 目標到達の直前で切ってしまう。
                var l_埋める長さ = l_継ぎ足した数 - p_k長;
                if (l_埋める長さ > l_最大長)
                {
                    continue;
                }

                if (l_埋める長さ >= l_最小長 && l_現在のkmer.AsSpan().SequenceEqual(l_目標kmer))
                {
                    l_見つかった経路.Add(Get_復元経路(l_節点, l_現在, l_埋める長さ));
                    if (l_見つかった経路.Count > 1)
                    {
                        // 2本見つかった時点で一意には定まらない。
                        p_判定 = ギャップ充填判定.一意でない;
                        return null;
                    }
                    continue;
                }

                if (l_節点.Count > ギャップあたりの状態数上限)
                {
                    p_判定 = ギャップ充填判定.一意でない;
                    return null;
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

            if (l_見つかった経路.Count == 1)
            {
                p_判定 = ギャップ充填判定.充填済み;
                return l_見つかった経路[0];
            }

            p_判定 = l_見つかった経路.Count > 1 ? ギャップ充填判定.一意でない : ギャップ充填判定.到達不能;
            return null;
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

        public static void V_出力_充填統計(ギャップ充填統計 p_統計)
        {
            if (p_統計.A_総ギャップ数 == 0)
            {
                Console.WriteLine("[Info] Gap filling: no gaps to fill.");
                return;
            }
            Console.WriteLine(
                $"[Info] Gap filling: {p_統計.A_埋めたギャップ数}/{p_統計.A_総ギャップ数} gap(s) closed with real sequence " +
                $"({p_統計.A_埋めた塩基数:N0}bp of N replaced). " +
                $"{p_統計.A_一意に定まらなかった数} left as N because more than one path fits, " +
                $"{p_統計.A_到達できなかった数} because no path through the graph connects the two sides.");
        }
    }
}
