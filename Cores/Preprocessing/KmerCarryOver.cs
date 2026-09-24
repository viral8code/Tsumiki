using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// 前段の k で組み上がった配列を、次の k の k-mer 集合へ引き継ぐ
    /// </summary>
    internal static class KmerCarryOver
    {
        #region 定数

        /// <summary>
        /// 引き継ぐ配列の最小長
        /// </summary>
        /// <remarks>
        /// 前段で短く切れた断片は連結の役に立たないうえ、エラー由来の残骸である可能性が相対的に高い
        /// </remarks>
        private const int 引き継ぐ配列の最小長 = 500;

        /// <summary>
        /// 持ち越しの裏付けを確かめる r-mer の長さ
        /// </summary>
        /// <remarks>
        /// k+1 では足りない。反復の別コピーが繋がった継ぎ目でも、短い窓なら片方のコピーを読んだリードだけで真になる<br/>
        /// リード長に収まり、偶然の一致がまず起きない長さとして 41 を使う
        /// </remarks>
        internal const int 持ち越し検証のr長 = 41;

        /// <summary>
        /// 持ち越さないと判断する、未観測の窓の連続数
        /// </summary>
        /// <remarks>
        /// 低カバレッジでは未観測の窓が散発するため、1 つでも外すと本来の穴埋めまで消える<br/>
        /// 周りのカバレッジとの比でも判定できそうに見えるが、内側が沈むのはカバレッジの谷だけでなく反復の別コピーへ乗り換えた偽の経路もそうなので、使えない
        /// </remarks>
        internal const int 未観測の連続の下限 = 14;

        /// <summary>
        /// 未観測の連続の位置を書き出すファイル名 (持ち越し元の FASTA と同じ場所に置く)
        /// </summary>
        private const string 未観測の連続の書き出し名 = "unobserved_runs.tsv";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 引き継ぎ元の配列とカバレッジを、その k の成果物から作る
        /// </summary>
        /// <param name="p_FASTAパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_分岐の継ぎ目">その k で分岐のある継ぎ目を通った辺 ((k+1)-mer、両向き)、無ければ null</param>
        /// <param name="p_検証器">渡すと、リードで観測されていない r-mer が続く範囲を調べ、そこを持ち越しから外す</param>
        /// <remarks>
        /// k-mer インデックスが生きているうちにしか作れない
        /// </remarks>
        /// <returns></returns>
        public static List<引き継ぎ配列> Get_引き継ぎ配列(string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長, HashSet<string>? p_分岐の継ぎ目 = null, RepeatRMerVerifier? p_検証器 = null)
        {
            List<引き継ぎ配列> l_結果 = [];
            var l_度数 = p_検証器 is null ? null : new 連続長の度数();
            using var l_読み込み = new FastaReader(p_FASTAパス);

            while (l_読み込み.Has続き())
            {
                var l_エントリ = l_読み込み.Get_次の配列();
                var l_配列 = l_エントリ.A_配列;
                if (l_配列.Length < Math.Max(引き継ぐ配列の最小長, p_k長))
                {
                    continue;
                }

                var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                var l_カバレッジ = new int[l_配列.Length - p_k長 + 1];
                for (var i = 0; i < l_カバレッジ.Length; i++)
                {
                    l_カバレッジ[i] = (int)Math.Min(int.MaxValue, p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長)));
                }
                var l_継ぎ目位置 = Get_継ぎ目位置(l_配列, p_分岐の継ぎ目, p_k長);
                l_結果.Add(new 引き継ぎ配列(l_配列, l_カバレッジ, p_k長, A_Is確定経路: true, A_分岐の継ぎ目位置: l_継ぎ目位置, A_未観測の連続範囲: Get_未観測の連続範囲(l_配列, p_検証器, l_継ぎ目位置, p_k長, l_度数, l_エントリ.A_ID.TrimStart('>'), l_カバレッジ)));
            }

            if (l_度数 is not null)
            {
                l_度数.V_出力(p_k長);
                l_度数.V_書き出し(Path.Combine(Path.GetDirectoryName(Path.GetFullPath(p_FASTAパス))!, 未観測の連続の書き出し名));
            }
            return l_結果;
        }

        /// <summary>
        /// 配列の中で、分岐のある継ぎ目の辺が現れる開始位置
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_分岐の継ぎ目">辺 ((k+1)-mer、両向き) の集合、無ければ null</param>
        /// <param name="p_k長">辺を作った k</param>
        /// <returns>1 つも無ければ null</returns>
        internal static IReadOnlyList<int>? Get_継ぎ目位置(string p_配列, HashSet<string>? p_分岐の継ぎ目, int p_k長)
        {
            if (p_分岐の継ぎ目 is not { Count: > 0 })
            {
                return null;
            }

            var l_辺長 = p_k長 + 1;
            var l_参照 = p_分岐の継ぎ目.GetAlternateLookup<ReadOnlySpan<char>>();
            List<int>? l_位置群 = null;
            for (var i = 0; i + l_辺長 <= p_配列.Length; i++)
            {
                if (l_参照.Contains(p_配列.AsSpan(i, l_辺長)))
                {
                    (l_位置群 ??= []).Add(i);
                }
            }
            return l_位置群;
        }

        /// <summary>
        /// リードで観測されていない r-mer が続く範囲
        /// </summary>
        /// <param name="p_配列">調べる配列</param>
        /// <param name="p_検証器">r-mer の検証器、無ければ null</param>
        /// <param name="p_継ぎ目位置">分岐のある継ぎ目の辺の開始位置、無ければ null</param>
        /// <param name="p_k長">継ぎ目の辺を作った k</param>
        /// <param name="p_度数">渡すと、下限未満も含めた全ての連続を数える</param>
        /// <param name="p_配列名">度数の明細に残す配列の名前</param>
        /// <param name="p_カバレッジ">k-mer ごとのカバレッジ、渡すと連続の内側と前後の中央値を明細に残す</param>
        /// <remarks>
        /// 下限未満の分布も見えるように、全連続を一度拾ってから下限で絞る
        /// </remarks>
        /// <returns>1 つも無ければ null</returns>
        internal static IReadOnlyList<(int A_開始, int A_終了)>? Get_未観測の連続範囲(string p_配列, RepeatRMerVerifier? p_検証器, IReadOnlyList<int>? p_継ぎ目位置 = null, int p_k長 = 0, 連続長の度数? p_度数 = null, string p_配列名 = "", int[]? p_カバレッジ = null)
        {
            if (p_検証器 is null)
            {
                return null;
            }

            List<(int A_開始, int A_終了)> l_全連続 = [];
            p_検証器.V_収集_未観測の連続範囲(p_配列, 1, l_全連続);

            List<(int A_開始, int A_終了)> l_範囲 = [];
            foreach (var l_連続 in l_全連続)
            {
                if (p_度数 is not null)
                {
                    var (l_内側, l_外側) = Get_内外のカバレッジ(l_連続, p_カバレッジ);
                    p_度数.V_加算(p_配列名, l_連続.A_開始, l_連続.A_終了, Is継ぎ目に掛かる(l_連続, p_継ぎ目位置, p_k長), l_内側, l_外側);
                }
                if (l_連続.A_終了 - l_連続.A_開始 + 1 >= 未観測の連続の下限)
                {
                    l_範囲.Add(l_連続);
                }
            }
            return l_範囲.Count > 0 ? l_範囲 : null;
        }

        /// <summary>
        /// 未観測の連続に掛かる k-mer と、その前後の k-mer のカバレッジの中央値
        /// </summary>
        /// <param name="p_連続">r-mer の窓の開始と終了</param>
        /// <param name="p_カバレッジ">k-mer ごとのカバレッジ、無ければ (0, 0)</param>
        /// <remarks>
        /// 未観測の連続がカバレッジの谷なのか、リードに無い経路なのかを後から調べるために残す
        /// </remarks>
        /// <returns></returns>
        private static (int A_内側, int A_外側) Get_内外のカバレッジ((int A_開始, int A_終了) p_連続, int[]? p_カバレッジ)
        {
            if (p_カバレッジ is not { Length: > 0 })
            {
                return (0, 0);
            }

            // k が r より長いと k-mer の窓は r-mer の窓より少ないので、連続の端が配列の外に出うる
            const int l_前後の幅 = 150;
            var l_長さ = p_カバレッジ.Length;
            var l_開始 = Math.Min(p_連続.A_開始, l_長さ - 1);
            var l_終端 = Math.Clamp(p_連続.A_終了 + 持ち越し検証のr長 - 1, l_開始, l_長さ - 1);
            var l_内側 = p_カバレッジ[l_開始..(l_終端 + 1)];
            var l_前 = p_カバレッジ[Math.Max(0, l_開始 - l_前後の幅)..l_開始];
            var l_後 = p_カバレッジ[(l_終端 + 1)..Math.Min(l_長さ, l_終端 + 1 + l_前後の幅)];
            return (Get_中央値(l_内側), Get_中央値([.. l_前, .. l_後]));
        }

        /// <summary>
        /// 中央値 (空なら 0)
        /// </summary>
        /// <param name="p_値群"></param>
        /// <returns></returns>
        private static int Get_中央値(int[] p_値群)
        {
            if (p_値群.Length == 0)
            {
                return 0;
            }
            var l_並び = p_値群.Order().ToArray();
            return l_並び[l_並び.Length / 2];
        }

        /// <summary>
        /// 未観測の連続が、分岐のある継ぎ目の辺と塩基で重なるか
        /// </summary>
        /// <param name="p_連続">r-mer の窓の開始と終了</param>
        /// <param name="p_継ぎ目位置">継ぎ目の辺 ((k+1)-mer) の開始位置</param>
        /// <param name="p_k長">継ぎ目の辺を作った k</param>
        /// <returns></returns>
        private static bool Is継ぎ目に掛かる((int A_開始, int A_終了) p_連続, IReadOnlyList<int>? p_継ぎ目位置, int p_k長)
        {
            if (p_継ぎ目位置 is null)
            {
                return false;
            }

            var l_塩基の終端 = p_連続.A_終了 + 持ち越し検証のr長 - 1;
            foreach (var l_位置 in p_継ぎ目位置)
            {
                if (l_位置 <= l_塩基の終端 && l_位置 + p_k長 >= p_連続.A_開始)
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// リードで観測されていない範囲に掛かる窓の最小カバレッジを 0 にして、足さないようにする
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_k長">足す先の k</param>
        /// <param name="p_最小値列">窓ごとの最小カバレッジ</param>
        /// <remarks>
        /// 範囲に少しでも掛かる k-mer を落とす。継ぎ目そのものを含む k-mer だけを外しても、その両隣から同じ経路が繋がる
        /// </remarks>
        internal static void V_除外_未観測の範囲(引き継ぎ配列 p_引き継ぎ, int p_k長, Span<int> p_最小値列)
        {
            if (p_引き継ぎ.A_未観測の連続範囲 is not { Count: > 0 } l_範囲群)
            {
                return;
            }

            foreach (var (l_開始, l_終了) in l_範囲群)
            {
                var l_下 = Math.Max(0, l_開始 - p_k長 + 1);
                var l_上 = Math.Min(p_最小値列.Length - 1, l_終了 + 持ち越し検証のr長 - 1);
                for (var s = l_下; s <= l_上; s++)
                {
                    p_最小値列[s] = 0;
                }
            }
        }

        /// <summary>
        /// 分岐のある継ぎ目の辺を丸ごと含む窓の最小カバレッジを 0 にして、足さないようにする
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_k長">足す先の k</param>
        /// <param name="p_最小値列">窓ごとの最小カバレッジ</param>
        internal static void V_除外_継ぎ目を含む窓(引き継ぎ配列 p_引き継ぎ, int p_k長, Span<int> p_最小値列)
        {
            if (p_引き継ぎ.A_分岐の継ぎ目位置 is not { Count: > 0 } l_位置群)
            {
                return;
            }

            var l_辺長 = p_引き継ぎ.A_k長 + 1;
            foreach (var l_位置 in l_位置群)
            {
                var l_終了 = Math.Min(l_位置, p_最小値列.Length - 1);
                for (var s = Math.Max(0, l_位置 + l_辺長 - p_k長); s <= l_終了; s++)
                {
                    p_最小値列[s] = 0;
                }
            }
        }

        /// <summary>
        /// 引き継ぎ配列のうち、この k の集合に無い k-mer を足す
        /// </summary>
        /// <param name="p_引き継ぎ配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 既にある k-mer は触らない (実際のリード由来の観測を優先する) <br/>
        /// 戻り値は足した k-mer の数
        /// </remarks>
        /// <returns></returns>
        public static int V_引き継ぎ(IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            using var l_計測 = new StageTimer($"carry-over k={p_k長}");
            return p_k長 <= TrustedKmerIndex.パック値のk上限
                ? Get_追加数_パック値(p_引き継ぎ配列, p_kmerインデックス, p_k長, p_リード長)
                : Get_追加数_塩基列(p_引き継ぎ配列, p_kmerインデックス, p_k長, p_リード長);
        }

        /// <summary>
        /// この k-mer に与えるカバレッジ
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_位置"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <returns></returns>
        public static ulong Get_引き継ぐカバレッジ(引き継ぎ配列 p_引き継ぎ, int p_位置, int p_k長, int? p_リード長)
        {
            var l_終端 = Math.Min(p_引き継ぎ.A_カバレッジ.Length - 1, p_位置 + p_k長 - p_引き継ぎ.A_k長);
            if (p_位置 > l_終端)
            {
                return 0UL;
            }

            if (p_引き継ぎ.A_分岐の継ぎ目位置 is { } l_継ぎ目 && l_継ぎ目.Any(x => p_位置 <= x && x + p_引き継ぎ.A_k長 + 1 <= p_位置 + p_k長))
            {
                return 0UL;
            }

            var l_最小 = int.MaxValue;
            for (var i = p_位置; i <= l_終端; i++)
            {
                l_最小 = Math.Min(l_最小, p_引き継ぎ.A_カバレッジ[i]);
            }
            return Get_換算カバレッジ(l_最小, p_引き継ぎ.A_k長, p_k長, p_リード長);
        }

        /// <summary>各窓の最小値を単調な待ち行列で求める</summary>
        /// <param name="p_値">元のカバレッジ</param>
        /// <param name="p_幅">窓の幅</param>
        /// <param name="p_結果">位置ごとの最小値</param>
        /// <param name="p_待ち行列">元のカバレッジ以上の長さの作業領域</param>
        internal static void V_計算_最小値列(int[] p_値, int p_幅, Span<int> p_結果, Span<int> p_待ち行列)
        {
            var l_先頭 = 0;
            var l_末尾 = 0;
            var l_読込位置 = 0;
            for (var i = 0; i < p_結果.Length; i++)
            {
                while (l_先頭 < l_末尾 && p_待ち行列[l_先頭] < i)
                {
                    l_先頭++;
                }
                var l_終端 = Math.Min(p_値.Length, i + p_幅);
                while (l_読込位置 < l_終端)
                {
                    while (l_先頭 < l_末尾 && p_値[p_待ち行列[l_末尾 - 1]] >= p_値[l_読込位置])
                    {
                        l_末尾--;
                    }
                    p_待ち行列[l_末尾++] = l_読込位置++;
                }
                p_結果[i] = p_幅 <= 0 || i >= p_値.Length || l_先頭 == l_末尾 ? 0 : p_値[p_待ち行列[l_先頭]];
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// パック値で窓を転がし、集合に無い k-mer を足す (k &lt;= 128)
        /// </summary>
        /// <param name="p_引き継ぎ配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 集合に有るかの判定は読むだけなので配列ごとに並列に行い、足すのは配列の順に 1 本のスレッドで行う<br/>
        /// 同じ k-mer を複数の配列が持つとき先に来た配列のカバレッジが残る点は、逐次に足した場合と変わらない
        /// </remarks>
        /// <returns>足した k-mer の数</returns>
        private static int Get_追加数_パック値(IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_バッチ長 = (int)Consts.進捗ログ間隔;
            var l_未登録群 = new List<(UInt128 A_上位, UInt128 A_下位, ulong A_カバレッジ)>?[Math.Min(l_バッチ長, p_引き継ぎ配列.Count)];
            var l_追加数 = 0;

            for (var l_開始 = 0; l_開始 < p_引き継ぎ配列.Count; l_開始 += l_バッチ長)
            {
                var l_件数 = Math.Min(l_バッチ長, p_引き継ぎ配列.Count - l_開始);
                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, () => (A_最小値列: Array.Empty<int>(), A_待ち行列: Array.Empty<int>()), (i, _, l_作業域) =>
                {
                    l_未登録群[i] = Get_未登録kmer(p_引き継ぎ配列[l_開始 + i], p_kmerインデックス, p_k長, p_リード長, ref l_作業域.A_最小値列, ref l_作業域.A_待ち行列);
                    return l_作業域;
                }, _ => { });

                for (var i = 0; i < l_件数; i++)
                {
                    if (l_未登録群[i] is not { } l_未登録)
                    {
                        continue;
                    }
                    foreach (var (l_上位, l_下位, l_カバレッジ) in l_未登録)
                    {
                        if (p_kmerインデックス.Try追加_信頼kmer_正規形(l_上位, l_下位, l_カバレッジ))
                        {
                            l_追加数++;
                        }
                    }
                    l_未登録群[i] = null;
                }

                // 合成リードを作ると引き継ぎ配列はリードペアの数まで増える
                // 全 k-mer を足し終えるまで無言だと止まったように見える
                if (l_件数 == l_バッチ長)
                {
                    Logger.V_出力(メッセージID.引き継ぎの統合進捗, l_開始 + l_件数, p_引き継ぎ配列.Count);
                }
            }
            return l_追加数;
        }

        /// <summary>
        /// 引き継ぎ配列 1 本のうち、集合に無く、与えるカバレッジが正の k-mer を配列の順に返す (k &lt;= 128)
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_最小値列">窓ごとの最小カバレッジの作業領域、足りなければ確保し直す</param>
        /// <param name="p_待ち行列">最小値を求める待ち行列の作業領域、足りなければ確保し直す</param>
        /// <returns>1 つも無ければ null</returns>
        private static List<(UInt128 A_上位, UInt128 A_下位, ulong A_カバレッジ)>? Get_未登録kmer(引き継ぎ配列 p_引き継ぎ, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長, ref int[] p_最小値列, ref int[] p_待ち行列)
        {
            if (p_引き継ぎ.A_配列.Length < p_k長)
            {
                return null;
            }

            var l_窓数 = p_引き継ぎ.A_配列.Length - p_k長 + 1;
            if (p_最小値列.Length < l_窓数)
            {
                p_最小値列 = new int[l_窓数];
            }
            if (p_待ち行列.Length < p_引き継ぎ.A_カバレッジ.Length)
            {
                p_待ち行列 = new int[p_引き継ぎ.A_カバレッジ.Length];
            }
            V_計算_最小値列(p_引き継ぎ.A_カバレッジ, p_k長 - p_引き継ぎ.A_k長 + 1, p_最小値列.AsSpan(0, l_窓数), p_待ち行列);
            V_除外_継ぎ目を含む窓(p_引き継ぎ, p_k長, p_最小値列.AsSpan(0, l_窓数));
            V_除外_未観測の範囲(p_引き継ぎ, p_k長, p_最小値列.AsSpan(0, l_窓数));

            List<(UInt128 A_上位, UInt128 A_下位, ulong A_カバレッジ)>? l_未登録 = null;
            var l_窓 = new RollingKmer(p_k長);
            for (var i = 0; i < p_引き継ぎ.A_配列.Length; i++)
            {
                if (!l_窓.Try追加(p_引き継ぎ.A_配列[i], out var l_キー) || p_kmerインデックス.Haskmer_正規形(l_キー.A_上位, l_キー.A_下位))
                {
                    continue;
                }
                var l_カバレッジ = Get_換算カバレッジ(p_最小値列[i - p_k長 + 1], p_引き継ぎ.A_k長, p_k長, p_リード長);
                if (l_カバレッジ > 0UL)
                {
                    (l_未登録 ??= []).Add((l_キー.A_上位, l_キー.A_下位, l_カバレッジ));
                }
            }
            return l_未登録;
        }

        /// <summary>
        /// 語の並びで窓を転がし、集合に無い k-mer を足す (k &gt; 128)
        /// </summary>
        /// <param name="p_引き継ぎ配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 足す順序は <see cref="Get_追加数_パック値"/> と同じく配列の順に 1 本のスレッドで保つ
        /// </remarks>
        /// <returns>足した k-mer の数</returns>
        private static int Get_追加数_塩基列(IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_バッチ長 = (int)Consts.進捗ログ間隔;
            var l_未登録群 = new List<(KmerKey A_キー, ulong A_カバレッジ)>?[Math.Min(l_バッチ長, p_引き継ぎ配列.Count)];
            var l_追加数 = 0;

            for (var l_開始 = 0; l_開始 < p_引き継ぎ配列.Count; l_開始 += l_バッチ長)
            {
                var l_件数 = Math.Min(l_バッチ長, p_引き継ぎ配列.Count - l_開始);
                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, () => (A_最小値列: Array.Empty<int>(), A_待ち行列: Array.Empty<int>()), (i, _, l_作業域) =>
                {
                    l_未登録群[i] = Get_未登録kmer_塩基列(p_引き継ぎ配列[l_開始 + i], p_kmerインデックス, p_k長, p_リード長, ref l_作業域.A_最小値列, ref l_作業域.A_待ち行列);
                    return l_作業域;
                }, _ => { });

                for (var i = 0; i < l_件数; i++)
                {
                    if (l_未登録群[i] is not { } l_未登録)
                    {
                        continue;
                    }
                    foreach (var (l_キー, l_カバレッジ) in l_未登録)
                    {
                        if (p_kmerインデックス.Try追加_信頼kmer_正規形(l_キー, l_カバレッジ))
                        {
                            l_追加数++;
                        }
                    }
                    l_未登録群[i] = null;
                }

                if (l_件数 == l_バッチ長)
                {
                    Logger.V_出力(メッセージID.引き継ぎの統合進捗, l_開始 + l_件数, p_引き継ぎ配列.Count);
                }
            }
            return l_追加数;
        }

        /// <summary>
        /// 引き継ぎ配列 1 本のうち、集合に無く、与えるカバレッジが正の k-mer を配列の順に返す (k &gt; 128)
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_最小値列">窓ごとの最小カバレッジの作業領域、足りなければ確保し直す</param>
        /// <param name="p_待ち行列">最小値を求める待ち行列の作業領域、足りなければ確保し直す</param>
        /// <returns>1 つも無ければ null</returns>
        private static List<(KmerKey A_キー, ulong A_カバレッジ)>? Get_未登録kmer_塩基列(引き継ぎ配列 p_引き継ぎ, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長, ref int[] p_最小値列, ref int[] p_待ち行列)
        {
            if (p_引き継ぎ.A_配列.Length < p_k長)
            {
                return null;
            }

            var l_窓数 = p_引き継ぎ.A_配列.Length - p_k長 + 1;
            if (p_最小値列.Length < l_窓数)
            {
                p_最小値列 = new int[l_窓数];
            }
            if (p_待ち行列.Length < p_引き継ぎ.A_カバレッジ.Length)
            {
                p_待ち行列 = new int[p_引き継ぎ.A_カバレッジ.Length];
            }
            V_計算_最小値列(p_引き継ぎ.A_カバレッジ, p_k長 - p_引き継ぎ.A_k長 + 1, p_最小値列.AsSpan(0, l_窓数), p_待ち行列);
            V_除外_継ぎ目を含む窓(p_引き継ぎ, p_k長, p_最小値列.AsSpan(0, l_窓数));
            V_除外_未観測の範囲(p_引き継ぎ, p_k長, p_最小値列.AsSpan(0, l_窓数));

            List<(KmerKey A_キー, ulong A_カバレッジ)>? l_未登録 = null;
            var l_窓 = new WideRollingKmer(p_k長);
            for (var i = 0; i < p_引き継ぎ.A_配列.Length; i++)
            {
                if (!l_窓.Try追加(p_引き継ぎ.A_配列[i], out var l_キー))
                {
                    continue;
                }

                var l_カバレッジ = Get_換算カバレッジ(p_最小値列[i - p_k長 + 1], p_引き継ぎ.A_k長, p_k長, p_リード長);
                if (l_カバレッジ > 0UL && !p_kmerインデックス.Haskmer_正規形(l_キー))
                {
                    (l_未登録 ??= []).Add((l_キー.Get_複製(), l_カバレッジ));
                }
            }
            return l_未登録;
        }

        /// <summary>最小カバレッジを次の k の観測本数へ換算する</summary>
        /// <param name="p_最小">最小カバレッジ</param>
        /// <param name="p_前k長">前段の k</param>
        /// <param name="p_k長">次段の k</param>
        /// <param name="p_リード長">リード長</param>
        /// <returns>換算したカバレッジ</returns>
        private static ulong Get_換算カバレッジ(int p_最小, int p_前k長, int p_k長, int? p_リード長)
        {
            if (p_最小 <= 0)
            {
                return 0UL;
            }

            if (p_リード長 is not { } l_リード長 || l_リード長 <= p_k長)
            {
                return (ulong)p_最小;
            }

            var l_前段の本数 = l_リード長 - p_前k長 + 1;
            var l_今の本数 = l_リード長 - p_k長 + 1;
            return l_前段の本数 <= 0 || l_今の本数 <= 0 ? (ulong)p_最小 : (ulong)Math.Max(1L, (long)Math.Round((double)p_最小 * l_今の本数 / l_前段の本数));
        }

        #endregion
    }
}
