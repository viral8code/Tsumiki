using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Model.Preprocessing;
using Tsumiki.Model.Scaffolding;
using Tsumiki.Utility;

namespace Tsumiki.Core.Preprocessing
{
    /// <summary>
    /// ペアエンドの 2 本を、間の未読区間ごと 1 本の合成リード (SuperRead) に統合する
    /// </summary>
    /// <remarks>
    /// MaSuRCA の SuperReads に相当する<br/>
    /// 重なりが一つに定まるペアはそれで繋ぎ、繋げないペアは read1 の末尾 k-mer から
    /// RC(read2) の先頭 k-mer まで信頼できる k-mer 集合の中で経路を探す<br/>
    /// 使う探索エンジンは GapFiller (スキャフォールドのギャップ埋め) と共通の
    /// ConstrainedPathFinder で、対象をスキャフォールドのギャップからリードペアに変えただけになる<br/>
    /// Tsumiki のマルチ k は自動選択の上限がリード長に縛られている (<see cref="Consts.マルチk上限のリード長比"/>)<br/>
    /// 合成リードは元のリードより長いため、KmerCarryOver 経由で次の k へ渡せばこの上限を実効的に外せる<br/>
    /// 統合に失敗したペアは単に対象から外れるだけで、元のリードは通常どおり k-mer カウントに使われ続ける<br/>
    /// 失敗が増えても元のペアのままに退化するだけで、悪化はしない
    /// </remarks>
    internal static class SuperReadJoiner
    {
        /// <summary>
        /// 橋渡しする長さの上限
        /// </summary>
        /// <remarks>
        /// 細菌ゲノムの一般的なライブラリではフラグメント長は高々 1 kb 程度に収まるため、
        /// これを大きく超える探索は時間をかけても一意に定まる見込みが薄い
        /// (GapFiller のギャップ長上限と同じ考え方)
        /// </remarks>
        private const int 橋渡し長の上限 = 500;

        /// <summary>
        /// 1 バッチあたりのペア数
        /// </summary>
        private const int バッチサイズ = 5000;

        /// <summary>
        /// 1 ペアあたりに展開してよい探索状態の上限
        /// </summary>
        /// <remarks>
        /// ギャップ充填は数千箇所だが橋渡しはリードペアの数だけ走るため、
        /// 1 件あたりの上限を絞らないと解けない少数のペアに全体の時間を持っていかれる
        /// </remarks>
        private const int 橋渡しの状態数上限 = 20_000;

        /// <summary>
        /// -i でインサートサイズが分かっているときに、そこから見積もった橋渡し長に掛ける許容比
        /// </summary>
        /// <remarks>
        /// 探索の深さは実行時間を直接決めるので、分かっている手掛かりで絞る<br/>
        /// 実際より短く見積もったペアは橋渡しに失敗するだけで、元のリードとしては通常どおり残る
        /// </remarks>
        private const double インサートサイズの許容比 = 1.5D;

        /// <summary>
        /// ペアの FASTQ を読み込み、統合できたペアを合成配列 (引き継ぎ配列) として返す
        /// </summary>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_統計">統合の内訳</param>
        /// <returns>統合できたペアの合成配列</returns>
        public static List<引き継ぎ配列> Get_合成リード(
            string p_リード1のパス, string p_リード2のパス,
            TrustedKmerIndex p_kmerインデックス, int p_k長, out SuperRead統計 p_統計)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_インサートサイズ = ConfigurationManager.A_実行時引数.A_インサートサイズ;

            var l_総ペア数 = 0;
            var l_統合数 = 0;
            var l_重なり結合数 = 0;
            var l_曖昧で捨てた数 = 0;
            ulong l_出力済みの区切り = 0UL;
            List<引き継ぎ配列> l_結果 = [];

            Logger.V_出力(メッセージID.SuperRead橋渡し開始);

            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);

            var l_配列1群 = new string[バッチサイズ];
            var l_配列2群 = new string[バッチサイズ];
            var l_統合結果群 = new (string? A_配列, bool A_重なりで結合したか, int A_曖昧で捨てた数)[バッチサイズ];

            while (l_読み込み1.Get_続きがあるか() && l_読み込み2.Get_続きがあるか())
            {
                var l_件数 = 0;
                while (l_件数 < バッチサイズ
                    && l_読み込み1.Get_続きがあるか() && l_読み込み2.Get_続きがあるか())
                {
                    l_配列1群[l_件数] = l_読み込み1.Get_次のリード_軽量().A_生リード;
                    l_配列2群[l_件数] = l_読み込み2.Get_次のリード_軽量().A_生リード;
                    l_件数++;
                }
                l_総ペア数 += l_件数;

                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                {
                    l_統合結果群[i] = Get_合成配列_内訳つき(
                        l_配列1群[i], l_配列2群[i], p_kmerインデックス, p_k長, l_インサートサイズ);
                });

                for (var i = 0; i < l_件数; i++)
                {
                    l_曖昧で捨てた数 += l_統合結果群[i].A_曖昧で捨てた数;
                    if (l_統合結果群[i].A_配列 is not { } l_配列)
                    {
                        continue;
                    }
                    l_統合数++;
                    if (l_統合結果群[i].A_重なりで結合したか)
                    {
                        l_重なり結合数++;
                    }
                    l_結果.Add(Get_引き継ぎ配列(l_配列, p_kmerインデックス, p_k長));
                }

                // 全ペアを走査するうえ 1 件ごとの探索も重く、無言のまま数十分経つことがあるため、進み具合が分かるようにする
                var l_区切り = (ulong)l_総ペア数 / Consts.進捗ログ間隔;
                if (l_区切り > l_出力済みの区切り)
                {
                    l_出力済みの区切り = l_区切り;
                    Logger.V_出力(メッセージID.SuperRead橋渡し進捗, l_総ペア数, l_統合数);
                }
            }

            p_統計 = new SuperRead統計(l_総ペア数, l_統合数, l_重なり結合数, l_曖昧で捨てた数);
            return l_結果;
        }

        /// <summary>
        /// 1 ペア分の統合を試みる
        /// </summary>
        /// <remarks>
        /// 副作用のない純粋関数
        /// </remarks>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>合成配列、統合できなかった場合は null</returns>
        internal static string? Get_合成配列(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_インサートサイズ = null)
        {
            return Get_合成配列_内訳つき(p_配列1, p_配列2, p_kmerインデックス, p_k長, p_インサートサイズ).A_配列;
        }

        /// <summary>
        /// 1 ペア分の統合を、どちらの手段で繋いだかと一緒に返す
        /// </summary>
        /// <remarks>
        /// 先にペアの重なりを試すのは、断片長がリード長の 2 倍を下回るライブラリでは
        /// read1 と RC(read2) が重なり、橋渡しに必要な長さが負になるため<br/>
        /// この場合グラフ探索は前向きにしか進めないので構造的に解けない
        /// </remarks>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>合成配列と、重なりで繋いだかどうか、重なりが曖昧で捨てた数</returns>
        internal static (string? A_配列, bool A_重なりで結合したか, int A_曖昧で捨てた数) Get_合成配列_内訳つき(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_インサートサイズ)
        {
            var l_曖昧 = 0;
            if (Get_重なりで結合(p_配列1, p_配列2, p_kmerインデックス, p_k長, ref l_曖昧) is { } l_重なり結合)
            {
                return (l_重なり結合, true, 0);
            }
            return (Get_橋渡しで結合(p_配列1, p_配列2, p_kmerインデックス, p_k長, p_インサートサイズ), false, l_曖昧);
        }

        /// <summary>
        /// read1 と RC(read2) の重なりから断片を復元する
        /// </summary>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_曖昧で捨てた数">重なりが一つに定まらず捨てた数、該当すれば 1 加算する</param>
        /// <returns>
        /// 復元した断片<br/>
        /// 重なりが見つからない (断片がリード長の 2 倍を超える) 場合と、繋いでも伸びない場合は null
        /// </returns>
        private static string? Get_重なりで結合(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, ref int p_曖昧で捨てた数)
        {
            var l_RC配列2 = Util.V_逆相補_曖昧塩基あり(p_配列2);
            var l_重なり = Preprocessor.Get_最適オーバーラップ(
                Util.V_変換_塩基列(p_配列1), Util.V_変換_塩基列(l_RC配列2),
                Consts.ペア結合の最小重なり長, Consts.ペア結合の許容不一致率,
                out var l_対抗馬があるか);
            if (l_重なり is not { } l_位置合わせ || l_位置合わせ.A_offset < 0)
            {
                return null;
            }

            // 反復配列の中では周期のぶんだけずれた位置も同じくらい良く合い、
            // 最良を 1 つ選べてしまうため、対抗馬の有無を見ないと別コピーを掴んだことに気づけない
            if (l_対抗馬があるか)
            {
                p_曖昧で捨てた数++;
                return null;
            }

            // 断片が read1 に収まっているならアダプタ読み抜けで前処理のトリムの領分になり、ここで作っても長さが伸びない
            var l_フラグメント長 = l_位置合わせ.A_offset + p_配列2.Length;
            if (l_フラグメント長 <= p_配列1.Length)
            {
                return null;
            }

            var l_合成 = p_配列1 + l_RC配列2[(p_配列1.Length - l_位置合わせ.A_offset)..];
            return l_合成.Length >= p_k長 && Get_継ぎ目が支持されているか(l_合成, p_配列1.Length, p_kmerインデックス, p_k長)
                ? l_合成
                : null;
        }

        /// <summary>
        /// 合成配列の継ぎ目 (read1 から RC(read2) へ切り替わる位置) を跨ぐ k-mer が、
        /// 信頼できる k-mer 集合に入っているか
        /// </summary>
        /// <remarks>
        /// 重なりに許容範囲内の不一致が残っていると、繋いだ配列にはどちらのリードとも違う塩基が入り、
        /// 継ぎ目を跨ぐ k-mer だけがどのリードにも存在しないものになる<br/>
        /// 存在しない配列をグラフへ持ち込まないための関門になる<br/>
        /// 重なりが完全一致なら継ぎ目の k-mer は read2 側の k-mer と一致するため、
        /// 反復配列の別コピーを掴む取り違えはここでは防げず、それは重なりの一意性で抑える
        /// </remarks>
        /// <param name="p_合成">継ぎ目を含む合成配列</param>
        /// <param name="p_継ぎ目">read1 から RC(read2) へ切り替わる位置</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <returns>継ぎ目を跨ぐ k-mer がすべて集合にあれば true</returns>
        private static bool Get_継ぎ目が支持されているか(string p_合成, int p_継ぎ目, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_開始 = Math.Max(0, p_継ぎ目 - p_k長 + 1);
            var l_終了 = Math.Min(p_継ぎ目, p_合成.Length - p_k長);
            for (var l_位置 = l_開始; l_位置 <= l_終了; l_位置++)
            {
                var l_kmer = Util.V_変換_塩基列(p_合成.Substring(l_位置, p_k長));
                if (Array.IndexOf(l_kmer, Consts.無効な塩基) >= 0
                    || !p_kmerインデックス.Get_含まれるか(l_kmer))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// read1 の末尾 k-mer から RC(read2) の先頭 k-mer まで、信頼できる k-mer 集合の中で経路を探して繋ぐ
        /// </summary>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>
        /// 繋いだ配列<br/>
        /// どちらかの端の k-mer が集合に無い、あるいは経路が一意に定まらない場合は null
        /// </returns>
        private static string? Get_橋渡しで結合(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_インサートサイズ)
        {
            if (p_配列1.Length < p_k長 || p_配列2.Length < p_k長)
            {
                return null;
            }

            var l_最大長 = Get_橋渡し長の上限(p_配列1.Length, p_配列2.Length, p_インサートサイズ);
            if (l_最大長 < 0)
            {
                return null;
            }

            // 端の k-mer が集合に無いペアが大半を占めるので、全長の変換と RC を先に作ると、その大半で捨てるだけの配列を確保することになる
            var l_左のkmer = Util.V_変換_塩基列(p_配列1[^p_k長..]);
            if (Array.IndexOf(l_左のkmer, Consts.無効な塩基) >= 0
                || !p_kmerインデックス.Get_含まれるか(l_左のkmer))
            {
                return null;
            }

            // RC(read2) の先頭 k-mer は、read2 の末尾 k 塩基の逆相補と一致する
            var l_目標kmer = Util.V_変換_塩基列(Util.V_逆相補_曖昧塩基あり(p_配列2[^p_k長..]));
            if (Array.IndexOf(l_目標kmer, Consts.無効な塩基) >= 0
                || !p_kmerインデックス.Get_含まれるか(l_目標kmer))
            {
                return null;
            }

            var (l_橋渡し配列, l_判定) = ConstrainedPathFinder.Get_経路(
                l_左のkmer, l_目標kmer, p_最小長: 0, p_最大長: l_最大長, p_kmerインデックス, p_k長,
                橋渡しの状態数上限);
            return l_判定 == ギャップ充填判定.充填済み
                ? p_配列1 + l_橋渡し配列 + Util.V_逆相補_曖昧塩基あり(p_配列2)
                : null;
        }

        /// <summary>
        /// このペアで探索してよい橋渡し長の上限
        /// </summary>
        /// <remarks>
        /// インサートサイズが分かっていれば、そこから見積もった長さまでに絞る
        /// </remarks>
        /// <param name="p_長さ1">read1 の長さ</param>
        /// <param name="p_長さ2">read2 の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>探索してよい橋渡し長の上限</returns>
        private static int Get_橋渡し長の上限(int p_長さ1, int p_長さ2, int? p_インサートサイズ)
        {
            if (p_インサートサイズ is not { } l_インサートサイズ)
            {
                return 橋渡し長の上限;
            }
            var l_見積もり = (int)(l_インサートサイズ * インサートサイズの許容比) - p_長さ1 - p_長さ2;
            return Math.Min(橋渡し長の上限, l_見積もり);
        }

        /// <summary>
        /// 合成配列を、位置ごとのカバレッジ付きの引き継ぎ配列にする
        /// </summary>
        /// <remarks>
        /// KmerCarryOver.Get_引き継ぎ配列 と同じ計算になる
        /// </remarks>
        /// <param name="p_配列">合成配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <returns>次の k へ渡す引き継ぎ配列</returns>
        private static 引き継ぎ配列 Get_引き継ぎ配列(string p_配列, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_塩基列 = p_配列.Select(Util.Get_塩基ID).ToArray();
            var l_カバレッジ = new int[p_配列.Length - p_k長 + 1];
            for (var i = 0; i < l_カバレッジ.Length; i++)
            {
                l_カバレッジ[i] = (int)Math.Min(int.MaxValue, p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長)));
            }
            return new 引き継ぎ配列(p_配列, l_カバレッジ, p_k長);
        }

        /// <summary>
        /// 統合の内訳をログへ出力する
        /// </summary>
        /// <param name="p_統計">統合の内訳</param>
        public static void V_出力_統計(SuperRead統計 p_統計)
        {
            Logger.V_出力(メッセージID.SuperRead統計, p_統計.A_統合数, p_統計.A_総ペア数);
            Logger.V_出力(メッセージID.SuperRead統計_重なり, p_統計.A_重なり結合数, p_統計.A_橋渡し数);
            if (p_統計.A_曖昧で捨てた数 > 0)
            {
                Logger.V_出力(メッセージID.SuperRead統計_曖昧, p_統計.A_曖昧で捨てた数);
            }
        }
    }
}
