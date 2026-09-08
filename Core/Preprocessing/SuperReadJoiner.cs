using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// ペアエンドの2本を、間の未読区間ごと1本の「合成リード(SuperRead)」に
    /// 統合する(MaSuRCA の SuperReads に相当)。read1 の末尾 k-mer から
    /// RC(read2) の先頭 k-mer まで、信頼できる k-mer 集合の中で経路が
    /// ちょうど1本に定まったペアだけを統合する。使う探索エンジンは
    /// GapFiller(スキャフォールドのギャップ埋め)と共通の ConstrainedPathFinder で、
    /// 対象をスキャフォールドのギャップからリードペアに変えただけになる。
    ///
    /// Tsumiki のマルチ k は自動選択の上限がリード長に縛られている
    /// (Consts.マルチk上限のリード長比)。合成リードは元のリードの2〜4倍の
    /// 長さを持つため、KmerCarryOver 経由で次の k へ渡せばこの上限を実効的に外せる。
    ///
    /// 統合に失敗したペアは単に対象から外れるだけで、元のリードは通常どおり
    /// k-mer カウントに使われ続ける。失敗が増えても「元のペアのまま」に
    /// 退化するだけで、悪化はしない。
    /// </summary>
    internal static class SuperReadJoiner
    {
        /// <summary>
        /// 橋渡しする長さの上限。細菌ゲノムの一般的なライブラリでは
        /// フラグメント長は高々1kb程度に収まるため、これを大きく超える
        /// 探索は時間をかけても一意に定まる見込みが薄い
        /// (GapFiller のギャップ長上限と同じ考え方)。
        /// </summary>
        private const int 橋渡し長の上限 = 500;

        private const int バッチサイズ = 5000;

        /// <summary>
        /// 1ペアあたりに展開してよい探索状態の上限。ギャップ充填は数千箇所だが
        /// 橋渡しはリードペアの数だけ走るため、1件あたりの上限を絞らないと
        /// 解けない少数のペアに全体の時間を持っていかれる。
        /// </summary>
        private const int 橋渡しの状態数上限 = 20_000;

        /// <summary>
        /// -i でインサートサイズが分かっているときに、そこから見積もった
        /// 橋渡し長に掛ける許容比。探索の深さは実行時間を直接決めるので、
        /// 分かっている手掛かりで絞る。実際より短く見積もったペアは
        /// 橋渡しに失敗するだけで、元のリードとしては通常どおり残る。
        /// </summary>
        private const double インサートサイズの許容比 = 1.5;

        /// <summary>
        /// ペアの FASTQ を読み込み、統合できたペアを合成配列(引き継ぎ配列)として返す。
        /// </summary>
        public static List<引き継ぎ配列> Get_合成リード(
            string p_リード1のパス, string p_リード2のパス,
            TrustedKmerIndex p_kmerインデックス, int p_k長, out SuperRead統計 p_統計)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_インサートサイズ = ConfigurationManager.A_実行時引数.A_インサートサイズ;

            var l_総ペア数 = 0;
            var l_統合数 = 0;
            ulong l_出力済みの区切り = 0;
            List<引き継ぎ配列> l_結果 = [];

            Logger.V_出力(メッセージID.SuperRead橋渡し開始);

            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);

            var l_配列1群 = new string[バッチサイズ];
            var l_配列2群 = new string[バッチサイズ];
            var l_統合結果群 = new string?[バッチサイズ];

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
                    l_統合結果群[i] = Get_合成配列(
                        l_配列1群[i], l_配列2群[i], p_kmerインデックス, p_k長, l_インサートサイズ);
                });

                for (var i = 0; i < l_件数; i++)
                {
                    if (l_統合結果群[i] is not { } l_配列)
                    {
                        continue;
                    }
                    l_統合数++;
                    l_結果.Add(Get_引き継ぎ配列(l_配列, p_kmerインデックス, p_k長));
                }

                // 全ペアを走査するうえ1件ごとの探索も重い。無言のまま数十分
                // 経つことがあるため、進み具合が分かるようにする。
                var l_区切り = (ulong)l_総ペア数 / Consts.進捗ログ間隔;
                if (l_区切り > l_出力済みの区切り)
                {
                    l_出力済みの区切り = l_区切り;
                    Logger.V_出力(メッセージID.SuperRead橋渡し進捗, l_総ペア数, l_統合数);
                }
            }

            p_統計 = new SuperRead統計(l_総ペア数, l_統合数);
            return l_結果;
        }

        /// <summary>
        /// 1ペア分の統合を試みる。副作用のない純粋関数。
        /// read1 の末尾 k-mer と RC(read2) の先頭 k-mer のどちらかが信頼できる
        /// k-mer 集合に無い、あるいは繋ぐ経路が一意に定まらない場合は null。
        /// </summary>
        internal static string? Get_合成配列(
            string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長,
            int? p_インサートサイズ = null)
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

            // 端の k-mer が集合に無いペアが大半を占める。全長の変換と RC を
            // 先に作ると、その大半で捨てるだけの配列を確保することになる。
            var l_左のkmer = Util.V_変換_塩基列(p_配列1[^p_k長..]);
            if (Array.IndexOf(l_左のkmer, Consts.無効な塩基) >= 0
                || !p_kmerインデックス.Get_含まれるか(l_左のkmer))
            {
                return null;
            }

            // RC(read2) の先頭 k-mer は、read2 の末尾 k 塩基の逆相補と一致する。
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
        /// このペアで探索してよい橋渡し長の上限。インサートサイズが分かって
        /// いれば、そこから見積もった長さまでに絞る。
        /// </summary>
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
        /// (KmerCarryOver.Get_引き継ぎ配列 と同じ計算)。
        /// </summary>
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

        public static void V_出力_統計(SuperRead統計 p_統計)
        {
            Logger.V_出力(メッセージID.SuperRead統計, p_統計.A_統合数, p_統計.A_総ペア数);
        }
    }
}
