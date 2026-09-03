using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-merの出現回数カウント(CountingDBによる厳密な外部マージソート)と、
    /// カットオフを通過した「信頼できるk-mer」の厳密な集合を保持するクラス。
    ///
    /// 以前はカットオフ後の膜(membership)判定に多重ハッシュのビット配列
    /// (Bloom filter)を使っていたが、複数ハッシュのうち1つ(shift=1)が
    /// 事実上塩基の並び順に依存しない合計値に退化しておりハッシュの
    /// 独立性が低いこと、そもそも近似判定である以上フォールスポジティブに
    /// よるグラフ構造の誤判定(誤った分岐点検出・誤った隣接判定)を
    /// 原理的に排除できないことが問題だった。
    ///
    /// 7Mbp程度のバクテリアゲノムであれば、カットオフを通過した信頼できる
    /// k-merの総数は現実的にせいぜい数千万件程度に収まり、厳密な
    /// HashSet&lt;KmerKey&gt; としてメモリに保持できる規模である。そのため
    /// カットオフ後は厳密な集合に置き換え、近似判定を完全に排除した。
    /// </summary>
    internal class TrustedKmerIndex : IDisposable
    {
        private readonly string _tempDirectory;

        // 並列読み込み用に、ワーカースレッドの数だけ CountingDB を用意する。
        // 各スレッドは自分専用のインスタンスにのみ書き込むため、
        // Add 呼び出し時のロックが不要になる。
        private CountingDB[]? _counters;

        // シャードごとのロック。k-mer をワーカー単位ではなくハッシュ値で
        // 振り分けるようにしたため、複数スレッドが同じシャードへ書きうる。
        private object[]? _counterLocks;

        // Cutoff() 実行後に確定する、カットオフを通過したk-merの厳密な集合
        // (常に正規化(Canonical)された形で保持する)。値はそのk-merの
        // 出現回数(カバレッジ)。GraphSimplifierの低カバレッジunitig除去
        // (bubble/tip相当のアーティファクトを、厳密なトポロジー判定なしに
        // カバレッジの相対的な低さで検出する)に使う。以降のグラフ探索
        // (IsFirstKmer/CountInEdges/UnitigMakerの伸長判定)は、値を見ない
        // 所属判定(ContainsKey相当)のみで行う。
        //
        // k<=32(実用上ほぼ常にこちら。デフォルトk=31)の場合は、KmerKeyの
        // ulong[]割り当てを経由せず、2bitパックしたulong1個で直接
        // Dictionary<ulong,ulong>を引く高速経路(_trustedKmersSmall)を使う。
        // ErrorCorrector は1リードあたり数百〜数千回 Contains を呼ぶため、
        // KmerKey経由(構築のたびにulong[]・byte[]を複数回ヒープ確保する)の
        // オーバーヘッドが無視できず、実データ規模で致命的に遅くなることが
        // 実測で判明したため導入した。k>32の場合のみ、従来通り厳密だが
        // 低速な Dictionary&lt;KmerKey,ulong&gt; にフォールバックする。
        private Dictionary<KmerKey, ulong>? _trustedKmers;

        private Dictionary<ulong, ulong>? _trustedKmersSmall;

        // 33 <= k <= 64 用。150bp リードで k=31 のままだと 31bp 以上の反復配列が
        // すべて潰れてしまい contig N50 が伸びないため、k を 63 前後まで上げられる
        // ことが品質上きわめて重要になる。UInt128 は値型なので、KmerKey
        // (ulong[] を毎回ヒープ確保する)と違い割り当てが一切発生せず、
        // k<=32 の ulong 経路とほぼ同等の速度で動く。
        private Dictionary<UInt128, ulong>? _trustedKmersMid;

        private static bool UseSmallPath => ConfigurationManager.Arguments.Kmer <= 32;

        private static bool UseMidPath => ConfigurationManager.Arguments.Kmer is > 32 and <= 64;

        public TrustedKmerIndex(string tempDirectory)
        {
            this._tempDirectory = tempDirectory;
            var workerCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);
            this._counters = new CountingDB[workerCount];
            this._counterLocks = new object[workerCount];
            for (var i = 0; i < workerCount; i++)
            {
                this._counters[i] = new CountingDB(tempDirectory, workerCount);
                this._counterLocks[i] = new object();
            }
        }

        /// <summary>
        /// 指定したワーカー番号(0 始まり)専用の CountingDB に登録する。
        /// 並列読み込み時、各スレッドは自分の workerIndex を固定して呼び出すことで
        /// スレッドセーフに(ロックなしで)k-mer を登録できる。
        /// </summary>
        /// <summary>
        /// 曖昧塩基(N など、候補が複数ある位置)を含む k-mer を、
        /// ありうる塩基の組み合わせすべてに展開して登録する。
        /// read[i] はその位置で取りうる塩基ID の一覧。
        ///
        /// 展開は塩基ID の空間で行い、1件ずつ通常の Add へ渡す。
        /// こうすることで正規化(順鎖・逆鎖の寄せ)とシャード振り分けが
        /// 通常経路とまったく同じ扱いになる。
        ///
        /// 以前は CountingDB 側でパック済みバイト列を組み立てており、
        /// (a) 正規化されない (b) 組み立て途中のバッファをそのまま辞書のキーと
        /// して格納しており、その後も書き換え続けるため既に登録済みのキーが
        /// 壊れる、という2つの問題があった。
        /// </summary>
        public void Add(Span<byte[]> read, int workerIndex)
        {
            if (this._counters is null)
            {
                return;
            }
            var kmer = new byte[read.Length];
            this.AddCombinations(read, 0, kmer, workerIndex);
        }

        private void AddCombinations(Span<byte[]> options, int position, byte[] kmer, int workerIndex)
        {
            if (position == options.Length)
            {
                this.Add(kmer.AsSpan(), workerIndex);
                return;
            }
            foreach (var id in options[position])
            {
                kmer[position] = id;
                this.AddCombinations(options, position + 1, kmer, workerIndex);
            }
        }

        /// <summary>
        /// k-mer を1件カウントする。
        ///
        /// 振り分けは workerIndex ではなく k-mer 自身のハッシュ値で行う。
        /// ワーカーごとに別の辞書へ入れていた頃は、同じ k-mer が最大で
        /// スレッド数ぶんの辞書に重複して載り、メモリも書き出し量も
        /// そのぶん膨らんでいた(実データ 100x・16スレッドでピーク12.5GB)。
        /// ハッシュで振り分ければ、ある k-mer は必ず1つのシャードにしか載らない。
        ///
        /// さらに、順鎖・逆鎖のうち辞書順で小さいほう(正規化形)に寄せてから
        /// 数える。以前は両向きを別キーとして数え、カットオフ時に合算していたため、
        /// エントリ数・書き出し量ともに2倍になっていた。
        /// </summary>
        public void Add(Span<byte> read, int workerIndex)
        {
            if (this._counters is not { } counters)
            {
                return;
            }

            var packed = PackCanonical(read);
            var shard = (int)(Hash(packed) % (uint)counters.Length);
            lock (this._counterLocks![shard])
            {
                counters[shard].AddPacked(packed);
            }
        }

        /// <summary>
        /// k-mer を、順鎖・逆鎖のうち塩基列として辞書順で小さいほうの向きで
        /// 2bit パックした byte 配列にする。
        ///
        /// パック後のバイト列の辞書順は塩基列の辞書順と一致する(先頭塩基が
        /// 上位ビットに来るため)ので、外部マージソートの順序とも整合する。
        /// </summary>
        private static byte[] PackCanonical(ReadOnlySpan<byte> kmer)
        {
            var useForward = IsCanonicalForward(kmer);
            var arr = new byte[(kmer.Length + 3) / 4];
            for (var i = 0; i < kmer.Length; i++)
            {
                // 逆鎖側を採用する場合は、末尾から相補塩基を取り出す。
                // 相補は A(1)<->T(4), C(2)<->G(3) なので 5 - x で得られる。
                var id = useForward ? kmer[i] : (byte)(5 - kmer[kmer.Length - 1 - i]);
                arr[i >> 2] |= (byte)((id - 1) << ((3 - (i & 3)) << 1));
            }
            return arr;
        }

        /// <summary>
        /// 順鎖側がその逆相補以下(辞書順)かどうか。確保なしで判定する。
        /// </summary>
        private static bool IsCanonicalForward(ReadOnlySpan<byte> kmer)
        {
            int i = 0, j = kmer.Length - 1;
            while (i <= j)
            {
                var forward = kmer[i];
                var reverse = (byte)(5 - kmer[j]);
                if (forward != reverse)
                {
                    return forward < reverse;
                }
                i++;
                j--;
            }
            // 回文(自身が逆相補と一致)。どちらでも同じなので順鎖扱い。
            return true;
        }

        /// <summary>パック済みキーの FNV-1a ハッシュ。シャードの振り分けに使う。</summary>
        private static uint Hash(byte[] packed)
        {
            var hash = 2166136261u;
            foreach (var b in packed)
            {
                hash ^= b;
                hash *= 16777619u;
            }
            return hash;
        }

        /// <summary>
        /// kmer(順鎖・逆鎖いずれの向きでもよい)がカットオフを通過した
        /// 信頼できるk-mer集合に含まれるかどうかを厳密に判定する。
        /// </summary>
        public bool Contains(Span<byte> kmer)
        {
            if (UseSmallPath)
            {
                return this._trustedKmersSmall!.ContainsKey(CanonicalSmall(kmer));
            }
            if (UseMidPath)
            {
                return this._trustedKmersMid!.ContainsKey(CanonicalMid(kmer));
            }
            return this._trustedKmers!.ContainsKey(new KmerKey(kmer).Canonical());
        }

        /// <summary>
        /// kmerの出現回数(カバレッジ)を返す。信頼できるk-mer集合に
        /// 含まれない場合は0を返す。
        /// </summary>
        public ulong GetCoverage(Span<byte> kmer)
        {
            if (UseSmallPath)
            {
                return this._trustedKmersSmall!.GetValueOrDefault(CanonicalSmall(kmer), 0UL);
            }
            if (UseMidPath)
            {
                return this._trustedKmersMid!.GetValueOrDefault(CanonicalMid(kmer), 0UL);
            }
            return this._trustedKmers!.GetValueOrDefault(new KmerKey(kmer).Canonical(), 0UL);
        }

        /// <summary>
        /// kmer(塩基ID 1-4、長さ32以下)を2bit/塩基でulong1個にパックする。
        /// kmer[0]が最上位側、kmer[^1]が最下位側に来る(空きビットは下位側に残る)。
        /// </summary>
        private static ulong PackSmall(ReadOnlySpan<byte> kmer)
        {
            var value = 0UL;
            foreach (var b in kmer)
            {
                value = (value << 2) | ((ulong)b - 1);
            }
            return value;
        }

        /// <summary>
        /// PackSmallでパックした値の逆相補を、ヒープ確保なしで直接計算する。
        /// 2bitコドンごとに comp = codon ^ 0b11 (A&lt;-&gt;T, C&lt;-&gt;G)で複製し、
        /// 下位から順に取り出しつつ上位へ積み直すことでコドン順序も反転させる。
        /// </summary>
        private static ulong ReverseComplementSmall(ulong packed, int length)
        {
            var temp = packed;
            var result = 0UL;
            for (var i = 0; i < length; i++)
            {
                var codon = temp & 0x3UL;
                result = (result << 2) | (codon ^ 0x3UL);
                temp >>= 2;
            }
            return result;
        }

        private static ulong CanonicalSmall(ReadOnlySpan<byte> kmer)
        {
            var packed = PackSmall(kmer);
            var rev = ReverseComplementSmall(packed, kmer.Length);
            return Math.Min(packed, rev);
        }

        /// <summary>PackSmallの逆変換。kmer[^1]が最下位ビット側にあるため、末尾から復元する。</summary>
        private static byte[] UnpackSmall(ulong packed, int length)
        {
            var bytes = new byte[length];
            for (var i = length - 1; i >= 0; i--)
            {
                bytes[i] = (byte)((packed & 0x3UL) + 1);
                packed >>= 2;
            }
            return bytes;
        }

        /// <summary>
        /// PackSmall の 128bit 版(k は 64 以下)。ビット配置の規約は PackSmall と
        /// 同じで、kmer の先頭塩基が最上位側、末尾塩基が最下位側に来る。
        /// </summary>
        private static UInt128 PackMid(ReadOnlySpan<byte> kmer)
        {
            UInt128 value = 0;
            foreach (var b in kmer)
            {
                value = (value << 2) | (UInt128)(b - 1);
            }
            return value;
        }

        /// <summary>ReverseComplementSmall の 128bit 版。</summary>
        private static UInt128 ReverseComplementMid(UInt128 packed, int length)
        {
            var temp = packed;
            UInt128 result = 0;
            for (var i = 0; i < length; i++)
            {
                var codon = temp & 3;
                result = (result << 2) | (codon ^ 3);
                temp >>= 2;
            }
            return result;
        }

        private static UInt128 CanonicalMid(ReadOnlySpan<byte> kmer)
        {
            var packed = PackMid(kmer);
            var rev = ReverseComplementMid(packed, kmer.Length);
            return packed < rev ? packed : rev;
        }

        /// <summary>PackMid の逆変換。</summary>
        private static byte[] UnpackMid(UInt128 packed, int length)
        {
            var bytes = new byte[length];
            for (var i = length - 1; i >= 0; i--)
            {
                bytes[i] = (byte)((ulong)(packed & 3) + 1);
                packed >>= 2;
            }
            return bytes;
        }

        /// <summary>
        /// カットオフを通過した信頼できるk-merを(正規化された、いずれかの向きの)
        /// byte配列として1件ずつ列挙する。GraphSimplifier のtip clipping等、
        /// 集合全体を舐めて再判定する処理から使う。
        /// </summary>
        public IEnumerable<byte[]> EnumerateTrustedKmers()
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            if (UseSmallPath)
            {
                foreach (var packed in this._trustedKmersSmall!.Keys)
                {
                    yield return UnpackSmall(packed, kmerLength);
                }
            }
            else if (UseMidPath)
            {
                foreach (var packed in this._trustedKmersMid!.Keys)
                {
                    yield return UnpackMid(packed, kmerLength);
                }
            }
            else
            {
                foreach (var key in this._trustedKmers!.Keys)
                {
                    yield return key.ToBytes(kmerLength);
                }
            }
        }

        /// <summary>
        /// kmerを信頼できるk-mer集合から除去する(順鎖・逆鎖どちらの向きで
        /// 渡してもよい)。GraphSimplifier がtipの構成k-merを取り除く際に使う。
        /// </summary>
        public void RemoveTrusted(ReadOnlySpan<byte> kmer)
        {
            if (UseSmallPath)
            {
                _ = this._trustedKmersSmall!.Remove(CanonicalSmall(kmer));
            }
            else if (UseMidPath)
            {
                _ = this._trustedKmersMid!.Remove(CanonicalMid(kmer));
            }
            else
            {
                _ = this._trustedKmers!.Remove(new KmerKey(kmer).Canonical());
            }
        }

        /// <summary>
        /// 現在の信頼できるk-mer集合を1回走査し、unitigの開始点となる
        /// 「入次数が1でない」k-merをすべて再検出する。ファイルの読み直しではなく
        /// インメモリの集合をそのまま使うため、tip clippingで集合を縮小した後の
        /// 再構築にも安価に使える。
        ///
        /// EnumerateTrustedKmersは各座位ごとに正規化された(canonicalな)
        /// 向きのk-merを1つだけ返すが、IsFirstKmerは向き依存(そのk-mer自身の
        /// 入次数を見る)の判定である。ある座位が「順鎖では分岐点の直後」でも
        /// 「逆鎖(=canonical側)では分岐点そのものではない」ことがありえるため、
        /// canonical側だけを調べると本来開始点であるべき向きを見逃す
        /// (小さなbubbleのテストケースで実際に見逃しを確認: 分岐解消後も
        /// 低カバレッジ側の枝が除去されないままになっていた)。
        /// そのため、各座位について両方の向き(kmerとその逆相補)を
        /// 個別に判定する。
        /// </summary>
        public List<byte[]> FindFirstKmers()
        {
            List<byte[]> kmers = [];
            foreach (var kmer in this.EnumerateTrustedKmers())
            {
                if (this.IsFirstKmer(kmer))
                {
                    kmers.Add(kmer);
                }

                var revComp = Util.ReverseComprement(kmer).ToArray();
                if (this.IsFirstKmer(revComp))
                {
                    kmers.Add(revComp);
                }
            }
            return kmers;
        }

        public List<byte[]> Cutoff(ulong bounds)
        {
            // 各ワーカーの CountingDB をそれぞれ MergeAll し、
            // 出来上がった複数のソート済みファイルをさらに1本にマージする。
            var mergedFiles = new List<string>();
            foreach (var counter in this._counters!)
            {
                mergedFiles.Add(counter.MergeAll());
                counter.Dispose();
            }
            this._counters = null;

            var filePath = CountingDB.MergeExternalFiles(this._tempDirectory, mergedFiles);

            var length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            var useSmallPath = UseSmallPath;
            var useMidPath = UseMidPath;
            var trustedKmers = useSmallPath || useMidPath ? null : new Dictionary<KmerKey, ulong>();
            var trustedKmersSmall = useSmallPath ? new Dictionary<ulong, ulong>() : null;
            var trustedKmersMid = useMidPath ? new Dictionary<UInt128, ulong>() : null;
            using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
            {
                ulong addedKmer = 0;
                ulong countKmer = 0;
                // count(k-merの出現回数) -> その回数を持つユニークk-merの種類数。
                // エラー由来の低頻度k-merと真のゲノム由来k-merを分ける「谷」を
                // 推定するために、カットオフ判定と同じこのループで集計する
                // (このファイルはこの後 File.Delete されるため、ここでしか見られない)。
                Dictionary<ulong, long> countHistogram = [];
                while (Util.HasNext(reader))
                {
                    var read = reader.ReadBytes(length);
                    var count = reader.ReadUInt64();
                    countKmer += 1;
                    countHistogram[count] = countHistogram.GetValueOrDefault(count, 0L) + 1;
                    if (count >= bounds)
                    {
                        addedKmer += 1;
                        List<byte> bytes = [];
                        foreach (var b in read)
                        {
                            bytes.AddRange(Util.ByteToNucleotideSequence(b));
                        }
                        var kmer = CollectionsMarshal.AsSpan(bytes)[..ConfigurationManager.Arguments.Kmer];
                        // カウント段階(ProcessRead)では正規化前の順鎖/逆鎖を
                        // 別々のキーとしてカウントしているため、同じ正規化k-merに
                        // 対応する2エントリ(順鎖側・逆鎖側)が別々にここへ来うる。
                        // 上書きではなく加算することで、両ストランド分の
                        // カバレッジを正しく合算する。
                        if (useSmallPath)
                        {
                            var canonical = CanonicalSmall(kmer);
                            trustedKmersSmall![canonical] = trustedKmersSmall.GetValueOrDefault(canonical, 0UL) + count;
                        }
                        else if (useMidPath)
                        {
                            var canonical = CanonicalMid(kmer);
                            trustedKmersMid![canonical] = trustedKmersMid.GetValueOrDefault(canonical, 0UL) + count;
                        }
                        else
                        {
                            var canonical = new KmerKey(kmer).Canonical();
                            trustedKmers![canonical] = trustedKmers.GetValueOrDefault(canonical, 0UL) + count;
                        }
                    }
                }
                Console.WriteLine("kmer count: " + countKmer);
                Console.WriteLine("good kmer: " + addedKmer);
                Console.WriteLine($"[Info] k-mer count histogram (count:#distinct kmers): {KmerHistogram.FormatSummary(countHistogram)}");
                var suggestedCutoff = KmerHistogram.SuggestCutoff(countHistogram);
                if (suggestedCutoff is { } suggestion)
                {
                    var note = suggestion == bounds ? " (matches the cutoff currently in effect)" : $" (currently using -kc {bounds})";
                    Console.WriteLine($"[Info] Suggested k-mer cutoff from histogram valley: {suggestion}{note}");
                }
                else
                {
                    Console.WriteLine("[Info] Could not identify a clear histogram valley to suggest a k-mer cutoff (spectrum may not be bimodal at this coverage).");
                }
            }
            File.Delete(filePath);
            this._trustedKmers = trustedKmers;
            this._trustedKmersSmall = trustedKmersSmall;
            this._trustedKmersMid = trustedKmersMid;

            Console.WriteLine("Search First k-mer");
            // 以前はここで一度カットオフ通過k-merをファイルへ書き出し、
            // 読み直して各k-merのIsFirstKmerを判定していた。厳密な集合を
            // インメモリで保持するようになった(Phase 1)ため、その集合を
            // 直接走査すれば同じ結果が得られ、ディスクI/Oを1往復省略できる。
            return this.FindFirstKmers();
        }

        /// <summary>
        /// 与えられた k-mer が unitig の開始点かどうかを判定する。
        ///
        /// 入次数が0または2個以上(=分岐点そのもの)であれば当然開始点になる。
        /// 入次数がちょうど1の場合でも、その唯一の予測元(predecessor)自身が
        /// 分岐点(出次数が1でない)であれば、UnitigMaker の前進walkは
        /// predecessorの時点で停止してしまいこのk-merへは到達しない
        /// (=このk-merは誰からも「walkで訪れてもらえない」)ため、
        /// 新たなunitigの開始点として別途扱う必要がある。
        /// これを見落とすと、分岐点の直後から始まる配列がunitig化されず
        /// 丸ごと欠落する(小さなテストケースで実際に発生を確認した)。
        /// </summary>
        public bool IsFirstKmer(Span<byte> kmer)
        {
            var inDegree = this.CountInEdges(kmer, out var uniquePredecessor);
            if (inDegree != 1)
            {
                return true;
            }
            return this.CountOutEdges(uniquePredecessor!) != 1;
        }

        /// <summary>
        /// kmer への入次数(前方に接続しうる異なる1塩基拡張の数)を数える。
        ///
        /// Core/UnitigMaker.cs の前進伸長規則は「kmerの先頭1文字を落とし、
        /// 末尾に候補塩基cを付加する」(successor = kmer[1..] + c)。
        /// この関係の逆(predecessor)を解くと、predecessor P は
        /// 「P[1..] + (kmerの末尾文字) == kmer」を満たす必要があり、
        /// P[1..] = kmer[..^1](kmerの末尾を落としたもの)、
        /// P[0] = 任意の候補塩基c、すなわち P = c + kmer[..^1] となる。
        ///
        /// 以前の実装は candidate = c + kmer[1..](kmerの"先頭"を落としたもの)
        /// を試しており、c = kmer[0] のとき candidate が kmer 自身と一致して
        /// しまう(=常に最低1回は自己ヒットする)退化バグがあった。これにより
        /// 真の入次数0(=配列の先頭)が絶対に検出できず、IsFirstKmer が
        /// 意図通りに開始点を拾えていなかった。
        /// </summary>
        public int CountInEdges(Span<byte> kmer)
        {
            return this.CountInEdges(kmer, out _);
        }

        /// <summary>
        /// CountInEdgesの本体。入次数がちょうど1だった場合、その唯一の
        /// predecessor(kmer長のbyte配列)も同時に返す(IsFirstKmerが使う)。
        /// </summary>
        private int CountInEdges(Span<byte> kmer, out byte[]? uniquePredecessor)
        {
            var candidate = new byte[kmer.Length];
            kmer[..^1].CopyTo(candidate.AsSpan(1));
            var count = 0;
            byte[]? match = null;
            for (byte i = 1; i <= 4; i++)
            {
                candidate[0] = i;
                if (this.Contains(candidate))
                {
                    count++;
                    match = count == 1 ? (byte[])candidate.Clone() : null;
                }
            }
            uniquePredecessor = count == 1 ? match : null;
            return count;
        }

        /// <summary>
        /// kmer からの出次数(後方に接続しうる異なる1塩基拡張の数)を数える。
        /// Core/UnitigMaker.cs の前進伸長規則(kmer[1..] + c)そのものを試す。
        /// GraphSimplifier がunitigの末尾端の次数(=tip判定)を見る際に使う。
        /// </summary>
        public int CountOutEdges(Span<byte> kmer)
        {
            var candidate = new byte[kmer.Length];
            kmer[1..].CopyTo(candidate.AsSpan(0, kmer.Length - 1));
            var count = 0;
            for (byte i = 1; i <= 4; i++)
            {
                candidate[^1] = i;
                if (this.Contains(candidate))
                {
                    count++;
                }
            }
            return count;
        }

        public void Dispose()
        {
            if (this._counters != null)
            {
                foreach (var counter in this._counters)
                {
                    counter.Dispose();
                }
            }
        }
    }
}
