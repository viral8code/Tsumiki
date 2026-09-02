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

        // Cutoff() 実行後に確定する、カットオフを通過したk-merの厳密な集合
        // (常に正規化(Canonical)された形で保持する)。以降のグラフ探索
        // (IsFirstKmer/CountInEdges/UnitigMakerの伸長判定)はすべて
        // この集合への厳密な所属判定のみで行う。
        private HashSet<KmerKey>? _trustedKmers;

        public TrustedKmerIndex(string tempDirectory)
        {
            this._tempDirectory = tempDirectory;
            var workerCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);
            this._counters = new CountingDB[workerCount];
            for (var i = 0; i < workerCount; i++)
            {
                this._counters[i] = new CountingDB(tempDirectory);
            }
        }

        /// <summary>
        /// 指定したワーカー番号(0 始まり)専用の CountingDB に登録する。
        /// 並列読み込み時、各スレッドは自分の workerIndex を固定して呼び出すことで
        /// スレッドセーフに(ロックなしで)k-mer を登録できる。
        /// </summary>
        public void Add(Span<byte[]> read, int workerIndex)
        {
            this._counters?[workerIndex].Add(read);
        }

        /// <summary>
        /// 指定したワーカー番号(0 始まり)専用の CountingDB に登録する。
        /// </summary>
        public void Add(Span<byte> read, int workerIndex)
        {
            this._counters?[workerIndex].Add(read);
        }

        /// <summary>
        /// kmer(順鎖・逆鎖いずれの向きでもよい)がカットオフを通過した
        /// 信頼できるk-mer集合に含まれるかどうかを厳密に判定する。
        /// </summary>
        public bool Contains(Span<byte> kmer)
        {
            return this._trustedKmers!.Contains(new KmerKey(kmer).Canonical());
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
            var kmerPath = Path.Combine(this._tempDirectory, Consts.KmerFileName);
            var trustedKmers = new HashSet<KmerKey>();
            using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
            {
                using var writer = new BinaryWriter(File.Open(kmerPath, FileMode.Create, FileAccess.Write));
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
                        _ = trustedKmers.Add(new KmerKey(kmer).Canonical());
                        writer.Write(kmer);
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

            Console.WriteLine("Search First k-mer");
            List<byte[]> kmers = [];
            using (var reader = new BinaryReader(File.Open(kmerPath, FileMode.Open, FileAccess.Read)))
            {
                while (Util.HasNext(reader))
                {
                    var read = reader.ReadBytes(ConfigurationManager.Arguments.Kmer);
                    if (this.IsFirstKmer(read))
                    {
                        kmers.Add(read);
                    }
                }
            }
            File.Delete(kmerPath);
            return kmers;
        }

        /// <summary>
        /// 与えられた k-mer が unitig の開始点(入次数が1でない = 0個または2個以上の
        /// prefix 拡張が存在する)かどうかを判定する。
        /// 入次数は「kmer の末尾 k-1 文字」の先頭に任意の1塩基を付加した k-mer が
        /// 信頼できるk-mer集合に存在するかで数える。
        /// </summary>
        private bool IsFirstKmer(Span<byte> kmer)
        {
            return this.CountInEdges(kmer) != 1;
        }

        /// <summary>
        /// kmer への入次数(前方に接続しうる異なる1塩基拡張の数)を数える。
        /// Contains が順鎖・逆鎖のどちらでも正しく判定するため、
        /// ここでは逆相補側への個別フォールバックは不要。
        /// </summary>
        private int CountInEdges(Span<byte> kmer)
        {
            var candidate = new byte[kmer.Length];
            kmer[1..].CopyTo(candidate.AsSpan(1));
            var count = 0;
            for (byte i = 1; i <= 4; i++)
            {
                candidate[0] = i;
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
