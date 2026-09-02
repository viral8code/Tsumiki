using System.Collections.Concurrent;
using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    internal class ContigMaker
    {
        // 同一 k-mer が複数の unitig にまたがって出現した(=リピート配列等に
        // 由来する曖昧な k-mer である)ことを示すセンチネル値。
        // unitig ID は 1 始まりの正数、逆鎖側はその負数を使うため int.MinValue と衝突しない。
        private const int AmbiguousKmer = int.MinValue;

        private readonly Dictionary<KmerKey, int> kmerDict;

        // unitig ID(1始まり) -> unitig の塩基長。ギャップ長推定で
        // 「unitig の末尾からリードのヒット位置までの残り長」を求めるのに使う。
        private readonly Dictionary<int, int> unitigLengths;

        private readonly string unitigFilePath;

        // 単一リード内で直接検出された隣接(=k-1塩基のオーバーラップで
        // 実際に結合できる可能性が高いエッジ)。UniteContigs はこちらのみを使う。
        private readonly Dictionary<(int, int), ulong> kmerPath;

        // ペアエンド情報(read1/read2 がそれぞれ別 unitig にマップされたこと)由来の
        // 隣接候補。キーは kmerPath と同じ (from, to) 形式(符号がunitigの向きを表す)。
        // 値は「観測されたペアの一覧」で、count(サポート数)に加えて各観測ごとの
        // ギャップ長推定値を保持し、Scaffolder 側で代表値(中央値)を計算できるようにする。
        private readonly Dictionary<(int, int), List<int>> pairPath;

        public ContigMaker(string unitigFilePath)
        {
            this.unitigFilePath = unitigFilePath;
            this.kmerDict = [];
            this.unitigLengths = [];
            this.kmerPath = [];
            this.pairPath = [];
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            using FastaReader reader = new(unitigFilePath);
            var id = 1;
            var skippedShortUnitigs = 0;
            var ambiguousCount = 0;
            while (reader.HasNext())
            {
                var unitig = reader.NextSequence();
                this.unitigLengths[id] = unitig.Seq.Length;
                if (unitig.Seq.Length < kmerLength)
                {
                    // k-mer 長より短い unitig は本来 k-mer を1つも持てないため、
                    // kmerDict に登録できずリードマッピングの対象から漏れる。
                    // (unitig 自体は MaximumUnitigCount との兼ね合いで一定数生じうる。)
                    // 完全な解決にはより短い k-mer での再マッピング等が必要だが、
                    // ここでは少なくとも「登録が1件もされないまま黙って id が進む」
                    // 状態を可視化するためカウントしておく。
                    skippedShortUnitigs++;
                    id++;
                    continue;
                }
                for (var i = kmerLength; i <= unitig.Seq.Length; i++)
                {
                    var key = new KmerKey(unitig.Seq.AsSpan(i - kmerLength, kmerLength));
                    var revKey = key.ReverseComprement();
                    ambiguousCount += RegisterKmer(this.kmerDict, key, id);
                    ambiguousCount += RegisterKmer(this.kmerDict, revKey, -id);
                }
                id++;
            }
            if (skippedShortUnitigs > 0)
            {
                Console.WriteLine($"[Warning] {skippedShortUnitigs} unitig(s) shorter than k-mer length were skipped in mapping.");
            }
            if (ambiguousCount > 0)
            {
                Console.WriteLine($"[Warning] {ambiguousCount} k-mer registration(s) were ambiguous (shared by multiple unitigs) and will be ignored during mapping.");
            }
        }

        /// <summary>
        /// kmerDict へ1件登録する。既に別の unitig ID が登録されていた場合、
        /// そのままでは後勝ちで上書きされてしまい、実際には異なる unitig 由来の
        /// リードが同じ ID にマップされたかのように誤って隣接関係を作ってしまう。
        /// これを防ぐため、衝突した k-mer は AmbiguousKmer としてマークし、
        /// マッピング時にはヒットとして扱わないようにする。
        /// 戻り値: 新たに曖昧マークを付けた場合は 1、そうでなければ 0。
        /// </summary>
        private static int RegisterKmer(Dictionary<KmerKey, int> dict, KmerKey key, int id)
        {
            if (dict.TryGetValue(key, out var existing))
            {
                if (existing == AmbiguousKmer || existing == id)
                {
                    return 0;
                }
                dict[key] = AmbiguousKmer;
                return 1;
            }
            dict[key] = id;
            return 0;
        }

        public void MappingRead(string readPath)
        {
            var threadCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);

            // kmerDict は構築後に変更されない読み取り専用データなので、
            // 複数スレッドから安全に参照できる。
            // kmerPath への書き込みはスレッドごとにローカルな辞書に集計し、
            // 最後にメインの kmerPath へマージすることでロックを避ける。
            var localPaths = new Dictionary<(int, int), ulong>[threadCount];
            for (var i = 0; i < threadCount; i++)
            {
                localPaths[i] = [];
            }

            using (var queue = new BlockingCollection<string>(boundedCapacity: threadCount * 256))
            {
                var workers = new Task[threadCount];
                for (var w = 0; w < threadCount; w++)
                {
                    var workerIndex = w;
                    workers[w] = Task.Run(() =>
                    {
                        var local = localPaths[workerIndex];
                        foreach (var read in queue.GetConsumingEnumerable())
                        {
                            this.MapSingleRead(read, local);
                        }
                    });
                }

                using (var reader = new FastqReader(readPath))
                {
                    while (reader.HasNext())
                    {
                        var read = reader.NextRead().RowRead;
                        queue.Add(read);
                    }
                }
                queue.CompleteAdding();

                Task.WaitAll(workers);
            }

            // 各ワーカーのローカル集計結果をメインの kmerPath へマージする。
            foreach (var local in localPaths)
            {
                foreach (var (pathKey, value) in local)
                {
                    this.kmerPath[pathKey] = this.kmerPath.TryGetValue(pathKey, out var existing) ? existing + value : value;
                }
            }
        }

        /// <summary>
        /// ペアエンドリードの情報を使って unitig 間の隣接関係を検出する。
        /// 単一リード内で複数 unitig をまたぐ場合(MappingRead と同様)に加え、
        /// read1 と read2 がそれぞれ別々の unitig に(単独で)マップされた場合も、
        /// 「インサートサイズ程度の距離で隣接している」という情報として pairPath に
        /// 追加登録する。これにより、unitig 長がリード長よりずっと長く
        /// 単一リードでは境界をまたげないケースでも隣接関係を検出できる。
        ///
        /// read1/read2 が本当にペアであることを、リード ID の対応(/1,/2 や
        /// Casava 1.8+ の "1:.../2:..." 記法)で検証する。対応が取れない場合は
        /// 警告を出し、そのペアはペアエンド由来の隣接検出をスキップする
        /// (単一リード内の隣接検出は通常どおり行う)。
        /// </summary>
        public void MappingPairedReads(string readPath1, string readPath2)
        {
            var threadCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);

            var localPaths = new Dictionary<(int, int), ulong>[threadCount];
            // localPairPaths: (from,to) -> このワーカーで観測した各ペアの
            // 「未読了長の合計(read1側残り長 + read2側残り長)」のリスト。
            // Scaffolder はこれと InsertSize から
            // gap = InsertSize - (read1側残り長 + read2側残り長) を計算する。
            var localPairPaths = new Dictionary<(int, int), List<int>>[threadCount];
            // インサートサイズ自動推定用: 両リードが同一unitigに単独マップされた
            // ペアについて、そのunitig内での距離をサンプリングする。
            // ライブラリの向きの組み合わせ(FR/RF/FF/RR)を決め打ちできないため、
            // 「符号が一致するヒット」と「符号が不一致のヒット」を別々に集計し、
            // 全体マージ後にサンプル数が多い方を採用する。
            var localSameOrientationSampleLists = new List<int>[threadCount];
            var localOppositeOrientationSampleLists = new List<int>[threadCount];
            for (var i = 0; i < threadCount; i++)
            {
                localPaths[i] = [];
                localPairPaths[i] = [];
                localSameOrientationSampleLists[i] = [];
                localOppositeOrientationSampleLists[i] = [];
            }

            using (var queue = new BlockingCollection<(string Read1, string Read2)>(boundedCapacity: threadCount * 256))
            {
                var workers = new Task[threadCount];
                for (var w = 0; w < threadCount; w++)
                {
                    var workerIndex = w;
                    workers[w] = Task.Run(() =>
                    {
                        var local = localPaths[workerIndex];
                        var localPair = localPairPaths[workerIndex];
                        var localSameOrientationSamples = localSameOrientationSampleLists[workerIndex];
                        var localOppositeOrientationSamples = localOppositeOrientationSampleLists[workerIndex];
                        foreach (var (read1, read2) in queue.GetConsumingEnumerable())
                        {
                            // 単一リード内の隣接検出は従来どおり両方に対して行う。
                            // (直接オーバーラップで結合できる可能性が高いエッジ。)
                            this.MapSingleRead(read1, local);
                            this.MapSingleRead(read2, local);

                            // ペアエンド情報による隣接検出: それぞれのリードが
                            // 単独でどの unitig にマップされるかを求め、
                            // 異なる unitig であれば「インサートサイズ程度の
                            // 距離で隣接している」という弱い証拠として記録する。
                            // こちらは直接のオーバーラップを保証しないため、
                            // kmerPath とは別の localPair に集計する。
                            var hit1 = this.FindDominantUnitig(read1);
                            var hit2 = this.FindDominantUnitig(read2);

                            if (hit1.UnitigId != 0 && hit2.UnitigId != 0 && Math.Abs(hit1.UnitigId) == Math.Abs(hit2.UnitigId))
                            {
                                // 両リードが同一unitigにマップされた場合、
                                // インサートサイズの実測サンプルとして使える。
                                // hit1/hit2 はいずれも「その unitig ID の符号が
                                // 示す向き」の座標系での LastMatchEndOffset を持つ。
                                //
                                // ペアエンドライブラリの向きの組み合わせ(FR/RF/FF/RR)は
                                // シーケンサ・ライブラリ調製方法に依存し、コード側で
                                // 一方を「正しい配置」と決め打ちすることはできない。
                                // そのため、符号が一致(同じ向き同士でヒット)する
                                // ケースと不一致(互いに逆向きでヒット)するケースの
                                // 両方についてサンプルを集めておき、実際にどちらが
                                // 多数派かを全ワーカー分集計してから判断する
                                // (このメソッドの最後で多数派側だけを採用する)。
                                var distance = Math.Abs(hit1.LastMatchEndOffset - hit2.LastMatchEndOffset);
                                if (distance > 0)
                                {
                                    if ((hit1.UnitigId > 0) == (hit2.UnitigId > 0))
                                    {
                                        localSameOrientationSamples.Add(distance);
                                    }
                                    else
                                    {
                                        localOppositeOrientationSamples.Add(distance);
                                    }
                                }
                            }
                            else if (hit1.UnitigId != 0 && hit2.UnitigId != 0 && hit1.UnitigId != hit2.UnitigId)
                            {
                                // read2 はリード分子の逆鎖側から読まれるため、
                                // read1 の向きに揃えるには read2 の逆相補鎖が
                                // 実際に「read1 の下流」に来る、という関係になる。
                                // read2 自体でヒットした unitig ID の符号を反転させ、
                                // read1 の向きに揃えたうえでペアを記録する。
                                var pathKey = (hit1.UnitigId, -hit2.UnitigId);

                                // ギャップ長推定に使う「未読了長」を計算する。
                                // unitig1 は hit1.UnitigId の符号の向きで見て、
                                // 読み取り済み末端(LastMatchEndOffset)から
                                // unitig の終端までの残り長(RemainingLength)が
                                // ギャップ側に残る未知区間の長さになる。
                                // unitig2 は pathKey.Item2(=-hit2.UnitigId)の
                                // 向きで見る必要があるが、hit2 は hit2.UnitigId の
                                // 向きで計算されているため、符号を反転させた
                                // 座標系に変換(=前後を入れ替える)する必要がある。
                                var remaining1 = hit1.RemainingLength;
                                var remaining2 = FlipOffsetToRemaining(hit2);

                                var pairLocal = localPair;
                                var totalRemaining = remaining1 + remaining2;
                                if (pairLocal.TryGetValue(pathKey, out var list))
                                {
                                    list.Add(totalRemaining);
                                }
                                else
                                {
                                    pairLocal[pathKey] = [totalRemaining];
                                }
                            }
                        }
                    });
                }

                using (var reader1 = new FastqReader(readPath1))
                using (var reader2 = new FastqReader(readPath2))
                {
                    var mismatchWarned = false;
                    while (reader1.HasNext() && reader2.HasNext())
                    {
                        var data1 = reader1.NextRead();
                        var data2 = reader2.NextRead();

                        var base1 = Util.GetPairedReadBaseId(data1.ID);
                        var base2 = Util.GetPairedReadBaseId(data2.ID);
                        if (base1 != base2)
                        {
                            if (!mismatchWarned)
                            {
                                Console.WriteLine($"[Warning] Paired read IDs do not match at this position (\"{data1.ID}\" vs \"{data2.ID}\"). " +
                                    "Paired-end adjacency detection may be unreliable for reads after this point; " +
                                    "single-read adjacency detection is unaffected.");
                                mismatchWarned = true;
                            }
                            // ペアが崩れている場合でも、単一リードとしての処理は続行する。
                            // (お互いを誤ってペアとして扱わないよう、キューには
                            //  「ペアなし」を示す空文字を積む。)
                            queue.Add((data1.RowRead, string.Empty));
                            queue.Add((data2.RowRead, string.Empty));
                            continue;
                        }

                        queue.Add((data1.RowRead, data2.RowRead));
                    }

                    // 片方のファイルだけ残っている場合は単一リードとして処理する。
                    while (reader1.HasNext())
                    {
                        queue.Add((reader1.NextRead().RowRead, string.Empty));
                    }
                    while (reader2.HasNext())
                    {
                        queue.Add((reader2.NextRead().RowRead, string.Empty));
                    }
                }
                queue.CompleteAdding();

                Task.WaitAll(workers);
            }

            foreach (var local in localPaths)
            {
                foreach (var (pathKey, value) in local)
                {
                    this.kmerPath[pathKey] = this.kmerPath.TryGetValue(pathKey, out var existing) ? existing + value : value;
                }
            }

            foreach (var localPair in localPairPaths)
            {
                foreach (var (pathKey, values) in localPair)
                {
                    if (this.pairPath.TryGetValue(pathKey, out var list))
                    {
                        list.AddRange(values);
                    }
                    else
                    {
                        this.pairPath[pathKey] = [.. values];
                    }
                }
            }

            // 「符号一致」「符号不一致」それぞれの総サンプル数を集計し、
            // 多数派の側だけを実際のライブラリ配置として InsertSizeSamples に採用する。
            // 少数派側は測定ノイズ・誤マッピング・稀な異常配置とみなして捨てる。
            var sameOrientationTotal = localSameOrientationSampleLists.Sum(l => l.Count);
            var oppositeOrientationTotal = localOppositeOrientationSampleLists.Sum(l => l.Count);

            IEnumerable<List<int>> chosenLists;
            string chosenLabel;
            if (sameOrientationTotal == 0 && oppositeOrientationTotal == 0)
            {
                chosenLists = [];
                chosenLabel = "none";
            }
            else if (sameOrientationTotal >= oppositeOrientationTotal)
            {
                chosenLists = localSameOrientationSampleLists;
                chosenLabel = "same-orientation";
            }
            else
            {
                chosenLists = localOppositeOrientationSampleLists;
                chosenLabel = "opposite-orientation";
            }

            foreach (var samples in chosenLists)
            {
                this.InsertSizeSamples.AddRange(samples);
            }

            var pairSupportCount = this.pairPath.Values.Sum(v => v.Count);
            Console.WriteLine($"[Info] Paired-end adjacency candidates detected: {this.pairPath.Count} edges ({pairSupportCount} supporting pairs total).");
            Console.WriteLine($"[Info] Same-unitig pair orientation counts: same-orientation={sameOrientationTotal}, opposite-orientation={oppositeOrientationTotal}. Using '{chosenLabel}' as the library's observed orientation for InsertSize estimation ({this.InsertSizeSamples.Count} samples).");
        }

        /// <summary>
        /// hit(ある unitig ID の符号が示す向きで計算された LastMatchEndOffset)を、
        /// その unitig を「逆向き」に見た座標系での RemainingLength に変換する。
        /// 元の向きでの RemainingLength(終端までの残り長)が、逆向きで見たときの
        /// 「先頭からの既知長」に相当するため、逆向きでの RemainingLength は
        /// 元の向きでの LastMatchEndOffset(先頭からの既知長)がそのまま使える。
        ///
        /// つまり: 向きを反転すると「先頭からの距離」と「末尾からの距離」が
        /// 入れ替わるので、反転後の RemainingLength = 元の LastMatchEndOffset。
        /// </summary>
        private static int FlipOffsetToRemaining(DominantUnitigHit hit)
        {
            return Math.Max(0, hit.LastMatchEndOffset);
        }

        /// <summary>
        /// InsertSize 自動推定用にサンプリングされた、単一unitig内での
        /// read1-read2 間距離のリスト。MappingPairedReads 実行後に
        /// Scaffolder から参照される。
        /// </summary>
        public List<int> InsertSizeSamples { get; } = [];

        /// <summary>
        /// ペアエンド由来の隣接候補。キーは (from, to) の unitig ID(符号は向き)、
        /// 値は各観測ペアのギャップ長推定用「未読了長の合計」のリスト。
        /// Scaffolder から参照される。
        /// </summary>
        public IReadOnlyDictionary<(int, int), List<int>> PairPath => this.pairPath;

        /// <summary>
        /// unitig ID(1始まり、符号なし)からその塩基長を引く。Scaffolder が
        /// contig 側の末端 unitig の長さを参照する際に使う。
        /// </summary>
        public IReadOnlyDictionary<int, int> UnitigLengths => this.unitigLengths;

        /// <summary>
        /// 1本のリードが「代表として」どの unitig にマップされるかを判定する。
        /// リード中で最も安定して(連続して)ヒットし続けた unitig ID に加え、
        /// スキャフォールディングのギャップ長推定に使うための最終ヒット位置も返す。
        /// どの unitig にもヒットしなかった場合は DominantUnitigHit.None を返す。
        /// ペアエンドの隣接検出でのみ使用する軽量な単方向スキャン。
        /// </summary>
        private DominantUnitigHit FindDominantUnitig(string read)
        {
            if (string.IsNullOrEmpty(read))
            {
                return DominantUnitigHit.None;
            }

            var kmerLength = ConfigurationManager.Arguments.Kmer;
            if (read.Length < kmerLength)
            {
                return DominantUnitigHit.None;
            }

            var counts = new Dictionary<int, int>();
            // 各候補 id ごとに、その id として最後にヒットしたk-merの
            // read内終端位置(exclusive, 0-based)を記録しておく。
            // read 中の相対位置は id の向きが unitig の順鎖(id>0)であれば
            // そのまま unitig 内オフセットとみなせる
            // (kmerDict は unitig の順鎖・逆鎖それぞれの k-mer をそのまま
            //  登録しているため、一致した時点の read 側の位置 = unitig 内の
            //  対応位置になる)。
            var lastReadEndOffset = new Dictionary<int, int>();
            var badBase = 0;
            for (var i = 0; i < kmerLength; i++)
            {
                if (Util.GetNucleotideIDs(read[i]).Count > 1)
                {
                    badBase++;
                }
            }
            for (var i = kmerLength; i <= read.Length; i++)
            {
                if (Util.GetNucleotideIDs(read[i - kmerLength]).Count > 1)
                {
                    badBase--;
                }
                if (badBase == 0)
                {
                    var key = new KmerKey(read.AsSpan(i - kmerLength, kmerLength));
                    if (this.kmerDict.TryGetValue(key, out var id) && id != AmbiguousKmer)
                    {
                        counts[id] = counts.GetValueOrDefault(id) + 1;
                        lastReadEndOffset[id] = i;
                    }
                }
            }

            if (counts.Count == 0)
            {
                return DominantUnitigHit.None;
            }

            var best = 0;
            var bestCount = 0;
            foreach (var (id, count) in counts)
            {
                if (count > bestCount)
                {
                    best = id;
                    bestCount = count;
                }
            }

            var unitigId = Math.Abs(best);
            var unitigLength = this.unitigLengths.GetValueOrDefault(unitigId, 0);
            var readEndOffset = lastReadEndOffset[best];

            // best > 0 の場合、read の k-mer は unitig の順鎖にそのまま一致して
            // いるため、read 内終端位置がそのまま unitig 内終端位置になる。
            // best < 0 の場合、read の k-mer は unitig の逆鎖(=unitig の
            // 逆相補鎖)に一致している。この場合、read を順方向に読み進めるほど
            // unitig の座標としては先頭側へ向かって進むことになるため、
            // 「unitig をその逆鎖の向きで見た座標系」での終端位置は
            // read の終端位置をそのまま使ってよい(逆鎖の kmerDict エントリは
            // 既に逆鎖の並びで登録されているため、座標系はその逆鎖基準になっている)。
            // つまりどちらの符号でも、read 内終端位置 = 「best の符号が示す
            // 向きで見た unitig 内終端位置」としてそのまま使える。
            var lastMatchEndOffset = readEndOffset;

            return new DominantUnitigHit(best, bestCount, lastMatchEndOffset, unitigLength);
        }

        /// <summary>
        /// 1リード分の k-mer マッピングを行い、隣接関係を(スレッドローカルな)
        /// localPath に集計する。並列化前の MappingRead の本体ロジックをそのまま
        /// 1リード単位の処理として切り出したもの。
        /// </summary>
        private void MapSingleRead(string read, Dictionary<(int, int), ulong> localPath)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            // FASTQ の生リードには N 等の曖昧塩基が混入しうるため、
            // A/C/G/T のみを前提とする ReverseComprement(string) ではなく
            // 曖昧塩基を許容する版を使う。曖昧塩基を含む区間の k-mer は
            // 後段の badBase/revBadBase によるスキップで除外される。
            var revRead = Util.ReverseComprementAllowAmbiguous(read);
            var bef = 0;
            var revBef = 0;
            var badBase = 0;
            var revBadBase = 0;
            for (var i = 0; i < kmerLength; i++)
            {
                if (Util.GetNucleotideIDs(read[i]).Count > 1)
                {
                    badBase++;
                }
                if (Util.GetNucleotideIDs(revRead[i]).Count > 1)
                {
                    revBadBase++;
                }
            }
            for (var i = kmerLength; i <= read.Length; i++)
            {
                if (Util.GetNucleotideIDs(read[i - kmerLength]).Count > 1)
                {
                    badBase--;
                }
                if (badBase == 0)
                {
                    var key = new KmerKey(read.AsSpan(i - kmerLength, kmerLength));
                    if (this.kmerDict.TryGetValue(key, out var id) && id != AmbiguousKmer)
                    {
                        if (bef == 0)
                        {
                            bef = id;
                        }
                        else if (bef != id)
                        {
                            var pathKey = (bef, id);
                            localPath[pathKey] = localPath.TryGetValue(pathKey, out var count) ? count + 1 : 1;
                            // 直前にヒットした unitig を更新する。これを怠ると、
                            // リード内で3つ以上の unitig にまたがった場合でも
                            // 常に「最初にヒットした unitig」との組しか記録されず、
                            // 実際の隣接関係(直前→直後)を反映できない。
                            bef = id;
                        }
                    }
                }
                if (Util.GetNucleotideIDs(revRead[i - kmerLength]).Count > 1)
                {
                    revBadBase--;
                }
                if (revBadBase == 0)
                {
                    var revKey = new KmerKey(revRead.AsSpan(i - kmerLength, kmerLength));
                    if (this.kmerDict.TryGetValue(revKey, out var revId) && revId != AmbiguousKmer)
                    {
                        if (revBef == 0)
                        {
                            revBef = revId;
                        }
                        else if (revBef != revId)
                        {
                            var pathKey = (revBef, revId);
                            localPath[pathKey] = localPath.TryGetValue(pathKey, out var count) ? count + 1 : 1;
                            // bef 側と同様、直前にヒットした unitig を更新する。
                            revBef = revId;
                        }
                    }
                }
            }
        }

        // unitig ID(1始まり) -> その unitig が最終的にどの contig の
        // どの位置に配置されたか。UniteContigs 実行後、Scaffolder から参照される。
        private readonly Dictionary<int, UnitigPlacement> unitigPlacements = [];

        public IReadOnlyDictionary<int, UnitigPlacement> UnitigPlacements => this.unitigPlacements;

        public void UniteContigs(string contigPath, decimal uniteThreshold, ulong countThreshold)
        {
            // リードマッピングによって得られた unitig 間の隣接情報(kmerPath)の
            // 規模を可視化する。ここが極端に少ない/空の場合、unitig の結合が
            // ほとんど起きず contigs.fasta が unitigs.fasta とほぼ同一になる。
            Console.WriteLine($"kmerPath entries (raw adjacency pairs found): {this.kmerPath.Count}");
            if (this.kmerPath.Count > 0)
            {
                var totalSupport = this.kmerPath.Values.Aggregate(0UL, (acc, v) => acc + v);
                var maxSupport = this.kmerPath.Values.Max();
                Console.WriteLine($"kmerPath total read support: {totalSupport}, max single-edge support: {maxSupport}");
            }

            List<string> unitigList = [string.Empty, string.Empty];
            var unitigCount = 0;
            using (FastaReader reader = new(this.unitigFilePath))
            {
                while (reader.HasNext())
                {
                    var unitig = reader.NextSequence().Seq;
                    unitigList.Add(unitig);
                    unitigList.Add(Util.ReverseComprement(unitig));
                    unitigCount++;
                }
            }
            List<List<(int, ulong)>> adjacencyList = [];
            for (var i = 0; i < unitigList.Count; i++)
            {
                adjacencyList.Add([]);
            }
            for (var i = 1; i <= unitigCount; i++)
            {
                for (var j = 1; j <= unitigCount; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    if (this.kmerPath.TryGetValue((i, j), out var count))
                    {
                        adjacencyList[i << 1].Add((j << 1, count));
                    }
                    if (this.kmerPath.TryGetValue((i, -j), out count))
                    {
                        adjacencyList[i << 1].Add((j << 1 | 1, count));
                    }
                    if (this.kmerPath.TryGetValue((-i, j), out count))
                    {
                        adjacencyList[i << 1 | 1].Add((j << 1, count));
                    }
                    if (this.kmerPath.TryGetValue((-i, -j), out count))
                    {
                        adjacencyList[i << 1 | 1].Add((j << 1 | 1, count));
                    }
                }
            }

            // FixPath 適用前の「候補を持つ頂点数」「合計候補数」を記録しておき、
            // countThreshold / uniteThreshold によってどれだけ絞り込まれた(あるいは
            // 消え去った)かを比較できるようにする。
            var verticesWithCandidatesBefore = 0;
            var totalCandidatesBefore = 0;
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                if (adjacencyList[i].Count > 0)
                {
                    verticesWithCandidatesBefore++;
                    totalCandidatesBefore += adjacencyList[i].Count;
                }
            }
            Console.WriteLine($"[Debug] Before FixPath: {verticesWithCandidatesBefore} vertices have candidate edges (total {totalCandidatesBefore} candidates).");

            for (var i = 2; i < adjacencyList.Count; i++)
            {
                FixPath(adjacencyList, i, uniteThreshold, countThreshold);
            }

            var verticesResolvedAfter = 0;
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                if (adjacencyList[i].Count == 1)
                {
                    verticesResolvedAfter++;
                }
            }
            Console.WriteLine($"[Debug] After FixPath: {verticesResolvedAfter} vertices resolved to exactly one edge (i.e. will actually be joined).");

            // InsertSize 自動推定: FixPath によって「実際に結合される」と確定した
            // unitig ペア(k-1 オーバーラップで直接連結される = ギャップなし)について、
            // pairPath に記録済みの「未読了長の合計」から
            // InsertSize = totalRemaining + (k-1) として逆算し、サンプルとして使う。
            // (pairPath 自体は MappingPairedReads 時点で「まだ結合されるかどうか
            //  わからない」候補として全件保持されているため、ここで FixPath の
            //  結果と突き合わせて実際に使うものだけを選び出す。)
            this.CollectInsertSizeSamplesFromResolvedEdges(adjacencyList);

            var enterCount = new int[adjacencyList.Count];
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                if (adjacencyList[i].Count == 1)
                {
                    enterCount[adjacencyList[i][0].Item1]++;
                }
            }

            // enterCount(その頂点を「唯一の行き先」として指している頂点の数)の分布。
            // 2以上になっている頂点が多い場合、複数の異なる unitig が同じ次の
            // unitig を指しており、WalkPath の visited 管理によって
            // 最初に辿り着いた1本しか結合されず、残りは孤立 unitig として
            // 個別出力される(=contig数が unitig 数を上回る一因になりうる)。
            var enterCountZero = 0;
            var enterCountOne = 0;
            var enterCountMulti = 0;
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                switch (enterCount[i])
                {
                    case 0:
                        enterCountZero++;
                        break;
                    case 1:
                        enterCountOne++;
                        break;
                    default:
                        enterCountMulti++;
                        break;
                }
            }
            Console.WriteLine($"[Debug] enterCount distribution: 0={enterCountZero}, 1={enterCountOne}, 2+={enterCountMulti}");

            var firstUnitig = new List<int>();
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                if (enterCount[i] == 0)
                {
                    firstUnitig.Add(i);
                }
            }
            Console.WriteLine($"[Debug] firstUnitig (walk start points) count: {firstUnitig.Count}");
            List<string> contigList = [];
            // 各 contig ごとの walk 順(vertexIndex のリスト)。vertexIndex は
            // adjacencyList の添字(unitig ID の << 1 / << 1|1 形式)であり、
            // unitig ID と向きの両方を含む。
            List<List<int>> walkOrders = [];
            var visited = new bool[adjacencyList.Count];
            foreach (var index in firstUnitig)
            {
                var walkOrder = new List<int>();
                var contig = WalkPath(unitigList, adjacencyList, index, visited, walkOrder);
                contigList.Add(contig);
                walkOrders.Add(walkOrder);
            }
            for (var i = 2; i < adjacencyList.Count; i += 2)
            {
                if (!visited[i] && !visited[i + 1])
                {
                    contigList.Add(unitigList[i]);
                    walkOrders.Add([i]);
                }
            }
            HashSet<string> set = [];
            using var writer = new FastaWriter(contigPath);
            var ID = 1;
            var genomeSize = 0L;
            for (var c = 0; c < contigList.Count; c++)
            {
                var contig = contigList[c];
                var walkOrder = walkOrders[c];
                if (set.Add(contig))
                {
                    var revContig = Util.ReverseComprement(contig);
                    if (contig != revContig)
                    {
                        if (!set.Add(revContig))
                        {
                            continue;
                        }
                    }
                    var isReverseComplemented = contig.CompareTo(revContig) > 0;
                    if (isReverseComplemented)
                    {
                        writer.Write($"NODE{ID}", revContig);
                    }
                    else
                    {
                        writer.Write($"NODE{ID}", contig);
                    }

                    // walkOrder に含まれる各頂点(unitig の向き付きインデックス)を
                    // unitigPlacements に記録する。walkOrder は「実際に配列へ
                    // 連結された順」なので、そのままこの contig 内での並び順になる。
                    // isReverseComplemented な場合、contigs.fasta 上の配列は
                    // walk 順と逆向きになっているため、位置(先頭/末尾)の解釈も
                    // 反転させる必要がある。ここでは「walk 順そのまま」の
                    // WalkOrderIndex を記録し、IsAtContigStart/End の判定はそのまま
                    // walk 順ベースで行い、Scaffolder 側で
                    // isContigReverseComplemented を見て解釈を反転させる。
                    for (var w = 0; w < walkOrder.Count; w++)
                    {
                        var vertexIndex = walkOrder[w];
                        var unitigId = vertexIndex >> 1;
                        var isReverseVertex = (vertexIndex & 1) == 1;
                        this.unitigPlacements[unitigId] = new UnitigPlacement(
                            contigId: ID,
                            isContigReverseComplemented: isReverseComplemented,
                            walkOrderIndex: w,
                            walkOrderCount: walkOrder.Count,
                            isUnitigReverseInWalk: isReverseVertex);
                    }

                    ID++;
                    genomeSize += contig.Length;
                }
            }
            Console.WriteLine("Total Length of contigs : " + genomeSize);
        }

        /// <summary>
        /// FixPath 後の adjacencyList を走査し、「頂点 v が唯一のエッジとして
        /// 頂点 next を指している」= v から next への直接結合が確定した、
        /// という関係を pairPath のキー形式(符号付き unitig ID のペア)に変換する。
        /// 変換後、該当する pairPath エントリの「未読了長合計」から
        /// InsertSize = totalRemaining + (k-1) を計算し、InsertSizeSamples に積む。
        ///
        /// adjacencyList の頂点インデックスは unitigId &lt;&lt; 1 (順鎖) /
        /// unitigId &lt;&lt; 1 | 1 (逆鎖) の形式。これを pairPath のキー形式
        /// (正の unitig ID = 順鎖, 負の unitig ID = 逆鎖)に変換する。
        /// </summary>
        private void CollectInsertSizeSamplesFromResolvedEdges(List<List<(int, ulong)>> adjacencyList)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            var overlap = kmerLength - 1;
            var collected = 0;

            for (var v = 2; v < adjacencyList.Count; v++)
            {
                if (adjacencyList[v].Count != 1)
                {
                    continue;
                }

                var next = adjacencyList[v][0].Item1;

                // 頂点インデックス -> 符号付き unitig ID。
                var fromUnitig = (v >> 1) * ((v & 1) == 0 ? 1 : -1);
                var toUnitig = (next >> 1) * ((next & 1) == 0 ? 1 : -1);

                if (!this.pairPath.TryGetValue((fromUnitig, toUnitig), out var totalRemainingSamples))
                {
                    continue;
                }

                foreach (var totalRemaining in totalRemainingSamples)
                {
                    var insertSize = totalRemaining + overlap;
                    if (insertSize > 0)
                    {
                        this.InsertSizeSamples.Add(insertSize);
                        collected++;
                    }
                }
            }

            Console.WriteLine($"[Info] InsertSize samples derived from resolved (actually-joined) unitig adjacency: {collected}.");
        }

        private static void FixPath(List<List<(int, ulong)>> adjacencyList, int index, decimal uniteThreshold, ulong countThreshold)
        {
            var pathList = adjacencyList[index];
            var sum = 0UL;
            for (var j = pathList.Count - 1; j >= 0; j--)
            {
                if (pathList[j].Item2 < countThreshold)
                {
                    pathList.RemoveAt(j);
                }
                else
                {
                    sum += pathList[j].Item2;
                }
            }
            (int, ulong)? path = null;
            var max = 0UL;
            foreach (var item in pathList)
            {
                if (max < item.Item2)
                {
                    max = item.Item2;
                    path = item;
                }
            }
            // uniteThreshold は「最多パスが全体の支持のうち何割を占めるか」の
            // 比率(例: 0.8 = 80%)として設計されている。
            // 以前は max(絶対リード数) と uniteThreshold(比率) をそのまま比較しており、
            // 実質 max >= 1 とほぼ同義になってしまっていた(countThreshold を
            // 通過した時点でほぼ常に真になる)。sum に対する比率で正しく判定する。
            adjacencyList[index] = sum > 0 && path != null && (decimal)max / sum >= uniteThreshold ? [((int, ulong))path!] : [];
        }

        private static string WalkPath(List<string> unitigList, List<List<(int, ulong)>> adjacencyList, int index, bool[] visited, List<int> walkOrder)
        {
            // de Bruijn グラフ上の unitig 同士は理論上ちょうど k-1 塩基だけ
            // オーバーラップするはずである。以前の実装は
            // Math.Min(sb.Length, unitig.Length) から降順に試し、
            // 最初に一致した長さ(理論値より短い偶然の一致でも)で結合していたため、
            // 誤結合のリスクがあった。k-1 を最優先で試し、
            // それで一致しない場合のみ他の長さへフォールバックする。
            var expectedOverlap = ConfigurationManager.Arguments.Kmer - 1;

            var sb = new StringBuilder(unitigList[index]);
            walkOrder.Add(index);
            while (adjacencyList[index].Count > 0 && !visited[index])
            {
                visited[index] = true;
                var next = adjacencyList[index][0].Item1;
                var unitig = unitigList[next];
                var flag = false;
                var maxLen = Math.Min(sb.Length, unitig.Length);

                // まず理論値(k-1)ちょうどのオーバーラップを試す。
                if (expectedOverlap > 0 && expectedOverlap <= maxLen &&
                    TryMatchOverlap(sb, unitig, expectedOverlap))
                {
                    _ = sb.Append(unitig[expectedOverlap..]);
                    flag = true;
                }
                else
                {
                    // 理論値で一致しなかった場合のフォールバックとして、
                    // 長い方から順に一致するオーバーラップを探す
                    // (元の実装と同じ挙動)。
                    for (var i = maxLen; i > 0; i--)
                    {
                        if (i == expectedOverlap)
                        {
                            // 上で既に試して失敗しているのでスキップ。
                            continue;
                        }
                        if (TryMatchOverlap(sb, unitig, i))
                        {
                            _ = sb.Append(unitig[i..]);
                            flag = true;
                            break;
                        }
                    }
                }

                if (!flag)
                {
                    break;
                }
                index = next;
                walkOrder.Add(index);
            }
            visited[index] = true;
            return sb.ToString();
        }

        /// <summary>
        /// sb の末尾 overlapLength 文字と unitig の先頭 overlapLength 文字が
        /// 一致するかどうかを判定する。
        /// </summary>
        private static bool TryMatchOverlap(StringBuilder sb, string unitig, int overlapLength)
        {
            var offset = sb.Length - overlapLength;
            for (var j = 0; j < overlapLength; j++)
            {
                if (sb[offset + j] != unitig[j])
                {
                    return false;
                }
            }
            return true;
        }
    }
}