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

        private readonly string unitigFilePath;

        // 単一リード内で直接検出された隣接(=k-1塩基のオーバーラップで
        // 実際に結合できる可能性が高いエッジ)。UniteContigs はこちらのみを使う。
        private readonly Dictionary<(int, int), ulong> kmerPath;

        // ペアエンド情報(read1/read2 がそれぞれ別 unitig にマップされたこと)由来の
        // 隣接。read1・read2 の間には既知の(が読まれていない)ギャップがあるため、
        // 単純な文字列オーバーラップでは結合できない。
        // 現時点ではスキャフォールディング機能が未実装のため UniteContigs には使わず、
        // 参考情報としてログに残すだけに留める。
        private readonly Dictionary<(int, int), ulong> pairPath;

        public ContigMaker(string unitigFilePath)
        {
            this.unitigFilePath = unitigFilePath;
            this.kmerDict = [];
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
        /// 「インサートサイズ程度の距離で隣接している」という情報として kmerPath に
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
            var localPairPaths = new Dictionary<(int, int), ulong>[threadCount];
            for (var i = 0; i < threadCount; i++)
            {
                localPaths[i] = [];
                localPairPaths[i] = [];
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
                            var id1 = this.FindDominantUnitig(read1);
                            var id2 = this.FindDominantUnitig(read2);
                            if (id1 != 0 && id2 != 0 && id1 != id2)
                            {
                                // read2 はリード分子の逆鎖側から読まれるため、
                                // read1 の向きに揃えるには read2 の逆相補鎖が
                                // 実際に「read1 の下流」に来る、という関係になる。
                                // read2 自体でヒットした unitig ID の符号を反転させ、
                                // read1 の向きに揃えたうえでペアを記録する。
                                var pathKey = (id1, -id2);
                                localPair[pathKey] = localPair.TryGetValue(pathKey, out var count) ? count + 1 : 1;
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
                foreach (var (pathKey, value) in localPair)
                {
                    this.pairPath[pathKey] = this.pairPath.TryGetValue(pathKey, out var existing) ? existing + value : value;
                }
            }

            // ペアエンド由来の隣接は現状 UniteContigs では使用しない
            // (スキャフォールディング機構が未実装のため)。参考情報として件数のみ出力する。
            Console.WriteLine($"[Info] Paired-end adjacency candidates detected: {this.pairPath.Count} (not used for direct unitig joining yet).");
        }

        /// <summary>
        /// 1本のリードが「代表として」どの unitig にマップされるかを判定する。
        /// リード中で最も安定して(連続して)ヒットし続けた unitig ID を返す。
        /// どの unitig にもヒットしなかった場合は 0 を返す。
        /// ペアエンドの隣接検出でのみ使用する軽量な単方向スキャン。
        /// </summary>
        private int FindDominantUnitig(string read)
        {
            if (string.IsNullOrEmpty(read))
            {
                return 0;
            }

            var kmerLength = ConfigurationManager.Arguments.Kmer;
            if (read.Length < kmerLength)
            {
                return 0;
            }

            var counts = new Dictionary<int, int>();
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
                    }
                }
            }

            if (counts.Count == 0)
            {
                return 0;
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
            return best;
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
            var visited = new bool[adjacencyList.Count];
            foreach (var index in firstUnitig)
            {
                var contig = WalkPath(unitigList, adjacencyList, index, visited);
                contigList.Add(contig);
            }
            for (var i = 2; i < adjacencyList.Count; i += 2)
            {
                if (!visited[i] && !visited[i + 1])
                {
                    contigList.Add(unitigList[i]);
                }
            }
            HashSet<string> set = [];
            using var writer = new FastaWriter(contigPath);
            var ID = 1;
            var genomeSize = 0L;
            foreach (var contig in contigList)
            {
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
                    if (contig.CompareTo(revContig) <= 0)
                    {
                        writer.Write($"NODE{ID}", contig);
                    }
                    else
                    {
                        writer.Write($"NODE{ID}", revContig);
                    }
                    ID++;
                    genomeSize += contig.Length;
                }
            }
            Console.WriteLine("Total Length of contigs : " + genomeSize);
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

        private static string WalkPath(List<string> unitigList, List<List<(int, ulong)>> adjacencyList, int index, bool[] visited)
        {
            // de Bruijn グラフ上の unitig 同士は理論上ちょうど k-1 塩基だけ
            // オーバーラップするはずである。以前の実装は
            // Math.Min(sb.Length, unitig.Length) から降順に試し、
            // 最初に一致した長さ(理論値より短い偶然の一致でも)で結合していたため、
            // 誤結合のリスクがあった。k-1 を最優先で試し、
            // それで一致しない場合のみ他の長さへフォールバックする。
            var expectedOverlap = ConfigurationManager.Arguments.Kmer - 1;

            var sb = new StringBuilder(unitigList[index]);
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