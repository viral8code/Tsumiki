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

        // 値は (符号付きunitig ID, そのunitig内でのk-mer開始位置(0始まり、
        // 符号が示す向きの座標系))。位置情報は FindDominantUnitig が
        // 「read内での最後のヒット位置」ではなく「unitig内での最後のヒット
        // 位置」を正しく求めるために必要(ギャップ長・インサートサイズ推定に使う)。
        private readonly Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict;

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
                    var startPos = i - kmerLength;
                    var key = new KmerKey(unitig.Seq.AsSpan(startPos, kmerLength));
                    var revKey = key.ReverseComprement();
                    // revKey は unitig 全体を逆相補した(=逆鎖の向きで読んだ)場合の
                    // 配列に対応する。区間 [startPos, startPos+kmerLength) を
                    // 長さ L の配列の逆側に写すと [L-i, L-startPos) になるため、
                    // 逆鎖側での開始位置は L-i。
                    var revStartPos = unitig.Seq.Length - i;
                    ambiguousCount += RegisterKmer(this.kmerDict, key, id, startPos);
                    ambiguousCount += RegisterKmer(this.kmerDict, revKey, -id, revStartPos);
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
        private static int RegisterKmer(Dictionary<KmerKey, (int, int)> dict, KmerKey key, int id, int position)
        {
            if (dict.TryGetValue(key, out var existing))
            {
                if (existing.Item1 == AmbiguousKmer || existing.Item1 == id)
                {
                    return 0;
                }
                dict[key] = (AmbiguousKmer, 0);
                return 1;
            }
            dict[key] = (id, position);
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
                                //
                                // hit1/hit2 の LastMatchEndOffset は、それぞれの
                                // ヒット自身の符号が示す向きの座標系(順鎖なら
                                // unitig先頭起点、逆鎖ならunitigを逆相補した
                                // 向きの起点)で測られている。符号が一致しない
                                // (=一方は順鎖、他方は逆鎖でヒットした)場合、
                                // 座標系が異なる2つの値をそのまま引き算しても
                                // 意味のある距離にはならない。ToForwardFrame で
                                // 両方を共通の(順鎖)座標系に変換してから
                                // 差を取る(順鎖同士・逆鎖同士の場合はこの変換は
                                // 距離を変えないため、常にこの経路で問題ない)。
                                //
                                // ペアエンドライブラリの向きの組み合わせ(FR/RF/FF/RR)は
                                // シーケンサ・ライブラリ調製方法に依存し、コード側で
                                // 一方を「正しい配置」と決め打ちすることはできない。
                                // そのため、符号が一致(同じ向き同士でヒット)する
                                // ケースと不一致(互いに逆向きでヒット)するケースの
                                // 両方についてサンプルを集めておき、実際にどちらが
                                // 多数派かを全ワーカー分集計してから判断する
                                // (このメソッドの最後で多数派側だけを採用する)。
                                // ToForwardFrame は「そのヒットの向きで見た既知長」を
                                // 順鎖座標へ写した値、すなわち順鎖から見たリードの
                                // 「内側の端」の座標になる。したがって2つの差は
                                // フラグメント長ではなく、2リードに挟まれた内側の
                                // 未読区間(inner distance)の長さである。
                                //
                                // 実データ(150bpリード・IS350ライブラリ)で
                                // この差の中央値が58になり、リード長150bpより
                                // 短いという物理的にありえない推定値になっていた。
                                // 内側距離 d と真のフラグメント長 F の関係は
                                // FR配置で F = d + len(read1) + len(read2) であり、
                                // 58 + 150 + 150 = 358 でライブラリ名(IS350)と一致する。
                                // ここで両リード長を足し戻し、以降の推定値が
                                // 一貫して「フラグメント長」の単位になるようにする。
                                if ((hit1.UnitigId > 0) == (hit2.UnitigId > 0))
                                {
                                    // 同じ向き同士(FF/RR相当)。両リードの内側の端は
                                    // どちらも同じ側を向いているため、差は
                                    // 「開始位置の差」に相当する。下流側リード1本分を
                                    // 足すとフラグメント長になる。
                                    var offsetDistance = Math.Abs(ToForwardFrame(hit1) - ToForwardFrame(hit2));
                                    var fragment = offsetDistance + Math.Max(read1.Length, read2.Length);
                                    if (fragment > 0)
                                    {
                                        localSameOrientationSamples.Add(fragment);
                                    }
                                }
                                else
                                {
                                    // 互いに逆向き(FR相当、Illuminaペアエンドの通常配置)。
                                    // 順鎖側ヒットのリードがフラグメントの左端、
                                    // 逆鎖側ヒットのリードが右端を占める。
                                    var hit1IsForward = hit1.UnitigId > 0;
                                    var forwardEnd = ToForwardFrame(hit1IsForward ? hit1 : hit2);
                                    var reverseStart = ToForwardFrame(hit1IsForward ? hit2 : hit1);
                                    var forwardReadLength = hit1IsForward ? read1.Length : read2.Length;
                                    var reverseReadLength = hit1IsForward ? read2.Length : read1.Length;

                                    // フラグメントの左端 = 順鎖リードの開始位置、
                                    // 右端 = 逆鎖リードの終了位置。
                                    var fragment = (reverseStart + reverseReadLength) - (forwardEnd - forwardReadLength);
                                    if (fragment > 0)
                                    {
                                        localOppositeOrientationSamples.Add(fragment);
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

                                // 記録するのは「フラグメントのうち、2つのunitigの
                                // 内側に既に見えている分の長さ」:
                                //   read1の長さ + unitig1末端までの残り
                                //   + unitig2先頭からの残り + read2の長さ
                                // 未知区間(ギャップ)長を G とすると
                                //   フラグメント長 = この値 + G
                                // という関係が常に成り立つ(直接k-1で結合された
                                // 場合は G = -(k-1))。
                                //
                                // 以前は read1/read2 の長さを含めない
                                // remaining1 + remaining2 だけを記録していたため、
                                // ここから逆算されるインサートサイズもギャップ長も
                                // 両リード長ぶん(実データで300bp)ずれていた。
                                var pairLocal = localPair;
                                var spannedLength = remaining1 + remaining2 + read1.Length + read2.Length;
                                if (pairLocal.TryGetValue(pathKey, out var list))
                                {
                                    list.Add(spannedLength);
                                }
                                else
                                {
                                    pairLocal[pathKey] = [spannedLength];
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

            var sameUnitigSamples = new List<int>();
            foreach (var samples in chosenLists)
            {
                sameUnitigSamples.AddRange(samples);
            }
            this.InsertSizeSamples.AddRange(sameUnitigSamples);
            this.SameUnitigInsertSizeSamples.AddRange(sameUnitigSamples);

            var pairSupportCount = this.pairPath.Values.Sum(v => v.Count);
            Console.WriteLine($"[Info] Paired-end adjacency candidates detected: {this.pairPath.Count} edges ({pairSupportCount} supporting pairs total).");
            Console.WriteLine($"[Info] Same-unitig pair orientation counts: same-orientation={sameOrientationTotal}, opposite-orientation={oppositeOrientationTotal}. Using '{chosenLabel}' as the library's observed orientation for InsertSize estimation ({sameUnitigSamples.Count} samples).");
            if (sameUnitigSamples.Count > 0)
            {
                // 同一unitig内サンプルは、unitig自体がフラグメント長より短い場合
                // 両端が同じunitig内に収まるペアしか観測できず、より短い
                // フラグメントに偏った標本になりやすい(unitigが短いほど顕著)。
                // resolved-edge由来の中央値(下のCollectInsertSizeSamplesFromMerges
                // が出力)と比較することで、このバイアスの有無を確認できる。
                Console.WriteLine($"[Info] Same-unitig fragment-length distribution: {FormatDistribution(sameUnitigSamples)}.");
                Console.WriteLine($"[Info] Same-unitig fragment-length median: {Median(sameUnitigSamples)} (from {sameUnitigSamples.Count} samples; read lengths added back to the inner distance, so this is a true fragment length. May still be biased short if unitigs are shorter than the true insert size).");
            }
        }

        private static int Median(List<int> values)
        {
            var sorted = values.OrderBy(x => x).ToList();
            var mid = sorted.Count / 2;
            return sorted.Count % 2 == 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
        }

        /// <summary>
        /// フラグメント長分布の分位点を要約する。中央値だけでは
        /// 「このライブラリがどれだけの長さのギャップを跨げるか」が分からない。
        /// リード長の2倍を超える分だけがスキャフォールディングで橋渡しできる
        /// 未知区間の長さなので、分布の裾(特に上側)が実際の橋渡し能力を決める。
        /// </summary>
        private static string FormatDistribution(List<int> values)
        {
            var sorted = values.OrderBy(x => x).ToList();
            int At(double q) => sorted[Math.Clamp((int)(q * (sorted.Count - 1)), 0, sorted.Count - 1)];
            return $"p1={At(0.01)}, p10={At(0.10)}, p25={At(0.25)}, p50={At(0.50)}, " +
                $"p75={At(0.75)}, p90={At(0.90)}, p99={At(0.99)}, max={sorted[^1]}";
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
        /// hit の LastMatchEndOffset(hit自身の符号が示す向きの座標系での値)を、
        /// 常に unitig の「順鎖」座標系での位置に変換する。順鎖ヒットは
        /// そのまま、逆鎖ヒットは UnitigLength - LastMatchEndOffset に変換する
        /// (順鎖・逆鎖どちらも同じ変換を経由するため、同じ符号同士を比較する
        /// 場合でも距離は変わらず、異符号同士の比較でも正しく意味を持つ)。
        /// 同一unitig上の2ヒット間の距離(インサートサイズ推定)を求める際、
        /// 座標系を揃えるために使う。
        /// </summary>
        private static int ToForwardFrame(DominantUnitigHit hit)
        {
            return hit.UnitigId > 0 ? hit.LastMatchEndOffset : hit.UnitigLength - hit.LastMatchEndOffset;
        }

        /// <summary>
        /// InsertSize 自動推定用にサンプリングされた距離のリスト
        /// (SameUnitigInsertSizeSamples と ResolvedEdgeInsertSizeSamples の結合)。
        /// 後方互換のため残しているが、Scaffolder はサンプリング元による
        /// バイアスの違いを考慮するため個別のリストを優先的に参照する。
        /// </summary>
        public List<int> InsertSizeSamples { get; } = [];

        /// <summary>
        /// 単一unitig内で両リードがヒットしたペアからのサンプル。
        /// unitig自体がフラグメント長より短い場合、両端が収まるペアしか
        /// 観測できないため、より短いフラグメントに偏りやすい
        /// (unitigが短いほど顕著)。
        /// </summary>
        public List<int> SameUnitigInsertSizeSamples { get; } = [];

        /// <summary>
        /// unitig同士がk-1オーバーラップで直接結合されたペアからのサンプル。
        /// SameUnitigInsertSizeSamples のような長さバイアスを受けない。
        /// </summary>
        public List<int> ResolvedEdgeInsertSizeSamples { get; } = [];

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
        internal DominantUnitigHit FindDominantUnitig(string read)
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
            // 「unitig内での」終端位置(exclusive, 0-based, idの符号が示す
            // 向きの座標系)を記録する。kmerDict が (unitigId, unitig内開始位置)
            // を保持するようになったため、read内での相対位置ではなく
            // kmerDict から得た本物のunitig内位置を使う
            // (以前はread内終端位置をそのままunitig内終端位置として誤用しており、
            //  unitigがread長より十分短い場合はたまたま近い値になり問題が
            //  表面化しにくかったが、unitigが長くなると全く違う値になっていた)。
            var lastUnitigEndOffset = new Dictionary<int, int>();
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
                    if (this.kmerDict.TryGetValue(key, out var entry) && entry.UnitigId != AmbiguousKmer)
                    {
                        var id = entry.UnitigId;
                        counts[id] = counts.GetValueOrDefault(id) + 1;
                        lastUnitigEndOffset[id] = entry.Position + kmerLength;
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
            var lastMatchEndOffset = lastUnitigEndOffset[best];

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
                    if (this.kmerDict.TryGetValue(key, out var entry) && entry.UnitigId != AmbiguousKmer)
                    {
                        var id = entry.UnitigId;
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
                    if (this.kmerDict.TryGetValue(revKey, out var revEntry) && revEntry.UnitigId != AmbiguousKmer)
                    {
                        var revId = revEntry.UnitigId;
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
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            var overlap = kmerLength - 1;

            List<string> unitigList = [string.Empty, string.Empty];
            using (FastaReader reader = new(this.unitigFilePath))
            {
                while (reader.HasNext())
                {
                    var unitig = reader.NextSequence().Seq;
                    unitigList.Add(unitig);
                    unitigList.Add(Util.ReverseComprement(unitig));
                }
            }

            // 隣接は de Bruijn グラフから厳密に導く(UnitigGraph の説明を参照)。
            // リードマッピング由来の kmerPath は「辺を作る」ためではなく、
            // 分岐点でどの辺を選ぶかの「重み」としてのみ使う。
            var graph = UnitigGraph.Build(unitigList, this.kmerDict, kmerLength, AmbiguousKmer);

            var edgeCount = 0;
            var branchingVertices = 0;
            for (var v = 2; v < graph.VertexCount; v++)
            {
                edgeCount += graph.OutEdges[v].Count;
                if (graph.OutEdges[v].Count > 1)
                {
                    branchingVertices++;
                }
            }
            Console.WriteLine($"[Debug] Exact de Bruijn unitig graph: {edgeCount} directed edge(s), {branchingVertices} branching vertex(es) out of {graph.VertexCount - 2}.");
            Console.WriteLine($"kmerPath entries (raw read-support pairs found): {this.kmerPath.Count}");

            // リード由来の支持数を逆鎖対称に集計する。辺 v→w と w^1→v^1 は
            // 同一の物理的な隣接を表すため、重みも同一でなければ順鎖側と
            // 逆鎖側で異なる経路が選ばれ、同じ領域が 2 通りに組み立てられてしまう。
            Dictionary<(int, int), ulong> support = [];
            foreach (var ((from, to), count) in this.kmerPath)
            {
                if (from == to)
                {
                    continue;
                }
                var v = VertexIndex(from);
                var w = VertexIndex(to);
                support[(v, w)] = support.GetValueOrDefault((v, w)) + count;
                support[(w ^ 1, v ^ 1)] = support.GetValueOrDefault((w ^ 1, v ^ 1)) + count;
            }

            // 単純バブルを潰してから辺を選ぶ。相互一意性を課す以上、
            // 再合流点の入次数が2以上のまま残っているとその経路全体が
            // 結合されなくなるため、先に枝を1本に絞っておく必要がある。
            var poppedBubbles = graph.PopSimpleBubbles(unitigList, support);
            if (poppedBubbles > 0)
            {
                Console.WriteLine($"[Debug] Popped {poppedBubbles} simple bubble branch(es) (kept as standalone contigs; only their graph edges were removed).");
            }

            // 各頂点について「出て行く先」を高々 1 つに絞る。
            var chosen = new int[graph.VertexCount];
            Array.Fill(chosen, -1);
            var unambiguous = 0;
            var resolvedByReads = 0;
            for (var v = 2; v < graph.VertexCount; v++)
            {
                var outs = graph.OutEdges[v];
                if (outs.Count == 0)
                {
                    continue;
                }
                if (outs.Count == 1)
                {
                    chosen[v] = outs[0];
                    unambiguous++;
                    continue;
                }

                var sum = 0UL;
                var best = -1;
                var bestCount = 0UL;
                foreach (var w in outs)
                {
                    var c = support.GetValueOrDefault((v, w));
                    sum += c;
                    if (c > bestCount)
                    {
                        bestCount = c;
                        best = w;
                    }
                }
                if (best >= 0 && bestCount >= countThreshold && sum > 0 && (decimal)bestCount / sum >= uniteThreshold)
                {
                    chosen[v] = best;
                    resolvedByReads++;
                }
            }

            // 相互一意(mutual unique)な辺だけを実際の結合として採用する。
            // v→w を結合してよいのは「v の唯一の行き先が w」であり、かつ
            // 「w の唯一の来訪元が v」であるときに限る。後者は逆鎖対称性より
            // chosen[w^1] == v^1 と同値。この条件を欠くと、複数の異なる
            // unitig が同じ次の unitig を指し(実データで 1550 頂点)、
            // 先着 1 本だけが結合されて残りが千切れる形になっていた。
            var merge = new int[graph.VertexCount];
            Array.Fill(merge, -1);
            var mergeCount = 0;
            for (var v = 2; v < graph.VertexCount; v++)
            {
                var w = chosen[v];
                if (w >= 0 && chosen[w ^ 1] == (v ^ 1))
                {
                    merge[v] = w;
                    mergeCount++;
                }
            }
            Console.WriteLine($"[Debug] Edge selection: {unambiguous} vertex(es) had a single out-edge, {resolvedByReads} branch(es) resolved by read support; {mergeCount} directed merge(s) survived the mutual-uniqueness check ({mergeCount / 2} undirected join(s)).");

            this.CollectInsertSizeSamplesFromMerges(merge);

            // 双子(v と v^1)は同一 unitig の裏表なので、unitig 単位で訪問済みを
            // 管理する。これを頂点単位でやっていたため、順鎖側の walk と逆鎖側の
            // walk が同じ unitig を別々に出力し、contig 総長が unitig 総長の
            // ちょうど 2 倍に膨れていた。
            var unitigCount = (unitigList.Count - 2) / 2;
            var unitigVisited = new bool[unitigCount + 1];

            List<string> contigList = [];
            List<List<int>> walkOrders = [];

            string Walk(int startVertex, List<int> walkOrder)
            {
                var sb = new StringBuilder(unitigList[startVertex]);
                walkOrder.Add(startVertex);
                unitigVisited[startVertex >> 1] = true;
                var cur = startVertex;
                while (true)
                {
                    var next = merge[cur];
                    if (next < 0 || unitigVisited[next >> 1])
                    {
                        break;
                    }
                    var seq = unitigList[next];
                    if (seq.Length < overlap || sb.Length < overlap)
                    {
                        break;
                    }
                    // 構築方法より k-1 のオーバーラップは保証されているが、
                    // 万一崩れていた場合に誤った配列を作らないよう検証する。
                    if (!TryMatchOverlap(sb, seq, overlap))
                    {
                        break;
                    }
                    _ = sb.Append(seq[overlap..]);
                    unitigVisited[next >> 1] = true;
                    walkOrder.Add(next);
                    cur = next;
                }
                return sb.ToString();
            }

            // 結合グラフ上で「入ってくる結合を持たない」頂点が経路の始点。
            // v への結合が存在することは、逆鎖対称性より merge[v^1] != -1 と同値。
            for (var v = 2; v < graph.VertexCount; v++)
            {
                if (merge[v ^ 1] != -1 || unitigVisited[v >> 1])
                {
                    continue;
                }
                List<int> walkOrder = [];
                contigList.Add(Walk(v, walkOrder));
                walkOrders.Add(walkOrder);
            }

            // 始点を持たない=循環している経路を拾う(環状ゲノム/プラスミド等)。
            for (var v = 2; v < graph.VertexCount; v += 2)
            {
                if (unitigVisited[v >> 1])
                {
                    continue;
                }
                List<int> walkOrder = [];
                contigList.Add(Walk(v, walkOrder));
                walkOrders.Add(walkOrder);
            }

            using var writer = new FastaWriter(contigPath);
            var ID = 1;
            var genomeSize = 0L;
            for (var c = 0; c < contigList.Count; c++)
            {
                var contig = contigList[c];
                var walkOrder = walkOrders[c];
                var revContig = Util.ReverseComprement(contig);
                var isReverseComplemented = string.CompareOrdinal(contig, revContig) > 0;
                writer.Write($"NODE{ID}", isReverseComplemented ? revContig : contig);

                // walkOrder に含まれる各頂点(unitig の向き付きインデックス)を
                // unitigPlacements に記録する。walkOrder は「実際に配列へ
                // 連結された順」なので、そのままこの contig 内での並び順になる。
                // isReverseComplemented な場合、contigs.fasta 上の配列は
                // walk 順と逆向きになっているため、位置(先頭/末尾)の解釈は
                // Scaffolder 側で isContigReverseComplemented を見て反転させる。
                for (var w = 0; w < walkOrder.Count; w++)
                {
                    var vertexIndex = walkOrder[w];
                    this.unitigPlacements[vertexIndex >> 1] = new UnitigPlacement(
                        contigId: ID,
                        isContigReverseComplemented: isReverseComplemented,
                        walkOrderIndex: w,
                        walkOrderCount: walkOrder.Count,
                        isUnitigReverseInWalk: (vertexIndex & 1) == 1);
                }

                ID++;
                genomeSize += contig.Length;
            }
            Console.WriteLine("Total Length of contigs : " + genomeSize);
        }

        /// <summary>
        /// 相互一意性の検査を通って実際に結合が確定した辺(merge[v] = next)を
        /// 走査し、その関係を pairPath のキー形式(符号付き unitig ID のペア)に
        /// 変換する。変換後、該当する pairPath エントリの「未読了長合計」から
        /// InsertSize = totalRemaining + (k-1) を計算し、InsertSizeSamples に積む。
        ///
        /// merge の頂点インデックスは unitigId &lt;&lt; 1 (順鎖) /
        /// unitigId &lt;&lt; 1 | 1 (逆鎖) の形式。これを pairPath のキー形式
        /// (正の unitig ID = 順鎖, 負の unitig ID = 逆鎖)に変換する。
        /// </summary>
        private void CollectInsertSizeSamplesFromMerges(int[] merge)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            var overlap = kmerLength - 1;
            List<int> resolvedEdgeSamples = [];

            for (var v = 2; v < merge.Length; v++)
            {
                var next = merge[v];
                if (next < 0)
                {
                    continue;
                }

                // 頂点インデックス -> 符号付き unitig ID。
                var fromUnitig = (v >> 1) * ((v & 1) == 0 ? 1 : -1);
                var toUnitig = (next >> 1) * ((next & 1) == 0 ? 1 : -1);

                if (!this.pairPath.TryGetValue((fromUnitig, toUnitig), out var spannedLengthSamples))
                {
                    continue;
                }

                foreach (var spannedLength in spannedLengthSamples)
                {
                    // 直接結合されたエッジでは2つのunitigがk-1塩基重なるので、
                    // 未知区間の長さは G = -(k-1)。よって
                    // フラグメント長 = spannedLength - (k-1)。
                    // (以前は符号を逆にした + (k-1) を使っており、
                    //  リード長ぶんのずれと合わせて推定値が大きく外れていた。)
                    var insertSize = spannedLength - overlap;
                    if (insertSize > 0)
                    {
                        resolvedEdgeSamples.Add(insertSize);
                    }
                }
            }

            this.InsertSizeSamples.AddRange(resolvedEdgeSamples);
            this.ResolvedEdgeInsertSizeSamples.AddRange(resolvedEdgeSamples);

            Console.WriteLine($"[Info] InsertSize samples derived from resolved (actually-joined) unitig adjacency: {resolvedEdgeSamples.Count}.");
            if (resolvedEdgeSamples.Count > 0)
            {
                // このプールは「unitig同士がk-1オーバーラップで直接結合された」
                // ペアのみを対象とするため、same-unitigサンプルのような
                // 「フラグメントが1つのunitigに収まる必要がある」制約が
                // なく、短いunitigによる短フラグメントへの偏りを受けにくい。
                Console.WriteLine($"[Info] Resolved-edge sample median: {Median(resolvedEdgeSamples)} (from {resolvedEdgeSamples.Count} samples; not subject to the same-unitig length bias).");
            }
        }

        /// <summary>
        /// 符号付き unitig ID(正=順鎖、負=逆鎖)を adjacencyList の頂点インデックスに変換する。
        /// </summary>
        internal static int VertexIndex(int signedUnitigId)
        {
            return (Math.Abs(signedUnitigId) << 1) | (signedUnitigId > 0 ? 0 : 1);
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