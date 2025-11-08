using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    internal class ContigMaker
    {
        private readonly Dictionary<KmerKey, int> kmerDict;

        private readonly string unitigFilePath;

        private readonly Dictionary<(int, int), ulong> kmerPath;

        public ContigMaker(string unitigFilePath)
        {
            this.unitigFilePath = unitigFilePath;
            this.kmerDict = [];
            this.kmerPath = [];
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            using FastaReader reader = new(unitigFilePath);
            var id = 1;
            while (reader.HasNext())
            {
                var unitig = reader.NextSequence();
                for (var i = kmerLength; i <= unitig.Seq.Length; i++)
                {
                    var key = new KmerKey(unitig.Seq.AsSpan(i - kmerLength, kmerLength));
                    var revKey = key.ReverseComprement();
                    this.kmerDict[key] = id;
                    this.kmerDict[revKey] = -id;
                }
                id++;
            }
        }

        public void MappingRead(string readPath)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            using FastqReader reader = new(readPath);
            while (reader.HasNext())
            {
                var read = reader.NextRead().RowRead;
                var revRead = Util.ReverseComprement(read);
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
                        if (this.kmerDict.TryGetValue(key, out var id))
                        {
                            if (bef == 0)
                            {
                                bef = id;
                            }
                            else if (bef != id)
                            {
                                var pathKey = (bef, id);
                                if (this.kmerPath.ContainsKey(pathKey))
                                {
                                    this.kmerPath[pathKey] += 1;
                                }
                                else
                                {
                                    this.kmerPath[pathKey] = 1;
                                }
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
                        if (this.kmerDict.TryGetValue(revKey, out var revId))
                        {
                            if (revBef == 0)
                            {
                                revBef = revId;
                            }
                            else if (revBef != revId)
                            {
                                var pathKey = (revBef, revId);
                                if (this.kmerPath.ContainsKey(pathKey))
                                {
                                    this.kmerPath[pathKey] += 1;
                                }
                                else
                                {
                                    this.kmerPath[pathKey] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        public void UniteContigs(string contigPath, decimal uniteThreshold, ulong countThreshold)
        {
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

            for (var i = 2; i < adjacencyList.Count; i++)
            {
                FixPath(adjacencyList, i, uniteThreshold, countThreshold);
            }

            var enterCount = new int[adjacencyList.Count];
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                if (adjacencyList[i].Count == 1)
                {
                    enterCount[adjacencyList[i][0].Item1]++;
                }
            }
            var firstUnitig = new List<int>();
            for (var i = 2; i < adjacencyList.Count; i++)
            {
                if (enterCount[i] == 0)
                {
                    firstUnitig.Add(i);
                }
            }
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
            if (max >= uniteThreshold)
            {
                adjacencyList[index] = [((int, ulong))path!];
            }
            else
            {
                adjacencyList[index] = [];
            }
        }

        private static string WalkPath(List<string> unitigList, List<List<(int, ulong)>> adjacencyList, int index, bool[] visited)
        {
            StringBuilder sb = new StringBuilder(unitigList[index]);
            while (adjacencyList[index].Count > 0 && !visited[index])
            {
                visited[index] = true;
                var next = adjacencyList[index][0].Item1;
                var unitig = unitigList[next];
                var flag = false;
                for (var i = Math.Min(sb.Length, unitig.Length); i > 0; i--)
                {
                    flag = true;
                    var offset = sb.Length - i;
                    for (var j = 0; j < i; j++)
                    {
                        if (sb[offset + j] != unitig[j])
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)
                    {
                        sb.Append(unitig[i..]);
                        break;
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
    }
}
