using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// unitig 間の隣接関係を、リードマッピングからの推測ではなく
    /// de Bruijn グラフそのものから厳密に構築する。
    ///
    /// これを導入した経緯(実データでの計測):
    /// 従来の UniteContigs は、リード上で連続して観測された unitig ペア
    /// (kmerPath)を隣接候補とし、結合時に「k-1 塩基のオーバーラップ」を
    /// 試し、一致しなければ長い方から順に任意長のオーバーラップを探す
    /// フォールバックを持っていた。しかし実データでは k-1(=30)での一致が
    /// ほぼ起こらず、フォールバックが平均 2.96 塩基という偶然の一致で
    /// unitig を接着していた(2135 箇所)。つまり結合のほぼ全てが誤結合で
    /// あり、かつ contig 総長は unitig 総長のちょうど 2.009 倍に膨れていた
    /// (順鎖・逆鎖の両方が別々の contig として出力されていた)。
    ///
    /// 隣接の唯一の正しい根拠は de Bruijn グラフの辺であり、それは
    /// 「unitig A の末尾 k-mer から 1 塩基伸ばした k-mer が unitig B の
    /// 先頭 k-mer に一致する」ことと同値である。この条件を満たす辺は
    /// 定義上ちょうど k-1 塩基のオーバーラップを持つため、結合時に
    /// オーバーラップ長を探索する必要がそもそも無くなる。
    ///
    /// 頂点は「符号付き向き」を持つ: unitig ID u に対し
    /// 頂点 2u(順鎖)と 2u+1(逆鎖)。ある頂点 v の逆鎖側の双子は v^1。
    /// 本クラスの構築方法により、辺 v→w が存在すれば必ず w^1→v^1 も
    /// 存在する(逆相補を取れば同じ重なりが成立するため)。
    /// </summary>
    internal sealed class UnitigGraph
    {
        /// <summary>頂点数。unitigList と同じ長さ(添字 0,1 は未使用のダミー)。</summary>
        public int VertexCount => this.OutEdges.Count;

        /// <summary>頂点ごとの出辺(行き先の頂点インデックス)。</summary>
        public List<List<int>> OutEdges { get; }

        private UnitigGraph(List<List<int>> outEdges)
        {
            this.OutEdges = outEdges;
        }

        /// <summary>頂点 v の入次数。辺の逆鎖対称性より、v の入次数は v^1 の出次数に等しい。</summary>
        public int InDegree(int vertex)
        {
            return this.OutEdges[vertex ^ 1].Count;
        }

        /// <summary>
        /// unitigList(添字 2u=順鎖, 2u+1=逆鎖の配列)と、k-mer から
        /// (符号付き unitig ID, その向きでの開始位置) への辞書から、
        /// 厳密な隣接グラフを構築する。
        ///
        /// 「先頭 k-mer である(position==0)」ことを要求するのが要点で、
        /// これにより結合が必ず k-1 オーバーラップの単純連結になる。
        /// 複数 unitig に跨る曖昧 k-mer(ambiguousKmerSentinel)は
        /// 行き先を一意に決められないため辺を張らない。
        /// </summary>
        public static UnitigGraph Build(
            List<string> unitigList,
            IReadOnlyDictionary<KmerKey, (int UnitigId, int Position)> kmerDict,
            int kmerLength,
            int ambiguousKmerSentinel)
        {
            List<List<int>> outEdges = [];
            for (var i = 0; i < unitigList.Count; i++)
            {
                outEdges.Add([]);
            }

            // 末尾 k-mer から 1 塩基伸ばした候補を組み立てるための作業バッファ。
            var candidate = new byte[kmerLength];

            for (var vertex = 2; vertex < unitigList.Count; vertex++)
            {
                var seq = unitigList[vertex];
                if (seq.Length < kmerLength)
                {
                    continue;
                }

                // 末尾 k-mer の 2 文字目以降(k-1 塩基)を候補の先頭に置く。
                var tailStart = seq.Length - kmerLength + 1;
                var hasInvalidBase = false;
                for (var i = 0; i < kmerLength - 1; i++)
                {
                    var id = Util.GetSimpleNucleotideID(seq[tailStart + i]);
                    if (id is < Consts.NucleotideID.A or > Consts.NucleotideID.T)
                    {
                        hasInvalidBase = true;
                        break;
                    }
                    candidate[i] = id;
                }
                if (hasInvalidBase)
                {
                    continue;
                }

                for (byte last = Consts.NucleotideID.A; last <= Consts.NucleotideID.T; last++)
                {
                    candidate[kmerLength - 1] = last;
                    if (!kmerDict.TryGetValue(new KmerKey(candidate.AsSpan()), out var hit))
                    {
                        continue;
                    }
                    if (hit.UnitigId == ambiguousKmerSentinel || hit.Position != 0)
                    {
                        // Position != 0 は「その k-mer が unitig の途中に現れる」
                        // ことを意味し、そこへ k-1 オーバーラップで連結することは
                        // できない(unitig 分割が正しければ本来起きないが、
                        // グラフ簡略化で k-mer を削った結果として起こりうる)。
                        continue;
                    }
                    var target = ContigMaker.VertexIndex(hit.UnitigId);
                    if (target == vertex)
                    {
                        // 自己ループは辿ると無限に伸びるため辺として持たない。
                        continue;
                    }
                    outEdges[vertex].Add(target);
                }
            }

            return new UnitigGraph(outEdges);
        }
    }
}
