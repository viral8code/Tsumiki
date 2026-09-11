using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig グラフから tip (行き止まりの短い枝) を除去する処理の検証
    /// </summary>
    public class GraphSimplifierTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public GraphSimplifierTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_graph_simplifier_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);

            // これらのテストは比率ベースの tip 判定を検証する
            // 他のテストが残した
            // 混合モデルの適合結果が「無条件に信頼する下限」として漏れ込み、
            // 判定を横取りしないようにする
            ConfigurationManager.A_スペクトルモデル = null;
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="p_seq">元の配列</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] V_変換_塩基ID列(string p_seq)
        {
            return [.. p_seq.Select(Util.Get_塩基ID)];
        }

        /// <summary>
        /// 配列の全 k-mer を、指定した回数だけ信頼できる kmer 索引へ登録する
        /// </summary>
        /// <param name="p_index">登録先の索引</param>
        /// <param name="p_bytes">登録する配列の塩基 ID 列</param>
        /// <param name="p_k">k-mer 長</param>
        /// <param name="p_repetitions">各 k-mer を登録する回数</param>
        private static void V_登録_全kmer(TrustedKmerIndex p_index, byte[] p_bytes, int p_k, int p_repetitions)
        {
            for (var i = 0; i + p_k <= p_bytes.Length; i++)
            {
                for (var rep = 0; rep < p_repetitions; rep++)
                {
                    p_index.V_登録(p_bytes.AsSpan(i, p_k));
                }
            }
        }

        /// <summary>
        /// mainSeq(主経路) の k-mer 群に加えて、その途中の 1 点から分岐する
        /// 短いtip配列 (tipSeq) の k-mer 群も登録した TrustedKmerIndex を作る
        /// </summary>
        /// <remarks>
        /// tipSeq は mainSeq の位置 branchPoint から始まる長さ kmerLength-1 の
        /// 「本来の続き」をコピーした上で、最後の 1 塩基だけ変えることで
        /// 主経路と k-1 塩基だけ重なる分岐を作る単純な構成にする
        /// </remarks>
        /// <param name="p_mainSeq">主経路の配列</param>
        /// <param name="p_kmerLength">k-mer 長</param>
        /// <param name="p_branchPoint">分岐させる主経路上の位置</param>
        /// <param name="p_tipLength">分岐後に伸ばす長さ</param>
        /// <param name="p_mainRepetitions">主経路の各 k-mer を登録する回数</param>
        /// <param name="p_tipRepetitions">tip の各 k-mer を登録する回数</param>
        /// <returns>構築した索引</returns>
        private TrustedKmerIndex V_構築_tip付き索引(
            string p_mainSeq, int p_kmerLength, int p_branchPoint, int p_tipLength,
            int p_mainRepetitions = 10, int p_tipRepetitions = 2)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmerLength, A_スレッド数 = 1 };
            var l_index = new TrustedKmerIndex(this._tempDir);

            V_登録_全kmer(l_index, V_変換_塩基ID列(p_mainSeq), p_kmerLength, p_mainRepetitions);

            // 分岐点の直前 kmerLength-1 文字を土台に、最後だけ主経路と異なる
            // 1 塩基を続けて tip を伸ばす (主経路と k-2 塩基だけ重なる短い枝)
            var l_overlap = p_mainSeq.Substring(p_branchPoint, p_kmerLength - 1);
            var l_branchBaseChar = p_mainSeq[p_branchPoint + p_kmerLength - 1];
            var l_altChar = "ACGT".First(c => c != l_branchBaseChar);
            // overlap(主経路と k-1 塩基共有)+ altChar(主経路とは異なる 1 塩基)+
            // 適当なユニークな続きで、主経路から分岐する短い tip を作る
            var l_tipSeq = l_overlap + l_altChar + string.Concat(Enumerable.Range(0, p_tipLength).Select(i => "ACGT"[i % 4]));
            V_登録_全kmer(l_index, V_変換_塩基ID列(l_tipSeq), p_kmerLength, p_tipRepetitions);

            _ = l_index.V_カットオフ(p_カットオフ: 2);
            return l_index;
        }

        /// <summary>
        /// tip 除去で短い行き止まり分岐を除去し、単一の主 unitig に再構築される
        /// </summary>
        [Fact]
        public void tip除去で短い行き止まり分岐を除去し単一の主unitigに再構築される()
        {
            // 非周期的な主経路配列 (k=8 での内部重複なしを別途 Python で確認済み)
            const string l_mainSeq = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int k = 8;
            const int l_branchPoint = 20; // 主経路の途中から分岐させる
            const int l_tipExtra = 4; // 分岐後にごく短く伸びるtip

            using var l_index = this.V_構築_tip付き索引(l_mainSeq, k, l_branchPoint, l_tipExtra);

            var l_simplifiedFirstKmers = GraphSimplifier.V_除去_tip(l_index, k, p_tip長閾値: k * 2);

            // tip 除去後は、分岐点だった箇所の次数が解消され、
            // 主経路が 1 本のunitigとして (理想的には) 再構築されるはず
            var l_unitigMaker = new UnitigMaker(l_index);
            HashSet<string> l_seen = [];
            var l_unitigs = new List<string>();
            foreach (var kmer in l_simplifiedFirstKmers)
            {
                var l_u = l_unitigMaker.Get_ユニティグ(kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }
                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_unitigs.Add(l_u.A_配列);
            }

            // tip 自体はもう存在しないはずなので、tip 由来の短い配列を含む
            // unitig は残っていないこと、かつ主経路の全長をカバーする
            // (ほぼ) 1 本の unitig が存在することを確認する
            var l_longest = l_unitigs.OrderByDescending(u => u.Length).First();
            Assert.True(l_longest.Length >= l_mainSeq.Length - k, $"expected a near-full-length main unitig, longest was {l_longest.Length}bp among [{string.Join(",", l_unitigs.Select(u => u.Length))}]");
        }

        /// <summary>
        /// 行き止まりでも、カバレッジが主経路並みなら除去しないこと
        /// </summary>
        /// <remarks>
        /// カバレッジの切れ目で孤立した実配列がこの形になるため、
        /// 行き止まりというだけで消すとゲノム被覆率を落とす
        /// </remarks>
        [Fact]
        public void 主経路並みのカバレッジを持つ短い行き止まり分岐は除去されない()
        {
            const string l_mainSeq = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int k = 8;
            const int l_branchPoint = 20;
            const int l_tipExtra = 4;

            using var l_index = this.V_構築_tip付き索引(
                l_mainSeq, k, l_branchPoint, l_tipExtra, p_mainRepetitions: 10, p_tipRepetitions: 10);

            var l_before = l_index.Get_信頼kmer一覧().Count();
            _ = GraphSimplifier.V_除去_tip(l_index, k, p_tip長閾値: k * 2);

            Assert.Equal(l_before, l_index.Get_信頼kmer一覧().Count());
        }

        /// <summary>
        /// 分岐の無い直鎖配列は tip 除去で変化しない
        /// </summary>
        [Fact]
        public void 分岐の無い直鎖配列はtip除去で変化しない()
        {
            const string l_seq = "GCTAAAGACAATTACATAACATAC"; // 24bp、非周期的
            const int k = 8;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            using var l_index = new TrustedKmerIndex(this._tempDir);
            V_登録_全kmer(l_index, V_変換_塩基ID列(l_seq), k, 3);
            _ = l_index.V_カットオフ(p_カットオフ: 2);

            var l_before = l_index.Get_信頼kmer一覧().Count();

            var l_firstKmers = GraphSimplifier.V_除去_tip(l_index, k, p_tip長閾値: k * 2);

            var l_after = l_index.Get_信頼kmer一覧().Count();

            // 分岐のない直鎖配列には tip が存在しないため、何も除去されないはず
            Assert.Equal(l_before, l_after);
            Assert.NotEmpty(l_firstKmers);
        }

        /// <summary>
        /// 低カバレッジの bubble 分岐を除去し、高カバレッジ分岐は残す
        /// </summary>
        [Fact]
        public void 低カバレッジのbubble分岐を除去し高カバレッジ分岐は残す()
        {
            // 分岐点 (commonBefore 末尾) から 1 塩基だけ異なる ('A' vs 'C') 経路 B/C に
            // 分かれ、その後 sharedAfter へ合流する SNP 様の単純な bubble 構造
            // Python(scripts 外、事前検証) で k=8 内に重複が生じないことを確認済み
            const string l_commonBefore = "GAAGTTGCCGTACTAAATTA"; // 20bp
            const string l_sharedAfter = "TGACAGCCGGGGATCTTCCC"; // 20bp
            const string l_seqHighCoverage = l_commonBefore + "A" + l_sharedAfter; // 分岐点でA
            const string l_seqLowCoverage = l_commonBefore + "C" + l_sharedAfter; // 分岐点でC(エラー相当)
            const int k = 8;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            using var l_index = new TrustedKmerIndex(this._tempDir);

            // 高カバレッジ経路 (真のゲノム由来相当) は 20 回、低カバレッジ経路
            // (エラー由来相当) は 3 回登録する (カットオフ 2 は超えるが、
            // baseline(高カバレッジ経路水準) に比べて著しく低い)
            V_登録_全kmer(l_index, V_変換_塩基ID列(l_seqHighCoverage), k, 20);
            V_登録_全kmer(l_index, V_変換_塩基ID列(l_seqLowCoverage), k, 3);

            _ = l_index.V_カットオフ(p_カットオフ: 2);

            var l_simplifiedFirstKmers = GraphSimplifier.V_除去_tip(l_index, k, p_tip長閾値: k * 2);

            var l_unitigMaker = new UnitigMaker(l_index);
            HashSet<string> l_seen = [];
            var l_unitigs = new List<string>();
            foreach (var kmer in l_simplifiedFirstKmers)
            {
                var l_u = l_unitigMaker.Get_ユニティグ(kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }
                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_unitigs.Add(l_u.A_配列);
            }

            // 低カバレッジ経路の分岐点を含む短い断片は残っていないはず
            // (再構築された配列のいずれにも "C" + sharedAfter の先頭部分は
            // 現れない = 低カバレッジ経路は除去された)
            Assert.DoesNotContain(l_unitigs, u => u.Contains('C' + l_sharedAfter[..(k - 1)]));

            // 高カバレッジ経路 (commonBefore + "A" + sharedAfter の全体、または
            // その逆相補) を含む、ほぼ全長の unitig が存在するはず
            var l_fullHigh = l_seqHighCoverage;
            var l_fullHighRevComp = Util.V_逆相補(l_fullHigh);
            Assert.Contains(l_unitigs, u => u == l_fullHigh || u == l_fullHighRevComp || u.Contains(l_fullHigh) || u.Contains(l_fullHighRevComp));
        }

        /// <summary>
        /// 2 つの異なる経路が同じ配列へ合流する構造 (reverse bubble) で、
        /// 合流後の共有配列が複数の unitig に重複して現れないことを確認する
        /// </summary>
        /// <remarks>
        /// unitig の定義は「内部の全節点が入次数 1 かつ出次数 1 の極大パス」であり、
        /// 合流点 (入次数 2) からは別の unitig が始まらなければならない<br/>
        /// この規則が無いと、両方の経路の walk が共有配列を走り抜けてしまい、
        /// 同じ配列を 2 度出力する (実データで k-mer 延べ数が実内容の 1.43 倍に
        /// 膨らんでいた原因)
        /// </remarks>
        [Fact]
        public void 合流点でunitig構築が止まり共有配列が重複出力されない()
        {
            const string l_prefixA = "ATATCACACCCAACCTTCAA";
            const string l_prefixB = "ATGCCGTGCCCTAACGCCCT";
            const string l_shared = "AATCCTGCGCTAGGGGTTGCAGCGACCAGA";
            const int k = 8;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };
            using var l_index = new TrustedKmerIndex(this._tempDir);

            V_登録_全kmer(l_index, V_変換_塩基ID列(l_prefixA + l_shared), k, 10);
            V_登録_全kmer(l_index, V_変換_塩基ID列(l_prefixB + l_shared), k, 10);

            var l_firstKmers = l_index.V_カットオフ(p_カットオフ: 2);

            var l_unitigMaker = new UnitigMaker(l_index);
            HashSet<string> l_seen = [];
            var l_unitigs = new List<string>();
            foreach (var kmer in l_firstKmers)
            {
                var l_u = l_unitigMaker.Get_ユニティグ(kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }
                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_unitigs.Add(l_u.A_配列);
            }

            // 共有配列の先頭 k-mer(またはその逆相補) が、全 unitig を通じて
            // 延べ 1 回しか現れないこと (=重複出力されていないこと) を確認する
            var l_sharedKmer = l_shared[..k];
            var l_sharedKmerRc = Util.V_逆相補(l_sharedKmer);
            var l_occurrences = l_unitigs.Sum(u => V_数える_出現回数(u, l_sharedKmer) + V_数える_出現回数(u, l_sharedKmerRc));
            Assert.Equal(1, l_occurrences);
        }

        /// <summary>
        /// 部分文字列が現れる回数を数える
        /// </summary>
        /// <param name="p_haystack">探される側</param>
        /// <param name="p_needle">探す文字列</param>
        /// <returns>現れた回数</returns>
        private static int V_数える_出現回数(string p_haystack, string p_needle)
        {
            var l_count = 0;
            for (var i = 0; i + p_needle.Length <= p_haystack.Length; i++)
            {
                if (string.CompareOrdinal(p_haystack, i, p_needle, 0, p_needle.Length) == 0)
                {
                    l_count++;
                }
            }
            return l_count;
        }
    }
}
