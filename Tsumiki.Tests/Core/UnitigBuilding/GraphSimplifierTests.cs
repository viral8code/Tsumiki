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
        #region 定数

        /// <summary>
        /// 項目 tsumiki graph simplifier tests
        /// </summary>
        private const string C_項目_tsumiki_graph_simplifier_tests = "tsumiki_graph_simplifier_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 塩基配列 GCTAAAGACAATTACATAACATAC GGCAATTGACCTGAAT
        /// </summary>
        private const string C_塩基配列_GCTAAAGACAATTACATAACATAC_GGCAATTGACCTGAAT = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";

        /// <summary>
        /// 区切り
        /// </summary>
        private const string C_区切り = ",";

        /// <summary>
        /// 塩基配列 GCTAAAGACAATTACATAACATAC
        /// </summary>
        private const string C_塩基配列_GCTAAAGACAATTACATAACATAC = "GCTAAAGACAATTACATAACATAC";

        /// <summary>
        /// 塩基配列 GAAGTTGCCGTACTAAATTA
        /// </summary>
        private const string C_塩基配列_GAAGTTGCCGTACTAAATTA = "GAAGTTGCCGTACTAAATTA";

        /// <summary>
        /// 塩基配列 TGACAGCCGGGGATCTTCCC
        /// </summary>
        private const string C_塩基配列_TGACAGCCGGGGATCTTCCC = "TGACAGCCGGGGATCTTCCC";

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        /// <summary>
        /// 項目 C
        /// </summary>
        private const string C_項目_C = "C";

        /// <summary>
        /// 塩基配列 ATATCACACCCAACCTTCAA
        /// </summary>
        private const string C_塩基配列_ATATCACACCCAACCTTCAA = "ATATCACACCCAACCTTCAA";

        /// <summary>
        /// 塩基配列 ATGCCGTGCCCTAACGCCCT
        /// </summary>
        private const string C_塩基配列_ATGCCGTGCCCTAACGCCCT = "ATGCCGTGCCCTAACGCCCT";

        /// <summary>
        /// 塩基配列 AATCCTGCGCTAGGGGTTGCAGCGACCAGA
        /// </summary>
        private const string C_塩基配列_AATCCTGCGCTAGGGGTTGCAGCGACCAGA = "AATCCTGCGCTAGGGGTTGCAGCGACCAGA";

        /// <summary>
        /// 塩基配列 ATGCCGTGCCCTAACGCCCTAATCCTGCGCTAGGGG
        /// </summary>
        private const string C_塩基配列_ATGCCGTGCCCTAACGCCCTAATCCTGCGCTAGGGG = "ATGCCGTGCCCTAACGCCCTAATCCTGCGCTAGGGG";

        /// <summary>
        /// 塩基配列 GGATCCTTCACGT
        /// </summary>
        private const string C_塩基配列_GGATCCTTCACGT = "GGATCCTTCACGT";

        /// <summary>
        /// 塩基配列 ACATAACAGTCA
        /// </summary>
        private const string C_塩基配列_ACATAACAGTCA = "ACATAACAGTCA";

        /// <summary>
        /// 塩基配列 ACATAACCTTGA
        /// </summary>
        private const string C_塩基配列_ACATAACCTTGA = "ACATAACCTTGA";

        /// <summary>
        /// 塩基配列 CAATTGACGATC
        /// </summary>
        private const string C_塩基配列_CAATTGACGATC = "CAATTGACGATC";

        /// <summary>
        /// 塩基配列 GACAATTTGCA
        /// </summary>
        private const string C_塩基配列_GACAATTTGCA = "GACAATTTGCA";

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public GraphSimplifierTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_graph_simplifier_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_スペクトルモデル = null;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// tip 除去で短い行き止まり分岐を除去し、単一の主 unitig に再構築される
        /// </summary>
        [Fact]
        public void V_tip除去で短い行き止まり分岐を除去し単一の主unitigに再構築される()
        {
            const string l_主配列 = C_塩基配列_GCTAAAGACAATTACATAACATAC_GGCAATTGACCTGAAT;
            const int l_k長 = 8;
            const int l_分岐点 = 20;
            const int l_tipExtra = 4;
            using var l_インデックス = this.V_構築_tip付き索引(l_主配列, l_k長, l_分岐点, l_tipExtra);
            var l_simplifiedFirstKmers = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);
            var l_unitig構築 = new UnitigMaker(l_インデックス);
            HashSet<string> l_seen = [];
            var l_unitig群 = new List<string>();
            foreach (var l_kmer in l_simplifiedFirstKmers)
            {
                var l_u = l_unitig構築.Get_Unitig(l_kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }

                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_unitig群.Add(l_u.A_配列);
            }

            var l_longest = l_unitig群.OrderByDescending(u => u.Length).First();
            Assert.True(l_longest.Length >= l_主配列.Length - l_k長, $"expected a near-full-length main unitig, longest was {l_longest.Length}bp among [{string.Join(C_区切り, l_unitig群.Select(u => u.Length))}]");
        }

        /// <summary>
        /// 行き止まりでも、カバレッジが主経路並みなら除去しないこと
        /// </summary>
        [Fact]
        public void V_主経路並みのカバレッジを持つ短い行き止まり分岐は除去されない()
        {
            const string l_主配列 = C_塩基配列_GCTAAAGACAATTACATAACATAC_GGCAATTGACCTGAAT;
            const int l_k長 = 8;
            const int l_分岐点 = 20;
            const int l_tipExtra = 4;
            using var l_インデックス = this.V_構築_tip付き索引(l_主配列, l_k長, l_分岐点, l_tipExtra, p_主経路反復回数: 10, p_末端反復回数: 10);
            var l_before = l_インデックス.Get_信頼kmer一覧().Count();
            _ = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);
            Assert.Equal(l_before, l_インデックス.Get_信頼kmer一覧().Count());
        }

        /// <summary>
        /// E2 では低カバレッジ端トリミングだけを独立に止められることを確かめる
        /// </summary>
        [Fact]
        public void V_E2は端トリミングを止めてもtip除去の判定を変えない()
        {
            const string l_主配列 = C_塩基配列_GCTAAAGACAATTACATAACATAC_GGCAATTGACCTGAAT;
            const int l_k長 = 8;
            using var l_インデックス = this.V_構築_tip付き索引(l_主配列, l_k長, 20, 4);
            var l_前 = l_インデックス.Get_信頼kmer一覧().Count();
            _ = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2, p_Is低カバレッジ端トリミング: false);
            Assert.True(l_インデックス.Get_信頼kmer一覧().Count() < l_前);
        }

        /// <summary>
        /// 分岐の無い直鎖配列は tip 除去で変化しない
        /// </summary>
        [Fact]
        public void V_分岐の無い直鎖配列はtip除去で変化しない()
        {
            const string l_配列 = C_塩基配列_GCTAAAGACAATTACATAACATAC;
            const int l_k長 = 8;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = 1
            };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_配列), l_k長, 3);
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            var l_before = l_インデックス.Get_信頼kmer一覧().Count();
            var l_firstKmers = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);
            var l_after = l_インデックス.Get_信頼kmer一覧().Count();
            Assert.Equal(l_before, l_after);
            Assert.NotEmpty(l_firstKmers);
        }

        /// <summary>
        /// 低カバレッジの bubble 分岐を除去し、高カバレッジ分岐は残す
        /// </summary>
        [Fact]
        public void V_低カバレッジのbubble分岐を除去し高カバレッジ分岐は残す()
        {
            const string l_commonBefore = C_塩基配列_GAAGTTGCCGTACTAAATTA;
            const string l_sharedAfter = C_塩基配列_TGACAGCCGGGGATCTTCCC;
            const string l_seqHighCoverage = l_commonBefore + C_項目_A + l_sharedAfter;
            const string l_seqLowCoverage = l_commonBefore + C_項目_C + l_sharedAfter;
            const int l_k長 = 8;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = 1
            };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_seqHighCoverage), l_k長, 20);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_seqLowCoverage), l_k長, 3);
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            var l_simplifiedFirstKmers = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);
            var l_unitig構築 = new UnitigMaker(l_インデックス);
            HashSet<string> l_seen = [];
            var l_unitig群 = new List<string>();
            foreach (var l_kmer in l_simplifiedFirstKmers)
            {
                var l_u = l_unitig構築.Get_Unitig(l_kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }

                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_unitig群.Add(l_u.A_配列);
            }

            Assert.DoesNotContain(l_unitig群, u => u.Contains('C' + l_sharedAfter[..(l_k長 - 1)]));
            var l_fullHigh = l_seqHighCoverage;
            var l_fullHighRevComp = Util.V_逆相補(l_fullHigh);
            Assert.Contains(l_unitig群, u => u == l_fullHigh || u == l_fullHighRevComp || u.Contains(l_fullHigh) || u.Contains(l_fullHighRevComp));
        }

        /// <summary>
        /// 2 つの異なる経路が同じ配列へ合流する構造 (reverse bubble) で、合流後の共有配列が複数の unitig に重複して現れないことを確認する
        /// </summary>
        [Fact]
        public void V_合流点でunitig構築が止まり共有配列が重複出力されない()
        {
            const string l_prefixA = C_塩基配列_ATATCACACCCAACCTTCAA;
            const string l_prefixB = C_塩基配列_ATGCCGTGCCCTAACGCCCT;
            const string l_shared = C_塩基配列_AATCCTGCGCTAGGGGTTGCAGCGACCAGA;
            const int l_k長 = 8;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = 1
            };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_prefixA + l_shared), l_k長, 10);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_prefixB + l_shared), l_k長, 10);
            var l_firstKmers = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            var l_unitig構築 = new UnitigMaker(l_インデックス);
            HashSet<string> l_seen = [];
            var l_unitig群 = new List<string>();
            foreach (var l_kmer in l_firstKmers)
            {
                var l_u = l_unitig構築.Get_Unitig(l_kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }

                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_unitig群.Add(l_u.A_配列);
            }

            var l_sharedKmer = l_shared[..l_k長];
            var l_sharedKmerRc = Util.V_逆相補(l_sharedKmer);
            var l_occurrences = l_unitig群.Sum(u => Get_出現回数(u, l_sharedKmer) + Get_出現回数(u, l_sharedKmerRc));
            Assert.Equal(1, l_occurrences);
        }

        /// <summary>
        /// 分岐へ入る行き止まりの短い枝は、競合する枝が無ければ全体の基準値より低くても除去しない
        /// </summary>
        [Fact]
        public void V_競合する枝の無い低カバレッジの行き止まりは除去されない()
        {
            const string l_主配列 = C_塩基配列_GCTAAAGACAATTACATAACATAC_GGCAATTGACCTGAAT;
            const string l_別の続き = C_塩基配列_ATGCCGTGCCCTAACGCCCTAATCCTGCGCTAGGGG;
            const int l_k長 = 8;
            var l_上流 = l_主配列[..12];
            var l_続き = l_主配列[12..];
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = 1
            };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_上流 + l_続き), l_k長, 4);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_上流 + l_別の続き), l_k長, 4);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_続き), l_k長, 36);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_別の続き), l_k長, 36);
            l_インデックス.V_適用_カットオフ(2UL);
            var l_前 = l_インデックス.Get_信頼kmer一覧().Count();
            _ = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);
            Assert.Equal(l_前, l_インデックス.Get_信頼kmer一覧().Count());
        }

        /// <summary>
        /// 後から足した k-mer の順番を変えても、整理の結果は変わらない
        /// </summary>
        [Fact]
        public void V_kmerを足す順番を変えても整理の結果は同じ()
        {
            const string l_主配列 = C_塩基配列_GCTAAAGACAATTACATAACATAC_GGCAATTGACCTGAAT;
            const int l_k長 = 8;
            string[] l_枝群 = [C_塩基配列_GGATCCTTCACGT, C_塩基配列_ACATAACAGTCA, C_塩基配列_ACATAACCTTGA, C_塩基配列_CAATTGACGATC, C_塩基配列_GACAATTTGCA];
            List<string>[] l_結果 = [[], []];
            for (var l_回 = 0; l_回 < 2; l_回++)
            {
                ConfigurationManager.A_実行時引数 = new Parameters
                {
                    A_k長 = l_k長,
                    A_スレッド数 = 4
                };
                var l_場所 = Directory.CreateDirectory(Path.Combine(this._作業ディレクトリ, l_回.ToString())).FullName;
                using var l_インデックス = new TrustedKmerIndex(l_場所);
                V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_主配列), l_k長, 10);
                _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
                List<(byte[] A_kmer, ulong A_出現回数)> l_足すkmer = [];
                for (var i = 0; i < l_枝群.Length; i++)
                {
                    var l_塩基列 = V_変換_塩基ID列(l_枝群[i]);
                    for (var j = 0; j + l_k長 <= l_塩基列.Length; j++)
                    {
                        l_足すkmer.Add((l_塩基列[j..(j + l_k長)], (ulong)(2 + i)));
                    }
                }

                if (l_回 == 1)
                {
                    l_足すkmer.Reverse();
                }

                foreach (var (l_kmer, l_出現回数)in l_足すkmer)
                {
                    _ = l_インデックス.Is成功_追加_信頼kmer(l_kmer, l_出現回数);
                }

                _ = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);
                l_結果[l_回] = [..l_インデックス.Get_信頼kmer一覧().Select(x => string.Concat(x.Select(Util.V_変換_塩基文字))).Order(StringComparer.Ordinal)];
            }

            Assert.Equal(l_結果[0], l_結果[1]);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] V_変換_塩基ID列(string p_配列)
        {
            return[..p_配列.Select(Util.Get_塩基ID)];
        }

        /// <summary>
        /// 配列の全 k-mer を、指定した回数だけ信頼できる kmer 索引へ登録する
        /// </summary>
        /// <param name="p_インデックス">登録先の索引</param>
        /// <param name="p_バイト列">登録する配列の塩基 ID 列</param>
        /// <param name="p_k長">k-mer 長</param>
        /// <param name="p_反復回数">各 k-mer を登録する回数</param>
        private static void V_登録_全kmer(TrustedKmerIndex p_インデックス, byte[] p_バイト列, int p_k長, int p_反復回数)
        {
            for (var i = 0; i + p_k長 <= p_バイト列.Length; i++)
            {
                for (var l_反復回数 = 0; l_反復回数 < p_反復回数; l_反復回数++)
                {
                    p_インデックス.V_登録(p_バイト列.AsSpan(i, p_k長));
                }
            }
        }

        /// <summary>
        /// mainSeq (主経路) の k-mer 群に加えて、その途中の 1 点から分岐する短い tip 配列 (tipSeq) の k-mer 群も登録した TrustedKmerIndex を作る
        /// </summary>
        /// <param name="p_主配列">主経路の配列</param>
        /// <param name="p_kmer長">k-mer 長</param>
        /// <param name="p_分岐点">分岐させる主経路上の位置</param>
        /// <param name="p_末端長">分岐後に伸ばす長さ</param>
        /// <param name="p_主経路反復回数">主経路の各 k-mer を登録する回数</param>
        /// <param name="p_末端反復回数">tip の各 k-mer を登録する回数</param>
        /// <returns>構築した索引</returns>
        private TrustedKmerIndex V_構築_tip付き索引(string p_主配列, int p_kmer長, int p_分岐点, int p_末端長, int p_主経路反復回数 = 10, int p_末端反復回数 = 2)
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = p_kmer長,
                A_スレッド数 = 1
            };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(p_主配列), p_kmer長, p_主経路反復回数);
            var l_overlap = p_主配列.Substring(p_分岐点, p_kmer長 - 1);
            var l_branchBaseChar = p_主配列[p_分岐点 + p_kmer長 - 1];
            var l_altChar = Consts.塩基文字.First(c => c != l_branchBaseChar);
            var l_tipSeq = l_overlap + l_altChar + string.Concat(Enumerable.Range(0, p_末端長).Select(i => Consts.塩基文字[i % 4]));
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_tipSeq), p_kmer長, p_末端反復回数);
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        /// <summary>
        /// 部分文字列が現れる回数を数える
        /// </summary>
        /// <param name="p_検索対象">探される側</param>
        /// <param name="p_検索配列">探す文字列</param>
        /// <returns>現れた回数</returns>
        private static int Get_出現回数(string p_検索対象, string p_検索配列)
        {
            var l_count = 0;
            for (var i = 0; i + p_検索配列.Length <= p_検索対象.Length; i++)
            {
                if (string.CompareOrdinal(p_検索対象, i, p_検索配列, 0, p_検索配列.Length) == 0)
                {
                    l_count++;
                }
            }

            return l_count;
        }

        #endregion
    }
}
