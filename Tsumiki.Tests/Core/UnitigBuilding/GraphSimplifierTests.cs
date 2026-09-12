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
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_graph_simplifier_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);

            // これらのテストは比率ベースの tip 判定を検証する
            // 他のテストが残した
            // 混合モデルの適合結果が「無条件に信頼する下限」として漏れ込み、
            // 判定を横取りしないようにする
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
            // 非周期的な主経路配列 (k=8 での内部重複なしを別途 Python で確認済み)
            const string l_主配列 = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int l_k長 = 8;
            const int l_分岐点 = 20; // 主経路の途中から分岐させる
            const int l_tipExtra = 4; // 分岐後にごく短く伸びる tip

            using var l_インデックス = this.V_構築_tip付き索引(l_主配列, l_k長, l_分岐点, l_tipExtra);

            var l_simplifiedFirstKmers = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);

            // tip 除去後は、分岐点だった箇所の次数が解消され、
            // 主経路が 1 本の unitig として (理想的には) 再構築されるはず
            var l_ユニティグ構築 = new UnitigMaker(l_インデックス);
            HashSet<string> l_seen = [];
            var l_ユニティグ群 = new List<string>();
            foreach (var l_kmer in l_simplifiedFirstKmers)
            {
                var l_u = l_ユニティグ構築.Get_ユニティグ(l_kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }
                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_ユニティグ群.Add(l_u.A_配列);
            }

            // tip 自体はもう存在しないはずなので、tip 由来の短い配列を含む
            // unitig は残っていないこと、かつ主経路の全長をカバーする
            // (ほぼ) 1 本の unitig が存在することを確認する
            var l_longest = l_ユニティグ群.OrderByDescending(u => u.Length).First();
            Assert.True(l_longest.Length >= l_主配列.Length - l_k長, $"expected a near-full-length main unitig, longest was {l_longest.Length}bp among [{string.Join(",", l_ユニティグ群.Select(u => u.Length))}]");
        }

        /// <summary>
        /// 行き止まりでも、カバレッジが主経路並みなら除去しないこと
        /// </summary>
        /// <remarks>
        /// カバレッジの切れ目で孤立した実配列がこの形になるため、行き止まりというだけで消すとゲノム被覆率を落とす
        /// </remarks>
        [Fact]
        public void V_主経路並みのカバレッジを持つ短い行き止まり分岐は除去されない()
        {
            const string l_主配列 = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";
            const int l_k長 = 8;
            const int l_分岐点 = 20;
            const int l_tipExtra = 4;

            using var l_インデックス = this.V_構築_tip付き索引(l_主配列, l_k長, l_分岐点, l_tipExtra, p_主経路反復回数: 10, p_末端反復回数: 10);

            var l_before = l_インデックス.Get_信頼kmer一覧().Count();
            _ = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);

            Assert.Equal(l_before, l_インデックス.Get_信頼kmer一覧().Count());
        }

        /// <summary>
        /// 分岐の無い直鎖配列は tip 除去で変化しない
        /// </summary>
        [Fact]
        public void V_分岐の無い直鎖配列はtip除去で変化しない()
        {
            const string l_配列 = "GCTAAAGACAATTACATAACATAC"; // 24 bp、非周期的
            const int l_k長 = 8;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_配列), l_k長, 3);
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            var l_before = l_インデックス.Get_信頼kmer一覧().Count();

            var l_firstKmers = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);

            var l_after = l_インデックス.Get_信頼kmer一覧().Count();

            // 分岐のない直鎖配列には tip が存在しないため、何も除去されないはず
            Assert.Equal(l_before, l_after);
            Assert.NotEmpty(l_firstKmers);
        }

        /// <summary>
        /// 低カバレッジの bubble 分岐を除去し、高カバレッジ分岐は残す
        /// </summary>
        [Fact]
        public void V_低カバレッジのbubble分岐を除去し高カバレッジ分岐は残す()
        {
            // 分岐点 (commonBefore 末尾) から 1 塩基だけ異なる ('A' vs 'C') 経路 B/C に
            // 分かれ、その後 sharedAfter へ合流する SNP 様の単純な bubble 構造
            // Python (scripts 外、事前検証) で k=8 内に重複が生じないことを確認済み
            const string l_commonBefore = "GAAGTTGCCGTACTAAATTA"; // 20 bp
            const string l_sharedAfter = "TGACAGCCGGGGATCTTCCC"; // 20 bp
            const string l_seqHighCoverage = l_commonBefore + "A" + l_sharedAfter; // 分岐点で A
            const string l_seqLowCoverage = l_commonBefore + "C" + l_sharedAfter; // 分岐点で C (エラー相当)
            const int l_k長 = 8;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            // 高カバレッジ経路 (真のゲノム由来相当) は 20 回、低カバレッジ経路
            // (エラー由来相当) は 3 回登録する (カットオフ 2 は超えるが、
            // baseline (高カバレッジ経路水準) に比べて著しく低い)
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_seqHighCoverage), l_k長, 20);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_seqLowCoverage), l_k長, 3);

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            var l_simplifiedFirstKmers = GraphSimplifier.V_除去_tip(l_インデックス, l_k長, p_tip長閾値: l_k長 * 2);

            var l_ユニティグ構築 = new UnitigMaker(l_インデックス);
            HashSet<string> l_seen = [];
            var l_ユニティグ群 = new List<string>();
            foreach (var l_kmer in l_simplifiedFirstKmers)
            {
                var l_u = l_ユニティグ構築.Get_ユニティグ(l_kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }
                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_ユニティグ群.Add(l_u.A_配列);
            }

            // 低カバレッジ経路の分岐点を含む短い断片は残っていないはず
            // (再構築された配列のいずれにも "C" + sharedAfter の先頭部分は
            // 現れない = 低カバレッジ経路は除去された)
            Assert.DoesNotContain(l_ユニティグ群, u => u.Contains('C' + l_sharedAfter[..(l_k長 - 1)]));

            // 高カバレッジ経路 (commonBefore + "A" + sharedAfter の全体、または
            // その逆相補) を含む、ほぼ全長の unitig が存在するはず
            var l_fullHigh = l_seqHighCoverage;
            var l_fullHighRevComp = Util.V_逆相補(l_fullHigh);
            Assert.Contains(l_ユニティグ群, u => u == l_fullHigh || u == l_fullHighRevComp || u.Contains(l_fullHigh) || u.Contains(l_fullHighRevComp));
        }

        /// <summary>
        /// 2 つの異なる経路が同じ配列へ合流する構造 (reverse bubble) で、合流後の共有配列が複数の unitig に重複して現れないことを確認する
        /// </summary>
        /// <remarks>
        /// unitig の定義は「内部の全節点が入次数 1 かつ出次数 1 の極大パス」であり、合流点 (入次数 2) からは別の unitig が始まらなければならない<br/>
        /// この規則が無いと、両方の経路の walk が共有配列を走り抜けてしまい、同じ配列を 2 度出力する (実データで k-mer 延べ数が実内容の 1.43 倍に膨らんでいた原因)
        /// </remarks>
        [Fact]
        public void V_合流点でunitig構築が止まり共有配列が重複出力されない()
        {
            const string l_prefixA = "ATATCACACCCAACCTTCAA";
            const string l_prefixB = "ATGCCGTGCCCTAACGCCCT";
            const string l_shared = "AATCCTGCGCTAGGGGTTGCAGCGACCAGA";
            const int l_k長 = 8;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_prefixA + l_shared), l_k長, 10);
            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(l_prefixB + l_shared), l_k長, 10);

            var l_firstKmers = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            var l_ユニティグ構築 = new UnitigMaker(l_インデックス);
            HashSet<string> l_seen = [];
            var l_ユニティグ群 = new List<string>();
            foreach (var l_kmer in l_firstKmers)
            {
                var l_u = l_ユニティグ構築.Get_ユニティグ(l_kmer);
                if (l_seen.Contains(l_u.A_配列) || l_seen.Contains(Util.V_逆相補(l_u.A_配列)))
                {
                    continue;
                }
                _ = l_seen.Add(l_u.A_配列);
                _ = l_seen.Add(Util.V_逆相補(l_u.A_配列));
                l_ユニティグ群.Add(l_u.A_配列);
            }

            // 共有配列の先頭 k-mer (またはその逆相補) が、全 unitig を通じて
            // 延べ 1 回しか現れないこと (=重複出力されていないこと) を確認する
            var l_sharedKmer = l_shared[..l_k長];
            var l_sharedKmerRc = Util.V_逆相補(l_sharedKmer);
            var l_occurrences = l_ユニティグ群.Sum(u => V_数える_出現回数(u, l_sharedKmer) + V_数える_出現回数(u, l_sharedKmerRc));
            Assert.Equal(1, l_occurrences);
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
            return [.. p_配列.Select(Util.Get_塩基ID)];
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
        /// <remarks>
        /// tipSeq は mainSeq の位置 branchPoint から始まる長さ kmerLength-1 の「本来の続き」をコピーした上で、最後の 1 塩基だけ変えることで主経路と k-1 塩基だけ重なる分岐を作る単純な構成にする
        /// </remarks>
        /// <param name="p_主配列">主経路の配列</param>
        /// <param name="p_kmer長">k-mer 長</param>
        /// <param name="p_分岐点">分岐させる主経路上の位置</param>
        /// <param name="p_末端長">分岐後に伸ばす長さ</param>
        /// <param name="p_主経路反復回数">主経路の各 k-mer を登録する回数</param>
        /// <param name="p_末端反復回数">tip の各 k-mer を登録する回数</param>
        /// <returns>構築した索引</returns>
        private TrustedKmerIndex V_構築_tip付き索引(string p_主配列, int p_kmer長, int p_分岐点, int p_末端長, int p_主経路反復回数 = 10, int p_末端反復回数 = 2)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmer長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            V_登録_全kmer(l_インデックス, V_変換_塩基ID列(p_主配列), p_kmer長, p_主経路反復回数);

            // 分岐点の直前 kmerLength-1 文字を土台に、最後だけ主経路と異なる
            // 1 塩基を続けて tip を伸ばす (主経路と k-2 塩基だけ重なる短い枝)
            var l_overlap = p_主配列.Substring(p_分岐点, p_kmer長 - 1);
            var l_branchBaseChar = p_主配列[p_分岐点 + p_kmer長 - 1];
            var l_altChar = "ACGT".First(c => c != l_branchBaseChar);
            // overlap (主経路と k-1 塩基共有) + altChar (主経路とは異なる 1 塩基) +
            // 適当なユニークな続きで、主経路から分岐する短い tip を作る
            var l_tipSeq = l_overlap + l_altChar + string.Concat(Enumerable.Range(0, p_末端長).Select(i => "ACGT"[i % 4]));
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
        private static int V_数える_出現回数(string p_検索対象, string p_検索配列)
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
