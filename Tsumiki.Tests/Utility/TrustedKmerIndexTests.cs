using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 信頼できる k-mer 集合 (TrustedKmerIndex) の登録・カットオフ・問い合わせの検証
    /// </summary>
    public class TrustedKmerIndexTests : IDisposable
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
        public TrustedKmerIndexTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_trusted_kmer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 大量の曖昧候補を展開せず、観測済み k-mer だけを残す
        /// </summary>
        [Fact]
        public void V_曖昧展開_上限超過を除外()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_候補 = Enumerable.Range(0, 8).Select(_ => new byte[] { 1, 2, 3, 4 }).ToArray();
            l_索引.V_登録_曖昧塩基あり(l_候補, 0);
            l_索引.V_登録(V_変換_塩基ID列("ACGTACGT"));
            _ = l_索引.V_カットオフ(1UL);
            Assert.Single(l_索引.Get_信頼kmer一覧());
        }

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
        /// 挿入した k-mer はどちらの向きで問い合わせても含まれ、未挿入の k-mer は含まれない
        /// </summary>
        [Fact]
        public void V_挿入したkmerはどちらの向きでも含まれ未挿入のkmerは含まれない()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_inserted = V_変換_塩基ID列("ACGTACGT");
            var l_revComp = V_変換_塩基ID列(Util.V_逆相補("ACGTACGT"));
            var l_neverInserted = V_変換_塩基ID列("TTTTTTTT");

            // カットオフ (2) を超えるよう複数回登録する
            for (var i = 0; i < 5; i++)
            {
                l_インデックス.V_登録(l_inserted.AsSpan());
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Assert.True(l_インデックス.Get_含まれるか(l_inserted));
            Assert.True(l_インデックス.Get_含まれるか(l_revComp)); // 正規化されるため逆鎖側からの問い合わせでもヒットする
            Assert.False(l_インデックス.Get_含まれるか(l_neverInserted));
        }

        /// <summary>
        /// カットオフ未満の k-mer は除外される
        /// </summary>
        [Fact]
        public void V_カットオフ未満のkmerは除外される()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_belowThreshold = V_変換_塩基ID列("GGGGCCCC");

            l_インデックス.V_登録(l_belowThreshold.AsSpan()); // 1 回だけ = カットオフ 2 未満

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Assert.False(l_インデックス.Get_含まれるか(l_belowThreshold));
        }

        /// <summary>
        /// カバレッジは、順鎖と逆鎖のカウントを合算する
        /// </summary>
        [Fact]
        public void V_カバレッジは順鎖と逆鎖のカウントを合算する()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_forward = V_変換_塩基ID列("ACGTACGT");
            var l_revComp = V_変換_塩基ID列(Util.V_逆相補("ACGTACGT"));

            // 順鎖を 3 回、逆鎖を 2 回登録する
            // カウント段階では別キー扱いだが、
            // カットオフ後の正規化されたエントリでは合算されているはず
            for (var i = 0; i < 3; i++)
            {
                l_インデックス.V_登録(l_forward.AsSpan());
            }
            for (var i = 0; i < 2; i++)
            {
                l_インデックス.V_登録(l_revComp.AsSpan());
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Assert.Equal(5UL, l_インデックス.Get_カバレッジ(l_forward));
            Assert.Equal(5UL, l_インデックス.Get_カバレッジ(l_revComp)); // 正規化されるため同じ値
        }

        /// <summary>
        /// 存在しない k-mer のカバレッジは 0 を返す
        /// </summary>
        [Fact]
        public void V_存在しないkmerのカバレッジは0を返す()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            l_インデックス.V_登録(V_変換_塩基ID列("AAAAAAAA").AsSpan());
            l_インデックス.V_登録(V_変換_塩基ID列("AAAAAAAA").AsSpan());

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Assert.Equal(0UL, l_インデックス.Get_カバレッジ(V_変換_塩基ID列("TTTTGGGG")));
        }

        /// <summary>
        /// 長い k での正規化とカバレッジ合算が正しく行われること
        /// </summary>
        /// <remarks>
        /// k が 32 を超え 64 以下のときの UInt128 経路と、64 を超えたときの KmerKey へのフォールバック経路をそれぞれ確かめる<br/>
        /// 正規化とは、順鎖と逆鎖のどちらから問い合わせても同じ結果になることをいう
        /// </remarks>
        /// <remarks>
        /// 150 bp リードでは k=31 のままだと 31 bp 以上の反復配列がすべて潰れ contig N50 が伸びないため、k=63 前後で正しく動くことは品質上重要
        /// </remarks>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(33)] // UInt128 経路の下限
        [InlineData(63)] // 150 bp リードでの実用値
        [InlineData(64)] // UInt128 経路の上限 (ちょうど 128 bit を使い切る)
        [InlineData(65)] // KmerKey フォールバック経路
        public void V_kmer長が32を超えても含有判定とカバレッジが正しく動く(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            // 逆相補と自己一致しないよう、非周期的な塩基列を決定的に生成する
            var l_配列 = string.Concat(Enumerable.Range(0, p_k長).Select(i => "ACGGTCATTGAC"[(i * 7) % 12]));
            var l_inserted = V_変換_塩基ID列(l_配列);
            var l_revComp = V_変換_塩基ID列(Util.V_逆相補(l_配列));
            var l_neverInserted = V_変換_塩基ID列(new string('T', p_k長));

            for (var i = 0; i < 3; i++)
            {
                l_インデックス.V_登録(l_inserted.AsSpan());
            }
            for (var i = 0; i < 2; i++)
            {
                l_インデックス.V_登録(l_revComp.AsSpan());
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Assert.True(l_インデックス.Get_含まれるか(l_inserted));
            Assert.True(l_インデックス.Get_含まれるか(l_revComp));
            Assert.False(l_インデックス.Get_含まれるか(l_neverInserted));

            // 順鎖 3 回 + 逆鎖 2 回 が同一の正規化キーへ合算されているはず
            Assert.Equal(5UL, l_インデックス.Get_カバレッジ(l_inserted));
            Assert.Equal(5UL, l_インデックス.Get_カバレッジ(l_revComp));
        }

        /// <summary>
        /// k=63 の直鎖配列で、EnumerateTrustedKmers が UInt128 経路でも正しく塩基列へ復元でき (UnpackMid)、隣接判定 (CountOutEdges) が成立することを確認する
        /// </summary>
        /// <remarks>
        /// パック/アンパックの往復が壊れていると unitig 構築が丸ごと機能しなくなるため、経路ごとに固定しておく
        /// </remarks>
        [Fact]
        public void V_UInt128経路でも列挙と次数判定が正しく往復する()
        {
            const int l_k長 = 63;
            // 70 塩基の非周期的な配列 (k=63 の k-mer が 8 個取れる)
            var l_配列 = string.Concat(Enumerable.Range(0, 70).Select(i => "ACGGTCATTGACCTA"[(i * 11) % 15]));

            using var l_インデックス = this.V_構築_直鎖索引(l_配列, l_k長);

            var l_kmers = l_インデックス.Get_信頼kmer一覧().ToList();
            Assert.NotEmpty(l_kmers);
            Assert.All(l_kmers, l_kmer => Assert.Equal(l_k長, l_kmer.Length));
            // 復元した k-mer は必ず集合に含まれていなければならない
            Assert.All(l_kmers, l_kmer => Assert.True(l_インデックス.Get_含まれるか(l_kmer)));

            var l_バイト列 = V_変換_塩基ID列(l_配列);
            Assert.Equal(1, l_インデックス.Get_出次数(l_バイト列.AsSpan(0, l_k長)));
            Assert.Equal(0, l_インデックス.Get_出次数(l_バイト列.AsSpan(l_バイト列.Length - l_k長, l_k長)));
            Assert.Equal(0, l_インデックス.Get_入次数(l_バイト列.AsSpan(0, l_k長)));
        }

        /// <summary>
        /// 直鎖配列の内部では出次数 1 で、末尾では 0 になる
        /// </summary>
        [Fact]
        public void V_直鎖配列の内部では出次数1で末尾では0になる()
        {
            const string l_配列 = "GCTAAAGACAATTACATAACATAC"; // 24 bp、非周期的 (内部に k=8 の重複なし)
            const int l_k長 = 8;
            using var l_インデックス = this.V_構築_直鎖索引(l_配列, l_k長);
            var l_バイト列 = V_変換_塩基ID列(l_配列);

            // 途中の k-mer: ちょうど 1 通りだけ後続がある
            Assert.Equal(1, l_インデックス.Get_出次数(l_バイト列.AsSpan(0, l_k長)));

            // 配列の末尾 k-mer: これ以上後続がない (out-degree 0)
            Assert.Equal(0, l_インデックス.Get_出次数(l_バイト列.AsSpan(l_バイト列.Length - l_k長, l_k長)));
        }

        /// <summary>
        /// 配列の先頭では入次数が 0 になる
        /// </summary>
        [Fact]
        public void V_配列の先頭では入次数が0になる()
        {
            const string l_配列 = "GCTAAAGACAATTACATAACATAC";
            const int l_k長 = 8;
            using var l_インデックス = this.V_構築_直鎖索引(l_配列, l_k長);
            var l_バイト列 = V_変換_塩基ID列(l_配列);

            Assert.Equal(0, l_インデックス.Get_入次数(l_バイト列.AsSpan(0, l_k長)));
        }

        /// <summary>
        /// 除去は、両方の向き (順鎖・逆鎖) を取り除く
        /// </summary>
        [Fact]
        public void V_除去は両方の向きを取り除く()
        {
            const string l_配列 = "GCTAAAGACAATTACATAACATAC";
            const int l_k長 = 8;
            using var l_インデックス = this.V_構築_直鎖索引(l_配列, l_k長);
            var l_kmer = V_変換_塩基ID列(l_配列[..l_k長]);
            var l_revComp = V_変換_塩基ID列(Util.V_逆相補(l_配列[..l_k長]));

            Assert.True(l_インデックス.Get_含まれるか(l_kmer));

            l_インデックス.V_除去(l_kmer);

            Assert.False(l_インデックス.Get_含まれるか(l_kmer));
            Assert.False(l_インデックス.Get_含まれるか(l_revComp));
        }

        /// <summary>
        /// 信頼できる k-mer 一覧は、期待した件数を返す
        /// </summary>
        [Fact]
        public void V_信頼kmer一覧は期待した件数を返す()
        {
            const string l_配列 = "GCTAAAGACAATTACATAACATAC"; // 24 bp、非周期的 (内部に k=8 の重複なし)
            const int l_k長 = 8;
            using var l_インデックス = this.V_構築_直鎖索引(l_配列, l_k長);

            // 24 bp ・ k=8 の非周期的な直鎖配列は 24-8+1=17 個のユニーク k-mer 位置を持ち、
            // 内部に重複 (自己一致・逆相補との一致) がないよう検証済みの配列なので、
            // 正規化後もちょうど 17 件になるはず
            var l_count = l_インデックス.Get_信頼kmer一覧().Count();
            Assert.Equal(17, l_count);
        }

        /// <summary>
        /// 直鎖配列の唯一の開始点を見つける
        /// </summary>
        [Fact]
        public void V_直鎖配列の唯一の開始点を見つける()
        {
            // "ACGT"の繰り返しだと逆相補と自己一致してしまい分岐点が
            // 複雑になるため、非周期的な配列を使う
            const string l_配列 = "GCTAAAGACAATTACATAACATAC"; // 非周期的
            const int l_k長 = 8;
            using var l_インデックス = this.V_構築_直鎖索引(l_配列, l_k長);
            var l_バイト列 = V_変換_塩基ID列(l_配列);

            var l_firstKmers = l_インデックス.Get_開始kmer一覧();

            // 開始 k-mer として、配列の先頭 (またはその正規化された逆鎖) が
            // 含まれているはず
            var l_startKmer = l_バイト列.AsSpan(0, l_k長).ToArray();
            var l_startRevComp = V_変換_塩基ID列(Util.V_逆相補(l_配列[..l_k長]));
            Assert.Contains(l_firstKmers, l_開始kmer => l_開始kmer.SequenceEqual(l_startKmer) || l_開始kmer.SequenceEqual(l_startRevComp));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列を塩基 ID 列へ変換する
        /// </summary>
        /// <param name="p_kmer">元の k-mer</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] V_変換_塩基ID列(string p_kmer)
        {
            return [.. p_kmer.Select(c => c switch { 'A' => Consts.塩基ID.A, 'C' => Consts.塩基ID.C, 'G' => Consts.塩基ID.G, 'T' => Consts.塩基ID.T, _ => throw new InvalidOperationException(), })];
        }

        /// <summary>
        /// 長さ len の直鎖配列 (分岐なし) の全 k-mer をカットオフ以上登録する
        /// </summary>
        /// <remarks>
        /// GraphSimplifier のテストとも共通で使える小さなヘルパー
        /// </remarks>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_kmer長">k-mer 長</param>
        /// <returns>構築した索引</returns>
        private TrustedKmerIndex V_構築_直鎖索引(string p_配列, int p_kmer長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmer長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_バイト列 = V_変換_塩基ID列(p_配列);
            for (var i = 0; i + p_kmer長 <= l_バイト列.Length; i++)
            {
                for (var l_反復回数 = 0; l_反復回数 < 3; l_反復回数++)
                {
                    l_インデックス.V_登録(l_バイト列.AsSpan(i, p_kmer長));
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        #endregion

    }
}
