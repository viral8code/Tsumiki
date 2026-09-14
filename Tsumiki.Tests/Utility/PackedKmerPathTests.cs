using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// パック値で k-mer を扱う高速経路が、塩基列で扱う参照経路と同じ結果を返すことの検証
    /// </summary>
    public class PackedKmerPathTests : IDisposable
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
        public PackedKmerPathTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_packed_path_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
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
        /// 右詰めのパック値とファイル上のバイト列は往復で一致する
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(1)]
        [InlineData(31)]
        [InlineData(64)]
        [InlineData(65)]
        [InlineData(97)]
        [InlineData(128)]
        public void V_パック値とバイト列は往復で一致する(int p_k長)
        {
            var l_乱数 = new Random(p_k長);
            var l_塩基 = Enumerable.Range(0, p_k長).Select(_ => (byte)(l_乱数.Next(4) + 1)).ToArray();
            var l_パック長 = (p_k長 + 3) / 4;
            var l_余り = (8 * l_パック長) - (2 * p_k長);
            var l_値 = TrustedKmerIndex.TryGet_パック_長(l_塩基);
            var l_バイト列 = new byte[l_パック長];

            CountingDB.V_変換_パック済みバイト列(l_値, p_k長, l_バイト列);

            Assert.Equal(l_値, TrustedKmerIndex.Get_読み替え_長(l_バイト列, l_余り));
        }

        /// <summary>
        /// 正規形のパック値は逆相補を渡しても同じで、ローリングで作った値とも一致する
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(21)]
        [InlineData(63)]
        [InlineData(65)]
        [InlineData(89)]
        [InlineData(128)]
        public void V_正規形のパック値はローリングと一致する(int p_k長)
        {
            var l_乱数 = new Random(p_k長 * 7);
            var l_配列 = new string([.. Enumerable.Range(0, p_k長 + 20).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_窓 = new RollingKmer(p_k長);
            for (var i = 0; i < l_配列.Length; i++)
            {
                if (!l_窓.Try追加(l_配列[i], out var l_キー))
                {
                    continue;
                }
                var l_塩基 = Util.V_変換_塩基列(l_配列.Substring(i - p_k長 + 1, p_k長));
                Assert.Equal((l_キー.A_上位, l_キー.A_下位), TrustedKmerIndex.Get_正規形_値(l_塩基));
                Assert.Equal(TrustedKmerIndex.Get_正規形_値(l_塩基), TrustedKmerIndex.Get_正規形_値(Util.V_逆相補(l_塩基.AsSpan()).ToArray()));
            }
        }

        /// <summary>
        /// k=65〜128 でも、登録回数どおりに数え、カットオフ後の所属判定・列挙・除去が逆相補を同一視して働く
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(65)]
        [InlineData(101)]
        [InlineData(128)]
        public void V_長いkでも数えて判定できる(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 3 };
            var l_乱数 = new Random(p_k長);
            var l_配列 = Enumerable.Range(0, p_k長 + 30).Select(_ => (byte)(l_乱数.Next(4) + 1)).ToArray();

            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            for (var i = 0; i + p_k長 <= l_配列.Length; i++)
            {
                for (var l_回数 = 0; l_回数 < 1 + (i % 3); l_回数++)
                {
                    l_索引.V_登録(l_配列.AsSpan(i, p_k長));
                }
            }
            l_索引.V_適用_カットオフ(2UL);

            for (var i = 0; i + p_k長 <= l_配列.Length; i++)
            {
                var l_kmer = l_配列.AsSpan(i, p_k長).ToArray();
                var l_期待回数 = (ulong)(1 + (i % 3));
                Assert.Equal(l_期待回数 >= 2UL, l_索引.Haskmer(l_kmer));
                Assert.Equal(l_期待回数 >= 2UL ? l_期待回数 : 0UL, l_索引.Get_カバレッジ(Util.V_逆相補(l_kmer.AsSpan()).ToArray()));
            }
            Assert.All(l_索引.Get_信頼kmer一覧(), x => Assert.True(l_索引.Haskmer(x)));

            var l_先頭 = l_配列.AsSpan(2, p_k長).ToArray();
            l_索引.V_除去(l_先頭);
            Assert.False(l_索引.Haskmer(l_先頭));
        }

        /// <summary>
        /// パック値で状態を持つ経路探索は、塩基列で状態を持つ参照実装と同じ経路と判定を返す
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(9)]
        [InlineData(33)]
        [InlineData(65)]
        [InlineData(80)]
        public void V_パック値の経路探索は参照実装と一致する(int p_k長)
        {
            var l_乱数 = new Random(p_k長 * 13);
            for (var l_試行 = 0; l_試行 < 40; l_試行++)
            {
                // 短い反復を混ぜた配列から局所集合を作り、一意・複数・到達不能のすべてが現れるようにする
                var l_片 = new string([.. Enumerable.Range(0, p_k長 + 5).Select(_ => "ACGT"[l_乱数.Next(4)])]);
                var l_配列 = string.Concat(Enumerable.Range(0, 6).Select(_ => l_乱数.Next(3) == 0 ? l_片 : new string([.. Enumerable.Range(0, l_乱数.Next(p_k長, 3 * p_k長)).Select(_ => "ACGT"[l_乱数.Next(4)])])));
                var l_塩基 = Util.V_変換_塩基列(l_配列);
                var l_集合 = new LocalKmerSet(p_k長);
                for (var i = 0; i + p_k長 <= l_塩基.Length; i++)
                {
                    l_集合.V_登録(l_塩基.AsSpan(i, p_k長));
                }

                var l_左位置 = l_乱数.Next(0, l_塩基.Length - p_k長);
                var l_右位置 = l_乱数.Next(0, l_塩基.Length - p_k長);
                var l_左 = l_塩基.AsSpan(l_左位置, p_k長).ToArray();
                var l_右 = l_塩基.AsSpan(l_右位置, p_k長).ToArray();
                var l_最大長 = l_乱数.Next(0, 4 * p_k長);

                var l_高速 = ConstrainedPathFinder.Get_経路(l_左, l_右, 0, l_最大長, l_集合, p_k長, 5_000);
                var l_参照 = ConstrainedPathFinder.Get_経路_参照(l_左, l_右, 0, l_最大長, l_集合, p_k長, 5_000);

                Assert.Equal(l_参照, l_高速);
            }
        }

        /// <summary>
        /// k=65〜128 のパック値の walk は、塩基列で walk する参照実装と同じ unitig を返す
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(65)]
        [InlineData(99)]
        public void V_長いkのwalkは参照実装と一致する(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };
            var l_乱数 = new Random(p_k長 * 17);
            var l_片 = new string([.. Enumerable.Range(0, p_k長 + 3).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_配列 = string.Concat(Enumerable.Range(0, 8).Select(_ => l_乱数.Next(3) == 0 ? l_片 : new string([.. Enumerable.Range(0, l_乱数.Next(p_k長, 3 * p_k長)).Select(_ => "ACGT"[l_乱数.Next(4)])])));
            var l_塩基 = Util.V_変換_塩基列(l_配列);

            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            for (var i = 0; i + p_k長 <= l_塩基.Length; i++)
            {
                l_索引.V_登録(l_塩基.AsSpan(i, p_k長));
                l_索引.V_登録(l_塩基.AsSpan(i, p_k長));
            }
            var l_開始点 = l_索引.V_カットオフ(2UL);

            var l_高速 = new UnitigWalk(l_索引, p_k長);
            var l_参照 = new UnitigMaker(l_索引);
            HashSet<(UInt128 A_上位, UInt128 A_下位)> l_訪問済み = [];
            foreach (var l_開始 in l_開始点)
            {
                var l_期待 = l_参照.Get_Unitig(l_開始).A_配列;
                var l_実際 = new string([.. l_高速.Get_塩基列_長(l_開始, l_訪問済み).Select(Util.Get_塩基文字)]);
                Assert.Equal(l_期待, l_実際);
            }
        }

        #endregion
    }
}
