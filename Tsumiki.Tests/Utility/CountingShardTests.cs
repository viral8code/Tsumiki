using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer のカウントが、シャード数 (スレッド数) に依らず正確であることを確認する
    /// </summary>
    public class CountingShardTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki shard tests
        /// </summary>
        private const string C_項目_tsumiki_shard_tests = "tsumiki_shard_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 塩基配列 ACGGTCATTGACCTAGGATCA
        /// </summary>
        private const string C_塩基配列_ACGGTCATTGACCTAGGATCA = "ACGGTCATTGACCTAGGATCA";

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
        public CountingShardTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_shard_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// シャード数 (スレッド数) を変えても、登録した回数どおりに数えられること
        /// </summary>
        /// <param name="p_スレッド数"></param>
        [Theory]
        [InlineData(1)]
        [InlineData(2)]
        [InlineData(8)]
        public void V_登録回数どおりに数えられる_シャード数によらず(int p_スレッド数)
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = p_スレッド数
            };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_配列 = C_塩基配列_ACGGTCATTGACCTAGGATCA;
            var l_kmer = l_配列.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i < 7; i++)
            {
                l_インデックス.V_登録(l_kmer.AsSpan());
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            Assert.Equal(7UL, l_インデックス.Get_カバレッジ(l_kmer));
        }

        /// <summary>
        /// 多数の異なる k-mer を、それぞれ異なる回数だけ登録しても正確に数えられること (シャード分割とマージの整合性)
        /// </summary>
        /// <param name="p_スレッド数"></param>
        [Theory]
        [InlineData(1)]
        [InlineData(8)]
        public void V_多数の異なるkmerでも正確に数えられる(int p_スレッド数)
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = p_スレッド数
            };
            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_乱数 = new Random(1_234);
            var l_配列 = string.Concat(Enumerable.Range(0, 500).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
            var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
            var l_期待値 = new Dictionary<int, ulong>();
            for (var i = 0; i + l_k長 <= l_塩基列.Length; i++)
            {
                var l_回数 = (ulong)((i % 5) + 2);
                l_期待値[i] = l_回数;
                for (ulong t = 0UL; t < l_回数; t++)
                {
                    l_インデックス.V_登録(l_塩基列.AsSpan(i, l_k長));
                }
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            foreach (var (l_位置, l_回数)in l_期待値)
            {
                Assert.Equal(l_回数, l_インデックス.Get_カバレッジ(l_塩基列.AsSpan(l_位置, l_k長)));
            }
        }

        #endregion
    }
}
