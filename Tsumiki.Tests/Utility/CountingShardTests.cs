using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer のカウントが、シャード数 (スレッド数) に依らず正確であることを確認する
    /// </summary>
    /// <remarks>
    /// k-mer をワーカー単位ではなくハッシュ値でシャードへ振り分けるようにした際、
    /// 実データでカウントがちょうど 2 倍になる不具合が出た (ヒストグラムが
    /// 偶数のカウントしか持たない、という形で表面化した)<br/>
    /// スレッド数を変えて同じ答えになることを固定しておく
    /// </remarks>
    public class CountingShardTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public CountingShardTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_shard_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
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
        /// シャード数(スレッド数)を変えても、登録した回数どおりに数えられること
        /// </summary>
        [Theory]
        [InlineData(1)]
        [InlineData(2)]
        [InlineData(8)]
        public void 登録回数どおりに数えられる_シャード数によらず(int p_スレッド数)
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = p_スレッド数 };

            using var index = new TrustedKmerIndex(this._tempDir);

            var l_配列 = "ACGGTCATTGACCTAGGATCA"; // 21 塩基
            var l_kmer = l_配列.Select(Util.Get_塩基ID).ToArray();

            // ちょうど 7 回登録する (奇数にして「2 倍になっていないか」を確実に見る)
            for (var i = 0; i < 7; i++)
            {
                index.V_登録(l_kmer.AsSpan(), p_ワーカー番号: i % p_スレッド数);
            }

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.Equal(7UL, index.Get_カバレッジ(l_kmer));
        }

        /// <summary>
        /// 多数の異なる k-mer を、それぞれ異なる回数だけ登録しても
        /// 正確に数えられること (シャード分割とマージの整合性)
        /// </summary>
        [Theory]
        [InlineData(1)]
        [InlineData(8)]
        public void 多数の異なるkmerでも正確に数えられる(int p_スレッド数)
        {
            const int k = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = p_スレッド数 };

            using var index = new TrustedKmerIndex(this._tempDir);

            var l_乱数 = new Random(1234);
            var l_配列 = string.Concat(Enumerable.Range(0, 500).Select(_ => "ACGT"[l_乱数.Next(4)]));
            var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();

            // 位置 i の k-mer を (i % 5) + 2 回登録する
            var l_期待値 = new Dictionary<int, ulong>();
            for (var i = 0; i + k <= l_塩基列.Length; i++)
            {
                var l_回数 = (ulong)((i % 5) + 2);
                l_期待値[i] = l_回数;
                for (ulong t = 0UL; t < l_回数; t++)
                {
                    index.V_登録(l_塩基列.AsSpan(i, k), p_ワーカー番号: (int)(t % (ulong)p_スレッド数));
                }
            }

            _ = index.V_カットオフ(p_カットオフ: 2);

            foreach (var (l_位置, l_回数) in l_期待値)
            {
                Assert.Equal(l_回数, index.Get_カバレッジ(l_塩基列.AsSpan(l_位置, k)));
            }
        }
    }
}
