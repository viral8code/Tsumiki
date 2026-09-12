using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Model
{
    /// <summary>
    /// k-mer をパックして正規形へ寄せる処理の検証
    /// </summary>
    public class KmerKeyTests
    {
        #region 公開メソッド

        /// <summary>
        /// 大域 k の変更が既存キーの意味を変えない
        /// </summary>
        [Fact]
        public void V_生成後のk変更に影響されない()
        {
            V_設定_k長(31);
            var l_配列 = "ACG" + new string('T', 64);
            var l_キー = new KmerKey(l_配列);
            V_設定_k長(93);
            Assert.Equal(new KmerKey(Util.V_逆相補(l_配列)), l_キー.Get_逆相補());
            Assert.NotEqual(new KmerKey("A"), new KmerKey("AA"));
            Assert.Equal(default(KmerKey), default(KmerKey));
        }

        /// <summary>
        /// バイト列版と char 版のコンストラクタが同じキーを作る
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(4)]   // 1 つの ulong に収まる短い k-mer
        [InlineData(31)]  // デフォルトの k-mer 長
        [InlineData(33)]  // 32 境界をまたぐ長さ (Data が複数 ulong になる)
        [InlineData(64)]  // ちょうど 2 ulong 分
        public void V_バイト列版とchar版で同じキーになる(int p_k長)
        {
            V_設定_k長(p_k長);
            var l_塩基列 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[..p_k長];
            var l_バイトkmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_バイトkmer[i] = l_塩基列[i] switch
                {
                    'A' => Consts.塩基ID.A,
                    'C' => Consts.塩基ID.C,
                    'G' => Consts.塩基ID.G,
                    'T' => Consts.塩基ID.T,
                    _ => throw new InvalidOperationException(),
                };
            }

            var l_char版 = new KmerKey(l_塩基列.AsSpan());
            var l_バイト版 = new KmerKey(l_バイトkmer);

            Assert.True(l_char版.Equals(l_バイト版));
            Assert.Equal(l_char版.GetHashCode(), l_バイト版.GetHashCode());
        }

        /// <summary>
        /// 逆相補が文字列版の逆相補と一致する
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(4)]
        [InlineData(31)]
        [InlineData(33)]
        [InlineData(64)]
        public void V_逆相補が文字列版と一致する(int p_k長)
        {
            V_設定_k長(p_k長);
            var l_フォワード = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[..p_k長];

            var l_期待値 = new KmerKey(Util.V_逆相補(l_フォワード).AsSpan());
            var l_実際 = new KmerKey(l_フォワード.AsSpan()).Get_逆相補();

            Assert.True(l_期待値.Equals(l_実際), $"expected Data=[{string.Join(",", l_期待値.A_パック済みデータ)}] actual Data=[{string.Join(",", l_実際.A_パック済みデータ)}]");
        }

        /// <summary>
        /// 正規形は kmer とその逆相補で同じになる
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(4)]
        [InlineData(31)]
        [InlineData(33)]
        public void V_正規形はkmerと逆相補で一致する(int p_k長)
        {
            V_設定_k長(p_k長);
            var l_フォワード = "ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA"[..p_k長];
            var l_逆相補 = Util.V_逆相補(l_フォワード);

            var l_フォワードキー = new KmerKey(l_フォワード.AsSpan());
            var l_逆相補キー = new KmerKey(l_逆相補.AsSpan());

            Assert.True(l_フォワードキー.Get_正規形().Equals(l_逆相補キー.Get_正規形()));
        }

        /// <summary>
        /// 正規形はべき等
        /// </summary>
        [Fact]
        public void V_正規形はべき等()
        {
            V_設定_k長(31);
            var l_キー = new KmerKey("ACGTGGCCTTAAACGTGGCCTTAAACGTG".PadRight(31, 'A').AsSpan());

            var l_正規形 = l_キー.Get_正規形();

            Assert.True(l_正規形.Equals(l_正規形.Get_正規形()));
        }

        /// <summary>
        /// 異なる kmer の正規形は区別される
        /// </summary>
        [Fact]
        public void V_異なるkmerの正規形は区別される()
        {
            V_設定_k長(4);
            var l_甲 = new KmerKey("ACGT".AsSpan()).Get_正規形();
            var l_乙 = new KmerKey("TTTT".AsSpan()).Get_正規形();

            Assert.False(l_甲.Equals(l_乙));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 実行時引数の k 長を差し替える
        /// </summary>
        /// <param name="p_k長">設定する k 長</param>
        private static void V_設定_k長(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長 };
        }

        #endregion

    }
}
