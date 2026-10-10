using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Model
{
    /// <summary>
    /// k-mer をパックして正規形へ寄せる処理の検証
    /// </summary>
    public class KmerKeyTests
    {
        #region 定数

        /// <summary>
        /// 塩基配列 ACG
        /// </summary>
        private const string C_塩基配列_ACG = "ACG";

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        /// <summary>
        /// 塩基配列 AA
        /// </summary>
        private const string C_塩基配列_AA = "AA";

        /// <summary>
        /// 塩基配列 ACGTACGTACGTACGTACGTACGT ACGTACGTACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGTACGTACGTACGTACGT_ACGTACGTACGTACGT = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        /// <summary>
        /// 区切り
        /// </summary>
        private const string C_区切り = ",";

        /// <summary>
        /// 塩基配列 ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA
        /// </summary>
        private const string C_塩基配列_ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA = "ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA";

        /// <summary>
        /// 塩基配列 ACGTGGCCTTAAACGTGGCCTTAAACGTG
        /// </summary>
        private const string C_塩基配列_ACGTGGCCTTAAACGTGGCCTTAAACGTG = "ACGTGGCCTTAAACGTGGCCTTAAACGTG";

        /// <summary>
        /// 塩基配列 TTTT
        /// </summary>
        private const string C_塩基配列_TTTT = "TTTT";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 大域 k の変更が既存キーの意味を変えない
        /// </summary>
        [Fact]
        public void V_生成後のk変更に影響されない()
        {
            V_設定_k長(31);
            var l_配列 = C_塩基配列_ACG + new string('T', 64);
            var l_キー = new KmerKey(l_配列);
            V_設定_k長(93);
            Assert.Equal(new KmerKey(Util.V_逆相補(l_配列)), l_キー.Get_逆相補());
            Assert.NotEqual(new KmerKey(C_項目_A), new KmerKey(C_塩基配列_AA));
            Assert.Equal(default(KmerKey), default(KmerKey));
        }

        /// <summary>
        /// バイト列版と char 版のコンストラクタが同じキーを作る
        /// </summary>
        /// <param name="p_k長"></param>
        [Theory]
        [InlineData(4)]
        [InlineData(31)]
        [InlineData(33)]
        [InlineData(64)]
        public void V_バイト列版とchar版で同じキーになる(int p_k長)
        {
            V_設定_k長(p_k長);
            var l_塩基列 = C_塩基配列_ACGTACGTACGTACGTACGTACGT_ACGTACGTACGTACGT[..p_k長];
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
            var l_フォワード = C_塩基配列_ACGTACGTACGTACGTACGTACGT_ACGTACGTACGTACGT[..p_k長];
            var l_期待値 = new KmerKey(Util.V_逆相補(l_フォワード).AsSpan());
            var l_実際 = new KmerKey(l_フォワード.AsSpan()).Get_逆相補();
            Assert.True(l_期待値.Equals(l_実際), $"expected Data=[{string.Join(C_区切り, l_期待値.A_パック済みデータ)}] actual Data=[{string.Join(C_区切り, l_実際.A_パック済みデータ)}]");
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
            var l_フォワード = C_塩基配列_ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA[..p_k長];
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
            var l_キー = new KmerKey(C_塩基配列_ACGTGGCCTTAAACGTGGCCTTAAACGTG.PadRight(31, 'A').AsSpan());
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
            var l_甲 = new KmerKey(Consts.塩基文字.AsSpan()).Get_正規形();
            var l_乙 = new KmerKey(C_塩基配列_TTTT.AsSpan()).Get_正規形();
            Assert.False(l_甲.Equals(l_乙));
        }

        /// <summary>
        /// 塩基 ID 列から直接作った正規形は、逆相補を組み立てて比べた正規形と同じ
        /// </summary>
        /// <param name="p_長さ">k-mer の長さ</param>
        [Theory]
        [InlineData(33)]
        [InlineData(64)]
        [InlineData(65)]
        [InlineData(135)]
        [InlineData(300)]
        public void Get_正規形_逆相補を組み立てた正規形と同じ(int p_長さ)
        {
            var l_乱数 = new Random(p_長さ);
            for (var l_回 = 0; l_回 < 200; l_回++)
            {
                var l_kmer = new byte[p_長さ];
                for (var i = 0; i < p_長さ; i++)
                {
                    l_kmer[i] = (byte)l_乱数.Next(1, 5);
                }

                if (l_回 % 4 == 0)
                {
                    for (var i = 0; i < p_長さ / 2; i++)
                    {
                        l_kmer[p_長さ - 1 - i] = (byte)(5 - l_kmer[i]);
                    }
                }

                Assert.Equal(new KmerKey(l_kmer).Get_正規形(), KmerKey.Get_正規形(l_kmer));
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 実行時引数の k 長を差し替える
        /// </summary>
        /// <param name="p_k長">設定する k 長</param>
        private static void V_設定_k長(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = p_k長
            };
        }

        #endregion
    }
}
