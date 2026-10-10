using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer スペクトルの谷からの -kc 自動選択を、実際に数えた k-mer から一気通貫で検証する
    /// </summary>
    public class KmerCutoffSelectorTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki cutoffsel tests
        /// </summary>
        private const string C_項目_tsumiki_cutoffsel_tests = "tsumiki_cutoffsel_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 21;

        /// <summary>
        /// 谷の位置
        /// </summary>
        private const ulong C_谷の位置 = 8UL;

        /// <summary>
        /// このスペクトルに対して選ばれるべきカットオフ
        /// </summary>
        private const ulong C_選ばれるべきカットオフ = 6UL;

        /// <summary>
        /// (出現回数, その回数を持たせる k-mer の種類数)
        /// </summary>
        private static readonly (ulong A_出現回数, int A_種類数)[] C_スペクトルの形 = [
            (1UL, 2_000),
            (2UL, 700),
            (3UL, 300),
            (4UL, 150),
            (5UL, 90),
            (6UL, 70),
            (7UL, 60),
            (8UL, 58),
            (9UL, 70),
            (10UL, 120),
            (11UL, 220),
            (12UL, 400),
            (13UL, 600),
            (14UL, 800),
            (15UL, 900),
            (16UL, 800),
            (17UL, 600),
            (18UL, 400),
            (19UL, 220),
            (20UL, 120),
        ];

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
        public KmerCutoffSelectorTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_cutoffsel_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// カットオフ未指定の場合、誤り由来の k-mer が支配しない最小のカットオフを選ぶ
        /// </summary>
        [Fact]
        public void V_カットオフ未指定なら誤りが支配しない最小のカットオフを選ぶ()
        {
            using var l_インデックス = this.V_構築_索引();
            var l_param = new Parameters();
            Assert.False(l_param.A_Iskmerカットオフ明示指定);
            KmerCutoffSelector.V_解決_kmerカットオフ(l_param, l_インデックス);
            Assert.Equal(C_選ばれるべきカットオフ, l_param.A_kmerカットオフ);
            Assert.True(l_param.A_kmerカットオフ < C_谷の位置);
            Assert.False(l_param.A_Iskmerカットオフ明示指定);
        }

        /// <summary>
        /// 明示指定はユーザーの判断なので、推定値で上書きしてはいけない
        /// </summary>
        [Fact]
        public void V_明示指定されたカットオフはそのまま残す()
        {
            using var l_インデックス = this.V_構築_索引();
            var l_param = new Parameters
            {
                A_kmerカットオフ = 2UL
            };
            KmerCutoffSelector.V_解決_kmerカットオフ(l_param, l_インデックス);
            Assert.Equal(2UL, l_param.A_kmerカットオフ);
        }

        /// <summary>
        /// 自動選択したカットオフをそのまま適用したとき、実際に残る k-mer が選択値と辻褄が合っていること (選択とカットオフの適用が同じヒストグラムを見ている、という一気通貫の確認)
        /// </summary>
        [Fact]
        public void V_選択したカットオフ以上のkmerだけがそのまま残る()
        {
            using var l_インデックス = this.V_構築_索引();
            var l_param = new Parameters();
            KmerCutoffSelector.V_解決_kmerカットオフ(l_param, l_インデックス);
            _ = l_インデックス.V_カットオフ(l_param.A_kmerカットオフ);
            var l_残るはずの種類数 = C_スペクトルの形.Where(x => x.A_出現回数 >= l_param.A_kmerカットオフ).Sum(x => x.A_種類数);
            Assert.Equal(l_残るはずの種類数, l_インデックス.Get_信頼kmer一覧().Count());
        }

        /// <summary>
        /// 谷がモデルのカットオフの 2 倍を超えるときだけ谷を採り、それ以外はモデルのカットオフを残す
        /// </summary>
        /// <param name="p_モデルのカットオフ"></param>
        /// <param name="p_谷"></param>
        /// <param name="p_期待"></param>
        [Theory]
        [InlineData(10UL, 42UL, 42UL)]
        [InlineData(13UL, 14UL, 13UL)]
        [InlineData(7UL, 9UL, 7UL)]
        [InlineData(10UL, 20UL, 10UL)]
        [InlineData(10UL, 21UL, 21UL)]
        [InlineData(19UL, 14UL, 19UL)]
        public void V_谷が明確に食い違うときだけ谷を採る(ulong p_モデルのカットオフ, ulong p_谷, ulong p_期待)
        {
            Assert.Equal(p_期待, KmerCutoffSelector.Get_谷で補正したカットオフ(p_モデルのカットオフ, p_谷));
        }

        /// <summary>
        /// 谷が求められなければモデルのカットオフを残す
        /// </summary>
        [Fact]
        public void V_谷が無ければモデルのカットオフを残す()
        {
            Assert.Equal(10UL, KmerCutoffSelector.Get_谷で補正したカットオフ(10UL, null));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 上のスペクトルの形どおりに k-mer を登録したインデックスを作る
        /// </summary>
        /// <returns></returns>
        private TrustedKmerIndex V_構築_索引()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 4
            };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_種類数の合計 = C_スペクトルの形.Sum(x => x.A_種類数);
            var l_乱数生成器 = new Random(20_260_904);
            var l_bases = string.Concat(Enumerable.Range(0, l_種類数の合計 + C_k長 - 1).Select(_ => Consts.塩基文字[l_乱数生成器.Next(4)])).Select(Util.Get_塩基ID).ToArray();
            var l_位置 = 0;
            foreach (var (l_出現回数, l_種類数)in C_スペクトルの形)
            {
                for (var i = 0; i < l_種類数; i++, l_位置++)
                {
                    for (var t = 0UL; t < l_出現回数; t++)
                    {
                        l_インデックス.V_登録(l_bases.AsSpan(l_位置, C_k長));
                    }
                }
            }

            return l_インデックス;
        }

        #endregion
    }
}
