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
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 21;

        /// <summary>
        /// 谷の位置
        /// </summary>
        private const ulong 谷の位置 = 8UL;

        /// <summary>
        /// このスペクトルに対して選ばれるべきカットオフ
        /// </summary>
        /// <remarks>
        /// 谷 (8) ではない<br/>
        /// V_解決_kmer カットオフ はまず 2 成分混合モデル (KmerSpectrumMixtureModel) の適合を試み、この形のスペクトルなら適合に成功して谷検出 (KmerHistogram) より低い 6 を返す (事後誤り確率が有意水準を下回る最小の出現回数) <br/>
        /// どちらの経路でも「谷までは上げない」という結論は変わらない
        /// </remarks>
        private const ulong 選ばれるべきカットオフ = 6UL;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        /// <summary>
        /// (出現回数, その回数を持たせる k-mer の種類数)
        /// </summary>
        /// <remarks>
        /// 出現回数 8 を底とする谷と、15 を頂点とする単一コピーの山を持つ、連続した二峰性スペクトルになるように組んである
        /// </remarks>
        private static readonly (ulong A_出現回数, int A_種類数)[] スペクトルの形 = [
            (1UL, 2_000), (2UL, 700), (3UL, 300), (4UL, 150), (5UL, 90),
            (6UL, 70), (7UL, 60), (8UL, 58), (9UL, 70), (10UL, 120),
            (11UL, 220), (12UL, 400), (13UL, 600), (14UL, 800), (15UL, 900),
            (16UL, 800), (17UL, 600), (18UL, 400), (19UL, 220), (20UL, 120),
        ];

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public KmerCutoffSelectorTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_cutoffsel_tests_" + Guid.NewGuid().ToString("N"));
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

            Assert.Equal(選ばれるべきカットオフ, l_param.A_kmerカットオフ);
            // 谷より上へは決して行かないこと
            // 谷で切ると本物の k-mer の左裾まで
            // 削れてグラフが切れる (実データで N50 が半分以下になった)
            Assert.True(l_param.A_kmerカットオフ < 谷の位置);
            // 自動適用は「明示指定された」扱いにしない
            Assert.False(l_param.A_Iskmerカットオフ明示指定);
        }

        /// <summary>
        /// 明示指定はユーザーの判断なので、推定値で上書きしてはいけない
        /// </summary>
        [Fact]
        public void V_明示指定されたカットオフはそのまま残す()
        {
            using var l_インデックス = this.V_構築_索引();
            var l_param = new Parameters { A_kmerカットオフ = 2UL };

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

            var l_残るはずの種類数 = スペクトルの形
                .Where(x => x.A_出現回数 >= l_param.A_kmerカットオフ)
                .Sum(x => x.A_種類数);

            Assert.Equal(l_残るはずの種類数, l_インデックス.Get_信頼kmer一覧().Count());
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 上のスペクトルの形どおりに k-mer を登録したインデックスを作る
        /// </summary>
        /// <remarks>
        /// 乱数配列から取った連続する k-mer は k=21 なら実質すべて相異なる
        /// </remarks>
        /// <returns></returns>
        private TrustedKmerIndex V_構築_索引()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            var l_種類数の合計 = スペクトルの形.Sum(x => x.A_種類数);
            var l_乱数生成器 = new Random(20_260_904);
            var l_bases = string.Concat(Enumerable.Range(0, l_種類数の合計 + k長 - 1).Select(_ => "ACGT"[l_乱数生成器.Next(4)]))
                .Select(Util.Get_塩基ID).ToArray();

            var l_位置 = 0;
            foreach (var (l_出現回数, l_種類数) in スペクトルの形)
            {
                for (var i = 0; i < l_種類数; i++, l_位置++)
                {
                    for (var t = 0UL; t < l_出現回数; t++)
                    {
                        l_インデックス.V_登録(l_bases.AsSpan(l_位置, k長));
                    }
                }
            }
            return l_インデックス;
        }

        #endregion

    }
}
