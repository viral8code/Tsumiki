using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer スペクトルの谷からの -kc 自動選択を、実際に数えた k-mer から
    /// 一気通貫で検証する
    /// </summary>
    /// <remarks>
    /// 既定値の 2 はどのカバレッジ帯にも合わない<br/>
    /// 実測では同じ検体でも
    /// 35 x で 4、100 x で 6〜11 が谷であり、2 のままだとエラー由来の k-mer が
    /// 大量に残る (35 x の実データで「良い k-mer」が 12.9 M と、ゲノムサイズの
    /// 倍に膨れていた)
    /// </remarks>
    public class KmerCutoffSelectorTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public KmerCutoffSelectorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_cutoffsel_tests_" + Guid.NewGuid().ToString("N"));
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
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 21;

        /// <summary>
        /// (出現回数, その回数を持たせる k-mer の種類数)
        /// </summary>
        /// <remarks>
        /// 出現回数 8 を底とする谷と、15 を頂点とする単一コピーの山を持つ、
        /// 連続した二峰性スペクトルになるように組んである
        /// </remarks>
        private static readonly (ulong A_出現回数, int A_種類数)[] スペクトルの形 = [
            (1, 2000), (2, 700), (3, 300), (4, 150), (5, 90),
            (6, 70), (7, 60), (8, 58), (9, 70), (10, 120),
            (11, 220), (12, 400), (13, 600), (14, 800), (15, 900),
            (16, 800), (17, 600), (18, 400), (19, 220), (20, 120),
        ];

        /// <summary>
        /// 谷の位置
        /// </summary>
        private const ulong 谷の位置 = 8UL;

        /// <summary>
        /// このスペクトルに対して選ばれるべきカットオフ
        /// </summary>
        /// <remarks>
        /// 谷 (8) ではない<br/>
        /// V_解決_kmerカットオフ はまず 2 成分混合モデル
        /// (KmerSpectrumMixtureModel) の適合を試み、この形のスペクトルなら
        /// 適合に成功して谷検出 (KmerHistogram) より低い 6 を返す
        /// (事後誤り確率が有意水準を下回る最小の出現回数)<br/>
        /// どちらの経路でも「谷までは上げない」という結論は変わらない
        /// </remarks>
        private const ulong 選ばれるべきカットオフ = 6UL;

        /// <summary>
        /// 上のスペクトルの形どおりに k-mer を登録したインデックスを作る
        /// </summary>
        /// <remarks>
        /// 乱数配列から取った連続する k-mer は k=21 なら実質すべて相異なる
        /// </remarks>
        private TrustedKmerIndex V_構築_索引()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };
            var l_index = new TrustedKmerIndex(this._tempDir);

            var l_種類数の合計 = スペクトルの形.Sum(x => x.A_種類数);
            var l_rng = new Random(20260904);
            var l_bases = string.Concat(Enumerable.Range(0, l_種類数の合計 + k長 - 1).Select(_ => "ACGT"[l_rng.Next(4)]))
                .Select(Util.Get_塩基ID).ToArray();

            var l_position = 0;
            foreach (var (l_出現回数, l_種類数) in スペクトルの形)
            {
                for (var i = 0; i < l_種類数; i++, l_position++)
                {
                    for (var t = 0UL; t < l_出現回数; t++)
                    {
                        l_index.V_登録(l_bases.AsSpan(l_position, k長), p_ワーカー番号: (int)(t % 4));
                    }
                }
            }
            return l_index;
        }

        /// <summary>
        /// カットオフ未指定の場合、誤り由来の k-mer が支配しない最小のカットオフを選ぶ
        /// </summary>
        [Fact]
        public void カットオフ未指定なら誤りが支配しない最小のカットオフを選ぶ()
        {
            using var l_index = this.V_構築_索引();
            var l_param = new Parameters();
            Assert.False(l_param.A_kmerカットオフが明示指定されたか);

            KmerCutoffSelector.V_解決_kmerカットオフ(l_param, l_index);

            Assert.Equal(選ばれるべきカットオフ, l_param.A_kmerカットオフ);
            // 谷より上へは決して行かないこと
            // 谷で切ると本物の k-mer の左裾まで
            // 削れてグラフが切れる (実データで N50 が半分以下になった)
            Assert.True(l_param.A_kmerカットオフ < 谷の位置);
            // 自動適用は「明示指定された」扱いにしない
            Assert.False(l_param.A_kmerカットオフが明示指定されたか);
        }

        /// <summary>
        /// 明示指定はユーザーの判断なので、推定値で上書きしてはいけない
        /// </summary>
        [Fact]
        public void 明示指定されたカットオフはそのまま残す()
        {
            using var l_index = this.V_構築_索引();
            var l_param = new Parameters { A_kmerカットオフ = 2 };

            KmerCutoffSelector.V_解決_kmerカットオフ(l_param, l_index);

            Assert.Equal(2UL, l_param.A_kmerカットオフ);
        }

        /// <summary>
        /// 自動選択したカットオフをそのまま適用したとき、実際に残る k-mer が
        /// 選択値と辻褄が合っていること (選択とカットオフの適用が同じ
        /// ヒストグラムを見ている、という一気通貫の確認)
        /// </summary>
        [Fact]
        public void 選択したカットオフ以上のkmerだけがそのまま残る()
        {
            using var l_index = this.V_構築_索引();
            var l_param = new Parameters();

            KmerCutoffSelector.V_解決_kmerカットオフ(l_param, l_index);
            _ = l_index.V_カットオフ(l_param.A_kmerカットオフ);

            var l_残るはずの種類数 = スペクトルの形
                .Where(x => x.A_出現回数 >= l_param.A_kmerカットオフ)
                .Sum(x => x.A_種類数);

            Assert.Equal(l_残るはずの種類数, l_index.Get_信頼kmer一覧().Count());
        }
    }
}
