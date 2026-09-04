using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer スペクトルの谷からの -kc 自動選択を、実際に数えた k-mer から
    /// 一気通貫で検証する。
    ///
    /// 既定値の 2 はどのカバレッジ帯にも合わない。実測では同じ検体でも
    /// 35x で 4、100x で 6〜11 が谷であり、2 のままだとエラー由来の k-mer が
    /// 大量に残る(35x の実データで「良い k-mer」が 12.9M と、ゲノムサイズの
    /// 倍に膨れていた)。
    /// </summary>
    public class KmerCutoffSelectorTests : IDisposable
    {
        private readonly string _tempDir;

        public KmerCutoffSelectorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_cutoffsel_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const int K = 21;

        /// <summary>
        /// (出現回数, その回数を持たせる k-mer の種類数)。
        /// 出現回数 8 を底とする谷と、15 を頂点とする単一コピーの山を持つ、
        /// 連続した二峰性スペクトルになるように組んである。
        /// </summary>
        private static readonly (ulong A_出現回数, int A_種類数)[] スペクトルの形 = [
            (1, 2000), (2, 700), (3, 300), (4, 150), (5, 90),
            (6, 70), (7, 60), (8, 58), (9, 70), (10, 120),
            (11, 220), (12, 400), (13, 600), (14, 800), (15, 900),
            (16, 800), (17, 600), (18, 400), (19, 220), (20, 120),
        ];

        private const ulong 谷の位置 = 8;

        /// <summary>
        /// このスペクトルに対して選ばれるべきカットオフ。
        /// 谷(8)ではない。出現回数3以上を残した時点で 5,978 種類となり、
        /// 推定ゲノムサイズ 5,252 の 1.14 倍に収まる(1.2 倍以内)ので、
        /// そこから先へ上げる理由が無い。谷まで上げると本物の k-mer を
        /// 670 種類も余計に削ることになる。
        /// </summary>
        private const ulong 選ばれるべきカットオフ = 3;

        /// <summary>
        /// 上のスペクトルの形どおりに k-mer を登録したインデックスを作る。
        /// 乱数配列から取った連続する k-mer は k=21 なら実質すべて相異なる。
        /// </summary>
        private TrustedKmerIndex BuildIndex()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 4 };
            var index = new TrustedKmerIndex(this._tempDir);

            var 種類数の合計 = スペクトルの形.Sum(x => x.A_種類数);
            var rng = new Random(20260904);
            var bases = string.Concat(Enumerable.Range(0, 種類数の合計 + K - 1).Select(_ => "ACGT"[rng.Next(4)]))
                .Select(Util.Get_塩基ID).ToArray();

            var position = 0;
            foreach (var (出現回数, 種類数) in スペクトルの形)
            {
                for (var i = 0; i < 種類数; i++, position++)
                {
                    for (var t = 0UL; t < 出現回数; t++)
                    {
                        index.V_登録(bases.AsSpan(position, K), p_ワーカー番号: (int)(t % 4));
                    }
                }
            }
            return index;
        }

        [Fact]
        public void Resolve_WhenCutoffWasNotGiven_PicksTheLowestCutoffThatKeepsErrorsFromDominating()
        {
            using var index = this.BuildIndex();
            var param = new Parameters();
            Assert.False(param.A_kmerカットオフが明示指定されたか);

            KmerCutoffSelector.V_解決_kmerカットオフ(param, index);

            Assert.Equal(選ばれるべきカットオフ, param.A_kmerカットオフ);
            // 谷より上へは決して行かないこと。谷で切ると本物の k-mer の左裾まで
            // 削れてグラフが切れる(実データで N50 が半分以下になった)。
            Assert.True(param.A_kmerカットオフ < 谷の位置);
            // 自動適用は「明示指定された」扱いにしない。
            Assert.False(param.A_kmerカットオフが明示指定されたか);
        }

        /// <summary>
        /// 明示指定はユーザーの判断なので、推定値で上書きしてはいけない。
        /// </summary>
        [Fact]
        public void Resolve_WhenCutoffWasGivenExplicitly_LeavesItAlone()
        {
            using var index = this.BuildIndex();
            var param = new Parameters { A_kmerカットオフ = 2 };

            KmerCutoffSelector.V_解決_kmerカットオフ(param, index);

            Assert.Equal(2UL, param.A_kmerカットオフ);
        }

        /// <summary>
        /// 自動選択したカットオフをそのまま適用したとき、実際に残る k-mer が
        /// 選択値と辻褄が合っていること(選択とカットオフの適用が同じ
        /// ヒストグラムを見ている、という一気通貫の確認)。
        /// </summary>
        [Fact]
        public void Resolve_ThenCutoff_KeepsExactlyTheKmersAtOrAboveTheSelectedCutoff()
        {
            using var index = this.BuildIndex();
            var param = new Parameters();

            KmerCutoffSelector.V_解決_kmerカットオフ(param, index);
            _ = index.V_カットオフ(param.A_kmerカットオフ);

            var 残るはずの種類数 = スペクトルの形
                .Where(x => x.A_出現回数 >= param.A_kmerカットオフ)
                .Sum(x => x.A_種類数);

            Assert.Equal(残るはずの種類数, index.Get_信頼kmer一覧().Count());
        }
    }
}
