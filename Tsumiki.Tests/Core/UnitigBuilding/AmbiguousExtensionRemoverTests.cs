using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 控えに同水準の兄弟がいる一本道を断つ処理の検証
    /// </summary>
    public class AmbiguousExtensionRemoverTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// k 長
        /// </summary>
        private const int k長 = 8;

        /// <summary>
        /// k=8 で内部に重複する k-mer を持たない配列
        /// </summary>
        private const string 主配列 = "GCTAAAGACAATTACATAACATACGGATCCTTAGGCAATTGACCTGAAT";

        /// <summary>
        /// 分岐が生じる位置 (この位置から始まる k-mer が競合する)
        /// </summary>
        private const int 分岐位置 = 20;

        /// <summary>
        /// カットオフ
        /// </summary>
        private const ulong カットオフ = 10UL;

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
        public AmbiguousExtensionRemoverTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_ambiguous_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
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
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 残った続きと同水準の兄弟が控えにいる一本道は断たれる
        /// </summary>
        [Fact]
        public void V_控えに同水準の兄弟がいる一本道は断つ()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_索引, 主配列, 12);
            // 分岐位置の k-mer と 1 塩基だけ違う兄弟を、カットオフのすぐ下かつほぼ同じ深さに置く
            V_登録(l_索引, Get_兄弟(分岐位置), 9);
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);

            var l_断つ前 = l_索引.Haskmer(Get_kmer(主配列, 分岐位置));
            var l_除去数 = AmbiguousExtensionRemover.Get_除去数(l_索引, k長);

            Assert.True(l_断つ前);
            Assert.Equal(1, l_除去数);
            Assert.False(l_索引.Haskmer(Get_kmer(主配列, 分岐位置)));
        }

        /// <summary>
        /// 残った続きが控えの兄弟より十分に優勢なら断たない
        /// </summary>
        [Fact]
        public void V_十分に優勢な一本道は断たない()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_索引, 主配列, 30);
            V_登録(l_索引, Get_兄弟(分岐位置), 9);
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);

            var l_除去数 = AmbiguousExtensionRemover.Get_除去数(l_索引, k長);

            Assert.Equal(0, l_除去数);
            Assert.True(l_索引.Haskmer(Get_kmer(主配列, 分岐位置)));
        }

        /// <summary>
        /// 控えが無ければ何も断たない
        /// </summary>
        [Fact]
        public void V_控えが無ければ断たない()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_索引, 主配列, 12);
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);

            var l_除去数 = AmbiguousExtensionRemover.Get_除去数(l_索引, k長);

            Assert.Equal(0, l_除去数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 指定位置の k-mer の末尾 1 塩基だけを別の塩基に替えたものを返す
        /// </summary>
        /// <param name="p_位置"></param>
        /// <returns></returns>
        private static byte[] Get_兄弟(int p_位置)
        {
            var l_kmer = Get_kmer(主配列, p_位置);
            l_kmer[^1] = l_kmer[^1] == Consts.塩基ID.A ? Consts.塩基ID.C : Consts.塩基ID.A;
            return l_kmer;
        }

        /// <summary>
        /// 配列の全 k-mer を同じ深さで登録する
        /// </summary>
        /// <param name="p_索引"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_深さ"></param>
        private static void V_登録_全kmer(TrustedKmerIndex p_索引, string p_配列, int p_深さ)
        {
            for (var i = 0; i + k長 <= p_配列.Length; i++)
            {
                V_登録(p_索引, Get_kmer(p_配列, i), p_深さ);
            }
        }

        /// <summary>
        /// k-mer を指定回数だけ登録する
        /// </summary>
        /// <param name="p_索引"></param>
        /// <param name="p_kmer"></param>
        /// <param name="p_回数"></param>
        private static void V_登録(TrustedKmerIndex p_索引, byte[] p_kmer, int p_回数)
        {
            for (var i = 0; i < p_回数; i++)
            {
                p_索引.V_登録(p_kmer);
            }
        }

        /// <summary>
        /// 配列の位置から k-mer を塩基 ID 列で取り出す
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_位置"></param>
        /// <returns></returns>
        private static byte[] Get_kmer(string p_配列, int p_位置)
        {
            return [.. p_配列.Substring(p_位置, k長).Select(Util.Get_塩基ID)];
        }

        #endregion
    }
}
