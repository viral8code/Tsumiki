using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// r-mer 検証器の、信頼できる k-mer による登録の絞り込みの検証
    /// </summary>
    public class RepeatRMerVerifierFilterTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki rmer filter tests
        /// </summary>
        private const string C_項目_tsumiki_rmer_filter_tests = "tsumiki_rmer_filter_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// 塩基配列 AAACCCGG
        /// </summary>
        private const string C_塩基配列_AAACCCGG = "AAACCCGG";

        /// <summary>
        /// 塩基配列 GATTACAG
        /// </summary>
        private const string C_塩基配列_GATTACAG = "GATTACAG";

        /// <summary>
        /// 塩基配列 CCTTGGAA
        /// </summary>
        private const string C_塩基配列_CCTTGGAA = "CCTTGGAA";

        /// <summary>
        /// r-mer 長
        /// </summary>
        private const int C_r長 = 18;

        /// <summary>
        /// アセンブリ側の k
        /// </summary>
        private const int C_アセンブリk長 = 8;

        /// <summary>
        /// 反復の手前にある配列
        /// </summary>
        private const string C_反復前配列 = "ACAGTTCGCGAGCCCTCCGTC";

        /// <summary>
        /// 手前と先に挟まれた反復配列
        /// </summary>
        private const string C_反復配列 = "CTCCGTCAGCTTGTTTGGAGCAGA";

        /// <summary>
        /// 反復の先にある配列
        /// </summary>
        private const string C_反復後配列 = "GAGCAGAGTCGTTCTGCGAGG";

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
        public RepeatRMerVerifierFilterTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_rmer_filter_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_アセンブリk長,
                A_スレッド数 = 2
            };
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
        /// グラフ上の経路の接合点支持は絞り込みの有無で変わらず、端の k-mer が集合に無い r-mer だけが登録から外れる
        /// </summary>
        [Fact]
        public void V_絞り込んでも経路の支持数は変わらずエラーを含むrMerだけが外れる()
        {
            var l_正しい配列 = C_反復前配列 + C_反復配列[(C_アセンブリk長 - 1)..] + C_反復後配列[(C_アセンブリk長 - 1)..];
            var l_誤りの反復 = C_反復配列[..12] + (C_反復配列[12] == 'A' ? 'C' : 'A') + C_反復配列[13..];
            var l_誤りの配列 = C_反復前配列 + l_誤りの反復[(C_アセンブリk長 - 1)..] + C_反復後配列[(C_アセンブリk長 - 1)..];
            var l_リードパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            File.WriteAllText(l_リードパス, $"@r1\n{l_正しい配列}\n+\n{new string('I', l_正しい配列.Length)}\n@r2\n{l_誤りの配列}\n+\n{new string('I', l_誤りの配列.Length)}\n");
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            var l_塩基列 = l_正しい配列.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + C_アセンブリk長 <= l_塩基列.Length; i++)
            {
                l_索引.V_登録(l_塩基列.AsSpan(i, C_アセンブリk長));
                l_索引.V_登録(l_塩基列.AsSpan(i, C_アセンブリk長));
            }

            l_索引.V_適用_カットオフ(2UL);
            var l_絞り込みなし = RepeatRMerVerifier.V_構築([l_リードパス], C_r長);
            var l_絞り込みあり = RepeatRMerVerifier.V_構築([l_リードパス], C_r長, l_索引, C_アセンブリk長);
            var l_正しい経路の支持 = l_絞り込みなし.Get_接合点の支持数(C_反復前配列, C_反復配列, C_反復後配列);
            Assert.True(l_正しい経路の支持 > 0);
            Assert.Equal(l_正しい経路の支持, l_絞り込みあり.Get_接合点の支持数(C_反復前配列, C_反復配列, C_反復後配列));
            Assert.True(l_絞り込みあり.Get_接合点の支持数(C_反復前配列, l_誤りの反復, C_反復後配列) < l_絞り込みなし.Get_接合点の支持数(C_反復前配列, l_誤りの反復, C_反復後配列));
        }

        /// <summary>
        /// カットオフ未満でも控えの下限以上の k-mer は信頼集合に入れずに控え、下限未満は捨てる
        /// </summary>
        [Fact]
        public void V_控えはカットオフ未満かつ下限以上のkmerだけを保持する()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録(l_索引, C_塩基配列_AAACCCGG, 5);
            V_登録(l_索引, C_塩基配列_GATTACAG, 3);
            V_登録(l_索引, C_塩基配列_CCTTGGAA, 1);
            l_索引.V_適用_カットオフ(5UL, 2UL);
            Assert.True(l_索引.Haskmer(Get_塩基列(C_塩基配列_AAACCCGG)));
            Assert.False(l_索引.Haskmer(Get_塩基列(C_塩基配列_GATTACAG)));
            Assert.Equal(3UL, l_索引.Get_控えカバレッジ(Get_塩基列(Util.V_逆相補(C_塩基配列_GATTACAG))));
            Assert.Equal(0UL, l_索引.Get_控えカバレッジ(Get_塩基列(C_塩基配列_CCTTGGAA)));
            Assert.Equal(1, l_索引.A_控えkmer数);
            l_索引.V_解放_控え();
            Assert.Equal(0, l_索引.A_控えkmer数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer を指定回数だけ登録する
        /// </summary>
        /// <param name="p_索引"></param>
        /// <param name="p_kmer"></param>
        /// <param name="p_回数"></param>
        private static void V_登録(TrustedKmerIndex p_索引, string p_kmer, int p_回数)
        {
            for (var i = 0; i < p_回数; i++)
            {
                p_索引.V_登録(Get_塩基列(p_kmer));
            }
        }

        /// <summary>
        /// 配列を塩基 ID 列にする
        /// </summary>
        /// <param name="p_配列"></param>
        /// <returns></returns>
        private static byte[] Get_塩基列(string p_配列)
        {
            return[..p_配列.Select(Util.Get_塩基ID)];
        }

        #endregion
    }
}
