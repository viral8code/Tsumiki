using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Models.Foundation;
using Tsumiki.Tests.Utility;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core.Preprocessing
{
    /// <summary>
    /// k &gt; 128 の数え上げを束でまとめても、窓ごとに 1 件ずつ登録したときと同じ数になることを検証する
    /// </summary>
    [Collection(中間データ置き場の集まり.C_名前)]
    public class KmerCountingWideTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki wide count tests
        /// </summary>
        private const string C_項目_tsumiki_wide_count_tests = "tsumiki_wide_count_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// FASTQ 品質区切り
        /// </summary>
        private const string C_FASTQ品質区切り = "+";

        /// <summary>
        /// 項目 batch
        /// </summary>
        private const string C_項目_batch = "batch";

        /// <summary>
        /// 項目 window
        /// </summary>
        private const string C_項目_window = "window";

        /// <summary>
        /// 確かめる k 長
        /// </summary>
        private const int C_k長 = 139;

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_wide_count_tests + Guid.NewGuid().ToString(C_GUID書式));

        /// <summary>
        /// 元の実行時引数
        /// </summary>
        private readonly Parameters _元の設定 = ConfigurationManager.A_実行時引数;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 作業ディレクトリと k 長を用意する
        /// </summary>
        public KmerCountingWideTests()
        {
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 4
            };
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 反復・回文・N・低品質を含むリードで、k-mer ごとの数が窓ごとの登録と一致する
        /// </summary>
        [Fact]
        public void V_束でまとめても窓ごとの登録と同じ数になる()
        {
            var l_乱数 = new Random(7);
            var l_ゲノム = new string([..Enumerable.Range(0, 3_000).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_回文 = l_ゲノム[..150] + Util.V_逆相補_曖昧塩基あり(l_ゲノム[..150]);
            List<(string A_配列, string A_クオリティ)> l_リード群 = [];
            for (var i = 0; i < 400; i++)
            {
                var l_開始 = l_乱数.Next(l_ゲノム.Length - 250);
                var l_配列 = (i % 7 == 0 ? l_回文 : l_ゲノム.Substring(l_開始, 250)).ToCharArray();
                var l_クオリティ = Enumerable.Repeat('I', l_配列.Length).ToArray();
                if (i % 5 == 0)
                {
                    l_配列[l_乱数.Next(l_配列.Length)] = 'N';
                }

                if (i % 3 == 0)
                {
                    l_クオリティ[l_乱数.Next(l_配列.Length)] = '!';
                }

                l_リード群.Add((new string(l_配列), new string(l_クオリティ)));
            }

            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            File.WriteAllLines(l_パス, l_リード群.SelectMany((x, i) => new[] { $"@r{i}", x.A_配列, C_FASTQ品質区切り, x.A_クオリティ }));
            var l_束の作業 = Directory.CreateDirectory(Path.Combine(this._作業ディレクトリ, C_項目_batch)).FullName;
            var l_窓の作業 = Directory.CreateDirectory(Path.Combine(this._作業ディレクトリ, C_項目_window)).FullName;
            using var l_束 = new TrustedKmerIndex(l_束の作業);
            using var l_窓 = new TrustedKmerIndex(l_窓の作業);
            KmerCounting.V_読込_リードファイル(l_パス, l_束, 33);
            foreach (var (l_配列, l_クオリティ)in l_リード群)
            {
                var l_塩基列 = Util.V_変換_塩基列(l_配列);
                for (var i = 0; i + C_k長 <= l_塩基列.Length; i++)
                {
                    if (Enumerable.Range(i, C_k長).All(j => l_塩基列[j] != Consts.無効な塩基 && l_クオリティ[j] - 33 >= ConfigurationManager.A_実行時引数.A_クオリティカットオフ))
                    {
                        l_窓.V_登録(l_塩基列.AsSpan(i, C_k長));
                    }
                }
            }

            l_束.V_適用_カットオフ(1UL);
            l_窓.V_適用_カットオフ(1UL);
            foreach (var (l_配列, _)in l_リード群)
            {
                var l_塩基列 = Util.V_変換_塩基列(l_配列);
                for (var i = 0; i + C_k長 <= l_塩基列.Length; i++)
                {
                    if (l_塩基列.AsSpan(i, C_k長).Contains(Consts.無効な塩基))
                    {
                        continue;
                    }

                    Assert.Equal(l_窓.Get_カバレッジ(l_塩基列.AsSpan(i, C_k長)), l_束.Get_カバレッジ(l_塩基列.AsSpan(i, C_k長)));
                }
            }
        }

        /// <summary>
        /// 作業ディレクトリと実行時引数を戻す
        /// </summary>
        public void Dispose()
        {
            ConfigurationManager.A_実行時引数 = this._元の設定;
            Directory.Delete(this._作業ディレクトリ, true);
            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
