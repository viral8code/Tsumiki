using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// リードの索引を使う検証器が、リードを全部流して作る検証器と同じ答えを返すことの検証
    /// </summary>
    [Collection(中間データ置き場の集まり.C_名前)]
    public class RepeatRMerVerifier索引Tests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki rmer index tests
        /// </summary>
        private const string C_項目_tsumiki_rmer_index_tests = "tsumiki_rmer_index_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// k 長
        /// </summary>
        private const int C_k長 = 21;

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_rmer_index_tests + Guid.NewGuid().ToString(C_GUID書式));

        /// <summary>
        /// 元の設定
        /// </summary>
        private readonly Parameters _元の設定 = ConfigurationManager.A_実行時引数;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// RepeatRMerVerifier 索引 Tests
        /// </summary>
        public RepeatRMerVerifier索引Tests()
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
        /// Dispose
        /// </summary>
        public void Dispose()
        {
            RepeatRMerVerifier.A_Is索引使用 = false;
            RepeatRMerVerifier.V_解放_共有索引();
            ConfigurationManager.A_実行時引数 = this._元の設定;
            Directory.Delete(this._作業ディレクトリ, true);
        }

        /// <summary>
        /// 索引を使っても答えが同じ
        /// </summary>
        /// <param name="p_r長"></param>
        [Theory]
        [InlineData(31)]
        [InlineData(35)]
        [InlineData(39)]
        [InlineData(40)]
        [InlineData(41)]
        [InlineData(47)]
        [InlineData(70)]
        public void V_索引を使っても答えが同じ(int p_r長)
        {
            var l_乱数 = new Random(p_r長);
            var l_反復 = new string([..Enumerable.Range(0, 80).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_ゲノム = string.Concat(Enumerable.Range(0, 40).Select(i => i % 9 == 4 ? l_反復 : new string([..Enumerable.Range(0, 120).Select(_ => Consts.塩基文字[l_乱数.Next(4)])])));
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            using (var l_書き込み = new StreamWriter(l_パス))
            {
                for (var i = 0; i < 3_000; i++)
                {
                    var l_長さ = l_乱数.Next(60, 151);
                    var l_リード = Get_変異(l_乱数, l_ゲノム.Substring(l_乱数.Next(0, l_ゲノム.Length - l_長さ), l_長さ), 0.01D);
                    l_リード = l_乱数.Next(2) == 0 ? l_リード : Get_逆相補(l_リード);
                    l_書き込み.WriteLine($"@r{i}\n{l_リード}\n+\n{new string('I', l_リード.Length)}");
                }
            }

            using var l_信頼 = new TrustedKmerIndex(this._作業ディレクトリ);
            Cores.Preprocessing.KmerCounting.V_読込_リードファイル(l_パス, l_信頼, 33);
            l_信頼.V_適用_カットオフ(3UL);
            var l_問い合わせ = Enumerable.Range(0, 40).Select(_ =>
            {
                var l_長さ = l_乱数.Next(p_r長 + 10, 400);
                var l_配列 = Get_変異(l_乱数, l_ゲノム.Substring(l_乱数.Next(0, l_ゲノム.Length - l_長さ), l_長さ), 0.005D);
                return l_乱数.Next(2) == 0 ? l_配列 : Get_逆相補(l_配列);
            }).ToList();
            foreach (var (l_絞り込み, l_候補)in new (TrustedKmerIndex? , IEnumerable<string>? )[]
            {
                (null, null),
                (l_信頼, null),
                (null, l_問い合わせ)
            }

            )
            {
                RepeatRMerVerifier.A_Is索引使用 = false;
                var l_旧 = RepeatRMerVerifier.V_構築([l_パス], p_r長, l_絞り込み, l_絞り込み is null ? 0 : C_k長, l_候補);
                RepeatRMerVerifier.A_Is索引使用 = true;
                var l_新 = RepeatRMerVerifier.V_構築([l_パス], p_r長, l_絞り込み, l_絞り込み is null ? 0 : C_k長, l_候補);
                foreach (var l_配列 in l_問い合わせ)
                {
                    List<(int, int)> l_旧の範囲 = [];
                    List<(int, int)> l_新の範囲 = [];
                    l_旧.V_収集_未観測の連続範囲(l_配列, 1, l_旧の範囲);
                    l_新.V_収集_未観測の連続範囲(l_配列, 1, l_新の範囲);
                    Assert.Equal(l_旧の範囲, l_新の範囲);
                    var l_区切り1 = l_配列.Length / 3;
                    var l_区切り2 = 2 * l_配列.Length / 3;
                    Assert.Equal(l_旧.Get_接合点の支持数(l_配列[..(l_区切り1 + C_k長 - 1)], l_配列[l_区切り1..l_区切り2], l_配列[(l_区切り2 - C_k長 + 1)..]), l_新.Get_接合点の支持数(l_配列[..(l_区切り1 + C_k長 - 1)], l_配列[l_区切り1..l_区切り2], l_配列[(l_区切り2 - C_k長 + 1)..]));
                }
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 逆相補
        /// </summary>
        /// <param name="p_配列"></param>
        /// <returns></returns>
        private static string Get_逆相補(string p_配列) => new([..p_配列.Reverse().Select(static x => x switch
        {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => 'N'
        })]);

        /// <summary>
        /// 変異
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_率"></param>
        /// <returns></returns>
        private static string Get_変異(Random p_乱数, string p_配列, double p_率)
        {
            var l_文字 = p_配列.ToCharArray();
            for (var i = 0; i < l_文字.Length; i++)
            {
                if (p_乱数.NextDouble() < p_率)
                {
                    l_文字[i] = p_乱数.Next(20) == 0 ? 'N' : Consts.塩基文字[p_乱数.Next(4)];
                }
            }

            return new string(l_文字);
        }

        #endregion
    }
}
