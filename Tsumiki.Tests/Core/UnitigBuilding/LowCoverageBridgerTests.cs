using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// カットオフ未満で控えた k-mer による行き止まりの架橋の検証
    /// </summary>
    public class LowCoverageBridgerTests : IDisposable
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
        /// カバレッジが落ち込む k-mer の開始位置の範囲 (両端を含む)
        /// </summary>
        private const int 谷の先頭 = 14;

        /// <summary>
        /// カバレッジが落ち込む k-mer の開始位置の範囲の終端 (含む)
        /// </summary>
        private const int 谷の末尾 = 22;

        /// <summary>
        /// 主経路の深さ
        /// </summary>
        private const int 主経路の深さ = 20;

        /// <summary>
        /// 谷の深さ
        /// </summary>
        private const int 谷の深さ = 4;

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
        public LowCoverageBridgerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_bridger_tests_" + Guid.NewGuid().ToString("N"));
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
        }

        /// <summary>
        /// カットオフを割る谷で途切れた実配列は、控えの k-mer で繋ぎ直される
        /// </summary>
        [Fact]
        public void V_カバレッジの谷で途切れた配列は架橋される()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_谷つき(l_索引, 主配列);
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);
            Assert.False(l_索引.Haskmer(Get_kmer(主配列, 谷の先頭)));

            var l_追加数 = LowCoverageBridger.Get_架橋kmer数(l_索引, k長, null);

            Assert.Equal(谷の末尾 - 谷の先頭 + 1, l_追加数);
            for (var i = 0; i + k長 <= 主配列.Length; i++)
            {
                Assert.True(l_索引.Haskmer(Get_kmer(主配列, i)), $"位置 {i} の k-mer が集合に無い");
            }
            Assert.Equal(0, l_索引.A_控えkmer数);
        }

        /// <summary>
        /// 行き止まりの先が控えの k-mer だけで終わり、信頼できる k-mer に合流しない枝は足さない
        /// </summary>
        [Fact]
        public void V_合流しない低カバレッジの続きは足さない()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_索引, 主配列, 主経路の深さ);
            V_登録_全kmer(l_索引, 主配列[^(k長 - 1)..] + "CGTCGA", 谷の深さ);
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);

            var l_追加数 = LowCoverageBridger.Get_架橋kmer数(l_索引, k長, null);

            Assert.Equal(0, l_追加数);
        }

        /// <summary>
        /// 同じくらいの出現回数の続きが 2 つあり、どちら向きからも一本道に決まらない谷は繋がない
        /// </summary>
        [Fact]
        public void V_続きが一本に決まらない谷は架橋しない()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_谷つき(l_索引, 主配列);

            // 谷の中の 1 塩基だけ違う対立配列を同じ深さで足し、谷の両端から見て分岐にする
            var l_対立配列 = 主配列[..21] + "G" + 主配列[22..];
            for (var i = 谷の先頭; i <= 21; i++)
            {
                V_登録(l_索引, Get_kmer(l_対立配列, i), 谷の深さ);
            }
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);

            var l_追加数 = LowCoverageBridger.Get_架橋kmer数(l_索引, k長, null);

            Assert.Equal(0, l_追加数);
        }

        /// <summary>
        /// 合流先の信頼できる k-mer に既に別の入口があるなら架橋しない
        /// </summary>
        /// <remarks>
        /// 薄いカバレッジで途切れた箇所なら合流先も行き止まりになる<br/>
        /// 別の入口があるのは、行き止まりから既存の配列へ新しく入る分岐で、短い反復を挟んだ近道にもなりうる
        /// </remarks>
        [Fact]
        public void V_合流先に別の入口があれば架橋しない()
        {
            using var l_索引 = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_谷つき(l_索引, 主配列);

            // 谷を抜けた最初の k-mer (位置 23) へ、谷とは別の塩基から入る信頼できる枝を足す
            var l_合流先 = 主配列.Substring(谷の末尾 + 1, k長);
            var l_枝 = "TGTTTGCA" + "G" + l_合流先[..^1];
            Assert.NotEqual(主配列[谷の末尾], l_枝[^k長]);
            V_登録_全kmer(l_索引, l_枝, 主経路の深さ);
            l_索引.V_適用_カットオフ(カットオフ, LowCoverageBridger.控えの最小出現回数);

            var l_追加数 = LowCoverageBridger.Get_架橋kmer数(l_索引, k長, null);

            Assert.Equal(0, l_追加数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 谷の範囲だけ浅く、それ以外を主経路の深さで全 k-mer を登録する
        /// </summary>
        /// <param name="p_索引"></param>
        /// <param name="p_配列"></param>
        private static void V_登録_谷つき(TrustedKmerIndex p_索引, string p_配列)
        {
            for (var i = 0; i + k長 <= p_配列.Length; i++)
            {
                V_登録(p_索引, Get_kmer(p_配列, i), i is >= 谷の先頭 and <= 谷の末尾 ? 谷の深さ : 主経路の深さ);
            }
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
