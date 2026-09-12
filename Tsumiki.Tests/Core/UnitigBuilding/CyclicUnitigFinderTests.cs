using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 分岐を 1 つも持たない閉路の回収を固定する
    /// </summary>
    /// <remarks>
    /// unitig の開始点は「入次数が 1 でない、または唯一の予測元が分岐している」k-mer として選ぶため、全頂点が入次数 1 ・出次数 1 の閉路は開始点を 1 つも持たない<br/>
    /// そのままだと、きれいな環状染色体や小さなプラスミドが出力から丸ごと消える
    /// </remarks>
    public class CyclicUnitigFinderTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// k 長
        /// </summary>
        private const int k長 = 21;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public CyclicUnitigFinderTests()
        {
            this._一時ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_cycle_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 分岐のない閉路は、通常の開始点判定 (開始 k-mer 一覧) では 1 つも拾えない
        /// </summary>
        [Fact]
        public void Get_閉路の開始kmer_分岐のない閉路は通常の開始点判定では拾えない()
        {
            var l_環状 = Get_乱数配列(300, 31);
            using var l_インデックス = this.Get_インデックス(l_環状, null);

            // この前提が崩れたら、以降のテストは意味を失う
            Assert.Empty(l_インデックス.Get_開始kmer一覧());
        }

        /// <summary>
        /// 分岐のない閉路を 1 つだけ回収し、1 周ぶんの配列に復元する
        /// </summary>
        [Fact]
        public void Get_閉路の開始kmer_閉路を1つだけ回収して1周ぶんの配列にする()
        {
            var l_環状 = Get_乱数配列(300, 31);
            using var l_インデックス = this.Get_インデックス(l_環状, null);

            var l_開始kmer = CyclicUnitigFinder.Get_閉路の開始kmer(l_インデックス, UnitigMaker.Get_walk結果(l_インデックス, l_インデックス.Get_開始kmer一覧()), k長);

            _ = Assert.Single(l_開始kmer);

            var l_配列 = Assert.Single(UnitigMaker.Get_walk結果(l_インデックス, l_開始kmer));

            // 環を 1 周し、次の k-mer で出発点に戻る手前まで伸びる
            Assert.Equal(l_環状.Length + k長 - 1, l_配列.Length);
            Assert.True(Get_環の1周か(l_配列, l_環状), "環の1周になっていない");
        }

        /// <summary>
        /// 線状の配列だけなら、閉路の開始 k-mer は 1 つも返さない
        /// </summary>
        [Fact]
        public void Get_閉路の開始kmer_線状の配列だけなら何も返さない()
        {
            var l_線状 = Get_乱数配列(300, 32);
            using var l_インデックス = this.Get_インデックス(null, l_線状);

            var l_通常の走査 = UnitigMaker.Get_walk結果(l_インデックス, l_インデックス.Get_開始kmer一覧());
            Assert.NotEmpty(l_通常の走査);

            Assert.Empty(CyclicUnitigFinder.Get_閉路の開始kmer(l_インデックス, l_通常の走査, k長));
        }

        /// <summary>
        /// 線状の配列と混ざっていても、閉路の開始 k-mer だけを拾う
        /// </summary>
        [Fact]
        public void Get_閉路の開始kmer_線状の配列と混ざっていても閉路だけを拾う()
        {
            var l_環状 = Get_乱数配列(300, 33);
            var l_線状 = Get_乱数配列(300, 34);
            using var l_インデックス = this.Get_インデックス(l_環状, l_線状);

            var l_通常の走査 = UnitigMaker.Get_walk結果(l_インデックス, l_インデックス.Get_開始kmer一覧());

            var l_開始kmer = CyclicUnitigFinder.Get_閉路の開始kmer(l_インデックス, l_通常の走査, k長);

            _ = Assert.Single(l_開始kmer);
            var l_配列 = Assert.Single(UnitigMaker.Get_walk結果(l_インデックス, l_開始kmer));
            Assert.True(Get_環の1周か(l_配列, l_環状), "拾ったのが環ではない");
        }

        /// <summary>
        /// 独立した閉路が 2 つあれば、それぞれ 1 つずつ開始 k-mer を返す
        /// </summary>
        [Fact]
        public void Get_閉路の開始kmer_閉路が2つあればそれぞれ1つずつ返す()
        {
            var l_環状1 = Get_乱数配列(300, 35);
            var l_環状2 = Get_乱数配列(250, 36);

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            using var l_インデックス = new TrustedKmerIndex(this._一時ディレクトリ);
            foreach (var l_環状 in new[] { l_環状1, l_環状2 })
            {
                var l_二周 = (l_環状 + l_環状).Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i < l_環状.Length; i++)
                {
                    V_登録(l_インデックス, l_二周.AsSpan(i, k長));
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            var l_開始kmer = CyclicUnitigFinder.Get_閉路の開始kmer(l_インデックス, [], k長);

            Assert.Equal(2, l_開始kmer.Count);
            var l_長さ = UnitigMaker.Get_walk結果(l_インデックス, l_開始kmer)
                .Select(x => x.Length).OrderBy(x => x).ToList();
            Assert.Equal([l_環状2.Length + k長 - 1, l_環状1.Length + k長 - 1], l_長さ);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string l_塩基 = "ACGT";
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => l_塩基[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 環状配列と線状配列から k-mer インデックスを作る
        /// </summary>
        /// <remarks>
        /// 環状側は末尾から先頭へ回り込む窓まで登録し、閉路そのものにする
        /// </remarks>
        /// <param name="p_環状配列"></param>
        /// <param name="p_線状配列"></param>
        /// <returns></returns>
        private TrustedKmerIndex Get_インデックス(string? p_環状配列, string? p_線状配列)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._一時ディレクトリ);

            if (p_環状配列 is { } l_環状)
            {
                var l_二周 = (l_環状 + l_環状).Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i < l_環状.Length; i++)
                {
                    V_登録(l_インデックス, l_二周.AsSpan(i, k長));
                }
            }
            if (p_線状配列 is { } l_線状)
            {
                var l_塩基列 = l_線状.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + k長 <= l_塩基列.Length; i++)
                {
                    V_登録(l_インデックス, l_塩基列.AsSpan(i, k長));
                }
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_インデックス">登録先の索引</param>
        /// <param name="p_kmer">登録する k-mer</param>
        private static void V_登録(TrustedKmerIndex p_インデックス, Span<byte> p_kmer)
        {
            for (var l_回 = 0; l_回 < 3; l_回++)
            {
                p_インデックス.V_登録(p_kmer);
            }
        }

        /// <summary>
        /// 配列が環状配列の回転 (順鎖・逆鎖のいずれか) になっているか
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_環状配列"></param>
        /// <returns></returns>
        private static bool Get_環の1周か(string p_配列, string p_環状配列)
        {
            var l_二周 = p_環状配列 + p_環状配列;
            return l_二周.Contains(p_配列, StringComparison.Ordinal)
                || l_二周.Contains(Util.V_逆相補(p_配列), StringComparison.Ordinal);
        }

        #endregion

    }
}
