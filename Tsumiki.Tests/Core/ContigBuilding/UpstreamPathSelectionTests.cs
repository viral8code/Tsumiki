using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 多コピーと推定された分岐元の行き先を、手前の単一コピー unitig から通り抜けたリードで決める検証
    /// </summary>
    public class UpstreamPathSelectionTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 8;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        /// <summary>
        /// 手前の単一コピー unitig
        /// </summary>
        private readonly string _U;

        /// <summary>
        /// 多コピーと推定される分岐元
        /// </summary>
        private readonly string _V;

        /// <summary>
        /// 分岐先の片方
        /// </summary>
        private readonly string _W1;

        /// <summary>
        /// 分岐先のもう片方
        /// </summary>
        private readonly string _W2;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// U → V → (W1 | W2) の構成を用意する
        /// </summary>
        public UpstreamPathSelectionTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_upstream_path_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };

            this._V = RepeatGraphFixture.Get_乱数配列(16, 4001);
            // 断片長の標本が無いとき起点に要る長さは k の 4 倍なので、それを超える長さにする
            this._U = RepeatGraphFixture.Get_乱数配列(30, 4002) + this._V[..(k長 - 1)];
            this._W1 = this._V[^(k長 - 1)..] + "G" + RepeatGraphFixture.Get_乱数配列(11, 4003);
            this._W2 = this._V[^(k長 - 1)..] + "T" + RepeatGraphFixture.Get_乱数配列(11, 4004);
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
        /// U から V を通って W1 へ抜けたリードがあれば、V が多コピーと推定されていても V → W1 を繋ぐ
        /// </summary>
        [Fact]
        public void V_手前の単一コピーから通り抜けたリードで多コピー推定の分岐を決める()
        {
            var l_UVW1 = this._U + this._V[(k長 - 1)..] + this._W1[(k長 - 1)..];
            var l_contig群 = this.Get_contig群([.. Get_窓群(l_UVW1, 30, this._U.Length - k長, this._U.Length + 16 - (2 * (k長 - 1)) + k長), .. this.Get_VW2の窓群()]);

            Assert.Equal(2, l_contig群.Count);
            Assert.Contains(l_contig群, x => x.Length == l_UVW1.Length);
            Assert.Contains(l_contig群, x => x.Length == this._W2.Length);
        }

        /// <summary>
        /// 通り抜けたリードが U に届いていなければ、多コピー推定の分岐はこれまでどおり決めない
        /// </summary>
        [Fact]
        public void V_手前の単一コピーに届かないリードでは多コピー推定の分岐を決めない()
        {
            var l_VW1 = this._V + this._W1[(k長 - 1)..];
            var l_contig群 = this.Get_contig群([.. Get_窓群(l_VW1, 20, 0, 17), .. this.Get_VW2の窓群()]);

            Assert.Equal(3, l_contig群.Count);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// V と W2 だけを跨ぐ読み取り窓
        /// </summary>
        /// <returns></returns>
        private List<string> Get_VW2の窓群()
        {
            return Get_窓群(this._V + this._W2[(k長 - 1)..], 20, 0, 17);
        }

        /// <summary>
        /// 開始位置が上限以下で、終了位置が下限以上の窓
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_長さ"></param>
        /// <param name="p_開始の上限"></param>
        /// <param name="p_終了の下限"></param>
        /// <returns></returns>
        private static List<string> Get_窓群(string p_配列, int p_長さ, int p_開始の上限, int p_終了の下限)
        {
            List<string> l_窓群 = [];
            for (var s = 0; s + p_長さ <= p_配列.Length; s++)
            {
                if (s <= p_開始の上限 && s + p_長さ >= p_終了の下限)
                {
                    l_窓群.Add(p_配列.Substring(s, p_長さ));
                }
            }
            return l_窓群;
        }

        /// <summary>
        /// V を 2 コピーと推定させて contig を組み、出力された配列を返す
        /// </summary>
        /// <param name="p_リード群"></param>
        /// <returns></returns>
        private List<string> Get_contig群(IReadOnlyList<string> p_リード群)
        {
            var l_unitigパス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            using (var l_書き込み = new FastaWriter(l_unitigパス))
            {
                l_書き込み.V_書き込み(1, this._U);
                l_書き込み.V_書き込み(2, this._V);
                l_書き込み.V_書き込み(3, this._W1);
                l_書き込み.V_書き込み(4, this._W2);
            }

            var l_リードパス = Path.Combine(this._作業ディレクトリ, "reads.fq");
            using (var l_書き込み = new StreamWriter(l_リードパス))
            {
                for (var i = 0; i < p_リード群.Count; i++)
                {
                    l_書き込み.WriteLine($"@r{i}\n{p_リード群[i]}\n+\n{new string('I', p_リード群[i].Length)}");
                }
            }

            var l_contig構築 = new ContigMaker(l_unitigパス);
            var l_グラフ = l_contig構築.Get_グラフ();
            Assert.Equal(1, l_グラフ.Get_入次数(ContigMaker.Get_頂点番号(2)));
            Assert.Equal(2, l_グラフ.A_出辺[ContigMaker.Get_頂点番号(2)].Count);

            l_contig構築.V_マッピング_リード(l_リードパス);
            var l_contigパス = Path.Combine(this._作業ディレクトリ, "contigs.fasta");
            Dictionary<int, int> l_コピー数 = new() { [1] = 1, [2] = 2, [3] = 1, [4] = 1 };
            l_contig構築.V_結合_Contig(l_contigパス, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_コピー数: l_コピー数);
            return [.. FastaReader.Get_全エントリ(l_contigパス).Select(x => x.A_配列)];
        }

        #endregion
    }
}
