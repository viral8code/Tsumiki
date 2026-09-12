using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドから推定される「インサートサイズ」が、リードに挟まれた内側の未読区間ではなく、真のフラグメント長 (左リードの 5'端から右リードの 3'端まで) の単位になっていることを検証する
    /// </summary>
    /// <remarks>
    /// 実データ (150 bp リード・ IS350 ライブラリ) で同一 unitig 由来サンプルの中央値が 58 と報告されていた<br/>
    /// リード長 150 bp より短いフラグメントは physically ありえないため、これは単位の取り違えを示していた<br/>
    /// 内側距離 58 に両リード長を足すと 358 となりライブラリ名と一致する<br/>
    /// この取り違えはギャップ長推定 (ギャップ = インサートサイズ - 既知長) にもそのまま伝播するため、単位を明示的に固定しておく
    /// </remarks>
    public class InsertSizeEstimationTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _一時ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public InsertSizeEstimationTests()
        {
            this._一時ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_insertsize_tests_" + Guid.NewGuid().ToString("N"));
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
        }

        /// <summary>
        /// V_SameUnitigSamples_MeasureFullFragmentLength_NotTheInnerDistance
        /// </summary>
        [Fact]
        public void V_SameUnitigSamples_MeasureFullFragmentLength_NotTheInnerDistance()
        {
            const int l_k長 = 21;
            const int l_ユニティグ長 = 600;
            const int l_リード長 = 50;
            const int l_断片開始位置 = 100;
            const int l_真の断片長 = 350;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };

            var l_ユニティグ配列 = V_生成_乱数配列(l_ユニティグ長, p_シード: 12345);
            var l_ユニティグパス = Path.Combine(this._一時ディレクトリ, "unitigs.fasta");
            File.WriteAllText(l_ユニティグパス, $">1\n{l_ユニティグ配列}\n");

            // FR 配置: read1 はフラグメント左端から順鎖方向、
            // read2 はフラグメント右端から逆鎖方向に読まれる
            var l_先行リード = l_ユニティグ配列.Substring(l_断片開始位置, l_リード長);
            var l_後続リード = Util.V_逆相補(l_ユニティグ配列.Substring(l_断片開始位置 + l_真の断片長 - l_リード長, l_リード長));

            var l_先行パス = Path.Combine(this._一時ディレクトリ, "r1.fq");
            var l_後続パス = Path.Combine(this._一時ディレクトリ, "r2.fq");
            // 中央値を安定させるため同一ペアを複数本入れる
            var l_ペア群 = Enumerable.Range(0, 5).ToList();
            V_書き込み_Fastq(l_先行パス, l_ペア群.Select(i => ($"pair{i}/1", l_先行リード)));
            V_書き込み_Fastq(l_後続パス, l_ペア群.Select(i => ($"pair{i}/2", l_後続リード)));

            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);
            l_コンティグ構築.V_マッピング_ペアリード(l_先行パス, l_後続パス);

            Assert.NotEmpty(l_コンティグ構築.A_同一ユニティグ標本);
            // 内側距離 (= 350 - 50 - 50 = 250) ではなく、フラグメント長 350 が
            // 得られなければならない
            Assert.All(l_コンティグ構築.A_同一ユニティグ標本, l_標本 => Assert.Equal(l_真の断片長, l_標本));
        }

        /// <summary>
        /// フラグメント長を変えたときに推定値が同じだけ動くこと (定数ぶんのずれではなく、単位そのものが一致していること) を確認する
        /// </summary>
        /// <param name="p_真の断片長"></param>
        [Theory]
        [InlineData(200)]
        [InlineData(350)]
        [InlineData(500)]
        public void V_SameUnitigSamples_TrackTheActualFragmentLength(int p_真の断片長)
        {
            const int l_k長 = 21;
            const int l_ユニティグ長 = 900;
            const int l_リード長 = 50;
            const int l_断片開始位置 = 120;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };

            var l_ユニティグ配列 = V_生成_乱数配列(l_ユニティグ長, p_シード: 777);
            var l_ユニティグパス = Path.Combine(this._一時ディレクトリ, $"unitigs_{p_真の断片長}.fasta");
            File.WriteAllText(l_ユニティグパス, $">1\n{l_ユニティグ配列}\n");

            var l_先行リード = l_ユニティグ配列.Substring(l_断片開始位置, l_リード長);
            var l_後続リード = Util.V_逆相補(l_ユニティグ配列.Substring(l_断片開始位置 + p_真の断片長 - l_リード長, l_リード長));

            var l_先行パス = Path.Combine(this._一時ディレクトリ, $"r1_{p_真の断片長}.fq");
            var l_後続パス = Path.Combine(this._一時ディレクトリ, $"r2_{p_真の断片長}.fq");
            V_書き込み_Fastq(l_先行パス, [("pair/1", l_先行リード)]);
            V_書き込み_Fastq(l_後続パス, [("pair/2", l_後続リード)]);

            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);
            l_コンティグ構築.V_マッピング_ペアリード(l_先行パス, l_後続パス);

            Assert.NotEmpty(l_コンティグ構築.A_同一ユニティグ標本);
            Assert.All(l_コンティグ構築.A_同一ユニティグ標本, l_標本 => Assert.Equal(p_真の断片長, l_標本));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 決定的な擬似乱数で非反復的な塩基配列を作る
        /// </summary>
        /// <remarks>
        /// k=21 ではこの長さの乱数配列に重複 k-mer が現れる確率は無視できる
        /// </remarks>
        /// <param name="p_長さ">生成する配列長</param>
        /// <param name="p_シード">乱数シード</param>
        /// <returns></returns>
        private static string V_生成_乱数配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        /// <param name="p_リード一覧">書き出すリード</param>
        private static void V_書き込み_Fastq(string p_パス, IEnumerable<(string A_ID, string A_配列)> p_リード一覧)
        {
            using var l_ライター = new StreamWriter(p_パス);
            foreach (var (l_ID, l_配列) in p_リード一覧)
            {
                l_ライター.WriteLine($"@{l_ID}");
                l_ライター.WriteLine(l_配列);
                l_ライター.WriteLine("+");
                l_ライター.WriteLine(new string('I', l_配列.Length)); // Q40 相当
            }
        }

        #endregion

    }
}
