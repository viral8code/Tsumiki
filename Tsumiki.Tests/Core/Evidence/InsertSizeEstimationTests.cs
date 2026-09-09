using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドから推定される「インサートサイズ」が、リードに挟まれた内側の
    /// 未読区間ではなく、真のフラグメント長 (左リードの 5'端から右リードの
    /// 3'端まで) の単位になっていることを検証する
    /// </summary>
    /// <remarks>
    /// 実データ (150 bp リード・IS350ライブラリ) で同一 unitig 由来サンプルの
    /// 中央値が 58 と報告されていた<br/>
    /// リード長150 bp より短いフラグメントは
    /// physically ありえないため、これは単位の取り違えを示していた<br/>
    /// 内側距離 58 に両リード長を足すと 358 となりライブラリ名と一致する<br/>
    /// この取り違えはギャップ長推定 (ギャップ = インサートサイズ - 既知長) にも
    /// そのまま伝播するため、単位を明示的に固定しておく
    /// </remarks>
    public class InsertSizeEstimationTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _一時ディレクトリ;

        public InsertSizeEstimationTests()
        {
            this._一時ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_insertsize_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

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
        /// 決定的な擬似乱数で非反復的な塩基配列を作る
        /// </summary>
        /// <remarks>
        /// k=21 では
        /// この長さの乱数配列に重複 k-mer が現れる確率は無視できる
        /// </remarks>
        /// <param name="p_長さ">生成する配列長</param>
        /// <param name="p_シード">乱数シード</param>
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
                l_ライター.WriteLine(new string('I', l_配列.Length)); // Q40相当
            }
        }

        [Fact]
        public void SameUnitigSamples_MeasureFullFragmentLength_NotTheInnerDistance()
        {
            const int k = 21;
            const int unitigLength = 600;
            const int readLength = 50;
            const int fragmentStart = 100;
            const int trueFragmentLength = 350;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var unitigSeq = V_生成_乱数配列(unitigLength, p_シード: 12345);
            var unitigsPath = Path.Combine(this._一時ディレクトリ, "unitigs.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            // FR 配置: read1 はフラグメント左端から順鎖方向、
            // read2 はフラグメント右端から逆鎖方向に読まれる
            var read1 = unitigSeq.Substring(fragmentStart, readLength);
            var read2 = Util.V_逆相補(
                unitigSeq.Substring(fragmentStart + trueFragmentLength - readLength, readLength));

            var path1 = Path.Combine(this._一時ディレクトリ, "r1.fq");
            var path2 = Path.Combine(this._一時ディレクトリ, "r2.fq");
            // 中央値を安定させるため同一ペアを複数本入れる
            var pairs = Enumerable.Range(0, 5).ToList();
            V_書き込み_Fastq(path1, pairs.Select(i => ($"pair{i}/1", read1)));
            V_書き込み_Fastq(path2, pairs.Select(i => ($"pair{i}/2", read2)));

            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.V_マッピング_ペアリード(path1, path2);

            Assert.NotEmpty(contigMaker.A_同一ユニティグ標本);
            // 内側距離 (= 350 - 50 - 50 = 250) ではなく、フラグメント長 350 が
            // 得られなければならない
            Assert.All(
                contigMaker.A_同一ユニティグ標本,
                sample => Assert.Equal(trueFragmentLength, sample));
        }

        /// <summary>
        /// フラグメント長を変えたときに推定値が同じだけ動くこと (定数ぶんの
        /// ずれではなく、単位そのものが一致していること) を確認する
        /// </summary>
        [Theory]
        [InlineData(200)]
        [InlineData(350)]
        [InlineData(500)]
        public void SameUnitigSamples_TrackTheActualFragmentLength(int trueFragmentLength)
        {
            const int k = 21;
            const int unitigLength = 900;
            const int readLength = 50;
            const int fragmentStart = 120;

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            var unitigSeq = V_生成_乱数配列(unitigLength, p_シード: 777);
            var unitigsPath = Path.Combine(this._一時ディレクトリ, $"unitigs_{trueFragmentLength}.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            var read1 = unitigSeq.Substring(fragmentStart, readLength);
            var read2 = Util.V_逆相補(
                unitigSeq.Substring(fragmentStart + trueFragmentLength - readLength, readLength));

            var path1 = Path.Combine(this._一時ディレクトリ, $"r1_{trueFragmentLength}.fq");
            var path2 = Path.Combine(this._一時ディレクトリ, $"r2_{trueFragmentLength}.fq");
            V_書き込み_Fastq(path1, [("pair/1", read1)]);
            V_書き込み_Fastq(path2, [("pair/2", read2)]);

            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.V_マッピング_ペアリード(path1, path2);

            Assert.NotEmpty(contigMaker.A_同一ユニティグ標本);
            Assert.All(
                contigMaker.A_同一ユニティグ標本,
                sample => Assert.Equal(trueFragmentLength, sample));
        }
    }
}
