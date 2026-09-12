using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// リード長が k より短いリードが混ざっていても処理が破綻しないことを固定する
    /// </summary>
    /// <remarks>
    /// トリミング済みのデータではリード長がばらつく<br/>
    /// GAGE-B の R. sphaeroides MiSeq (trimmed) では 755,847 本のうち 8% 以上が k=63 未満で、最短は 19 bp だった<br/>
    /// マッピング側に長さの判定が無く、19 bp のリードに対して添字 62 までアクセスして例外になっていた<br/>
    /// しかもその例外はワーカースレッドの中で起き、キューが満杯になったプロデューサーが永久に待ち続けたため、ログも例外も出ないまま 2 時間以上プロセスが停止した<br/>
    /// 長さの判定と、ワーカーの例外を伝える仕組み (ReadPipelineTests) の両方が要る
    /// </remarks>
    public class VariableLengthReadTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 31;

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
        public VariableLengthReadTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_varlen_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
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
        /// k より短いリードと十分長いリードが混ざったペアエンド入力
        /// </summary>
        /// <remarks>
        /// 短いリードは黙って読み飛ばされ、長いリード由来の隣接だけが残ること
        /// </remarks>
        [Fact]
        public void V_k未満のペアリードは失敗せずに読み飛ばされる()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };

            var l_ユニティグ配列 = V_生成_ランダム配列(600, p_シード: 987);
            var l_ユニティグパス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            File.WriteAllText(l_ユニティグパス, $">1\n{l_ユニティグ配列}\n");

            // 19 bp (最短の実例と同じ長さ) から 200 bp まで、k をまたぐ長さを混ぜる
            var l_長さ一覧 = new[] { 19, 30, k長 - 1, k長, k長 + 1, 120, 200 };
            var l_先行リード群 = new List<(string, string)>();
            var l_後続リード群 = new List<(string, string)>();
            for (var i = 0; i < l_長さ一覧.Length; i++)
            {
                var l_長さ = l_長さ一覧[i];
                l_先行リード群.Add(($"pair{i}/1", l_ユニティグ配列[..l_長さ]));
                l_後続リード群.Add(($"pair{i}/2", Util.V_逆相補(l_ユニティグ配列[^l_長さ..])));
            }

            var l_先行パス = this.V_書き出し_FASTQ("short.1.fq", l_先行リード群);
            var l_後続パス = this.V_書き出し_FASTQ("short.2.fq", l_後続リード群);

            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);

            // 例外を投げずに完走すること
            // 対策前はここで
            // IndexOutOfRangeException がワーカー内で起き、
            // そのままハングしていた
            l_コンティグ構築.V_マッピング_ペアリード(l_先行パス, l_後続パス);

            // k 以上のリードからは標本が取れていること
            // (短いリードのせいで全部落ちてしまっていないことの確認)
            Assert.NotEmpty(l_コンティグ構築.A_インサートサイズ標本);
        }

        /// <summary>
        /// k 未満のリードが混ざっていても、失敗せずに読み飛ばされる
        /// </summary>
        [Fact]
        public void V_k未満の単一リードは失敗せずに読み飛ばされる()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };

            var l_ユニティグ配列 = V_生成_ランダム配列(400, p_シード: 654);
            var l_ユニティグパス = Path.Combine(this._作業ディレクトリ, "unitigs_single.fasta");
            File.WriteAllText(l_ユニティグパス, $">1\n{l_ユニティグ配列}\n");

            var l_リード群 = new List<(string, string)>();
            for (var i = 0; i < 50; i++)
            {
                // 半分を k 未満にする
                var l_長さ = i % 2 == 0 ? 19 : 150;
                l_リード群.Add(($"read{i}", l_ユニティグ配列[..l_長さ]));
            }
            var l_パス = this.V_書き出し_FASTQ("short_single.fq", l_リード群);

            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);
            l_コンティグ構築.V_マッピング_リード(l_パス);
        }

        /// <summary>
        /// すべてのリードが k 未満でも、例外にならず単に何も得られないこと
        /// </summary>
        [Fact]
        public void V_すべてのリードがk未満のときは例外にならず標本が空になる()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };

            var l_ユニティグ配列 = V_生成_ランダム配列(400, p_シード: 321);
            var l_ユニティグパス = Path.Combine(this._作業ディレクトリ, "unitigs_allshort.fasta");
            File.WriteAllText(l_ユニティグパス, $">1\n{l_ユニティグ配列}\n");

            var l_先行リード群 = Enumerable.Range(0, 40).Select(i => ($"pair{i}/1", l_ユニティグ配列[..19]));
            var l_後続リード群 = Enumerable.Range(0, 40).Select(i => ($"pair{i}/2", l_ユニティグ配列[..20]));

            var l_先行パス = this.V_書き出し_FASTQ("allshort.1.fq", l_先行リード群);
            var l_後続パス = this.V_書き出し_FASTQ("allshort.2.fq", l_後続リード群);

            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);
            l_コンティグ構築.V_マッピング_ペアリード(l_先行パス, l_後続パス);

            Assert.Empty(l_コンティグ構築.A_インサートサイズ標本);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_ランダム配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_リード列">書き出すリード</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_FASTQ(string p_ファイル名, IEnumerable<(string A_ID, string A_配列)> p_リード列)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new StreamWriter(l_パス);
            foreach (var (l_ID, l_配列) in p_リード列)
            {
                l_書き込み.WriteLine($"@{l_ID}");
                l_書き込み.WriteLine(l_配列);
                l_書き込み.WriteLine("+");
                l_書き込み.WriteLine(new string('I', l_配列.Length)); // Q40 相当
            }
            return l_パス;
        }

        #endregion

    }
}
