using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// GapFiller が埋められなかったスキャフォールドのギャップを、その両端に
    /// 実際にマップされた局所リードだけで再アセンブリして埋める処理の検証
    /// </summary>
    /// <remarks>
    /// AssemblyMerger(-mg) と違い、他の k の「既に確定した結論」を持ち込むの
    /// ではなく、生リードから新しく証拠を集める<br/>
    /// ここでは生の FASTQ を
    /// 直接与えて、局所アセンブリだけでギャップが埋まる/埋まらないことを
    /// 検証する (GapFiller 側は経由しない)
    /// </remarks>
    public class LocalAssemblerTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public LocalAssemblerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_localasm_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
            // LocalAssembler は内部で TrustedKmerIndex/KmerKey を使うため、
            // 現在の実行時引数の k 長を、テストで使う k (21) に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 1 };
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

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
        /// スキャフォールドを FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_スキャフォールド(string p_ファイル名, string p_配列)
        {
            var l_パス = Path.Combine(this._tempDir, p_ファイル名);
            using var l_writer = new FastaWriter(l_パス);
            l_writer.V_書き込み("SCAFFOLD1", p_配列);
            return l_パス;
        }

        /// <summary>
        /// ゲノムを覆うリードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_ゲノム群">元になるゲノム</param>
        /// <param name="p_リード長">リード長</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_リード(string p_ファイル名, IEnumerable<string> p_ゲノム群, int p_リード長)
        {
            var l_パス = Path.Combine(this._tempDir, p_ファイル名);
            using var l_writer = new StreamWriter(l_パス);
            var l_連番 = 0;
            foreach (var l_ゲノム in p_ゲノム群)
            {
                for (var i = 0; i + p_リード長 <= l_ゲノム.Length; i++)
                {
                    l_writer.WriteLine($"@r{l_連番}");
                    l_writer.WriteLine(l_ゲノム.Substring(i, p_リード長));
                    l_writer.WriteLine("+");
                    l_writer.WriteLine(new string('I', p_リード長));
                    l_連番++;
                }
            }
            return l_パス;
        }

        /// <summary>
        /// FASTA から 1 本だけの配列を読み込む
        /// </summary>
        /// <param name="p_パス">読み込むパス</param>
        /// <returns>読み込んだ配列</returns>
        private static string Get_単一配列(string p_パス)
        {
            using var l_reader = new FastaReader(p_パス);
            Assert.True(l_reader.Get_続きがあるか());
            return l_reader.Get_次の配列().A_配列;
        }

        /// <summary>
        /// リードがギャップを実際に跨いでいるとき、真の配列を復元することを検証する
        /// </summary>
        [Fact]
        public void リードがギャップを実際に跨ぐときは真の配列を復元する()
        {
            const int k = 21;
            var prefix = V_生成_ランダム配列(300, p_シード: 1);
            var fill = V_生成_ランダム配列(60, p_シード: 2);
            var suffix = V_生成_ランダム配列(300, p_シード: 3);
            var truth = prefix + fill + suffix;

            var scaffoldPath = this.V_書き出し_スキャフォールド("scaffold.fasta", prefix + new string('N', fill.Length) + suffix);
            var readsPath = this.V_書き出し_リード("reads.fq", [truth], p_リード長: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(1, stats.A_埋めたギャップ数);
            Assert.Equal(fill.Length, stats.A_埋めた塩基数);
            Assert.Equal(truth, Get_単一配列(scaffoldPath));
        }

        /// <summary>
        /// どちらのアンカー付近にも局所リードが集まらないとき、ギャップを N のまま残すことを検証する
        /// </summary>
        [Fact]
        public void どちらのアンカー付近にも局所リードが無いときはギャップをNのまま残す()
        {
            const int k = 21;
            var prefix = V_生成_ランダム配列(300, p_シード: 10);
            var fill = V_生成_ランダム配列(60, p_シード: 11);
            var suffix = V_生成_ランダム配列(300, p_シード: 12);

            var scaffoldPath = this.V_書き出し_スキャフォールド("scaffold_noreads.fasta", prefix + new string('N', fill.Length) + suffix);
            // まったく無関係な配列からリードを取る
            var unrelated = V_生成_ランダム配列(500, p_シード: 999);
            var readsPath = this.V_書き出し_リード("reads_unrelated.fq", [unrelated], p_リード長: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(0, stats.A_埋めたギャップ数);
            Assert.Equal(1, stats.A_局所リードが集まらなかった数);
            Assert.Contains('N', Get_単一配列(scaffoldPath));
        }

        /// <summary>
        /// アンカー付近のリードはあるが、両者を橋渡しする配列 (ギャップの中身) を
        /// 読んだリードが無い場合、経路が繋がらないので埋められない
        /// </summary>
        [Fact]
        public void 局所リードがギャップを橋渡ししないときは経路無しとしてNのまま残す()
        {
            const int k = 21;
            var prefix = V_生成_ランダム配列(300, p_シード: 20);
            var fill = V_生成_ランダム配列(60, p_シード: 21);
            var suffix = V_生成_ランダム配列(300, p_シード: 22);

            var scaffoldPath = this.V_書き出し_スキャフォールド("scaffold_nopath.fasta", prefix + new string('N', fill.Length) + suffix);
            // prefix と suffix それぞれの内部だけを読んだリード (橋渡しは無い)
            var readsPath = this.V_書き出し_リード("reads_nopath.fq", [prefix, suffix], p_リード長: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(0, stats.A_埋めたギャップ数);
            Assert.Equal(1, stats.A_到達できなかった数);
            Assert.Contains('N', Get_単一配列(scaffoldPath));
        }

        /// <summary>
        /// 橋渡しの配列が 2 通りとも読まれている場合、どちらが正しいか
        /// 決められないので N のまま残す (誤った配列で埋めるより安全)
        /// </summary>
        [Fact]
        public void 橋渡しの配列が2通り同程度に支持されるときは推測せずNのまま残す()
        {
            const int k = 21;
            var prefix = V_生成_ランダム配列(300, p_シード: 30);
            var suffix = V_生成_ランダム配列(300, p_シード: 31);
            var fillA = V_生成_ランダム配列(60, p_シード: 32);
            var fillB = V_生成_ランダム配列(60, p_シード: 33);

            var scaffoldPath = this.V_書き出し_スキャフォールド("scaffold_ambiguous.fasta", prefix + new string('N', fillA.Length) + suffix);
            var readsPath = this.V_書き出し_リード(
                "reads_ambiguous.fq", [prefix + fillA + suffix, prefix + fillB + suffix], p_リード長: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(0, stats.A_埋めたギャップ数);
            Assert.Equal(1, stats.A_一意に定まらなかった数);
            Assert.Contains('N', Get_単一配列(scaffoldPath));
        }

        /// <summary>
        /// スキャフォールドにギャップが無いときは、対象ギャップ数をゼロとして報告することを検証する
        /// </summary>
        [Fact]
        public void スキャフォールドにギャップが無いときは対象ギャップ数をゼロとして報告する()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(300, p_シード: 40);
            var scaffoldPath = this.V_書き出し_スキャフォールド("scaffold_nogap.fasta", truth);
            var readsPath = this.V_書き出し_リード("reads_nogap.fq", [truth], p_リード長: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(0, stats.A_対象ギャップ数);
            Assert.Equal(truth, Get_単一配列(scaffoldPath));
        }
    }
}
