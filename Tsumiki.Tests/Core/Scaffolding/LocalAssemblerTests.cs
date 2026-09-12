using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// GapFiller が埋められなかった scaffold のギャップを、その両端に実際にマップされた局所リードだけで再アセンブリして埋める処理の検証
    /// </summary>
    /// <remarks>
    /// AssemblyMerger (-mg) と違い、他の k の「既に確定した結論」を持ち込むのではなく、生リードから新しく証拠を集める<br/>
    /// ここでは生の FASTQ を直接与えて、局所アセンブリだけでギャップが埋まる/埋まらないことを検証する (GapFiller 側は経由しない)
    /// </remarks>
    public class LocalAssemblerTests : IDisposable
    {
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
        public LocalAssemblerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_localasm_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            // LocalAssembler は内部で TrustedKmerIndex/KmerKey を使うため、
            // 現在の実行時引数の k 長を、テストで使う k (21) に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 1 };
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 短い文脈で曖昧な反復を長いリード文脈で識別する
        /// </summary>
        [Fact]
        public void V_局所文脈_短い反復を識別()
        {
            var l_左 = V_生成_ランダム配列(300, p_シード: 91);
            var l_右 = V_生成_ランダム配列(300, p_シード: 92);
            var l_反復 = V_生成_ランダム配列(25, p_シード: 93);
            var l_間 = V_生成_ランダム配列(5, p_シード: 94);
            var l_充填 = l_反復 + l_間 + l_反復;
            var l_正解 = l_左 + l_充填 + l_右;
            var l_パス = this.V_書き出し_スキャフォールド("adaptive.fasta", l_左 + new string('N', l_充填.Length) + l_右);
            var l_リード = this.V_書き出し_リード("adaptive.fq", [l_正解], p_リード長: 150);
            var l_固定結果 = LocalAssembler.Get_固定kの局所結果(new 局所ギャップ(0, 300, l_充填.Length, l_左, l_右), FastqReader.Get_生リード列(l_リード, null).ToList(), 21, out var l_固定判定);
            Assert.Null(l_固定結果);
            Assert.Equal(ギャップ充填判定.一意でない, l_固定判定);
            var l_結果 = LocalAssembler.V_充填_ギャップ(l_パス, l_リード, string.Empty, 21);
            Assert.Equal(1, l_結果.A_埋めたギャップ数);
            Assert.Equal(l_正解, Get_単一配列(l_パス));
        }

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
        /// リードがギャップを実際に跨いでいるとき、真の配列を復元することを検証する
        /// </summary>
        [Fact]
        public void V_リードがギャップを実際に跨ぐときは真の配列を復元する()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 1);
            var l_充填配列 = V_生成_ランダム配列(60, p_シード: 2);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 3);
            var l_正解 = l_接頭配列 + l_充填配列 + l_接尾配列;

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold.fasta", l_接頭配列 + new string('N', l_充填配列.Length) + l_接尾配列);
            var l_リードパス = this.V_書き出し_リード("reads.fq", [l_正解], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(1, l_統計.A_埋めたギャップ数);
            Assert.Equal(l_充填配列.Length, l_統計.A_埋めた塩基数);
            Assert.Equal(l_正解, Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// どちらのアンカー付近にも局所リードが集まらないとき、ギャップを N のまま残すことを検証する
        /// </summary>
        [Fact]
        public void V_どちらのアンカー付近にも局所リードが無いときはギャップをNのまま残す()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 10);
            var l_充填配列 = V_生成_ランダム配列(60, p_シード: 11);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 12);

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_noreads.fasta", l_接頭配列 + new string('N', l_充填配列.Length) + l_接尾配列);
            // まったく無関係な配列からリードを取る
            var l_無関係配列 = V_生成_ランダム配列(500, p_シード: 999);
            var l_リードパス = this.V_書き出し_リード("reads_unrelated.fq", [l_無関係配列], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_局所リードが集まらなかった数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// アンカー付近のリードはあるが、両者を橋渡しする配列 (ギャップの中身) を読んだリードが無い場合、経路が繋がらないので埋められない
        /// </summary>
        [Fact]
        public void V_局所リードがギャップを橋渡ししないときは経路無しとしてNのまま残す()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 20);
            var l_充填配列 = V_生成_ランダム配列(60, p_シード: 21);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 22);

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_nopath.fasta", l_接頭配列 + new string('N', l_充填配列.Length) + l_接尾配列);
            // prefix と suffix それぞれの内部だけを読んだリード (橋渡しは無い)
            var l_リードパス = this.V_書き出し_リード("reads_nopath.fq", [l_接頭配列, l_接尾配列], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_到達できなかった数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// 橋渡しの配列が 2 通りとも読まれている場合、どちらが正しいか決められないので N のまま残す (誤った配列で埋めるより安全)
        /// </summary>
        [Fact]
        public void V_橋渡しの配列が2通り同程度に支持されるときは推測せずNのまま残す()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 30);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 31);
            var l_主充填配列 = V_生成_ランダム配列(60, p_シード: 32);
            var l_代替充填配列 = V_生成_ランダム配列(60, p_シード: 33);

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_ambiguous.fasta", l_接頭配列 + new string('N', l_主充填配列.Length) + l_接尾配列);
            var l_リードパス = this.V_書き出し_リード("reads_ambiguous.fq", [l_接頭配列 + l_主充填配列 + l_接尾配列, l_接頭配列 + l_代替充填配列 + l_接尾配列], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_一意に定まらなかった数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// scaffold にギャップが無いときは、対象ギャップ数をゼロとして報告することを検証する
        /// </summary>
        [Fact]
        public void V_スキャフォールドにギャップが無いときは対象ギャップ数をゼロとして報告する()
        {
            const int l_k長 = 21;
            var l_正解 = V_生成_ランダム配列(300, p_シード: 40);
            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_nogap.fasta", l_正解);
            var l_リードパス = this.V_書き出し_リード("reads_nogap.fq", [l_正解], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(0, l_統計.A_対象ギャップ数);
            Assert.Equal(l_正解, Get_単一配列(l_スキャフォールドパス));
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
        /// scaffold を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_スキャフォールド(string p_ファイル名, string p_配列)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new FastaWriter(l_パス);
            l_書き込み.V_書き込み("SCAFFOLD1", p_配列);
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
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_連番 = 0;
            foreach (var l_ゲノム in p_ゲノム群)
            {
                for (var i = 0; i + p_リード長 <= l_ゲノム.Length; i++)
                {
                    l_書き込み.WriteLine($"@r{l_連番}");
                    l_書き込み.WriteLine(l_ゲノム.Substring(i, p_リード長));
                    l_書き込み.WriteLine("+");
                    l_書き込み.WriteLine(new string('I', p_リード長));
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
            using var l_読み込み = new FastaReader(p_パス);
            Assert.True(l_読み込み.Has続き());
            return l_読み込み.Get_次の配列().A_配列;
        }

        #endregion

    }
}
