using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// スキャフォールドのギャップ (N の連続) を、de Bruijn グラフ上で
    /// 両端を繋ぐ経路を探して実配列に置き換える処理の検証
    /// </summary>
    /// <remarks>
    /// contig が途切れるのは配列が存在しないからではなく、分岐でどちらへ
    /// 進むか決められなかったからであることが多い<br/>
    /// その場合ギャップを埋める
    /// 配列は k-mer 集合の中に実在しており、両端から辿れば復元できる
    /// </remarks>
    public class GapFillerTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        /// <summary>
        /// 一時ディレクトリを作る
        /// </summary>
        public GapFillerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_gapfiller_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
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
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_kmer長">k 長</param>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_索引(int p_kmer長, params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmer長, A_スレッド数 = 1 };
            var l_索引 = new TrustedKmerIndex(this._tempDir);
            foreach (var l_配列 in p_配列群)
            {
                var l_バイト列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + p_kmer長 <= l_バイト列.Length; i++)
                {
                    for (var l_繰り返し = 0; l_繰り返し < 3; l_繰り返し++)
                    {
                        l_索引.V_登録(l_バイト列.AsSpan(i, p_kmer長));
                    }
                }
            }
            _ = l_索引.V_カットオフ(p_カットオフ: 2);
            return l_索引;
        }

        /// <summary>
        /// スキャフォールドを FASTA として書き出す
        /// </summary>
        /// <param name="p_名前">ファイル名</param>
        /// <param name="p_配列">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き込み_スキャフォールド(string p_名前, string p_配列)
        {
            var l_パス = Path.Combine(this._tempDir, p_名前);
            using (var l_ライター = new FastaWriter(l_パス))
            {
                l_ライター.V_書き込み("SCAFFOLD1", p_配列);
            }
            return l_パス;
        }

        /// <summary>
        /// FASTA から 1 本だけの配列を読み込む
        /// </summary>
        /// <param name="p_パス">読み込むパス</param>
        /// <returns>読み込んだ配列</returns>
        private static string V_読み込み_単一配列(string p_パス)
        {
            using var l_リーダー = new FastaReader(p_パス);
            Assert.True(l_リーダー.Get_続きがあるか());
            return l_リーダー.Get_次の配列().A_配列;
        }

        /// <summary>
        /// 唯一の経路がグラフ上にあれば真の配列に復元される
        /// </summary>
        [Fact]
        public void 唯一経路の場合は真の配列に復元される()
        {
            const int l_k長 = 21;
            // 200 bp の非反復的な配列
            // k=21 なので偶然の重複はまず起きない
            var l_正解配列 = V_生成_ランダム配列(200, p_シード: 20260903);

            using var l_索引 = this.V_構築_索引(l_k長, l_正解配列);

            // 真ん中 40 bp を N に置き換えたスキャフォールドを作る
            const int l_ギャップ開始 = 80;
            const int l_ギャップ長 = 40;
            var l_ギャップ入り配列 = l_正解配列[..l_ギャップ開始] + new string('N', l_ギャップ長) + l_正解配列[(l_ギャップ開始 + l_ギャップ長)..];
            var l_パス = this.V_書き込み_スキャフォールド("scaffolds.fasta", l_ギャップ入り配列);

            var l_統計 = GapFiller.V_充填_ギャップ(l_パス, l_索引, l_k長);

            Assert.Equal(1, l_統計.A_総ギャップ数);
            Assert.Equal(1, l_統計.A_埋めたギャップ数);
            Assert.Equal(l_ギャップ長, l_統計.A_埋めた塩基数);

            // 埋めた結果は元の配列そのものに戻っていなければならない
            Assert.Equal(l_正解配列, V_読み込み_単一配列(l_パス));
        }

        /// <summary>
        /// ギャップ長推定が多少ずれていてもマージンで埋まる
        /// </summary>
        [Fact]
        public void ギャップ長推定が多少ずれていてもマージンで埋まる()
        {
            const int l_k長 = 21;
            var l_正解配列 = V_生成_ランダム配列(200, p_シード: 7);

            using var l_索引 = this.V_構築_索引(l_k長, l_正解配列);

            // 実際の欠損は 40 bp だが、推定を誤って 30 個の N になっている状況
            // ギャップ長推定はインサートサイズ推定のばらつきを引き継ぐため、
            // ぴったりの長さしか探さないと現実にはまず埋まらない
            const int l_ギャップ開始 = 80;
            const int l_実際の欠損 = 40;
            var l_ギャップ入り配列 = l_正解配列[..l_ギャップ開始] + new string('N', 30) + l_正解配列[(l_ギャップ開始 + l_実際の欠損)..];
            var l_パス = this.V_書き込み_スキャフォールド("scaffolds_off.fasta", l_ギャップ入り配列);

            var l_統計 = GapFiller.V_充填_ギャップ(l_パス, l_索引, l_k長);

            Assert.Equal(1, l_統計.A_埋めたギャップ数);
            Assert.Equal(l_正解配列, V_読み込み_単一配列(l_パス));
        }

        /// <summary>
        /// ギャップを埋める経路が複数ある場合、どれが正しいか決められない
        /// </summary>
        /// <remarks>
        /// 誤った配列で埋めるより N のまま残すほうが下流の解析にとって安全
        /// </remarks>
        [Fact]
        public void 経路が複数ある場合はNのまま残す()
        {
            const int l_k長 = 21;
            var l_前半 = V_生成_ランダム配列(80, p_シード: 11);
            var l_後半 = V_生成_ランダム配列(80, p_シード: 12);
            // 同じ長さで中身だけ違う 2 通りの中間配列を、どちらも k-mer 集合に入れる
            var l_中間A = V_生成_ランダム配列(40, p_シード: 13);
            var l_中間B = V_生成_ランダム配列(40, p_シード: 14);

            using var l_索引 = this.V_構築_索引(l_k長, l_前半 + l_中間A + l_後半, l_前半 + l_中間B + l_後半);

            var l_ギャップ入り配列 = l_前半 + new string('N', 40) + l_後半;
            var l_パス = this.V_書き込み_スキャフォールド("scaffolds_ambiguous.fasta", l_ギャップ入り配列);

            var l_統計 = GapFiller.V_充填_ギャップ(l_パス, l_索引, l_k長);

            Assert.Equal(1, l_統計.A_総ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_一意に定まらなかった数);
            // N はそのまま残っていること
            Assert.Contains('N', V_読み込み_単一配列(l_パス));
        }

        /// <summary>
        /// 両端を繋ぐ経路がグラフ上に存在しない (本当に配列が無い) 場合は、
        /// 当然埋められない
        /// </summary>
        [Fact]
        public void 経路が存在しない場合はNのまま残す()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(80, p_シード: 21);
            var l_右 = V_生成_ランダム配列(80, p_シード: 22);

            // 左右それぞれの k-mer は入れるが、両者を繋ぐ配列は入れない
            using var l_索引 = this.V_構築_索引(l_k長, l_左, l_右);

            var l_ギャップ入り配列 = l_左 + new string('N', 40) + l_右;
            var l_パス = this.V_書き込み_スキャフォールド("scaffolds_unreachable.fasta", l_ギャップ入り配列);

            var l_統計 = GapFiller.V_充填_ギャップ(l_パス, l_索引, l_k長);

            Assert.Equal(1, l_統計.A_総ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_到達できなかった数);
            Assert.Contains('N', V_読み込み_単一配列(l_パス));
        }

        /// <summary>
        /// ギャップが無ければ配列はそのまま保たれる
        /// </summary>
        [Fact]
        public void ギャップが無ければ配列はそのまま()
        {
            const int l_k長 = 21;
            var l_正解配列 = V_生成_ランダム配列(150, p_シード: 31);
            using var l_索引 = this.V_構築_索引(l_k長, l_正解配列);

            var l_パス = this.V_書き込み_スキャフォールド("scaffolds_nogap.fasta", l_正解配列);
            var l_統計 = GapFiller.V_充填_ギャップ(l_パス, l_索引, l_k長);

            Assert.Equal(0, l_統計.A_総ギャップ数);
            Assert.Equal(l_正解配列, V_読み込み_単一配列(l_パス));
        }
    }
}
