using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 複数の k のアセンブリを統合する処理の検証
    /// </summary>
    /// <remarks>
    /// 統合は誤った連結を持ち込みうる操作なので、繋ぐべきときに繋ぐことと同じくらい、根拠が無いときに繋がないことを固定しておく必要がある
    /// </remarks>
    public class AssemblyMergerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// アンカー k 長
        /// </summary>
        private const int アンカーk長 = 31;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 一時ディレクトリを作る
        /// </summary>
        public AssemblyMergerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_merger_tests_" + Guid.NewGuid().ToString("N"));
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
        /// 骨格が途切れている箇所を、別の k の配列が跨いでいる場合
        /// </summary>
        /// <remarks>
        /// 繋いだ結果が元のゲノムそのものに戻ること
        /// </remarks>
        [Fact]
        public void V_別のkが骨格の切れ目を跨ぐと繋いで元の配列に戻る()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 601);
            var l_中間 = V_生成_乱数配列(400, p_シード: 602);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 603);
            var l_真の配列 = l_左 + l_中間 + l_右;

            // 骨格は中間で切れている
            var l_骨格 = this.V_書き込み_アセンブリ("backbone.fasta", 63, l_左, l_右);
            // 別の k は切れ目を跨いでいる (両端に十分なアンカーを持つ)
            var l_他 = this.V_書き込み_アセンブリ("other.fasta", 31, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged.fasta");
            var l_繋いだか = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            Assert.True(l_繋いだか);
            var l_結果 = V_読み込み_配列群(l_出力);
            _ = Assert.Single(l_結果);
            Assert.Equal(l_真の配列, l_結果[0]);
        }

        /// <summary>
        /// 骨格側の片方が逆向きに出力されていても、向きを揃えて繋げること
        /// </summary>
        [Fact]
        public void V_骨格の片方が逆相補でも正しく繋がる()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 611);
            var l_中間 = V_生成_乱数配列(300, p_シード: 612);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 613);
            var l_真の配列 = l_左 + l_中間 + l_右;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_rc.fasta", 63, l_左, Util.V_逆相補(l_右));
            var l_他 = this.V_書き込み_アセンブリ("other_rc.fasta", 31, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_rc.fasta");
            var l_繋いだか = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            Assert.True(l_繋いだか);
            var l_結果 = V_読み込み_配列群(l_出力);
            _ = Assert.Single(l_結果);
            Assert.True(l_結果[0] == l_真の配列 || l_結果[0] == Util.V_逆相補(l_真の配列), "merged sequence should be the truth in one orientation or the other");
        }

        /// <summary>
        /// 跨いでいる配列が無ければ何もしないこと
        /// </summary>
        /// <remarks>
        /// 根拠が無いのに繋ぐのが最も避けたい失敗
        /// </remarks>
        [Fact]
        public void V_跨ぐ配列が無ければ何もしない()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 621);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 622);

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_none.fasta", 63, l_左, l_右);
            // 別の k も同じところで切れている
            var l_他 = this.V_書き込み_アセンブリ("other_none.fasta", 31, l_左, l_右);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_none.fasta");
            var l_繋いだか = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            Assert.False(l_繋いだか);
            Assert.False(File.Exists(l_出力));
        }

        /// <summary>
        /// 反復配列のせいで行き先が 2 つある場合は繋がないこと
        /// </summary>
        /// <remarks>
        /// 片方を選ぶ根拠が無く、選べば誤アセンブリになる
        /// </remarks>
        [Fact]
        public void V_行き先が2つ有り得るときは繋がない()
        {
            var l_共通の左 = V_生成_乱数配列(5_000, p_シード: 631);
            var l_右候補1 = V_生成_乱数配列(5_000, p_シード: 632);
            var l_右候補2 = V_生成_乱数配列(5_000, p_シード: 633);
            var l_中間 = V_生成_乱数配列(200, p_シード: 634);

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_amb.fasta", 63, l_共通の左, l_右候補1, l_右候補2);
            // 同じ左から 2 つの異なる右へ繋がる証拠が両方ある
            var l_他 = this.V_書き込み_アセンブリ("other_amb.fasta", 31, l_共通の左 + l_中間 + l_右候補1, l_共通の左 + l_中間 + l_右候補2);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_amb.fasta");
            var l_繋いだか = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            Assert.False(l_繋いだか);
        }

        /// <summary>
        /// 3 本を 2 箇所で繋ぐ連鎖
        /// </summary>
        /// <remarks>
        /// 1 回の統合で最後まで繋がること
        /// </remarks>
        [Fact]
        public void V_連鎖する3本を1回の統合で全て繋ぐ()
        {
            var l_先頭配列 = V_生成_乱数配列(4_000, p_シード: 641);
            var l_g1 = V_生成_乱数配列(200, p_シード: 642);
            var l_中間配列 = V_生成_乱数配列(4_000, p_シード: 643);
            var l_g2 = V_生成_乱数配列(200, p_シード: 644);
            var l_末尾配列 = V_生成_乱数配列(4_000, p_シード: 645);
            var l_真の配列 = l_先頭配列 + l_g1 + l_中間配列 + l_g2 + l_末尾配列;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_chain.fasta", 63, l_先頭配列, l_中間配列, l_末尾配列);
            var l_他 = this.V_書き込み_アセンブリ("other_chain.fasta", 31, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_chain.fasta");
            var l_繋いだか = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            Assert.True(l_繋いだか);
            var l_結果 = V_読み込み_配列群(l_出力);
            _ = Assert.Single(l_結果);
            Assert.Equal(l_真の配列, l_結果[0]);
        }

        /// <summary>
        /// 繋がらなかった骨格配列も、統合結果から失われないこと
        /// </summary>
        [Fact]
        public void V_繋がらなかった骨格配列も出力に残る()
        {
            var l_左 = V_生成_乱数配列(4_000, p_シード: 651);
            var l_中間 = V_生成_乱数配列(200, p_シード: 652);
            var l_右 = V_生成_乱数配列(4_000, p_シード: 653);
            var l_孤立 = V_生成_乱数配列(3_000, p_シード: 654);

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_iso.fasta", 63, l_左, l_右, l_孤立);
            var l_他 = this.V_書き込み_アセンブリ("other_iso.fasta", 31, l_左 + l_中間 + l_右);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_iso.fasta");
            var l_繋いだか = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            Assert.True(l_繋いだか);
            var l_結果 = V_読み込み_配列群(l_出力);
            Assert.Equal(2, l_結果.Count);
            Assert.Contains(l_結果, x => x == l_左 + l_中間 + l_右);
            Assert.Contains(l_結果, x => x == l_孤立 || x == Util.V_逆相補(l_孤立));
        }

        /// <summary>
        /// 既定では、1 つの k だけが主張する隣接は採らないこと
        /// </summary>
        [Fact]
        public void V_支持するkが1つだけの接合は既定では採用しない()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 671);
            var l_中間 = V_生成_乱数配列(300, p_シード: 672);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 673);

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_sup.fasta", 63, l_左, l_右);
            var l_他 = this.V_書き込み_アセンブリ("other_sup.fasta", 31, l_左 + l_中間 + l_右);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_sup.fasta");

            Assert.False(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力));
        }

        /// <summary>
        /// 2 つの k が同じ隣接を主張していれば採ること
        /// </summary>
        [Fact]
        public void V_独立した2つのkが一致すれば採用する()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 681);
            var l_中間 = V_生成_乱数配列(300, p_シード: 682);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 683);
            var l_真の配列 = l_左 + l_中間 + l_右;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_two.fasta", 63, l_左, l_右);
            var l_他1 = this.V_書き込み_アセンブリ("other_two_a.fasta", 31, l_真の配列);
            var l_他2 = this.V_書き込み_アセンブリ("other_two_b.fasta", 41, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_two.fasta");

            Assert.True(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他1, l_他2], アンカーk長, l_出力));
            Assert.Equal(l_真の配列, V_読み込み_配列群(l_出力)[0]);
        }

        /// <summary>
        /// 統合の総延長が、骨格の総延長を下回らないこと
        /// </summary>
        /// <remarks>
        /// 配列を落とすなら統合しないほうがましなので、これは不変条件
        /// </remarks>
        [Fact]
        public void Try統合結果の総延長は骨格を下回らない()
        {
            var l_先頭配列 = V_生成_乱数配列(4_000, p_シード: 661);
            var l_g = V_生成_乱数配列(150, p_シード: 662);
            var l_中間配列 = V_生成_乱数配列(4_000, p_シード: 663);
            var l_孤立 = V_生成_乱数配列(2_000, p_シード: 664);

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_len.fasta", 63, l_先頭配列, l_中間配列, l_孤立);
            var l_他 = this.V_書き込み_アセンブリ("other_len.fasta", 31, l_先頭配列 + l_g + l_中間配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_len.fasta");
            _ = AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1);

            var l_骨格の総延長 = l_先頭配列.Length + l_中間配列.Length + l_孤立.Length;
            Assert.True(V_読み込み_配列群(l_出力).Sum(x => x.Length) >= l_骨格の総延長);
        }

        #endregion

        #region 内部メソッド

        // 既定では 2 つ以上の k による裏付けを求める
        // 以下の多くのテストは
        // 証拠源が 1 つの状況を見たいので、明示的に 1 を渡している
        // 照合そのものは専用のテストで確かめる

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 配列を FASTA として書き出し、アセンブリの実行結果として返す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_配列群">書き出す配列</param>
        /// <returns>アセンブリの実行結果</returns>
        private アセンブリ実行結果 V_書き込み_アセンブリ(string p_ファイル名, int p_k長, params string[] p_配列群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using (var l_ライター = new FastaWriter(l_パス))
            {
                var l_通し番号 = 1;
                foreach (var l_配列 in p_配列群)
                {
                    l_ライター.V_書き込み($"NODE{l_通し番号++}", l_配列);
                }
            }
            return new アセンブリ実行結果(p_k長, l_パス, l_パス, null, 2UL, 20.0D);
        }

        /// <summary>
        /// FASTA から配列をすべて読み込む
        /// </summary>
        /// <param name="p_パス">読み込むパス</param>
        /// <returns>読み込んだ配列</returns>
        private static List<string> V_読み込み_配列群(string p_パス)
        {
            List<string> l_結果 = [];
            using var l_リーダー = new FastaReader(p_パス);
            while (l_リーダー.Has続き())
            {
                l_結果.Add(l_リーダー.Get_次の配列().A_配列);
            }
            return l_結果;
        }

        #endregion

    }
}
