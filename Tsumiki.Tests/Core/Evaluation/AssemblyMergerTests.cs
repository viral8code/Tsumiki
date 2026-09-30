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
        /// 骨格の末端が他の骨格にも現れる反復で終わり、一意なアンカーが末端より内側に下がっても、その区間を二重に挟まずに繋ぐこと
        /// </summary>
        [Fact]
        public void V_末端の反復でアンカーが内側に下がっても区間を重複させない()
        {
            var l_左の固有 = V_生成_乱数配列(4_000, p_シード: 691);
            var l_反復 = V_生成_乱数配列(300, p_シード: 692);
            var l_中間 = V_生成_乱数配列(250, p_シード: 693);
            var l_右 = V_生成_乱数配列(4_000, p_シード: 694);
            var l_別の場所 = V_生成_乱数配列(3_000, p_シード: 695);
            var l_真の配列 = l_左の固有 + l_反復 + l_中間 + l_右;

            // 反復の k-mer は 2 本の骨格に現れるので索引から外れる
            var l_骨格 = this.V_書き込み_アセンブリ("backbone_rep.fasta", 63, l_左の固有 + l_反復, l_右, l_反復 + l_別の場所);
            var l_他 = this.V_書き込み_アセンブリ("other_rep.fasta", 31, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_rep.fasta");
            Assert.True(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1));

            var l_結果 = V_読み込み_配列群(l_出力);
            Assert.Contains(l_結果, x => x == l_真の配列 || x == Util.V_逆相補(l_真の配列));
        }

        /// <summary>
        /// 2 本の骨格が末端と先頭で重なっているときは、重なりを 1 回だけにして繋ぐこと
        /// </summary>
        [Fact]
        public void V_重なって隣接する骨格は重なりを1回だけにして繋ぐ()
        {
            var l_左の固有 = V_生成_乱数配列(4_000, p_シード: 701);
            var l_重なり = V_生成_乱数配列(120, p_シード: 702);
            var l_右の固有 = V_生成_乱数配列(4_000, p_シード: 703);
            var l_真の配列 = l_左の固有 + l_重なり + l_右の固有;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_ovl.fasta", 63, l_左の固有 + l_重なり, l_重なり + l_右の固有);
            var l_他 = this.V_書き込み_アセンブリ("other_ovl.fasta", 31, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_ovl.fasta");
            Assert.True(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1));

            var l_結果 = V_読み込み_配列群(l_出力);
            _ = Assert.Single(l_結果);
            Assert.True(l_結果[0] == l_真の配列 || l_結果[0] == Util.V_逆相補(l_真の配列));
        }

        /// <summary>
        /// 長い骨格配列の末端と同じ配列を持つ短い骨格配列があっても、その長い配列を跨いで重複させないこと
        /// </summary>
        [Fact]
        public void V_短い骨格配列が末端を隠しても長い配列を跨いで重複させない()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 711);
            var l_隙間1 = V_生成_乱数配列(200, p_シード: 712);
            var l_中央 = V_生成_乱数配列(5_000, p_シード: 713);
            var l_隙間2 = V_生成_乱数配列(200, p_シード: 714);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 715);
            var l_真の配列 = l_左 + l_隙間1 + l_中央 + l_隙間2 + l_右;

            // 中央の両端と同じ配列を持つ短い断片が、中央の末端の k-mer を重複扱いにする
            var l_骨格 = this.V_書き込み_アセンブリ("backbone_hidden.fasta", 63, l_左, l_中央, l_右, l_中央[..150], l_中央[^150..]);
            var l_他 = this.V_書き込み_アセンブリ("other_hidden.fasta", 31, l_真の配列);

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_hidden.fasta");
            Assert.True(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他], アンカーk長, l_出力, p_必要な独立支持数: 1));

            var l_結果 = V_読み込み_配列群(l_出力);
            Assert.Contains(l_結果, x => x == l_真の配列 || x == Util.V_逆相補(l_真の配列));
            Assert.Equal(l_真の配列.Length + 300, l_結果.Sum(x => x.Length));
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

        /// <summary>
        /// 2 コピーの反復のうち 1 コピーだけが骨格に置かれていれば、もう 1 か所の繋ぎ目に反復を置いて繋ぐこと
        /// </summary>
        /// <remarks>
        /// コピー数が分からなければ、骨格に既にある配列を含む繋ぎ目は重複の持ち込みと区別できないので使わない
        /// </remarks>
        [Fact]
        public void V_コピー数に余裕のある反復を挟む繋ぎ目は繋ぐ()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 721);
            var l_反復 = V_生成_乱数配列(250, p_シード: 722);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 723);
            var l_別の場所 = V_生成_乱数配列(3_000, p_シード: 724) + l_反復 + V_生成_乱数配列(3_000, p_シード: 725);
            var l_真の配列 = l_左 + l_反復 + l_右;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_cn.fasta", 63, l_左, l_右, l_別の場所);
            var l_他1 = this.V_書き込み_アセンブリ("other_cn_a.fasta", 31, l_真の配列);
            var l_他2 = this.V_書き込み_アセンブリ("other_cn_b.fasta", 41, l_真の配列);
            var l_期待コピー数 = Get_期待コピー数(l_反復, 2);

            Assert.False(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他1, l_他2], アンカーk長, Path.Combine(this._作業ディレクトリ, "merged_cn_none.fasta")));

            var l_出力 = Path.Combine(this._作業ディレクトリ, "merged_cn.fasta");
            Assert.True(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他1, l_他2], アンカーk長, l_出力, l_期待コピー数));
            Assert.Contains(V_読み込み_配列群(l_出力), x => x == l_真の配列 || x == Util.V_逆相補(l_真の配列));
        }

        /// <summary>
        /// 反復のコピーが骨格に既に期待コピー数だけ置かれていれば、それ以上は置かないこと
        /// </summary>
        [Fact]
        public void V_コピー数を使い切った反復を挟む繋ぎ目は繋がない()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 731);
            var l_反復 = V_生成_乱数配列(250, p_シード: 732);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 733);
            var l_別の場所1 = V_生成_乱数配列(3_000, p_シード: 734) + l_反復 + V_生成_乱数配列(3_000, p_シード: 735);
            var l_別の場所2 = V_生成_乱数配列(3_000, p_シード: 736) + l_反復 + V_生成_乱数配列(3_000, p_シード: 737);
            var l_真の配列 = l_左 + l_反復 + l_右;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_cnfull.fasta", 63, l_左, l_右, l_別の場所1, l_別の場所2);
            var l_他1 = this.V_書き込み_アセンブリ("other_cnfull_a.fasta", 31, l_真の配列);
            var l_他2 = this.V_書き込み_アセンブリ("other_cnfull_b.fasta", 41, l_真の配列);

            Assert.False(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他1, l_他2], アンカーk長, Path.Combine(this._作業ディレクトリ, "merged_cnfull.fasta"), Get_期待コピー数(l_反復, 2)));
        }

        /// <summary>
        /// 反復を挟む繋ぎ目で、k によって繋ぎ長が食い違うときは繋がないこと
        /// </summary>
        /// <remarks>
        /// 別々の k が反復の別々のコピーを通って同じ 2 本を結んでいる形で、どちらの長さも正しいとは限らない
        /// </remarks>
        [Fact]
        public void V_kの間で繋ぎ長が揃わない反復の繋ぎ目は繋がない()
        {
            var l_左 = V_生成_乱数配列(5_000, p_シード: 741);
            var l_反復 = V_生成_乱数配列(250, p_シード: 742);
            var l_挿入 = V_生成_乱数配列(60, p_シード: 743);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 744);
            var l_別の場所 = V_生成_乱数配列(3_000, p_シード: 745) + l_反復 + V_生成_乱数配列(3_000, p_シード: 746);

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_len.fasta", 63, l_左, l_右, l_別の場所);
            var l_他1 = this.V_書き込み_アセンブリ("other_lenmix_a.fasta", 31, l_左 + l_反復 + l_右);
            var l_他2 = this.V_書き込み_アセンブリ("other_lenmix_b.fasta", 41, l_左 + l_反復 + l_挿入 + l_右);

            Assert.False(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他1, l_他2], アンカーk長, Path.Combine(this._作業ディレクトリ, "merged_lenmix.fasta"), Get_期待コピー数(l_反復, 2)));
        }

        /// <summary>
        /// 骨格配列の末端がゲノム中では複数コピーあるのに骨格に 1 つしか無い (畳まれた反復) なら、そこからは繋がないこと
        /// </summary>
        /// <remarks>
        /// 畳まれた反復の先はコピーごとに違い、別の k が通ったコピーが正しい行き先とは限らない
        /// </remarks>
        [Fact]
        public void V_畳まれた反復で終わる骨格配列からは反復を挟んで繋がない()
        {
            var l_左の固有 = V_生成_乱数配列(4_000, p_シード: 751);
            var l_畳まれた反復 = V_生成_乱数配列(800, p_シード: 752);
            var l_反復 = V_生成_乱数配列(250, p_シード: 753);
            var l_右 = V_生成_乱数配列(5_000, p_シード: 754);
            var l_別の場所 = V_生成_乱数配列(3_000, p_シード: 755) + l_反復 + V_生成_乱数配列(3_000, p_シード: 756);
            var l_真の配列 = l_左の固有 + l_畳まれた反復 + l_反復 + l_右;

            var l_骨格 = this.V_書き込み_アセンブリ("backbone_fold.fasta", 63, l_左の固有 + l_畳まれた反復, l_右, l_別の場所);
            var l_他1 = this.V_書き込み_アセンブリ("other_fold_a.fasta", 31, l_真の配列);
            var l_他2 = this.V_書き込み_アセンブリ("other_fold_b.fasta", 41, l_真の配列);
            var l_反復のkmer = Get_kmer集合(l_反復);
            var l_畳まれたkmer = Get_kmer集合(l_畳まれた反復);

            int l_期待コピー数(string p_配列, int p_位置)
            {
                var l_kmer = p_配列.Substring(p_位置, アンカーk長);
                return l_反復のkmer.Contains(l_kmer) || l_畳まれたkmer.Contains(l_kmer) ? 2 : 1;
            }

            Assert.False(AssemblyMerger.Try統合(l_骨格, [l_骨格, l_他1, l_他2], アンカーk長, Path.Combine(this._作業ディレクトリ, "merged_fold.fasta"), l_期待コピー数));
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
        /// 配列の k-mer を両方の向きで集める
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <returns>k-mer の集合</returns>
        private static HashSet<string> Get_kmer集合(string p_配列)
        {
            HashSet<string> l_集合 = [];
            foreach (var l_向き in new[] { p_配列, Util.V_逆相補(p_配列) })
            {
                for (var i = 0; i + アンカーk長 <= l_向き.Length; i++)
                {
                    _ = l_集合.Add(l_向き.Substring(i, アンカーk長));
                }
            }

            return l_集合;
        }

        /// <summary>
        /// 反復の k-mer だけを指定のコピー数、それ以外を 1 コピーと答える期待コピー数
        /// </summary>
        /// <param name="p_反復">複数コピーの配列</param>
        /// <param name="p_コピー数">反復のコピー数</param>
        /// <returns>期待コピー数を返す関数</returns>
        private static Func<string, int, int> Get_期待コピー数(string p_反復, int p_コピー数)
        {
            var l_反復のkmer = Get_kmer集合(p_反復);
            return (p_配列, p_位置) => l_反復のkmer.Contains(p_配列.Substring(p_位置, アンカーk長)) ? p_コピー数 : 1;
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
