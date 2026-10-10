using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// リファレンス無しでアセンブリの良さを測る評価器の検証
    /// </summary>
    public class AssemblyScorerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki scorer tests
        /// </summary>
        private const string C_項目_tsumiki_scorer_tests = "tsumiki_scorer_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 perfect fasta
        /// </summary>
        private const string C_ファイル名_perfect_fasta = "perfect.fasta";

        /// <summary>
        /// ファイル名 whole fasta
        /// </summary>
        private const string C_ファイル名_whole_fasta = "whole.fasta";

        /// <summary>
        /// ファイル名 fragmented fasta
        /// </summary>
        private const string C_ファイル名_fragmented_fasta = "fragmented.fasta";

        /// <summary>
        /// ファイル名 honest fasta
        /// </summary>
        private const string C_ファイル名_honest_fasta = "honest.fasta";

        /// <summary>
        /// ファイル名 chimera fasta
        /// </summary>
        private const string C_ファイル名_chimera_fasta = "chimera.fasta";

        /// <summary>
        /// 項目 this test is only meanin contiguity alone
        /// </summary>
        private const string C_項目_this_test_is_only_meanin_contiguity_alone = "this test is only meaningful while the chimera looks better on contiguity alone";

        /// <summary>
        /// ファイル名 single fasta
        /// </summary>
        private const string C_ファイル名_single_fasta = "single.fasta";

        /// <summary>
        /// ファイル名 doubled fasta
        /// </summary>
        private const string C_ファイル名_doubled_fasta = "doubled.fasta";

        /// <summary>
        /// ファイル名 repeat fasta
        /// </summary>
        private const string C_ファイル名_repeat_fasta = "repeat.fasta";

        /// <summary>
        /// ファイル名 partial fasta
        /// </summary>
        private const string C_ファイル名_partial_fasta = "partial.fasta";

        /// <summary>
        /// ファイル名 circular fasta
        /// </summary>
        private const string C_ファイル名_circular_fasta = "circular.fasta";

        /// <summary>
        /// 項目 NODE1 circular
        /// </summary>
        private const string C_項目_NODE1_circular = "NODE1_circular";

        /// <summary>
        /// ファイル名 linear fasta
        /// </summary>
        private const string C_ファイル名_linear_fasta = "linear.fasta";

        /// <summary>
        /// 項目 NODE1
        /// </summary>
        private const string C_項目_NODE1 = "NODE1";

        /// <summary>
        /// ファイル名 with plasmid fasta
        /// </summary>
        private const string C_ファイル名_with_plasmid_fasta = "with_plasmid.fasta";

        /// <summary>
        /// 項目 NODE2 circular
        /// </summary>
        private const string C_項目_NODE2_circular = "NODE2_circular";

        /// <summary>
        /// ファイル名 with artefacts fasta
        /// </summary>
        private const string C_ファイル名_with_artefacts_fasta = "with_artefacts.fasta";

        /// <summary>
        /// 項目 NODE3 circular
        /// </summary>
        private const string C_項目_NODE3_circular = "NODE3_circular";

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 21;

        /// <summary>
        /// 深さ
        /// </summary>
        private const int C_深さ = 20;

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
        public AssemblyScorerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_scorer_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// 完全なアセンブリは完全性も正確性も満点になることを確かめる
        /// </summary>
        [Fact]
        public void V_完全なアセンブリは完全性も正確性も満点になる()
        {
            var l_正解 = V_乱数配列(20_000, p_シード: 501);
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_パス = this.V_書き出し_Fasta(C_ファイル名_perfect_fasta, l_正解);
            var l_評価 = AssemblyScorer.Get_評価(l_パス, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価);
            Assert.Equal(0L, l_評価.A_欠損延べ数);
            Assert.Equal(0L, l_評価.A_過剰延べ数);
            Assert.Equal(1.0D, l_評価.A_完全性, 6);
            Assert.Equal(1.0D, l_評価.A_正確性, 6);
        }

        /// <summary>
        /// 断片化しているが取りこぼしの無いアセンブリは連続性だけが落ちることを確かめる
        /// </summary>
        [Fact]
        public void V_断片化しているが取りこぼしの無いアセンブリは連続性だけが落ちる()
        {
            var l_正解 = V_乱数配列(20_000, p_シード: 502);
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_断片 = new List<string>();
            for (var i = 0; i < l_正解.Length; i += 4_000)
            {
                var l_終端 = Math.Min(l_正解.Length, i + 4_000 + C_k長 - 1);
                l_断片.Add(l_正解[i..l_終端]);
            }

            var l_一本 = this.V_書き出し_Fasta(C_ファイル名_whole_fasta, l_正解);
            var l_断片化 = this.V_書き出し_Fasta(C_ファイル名_fragmented_fasta, [..l_断片]);
            var l_評価一本 = AssemblyScorer.Get_評価(l_一本, l_インデックス, C_k長, C_深さ, l_正解.Length);
            var l_評価断片 = AssemblyScorer.Get_評価(l_断片化, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価一本);
            Assert.NotNull(l_評価断片);
            Assert.Equal(0L, l_評価断片.A_欠損延べ数);
            Assert.True(l_評価断片.A_NG50 < l_評価一本.A_NG50);
        }

        /// <summary>
        /// 配列を飛ばしたキメラは断片化しているが正直なアセンブリより低く評価されることを確かめる
        /// </summary>
        [Fact]
        public void V_配列を飛ばしたキメラは断片化しているが正直なアセンブリより低く評価される()
        {
            var l_先頭配列 = V_乱数配列(8_000, p_シード: 511);
            var l_反復配列 = V_乱数配列(300, p_シード: 512);
            var l_中間配列 = V_乱数配列(8_000, p_シード: 513);
            var l_末尾配列 = V_乱数配列(8_000, p_シード: 514);
            var l_正解 = l_先頭配列 + l_反復配列 + l_中間配列 + l_反復配列 + l_末尾配列;
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_正直 = this.V_書き出し_Fasta(C_ファイル名_honest_fasta, l_先頭配列 + l_反復配列, l_反復配列 + l_中間配列 + l_反復配列, l_反復配列 + l_末尾配列);
            var l_キメラ = this.V_書き出し_Fasta(C_ファイル名_chimera_fasta, l_先頭配列 + l_反復配列 + l_末尾配列);
            var l_評価正直 = AssemblyScorer.Get_評価(l_正直, l_インデックス, C_k長, C_深さ, l_正解.Length);
            var l_評価キメラ = AssemblyScorer.Get_評価(l_キメラ, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価正直);
            Assert.NotNull(l_評価キメラ);
            Assert.True(l_評価キメラ.A_欠損延べ数 > l_評価正直.A_欠損延べ数, $"chimera missing={l_評価キメラ.A_欠損延べ数}, honest missing={l_評価正直.A_欠損延べ数}");
            Assert.True(l_評価キメラ.A_完全性 < l_評価正直.A_完全性);
            Assert.True(l_評価正直.A_完全性 - l_評価キメラ.A_完全性 > AssemblySelector.C_完全性の許容差, $"chimera completeness={l_評価キメラ.A_完全性:F3}, honest={l_評価正直.A_完全性:F3}");
            Assert.True(l_評価キメラ.A_NG50 > l_評価正直.A_NG50, C_項目_this_test_is_only_meanin_contiguity_alone);
        }

        /// <summary>
        /// 配列を 2 回出した水増しは正確性が落ちて総合点が下がることを確かめる
        /// </summary>
        [Fact]
        public void V_配列を2回出した水増しは正確性が落ちて総合点が下がる()
        {
            var l_正解 = V_乱数配列(20_000, p_シード: 521);
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_正常 = this.V_書き出し_Fasta(C_ファイル名_single_fasta, l_正解);
            var l_水増し = this.V_書き出し_Fasta(C_ファイル名_doubled_fasta, l_正解, l_正解);
            var l_評価正常 = AssemblyScorer.Get_評価(l_正常, l_インデックス, C_k長, C_深さ, l_正解.Length);
            var l_評価水増し = AssemblyScorer.Get_評価(l_水増し, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価正常);
            Assert.NotNull(l_評価水増し);
            Assert.True(l_評価水増し.A_過剰延べ数 > 0L);
            Assert.True(l_評価水増し.A_正確性 < l_評価正常.A_正確性);
        }

        /// <summary>
        /// 2 コピーの反復配列を 2 回出しても減点されないことを確かめる
        /// </summary>
        [Fact]
        public void V_二コピーの反復配列を2回出しても減点されない()
        {
            var l_単一 = V_乱数配列(10_000, p_シード: 531);
            var l_反復 = V_乱数配列(500, p_シード: 532);
            var l_正解 = l_単一 + l_反復 + V_乱数配列(5_000, p_シード: 533) + l_反復;
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_パス = this.V_書き出し_Fasta(C_ファイル名_repeat_fasta, l_正解);
            var l_評価 = AssemblyScorer.Get_評価(l_パス, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価);
            Assert.Equal(0L, l_評価.A_過剰延べ数);
            Assert.Equal(0L, l_評価.A_欠損延べ数);
        }

        /// <summary>
        /// NG50 は自分の総延長ではなく推定ゲノムサイズを分母にすることを確かめる
        /// </summary>
        [Fact]
        public void V_NG50は自分の総延長ではなく推定ゲノムサイズを分母にする()
        {
            var l_正解 = V_乱数配列(20_000, p_シード: 541);
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_一部 = this.V_書き出し_Fasta(C_ファイル名_partial_fasta, l_正解[..8_000]);
            var l_評価 = AssemblyScorer.Get_評価(l_一部, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価);
            Assert.Equal(0L, l_評価.A_NG50);
            Assert.InRange(l_評価.A_完全性, 0.35D, 0.45D);
        }

        /// <summary>
        /// 環状に閉じた complicon (ContigMaker が名前に "circular" を付けたもの) は本数・総延長に数えられることを確かめる
        /// </summary>
        [Fact]
        public void V_環状contigは環状本数に数えられる()
        {
            var l_正解 = V_乱数配列(20_000, p_シード: 601);
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_パス = this.V_書き出し_Fasta_名前付き(C_ファイル名_circular_fasta, (C_項目_NODE1_circular, l_正解));
            var l_評価 = AssemblyScorer.Get_評価(l_パス, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価);
            Assert.Equal(1, l_評価.A_環状本数);
            Assert.Equal(1.0D, l_評価.A_環状化率, 6);
        }

        /// <summary>
        /// 線状の contig は環状として数えられないことを確かめる
        /// </summary>
        [Fact]
        public void V_線状contigは環状として数えられない()
        {
            var l_正解 = V_乱数配列(20_000, p_シード: 602);
            using var l_インデックス = this.V_構築_索引(l_正解);
            var l_パス = this.V_書き出し_Fasta_名前付き(C_ファイル名_linear_fasta, (C_項目_NODE1, l_正解));
            var l_評価 = AssemblyScorer.Get_評価(l_パス, l_インデックス, C_k長, C_深さ, l_正解.Length);
            Assert.NotNull(l_評価);
            Assert.Equal(0, l_評価.A_環状本数);
            Assert.Equal(0.0D, l_評価.A_環状化率);
        }

        /// <summary>
        /// 染色体よりはるかに小さいプラスミドでも、複製単位として数えられる長さがあれば環状化率に反映されることを確かめる
        /// </summary>
        [Fact]
        public void V_小さい環状プラスミドも環状化率に反映される()
        {
            var l_染色体 = V_乱数配列(20_000, p_シード: 603);
            var l_プラスミド = V_乱数配列(2_000, p_シード: 604);
            using var l_インデックス = this.V_構築_索引(l_染色体, l_プラスミド);
            var l_genomeSize = l_染色体.Length + l_プラスミド.Length;
            var l_パス = this.V_書き出し_Fasta_名前付き(C_ファイル名_with_plasmid_fasta, (C_項目_NODE1_circular, l_染色体), (C_項目_NODE2_circular, l_プラスミド));
            var l_評価 = AssemblyScorer.Get_評価(l_パス, l_インデックス, C_k長, C_深さ, l_genomeSize);
            Assert.NotNull(l_評価);
            Assert.Equal(2, l_評価.A_環状本数);
            Assert.InRange(l_評価.A_環状化率, 0.99D, 1.0D);
        }

        /// <summary>
        /// 評価に含める最小長 (500 bp) を下回る配列は、環状の目印が付いていても数えられないことを確かめる
        /// </summary>
        [Fact]
        public void V_評価に含める最小長を下回る配列は数えられない()
        {
            var l_染色体 = V_乱数配列(20_000, p_シード: 605);
            using var l_インデックス = this.V_構築_索引(l_染色体);
            var l_パス = this.V_書き出し_Fasta_名前付き(C_ファイル名_with_artefacts_fasta, (C_項目_NODE1_circular, l_染色体), (C_項目_NODE2_circular, l_染色体[..30]), (C_項目_NODE3_circular, l_染色体[100..106]));
            var l_評価 = AssemblyScorer.Get_評価(l_パス, l_インデックス, C_k長, C_深さ, l_染色体.Length);
            Assert.NotNull(l_評価);
            Assert.Equal(1, l_評価.A_環状本数);
            Assert.Equal(1, l_評価.A_本数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_乱数配列(int p_長さ, int p_シード)
        {
            var l_乱数生成器 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数生成器.Next(4)]));
        }

        /// <summary>
        /// 与えた配列群から k-mer インデックスを作る
        /// </summary>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_索引(params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 1
            };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            foreach (var l_配列 in p_配列群)
            {
                var l_バイト列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + C_k長 <= l_バイト列.Length; i++)
                {
                    for (var l_反復回数 = 0; l_反復回数 < C_深さ; l_反復回数++)
                    {
                        l_インデックス.V_登録(l_バイト列.AsSpan(i, C_k長));
                    }
                }
            }

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列群">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_Fasta(string p_ファイル名, params string[] p_配列群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new FastaWriter(l_パス);
            var l_ID = 1;
            foreach (var l_配列 in p_配列群)
            {
                l_書き込み.V_書き込み($"NODE{l_ID++}", l_配列);
            }

            return l_パス;
        }

        /// <summary>
        /// 名前を付けた配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_エントリ群">書き出す名前と配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_Fasta_名前付き(string p_ファイル名, params (string A_名前, string A_配列)[] p_エントリ群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new FastaWriter(l_パス);
            foreach (var (l_名前, l_配列)in p_エントリ群)
            {
                l_書き込み.V_書き込み(l_名前, l_配列);
            }

            return l_パス;
        }

        #endregion
    }
}
