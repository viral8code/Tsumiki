using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 前段の k の配列を次の k へ引き継ぐ処理の検証
    /// </summary>
    /// <remarks>
    /// 引き継ぎで最も壊れやすいのはカバレッジである<br/>
    /// 名目値で埋めるとコピー数推定・低カバレッジ端のトリミング・自己検査がまとめて狂う<br/>
    /// 連結が保たれることと同じくらい、カバレッジが保たれることを固定する
    /// </remarks>
    public class KmerCarryOverTests : IDisposable
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
        public KmerCarryOverTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_carryover_tests_" + Guid.NewGuid().ToString("N"));
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
        /// 配列中のすべての k-mer についてカバレッジを記録することを検証する
        /// </summary>
        [Fact]
        public void V_配列中のすべてのkmerのカバレッジを記録する()
        {
            const int l_k長 = 21;
            var l_正解 = V_生成_ランダム配列(2_000, p_シード: 701);
            using var l_インデックス = this.V_構築_インデックス(l_k長, p_深さ: 17, l_正解);

            var l_パス = this.V_書き出し_FASTA("prep.fasta", l_正解);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(l_パス, l_インデックス, l_k長);

            var l_項目 = Assert.Single(l_引き継ぎ);
            Assert.Equal(l_正解, l_項目.A_配列);
            Assert.Equal(l_k長, l_項目.A_k長);
            Assert.Equal(l_正解.Length - l_k長 + 1, l_項目.A_カバレッジ.Length);
            Assert.All(l_項目.A_カバレッジ, x => Assert.Equal(17, x));
        }

        /// <summary>
        /// 引き継いだ k-mer のカバレッジが、k の差ぶんスケールされること
        /// </summary>
        /// <remarks>
        /// スケールしないと、引き継いだ領域だけカバレッジが高く見えてコピー数を過大に推定する
        /// </remarks>
        [Fact]
        public void V_引き継いだカバレッジは1リードあたりのkmer数の比でスケールされる()
        {
            // 前段 k=21 でカバレッジ 100、リード長 150
            // 次段 k=101 では 1 リードあたり 130 本から 50 本へ減るので、
            // 100 * 50 / 130 = 38
            var l_引き継ぎ = new 引き継ぎ配列(V_生成_ランダム配列(300, p_シード: 702), [.. Enumerable.Repeat(100, 280)], 21);

            var l_カバレッジ = KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, p_位置: 0, p_k長: 101, p_リード長: 150);

            Assert.Equal((ulong)Math.Round(100.0D * 50D / 130D), l_カバレッジ);
        }

        /// <summary>
        /// 引き継ぐ k-mer は、それを構成する前段の k-mer の最小値を超えないこと
        /// </summary>
        /// <remarks>
        /// 長い k-mer は短い k-mer をすべて含むので、最も弱い部分より強くはなれない
        /// </remarks>
        [Fact]
        public void V_引き継ぐkmerは構成する前段kmerの最小値を超えない()
        {
            var l_カバレッジ = new int[280];
            Array.Fill(l_カバレッジ, 100);
            l_カバレッジ[10] = 7;

            var l_引き継ぎ = new 引き継ぎ配列(V_生成_ランダム配列(300, p_シード: 703), l_カバレッジ, 21);

            // 位置 0 から k=41 の窓は前段の位置 0..20 を含むので、7 が効く
            var l_弱い部分を含む = KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, p_位置: 0, p_k長: 41, p_リード長: null);
            // 位置 30 の窓は 30..50 なので 7 を含まない
            var l_含まない = KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, p_位置: 30, p_k長: 41, p_リード長: null);

            Assert.Equal(7UL, l_弱い部分を含む);
            Assert.Equal(100UL, l_含まない);
        }

        /// <summary>
        /// 引き継ぎの本題
        /// </summary>
        /// <remarks>
        /// カバレッジが薄くて次の k では観測されなかった領域が、前段の配列から復元されること
        /// </remarks>
        [Fact]
        public void V_次のkでは観測されなかったkmerを前段の配列から復元する()
        {
            const int l_前段のk = 21;
            const int l_次のk = 41;
            var l_正解 = V_生成_ランダム配列(3_000, p_シード: 711);

            // 次の k では中央の領域が観測されていない状況を作る
            var l_左 = l_正解[..1_200];
            var l_右 = l_正解[1_800..];
            using var l_次段 = this.V_構築_インデックス(l_次のk, p_深さ: 20, l_左, l_右);

            // 前段は全体を観測している
            using var l_前段 = this.V_構築_インデックス(l_前段のk, p_深さ: 20, l_正解);
            var l_パス = this.V_書き出し_FASTA("carry.fasta", l_正解);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, l_前段のk);

            var l_塩基列 = l_正解.Select(Util.Get_塩基ID).ToArray();
            var l_中央 = l_塩基列.AsSpan(1_400, l_次のk);
            Assert.False(l_次段.Haskmer(l_中央));

            var l_追加数 = KmerCarryOver.V_引き継ぎ(l_引き継ぎ, l_次段, l_次のk, p_リード長: null);

            Assert.True(l_追加数 > 0);
            Assert.True(l_次段.Haskmer(l_中央));
            Assert.True(l_次段.Get_カバレッジ(l_中央) > 0UL);
        }

        /// <summary>
        /// 既に観測されている k-mer のカバレッジは書き換えないこと
        /// </summary>
        /// <remarks>
        /// 実際のリードによる観測のほうが、前段からの推定より確かである
        /// </remarks>
        [Fact]
        public void V_既に観測済みのkmerのカバレッジは書き換えない()
        {
            const int l_前段のk = 21;
            const int l_次のk = 31;
            var l_正解 = V_生成_ランダム配列(2_000, p_シード: 721);

            using var l_次段 = this.V_構築_インデックス(l_次のk, p_深さ: 40, l_正解);
            using var l_前段 = this.V_構築_インデックス(l_前段のk, p_深さ: 5, l_正解);

            var l_塩基列 = l_正解.Select(Util.Get_塩基ID).ToArray();
            var l_観測済み = l_塩基列.AsSpan(100, l_次のk);
            var l_元のカバレッジ = l_次段.Get_カバレッジ(l_観測済み);

            var l_パス = this.V_書き出し_FASTA("nooverwrite.fasta", l_正解);
            _ = KmerCarryOver.V_引き継ぎ(KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, l_前段のk), l_次段, l_次のk, p_リード長: null);

            Assert.Equal(l_元のカバレッジ, l_次段.Get_カバレッジ(l_観測済み));
        }

        /// <summary>
        /// 2 コピーの反復配列は、引き継いでも 2 コピー相当のカバレッジを保つこと
        /// </summary>
        /// <remarks>
        /// ここが崩れるとコピー数推定が壊れ、反復配列の扱いが総崩れになる
        /// </remarks>
        [Fact]
        public void V_反復配列の相対的なカバレッジ比を引き継ぎ後も保つ()
        {
            const int l_前段のk = 21;
            const int l_次のk = 41;
            var l_単一 = V_生成_ランダム配列(1_500, p_シード: 731);
            var l_反復 = V_生成_ランダム配列(1_500, p_シード: 732);

            // 反復側だけ倍の深さで観測されている状況
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_前段のk, A_スレッド数 = 1 };
            var l_前段の作業 = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_前段の作業);
            using var l_前段 = new TrustedKmerIndex(l_前段の作業);
            foreach (var (l_配列, l_深さ) in new[] { (l_単一, 20), (l_反復, 40) })
            {
                var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + l_前段のk <= l_塩基列.Length; i++)
                {
                    for (var l_回 = 0; l_回 < l_深さ; l_回++)
                    {
                        l_前段.V_登録(l_塩基列.AsSpan(i, l_前段のk));
                    }
                }
            }
            _ = l_前段.V_カットオフ(p_カットオフ: 2UL);

            // 次の k は何も観測していない
            using var l_次段 = this.V_構築_インデックス(l_次のk, p_深さ: 20, V_生成_ランダム配列(2_000, p_シード: 733));

            var l_パス = this.V_書き出し_FASTA("repeat_carry.fasta", l_単一, l_反復);
            _ = KmerCarryOver.V_引き継ぎ(KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, l_前段のk), l_次段, l_次のk, p_リード長: null);

            var l_単一の位置 = l_単一.Select(Util.Get_塩基ID).ToArray().AsSpan(50, l_次のk);
            var l_反復の位置 = l_反復.Select(Util.Get_塩基ID).ToArray().AsSpan(50, l_次のk);

            // 2 倍の関係が保たれていること
            Assert.Equal(2UL * l_次段.Get_カバレッジ(l_単一の位置), l_次段.Get_カバレッジ(l_反復の位置));
        }

        /// <summary>
        /// 短い断片は引き継がないこと
        /// </summary>
        /// <remarks>
        /// 連結の役に立たないうえ、エラー由来の残骸である可能性が相対的に高い
        /// </remarks>
        [Fact]
        public void V_短い断片は引き継がない()
        {
            const int l_k長 = 21;
            var l_長い = V_生成_ランダム配列(2_000, p_シード: 741);
            var l_短い = V_生成_ランダム配列(120, p_シード: 742);
            using var l_インデックス = this.V_構築_インデックス(l_k長, p_深さ: 20, l_長い, l_短い);

            var l_パス = this.V_書き出し_FASTA("short.fasta", l_長い, l_短い);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(l_パス, l_インデックス, l_k長);

            var l_項目 = Assert.Single(l_引き継ぎ);
            Assert.Equal(l_長い, l_項目.A_配列);
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
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_深さ">登録する深さ</param>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_インデックス(int p_k長, int p_深さ, params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };
            // インデックスごとに作業ディレクトリを分ける (同じ場所を使うと
            // 一時ファイルの後始末が互いに干渉する)
            var l_作業 = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業);
            var l_インデックス = new TrustedKmerIndex(l_作業);
            foreach (var l_配列 in p_配列群)
            {
                var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    for (var l_回 = 0; l_回 < p_深さ; l_回++)
                    {
                        l_インデックス.V_登録(l_塩基列.AsSpan(i, p_k長));
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
        private string V_書き出し_FASTA(string p_ファイル名, params string[] p_配列群)
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

        #endregion

    }
}
