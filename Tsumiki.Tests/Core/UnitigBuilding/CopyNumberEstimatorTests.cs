using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// カバレッジから unitig のコピー数を推定する処理の検証
    /// </summary>
    /// <remarks>
    /// ゲノム中に 1 回しか現れない領域のカバレッジを基準値とすると、n 回現れる反復配列にはリードが n 倍集まる<br/>
    /// したがってカバレッジ比を丸めればコピー数になる<br/>
    /// これが分かると、反復配列かどうかをグラフの形ではなく量的な根拠で判定でき、経路探索では「何回まで使ってよいか」の予算になる
    /// </remarks>
    public class CopyNumberEstimatorTests : IDisposable
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
        public CopyNumberEstimatorTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_copynumber_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);

            // これらのテストは長さ加重中央値のフォールバック経路を検証する
            // 他のテスト (KmerCutoffSelectorTests 等) が残した混合モデルの
            // 適合結果が ConfigurationManager 経由で漏れ込まないようにする
            ConfigurationManager.A_スペクトルモデル = null;
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
        /// 単一コピーと 2 倍・ 4 倍コピーの配列をカバレッジ比から分離できる
        /// </summary>
        [Fact]
        public void V_単一コピーと2倍_4倍コピーの配列を分離できる()
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };

            // 単一コピー相当を 2 本 (長さで基準値を支配させる)、
            // 2 倍・ 4 倍のカバレッジで登録する配列を 1 本ずつ用意する
            var l_先行単一配列 = V_生成_乱数配列(400, p_乱数種: 1);
            var l_後続単一配列 = V_生成_乱数配列(400, p_乱数種: 2);
            var l_二重配列 = V_生成_乱数配列(120, p_乱数種: 3);
            var l_四重配列 = V_生成_乱数配列(120, p_乱数種: 4);

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            V_登録_全kmer(l_インデックス, l_先行単一配列, 20, l_k長);
            V_登録_全kmer(l_インデックス, l_後続単一配列, 20, l_k長);
            V_登録_全kmer(l_インデックス, l_二重配列, 40, l_k長);
            V_登録_全kmer(l_インデックス, l_四重配列, 80, l_k長);

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Dictionary<int, string> l_ユニティグ群 = new()
            {
                [1] = l_先行単一配列,
                [2] = l_後続単一配列,
                [3] = l_二重配列,
                [4] = l_四重配列,
            };
            var l_長さ一覧 = l_ユニティグ群.ToDictionary(l_組 => l_組.Key, l_組 => l_組.Value.Length);

            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_インデックス, l_ユニティグ群, l_k長);
            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧);

            // 基準値は長さ加重中央値なので、長い単一コピー配列の水準になるはず
            Assert.InRange(l_結果.A_単一コピー基準値, 15D, 25D);

            Assert.Equal(1, l_結果.A_コピー数[1]);
            Assert.Equal(1, l_結果.A_コピー数[2]);
            Assert.Equal(2, l_結果.A_コピー数[3]);
            Assert.Equal(4, l_結果.A_コピー数[4]);
        }

        /// <summary>
        /// カバレッジがわずかに高いだけの配列を反復と誤判定してはいけない
        /// </summary>
        /// <remarks>
        /// 実データのカバレッジは領域ごとにかなりばらつくため、1.5 倍未満は単一コピーとして扱う
        /// </remarks>
        [Fact]
        public void V_わずかに高いカバレッジの配列を単一コピーとして扱う()
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };

            var l_基準配列 = V_生成_乱数配列(400, p_乱数種: 5);
            var l_微増配列 = V_生成_乱数配列(120, p_乱数種: 6);

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            V_登録_全kmer(l_インデックス, l_基準配列, 20, l_k長);
            V_登録_全kmer(l_インデックス, l_微増配列, 26, l_k長); // 1.3倍

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Dictionary<int, string> l_ユニティグ群 = new() { [1] = l_基準配列, [2] = l_微増配列 };
            var l_長さ一覧 = l_ユニティグ群.ToDictionary(l_組 => l_組.Key, l_組 => l_組.Value.Length);

            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_インデックス, l_ユニティグ群, l_k長);
            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧);

            Assert.Equal(1, l_結果.A_コピー数[2]);
        }

        /// <summary>
        /// k-mer 長より短い unitig はカバレッジを測れないが、コピー数 0 にして経路から締め出してはいけない (配列自体は存在する)
        /// </summary>
        [Fact]
        public void V_kmer長より短いunitigはコピー数0でなく1になる()
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };

            var l_通常配列 = V_生成_乱数配列(300, p_乱数種: 8);
            var l_短すぎる配列 = V_生成_乱数配列(10, p_乱数種: 9);

            using var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);
            V_登録_全kmer(l_インデックス, l_通常配列, 20, l_k長);
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Dictionary<int, string> l_ユニティグ群 = new() { [1] = l_通常配列, [2] = l_短すぎる配列 };
            var l_長さ一覧 = l_ユニティグ群.ToDictionary(l_組 => l_組.Key, l_組 => l_組.Value.Length);

            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_インデックス, l_ユニティグ群, l_k長);
            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧);

            Assert.Equal(0D, l_カバレッジ[2]);
            Assert.Equal(1, l_結果.A_コピー数[2]);
        }

        /// <summary>
        /// プラスミドのように染色体とは異なるカバレッジ水準を持つ領域は、大域基準値との比だけで見ると多コピーの反復に見える
        /// </summary>
        /// <remarks>
        /// しかしその単一コピー領域同士は分岐の無い (排他的な) 鎖で繋がっているため、接続構造を使えば「大域とは水準が違うだけの単一コピー」だと分かる<br/>
        /// unicycler の copy depth propagation が解決する問題そのもの
        /// </remarks>
        [Fact]
        public void V_グラフを使うと高カバレッジのプラスミド骨格を単一コピーと認識できる()
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            // 染色体相当 (長さで大域基準値=60 を支配する)
            var l_染色体 = V_生成_乱数配列(400, p_乱数種: 10);

            // プラスミド相当
            // 1 本の配列を k-1 (=20) ずつ重ねて 3 本に切り出し、
            // 分岐の無い鎖 plasmid1 -> plasmid2 -> plasmid3 を作る
            var l_プラスミド全配列 = V_生成_乱数配列(90, p_乱数種: 20);
            var l_プラスミド先頭 = l_プラスミド全配列[..40];
            var l_プラスミド中間 = l_プラスミド全配列[20..60];
            var l_プラスミド末尾 = l_プラスミド全配列[40..90];

            var l_FASTAパス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            using (var l_書き込み = new FastaWriter(l_FASTAパス))
            {
                l_書き込み.V_書き込み(1, l_染色体);
                l_書き込み.V_書き込み(2, l_プラスミド先頭);
                l_書き込み.V_書き込み(3, l_プラスミド中間);
                l_書き込み.V_書き込み(4, l_プラスミド末尾);
            }

            var l_コンティグ構築 = new ContigMaker(l_FASTAパス);
            var l_グラフ = l_コンティグ構築.Get_グラフ();

            Dictionary<int, double> l_カバレッジ = new()
            {
                [1] = 60.0D,
                [2] = 300.0D, // 大域基準値 (60) との比は 5 倍 -> 単独では多コピー判定
                [3] = 300.0D,
                [4] = 300.0D,
            };
            Dictionary<int, int> l_長さ一覧 = new()
            {
                [1] = l_染色体.Length,
                [2] = l_プラスミド先頭.Length,
                [3] = l_プラスミド中間.Length,
                [4] = l_プラスミド末尾.Length,
            };

            var l_グラフなし = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧);
            Assert.Equal(5, l_グラフなし.A_コピー数[2]);
            Assert.Equal(5, l_グラフなし.A_コピー数[3]);
            Assert.Equal(5, l_グラフなし.A_コピー数[4]);

            var l_グラフあり = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧, l_グラフ);
            Assert.Equal(1, l_グラフあり.A_コピー数[2]);
            Assert.Equal(1, l_グラフあり.A_コピー数[3]);
            Assert.Equal(1, l_グラフあり.A_コピー数[4]);
            // 染色体側は元々単一コピー判定であり、接続補正の対象にもならない
            Assert.Equal(1, l_グラフあり.A_コピー数[1]);
        }

        /// <summary>
        /// 排他的に繋がる相手がいない (孤立した) 高カバレッジ unitig は、比較材料が無いため接続補正の対象にせず、大域基準値との比のまま残す
        /// </summary>
        [Fact]
        public void V_グラフを使っても孤立した高カバレッジunitigは補正されない()
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            var l_染色体 = V_生成_乱数配列(400, p_乱数種: 11);
            var l_孤立反復 = V_生成_乱数配列(50, p_乱数種: 21); // 他のどれとも重ならない

            var l_FASTAパス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            using (var l_書き込み = new FastaWriter(l_FASTAパス))
            {
                l_書き込み.V_書き込み(1, l_染色体);
                l_書き込み.V_書き込み(2, l_孤立反復);
            }

            var l_コンティグ構築 = new ContigMaker(l_FASTAパス);
            var l_グラフ = l_コンティグ構築.Get_グラフ();

            Dictionary<int, double> l_カバレッジ = new() { [1] = 60.0D, [2] = 300.0D };
            Dictionary<int, int> l_長さ一覧 = new() { [1] = l_染色体.Length, [2] = l_孤立反復.Length };

            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧, l_グラフ);

            Assert.Equal(5, l_結果.A_コピー数[2]);
        }

        /// <summary>
        /// 小さなプラスミドが分岐無しの 1 本の unitig にきれいに閉じた、もっとも典型的なケース
        /// </summary>
        /// <remarks>
        /// 染色体側の成分とは一切繋がりが無い、十分な長さを持つ「島」なので、大域基準値との比が高くても単一コピーとみなしてよい (高コピープラスミド自身の水準で 1 コピー)
        /// </remarks>
        [Fact]
        public void V_グラフを使うと孤立した長いunitigを単独の単一コピーレプリコンと認識できる()
        {
            const int l_k長 = 21;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = l_k長, A_スレッド数 = 1 };
            ConfigurationManager.A_スペクトルモデル = null;

            var l_染色体 = V_生成_乱数配列(400, p_乱数種: 12);
            var l_プラスミド = V_生成_乱数配列(600, p_乱数種: 22); // 染色体とは無関係、500 bp 超

            var l_FASTAパス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            using (var l_書き込み = new FastaWriter(l_FASTAパス))
            {
                l_書き込み.V_書き込み(1, l_染色体);
                l_書き込み.V_書き込み(2, l_プラスミド);
            }

            var l_コンティグ構築 = new ContigMaker(l_FASTAパス);
            var l_グラフ = l_コンティグ構築.Get_グラフ();

            Dictionary<int, double> l_カバレッジ = new() { [1] = 60.0D, [2] = 300.0D };
            Dictionary<int, int> l_長さ一覧 = new() { [1] = l_染色体.Length, [2] = l_プラスミド.Length };

            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧, l_グラフ);

            Assert.Equal(1, l_結果.A_コピー数[2]);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_乱数種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_長さ, int p_乱数種)
        {
            var l_乱数生成器 = new Random(p_乱数種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数生成器.Next(4)]));
        }

        /// <summary>
        /// 配列の全 k-mer を、指定した深さだけ信頼できる kmer 索引へ登録する
        /// </summary>
        /// <param name="p_インデックス">登録先の索引</param>
        /// <param name="p_配列">登録する配列</param>
        /// <param name="p_深度">各 k-mer を登録する回数</param>
        /// <param name="p_k長">k-mer 長</param>
        private static void V_登録_全kmer(TrustedKmerIndex p_インデックス, string p_配列, int p_深度, int p_k長)
        {
            var l_バイト列 = p_配列.Select(Util.Get_塩基ID).ToArray();
            for (var i = 0; i + p_k長 <= l_バイト列.Length; i++)
            {
                for (var l_反復回数 = 0; l_反復回数 < p_深度; l_反復回数++)
                {
                    p_インデックス.V_登録(l_バイト列.AsSpan(i, p_k長));
                }
            }
        }

        #endregion

    }
}
