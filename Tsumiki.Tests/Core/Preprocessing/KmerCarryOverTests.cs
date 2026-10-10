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
    public class KmerCarryOverTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki carryover tests
        /// </summary>
        private const string C_項目_tsumiki_carryover_tests = "tsumiki_carryover_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 prep fasta
        /// </summary>
        private const string C_ファイル名_prep_fasta = "prep.fasta";

        /// <summary>
        /// ファイル名 carry fasta
        /// </summary>
        private const string C_ファイル名_carry_fasta = "carry.fasta";

        /// <summary>
        /// ファイル名 nooverwrite fasta
        /// </summary>
        private const string C_ファイル名_nooverwrite_fasta = "nooverwrite.fasta";

        /// <summary>
        /// ファイル名 repeat carry fasta
        /// </summary>
        private const string C_ファイル名_repeat_carry_fasta = "repeat_carry.fasta";

        /// <summary>
        /// ファイル名 junction fasta
        /// </summary>
        private const string C_ファイル名_junction_fasta = "junction.fasta";

        /// <summary>
        /// ファイル名 short fasta
        /// </summary>
        private const string C_ファイル名_short_fasta = "short.fasta";

        /// <summary>
        /// ファイル名 cover fq
        /// </summary>
        private const string C_ファイル名_cover_fq = "cover.fq";

        /// <summary>
        /// ファイル名 tail fq
        /// </summary>
        private const string C_ファイル名_tail_fq = "tail.fq";

        /// <summary>
        /// 書式 r0
        /// </summary>
        private const string C_書式_r0 = "@r0\n";

        /// <summary>
        /// 書式
        /// </summary>
        private const string C_書式 = "\n+\n";

        /// <summary>
        /// 改行
        /// </summary>
        private const string C_改行 = "\n";

        /// <summary>
        /// 項目 s
        /// </summary>
        private const string C_項目_s = "s";

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
        public KmerCarryOverTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_carryover_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// リードで観測されていない範囲に掛かる k-mer は引き継がない
        /// </summary>
        [Fact]
        public void V_除外_未観測の範囲_範囲に掛かる窓を外す()
        {
            const int l_前k長 = 21;
            const int l_k長 = 31;
            var l_カバレッジ = Enumerable.Repeat(40, 200).ToArray();
            var l_引き継ぎ = new 引き継ぎ配列(new string('A', 220), l_カバレッジ, l_前k長, A_未観測の連続範囲: [(100, 113)]);
            var l_最小値列 = Enumerable.Repeat(40, 190).ToArray();
            KmerCarryOver.V_除外_未観測の範囲(l_引き継ぎ, l_k長, l_最小値列);
            Assert.Equal(0, l_最小値列[100 - l_k長 + 1]);
            Assert.Equal(0, l_最小値列[113]);
            Assert.Equal(40, l_最小値列[100 - l_k長]);
            Assert.Equal(40, l_最小値列[113 + KmerCarryOver.C_持ち越し検証のr長]);
        }

        /// <summary>
        /// 未観測の範囲が無ければカバレッジは変わらない
        /// </summary>
        [Fact]
        public void V_除外_未観測の範囲_範囲が無ければ変えない()
        {
            var l_引き継ぎ = new 引き継ぎ配列(new string('A', 120), Enumerable.Repeat(40, 100).ToArray(), 21);
            var l_最小値列 = Enumerable.Repeat(40, 90).ToArray();
            KmerCarryOver.V_除外_未観測の範囲(l_引き継ぎ, 31, l_最小値列);
            Assert.All(l_最小値列, x => Assert.Equal(40, x));
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
            var l_パス = this.V_書き出し_FASTA(C_ファイル名_prep_fasta, l_正解);
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
        [Fact]
        public void V_引き継いだカバレッジは1リードあたりのkmer数の比でスケールされる()
        {
            var l_引き継ぎ = new 引き継ぎ配列(V_生成_ランダム配列(300, p_シード: 702), [..Enumerable.Repeat(100, 280)], 21);
            var l_カバレッジ = KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, p_位置: 0, p_k長: 101, p_リード長: 150);
            Assert.Equal((ulong)Math.Round(100.0D * 50D / 130D), l_カバレッジ);
        }

        /// <summary>
        /// 引き継ぐ k-mer は、それを構成する前段の k-mer の最小値を超えないこと
        /// </summary>
        [Fact]
        public void V_引き継ぐkmerは構成する前段kmerの最小値を超えない()
        {
            var l_カバレッジ = new int[280];
            Array.Fill(l_カバレッジ, 100);
            l_カバレッジ[10] = 7;
            var l_引き継ぎ = new 引き継ぎ配列(V_生成_ランダム配列(300, p_シード: 703), l_カバレッジ, 21);
            var l_弱い部分を含む = KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, p_位置: 0, p_k長: 41, p_リード長: null);
            var l_含まない = KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, p_位置: 30, p_k長: 41, p_リード長: null);
            Assert.Equal(7UL, l_弱い部分を含む);
            Assert.Equal(100UL, l_含まない);
        }

        /// <summary>
        /// 引き継ぎの本題
        /// </summary>
        [Fact]
        public void V_次のkでは観測されなかったkmerを前段の配列から復元する()
        {
            const int l_前段のk = 21;
            const int l_次のk = 41;
            var l_正解 = V_生成_ランダム配列(3_000, p_シード: 711);
            var l_左 = l_正解[..1_200];
            var l_右 = l_正解[1_800..];
            using var l_次段 = this.V_構築_インデックス(l_次のk, p_深さ: 20, l_左, l_右);
            using var l_前段 = this.V_構築_インデックス(l_前段のk, p_深さ: 20, l_正解);
            var l_パス = this.V_書き出し_FASTA(C_ファイル名_carry_fasta, l_正解);
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
            var l_パス = this.V_書き出し_FASTA(C_ファイル名_nooverwrite_fasta, l_正解);
            _ = KmerCarryOver.V_引き継ぎ(KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, l_前段のk), l_次段, l_次のk, p_リード長: null);
            Assert.Equal(l_元のカバレッジ, l_次段.Get_カバレッジ(l_観測済み));
        }

        /// <summary>
        /// 2 コピーの反復配列は、引き継いでも 2 コピー相当のカバレッジを保つこと
        /// </summary>
        [Fact]
        public void V_反復配列の相対的なカバレッジ比を引き継ぎ後も保つ()
        {
            const int l_前段のk = 21;
            const int l_次のk = 41;
            var l_単一 = V_生成_ランダム配列(1_500, p_シード: 731);
            var l_反復 = V_生成_ランダム配列(1_500, p_シード: 732);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_前段のk,
                A_スレッド数 = 1
            };
            var l_前段の作業 = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(l_前段の作業);
            using var l_前段 = new TrustedKmerIndex(l_前段の作業);
            foreach (var (l_配列, l_深さ)in new[]
            {
                (l_単一, 20),
                (l_反復, 40)
            }

            )
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
            using var l_次段 = this.V_構築_インデックス(l_次のk, p_深さ: 20, V_生成_ランダム配列(2_000, p_シード: 733));
            var l_パス = this.V_書き出し_FASTA(C_ファイル名_repeat_carry_fasta, l_単一, l_反復);
            _ = KmerCarryOver.V_引き継ぎ(KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, l_前段のk), l_次段, l_次のk, p_リード長: null);
            var l_単一の位置 = l_単一.Select(Util.Get_塩基ID).ToArray().AsSpan(50, l_次のk);
            var l_反復の位置 = l_反復.Select(Util.Get_塩基ID).ToArray().AsSpan(50, l_次のk);
            Assert.Equal(2UL * l_次段.Get_カバレッジ(l_単一の位置), l_次段.Get_カバレッジ(l_反復の位置));
        }

        /// <summary>
        /// 前段 k の分岐のある継ぎ目の辺を丸ごと含む k-mer は引き継がず、含まない k-mer は引き継ぐこと
        /// </summary>
        [Fact]
        public void V_分岐の継ぎ目を丸ごと含むkmerは引き継がない()
        {
            const int l_前段のk = 21;
            const int l_次のk = 41;
            const int l_継ぎ目の位置 = 1_500;
            var l_正解 = V_生成_ランダム配列(3_000, p_シード: 751);
            using var l_次段 = this.V_構築_インデックス(l_次のk, p_深さ: 20, l_正解[..1_200], l_正解[1_800..]);
            using var l_前段 = this.V_構築_インデックス(l_前段のk, p_深さ: 20, l_正解);
            var l_辺 = l_正解.Substring(l_継ぎ目の位置, l_前段のk + 1);
            HashSet<string> l_継ぎ目 = new(StringComparer.Ordinal)
            {
                l_辺,
                Util.V_逆相補(l_辺)
            };
            var l_パス = this.V_書き出し_FASTA(C_ファイル名_junction_fasta, l_正解);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(l_パス, l_前段, l_前段のk, l_継ぎ目);
            Assert.Equal([l_継ぎ目の位置], Assert.Single(l_引き継ぎ).A_分岐の継ぎ目位置);
            _ = KmerCarryOver.V_引き継ぎ(l_引き継ぎ, l_次段, l_次のk, p_リード長: null);
            var l_塩基列 = l_正解.Select(Util.Get_塩基ID).ToArray();
            var l_辺を含む最初 = l_継ぎ目の位置 + l_前段のk + 1 - l_次のk;
            Assert.True(l_次段.Haskmer(l_塩基列.AsSpan(l_辺を含む最初 - 1, l_次のk)));
            Assert.False(l_次段.Haskmer(l_塩基列.AsSpan(l_辺を含む最初, l_次のk)));
            Assert.False(l_次段.Haskmer(l_塩基列.AsSpan(l_継ぎ目の位置, l_次のk)));
            Assert.True(l_次段.Haskmer(l_塩基列.AsSpan(l_継ぎ目の位置 + 1, l_次のk)));
        }

        /// <summary>
        /// 短い断片は引き継がないこと
        /// </summary>
        [Fact]
        public void V_短い断片は引き継がない()
        {
            const int l_k長 = 21;
            var l_長い = V_生成_ランダム配列(2_000, p_シード: 741);
            var l_短い = V_生成_ランダム配列(120, p_シード: 742);
            using var l_インデックス = this.V_構築_インデックス(l_k長, p_深さ: 20, l_長い, l_短い);
            var l_パス = this.V_書き出し_FASTA(C_ファイル名_short_fasta, l_長い, l_短い);
            var l_引き継ぎ = KmerCarryOver.Get_引き継ぎ配列(l_パス, l_インデックス, l_k長);
            var l_項目 = Assert.Single(l_引き継ぎ);
            Assert.Equal(l_長い, l_項目.A_配列);
        }

        /// <summary>
        /// 全連続を拾ってから下限で絞っても、下限を検証器へ直接渡したときと同じ範囲になること
        /// </summary>
        [Fact]
        public void Get_未観測の連続範囲_下限で絞った範囲は検証器へ直接渡した場合と同じ()
        {
            const int l_r長 = KmerCarryOver.C_持ち越し検証のr長;
            var l_配列 = V_生成_ランダム配列(500, p_シード: 2301);
            string[] l_リード群 = [l_配列[..150], l_配列[115..300], l_配列[280..]];
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_cover_fq);
            File.WriteAllText(l_パス, string.Concat(l_リード群.Select((x, i) => $"@r{i}\n{x}\n+\n{new string('I', x.Length)}\n")));
            var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], l_r長);
            List<(int A_開始, int A_終了)> l_全連続 = [];
            l_検証器.V_収集_未観測の連続範囲(l_配列, 1, l_全連続);
            Assert.Equal([5, 20], l_全連続.Select(x => x.A_終了 - x.A_開始 + 1));
            List<(int A_開始, int A_終了)> l_直接 = [];
            l_検証器.V_収集_未観測の連続範囲(l_配列, KmerCarryOver.C_未観測の連続の下限, l_直接);
            var l_絞った = KmerCarryOver.Get_未観測の連続範囲(l_配列, l_検証器, p_度数: new 連続長の度数());
            Assert.NotNull(l_絞った);
            Assert.Equal(l_直接, l_絞った);
        }

        /// <summary>
        /// k が r より長く、未観測の連続が配列の末尾に掛かっても落ちないこと
        /// </summary>
        [Fact]
        public void Get_未観測の連続範囲_kがrより長く末尾に連続があっても落ちない()
        {
            const int l_k長 = 67;
            var l_配列 = V_生成_ランダム配列(400, p_シード: 2302);
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_tail_fq);
            var l_リード = l_配列[..380];
            File.WriteAllText(l_パス, C_書式_r0 + l_リード + C_書式 + new string('I', l_リード.Length) + C_改行);
            var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], KmerCarryOver.C_持ち越し検証のr長);
            var l_カバレッジ = Enumerable.Repeat(10, l_配列.Length - l_k長 + 1).ToArray();
            var l_度数 = new 連続長の度数();
            var l_範囲 = KmerCarryOver.Get_未観測の連続範囲(l_配列, l_検証器, null, l_k長, l_度数, C_項目_s, l_カバレッジ);
            var l_末尾 = Assert.Single(l_範囲!);
            Assert.Equal(l_配列.Length - KmerCarryOver.C_持ち越し検証のr長, l_末尾.A_終了);
            Assert.True(l_末尾.A_開始 >= l_カバレッジ.Length);
        }

        /// <summary>
        /// 窓の差分更新後も元の最小カバレッジと既存キーの優先を保つ
        /// </summary>
        /// <param name="p_k長">引き継ぎ先の k</param>
        [Theory]
        [InlineData(29)]
        [InlineData(43)]
        [InlineData(63)]
        [InlineData(93)]
        [InlineData(135)]
        [InlineData(160)]
        [InlineData(161)]
        public void V_高速引き継ぎを窓ごとの計算と比較(int p_k長)
        {
            var l_配列 = V_生成_ランダム配列(500, 1901) + C_GUID書式 + V_生成_ランダム配列(300, 1902);
            var l_カバレッジ = Enumerable.Range(0, l_配列.Length - 20).Select(i => 20 + i % 29).ToArray();
            l_カバレッジ[180] = 0;
            var l_引き継ぎ = new 引き継ぎ配列(l_配列, l_カバレッジ, 21);
            using var l_索引 = this.V_構築_インデックス(p_k長, 40, l_配列[..p_k長]);
            _ = KmerCarryOver.V_引き継ぎ([l_引き継ぎ], l_索引, p_k長, 150);
            var l_塩基 = Util.V_変換_塩基列(l_配列);
            for (var i = 0; i + p_k長 <= l_塩基.Length; i++)
            {
                if (l_配列.Substring(i, p_k長).Contains('N'))
                {
                    continue;
                }

                var l_期待 = i == 0 ? 40UL : KmerCarryOver.Get_引き継ぐカバレッジ(l_引き継ぎ, i, p_k長, 150);
                Assert.Equal(l_期待, l_索引.Get_カバレッジ(l_塩基.AsSpan(i, p_k長)));
            }
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
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
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
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = p_k長,
                A_スレッド数 = 1
            };
            var l_作業 = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString(C_GUID書式));
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
