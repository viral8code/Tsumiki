using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 置き場の有効・無効を切り替えるテストを、他のテストと並べて走らせないための集まり
    /// </summary>
    /// <remarks>
    /// 切り替えはプロセス全体に効くので、並列に走る他のテストの中間ファイルまでメモリへ回ってしまう
    /// </remarks>
    [CollectionDefinition(名前, DisableParallelization = true)]
    public class 中間データ置き場の集まり
    {
        /// <summary>
        /// 集まりの名前
        /// </summary>
        public const string 名前 = "中間データ置き場";
    }

    /// <summary>
    /// 中間データをメモリに置いても、ディスクに置いたときと同じ中身を読み書きできることの検証
    /// </summary>
    [Collection(中間データ置き場の集まり.名前)]
    public class 中間データ置き場Tests : IDisposable
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
        public 中間データ置き場Tests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_intermediate_store_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 片付ける
        /// </summary>
        public void Dispose()
        {
            中間データ置き場.A_Is有効 = false;
            中間データ置き場.V_消去();
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 塊の境目をまたぐ大きさでも、書いた中身がそのまま読め、位置と長さで終端を判定できること
        /// </summary>
        /// <remarks>
        /// k-mer 計数の統合は「位置 &lt; 長さ」で終端を判定するので、圧縮しても展開後の長さと位置を答える必要がある
        /// </remarks>
        [Fact]
        public void 塊をまたぐ大きさでも同じ中身が読める()
        {
            中間データ置き場.A_Is有効 = true;
            var l_パス = Path.Combine(this._作業ディレクトリ, "big.bin");
            var l_中身 = new byte[(40 << 20) + 12345];
            new Random(2401).NextBytes(l_中身);

            using (var l_書き込み = 中間データ置き場.Get_書込ストリーム(l_パス))
            {
                l_書き込み.Write(l_中身, 0, 1000);
                l_書き込み.Write(l_中身, 1000, l_中身.Length - 1000);
            }

            Assert.False(File.Exists(l_パス));
            Assert.True(中間データ置き場.Is存在(l_パス));

            using var l_読み込み = new BinaryReader(中間データ置き場.Get_読込ストリーム(l_パス));
            Assert.Equal(l_中身.Length, l_読み込み.BaseStream.Length);
            var l_読んだ = new List<byte>();
            while (Util.Has続き(l_読み込み))
            {
                l_読んだ.AddRange(l_読み込み.ReadBytes(7777));
            }
            Assert.Equal(l_中身, l_読んだ.ToArray());
        }

        /// <summary>
        /// 列挙・削除が置き場とディスクの両方に効くこと
        /// </summary>
        [Fact]
        public void 列挙と削除が置き場とディスクの両方に効く()
        {
            var l_ディスク = Path.Combine(this._作業ディレクトリ, "preprocessed.1.fq");
            File.WriteAllText(l_ディスク, "x");
            中間データ置き場.A_Is有効 = true;
            var l_メモリ = Path.Combine(this._作業ディレクトリ, "preprocessed.lib2.1.fq");
            using (var l_書き込み = new FastqWriter(l_メモリ))
            {
                l_書き込み.V_書き込み("r1", "ACGT", "IIII");
            }
            var l_対象外 = Path.Combine(this._作業ディレクトリ, "corrected.1.fq");
            using (var l_書き込み = new FastqWriter(l_対象外))
            {
                l_書き込み.V_書き込み("r1", "ACGT", "IIII");
            }

            var l_一覧 = 中間データ置き場.Get_一覧(this._作業ディレクトリ, "preprocessed", ".fq");
            Assert.Equal(2, l_一覧.Count);

            foreach (var l_パス in l_一覧)
            {
                中間データ置き場.V_削除(l_パス);
            }
            Assert.False(中間データ置き場.Is存在(l_ディスク));
            Assert.False(中間データ置き場.Is存在(l_メモリ));
            Assert.True(中間データ置き場.Is存在(l_対象外));
        }

        /// <summary>
        /// FASTQ を書いて読み直した結果が、ディスクに置いたときとメモリに置いたときで同じになること
        /// </summary>
        [Fact]
        public void FASTQの読み書きがディスクと同じ()
        {
            var l_リード群 = Enumerable.Range(0, 5000).Select(i => ($"read{i}", Get_ランダム配列(150, i), new string((char)('!' + (i % 40)), 150))).ToList();

            var l_ディスク = V_書いて読む(false, "disk.fq", l_リード群);
            var l_メモリ = V_書いて読む(true, "memory.fq", l_リード群);

            Assert.Equal(l_ディスク, l_メモリ);
            Assert.Equal(l_リード群.Count, l_メモリ.Count);
        }

        /// <summary>
        /// 計数の途中結果を外へ逃がすほど予算を絞っても、メモリに置いたときとディスクに置いたときで同じ k-mer と回数になること
        /// </summary>
        /// <remarks>
        /// 途中結果の書き出し・統合・削除の全部が置き場を通るかを確かめる。どこか 1 つでもディスクを直接触ると、ファイルが見つからずに落ちるか数がずれる
        /// </remarks>
        [Fact]
        public void 計数の途中結果をメモリに置いてもディスクと同じ結果になる()
        {
            var l_ディスク = this.Get_計数結果(false);
            var l_メモリ = this.Get_計数結果(true);

            Assert.NotEmpty(l_ディスク.A_信頼kmer);
            Assert.Equal(l_ディスク.A_信頼kmer, l_メモリ.A_信頼kmer);
            Assert.Equal(l_ディスク.A_ヒストグラム, l_メモリ.A_ヒストグラム);
            Assert.DoesNotContain(Directory.EnumerateFiles(this._作業ディレクトリ, "*", SearchOption.AllDirectories), x => x.Contains("memory_run", StringComparison.Ordinal));
        }

        /// <summary>
        /// 取り込んだ入力はメモリから読め、置き場から消してもディスク上の入力は残ること
        /// </summary>
        /// <remarks>
        /// 入力リードは工程ごとに読み直されるので、メモリに置くときは取り込んでおく。取り込みが入力そのものを壊してはいけない
        /// </remarks>
        [Fact]
        public void 取り込んだ入力はメモリから読めディスクの入力は消えない()
        {
            var l_リード群 = Enumerable.Range(0, 300).Select(i => ($"read{i}", Get_ランダム配列(100, i), new string('I', 100))).ToList();
            var l_ディスク = V_書いて読む(false, "input.fq", l_リード群);
            var l_パス = Path.Combine(this._作業ディレクトリ, "input.fq");

            中間データ置き場.A_Is有効 = true;
            中間データ置き場.V_取り込み(l_パス);
            File.Move(l_パス, l_パス + ".moved");

            List<(string, string, string)> l_メモリ = [];
            using (var l_読み込み = new FastqReader(l_パス))
            {
                while (l_読み込み.Has続き())
                {
                    l_メモリ.Add(l_読み込み.Get_次のレコード());
                }
            }
            Assert.Equal(l_ディスク, l_メモリ);

            File.Move(l_パス + ".moved", l_パス);
            中間データ置き場.V_削除(l_パス);
            中間データ置き場.V_削除(l_パス);
            Assert.True(File.Exists(l_パス));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// FASTQ を書いて読み直す
        /// </summary>
        /// <param name="p_Isメモリ">置き場を有効にするか</param>
        /// <param name="p_名前">ファイル名</param>
        /// <param name="p_リード群">書くリード</param>
        /// <returns>読み直したリード</returns>
        private List<(string, string, string)> V_書いて読む(bool p_Isメモリ, string p_名前, List<(string A_ID, string A_配列, string A_品質)> p_リード群)
        {
            中間データ置き場.A_Is有効 = p_Isメモリ;
            var l_パス = Path.Combine(this._作業ディレクトリ, p_名前);
            using (var l_書き込み = new FastqWriter(l_パス))
            {
                foreach (var (l_ID, l_配列, l_品質) in p_リード群)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列, l_品質);
                }
            }
            Assert.Equal(!p_Isメモリ, File.Exists(l_パス));

            List<(string, string, string)> l_結果 = [];
            using var l_読み込み = new FastqReader(l_パス);
            while (l_読み込み.Has続き())
            {
                l_結果.Add(l_読み込み.Get_次のレコード());
            }
            中間データ置き場.A_Is有効 = false;
            return l_結果;
        }

        /// <summary>
        /// 予算を絞って k-mer を数え、信頼 k-mer とヒストグラムを返す
        /// </summary>
        /// <param name="p_Isメモリ">置き場を有効にするか</param>
        /// <returns></returns>
        private (List<string> A_信頼kmer, List<(ulong, long)> A_ヒストグラム) Get_計数結果(bool p_Isメモリ)
        {
            中間データ置き場.A_Is有効 = p_Isメモリ;
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 2, A_メモリ予算 = "1M" };
            var l_作業 = Path.Combine(this._作業ディレクトリ, p_Isメモリ ? "memory_run" : "disk_run");
            _ = Directory.CreateDirectory(l_作業);

            var l_ゲノム = Get_ランダム配列(20_000, 2402);
            using var l_索引 = new TrustedKmerIndex(l_作業);
            var l_塩基 = l_ゲノム.Select(Util.Get_塩基ID).ToArray();
            for (var l_回 = 0; l_回 < 3; l_回++)
            {
                for (var i = 0; i + 21 <= l_塩基.Length; i++)
                {
                    l_索引.V_登録(l_塩基.AsSpan(i, 21));
                }
            }

            // 予算が足りずに途中結果を外へ逃がしていなければ、置き場を通る経路を確かめたことにならない
            if (p_Isメモリ)
            {
                Assert.True(中間データ置き場.A_使用量 > 0);
            }
            else
            {
                Assert.NotEmpty(Directory.EnumerateFiles(l_作業));
            }
            _ = l_索引.V_カットオフ(2UL);

            var l_信頼kmer = l_索引.Get_信頼kmer一覧().Select(Convert.ToHexString).Order(StringComparer.Ordinal).ToList();
            var l_ヒストグラム = l_索引.A_出現回数ヒストグラム.OrderBy(x => x.Key).Select(x => (x.Key, x.Value)).ToList();
            中間データ置き場.A_Is有効 = false;
            return (l_信頼kmer, l_ヒストグラム);
        }

        /// <summary>
        /// 決まった種から作る塩基列
        /// </summary>
        /// <param name="p_長さ"></param>
        /// <param name="p_種"></param>
        /// <returns></returns>
        private static string Get_ランダム配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return new string([.. Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)])]);
        }

        #endregion
    }
}
