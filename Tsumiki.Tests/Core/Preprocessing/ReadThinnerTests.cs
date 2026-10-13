using System.IO.Compression;
using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 深すぎる入力の間引き (採用割合、採否の決め方、ファイルの書き出し、総塩基数) の検証
    /// </summary>
    public class ReadThinnerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki read thinner tests
        /// </summary>
        private const string C_項目_tsumiki_read_thinner_tests = "tsumiki_read_thinner_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 塩基配列 1 文字ぶん
        /// </summary>
        private const string C_塩基 = "A";

        /// <summary>
        /// クオリティ 1 文字ぶん
        /// </summary>
        private const string C_品質 = "I";

        /// <summary>
        /// 比較の許容誤差
        /// </summary>
        private const int C_小数の桁数 = 10;

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業用ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の作業ディレクトリを用意する
        /// </summary>
        public ReadThinnerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_read_thinner_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 目標の 1.2 倍以下の深度なら間引かない (null)
        /// </summary>
        [Fact]
        public void V_目標の1_2倍以下なら間引かない()
        {
            Assert.Null(ReadThinner.Get_採用割合(11000, 100, 100));
        }

        /// <summary>
        /// 深度 500x を目標 100x へ間引くと採用割合は 0.2
        /// </summary>
        [Fact]
        public void V_深度500xを100xへ間引くと採用割合は0_2()
        {
            var l_割合 = ReadThinner.Get_採用割合(50000, 100, 100);
            Assert.NotNull(l_割合);
            Assert.Equal(0.2, l_割合!.Value, C_小数の桁数);
        }

        /// <summary>
        /// 深度 121x を目標 100x へ間引くと採用割合は 100/121
        /// </summary>
        [Fact]
        public void V_深度121xを100xへ間引くと採用割合は100割る121()
        {
            var l_割合 = ReadThinner.Get_採用割合(12100, 100, 100);
            Assert.NotNull(l_割合);
            Assert.Equal(100.0 / 121.0, l_割合!.Value, C_小数の桁数);
        }

        /// <summary>
        /// 目標 0、ゲノムサイズ 0、総塩基数 0 のいずれかなら間引かない (null)
        /// </summary>
        [Fact]
        public void V_目標0とゲノム0と総塩基0は間引かない()
        {
            Assert.Null(ReadThinner.Get_採用割合(50000, 100, 0));
            Assert.Null(ReadThinner.Get_採用割合(50000, 0, 100));
            Assert.Null(ReadThinner.Get_採用割合(0, 100, 100));
        }

        /// <summary>
        /// 割合 0.25 では番号 3 と 7 が採られる
        /// </summary>
        [Fact]
        public void V_割合0_25では3と7が採られる()
        {
            var l_採用 = Enumerable.Range(0, 8).Where(x => ReadThinner.Is採用(x, 0.25)).Select(x => (long)x).ToArray();
            Assert.Equal([3L, 7L], l_採用);
        }

        /// <summary>
        /// 割合 1.0 では全部採られる
        /// </summary>
        [Fact]
        public void V_割合1_0では全部採られる()
        {
            var l_採用 = Enumerable.Range(0, 8).Where(x => ReadThinner.Is採用(x, 1.0)).ToArray();
            Assert.Equal(Enumerable.Range(0, 8).ToArray(), l_採用);
        }

        /// <summary>
        /// 割合 0.5 では番号 1, 3, 5 が採られる
        /// </summary>
        [Fact]
        public void V_割合0_5では1と3と5が採られる()
        {
            var l_採用 = Enumerable.Range(0, 6).Where(x => ReadThinner.Is採用(x, 0.5)).Select(x => (long)x).ToArray();
            Assert.Equal([1L, 3L, 5L], l_採用);
        }

        /// <summary>
        /// 8 組のペアを割合 0.25 で間引くと 2 組が残り、番号 3 と 7 のリードが書かれる
        /// </summary>
        [Fact]
        public void V_ペアを間引くと採られた組だけが書かれる()
        {
            var l_順 = this.Get_パス("r1.fq");
            var l_逆 = this.Get_パス("r2.fq");
            var l_順件数 = Enumerable.Range(0, 8).Select(x => $"r{x}").ToArray();
            var l_逆件数 = Enumerable.Range(0, 8).Select(x => $"s{x}").ToArray();
            this.V_書込_FASTQ(l_順, l_順件数);
            this.V_書込_FASTQ(l_逆, l_逆件数);
            var l_出力1 = this.Get_パス("out1.fq");
            var l_出力2 = this.Get_パス("out2.fq");

            var (l_元, l_採用) = ReadThinner.V_間引き_ファイル(l_順, l_逆, l_出力1, l_出力2, 0.25);

            Assert.Equal(8, l_元);
            Assert.Equal(2, l_採用);
            var l_行1 = File.ReadAllLines(l_出力1);
            var l_行2 = File.ReadAllLines(l_出力2);
            Assert.Equal(8, l_行1.Length);
            Assert.Equal(8, l_行2.Length);
            Assert.Contains("r3", l_行1[0], StringComparison.Ordinal);
            Assert.Contains("r7", l_行1[4], StringComparison.Ordinal);
            Assert.Contains("s3", l_行2[0], StringComparison.Ordinal);
            Assert.Contains("s7", l_行2[4], StringComparison.Ordinal);
        }

        /// <summary>
        /// 単一リード 4 件を割合 0.5 で間引くと 2 件が残り、番号 1 と 3 が書かれる
        /// </summary>
        [Fact]
        public void V_単一リードを間引くと採られた件だけが書かれる()
        {
            var l_入力 = this.Get_パス("single.fq");
            this.V_書込_FASTQ(l_入力, ["r0", "r1", "r2", "r3"]);
            var l_出力 = this.Get_パス("single_out.fq");

            var (l_元, l_採用) = ReadThinner.V_間引き_ファイル(l_入力, null, l_出力, null, 0.5);

            Assert.Equal(4, l_元);
            Assert.Equal(2, l_採用);
            var l_行 = File.ReadAllLines(l_出力);
            Assert.Equal(8, l_行.Length);
            Assert.Contains("r1", l_行[0], StringComparison.Ordinal);
            Assert.Contains("r3", l_行[4], StringComparison.Ordinal);
        }

        /// <summary>
        /// 圧縮された入力でも、同じ件数と同じ採用結果になる
        /// </summary>
        [Fact]
        public void V_圧縮入力でも同じ結果になる()
        {
            var l_平文 = this.Get_パス("plain_src.fq");
            this.V_書込_FASTQ(l_平文, ["r0", "r1", "r2", "r3"]);
            var l_入力 = this.Get_パス("gz_src.fq.gz");
            using (var l_読 = File.OpenRead(l_平文))
            using (var l_書 = File.Create(l_入力))
            using (var l_圧縮 = new GZipStream(l_書, CompressionMode.Compress))
            {
                l_読.CopyTo(l_圧縮);
            }

            var l_出力 = this.Get_パス("gz_out.fq");
            var (l_元, l_採用) = ReadThinner.V_間引き_ファイル(l_入力, null, l_出力, null, 0.5);

            Assert.Equal(4, l_元);
            Assert.Equal(2, l_採用);
            var l_行 = File.ReadAllLines(l_出力);
            Assert.Contains("r1", l_行[0], StringComparison.Ordinal);
            Assert.Contains("r3", l_行[4], StringComparison.Ordinal);
        }

        /// <summary>
        /// 総塩基数は全ライブラリの全リードの塩基数の和になる (ペア 10 と 12、単一 5 で 27)
        /// </summary>
        [Fact]
        public void V_総塩基数は全リードの塩基数の和になる()
        {
            var l_順 = this.Get_パス("len_r1.fq");
            var l_逆 = this.Get_パス("len_r2.fq");
            var l_単一 = this.Get_パス("len_single.fq");
            this.V_書込_長さ指定(l_順, 10);
            this.V_書込_長さ指定(l_逆, 12);
            this.V_書込_長さ指定(l_単一, 5);
            var l_引数 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_順]);
            l_引数.Set_ライブラリ群([(l_順, l_逆), (l_単一, string.Empty)]);

            Assert.Equal(27, ReadThinner.Get_総塩基数(l_引数));
        }

        /// <summary>
        /// 作業ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 採用割合 0.3456789 を thinned.info に書いて読むと、同じ値が返る
        /// </summary>
        [Fact]
        public void V_採用割合を書いて読むと同じ値になり_無いファイルはnull()
        {
            var l_パス = Path.Combine(Path.GetTempPath(), $"{C_項目_tsumiki_read_thinner_tests}info_{Guid.NewGuid():N}.info");
            try
            {
                Assert.Null(ReadThinner.Get_読込_採用割合(l_パス));
                ReadThinner.V_書き出し_採用割合(l_パス, 0.3456789);
                Assert.Equal(0.3456789, ReadThinner.Get_読込_採用割合(l_パス));
            }
            finally
            {
                if (File.Exists(l_パス))
                {
                    File.Delete(l_パス);
                }
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 作業ディレクトリ内のパスを返す
        /// </summary>
        /// <param name="p_名前">ファイル名</param>
        /// <returns>パス</returns>
        private string Get_パス(string p_名前)
        {
            return Path.Combine(this._作業ディレクトリ, p_名前);
        }

        /// <summary>
        /// 指定の ID 群で、各 4 行の FASTQ を書く
        /// </summary>
        /// <param name="p_パス">書く先</param>
        /// <param name="p_ID群">リードの ID</param>
        private void V_書込_FASTQ(string p_パス, IEnumerable<string> p_ID群)
        {
            using var l_書込 = new FastqWriter(p_パス);
            foreach (var l_ID in p_ID群)
            {
                l_書込.V_書き込み(l_ID, C_塩基, C_品質);
            }
        }

        /// <summary>
        /// 指定の長さの塩基列を 1 件だけ持つ FASTQ を書く
        /// </summary>
        /// <param name="p_パス">書く先</param>
        /// <param name="p_長さ">塩基列の長さ</param>
        private void V_書込_長さ指定(string p_パス, int p_長さ)
        {
            using var l_書込 = new FastqWriter(p_パス);
            l_書込.V_書き込み("r", new string(C_塩基[0], p_長さ), new string(C_品質[0], p_長さ));
        }

        #endregion
    }
}
