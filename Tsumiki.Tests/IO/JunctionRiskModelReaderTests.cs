using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 継ぎ目の誤り確率モデルの既定の値、ファイルの読み込み、-jm と実行ファイルの隣のファイルでの差し替えの検証
    /// </summary>
    public class JunctionRiskModelReaderTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki junction model tests
        /// </summary>
        private const string C_項目_tsumiki_junction_model_tests = "tsumiki_junction_model_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 組み込みのモデルの版の名前
        /// </summary>
        private const string C_組み込みの名前 = "tsumiki-junction-v1";

        /// <summary>
        /// テストで渡す出所
        /// </summary>
        private const string C_テストの出所 = "test.tsv";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// ファイル名 none tsv
        /// </summary>
        private const string C_ファイル名_none_tsv = "none.tsv";

        /// <summary>
        /// リポジトリの根を見分けるファイル
        /// </summary>
        private const string C_リポジトリの印 = "Tsumiki.csproj";

        /// <summary>
        /// リポジトリの中の雛形の置き場
        /// </summary>
        private const string C_雛形のディレクトリ = "Cores/Evaluation";

        /// <summary>
        /// 注釈の行
        /// </summary>
        private const string C_注釈の行 = "# comment";

        /// <summary>
        /// LF
        /// </summary>
        private const string C_LF = "\n";

        /// <summary>
        /// CRLF
        /// </summary>
        private const string C_CRLF = "\r\n";

        /// <summary>
        /// タブ
        /// </summary>
        private const string C_タブ = "\t";

        /// <summary>
        /// 空白 2 つ
        /// </summary>
        private const string C_空白2つ = "  ";

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_junction_model_tests + Guid.NewGuid().ToString(C_GUID書式));

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 作業ディレクトリを作る
        /// </summary>
        public JunctionRiskModelReaderTests()
        {
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 組み込みのモデルは、版の名前と係数の数が決まったとおりで、ファイルに依らない
        /// </summary>
        [Fact]
        public void A_既定_版の名前と係数の数()
        {
            var l_モデル = JunctionRiskModel.A_既定;
            Assert.Equal(C_組み込みの名前, l_モデル.A_名前);
            Assert.Equal(継ぎ目モデル.C_組み込みの出所, l_モデル.A_出所);
            Assert.Equal(64, l_モデル.A_SHA256.Length);
            int[] l_期待する数 = [13, 13, 13, 4, 2, 4];
            int[] l_数 = [l_モデル.A_得点の係数.Length, l_モデル.A_得点の平均.Length, l_モデル.A_得点の幅.Length, l_モデル.A_ギャップの較正.Length, l_モデル.A_跨げる組の較正.Length, l_モデル.A_跨げない組の較正.Length];
            Assert.Equal(l_期待する数, l_数);
            Assert.Equal(-4.352129400596203D, l_モデル.A_得点の係数[0]);
            Assert.Equal(2.617854620679609D, l_モデル.A_跨げる組の較正[1]);
        }

        /// <summary>
        /// リポジトリの雛形 (Cores/Evaluation/JunctionRiskModel.tsv) は、コードの既定と同じ名前・同じ値・同じ SHA-256
        /// </summary>
        [Fact]
        public void Get_読込_雛形はコードの既定と同じ()
        {
            var l_雛形 = JunctionRiskModelReader.Get_読込(Path.Combine(Get_リポジトリの根(), C_雛形のディレクトリ, Consts.継ぎ目モデルのファイル名));
            var l_既定 = JunctionRiskModel.A_既定;
            Assert.Equal(l_既定.A_名前, l_雛形.A_名前);
            foreach (var ((l_項目, l_既定の値), (_, l_雛形の値)) in l_既定.Get_数の項目群().Zip(l_雛形.Get_数の項目群()))
            {
                Assert.True(l_既定の値.SequenceEqual(l_雛形の値), l_項目);
            }

            Assert.Equal(l_既定.A_SHA256, l_雛形.A_SHA256);
        }

        /// <summary>
        /// 改行 (LF・CRLF)、注釈の行、区切り (タブ・空白) の違いで、係数と SHA-256 は変わらない
        /// </summary>
        [Fact]
        public void Get_解釈_書き方の違いで係数とSHA256は変わらない()
        {
            var l_本文 = JunctionRiskModel.A_既定.Get_本文();
            var l_書き方群 = new[]
            {
                l_本文,
                l_本文.Replace(C_LF, C_CRLF, StringComparison.Ordinal),
                C_注釈の行 + C_LF + l_本文 + C_LF + C_注釈の行,
                l_本文.Replace(C_タブ, C_空白2つ, StringComparison.Ordinal)
            };
            foreach (var l_書き方 in l_書き方群)
            {
                var l_モデル = JunctionRiskModelReader.Get_解釈(l_書き方, C_テストの出所);
                Assert.Equal(JunctionRiskModel.A_既定.A_SHA256, l_モデル.A_SHA256);
                Assert.Equal(JunctionRiskModel.A_既定.A_跨げない組の較正, l_モデル.A_跨げない組の較正);
            }
        }

        /// <summary>
        /// 誤った本文は、行番号か項目名の入ったメッセージで止める
        /// </summary>
        /// <param name="p_置き換える行">置き換える行の頭 (項目名とタブ)</param>
        /// <param name="p_新しい行">置き換えたあとの行 (空なら行を消す)</param>
        /// <param name="p_メッセージの一部">メッセージに含まれるはずの文字列</param>
        [Theory]
        [InlineData("span\t", "span\t0.6", "\"span\" needs 2 values, found 1")]
        [InlineData("scale\t", "scale\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1", "must be positive")]
        [InlineData("gap\t", "gap\t1\tnan\t1\t1", "is not a finite number")]
        [InlineData("nospan\t", "nospan\t1\t1\t1\t1\nfoo\t1", "unknown item \"foo\"")]
        [InlineData("span\t", "", "missing item \"span\"")]
        [InlineData("name\t", "name\ta\nname\tb", "\"name\" appears more than once")]
        public void Get_解釈_誤った本文は止める(string p_置き換える行, string p_新しい行, string p_メッセージの一部)
        {
            var l_行群 = JunctionRiskModel.A_既定.Get_本文().Split(C_LF).Select(x => x.StartsWith(p_置き換える行, StringComparison.Ordinal) ? p_新しい行 : x);
            var l_例外 = Assert.Throws<InvalidDataException>(() => JunctionRiskModelReader.Get_解釈(string.Join(C_LF, l_行群), C_テストの出所));
            Assert.Contains(p_メッセージの一部, l_例外.Message, StringComparison.Ordinal);
            Assert.StartsWith(C_テストの出所, l_例外.Message, StringComparison.Ordinal);
        }

        /// <summary>
        /// 探す場所にファイルが無ければ null、あれば読んだモデル (出所はそのパス) を返す
        /// </summary>
        [Fact]
        public void Get_読込_あれば_無ければnullであれば読む()
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, Consts.継ぎ目モデルのファイル名);
            Assert.Null(JunctionRiskModelReader.Get_読込_あれば(l_パス));
            File.WriteAllText(l_パス, JunctionRiskModel.A_既定.Get_本文());
            var l_モデル = JunctionRiskModelReader.Get_読込_あれば(l_パス);
            Assert.NotNull(l_モデル);
            Assert.Equal(Path.GetFullPath(l_パス), l_モデル.A_出所);
            Assert.Equal(JunctionRiskModel.A_既定.A_SHA256, l_モデル.A_SHA256);
        }

        /// <summary>
        /// 実行ファイルの隣にファイルが無ければ、-jm を付けないときは組み込みのモデルを使う。-jm で渡したファイルはそのファイルの値を使い、無いファイルは止める
        /// </summary>
        [Fact]
        public void Get_実行時引数_jmで渡したファイルのモデルを使う()
        {
            Assert.False(File.Exists(Path.Combine(AppContext.BaseDirectory, Consts.継ぎ目モデルのファイル名)));
            var l_リード = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            var l_モデルのパス = Path.Combine(this._作業ディレクトリ, Consts.継ぎ目モデルのファイル名);
            File.WriteAllText(l_リード, string.Empty);
            File.WriteAllText(l_モデルのパス, JunctionRiskModel.A_既定.Get_本文());
            Assert.Same(JunctionRiskModel.A_既定, ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_リード]).A_継ぎ目モデル);
            var l_引数 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_リード, Consts.引数キー.継ぎ目モデル, l_モデルのパス]);
            var l_既定 = JunctionRiskModel.A_既定;
            Assert.Equal(Path.GetFullPath(l_モデルのパス), l_引数.A_継ぎ目モデル.A_出所);
            Assert.Equal(l_既定.A_名前, l_引数.A_継ぎ目モデル.A_名前);
            Assert.Equal(l_既定.A_SHA256, l_引数.A_継ぎ目モデル.A_SHA256);
            Assert.Equal(l_既定.A_得点の係数, l_引数.A_継ぎ目モデル.A_得点の係数);
            Assert.Equal(l_既定.A_跨げない組の較正, l_引数.A_継ぎ目モデル.A_跨げない組の較正);
            _ = Assert.ThrowsAny<IOException>(() => ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_リード, Consts.引数キー.継ぎ目モデル, Path.Combine(this._作業ディレクトリ, C_ファイル名_none_tsv)]));
        }

        /// <summary>
        /// 作業ディレクトリを消す
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }

            GC.SuppressFinalize(this);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// テストの実行ファイルの場所からさかのぼって、リポジトリの根 (Tsumiki.csproj のあるディレクトリ) を探す
        /// </summary>
        /// <returns>リポジトリの根</returns>
        private static string Get_リポジトリの根()
        {
            for (var l_場所 = new DirectoryInfo(AppContext.BaseDirectory); l_場所 is not null; l_場所 = l_場所.Parent)
            {
                if (File.Exists(Path.Combine(l_場所.FullName, C_リポジトリの印)))
                {
                    return l_場所.FullName;
                }
            }

            throw new DirectoryNotFoundException(C_リポジトリの印);
        }

        #endregion
    }
}
