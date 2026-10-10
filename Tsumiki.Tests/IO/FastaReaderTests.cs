using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// FASTA の読み込み (折り返した配列・空の配列・途中で終わるファイル) の確認
    /// </summary>
    public class FastaReaderTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 書式 a desc ACGT  b AAAA CCCC  GG   c TTT
        /// </summary>
        private const string C_書式_a_desc_ACGT__b_AAAA_CCCC__GG___c_TTT = ">a desc\nACGT\n>b\nAAAA\nCCCC\r\nGG\n\n>c\nTTT\n\n\n";

        /// <summary>
        /// 項目 a desc
        /// </summary>
        private const string C_項目_a_desc = "a desc";

        /// <summary>
        /// 項目 b
        /// </summary>
        private const string C_項目_b = "b";

        /// <summary>
        /// 塩基配列 AAAACCCCGG
        /// </summary>
        private const string C_塩基配列_AAAACCCCGG = "AAAACCCCGG";

        /// <summary>
        /// 項目 c
        /// </summary>
        private const string C_項目_c = "c";

        /// <summary>
        /// 塩基配列 TTT
        /// </summary>
        private const string C_塩基配列_TTT = "TTT";

        /// <summary>
        /// 書式 a  b AC  c
        /// </summary>
        private const string C_書式_a__b_AC__c = ">a\n>b\nAC\n>c\n";

        /// <summary>
        /// 項目 a
        /// </summary>
        private const string C_項目_a = "a";

        /// <summary>
        /// 塩基配列 AC
        /// </summary>
        private const string C_塩基配列_AC = "AC";

        /// <summary>
        /// 書式
        /// </summary>
        private const string C_書式 = "\n\n";

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ファイルのパス
        /// </summary>
        private readonly string _パス = Path.GetTempFileName();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 行の配列も、複数行に折り返した配列も、1 本につなげて読む
        /// </summary>
        [Fact]
        public void V_折り返した配列もつなげて読む()
        {
            File.WriteAllText(this._パス, C_書式_a_desc_ACGT__b_AAAA_CCCC__GG___c_TTT);
            var l_全件 = FastaReader.Get_全エントリ(this._パス);
            Assert.Equal([(C_項目_a_desc, Consts.塩基文字), (C_項目_b, C_塩基配列_AAAACCCCGG), (C_項目_c, C_塩基配列_TTT)], l_全件);
        }

        /// <summary>
        /// 配列の無い見出しは空の配列として読む
        /// </summary>
        [Fact]
        public void V_配列の無い見出しは空の配列()
        {
            File.WriteAllText(this._パス, C_書式_a__b_AC__c);
            Assert.Equal([(C_項目_a, string.Empty), (C_項目_b, C_塩基配列_AC), (C_項目_c, string.Empty)], FastaReader.Get_全エントリ(this._パス));
        }

        /// <summary>
        /// 空行しか無いファイルは、止まらずに例外にする
        /// </summary>
        [Fact]
        public void V_空行しか無いファイルは止まらずに例外()
        {
            File.WriteAllText(this._パス, C_書式);
            _ = Assert.Throws<InvalidDataException>(() => FastaReader.Get_全エントリ(this._パス));
        }

        /// <summary>
        /// 一時ファイルを片付ける
        /// </summary>
        public void Dispose()
        {
            File.Delete(this._パス);
            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
