using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ を 1 リードずつ読み込む
    /// </summary>
    internal class FastqReader(string p_パス) : SequenceFileReaderBase(p_パス)
    {
        #region 公開メソッド

        /// <summary>
        /// (オーバーライド) 次の 1 行を読み込んで返す
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// FASTQ は 4 行 1 組の固定構造なので、空行に見えても実は EOF という
        /// ケースを区別しないと 4 行の途中で切れたファイルで無限に回り続ける
        /// </remarks>
        protected override string Get_次の行()
        {
            var l_行 = this.Get_次の行_生();
            while (string.IsNullOrWhiteSpace(l_行))
            {
                if (l_行 is null && !this.Get_続きがあるか())
                {
                    throw new InvalidDataException($"{this.A_ファイルパス}: FASTQ が4行の途中で終わっている。");
                }
                l_行 = this.Get_次の行_生();
            }
            return l_行;
        }

        /// <summary>
        /// 指定したファイル群のリードを、塩基列だけを取り出して順に流す
        /// </summary>
        /// <param name="p_パス群"></param>
        /// <remarks>
        /// 最終成果物へリードを貼り直す処理 (ポリッシュ・閉じ目の検証) のように、
        /// ID もクオリティも要らない全走査のための入口
        /// </remarks>
        public static IEnumerable<string> Get_生リード列(params string?[] p_パス群)
        {
            foreach (var l_パス in p_パス群)
            {
                if (string.IsNullOrWhiteSpace(l_パス))
                {
                    continue;
                }
                using var l_読み込み = new FastqReader(l_パス);
                while (l_読み込み.Get_続きがあるか())
                {
                    yield return l_読み込み.Get_次のレコード().A_配列;
                }
            }
        }

        /// <summary>
        /// 次の 1 リードを読み込んで返す
        /// </summary>
        /// <returns>読み込んだリード</returns>
        public リードデータ Get_次のリード()
        {
            try
            {
                var (l_ID, l_配列, l_クオリティ) = this.Get_次のレコード();
                return new リードデータ()
                {
                    A_ID = l_ID,
                    A_塩基候補列 = Util.V_変換_塩基候補列(l_配列),
                    A_生リード = l_配列,
                    A_クオリティ = l_クオリティ,
                };
            }
            catch (Exception ex)
            {
                Logger.V_出力_警告(Logger.Get_メソッド名(), ex);
                throw;
            }
        }

        /// <summary>
        /// 曖昧塩基を無視する経路向けの軽量版
        /// </summary>
        /// <remarks>
        /// A_塩基候補列(List&lt;byte[]&gt;)の
        /// 代わりに A_塩基列(byte[])のみを構築する<br/>
        /// KmerCounting.V_読込_リードファイル から使用する
        /// </remarks>
        public リードデータ Get_次のリード_軽量()
        {
            try
            {
                var (l_ID, l_配列, l_クオリティ) = this.Get_次のレコード();
                return new リードデータ()
                {
                    A_ID = l_ID,
                    A_塩基列 = Util.V_変換_塩基列(l_配列),
                    A_生リード = l_配列,
                    A_クオリティ = l_クオリティ,
                };
            }
            catch (Exception ex)
            {
                Logger.V_出力_警告(Logger.Get_メソッド名(), ex);
                throw;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列とクオリティの長さが合わない FASTQ は、そのまま進めると
        /// 品質判定が配列の範囲外を触って落ちる
        /// </summary>
        /// <param name="p_ID"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_クオリティ"></param>
        /// <remarks>
        /// どのリードが不正かを言って止める
        /// </remarks>
        private void V_検査(string p_ID, string p_配列, string p_クオリティ)
        {
            if (p_配列.Length != p_クオリティ.Length)
            {
                throw new InvalidDataException($"{this.A_ファイルパス}: リード {p_ID} の塩基列({p_配列.Length}文字)とクオリティ({p_クオリティ.Length}文字)の長さが一致しない。");
            }
        }

        /// <summary>
        /// ID・配列・クオリティの 1 レコードを読み込んで返す
        /// </summary>
        /// <returns></returns>
        private (string A_ID, string A_配列, string A_クオリティ) Get_次のレコード()
        {
            var l_ID = this.Get_次の行();
            var l_配列 = this.Get_次の行();
            _ = this.Get_次の行();
            var l_クオリティ = this.Get_次の行();
            this.V_検査(l_ID, l_配列, l_クオリティ);
            return (l_ID, l_配列, l_クオリティ);
        }

        #endregion
    }
}
