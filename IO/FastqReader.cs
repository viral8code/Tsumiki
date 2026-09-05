using System.IO.Compression;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.IO
{
    internal class FastqReader : IDisposable
    {
        public string A_ファイルパス { get; private set; }

        private readonly StreamReader _読み込み;

        private const int バッファサイズ = 1 << 25;

        public FastqReader(string p_パス)
        {
            this.A_ファイルパス = p_パス;
            var l_入力ストリーム = new FileStream(p_パス, FileMode.Open, FileAccess.Read);
            if (Path.GetExtension(p_パス)?.ToLower() == ".gz")
            {
                var l_展開ストリーム = new GZipStream(l_入力ストリーム, CompressionMode.Decompress);
                this._読み込み = new(l_展開ストリーム, bufferSize: バッファサイズ);
            }
            else
            {
                this._読み込み = new(l_入力ストリーム, bufferSize: バッファサイズ);
            }
        }

        public bool Get_続きがあるか()
        {
            return !this._読み込み.EndOfStream;
        }

        private string Get_次の行()
        {
            var l_行 = this._読み込み.ReadLine();
            while (string.IsNullOrWhiteSpace(l_行))
            {
                // ReadLine は EOF でも null を返す。空行の読み飛ばしと区別しないと
                // 途中で切れたファイルで無限に回り続ける。
                if (l_行 is null && this._読み込み.EndOfStream)
                {
                    throw new InvalidDataException(
                        $"{this.A_ファイルパス}: FASTQ が4行の途中で終わっている。");
                }
                l_行 = this._読み込み.ReadLine();
            }
            return l_行;
        }

        /// <summary>
        /// 配列とクオリティの長さが合わない FASTQ は、そのまま進めると
        /// 品質判定が配列の範囲外を触って落ちる。どのリードが不正かを言って止める。
        /// </summary>
        private void V_検査(string p_ID, string p_配列, string p_クオリティ)
        {
            if (p_配列.Length != p_クオリティ.Length)
            {
                throw new InvalidDataException(
                    $"{this.A_ファイルパス}: リード {p_ID} の塩基列({p_配列.Length}文字)と" +
                    $"クオリティ({p_クオリティ.Length}文字)の長さが一致しない。");
            }
        }

        public リードデータ Get_次のリード()
        {
            try
            {
                var l_ID = this.Get_次の行();
                var l_配列 = this.Get_次の行();
                _ = this.Get_次の行();
                var l_クオリティ = this.Get_次の行();
                this.V_検査(l_ID, l_配列, l_クオリティ);

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
        /// 曖昧塩基を無視する経路向けの軽量版。A_塩基候補列(List&lt;byte[]&gt;)の
        /// 代わりに A_塩基列(byte[])のみを構築する。
        /// KmerCounting.V_読込_リードファイル から使用する。
        /// </summary>
        public リードデータ Get_次のリード_軽量()
        {
            try
            {
                var l_ID = this.Get_次の行();
                var l_配列 = this.Get_次の行();
                _ = this.Get_次の行();
                var l_クオリティ = this.Get_次の行();
                this.V_検査(l_ID, l_配列, l_クオリティ);

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

        public void Dispose()
        {
            this._読み込み.Dispose();
        }
    }
}
