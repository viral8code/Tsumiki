using System.Buffers.Binary;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer と出現回数の並んだファイルを、確保をせずに 1 件ずつ読む
    /// </summary>
    internal sealed class エントリ読み込み : IDisposable
    {
        #region 定数

        /// <summary>
        /// 1 回に読み込む件数
        /// </summary>
        private const int C_読み込む件数 = 4_096;

        #endregion

        #region 内部変数

        /// <summary>
        /// 読み込み元
        /// </summary>
        private readonly Stream _流れ;

        /// <summary>
        /// パック済みの k-mer の長さ
        /// </summary>
        private readonly int _パック長;

        /// <summary>
        /// 1 件の長さ
        /// </summary>
        private readonly int _エントリ長;

        /// <summary>
        /// 読み込んだ件の置き場
        /// </summary>
        private readonly byte[] _バッファ;

        /// <summary>
        /// バッファ内の有効な長さ
        /// </summary>
        private int _有効長;

        /// <summary>
        /// 今の件の位置
        /// </summary>
        private int _位置;

        #endregion

        #region プロパティ

        /// <summary>
        /// 今の件があるか
        /// </summary>
        public bool Has項目 => this._位置 + this._エントリ長 <= this._有効長;

        /// <summary>
        /// 今の件の k-mer
        /// </summary>
        public ReadOnlySpan<byte> A_キー => this._バッファ.AsSpan(this._位置, this._パック長);

        /// <summary>
        /// 今の件の出現回数
        /// </summary>
        public ulong A_出現回数 => BinaryPrimitives.ReadUInt64LittleEndian(this._バッファ.AsSpan(this._位置 + this._パック長, sizeof(ulong)));

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 読み込み元を開いて最初の件を読む
        /// </summary>
        /// <param name="p_流れ">読み込み元</param>
        /// <param name="p_パック長">パック済みの k-mer の長さ</param>
        public エントリ読み込み(Stream p_流れ, int p_パック長)
        {
            this._流れ = p_流れ;
            this._パック長 = p_パック長;
            this._エントリ長 = p_パック長 + sizeof(ulong);
            this._バッファ = new byte[this._エントリ長 * C_読み込む件数];
            this.V_補充();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 次の件へ進む
        /// </summary>
        public void V_進む()
        {
            this._位置 += this._エントリ長;
            if (!this.Has項目)
            {
                this.V_補充();
            }
        }

        /// <summary>
        /// 読み込み元を閉じる
        /// </summary>
        public void Dispose()
        {
            this._流れ.Dispose();
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 読み残しを先頭へ寄せ、空いた分を読み込む
        /// </summary>
        private void V_補充()
        {
            var l_残り = this._有効長 - this._位置;
            if (l_残り > 0)
            {
                Array.Copy(this._バッファ, this._位置, this._バッファ, 0, l_残り);
            }

            this._位置 = 0;
            this._有効長 = l_残り + this._流れ.ReadAtLeast(this._バッファ.AsSpan(l_残り), this._バッファ.Length - l_残り, throwOnEndOfStream: false);
        }

        #endregion
    }
}
