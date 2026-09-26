namespace Tsumiki.Utilities
{
    /// <summary>
    /// 塊の列を先頭から順に読む
    /// </summary>
    /// <param name="p_塊群">読む塊の列</param>
    internal sealed class 塊読込ストリーム(List<byte[]> p_塊群) : Stream
    {
        #region 内部変数

        /// <summary>
        /// 読んでいる塊の番号
        /// </summary>
        private int _塊番号;

        /// <summary>
        /// 読んでいる塊の中の位置
        /// </summary>
        private int _位置;

        #endregion

        #region 継承メソッド

        /// <summary>
        /// (オーバーライド) 読み込み可能かを返す
        /// </summary>
        public override bool CanRead => true;

        /// <summary>
        /// (オーバーライド) 位置を移動できるかを返す
        /// </summary>
        public override bool CanSeek => false;

        /// <summary>
        /// (オーバーライド) 書き込み可能かを返す
        /// </summary>
        public override bool CanWrite => false;

        /// <summary>
        /// (オーバーライド) ストリームの長さを返す
        /// </summary>
        public override long Length => throw new NotSupportedException();

        /// <summary>
        /// (オーバーライド) 現在位置を返す
        /// </summary>
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }

        /// <summary>
        /// (オーバーライド) バッファへ読み込む
        /// </summary>
        /// <param name="p_バッファ"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        public override int Read(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            return this.Read(p_バッファ.AsSpan(p_開始, p_長さ));
        }

        /// <summary>
        /// (オーバーライド) バッファへ読み込む
        /// </summary>
        /// <param name="p_バッファ"></param>
        /// <returns></returns>
        public override int Read(Span<byte> p_バッファ)
        {
            var l_読んだ = 0;
            while (!p_バッファ.IsEmpty && this._塊番号 < p_塊群.Count)
            {
                var l_塊 = p_塊群[this._塊番号];
                if (this._位置 == l_塊.Length)
                {
                    this._塊番号++;
                    this._位置 = 0;
                    continue;
                }

                var l_読む分 = Math.Min(p_バッファ.Length, l_塊.Length - this._位置);
                l_塊.AsSpan(this._位置, l_読む分).CopyTo(p_バッファ);
                this._位置 += l_読む分;
                l_読んだ += l_読む分;
                p_バッファ = p_バッファ[l_読む分..];
            }

            return l_読んだ;
        }

        /// <summary>
        /// (オーバーライド) 書き込み内容を確定する
        /// </summary>
        public override void Flush()
        {
        }

        /// <summary>
        /// (オーバーライド) 読み書き位置を移動する
        /// </summary>
        /// <param name="p_位置"></param>
        /// <param name="p_起点"></param>
        /// <returns></returns>
        public override long Seek(long p_位置, SeekOrigin p_起点)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        /// (オーバーライド) ストリームの長さを変更する
        /// </summary>
        /// <param name="p_長さ"></param>
        public override void SetLength(long p_長さ)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        /// (オーバーライド) バッファを書き込む
        /// </summary>
        /// <param name="p_バッファ"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_長さ"></param>
        public override void Write(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            throw new NotSupportedException();
        }

        #endregion
    }
}
