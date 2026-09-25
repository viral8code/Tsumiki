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

        /// <inheritdoc/>
        public override bool CanRead => true;

        /// <inheritdoc/>
        public override bool CanSeek => false;

        /// <inheritdoc/>
        public override bool CanWrite => false;

        /// <inheritdoc/>
        public override long Length => throw new NotSupportedException();

        /// <inheritdoc/>
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }

        /// <inheritdoc/>
        public override int Read(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            return this.Read(p_バッファ.AsSpan(p_開始, p_長さ));
        }

        /// <inheritdoc/>
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

        /// <inheritdoc/>
        public override void Flush()
        {
        }

        /// <inheritdoc/>
        public override long Seek(long p_位置, SeekOrigin p_起点)
        {
            throw new NotSupportedException();
        }

        /// <inheritdoc/>
        public override void SetLength(long p_長さ)
        {
            throw new NotSupportedException();
        }

        /// <inheritdoc/>
        public override void Write(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            throw new NotSupportedException();
        }

        #endregion
    }
}
