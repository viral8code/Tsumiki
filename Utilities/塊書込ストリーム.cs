namespace Tsumiki.Utilities
{
    /// <summary>
    /// 塊の列へ追記する
    /// </summary>
    internal sealed class 塊書込ストリーム : Stream
    {
        #region 内部変数

        /// <summary>
        /// 書きかけの塊
        /// </summary>
        private byte[] _今の塊 = new byte[中間データ置き場.C_塊の大きさ];

        /// <summary>
        /// 書きかけの塊の使用量
        /// </summary>
        private int _今の位置;

        /// <summary>
        /// 閉じたか
        /// </summary>
        private bool _Is確定済み;

        #endregion

        #region プロパティ

        /// <summary>
        /// 書き終えた塊 (閉じた後は書きかけの塊も含む)
        /// </summary>
        public List<byte[]> A_塊群 { get; } = [];

        #endregion

        #region 継承メソッド

        /// <inheritdoc/>
        public override bool CanRead => false;

        /// <inheritdoc/>
        public override bool CanSeek => false;

        /// <inheritdoc/>
        public override bool CanWrite => true;

        /// <inheritdoc/>
        public override long Length => throw new NotSupportedException();

        /// <inheritdoc/>
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }

        /// <inheritdoc/>
        public override void Write(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            this.Write(p_バッファ.AsSpan(p_開始, p_長さ));
        }

        /// <inheritdoc/>
        public override void Write(ReadOnlySpan<byte> p_バッファ)
        {
            while (!p_バッファ.IsEmpty)
            {
                if (this._今の位置 == this._今の塊.Length)
                {
                    this.A_塊群.Add(this._今の塊);
                    this._今の塊 = new byte[中間データ置き場.C_塊の大きさ];
                    this._今の位置 = 0;
                }

                var l_書く分 = Math.Min(p_バッファ.Length, this._今の塊.Length - this._今の位置);
                p_バッファ[..l_書く分].CopyTo(this._今の塊.AsSpan(this._今の位置));
                this._今の位置 += l_書く分;
                p_バッファ = p_バッファ[l_書く分..];
            }
        }

        /// <inheritdoc/>
        public override void Flush()
        {
        }

        /// <inheritdoc/>
        public override int Read(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            throw new NotSupportedException();
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

        /// <summary>
        /// 書きかけの塊を、使った分だけに詰めて塊の列へ加える
        /// </summary>
        /// <param name="p_Is明示">Dispose から呼ばれたか</param>
        protected override void Dispose(bool p_Is明示)
        {
            if (!this._Is確定済み)
            {
                this._Is確定済み = true;
                if (this._今の位置 > 0)
                {
                    this.A_塊群.Add(this._今の塊.AsSpan(0, this._今の位置).ToArray());
                }
            }

            base.Dispose(p_Is明示);
        }

        #endregion
    }
}
