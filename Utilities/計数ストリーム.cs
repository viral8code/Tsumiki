namespace Tsumiki.Utilities
{
    /// <summary>
    /// 展開後の読み書きの量を数え、長さと位置を答える
    /// </summary>
    internal sealed class 計数ストリーム : Stream
    {
        #region 内部変数

        /// <summary>
        /// 包んでいるストリーム
        /// </summary>
        private readonly Stream _中身;

        /// <summary>
        /// 書き終えたときに展開後の長さを渡す先 (読み込み用なら null)
        /// </summary>
        private readonly Action<long>? _確定;

        /// <summary>
        /// 展開後の長さ (読み込み用のときだけ分かる)
        /// </summary>
        private readonly long _長さ;

        /// <summary>
        /// ここまでに読み書きした展開後の量
        /// </summary>
        private long _位置;

        /// <summary>
        /// 閉じたか
        /// </summary>
        private bool _Is確定済み;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 書き込み用
        /// </summary>
        /// <param name="p_中身">包む圧縮ストリーム</param>
        /// <param name="p_確定">閉じたときに展開後の長さを渡す先</param>
        public 計数ストリーム(Stream p_中身, Action<long> p_確定)
        {
            this._中身 = p_中身;
            this._確定 = p_確定;
        }

        /// <summary>
        /// 読み込み用
        /// </summary>
        /// <param name="p_中身">包む展開ストリーム</param>
        /// <param name="p_長さ">展開後の長さ</param>
        public 計数ストリーム(Stream p_中身, long p_長さ)
        {
            this._中身 = p_中身;
            this._長さ = p_長さ;
        }

        #endregion

        #region 継承メソッド

        /// <summary>
        /// (オーバーライド) 読み込み可能かを返す
        /// </summary>
        public override bool CanRead => this._確定 is null;

        /// <summary>
        /// (オーバーライド) 位置を移動できるかを返す
        /// </summary>
        public override bool CanSeek => false;

        /// <summary>
        /// (オーバーライド) 書き込み可能かを返す
        /// </summary>
        public override bool CanWrite => this._確定 is not null;

        /// <summary>
        /// (オーバーライド) ストリームの長さを返す
        /// </summary>
        public override long Length => this._確定 is null ? this._長さ : this._位置;

        /// <summary>
        /// (オーバーライド) 現在位置を返す
        /// </summary>
        public override long Position { get => this._位置; set => throw new NotSupportedException(); }

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
            var l_読んだ = this._中身.Read(p_バッファ);
            this._位置 += l_読んだ;
            return l_読んだ;
        }

        /// <summary>
        /// (オーバーライド) バッファを書き込む
        /// </summary>
        /// <param name="p_バッファ"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_長さ"></param>
        public override void Write(byte[] p_バッファ, int p_開始, int p_長さ)
        {
            this.Write(p_バッファ.AsSpan(p_開始, p_長さ));
        }

        /// <summary>
        /// (オーバーライド) バッファを書き込む
        /// </summary>
        /// <param name="p_バッファ"></param>
        public override void Write(ReadOnlySpan<byte> p_バッファ)
        {
            this._中身.Write(p_バッファ);
            this._位置 += p_バッファ.Length;
        }

        /// <summary>
        /// (オーバーライド) 書き込み内容を確定する
        /// </summary>
        public override void Flush()
        {
            this._中身.Flush();
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
        /// (オーバーライド) 包んでいるストリームを閉じてから、展開後の長さを渡す
        /// </summary>
        /// <param name="p_Is明示">Dispose から呼ばれたか</param>
        protected override void Dispose(bool p_Is明示)
        {
            if (!this._Is確定済み)
            {
                this._Is確定済み = true;
                this._中身.Dispose();
                this._確定?.Invoke(this._位置);
            }

            base.Dispose(p_Is明示);
        }

        #endregion
    }
}
