namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// ContigMaker.Get_代表ユニティグ の結果
    /// </summary>
    /// <remarks>
    /// リードが代表としてマップされた unitig の ID (符号は向きを表す<br/>
    /// 正=順鎖、負=逆鎖として一致) に加え、スキャフォールディングのギャップ長推定に使うためのオフセット情報を持つ<br/>
    /// Read1 側で使う場合: unitig の「順方向」で見たときに、リードが最後にヒットした k-mer の終端位置 (unitig 内 0-based, 末尾からの残り塩基数を求めるのに使う) <br/>
    /// Read2 側で使う場合も同様に、read2 の逆相補鎖を unitig の順方向に揃えたうえでの最後のヒット位置を保持する<br/>
    /// これにより「read が unitig の末尾からどれだけ内側で止まっているか」を双方について計算し、ギャップ長 = インサートサイズ - 内側距離 1 - 内側距離 2 という見積りに使える
    /// </remarks>
    /// <param name="p_ユニティグID"></param>
    /// <param name="p_一致kmer数"></param>
    /// <param name="p_最終一致終端位置"></param>
    /// <param name="p_ユニティグ長"></param>
    internal readonly struct 代表ユニティグヒット(int p_ユニティグID, int p_一致kmer数, int p_最終一致終端位置, int p_ユニティグ長)
    {
        #region 定数

        /// <summary>
        /// ヒットなし
        /// </summary>
        public static readonly 代表ユニティグヒット A_ヒットなし = new(0, 0, 0, 0);

        #endregion

        #region 内部変数

        /// <summary>
        /// マップ先 unitig ID
        /// </summary>
        /// <remarks>
        /// 正=unitig の順鎖として一致、負=逆鎖として一致<br/>
        /// 0 はヒットなしを表す
        /// </remarks>
        public readonly int A_ユニティグID = p_ユニティグID;

        /// <summary>
        /// 採用された (最多得票の) unitig に対する一致 k-mer 数
        /// </summary>
        public readonly int A_一致kmer数 = p_一致kmer数;

        /// <summary>
        /// ユニティグの向きに揃えた最後の一致 k-mer の終端位置 (末尾の添字 + 1)
        /// </summary>
        /// <remarks>
        /// この値をユニティグ先頭からの既知長として使う
        /// </remarks>
        public readonly int A_最終一致終端位置 = p_最終一致終端位置;

        /// <summary>
        /// マップ先 unitig の全長 (向きに依存せず同じ)
        /// </summary>
        public readonly int A_ユニティグ長 = p_ユニティグ長;

        #endregion

        #region プロパティ

        /// <summary>
        /// unitig の末尾から、リードが最後にヒットした位置までの残り塩基数
        /// </summary>
        /// <remarks>
        /// この値が小さいほど、リードは unitig の末端近くまで到達している (＝ペアのもう一方までの未知区間が長くなる可能性が高い) ことを示す
        /// </remarks>
        public int A_末尾までの残り長 => Math.Max(0, this.A_ユニティグ長 - this.A_最終一致終端位置);

        #endregion
    }
}
