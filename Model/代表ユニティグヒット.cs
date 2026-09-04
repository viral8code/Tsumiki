namespace Tsumiki.Model
{
    /// <summary>
    /// ContigMaker.Get_代表ユニティグ の結果。リードが代表としてマップされた
    /// unitig の ID(符号は向きを表す。正=順鎖、負=逆鎖として一致)に加え、
    /// スキャフォールディングのギャップ長推定に使うためのオフセット情報を持つ。
    ///
    /// Read1 側で使う場合: unitig の「順方向」で見たときに、リードが最後に
    /// ヒットしたk-merの終端位置(unitig内 0-based, 末尾からの残り塩基数を
    /// 求めるのに使う)。
    /// Read2 側で使う場合も同様に、read2 の逆相補鎖を unitig の順方向に
    /// 揃えたうえでの最後のヒット位置を保持する。
    /// これにより「read が unitig の末尾からどれだけ内側で止まっているか」を
    /// 双方について計算し、ギャップ長 = インサートサイズ - 内側距離1 - 内側距離2 という
    /// 見積りに使える。
    /// </summary>
    internal readonly struct 代表ユニティグヒット(int p_ユニティグID, int p_一致kmer数, int p_最終一致終端位置, int p_ユニティグ長)
    {
        /// <summary>
        /// マップ先 unitig ID。正=unitig の順鎖として一致、負=逆鎖として一致。
        /// 0はヒットなしを表す。
        /// </summary>
        public readonly int A_ユニティグID = p_ユニティグID;

        /// <summary>
        /// 採用された(最多得票の) unitig に対する一致k-mer数。
        /// </summary>
        public readonly int A_一致kmer数 = p_一致kmer数;

        /// <summary>
        /// unitig を A_ユニティグID の符号が示す向きに揃えたときの座標系で、
        /// リードが最後にヒットしたk-merの終端位置(0-based, 末尾側の
        /// インデックス+1。つまりこの値がそのまま「先頭からの既知長」になる)。
        /// </summary>
        public readonly int A_最終一致終端位置 = p_最終一致終端位置;

        /// <summary>
        /// マップ先 unitig の全長(向きに依存せず同じ)。
        /// </summary>
        public readonly int A_ユニティグ長 = p_ユニティグ長;

        public static readonly 代表ユニティグヒット A_ヒットなし = new(0, 0, 0, 0);

        /// <summary>
        /// unitig の末尾から、リードが最後にヒットした位置までの残り塩基数。
        /// この値が小さいほど、リードは unitig の末端近くまで到達している
        /// (＝ペアのもう一方までの未知区間が長くなる可能性が高い)ことを示す。
        /// </summary>
        public int A_末尾までの残り長 => Math.Max(0, this.A_ユニティグ長 - this.A_最終一致終端位置);
    }
}
