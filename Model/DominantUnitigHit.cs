namespace Tsumiki.Model
{
    /// <summary>
    /// ContigMaker.FindDominantUnitig の結果。リードが代表としてマップされた
    /// unitig の ID(符号は向きを表す。正=順鎖、負=逆鎖として一致)に加え、
    /// スキャフォールディングのギャップ長推定に使うためのオフセット情報を持つ。
    ///
    /// Read1 側で使う場合: unitig の「順方向」で見たときに、リードが最後に
    /// ヒットしたk-merの終端位置(unitig内 0-based, 末尾からの残り塩基数を
    /// 求めるのに使う)。
    /// Read2 側で使う場合も同様に、read2 の逆相補鎖を unitig の順方向に
    /// 揃えたうえでの最後のヒット位置を保持する。
    /// これにより「read が unitig の末尾からどれだけ内側で止まっているか」を
    /// 双方について計算し、ギャップ長 = InsertSize - 内側距離1 - 内側距離2 という
    /// 見積りに使える。
    /// </summary>
    internal readonly struct DominantUnitigHit(int unitigId, int hitCount, int lastMatchEndOffset, int unitigLength)
    {
        /// <summary>
        /// マップ先 unitig ID。正=unitig の順鎖として一致、負=逆鎖として一致。
        /// 0はヒットなしを表す。
        /// </summary>
        public readonly int UnitigId = unitigId;

        /// <summary>
        /// 採用された(最多得票の) unitig に対する一致k-mer数。
        /// </summary>
        public readonly int HitCount = hitCount;

        /// <summary>
        /// unitig を UnitigId の符号が示す向きに揃えたときの座標系で、
        /// リードが最後にヒットしたk-merの終端位置(0-based, 末尾側の
        /// インデックス+1。つまりこの値がそのまま「先頭からの既知長」になる)。
        /// </summary>
        public readonly int LastMatchEndOffset = lastMatchEndOffset;

        /// <summary>
        /// マップ先 unitig の全長(向きに依存せず同じ)。
        /// </summary>
        public readonly int UnitigLength = unitigLength;

        public static readonly DominantUnitigHit None = new(0, 0, 0, 0);

        /// <summary>
        /// unitig の末尾から、リードが最後にヒットした位置までの残り塩基数。
        /// この値が小さいほど、リードは unitig の末端近くまで到達している
        /// (＝ペアのもう一方までの未知区間が長くなる可能性が高い)ことを示す。
        /// </summary>
        public int RemainingLength => Math.Max(0, this.UnitigLength - this.LastMatchEndOffset);
    }
}
