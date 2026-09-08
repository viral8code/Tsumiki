namespace Tsumiki.Model
{
    /// <summary>
    /// V_結合_コンティグ が unitig を結合して contig を作る際、各 unitig が
    /// 最終的にどの contig の中に、どの向きで、どこに位置したかを表す。
    /// Scaffolder はこれを使って「unitig 単位のペアエンド隣接候補」を
    /// 「contig 単位のスキャフォールディング候補」へ変換する。
    ///
    /// contigs.fasta に書き出される配列は、内部的に walk した向き
    /// (Forward)そのままの場合と、辞書順で正規化するために逆相補を
    /// 取った場合(Reverse)があるため、その正規化情報も保持する。
    /// </summary>
    internal readonly struct ユニティグ配置(
        int p_コンティグID,
        bool p_コンティグが逆相補か,
        int p_walk順の位置,
        int p_walk順の総数,
        bool p_walk中で逆鎖か)
    {
        /// <summary>
        /// この unitig が属する contig の ID(FastaWriter に書き出した ID、1始まり)。
        /// </summary>
        public readonly int A_コンティグID = p_コンティグID;

        /// <summary>
        /// contigs.fasta に書き出す際、walk した配列そのものではなく
        /// その逆相補を採用した(=辞書順で正規化した)場合 true。
        /// </summary>
        public readonly bool A_コンティグが逆相補か = p_コンティグが逆相補か;

        /// <summary>
        /// この unitig が contig の walk 順で何番目(0始まり)だったか。
        /// </summary>
        public readonly int A_walk順の位置 = p_walk順の位置;

        /// <summary>
        /// この unitig が属する contig を構成する unitig の総数。
        /// </summary>
        public readonly int A_walk順の総数 = p_walk順の総数;

        /// <summary>
        /// walk の過程でこの unitig が(元の unitigs.fasta 上の向きに対して)
        /// 逆鎖として使われた場合 true。
        /// </summary>
        public readonly bool A_walk中で逆鎖か = p_walk中で逆鎖か;

        /// <summary>
        /// この unitig が contig の先頭(5'端)に位置するか。
        /// </summary>
        public bool A_コンティグ先頭か => this.A_walk順の位置 == 0;

        /// <summary>
        /// この unitig が contig の末尾(3'端)に位置するか。
        /// </summary>
        public bool A_コンティグ末尾か => this.A_walk順の位置 == this.A_walk順の総数 - 1;
    }
}
