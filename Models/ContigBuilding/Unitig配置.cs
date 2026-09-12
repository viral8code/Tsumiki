namespace Tsumiki.Models.ContigBuilding
{
    /// <summary>
    /// V_結合_contig が unitig を結合して contig を作る際、各 unitig が最終的にどの contig の中に、どの向きで、どこに位置したかを表す
    /// </summary>
    /// <remarks>
    /// Scaffolder はこれを使って「unitig 単位のペアエンド隣接候補」を「contig 単位の scaffolding 候補」へ変換する<br/>
    /// contigs.fasta に書き出される配列は、内部的に walk した向き (Forward) そのままの場合と、辞書順で正規化するために逆相補を取った場合 (Reverse) があるため、その正規化情報も保持する
    /// </remarks>
    /// <param name="p_contigID"></param>
    /// <param name="p_IsContig逆相補"></param>
    /// <param name="p_walk順の位置"></param>
    /// <param name="p_walk順の総数"></param>
    /// <param name="p_Iswalk中逆鎖"></param>
    internal readonly struct Unitig配置(int p_contigID, bool p_IsContig逆相補, int p_walk順の位置, int p_walk順の総数, bool p_Iswalk中逆鎖)
    {
        #region 内部変数

        /// <summary>
        /// この unitig が属する contig の ID (FastaWriter に書き出した ID、1 始まり)
        /// </summary>
        public readonly int A_ContigID = p_contigID;

        /// <summary>
        /// contigs.fasta に書き出す際、walk した配列そのものではなくその逆相補を採用した (= 辞書順で正規化した) 場合 true
        /// </summary>
        public readonly bool A_IsContig逆相補 = p_IsContig逆相補;

        /// <summary>
        /// この unitig が contig の walk 順で何番目 (0 始まり) だったか
        /// </summary>
        public readonly int A_walk順の位置 = p_walk順の位置;

        /// <summary>
        /// この unitig が属する contig を構成する unitig の総数
        /// </summary>
        public readonly int A_walk順の総数 = p_walk順の総数;

        /// <summary>
        /// walk の過程でこの unitig が (元の unitigs.fasta 上の向きに対して) 逆鎖として使われた場合 true
        /// </summary>
        public readonly bool A_Iswalk中逆鎖 = p_Iswalk中逆鎖;

        #endregion

        #region プロパティ

        /// <summary>
        /// この unitig が contig の先頭 (5' 端) に位置するか
        /// </summary>
        public bool A_IsContig先頭 => this.A_walk順の位置 == 0;

        /// <summary>
        /// この unitig が contig の末尾 (3' 端) に位置するか
        /// </summary>
        public bool A_IsContig末尾 => this.A_walk順の位置 == this.A_walk順の総数 - 1;

        #endregion
    }
}
