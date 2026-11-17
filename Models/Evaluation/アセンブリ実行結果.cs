using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// ある k 長で 1 回アセンブリを走らせた結果
    /// </summary>
    /// <param name="A_k長"></param>
    /// <param name="A_unitigパス"></param>
    /// <param name="A_contigパス"></param>
    /// <param name="A_scaffoldパス"></param>
    /// <param name="A_kmerカットオフ"></param>
    /// <param name="A_単一コピー基準値"></param>
    /// <param name="A_GFAパス"></param>
    /// <param name="A_整合性検査"></param>
    /// <param name="A_実際のコピー数基準"></param>
    /// <param name="A_固定アンカー評価">採用した k・コピー数基準に関わらず、常に同じ固定アンカー k-mer 集合で測った評価<br/> アンカースペクトルが二峰でない等で測れなかった場合は null</param>
    internal record アセンブリ実行結果(int A_k長, string A_unitigパス, string A_contigパス, string? A_scaffoldパス, ulong A_kmerカットオフ, double A_単一コピー基準値, string? A_GFAパス = null, 整合性検査結果? A_整合性検査 = null, コピー数基準の出所 A_実際のコピー数基準 = コピー数基準の出所.Spectrum, アセンブリ評価? A_固定アンカー評価 = null)
    {
        #region カスタムプロパティ

        /// <summary>
        /// ペアエンドなら scaffold、シングルエンドなら contig
        /// </summary>
        public string A_最終パス => this.A_scaffoldパス ?? this.A_contigパス;

        /// <summary>
        /// マルチ k で採用しなかった k の contig の FASTA パス (最終段で、残った N の列を埋めるのに使う)<br/>
        /// 単一 k のときは空
        /// </summary>
        public IReadOnlyList<string> A_ほかのkのcontigパス群 { get; init; } = [];

        /// <summary>
        /// マルチ k で採用しなかった k の scaffold の FASTA パス (最終段で、縮んだ縦の反復を差し替えるのに contig と合わせて使う)<br/>
        /// 単一 k のときと、scaffold を作らなかった k は空
        /// </summary>
        public IReadOnlyList<string> A_ほかのkのscaffoldパス群 { get; init; } = [];

        #endregion
    }
}
