namespace Tsumiki.Model.Evaluation
{
    /// <summary>
    /// ある k 長で1回アセンブリを走らせた結果。
    /// </summary>
    internal record アセンブリ実行結果(
        int A_k長,
        string A_ユニティグパス,
        string A_コンティグパス,
        string? A_スキャフォールドパス,
        ulong A_kmerカットオフ,
        double A_単一コピー基準値,
        string? A_GFAパス = null,
        整合性検査結果? A_整合性検査 = null)
    {
        /// <summary>
        /// ペアエンドならスキャフォールド、シングルエンドならコンティグ。
        /// </summary>
        public string A_最終パス => this.A_スキャフォールドパス ?? this.A_コンティグパス;
    }
}
