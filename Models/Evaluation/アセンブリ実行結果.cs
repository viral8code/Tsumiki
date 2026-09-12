namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// ある k 長で 1 回アセンブリを走らせた結果
    /// </summary>
    /// <param name="A_k長"></param>
    /// <param name="A_ユニティグパス"></param>
    /// <param name="A_コンティグパス"></param>
    /// <param name="A_スキャフォールドパス"></param>
    /// <param name="A_kmerカットオフ"></param>
    /// <param name="A_単一コピー基準値"></param>
    /// <param name="A_GFAパス"></param>
    /// <param name="A_整合性検査"></param>
    internal record アセンブリ実行結果(int A_k長, string A_ユニティグパス, string A_コンティグパス, string? A_スキャフォールドパス, ulong A_kmerカットオフ, double A_単一コピー基準値, string? A_GFAパス = null, 整合性検査結果? A_整合性検査 = null)
    {
        #region カスタムプロパティ

        /// <summary>
        /// ペアエンドならスキャフォールド、シングルエンドならコンティグ
        /// </summary>
        public string A_最終パス => this.A_スキャフォールドパス ?? this.A_コンティグパス;

        #endregion
    }
}
