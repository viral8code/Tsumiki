namespace Tsumiki.Model
{
    /// <summary>
    /// ある k 長で1回アセンブリを走らせた結果。
    /// multi-k では k ごとにこれが1つ得られ、比較・マージの対象になる。
    /// </summary>
    internal record アセンブリ実行結果(
        int A_k長,
        string A_ユニティグパス,
        string A_コンティグパス,
        string? A_スキャフォールドパス,
        ulong A_kmerカットオフ,
        double A_単一コピー基準値)
    {
        /// <summary>
        /// このkでの最終成果物。ペアエンドならスキャフォールド、
        /// シングルエンドならコンティグになる。
        /// </summary>
        public string A_最終パス => this.A_スキャフォールドパス ?? this.A_コンティグパス;
    }
}
