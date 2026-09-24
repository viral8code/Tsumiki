namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// ビームサーチの部分経路 1 本ぶんの状態
    /// </summary>
    internal sealed class 先読み探索状態
    {
        #region プロパティ

        /// <summary>
        /// いま到達している頂点
        /// </summary>
        public required int A_現在の頂点 { get; init; }

        /// <summary>
        /// この経路が分岐点から最初に踏んだ頂点 (どの枝を選んだか)
        /// </summary>
        public required int A_最初の1歩 { get; init; }

        /// <summary>
        /// ここまでに積算したペアエンドの支持 (期待本数との比、較正器が使えない場合は生カウントそのもの)
        /// </summary>
        public required double A_スコア { get; init; }

        /// <summary>
        /// ここまでに積算した生カウント
        /// </summary>
        public required long A_生スコア { get; init; }

        /// <summary>
        /// ここまでに進んだ塩基数
        /// </summary>
        public required int A_進んだ長さ { get; init; }

        /// <summary>
        /// unitig ID -> この経路で何回通ったか
        /// </summary>
        public required Dictionary<int, int> A_使用回数 { get; init; }

        #endregion
    }
}
