namespace Tsumiki.Model
{
    /// <summary>
    /// ビームサーチの部分経路1本ぶんの状態。
    /// </summary>
    internal sealed class 先読み探索状態
    {
        /// <summary>いま到達している頂点。</summary>
        public required int A_現在の頂点 { get; init; }

        /// <summary>この経路が分岐点から最初に踏んだ頂点(どの枝を選んだか)。</summary>
        public required int A_最初の1歩 { get; init; }

        /// <summary>ここまでに積算したペアエンドの支持。</summary>
        public required long A_スコア { get; init; }

        /// <summary>ここまでに進んだ塩基数。先読みの打ち切り判定に使う。</summary>
        public required int A_進んだ長さ { get; init; }

        /// <summary>unitig ID -> この経路で何回通ったか。コピー数の予算管理に使う。</summary>
        public required Dictionary<int, int> A_使用回数 { get; init; }
    }
}
