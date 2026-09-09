namespace Tsumiki.Model.UnitigBuilding
{
    /// <summary>
    /// ビームサーチの部分経路 1 本ぶんの状態
    /// </summary>
    internal sealed class 先読み探索状態
    {
        /// <summary>
        /// いま到達している頂点
        /// </summary>
        public required int A_現在の頂点 { get; init; }

        /// <summary>
        /// この経路が分岐点から最初に踏んだ頂点 (どの枝を選んだか)
        /// </summary>
        public required int A_最初の1歩 { get; init; }

        /// <summary>
        /// ここまでに積算したペアエンドの支持 (期待本数との比、較正器が
        /// 使えない場合は生カウントそのもの)<br/>
        /// ビームの絞り込み・最終的な
        /// 優勢判定に使う
        /// </summary>
        public required double A_スコア { get; init; }

        /// <summary>
        /// ここまでに積算した生カウント<br/>
        /// 較正器の有無に関わらず、
        /// 「証拠が最低限あるか」の足切り判定にだけ使う
        /// (期待本数が極端に小さい場所では比が実態以上に跳ね上がりうるため)
        /// </summary>
        public required long A_生スコア { get; init; }

        /// <summary>
        /// ここまでに進んだ塩基数<br/>
        /// 先読みの打ち切り判定に使う
        /// </summary>
        public required int A_進んだ長さ { get; init; }

        /// <summary>
        /// unitig ID -> この経路で何回通ったか<br/>
        /// コピー数の予算管理に使う
        /// </summary>
        public required Dictionary<int, int> A_使用回数 { get; init; }
    }
}
