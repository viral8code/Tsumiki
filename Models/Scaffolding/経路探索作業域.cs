namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// パック値で状態を持つ経路探索の作業領域を、スレッドごとに使い回して持つ
    /// </summary>
    /// <remarks>
    /// リードペアの橋渡しは探索をペアの数だけ繰り返すため、呼び出しごとに確保すると確保と回収が支配的になる
    /// </remarks>
    internal sealed class 経路探索作業域
    {
        #region 定数

        /// <summary>
        /// 使い回さずに作り直す状態数
        /// </summary>
        /// <remarks>
        /// 辞書の消去は容量に比例するため、大きく育った作業域を消し続けると小さな探索まで遅くなる
        /// </remarks>
        public const int 作り直す状態数 = 8_192;

        #endregion

        #region プロパティ

        /// <summary>
        /// 状態ごとの親と、親から伸ばした 1 塩基
        /// </summary>
        public List<(int A_親, byte A_塩基)> A_節点 { get; } = new(1_024);

        /// <summary>
        /// 状態ごとの順鎖・逆鎖のパック値
        /// </summary>
        public List<(UInt128 A_順上, UInt128 A_順下, UInt128 A_逆上, UInt128 A_逆下)> A_状態群 { get; } = new(1_024);

        /// <summary>
        /// 状態ごとの継ぎ足した塩基数
        /// </summary>
        public List<int> A_深さ群 { get; } = new(1_024);

        /// <summary>
        /// 順鎖のパック値と深さから、最初に着いた状態を引く
        /// </summary>
        public Dictionary<(UInt128, UInt128, int), int> A_到達済み { get; } = new(1_024);

        /// <summary>
        /// 状態ごとの、別経路からも到達されたかの印
        /// </summary>
        public List<bool> A_多重到達 { get; } = new(1_024);

        /// <summary>
        /// 展開を待つ状態
        /// </summary>
        public Queue<int> A_キュー { get; } = new();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 次の探索に備えて空にする
        /// </summary>
        public void V_初期化()
        {
            this.A_節点.Clear();
            this.A_状態群.Clear();
            this.A_深さ群.Clear();
            this.A_到達済み.Clear();
            this.A_多重到達.Clear();
            this.A_キュー.Clear();
        }

        #endregion
    }
}
