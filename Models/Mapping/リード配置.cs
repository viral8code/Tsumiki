namespace Tsumiki.Models.Mapping
{
    /// <summary>
    /// リードを参照配列へ配置した結果
    /// </summary>
    internal readonly record struct リード配置(int A_配列番号, bool A_逆鎖か, int A_スコア, int A_信頼度, IReadOnlyList<整列位置> A_整列位置群)
    {
        #region 定数

        /// <summary>
        /// 配置できなかった結果
        /// </summary>
        public static readonly リード配置 配置なし = new(-1, false, 0, 0, []);

        #endregion
    }
}
