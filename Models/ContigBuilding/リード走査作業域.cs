namespace Tsumiki.Models.ContigBuilding
{
    /// <summary>
    /// リード 1 本を unitig の k-mer 索引で走査した結果を、ワーカーごとに使い回して持つ
    /// </summary>
    internal sealed class リード走査作業域
    {
        #region プロパティ

        /// <summary>
        /// 区間ごとの符号付き unitig ID (隣り合う区間は必ず異なる)
        /// </summary>
        public List<int> A_経路 { get; } = [];

        /// <summary>
        /// 区間ごとのヒット k-mer 数
        /// </summary>
        public List<int> A_票数 { get; } = [];

        /// <summary>
        /// 区間ごとの、unitig の向きに揃えた最後のヒットの終端位置
        /// </summary>
        public List<int> A_終端位置 { get; } = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 次のリードの走査に備えて空にする
        /// </summary>
        public void V_初期化()
        {
            this.A_経路.Clear();
            this.A_票数.Clear();
            this.A_終端位置.Clear();
        }

        #endregion
    }
}
