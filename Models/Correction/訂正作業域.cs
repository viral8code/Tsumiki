namespace Tsumiki.Models.Correction
{
    /// <summary>
    /// 1 本の訂正で使う配列を、リードをまたいで使い回す置き場
    /// </summary>
    internal sealed class 訂正作業域
    {
        #region プロパティ

        /// <summary>
        /// 窓ごとの順鎖のパック値
        /// </summary>
        public UInt128[] A_パック { get; private set; } = [];

        /// <summary>
        /// 窓ごとの逆相補のパック値
        /// </summary>
        public UInt128[] A_逆相補 { get; private set; } = [];

        /// <summary>
        /// 窓ごとの無効な塩基の数
        /// </summary>
        public int[] A_無効数 { get; private set; } = [];

        /// <summary>
        /// 窓ごとに信頼できるか
        /// </summary>
        public bool[] A_信頼状況 { get; private set; } = [];

        /// <summary>
        /// 信頼できない窓の数の累積
        /// </summary>
        public int[] A_未信頼累積 { get; private set; } = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 窓数ぶんの長さを確保する
        /// </summary>
        /// <param name="p_窓数">窓の数</param>
        public void V_確保(int p_窓数)
        {
            if (this.A_パック.Length >= p_窓数)
            {
                return;
            }

            this.A_パック = new UInt128[p_窓数];
            this.A_逆相補 = new UInt128[p_窓数];
            this.A_無効数 = new int[p_窓数];
            this.A_信頼状況 = new bool[p_窓数];
            this.A_未信頼累積 = new int[p_窓数 + 1];
        }

        #endregion
    }
}
