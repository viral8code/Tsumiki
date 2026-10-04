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

        /// <summary>
        /// 位置と置き換える塩基ごとの、信頼できる窓の純増数の控え (確かな値か、その値以下と分かっている上界)
        /// </summary>
        public int[] A_改善数の控え { get; private set; } = [];

        /// <summary>
        /// 改善数の控えの種類 (0 = 無し、1 = 確かな値、2 = 上界)
        /// </summary>
        public byte[] A_控えの種類 { get; private set; } = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 窓数と塩基数ぶんの長さを確保する
        /// </summary>
        /// <param name="p_窓数">窓の数</param>
        /// <param name="p_塩基数">リードの塩基数</param>
        public void V_確保(int p_窓数, int p_塩基数)
        {
            if (this.A_パック.Length < p_窓数)
            {
                this.A_パック = new UInt128[p_窓数];
                this.A_逆相補 = new UInt128[p_窓数];
                this.A_無効数 = new int[p_窓数];
                this.A_信頼状況 = new bool[p_窓数];
                this.A_未信頼累積 = new int[p_窓数 + 1];
            }

            if (this.A_改善数の控え.Length < 4 * p_塩基数)
            {
                this.A_改善数の控え = new int[4 * p_塩基数];
                this.A_控えの種類 = new byte[4 * p_塩基数];
            }
        }

        #endregion
    }
}
