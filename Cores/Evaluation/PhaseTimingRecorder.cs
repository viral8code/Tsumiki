using Tsumiki.Models.Reporting;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// StageTimer が計測した工程ごとの資源使用量を集めておく収集器
    /// </summary>
    internal static class PhaseTimingRecorder
    {
        #region 内部変数

        /// <summary>
        /// 錠
        /// </summary>
        private static readonly Lock _錠 = new();

        /// <summary>
        /// 記録した計測の一覧
        /// </summary>
        private static readonly List<フェーズ計測> _記録 = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 工程ぶんの計測を書き留める
        /// </summary>
        /// <param name="p_計測"></param>
        public static void V_記録(フェーズ計測 p_計測)
        {
            lock (_錠)
            {
                _記録.Add(p_計測);
            }
        }

        /// <summary>
        /// これまでに書き留めた計測の一覧
        /// </summary>
        /// <returns></returns>
        public static IReadOnlyList<フェーズ計測> Get_記録()
        {
            lock (_錠)
            {
                return [.. _記録];
            }
        }

        #endregion
    }
}
