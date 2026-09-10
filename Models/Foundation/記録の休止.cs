using Tsumiki.Commons;

namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// Logger の記録を一時的に止めるための解放用ハンドル
    /// </summary>
    internal sealed class 記録の休止 : IDisposable
    {
        #region 公開メソッド

        public 記録の休止()
        {
            lock (Logger._錠)
            {
                Logger._休止の深さ++;
            }
        }

        #endregion

        #region 継承メソッド

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            lock (Logger._錠)
            {
                Logger._休止の深さ--;
            }
        }

        #endregion
    }
}
