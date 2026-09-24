using System.Globalization;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Commons
{
    /// <summary>
    /// メッセージ ID から表示用の文言を作る
    /// </summary>
    internal static class Messages
    {
        #region プロパティ

        /// <summary>
        /// 表示に使う言語
        /// </summary>
        public static 言語 A_言語 { get; set; } = 言語.日本語;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_ID の文言に p_引数 を差し込んで返す
        /// </summary>
        /// <param name="p_ID"></param>
        /// <param name="p_引数"></param>
        /// <returns></returns>
        public static string Get_文言(メッセージID p_ID, params object?[] p_引数)
        {
            var l_書式 = MessageCatalog.Get_書式(A_言語, p_ID);
            return p_引数.Length == 0 ? l_書式 : string.Format(CultureInfo.InvariantCulture, l_書式, p_引数);
        }

        #endregion
    }
}
