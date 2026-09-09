using System.Globalization;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Common
{
    /// <summary>
    /// メッセージID から表示用の文言を作る<br/>
    /// 文言を呼び出し側に直接書かないのは、出す場所と訳す場所を分けるため<br/>
    /// 書式は <see cref="MessageCatalog"/> が言語ごとに持ち、ここは
    /// 現在の言語を選んで差し込むだけにする
    /// </summary>
    internal static class Messages
    {
        /// <summary>
        /// 表示に使う言語
        /// </summary>
        public static 言語 A_言語 { get; set; } = 言語.日本語;

        /// <summary>
        /// p_ID の文言に p_引数 を差し込んで返す
        /// </summary>
        public static string Get_文言(メッセージID p_ID, params object?[] p_引数)
        {
            var l_書式 = MessageCatalog.Get_書式(A_言語, p_ID);
            return p_引数.Length == 0
                ? l_書式
                : string.Format(CultureInfo.InvariantCulture, l_書式, p_引数);
        }
    }
}
