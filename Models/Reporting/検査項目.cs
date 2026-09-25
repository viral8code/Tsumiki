using Tsumiki.Models.Foundation;

namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 検査 1 項目
    /// </summary>
    /// <param name="A_キー"></param>
    /// <param name="A_見出し"></param>
    /// <param name="A_判定"></param>
    /// <param name="A_内訳"></param>
    internal readonly record struct 検査項目(string A_キー, メッセージID A_見出し, 検査判定 A_判定, string A_内訳);
}
