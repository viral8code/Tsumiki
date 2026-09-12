using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Xunit;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// 文言カタログの検査
    /// </summary>
    /// <remarks>
    /// ID を増やしたときに訳を入れ忘れる、あるいは訳の差し込み位置が原文とずれる、といった取りこぼしを防ぐ
    /// </remarks>
    public class MessageCatalogTests
    {
        #region 公開メソッド

        /// <summary>
        /// 全てのメッセージ ID に全ての言語の訳がある
        /// </summary>
        [Fact]
        public void V_全てのIDに全ての言語の文言がある()
        {
            var l_欠け = (from l_言語 in Enum.GetValues<言語>()
                          from l_ID in Get_全ID()
                          where !MessageCatalog.Has訳(l_言語, l_ID)
                          select $"{l_言語}/{l_ID}").ToList();

            Assert.Empty(l_欠け);
        }

        /// <summary>
        /// 各言語の訳の差し込み位置が英語の書式と一致する
        /// </summary>
        [Fact]
        public void V_訳の差し込み位置が英語と一致する()
        {
            foreach (var l_ID in Get_全ID())
            {
                var l_英語 = Get_差し込み位置(MessageCatalog.Get_書式(言語.英語, l_ID));
                foreach (var l_言語 in Enum.GetValues<言語>())
                {
                    Assert.Equal(l_英語, Get_差し込み位置(MessageCatalog.Get_書式(l_言語, l_ID)));
                }
            }
        }

        /// <summary>
        /// 文言に引数を差し込める
        /// </summary>
        [Fact]
        public void V_引数を差し込める()
        {
            foreach (var l_言語 in Enum.GetValues<言語>())
            {
                Messages.A_言語 = l_言語;
                var l_文言 = Messages.Get_文言(メッセージID.採用したk, 63);

                Assert.Contains("63", l_文言);
            }
            Messages.A_言語 = 言語.日本語;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// メッセージ ID をすべて返す
        /// </summary>
        /// <returns>メッセージ ID</returns>
        private static IEnumerable<メッセージID> Get_全ID()
        {
            return Enum.GetValues<メッセージID>();
        }

        /// <summary>
        /// 差し込み位置 ({0} など) の並び
        /// </summary>
        /// <remarks>
        /// 書式指定は無視する
        /// </remarks>
        /// <param name="p_書式"></param>
        /// <returns></returns>
        private static List<int> Get_差し込み位置(string p_書式)
        {
            return [.. Regex.Matches(p_書式, @"(?<!\{)\{(\d+)[^}]*\}")
                .Select(x => int.Parse(x.Groups[1].Value))
                .Distinct()
                .Order()];
        }

        #endregion

    }
}
