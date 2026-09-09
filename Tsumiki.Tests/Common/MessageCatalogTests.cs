using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System;
using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Xunit;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// 文言カタログの検査<br/>
    /// ID を増やしたときに訳を入れ忘れる、あるいは
    /// 訳の差し込み位置が原文とずれる、といった取りこぼしを防ぐ
    /// </summary>
    public class MessageCatalogTests
    {
        private static IEnumerable<メッセージID> Get_全ID()
        {
            return Enum.GetValues<メッセージID>();
        }

        /// <summary>
        /// 差し込み位置({0} など)の並び<br/>
        /// 書式指定は無視する
        /// </summary>
        private static List<int> Get_差し込み位置(string p_書式)
        {
            return [.. Regex.Matches(p_書式, @"(?<!\{)\{(\d+)[^}]*\}")
                .Select(x => int.Parse(x.Groups[1].Value))
                .Distinct()
                .Order()];
        }

        [Fact]
        public void 全てのIDに全ての言語の文言がある()
        {
            var l_欠け = (from l_言語 in Enum.GetValues<言語>()
                          from l_ID in Get_全ID()
                          where !MessageCatalog.Get_訳があるか(l_言語, l_ID)
                          select $"{l_言語}/{l_ID}").ToList();

            Assert.Empty(l_欠け);
        }

        [Fact]
        public void 訳の差し込み位置が英語と一致する()
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

        [Fact]
        public void 引数を差し込める()
        {
            foreach (var l_言語 in Enum.GetValues<言語>())
            {
                Messages.A_言語 = l_言語;
                var l_文言 = Messages.Get_文言(メッセージID.採用したk, 63);

                Assert.Contains("63", l_文言);
            }
            Messages.A_言語 = 言語.日本語;
        }
    }
}
