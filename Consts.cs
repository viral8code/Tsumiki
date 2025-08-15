using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Tsumiki
{
    internal class Consts
    {
        public const string VERSION = "1.0";

        public static readonly List<string> AUTHOR_LIST = [
            "viral",
            ];

        public static readonly string DETAILS_TEXT = $"""
            This Program 'Tsumiki' is a genome assebly.
            author: {string.Join(", ", AUTHOR_LIST)}
            version: {VERSION}
            """;
    }
}
