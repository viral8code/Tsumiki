using System.Globalization;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.IO
{
    /// <summary>
    /// 継ぎ目の誤り確率モデルのファイル (1 行 1 項目、空白区切り) を読む
    /// </summary>
    internal static class JunctionRiskModelReader
    {
        #region 定数

        /// <summary>
        /// 注釈の行の印
        /// </summary>
        private const char C_注釈の印 = '#';

        /// <summary>
        /// メッセージの中の一覧の区切り
        /// </summary>
        private const string C_一覧の区切り = ", ";

        /// <summary>
        /// 項目の区切り
        /// </summary>
        private static readonly char[] C_区切り = [' ', '\t'];

        /// <summary>
        /// 数の項目ごとの値の数
        /// </summary>
        private static readonly Dictionary<string, int> C_値の数 = new()
        {
            [継ぎ目モデル.C_項目_得点の係数] = 継ぎ目モデル.C_得点の項の数,
            [継ぎ目モデル.C_項目_得点の平均] = 継ぎ目モデル.C_得点の項の数,
            [継ぎ目モデル.C_項目_得点の幅] = 継ぎ目モデル.C_得点の項の数,
            [継ぎ目モデル.C_項目_ギャップの較正] = 継ぎ目モデル.C_深さつき較正の項の数,
            [継ぎ目モデル.C_項目_跨げる組の較正] = 継ぎ目モデル.C_跨げる組の較正の項の数,
            [継ぎ目モデル.C_項目_跨げない組の較正] = 継ぎ目モデル.C_深さつき較正の項の数
        };

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ファイルからモデルを読む
        /// </summary>
        /// <param name="p_パス">モデルのファイル</param>
        /// <returns>読んだモデル</returns>
        public static 継ぎ目モデル Get_読込(string p_パス)
        {
            var l_出所 = Path.GetFullPath(p_パス);
            return Get_解釈(File.ReadAllText(l_出所), l_出所);
        }

        /// <summary>
        /// ファイルがあればモデルを読む (無ければ null、あって壊れていれば例外)
        /// </summary>
        /// <param name="p_パス">探すファイル</param>
        /// <returns>読んだモデル、ファイルが無ければ null</returns>
        public static 継ぎ目モデル? Get_読込_あれば(string p_パス)
        {
            return File.Exists(p_パス) ? Get_読込(p_パス) : null;
        }

        /// <summary>
        /// モデルの本文を解釈する
        /// </summary>
        /// <param name="p_本文">モデルの本文</param>
        /// <param name="p_出所">メッセージと記録に使う出所</param>
        /// <returns>解釈したモデル</returns>
        public static 継ぎ目モデル Get_解釈(string p_本文, string p_出所)
        {
            string? l_名前 = null;
            Dictionary<string, double[]> l_値群 = [];
            var l_行群 = p_本文.TrimStart('﻿').ReplaceLineEndings("\n").Split('\n');
            for (var i = 0; i < l_行群.Length; i++)
            {
                var l_行 = l_行群[i].Trim();
                if (l_行.Length == 0 || l_行[0] == C_注釈の印)
                {
                    continue;
                }

                var l_欄 = l_行.Split(C_区切り, StringSplitOptions.RemoveEmptyEntries);
                var l_項目 = l_欄[0];
                if (l_項目 == 継ぎ目モデル.C_項目_名前)
                {
                    if (l_名前 is not null)
                    {
                        throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"\"{継ぎ目モデル.C_項目_名前}\" appears more than once"));
                    }

                    if (l_欄.Length != 2)
                    {
                        throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"\"{継ぎ目モデル.C_項目_名前}\" needs exactly one word"));
                    }

                    l_名前 = l_欄[1];
                    continue;
                }

                if (!C_値の数.TryGetValue(l_項目, out var l_数))
                {
                    throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"unknown item \"{l_項目}\"; expected one of {継ぎ目モデル.C_項目_名前}{C_一覧の区切り}{string.Join(C_一覧の区切り, C_値の数.Keys)}"));
                }

                if (l_値群.ContainsKey(l_項目))
                {
                    throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"\"{l_項目}\" appears more than once"));
                }

                if (l_欄.Length - 1 != l_数)
                {
                    throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"\"{l_項目}\" needs {l_数} values, found {l_欄.Length - 1}"));
                }

                var l_値 = new double[l_数];
                for (var j = 0; j < l_数; j++)
                {
                    if (!double.TryParse(l_欄[j + 1], NumberStyles.Float, CultureInfo.InvariantCulture, out l_値[j]) || !double.IsFinite(l_値[j]))
                    {
                        throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"\"{l_欄[j + 1]}\" in \"{l_項目}\" is not a finite number"));
                    }
                }

                if (l_項目 == 継ぎ目モデル.C_項目_得点の幅 && l_値.Any(static x => x <= 0D))
                {
                    throw new InvalidDataException(Get_行のメッセージ(p_出所, i, $"every value in \"{継ぎ目モデル.C_項目_得点の幅}\" must be positive"));
                }

                l_値群[l_項目] = l_値;
            }

            if (l_名前 is null)
            {
                throw new InvalidDataException($"{p_出所}: missing item \"{継ぎ目モデル.C_項目_名前}\"");
            }

            if (C_値の数.Keys.FirstOrDefault(x => !l_値群.ContainsKey(x)) is { } l_欠けた項目)
            {
                throw new InvalidDataException($"{p_出所}: missing item \"{l_欠けた項目}\"");
            }

            return new 継ぎ目モデル(l_名前, p_出所, l_値群[継ぎ目モデル.C_項目_得点の係数], l_値群[継ぎ目モデル.C_項目_得点の平均], l_値群[継ぎ目モデル.C_項目_得点の幅], l_値群[継ぎ目モデル.C_項目_ギャップの較正], l_値群[継ぎ目モデル.C_項目_跨げる組の較正], l_値群[継ぎ目モデル.C_項目_跨げない組の較正]);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 行番号つきのメッセージを作る
        /// </summary>
        /// <param name="p_出所">モデルの出所</param>
        /// <param name="p_行番号">0 始まりの行番号</param>
        /// <param name="p_内容">誤りの内容</param>
        /// <returns>メッセージ</returns>
        private static string Get_行のメッセージ(string p_出所, int p_行番号, string p_内容)
        {
            return $"{p_出所} line {p_行番号 + 1}: {p_内容}";
        }

        #endregion
    }
}
