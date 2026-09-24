using System.Text;
using Tsumiki.Commons;

namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 持ち越し検証で見つかった、未観測の r-mer が続く長さの度数
    /// </summary>
    /// <remarks>
    /// 継ぎ目に掛かるものとそれ以外を分けて数えるのは、下限が両者の分布の谷にあるかを見るため
    /// </remarks>
    internal sealed class 連続長の度数
    {
        #region 定数

        /// <summary>
        /// 1 窓ずつ数える上限。これを超える長さは幅を持った区間にまとめる
        /// </summary>
        private const int 個別に数える上限 = 20;

        /// <summary>
        /// 個別に数える上限より長い連続をまとめる区間の上端
        /// </summary>
        private static readonly int[] まとめる区間の上端 = [30, 50, 100, int.MaxValue];

        #endregion

        #region 内部変数

        /// <summary>
        /// 継ぎ目に掛からない連続の度数
        /// </summary>
        private readonly long[] _継ぎ目なし = new long[個別に数える上限 + まとめる区間の上端.Length];

        /// <summary>
        /// 継ぎ目に掛かる連続の度数
        /// </summary>
        private readonly long[] _継ぎ目あり = new long[個別に数える上限 + まとめる区間の上端.Length];

        /// <summary>
        /// 連続 1 本ごとの位置 (配列名・窓の開始・窓の終了・継ぎ目に掛かるか)
        /// </summary>
        private readonly List<(string A_配列名, int A_開始, int A_終了, bool A_Is継ぎ目, int A_内側, int A_外側)> _明細 = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 連続 1 本を数える
        /// </summary>
        /// <param name="p_配列名">連続が見つかった配列の名前</param>
        /// <param name="p_開始">窓の開始位置</param>
        /// <param name="p_終了">窓の終了位置 (含む)</param>
        /// <param name="p_Is継ぎ目">分岐のある継ぎ目に掛かるか</param>
        /// <param name="p_内側">連続に掛かる k-mer のカバレッジの中央値</param>
        /// <param name="p_外側">連続の前後の k-mer のカバレッジの中央値</param>
        public void V_加算(string p_配列名, int p_開始, int p_終了, bool p_Is継ぎ目, int p_内側 = 0, int p_外側 = 0)
        {
            var l_表 = p_Is継ぎ目 ? this._継ぎ目あり : this._継ぎ目なし;
            l_表[Get_区間番号(p_終了 - p_開始 + 1)]++;
            this._明細.Add((p_配列名, p_開始, p_終了, p_Is継ぎ目, p_内側, p_外側));
        }

        /// <summary>
        /// 連続ごとの位置を TSV に書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        /// <remarks>
        /// 度数だけでは継ぎ目に掛からない長い連続が正しい配列か誤結合かを決められないので、参照へ当てられるよう位置を残す
        /// </remarks>
        public void V_書き出し(string p_パス)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            l_書き込み.WriteLine("name\tstart\tend\twindows\tjunction\tcov_in\tcov_flank");
            foreach (var (l_名前, l_開始, l_終了, l_Is継ぎ目, l_内側, l_外側) in this._明細)
            {
                l_書き込み.WriteLine(FormattableString.Invariant($"{l_名前}\t{l_開始}\t{l_終了}\t{l_終了 - l_開始 + 1}\t{(l_Is継ぎ目 ? 1 : 0)}\t{l_内側}\t{l_外側}"));
            }
        }

        /// <summary>
        /// 度数をログへ出す
        /// </summary>
        /// <param name="p_k長">持ち越し元の k</param>
        public void V_出力(int p_k長)
        {
            Logger.V_出力_そのまま(FormattableString.Invariant($"[Info] 未観測 r-mer の連続長 k={p_k長} 継ぎ目なし: {Get_文字列(this._継ぎ目なし)}"));
            Logger.V_出力_そのまま(FormattableString.Invariant($"[Info] 未観測 r-mer の連続長 k={p_k長} 継ぎ目あり: {Get_文字列(this._継ぎ目あり)}"));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 長さが入る区間の番号
        /// </summary>
        /// <param name="p_長さ">窓の数 (1 以上)</param>
        /// <returns></returns>
        private static int Get_区間番号(int p_長さ)
        {
            if (p_長さ <= 個別に数える上限)
            {
                return p_長さ - 1;
            }

            for (var i = 0; i < まとめる区間の上端.Length; i++)
            {
                if (p_長さ <= まとめる区間の上端[i])
                {
                    return 個別に数える上限 + i;
                }
            }
            return 個別に数える上限 + まとめる区間の上端.Length - 1;
        }

        /// <summary>
        /// 度数を「長さ:本数」の並びにする (0 本の区間は省く)
        /// </summary>
        /// <param name="p_表"></param>
        /// <returns></returns>
        private static string Get_文字列(long[] p_表)
        {
            var l_文 = new StringBuilder();
            for (var i = 0; i < p_表.Length; i++)
            {
                if (p_表[i] == 0)
                {
                    continue;
                }

                string l_ラベル;
                if (i < 個別に数える上限)
                {
                    l_ラベル = FormattableString.Invariant($"{i + 1}");
                }
                else
                {
                    var j = i - 個別に数える上限;
                    var l_下端 = (j == 0 ? 個別に数える上限 : まとめる区間の上端[j - 1]) + 1;
                    l_ラベル = まとめる区間の上端[j] == int.MaxValue
                        ? FormattableString.Invariant($"{l_下端}+")
                        : FormattableString.Invariant($"{l_下端}-{まとめる区間の上端[j]}");
                }
                _ = l_文.Append(FormattableString.Invariant($"{l_ラベル}:{p_表[i]} "));
            }
            return l_文.Length == 0 ? "なし" : l_文.ToString().TrimEnd();
        }

        #endregion
    }
}
