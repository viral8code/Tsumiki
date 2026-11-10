using System.Globalization;
using System.Security.Cryptography;
using System.Text;

namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// 継ぎ目の誤りの確率を出すモデルの係数と、その出所
    /// </summary>
    /// <param name="A_名前">モデルの版の名前</param>
    /// <param name="A_出所">読み込んだファイルのパス、組み込みなら <see cref="C_組み込みの出所"/></param>
    /// <param name="A_得点の係数">得点の説明変数ごとの係数 (<see cref="C_得点の項の数"/> 個)</param>
    /// <param name="A_得点の平均">得点の説明変数ごとの標準化の平均 (<see cref="C_得点の項の数"/> 個)</param>
    /// <param name="A_得点の幅">得点の説明変数ごとの標準化の幅 (<see cref="C_得点の項の数"/> 個、正)</param>
    /// <param name="A_ギャップの較正">ギャップの組の較正 (切片、得点、低い側の深さ比、切れ端)</param>
    /// <param name="A_跨げる組の較正">組で跨げる組の較正 (切片、得点)</param>
    /// <param name="A_跨げない組の較正">組で跨げない組の較正 (切片、得点、低い側の深さ比、切れ端)</param>
    internal sealed record 継ぎ目モデル(string A_名前, string A_出所, double[] A_得点の係数, double[] A_得点の平均, double[] A_得点の幅, double[] A_ギャップの較正, double[] A_跨げる組の較正, double[] A_跨げない組の較正)
    {
        #region 定数

        /// <summary>
        /// 組み込みのモデルの出所
        /// </summary>
        public const string C_組み込みの出所 = "built-in";

        /// <summary>
        /// 得点の説明変数の数 (切片を含む)
        /// </summary>
        public const int C_得点の項の数 = 13;

        /// <summary>
        /// ギャップの組と、組で跨げない組の較正の項の数
        /// </summary>
        public const int C_深さつき較正の項の数 = 4;

        /// <summary>
        /// 組で跨げる組の較正の項の数
        /// </summary>
        public const int C_跨げる組の較正の項の数 = 2;

        /// <summary>
        /// 版の名前の項目
        /// </summary>
        public const string C_項目_名前 = "name";

        /// <summary>
        /// 得点の係数の項目
        /// </summary>
        public const string C_項目_得点の係数 = "score";

        /// <summary>
        /// 標準化の平均の項目
        /// </summary>
        public const string C_項目_得点の平均 = "mean";

        /// <summary>
        /// 標準化の幅の項目
        /// </summary>
        public const string C_項目_得点の幅 = "scale";

        /// <summary>
        /// ギャップの組の較正の項目
        /// </summary>
        public const string C_項目_ギャップの較正 = "gap";

        /// <summary>
        /// 組で跨げる組の較正の項目
        /// </summary>
        public const string C_項目_跨げる組の較正 = "span";

        /// <summary>
        /// 組で跨げない組の較正の項目
        /// </summary>
        public const string C_項目_跨げない組の較正 = "nospan";

        /// <summary>
        /// 正規の本文の値の区切り
        /// </summary>
        private const char C_値の区切り = '\t';

        /// <summary>
        /// 正規の本文の行の終わり
        /// </summary>
        private const char C_行の終わり = '\n';

        /// <summary>
        /// 値を元に戻せる数の書式
        /// </summary>
        private const string C_数の書式 = "R";

        #endregion

        #region プロパティ

        /// <summary>
        /// 正規の本文 (<see cref="Get_本文"/>) の SHA-256 (16 進の小文字、注釈・空白・改行の違いに依らない)
        /// </summary>
        public string A_SHA256 => Convert.ToHexStringLower(SHA256.HashData(Encoding.UTF8.GetBytes(this.Get_本文())));

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 数の項目を、ファイルに書く順に並べる
        /// </summary>
        /// <returns>項目名と値</returns>
        public (string A_項目, double[] A_値)[] Get_数の項目群()
        {
            return [(C_項目_得点の係数, this.A_得点の係数), (C_項目_得点の平均, this.A_得点の平均), (C_項目_得点の幅, this.A_得点の幅), (C_項目_ギャップの較正, this.A_ギャップの較正), (C_項目_跨げる組の較正, this.A_跨げる組の較正), (C_項目_跨げない組の較正, this.A_跨げない組の較正)];
        }

        /// <summary>
        /// 係数を、モデルのファイルの書式 (注釈なし、タブ区切り、値は元に戻せる書式) の本文にする
        /// </summary>
        /// <returns>正規の本文</returns>
        public string Get_本文()
        {
            var l_本文 = new StringBuilder().Append(C_項目_名前).Append(C_値の区切り).Append(this.A_名前).Append(C_行の終わり);
            foreach (var (l_項目, l_値) in this.Get_数の項目群())
            {
                _ = l_本文.Append(l_項目).Append(C_値の区切り).Append(string.Join(C_値の区切り, l_値.Select(static x => x.ToString(C_数の書式, CultureInfo.InvariantCulture)))).Append(C_行の終わり);
            }

            return l_本文.ToString();
        }

        #endregion
    }
}
