using System.Globalization;

namespace Tsumiki.Models.Correction
{
    /// <summary>
    /// エラー訂正の前後で、リードの中の誤りの無い区間 (信頼できる窓が続く範囲) の長さを数え、各 k のカバレッジを見積もる
    /// </summary>
    /// <param name="p_k長">数えた窓の k 長</param>
    internal sealed class 無誤り区間の度数(int p_k長)
    {
        #region 定数

        /// <summary>
        /// 項目 R
        /// </summary>
        private const string C_項目_R = "R";

        /// <summary>
        /// 数える区間の長さの上限 (超えたものはここに入れる)
        /// </summary>
        private const int C_最大長 = 4_095;

        /// <summary>
        /// 保存ファイルの列の区切り
        /// </summary>
        private const char C_列区切り = '\t';

        #endregion

        #region 内部変数

        /// <summary>
        /// 長さごとの訂正前の区間の数
        /// </summary>
        private readonly long[] _訂正前 = new long[C_最大長 + 1];

        /// <summary>
        /// 長さごとの訂正後の区間の数
        /// </summary>
        private readonly long[] _訂正後 = new long[C_最大長 + 1];

        #endregion

        #region プロパティ

        /// <summary>
        /// 数えた窓の k 長
        /// </summary>
        public int A_k長 { get; } = p_k長;

        /// <summary>
        /// 訂正前のリードの k-mer スペクトルで見た単一コピーのカバレッジ
        /// </summary>
        public double A_単一コピー平均 { get; set; }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 リードの窓の信頼状況から、誤りの無い区間を数える
        /// </summary>
        /// <param name="p_信頼状況">窓ごとの信頼状況</param>
        /// <param name="p_Is訂正後">訂正後の状況か</param>
        public void V_追加(ReadOnlySpan<bool> p_信頼状況, bool p_Is訂正後)
        {
            var l_度数 = p_Is訂正後 ? this._訂正後 : this._訂正前;
            var l_連続 = 0;
            for (var w = 0; w <= p_信頼状況.Length; w++)
            {
                if (w < p_信頼状況.Length && p_信頼状況[w])
                {
                    l_連続++;
                    continue;
                }

                if (l_連続 > 0)
                {
                    _ = Interlocked.Increment(ref l_度数[Math.Min(l_連続 + this.A_k長 - 1, C_最大長)]);
                    l_連続 = 0;
                }
            }
        }

        /// <summary>
        /// 訂正後のリードでの k のカバレッジを見積もる (訂正前の単一コピーのカバレッジに、取れる k-mer の数の比を掛ける)
        /// </summary>
        /// <param name="p_k長">見積もる k 長 (数えた窓の k 長以上)</param>
        /// <returns>見積もったカバレッジ、訂正前に k-mer が無ければ 0</returns>
        public double Get_予測カバレッジ(int p_k長)
        {
            var l_訂正前の数 = Get_kmer数(this._訂正前, this.A_k長);
            return l_訂正前の数 <= 0D ? 0D : this.A_単一コピー平均 * Get_kmer数(this._訂正後, p_k長) / l_訂正前の数;
        }

        /// <summary>
        /// 書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        public void V_書き出し(string p_パス)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            l_書き込み.WriteLine(string.Join(C_列区切り, this.A_k長.ToString(CultureInfo.InvariantCulture), this.A_単一コピー平均.ToString(C_項目_R, CultureInfo.InvariantCulture)));
            for (var l_長さ = 0; l_長さ <= C_最大長; l_長さ++)
            {
                if (this._訂正前[l_長さ] > 0L || this._訂正後[l_長さ] > 0L)
                {
                    l_書き込み.WriteLine(string.Join(C_列区切り, l_長さ, this._訂正前[l_長さ], this._訂正後[l_長さ]));
                }
            }
        }

        /// <summary>
        /// 読み込む
        /// </summary>
        /// <param name="p_パス">読み込み元</param>
        /// <returns>読み込んだ度数、ファイルが無ければ null</returns>
        public static 無誤り区間の度数? Get_読込(string p_パス)
        {
            if (!File.Exists(p_パス))
            {
                return null;
            }

            var l_行群 = File.ReadAllLines(p_パス);
            var l_見出し = l_行群[0].Split(C_列区切り);
            var l_度数 = new 無誤り区間の度数(int.Parse(l_見出し[0], CultureInfo.InvariantCulture))
            {
                A_単一コピー平均 = double.Parse(l_見出し[1], CultureInfo.InvariantCulture),
            };
            foreach (var l_行 in l_行群.Skip(1))
            {
                var l_列 = l_行.Split(C_列区切り);
                var l_長さ = int.Parse(l_列[0], CultureInfo.InvariantCulture);
                l_度数._訂正前[l_長さ] = long.Parse(l_列[1], CultureInfo.InvariantCulture);
                l_度数._訂正後[l_長さ] = long.Parse(l_列[2], CultureInfo.InvariantCulture);
            }

            return l_度数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 区間の長さの度数から、長さ k の k-mer がいくつ取れるかを数える
        /// </summary>
        /// <param name="p_度数">長さごとの区間の数</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>取れる k-mer の数</returns>
        private static double Get_kmer数(long[] p_度数, int p_k長)
        {
            var l_合計 = 0D;
            for (var l_長さ = p_k長; l_長さ < p_度数.Length; l_長さ++)
            {
                l_合計 += (double)p_度数[l_長さ] * (l_長さ - p_k長 + 1);
            }

            return l_合計;
        }

        #endregion
    }
}
