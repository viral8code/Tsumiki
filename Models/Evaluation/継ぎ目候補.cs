namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// 反復の区間を挟んで両側の固有配列を繋いだ箇所 (または N のギャップ) と、それを裏付けるリードの証拠
    /// </summary>
    /// <param name="A_配列名"></param>
    /// <param name="A_開始">反復の区間の始まり (0 始まり)</param>
    /// <param name="A_終了">反復の区間の終わり (含まない)</param>
    /// <param name="A_Isギャップ">区間に N を含むか</param>
    /// <param name="A_コピー数">区間の深さの中央値を、固有配列の深さの中央値で割った値</param>
    /// <param name="A_跨ぐ読み">区間の前後を 1 本で読んだリードの数</param>
    /// <param name="A_跨ぐ組">断片が区間の前後を覆う、向きと長さの正しいペアの数</param>
    /// <param name="A_期待の組">断片長の分布と両側の深さから見込む、跨ぐ組の数</param>
    /// <param name="A_外れ錨">区間の手前にあって区間の側を向くのに、相方が正しく組めていない一意なリードの数</param>
    /// <param name="A_切れ端">区間の両端で途中から整列が切れたリードの数を、深さの中央値で割った値</param>
    /// <param name="A_左の深さ比">区間の左 100 塩基の一意な深さを、深さの中央値で割った値</param>
    /// <param name="A_右の深さ比">区間の右 100 塩基の一意な深さを、深さの中央値で割った値</param>
    internal readonly record struct 継ぎ目候補(string A_配列名, int A_開始, int A_終了, bool A_Isギャップ, double A_コピー数, int A_跨ぐ読み, int A_跨ぐ組, double A_期待の組, int A_外れ錨, double A_切れ端, double A_左の深さ比, double A_右の深さ比)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 反復の区間の長さ
        /// </summary>
        public int A_反復長 => this.A_終了 - this.A_開始;

        /// <summary>
        /// 跨ぐ組を期待で割った値、期待が 0 なら -1
        /// </summary>
        public double A_組の比 => this.A_期待の組 > 0D ? this.A_跨ぐ組 / this.A_期待の組 : -1D;

        /// <summary>
        /// 外れ錨が、外れ錨と跨ぐ組の証拠の合計に占める割合
        /// </summary>
        public double A_外れ割合 => this.A_外れ錨 / (this.A_外れ錨 + (2D * this.A_跨ぐ組) + 1D);

        /// <summary>
        /// 両側の深さ比の低いほう
        /// </summary>
        public double A_低い側の深さ比 => Math.Min(this.A_左の深さ比, this.A_右の深さ比);

        #endregion
    }
}
