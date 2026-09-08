namespace Tsumiki.Model.Polishing
{
    /// <summary>
    /// 最終配列へのリード再マッピングと、その多数決による置換訂正の集計。
    /// 位置ごとの深度もここで得られるため、カバレッジの連続性の判定にも使う。
    /// </summary>
    internal readonly record struct ポリッシュ統計(
        int A_配列数,
        long A_総延長,
        long A_マップされたリード数,
        long A_棄却されたリード数,
        long A_訂正した塩基数,
        long A_深度不足の位置数,
        long A_評価できた位置数,
        double A_深度の中央値)
    {
        /// <summary>
        /// 深度が期待から大きく落ち込んだ位置の割合。連結の裏付けが
        /// 無い箇所は、その接合点の前後で深度が不連続になる。
        /// </summary>
        public double A_深度不足率 => this.A_評価できた位置数 == 0
            ? 0
            : (double)this.A_深度不足の位置数 / this.A_評価できた位置数;
    }
}
