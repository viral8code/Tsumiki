namespace Tsumiki.Cores.Evidence
{
    /// <summary>
    /// ペアエンドの隣接証拠を、生の観測本数ではなく期待本数との比で測るための較正器
    /// </summary>
    internal sealed class 証拠較正器
    {
        #region 定数

        /// <summary>
        /// 支持本数を 0〜1 の確信度へ潰すときの時定数
        /// </summary>
        public const double 飽和の時定数 = 3.0D;

        #endregion

        #region 内部変数

        /// <summary>
        /// モデル
        /// </summary>
        private readonly PairedDistanceModel? _モデル;

        /// <summary>
        /// 密度
        /// </summary>
        private readonly double _密度;

        #endregion

        #region プロパティ

        /// <summary>
        /// モデルが構築できたか
        /// </summary>
        public bool A_Is使用可能 => this._モデル is not null;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// モデルと密度から較正器を組み立てる
        /// </summary>
        /// <param name="p_モデル">距離モデル、使えない場合は null</param>
        /// <param name="p_密度">観測密度</param>
        private 証拠較正器(PairedDistanceModel? p_モデル, double p_密度)
        {
            this._モデル = p_モデル;
            this._密度 = p_密度;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 独立な支持本数を 0〜1 の確信度に変換する
        /// </summary>
        /// <param name="p_独立支持数">独立とみなせる支持本数</param>
        /// <returns></returns>
        public static double Get_飽和支持(double p_独立支持数)
        {
            return p_独立支持数 <= 0D ? 0D : 1D - Math.Exp(-p_独立支持数 / 飽和の時定数);
        }

        /// <summary>
        /// 同一 unitig 内標本とリード長から較正器を作る
        /// </summary>
        /// <param name="p_同一unitig標本">同一 unitig 内で観測された距離の標本</param>
        /// <param name="p_リード長">リード長、不明な場合は null</param>
        /// <param name="p_unitig長一覧">unitig ごとの長さ</param>
        /// <returns>較正器、モデルが使えない場合も返り値自体は null にならない</returns>
        public static 証拠較正器 Get_較正器(IReadOnlyList<int> p_同一unitig標本, int? p_リード長, IEnumerable<long> p_unitig長一覧)
        {
            if (p_リード長 is not { } l_リード長 || p_同一unitig標本.Count == 0)
            {
                return new 証拠較正器(null, 0D);
            }

            var l_モデル = new PairedDistanceModel(p_同一unitig標本, l_リード長);
            if (!l_モデル.A_Is使用可能)
            {
                return new 証拠較正器(null, 0D);
            }

            var l_期待位置数合計 = 0D;
            foreach (var l_長さ in p_unitig長一覧)
            {
                l_期待位置数合計 += l_モデル.Get_期待位置数_単一(l_長さ);
            }
            if (l_期待位置数合計 <= 0D)
            {
                return new 証拠較正器(null, 0D);
            }

            var l_密度 = p_同一unitig標本.Count / l_期待位置数合計;
            return l_密度 > 0D ? new 証拠較正器(l_モデル, l_密度) : new 証拠較正器(null, 0D);
        }

        /// <summary>
        /// 観測本数を期待本数 (密度 x 期待位置数) で割った比
        /// </summary>
        /// <param name="p_観測本数">観測された本数</param>
        /// <param name="p_長さ1">片側の長さ</param>
        /// <param name="p_長さ2">もう片側の長さ</param>
        /// <param name="p_ギャップ長">両者の間のギャップ長</param>
        /// <returns>正規化済みの支持</returns>
        public double Get_正規化済み支持(ulong p_観測本数, long p_長さ1, long p_長さ2, int p_ギャップ長)
        {
            var l_期待 = this.Get_期待本数(p_長さ1, p_長さ2, p_ギャップ長);
            return l_期待 > 0D ? p_観測本数 / l_期待 : 0D;
        }

        /// <summary>
        /// 2 本の端の間に、このライブラリのペアが何本跨ぐはずか
        /// </summary>
        /// <param name="p_長さ1">片側の長さ</param>
        /// <param name="p_長さ2">もう片側の長さ</param>
        /// <param name="p_ギャップ長">両者の間のギャップ長</param>
        /// <returns>期待本数、モデルが使えなければ 0</returns>
        public double Get_期待本数(long p_長さ1, long p_長さ2, int p_ギャップ長)
        {
            return this._モデル is { } l_モデル ? this._密度 * l_モデル.Get_期待位置数(p_長さ1, p_長さ2, p_ギャップ長) : 0D;
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// 観測された既知長の一覧から、独立とみなせる支持本数を数える
        /// </summary>
        /// <param name="p_既知長一覧">観測された既知長の一覧</param>
        /// <returns></returns>
        public static int Get_独立支持数(IReadOnlyList<int> p_既知長一覧)
        {
            if (p_既知長一覧.Count == 0)
            {
                return 0;
            }

            HashSet<int> l_相異なる距離 = [];
            foreach (var l_距離 in p_既知長一覧)
            {
                _ = l_相異なる距離.Add(l_距離);
            }
            return l_相異なる距離.Count;
        }

        #endregion
    }
}
