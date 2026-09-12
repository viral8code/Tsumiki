namespace Tsumiki.Cores.Evidence
{
    /// <summary>
    /// ペアエンドの隣接証拠を、観測本数ではなく期待本数との比で測るためのモデル
    /// </summary>
    /// <remarks>
    /// 観測本数をそのまま固定の下限と比べると、辺が長く距離が近いほど多く観測されるという幾何的な偏りをそのまま拾う<br/>
    /// 期待本数が数本しかない場所と数百本ある場所に同じ下限を課しても、どこに線を引いても正しくならない
    /// </remarks>
    internal sealed class PairedDistanceModel
    {
        #region 定数

        /// <summary>
        /// フラグメント長の経験分布を刻むビン幅
        /// </summary>
        private const int フラグメント長のビン幅 = 5;

        /// <summary>
        /// 既知長のばらつき幅の下限
        /// </summary>
        private const int 既知長のばらつき幅の下限 = 25;

        #endregion

        #region 内部変数

        /// <summary>
        /// フラグメント長の経験分布
        /// </summary>
        /// <remarks>
        /// 裾は誤マップなので両端を落とす
        /// </remarks>
        private readonly (int A_長さ, double A_確率)[] _分布;

        /// <summary>
        /// リード長
        /// </summary>
        private readonly int _リード長;

        /// <summary>
        /// 同一のギャップ長から出たとみなす既知長のばらつきの幅
        /// </summary>
        private readonly int _窓幅;

        /// <summary>
        /// 中央フラグメント長
        /// </summary>
        private readonly int _中央フラグメント長;

        #endregion

        #region プロパティ

        /// <summary>
        /// 使えるか
        /// </summary>
        public bool A_Is使用可能 => this._分布.Length > 0;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// フラグメント長標本から経験分布を組み立てる
        /// </summary>
        /// <param name="p_フラグメント長標本">フラグメント長の標本</param>
        /// <param name="p_リード長">リード長</param>
        public PairedDistanceModel(IReadOnlyList<int> p_フラグメント長標本, int p_リード長)
        {
            this._リード長 = Math.Max(1, p_リード長);

            var l_並び = p_フラグメント長標本.Where(x => x > 0).OrderBy(x => x).ToArray();
            if (l_並び.Length == 0)
            {
                this._分布 = [];
                this._窓幅 = 1;
                this._中央フラグメント長 = 0;
                return;
            }

            var l_下限 = Get_分位(l_並び, 0.01D);
            var l_上限 = Get_分位(l_並び, 0.99D);
            this._中央フラグメント長 = Get_分位(l_並び, 0.5D);
            this._窓幅 = Math.Max((l_上限 - l_下限) / 2, 既知長のばらつき幅の下限);

            var l_件数 = new Dictionary<int, int>();
            var l_総数 = 0;
            foreach (var l_長さ in l_並び)
            {
                if (l_長さ < l_下限 || l_長さ > l_上限)
                {
                    continue;
                }
                var l_ビン = l_長さ / フラグメント長のビン幅 * フラグメント長のビン幅;
                l_件数[l_ビン] = l_件数.GetValueOrDefault(l_ビン) + 1;
                l_総数++;
            }

            this._分布 = l_総数 == 0 ? [] : [.. l_件数.OrderBy(x => x.Key).Select(x => (x.Key, (double)x.Value / l_総数))];
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 既知長標本の中で最も密集した窓を選び、その本数と、そこから導かれるギャップ長を返す
        /// </summary>
        /// <param name="p_既知長標本">既知長の標本</param>
        /// <returns>最も密集した窓の本数と、そこから導かれるギャップ長</returns>
        /// <remarks>
        /// 中央値ではなく密集した窓を採るのは、反復の別コピーへ誤マップしたペアが長い裾を作るため<br/>
        /// 裾が過半を占めても、同じ隣接から出たペアは狭い範囲に集まるので峰は残る
        /// </remarks>
        public (int A_本数, int A_ギャップ長) Get_一貫した支持(IReadOnlyList<int> p_既知長標本)
        {
            if (p_既知長標本.Count == 0)
            {
                return (0, 0);
            }

            var l_並び = p_既知長標本.OrderBy(x => x).ToArray();
            var l_幅 = 2 * this._窓幅;
            var l_最良開始 = 0;
            var l_最良数 = 0;
            var l_右 = 0;
            for (var l_左 = 0; l_左 < l_並び.Length; l_左++)
            {
                if (l_右 < l_左)
                {
                    l_右 = l_左;
                }
                while (l_右 < l_並び.Length && l_並び[l_右] - l_並び[l_左] <= l_幅)
                {
                    l_右++;
                }
                if (l_右 - l_左 > l_最良数)
                {
                    l_最良数 = l_右 - l_左;
                    l_最良開始 = l_左;
                }
            }

            var l_代表 = l_並び[l_最良開始 + (l_最良数 / 2)];
            return (l_最良数, this._中央フラグメント長 - l_代表);
        }

        /// <summary>
        /// ギャップ長 p_ギャップ長 で隣り合う長さ p_長さ1・p_長さ2 の配列に跨がりうるフラグメントの開始位置の総数
        /// </summary>
        /// <param name="p_長さ1">片側の長さ</param>
        /// <param name="p_長さ2">もう片側の長さ</param>
        /// <param name="p_ギャップ長">両者の間のギャップ長</param>
        /// <returns>フラグメントの開始位置の総数</returns>
        /// <remarks>
        /// フラグメント開始位置の密度を掛けると期待ペア数になる<br/>
        /// 接合点を原点とし、フラグメント長 x の開始位置 s について左リードが左側に収まる条件 s &lt;= -リード長 かつ s &gt;= -長さ 1、右リードが右側に収まる条件 s &gt;= ギャップ長 + リード長 - x かつ s &lt;= ギャップ長 + 長さ 2 - x を満たす s の個数を数える
        /// </remarks>
        public double Get_期待位置数(long p_長さ1, long p_長さ2, int p_ギャップ長)
        {
            var l_合計 = 0D;
            foreach (var (l_長さ, l_確率) in this._分布)
            {
                var l_下 = Math.Max(-p_長さ1, p_ギャップ長 + this._リード長 - l_長さ);
                var l_上 = Math.Min(-this._リード長, p_ギャップ長 + p_長さ2 - l_長さ);
                if (l_上 >= l_下)
                {
                    l_合計 += l_確率 * (l_上 - l_下 + 1L);
                }
            }
            return l_合計;
        }

        /// <summary>
        /// 1 本の配列の内側に両端が収まるフラグメントの開始位置の総数
        /// </summary>
        /// <param name="p_長さ">配列の長さ</param>
        /// <returns>フラグメントの開始位置の総数</returns>
        /// <remarks>
        /// 同一 unitig 内の観測数からフラグメント開始位置の密度を較正するのに使う
        /// </remarks>
        public double Get_期待位置数_単一(long p_長さ)
        {
            var l_合計 = 0D;
            foreach (var (l_長さ, l_確率) in this._分布)
            {
                var l_個数 = p_長さ - l_長さ + 1L;
                if (l_個数 > 0L)
                {
                    l_合計 += l_確率 * l_個数;
                }
            }
            return l_合計;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 昇順に並んだ値から分位点を返す
        /// </summary>
        /// <param name="p_昇順">昇順に並んだ値</param>
        /// <param name="p_位置">求める分位 (0 から 1) </param>
        /// <returns>分位点</returns>
        private static int Get_分位(int[] p_昇順, double p_位置)
        {
            var l_i = (int)(p_昇順.Length * p_位置);
            return p_昇順[Math.Clamp(l_i, 0, p_昇順.Length - 1)];
        }

        #endregion
    }
}
