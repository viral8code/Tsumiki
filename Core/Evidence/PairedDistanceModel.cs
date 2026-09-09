using Tsumiki.Common;

namespace Tsumiki.Core.Evidence
{
    /// <summary>
    /// ペアエンドの隣接証拠を、観測本数ではなく期待本数との比で測るためのモデル<br/>
    /// 観測本数をそのまま固定の下限と比べると、辺が長く距離が近いほど多く
    /// 観測されるという幾何的な偏りをそのまま拾う<br/>
    /// 期待本数が数本しかない場所と
    /// 数百本ある場所に同じ下限を課しても、どこに線を引いても正しくならない
    /// </summary>
    internal sealed class PairedDistanceModel
    {
        /// <summary>
        /// フラグメント長の経験分布<br/>
        /// 裾は誤マップなので両端を落とす
        /// </summary>
        private readonly (int A_長さ, double A_確率)[] _分布;

        private readonly int _リード長;

        /// <summary>
        /// 同一のギャップ長から出たとみなす既知長のばらつきの幅
        /// </summary>
        private readonly int _窓幅;

        private readonly int _中央フラグメント長;

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

            var l_下限 = Get_分位(l_並び, 0.01);
            var l_上限 = Get_分位(l_並び, 0.99);
            this._中央フラグメント長 = Get_分位(l_並び, 0.5);
            this._窓幅 = Math.Max((l_上限 - l_下限) / 2, Consts.既知長のばらつき幅の下限);

            var l_件数 = new Dictionary<int, int>();
            var l_総数 = 0;
            foreach (var l_長さ in l_並び)
            {
                if (l_長さ < l_下限 || l_長さ > l_上限)
                {
                    continue;
                }
                var l_ビン = l_長さ / Consts.フラグメント長のビン幅 * Consts.フラグメント長のビン幅;
                l_件数[l_ビン] = l_件数.GetValueOrDefault(l_ビン) + 1;
                l_総数++;
            }

            this._分布 = l_総数 == 0
                ? []
                : [.. l_件数.OrderBy(x => x.Key).Select(x => (x.Key, (double)x.Value / l_総数))];
        }

        public bool A_使えるか => this._分布.Length > 0;

        private static int Get_分位(int[] p_昇順, double p_位置)
        {
            var i = (int)(p_昇順.Length * p_位置);
            return p_昇順[Math.Clamp(i, 0, p_昇順.Length - 1)];
        }

        /// <summary>
        /// 既知長標本の中で最も密集した窓を選び、その本数と、そこから導かれる
        /// ギャップ長を返す<br/>
        /// 中央値ではなく密集した窓を採るのは、反復の別コピーへ誤マップした
        /// ペアが長い裾を作るため<br/>
        /// 裾が過半を占めても、同じ隣接から出た
        /// ペアは狭い範囲に集まるので峰は残る
        /// </summary>
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
        /// ギャップ長 p_ギャップ長 で隣り合う長さ p_長さ1・p_長さ2 の配列に跨がりうる
        /// フラグメントの開始位置の総数<br/>
        /// フラグメント開始位置の密度を掛けると
        /// 期待ペア数になる<br/>
        /// 接合点を原点とし、フラグメント長 x の開始位置 s について
        /// 左リードが左側に収まる条件 s &lt;= -リード長 かつ s &gt;= -長さ 1、
        /// 右リードが右側に収まる条件 s &gt;= ギャップ長 + リード長 - x かつ
        /// s &lt;= ギャップ長 + 長さ 2 - x を満たす s の個数を数える
        /// </summary>
        public double Get_期待位置数(long p_長さ1, long p_長さ2, int p_ギャップ長)
        {
            double l_合計 = 0D;
            foreach (var (l_長さ, l_確率) in this._分布)
            {
                var l_下 = Math.Max(-p_長さ1, p_ギャップ長 + this._リード長 - l_長さ);
                var l_上 = Math.Min(-this._リード長, p_ギャップ長 + p_長さ2 - l_長さ);
                if (l_上 >= l_下)
                {
                    l_合計 += l_確率 * (l_上 - l_下 + 1);
                }
            }
            return l_合計;
        }

        /// <summary>
        /// 1 本の配列の内側に両端が収まるフラグメントの開始位置の総数<br/>
        /// 同一 unitig 内の観測数からフラグメント開始位置の密度を較正するのに使う
        /// </summary>
        public double Get_期待位置数_単一(long p_長さ)
        {
            double l_合計 = 0D;
            foreach (var (l_長さ, l_確率) in this._分布)
            {
                var l_個数 = p_長さ - l_長さ + 1;
                if (l_個数 > 0)
                {
                    l_合計 += l_確率 * l_個数;
                }
            }
            return l_合計;
        }
    }
}
