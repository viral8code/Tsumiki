namespace Tsumiki.Core
{
    /// <summary>
    /// ペアエンドの隣接証拠を、生の観測本数ではなく期待本数との比で測るための
    /// 較正器。ある辺で観測される本数は、その辺の両側の unitig 長と
    /// フラグメント長分布に応じて期待値そのものが桁で変わるため、
    /// 観測本数をそのまま固定閾値と比べると、短い辺には厳しすぎ、
    /// 長い辺には緩すぎる基準になる。
    ///
    /// Scaffolder(スキャフォールド辺の選択)・ContigMaker(分岐選択)・
    /// BeamSearchExtender(先読みスコア)が同じ較正を共有するために括り出した。
    /// </summary>
    internal sealed class 証拠較正器
    {
        private readonly PairedDistanceModel? _モデル;
        private readonly double _密度;

        private 証拠較正器(PairedDistanceModel? p_モデル, double p_密度)
        {
            this._モデル = p_モデル;
            this._密度 = p_密度;
        }

        /// <summary>モデルが構築できたか。false の場合、呼び出し側は生カウント方式に自分でフォールバックする。</summary>
        public bool A_使えるか => this._モデル is not null;

        /// <summary>
        /// 同一 unitig 内標本とリード長から較正器を作る。標本が無い・リード長が
        /// 不明・期待位置数の合計が 0(すべての unitig がフラグメント長より
        /// 短い等)のいずれかならモデルは使えないものとして返す。
        /// </summary>
        public static 証拠較正器 Get_較正器(
            IReadOnlyList<int> p_同一ユニティグ標本,
            int? p_リード長,
            IEnumerable<long> p_ユニティグ長一覧)
        {
            if (p_リード長 is not { } l_リード長 || p_同一ユニティグ標本.Count == 0)
            {
                return new 証拠較正器(null, 0);
            }

            var l_モデル = new PairedDistanceModel(p_同一ユニティグ標本, l_リード長);
            if (!l_モデル.A_使えるか)
            {
                return new 証拠較正器(null, 0);
            }

            double l_期待位置数合計 = 0;
            foreach (var l_長さ in p_ユニティグ長一覧)
            {
                l_期待位置数合計 += l_モデル.Get_期待位置数_単一(l_長さ);
            }
            if (l_期待位置数合計 <= 0)
            {
                return new 証拠較正器(null, 0);
            }

            var l_密度 = p_同一ユニティグ標本.Count / l_期待位置数合計;
            return l_密度 > 0 ? new 証拠較正器(l_モデル, l_密度) : new 証拠較正器(null, 0);
        }

        /// <summary>
        /// 観測本数を期待本数(密度 x 期待位置数)で割った比。
        /// モデルが使えない場合は 0 を返す(=証拠なしとして扱う)。
        /// 生カウントへのフォールバックが必要な呼び出し側は
        /// A_使えるか を先に見て自分で分岐すること。
        /// </summary>
        public double Get_正規化済み支持(ulong p_観測本数, long p_長さ1, long p_長さ2, int p_ギャップ長)
        {
            if (this._モデル is not { } l_モデル)
            {
                return 0;
            }
            var l_期待 = this._密度 * l_モデル.Get_期待位置数(p_長さ1, p_長さ2, p_ギャップ長);
            return l_期待 > 0 ? p_観測本数 / l_期待 : 0;
        }
    }
}
