namespace Tsumiki.Core.Evidence
{
    /// <summary>
    /// ペアエンドの隣接証拠を、生の観測本数ではなく期待本数との比で測るための
    /// 較正器<br/>
    /// ある辺で観測される本数は、その辺の両側の unitig 長と
    /// フラグメント長分布に応じて期待値そのものが桁で変わるため、
    /// 観測本数をそのまま固定閾値と比べると、短い辺には厳しすぎ、
    /// 長い辺には緩すぎる基準になる<br/>
    /// Scaffolder(スキャフォールド辺の選択)・ContigMaker(分岐選択)・
    /// BeamSearchExtender(先読みスコア)が同じ較正を共有するために括り出した
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

        /// <summary>
        /// 支持本数を 0〜1 の確信度へ潰すときの時定数<br/>
        /// この本数あたりで
        /// 6 割強に達し、以降は増やしてもほとんど動かなくなる
        /// </summary>
        public const double 飽和の時定数 = 3.0;

        /// <summary>
        /// 独立な支持本数を 0〜1 の確信度に変換する<br/>
        /// 本数をそのまま足し合わせると、ほぼ同じ条件のペアが 1000 本あるだけで
        /// 種類の違う証拠 1 本を完全に押し流してしまう<br/>
        /// 頭打ちにすることで、
        /// 「同じ話を何度も聞いた」ことが「別の裏付けを得た」ことに化けるのを防ぐ
        /// </summary>
        public static double Get_飽和支持(double p_独立支持数)
        {
            return p_独立支持数 <= 0 ? 0 : 1 - Math.Exp(-p_独立支持数 / 飽和の時定数);
        }

        /// <summary>
        /// 観測された既知長の一覧から、独立とみなせる支持本数を数える<br/>
        /// 同じ辺に対して全く同じ距離を示す観測は、PCR 重複か同一断片の
        /// 読み直しである可能性が高く、別々の分子から得た裏付けとは言えない<br/>
        /// 相異なる距離の個数を独立な証拠の数とみなす
        /// </summary>
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

        /// <summary>
        /// モデルが構築できたか<br/>
        /// false の場合、呼び出し側は生カウント方式に自分でフォールバックする
        /// </summary>
        public bool A_使えるか => this._モデル is not null;

        /// <summary>
        /// 同一 unitig 内標本とリード長から較正器を作る<br/>
        /// 標本が無い・リード長が
        /// 不明・期待位置数の合計が 0(すべての unitig がフラグメント長より
        /// 短い等)のいずれかならモデルは使えないものとして返す
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
        /// 観測本数を期待本数(密度 x 期待位置数)で割った比<br/>
        /// モデルが使えない場合は 0 を返す(=証拠なしとして扱う)<br/>
        /// 生カウントへのフォールバックが必要な呼び出し側は
        /// A_使えるか を先に見て自分で分岐すること
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
