using System.Numerics;
using Tsumiki.Commons;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 塩基列を 2 bit / 塩基 で 64 bit 語に詰めたもの
    /// </summary>
    /// <remarks>
    /// 重なりの探索はずらしながら 2 本を突き合わせる形をしており、
    /// 1 塩基ずつ比べると 1 ペアあたり数万回の比較になる<br/>
    /// 語単位で XOR して不一致のレーンを数えれば、同じことが 32 塩基まとめて片付く<br/>
    /// 曖昧塩基は 2 bit に収まらないため、含む配列は作れない
    /// </remarks>
    internal sealed class PackedBases
    {
        /// <summary>
        /// 1 語に詰まる塩基数
        /// </summary>
        public const int 語あたりの塩基数 = 32;

        /// <summary>
        /// 2 bit レーンの下位ビットだけを立てたマスク
        /// </summary>
        private const ulong 下位ビット = 0x5555555555555555UL;

        /// <summary>
        /// 詰めた語の並び
        /// </summary>
        private readonly ulong[] _語;

        /// <summary>
        /// 詰めた塩基の数
        /// </summary>
        public int A_長さ { get; }

        /// <summary>
        /// コンストラクタ
        /// </summary>
        /// <param name="p_語">詰めた語の並び</param>
        /// <param name="p_長さ">詰めた塩基の数</param>
        private PackedBases(ulong[] p_語, int p_長さ)
        {
            this._語 = p_語;
            this.A_長さ = p_長さ;
        }

        /// <summary>
        /// 塩基 ID 列を詰める
        /// </summary>
        /// <remarks>
        /// 窓の取り出しで次の語を無条件に読めるよう、末尾に空き語を 1 つ足す
        /// </remarks>
        /// <param name="p_塩基列">塩基 ID 列 (A = 1 .. T = 4)</param>
        /// <returns>詰めた結果、曖昧塩基を含む場合は null</returns>
        public static PackedBases? Get_作る(ReadOnlySpan<byte> p_塩基列)
        {
            var l_語 = new ulong[(p_塩基列.Length / 語あたりの塩基数) + 2];
            for (var i = 0; i < p_塩基列.Length; i++)
            {
                var l_塩基ID = p_塩基列[i];
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    return null;
                }
                var l_語番号 = i / 語あたりの塩基数;
                var l_ずらし = 62 - (2 * (i % 語あたりの塩基数));
                l_語[l_語番号] |= (ulong)(l_塩基ID - 1) << l_ずらし;
            }
            return new PackedBases(l_語, p_塩基列.Length);
        }

        /// <summary>
        /// 指定位置から 32 塩基ぶんを 1 語として取り出す
        /// </summary>
        /// <param name="p_位置">取り出しを始める塩基の位置</param>
        /// <returns>取り出した語、配列の末尾を越える分は 0</returns>
        public ulong Get_窓(int p_位置)
        {
            var l_語番号 = p_位置 / 語あたりの塩基数;
            var l_ずらし = 2 * (p_位置 % 語あたりの塩基数);
            var l_先頭 = this._語[l_語番号];
            return l_ずらし == 0
                ? l_先頭
                : (l_先頭 << l_ずらし) | (this._語[l_語番号 + 1] >> (64 - l_ずらし));
        }

        /// <summary>
        /// 2 つの窓の先頭から数えて、一致しない塩基の数を返す
        /// </summary>
        /// <param name="p_窓1">比べる語</param>
        /// <param name="p_窓2">比べる語</param>
        /// <param name="p_塩基数">先頭から何塩基を比べるか</param>
        /// <returns>一致しない塩基の数</returns>
        public static int Get_不一致数(ulong p_窓1, ulong p_窓2, int p_塩基数)
        {
            var l_差 = p_窓1 ^ p_窓2;

            // 2 bit のどちらかが立っていれば不一致なので、レーンごとに 1 ビットへ畳む
            var l_レーン = (l_差 | (l_差 >> 1)) & 下位ビット;

            if (p_塩基数 < 語あたりの塩基数)
            {
                l_レーン &= ulong.MaxValue << (64 - (2 * p_塩基数));
            }
            return BitOperations.PopCount(l_レーン);
        }
    }
}
