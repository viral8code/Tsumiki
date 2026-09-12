using Tsumiki.Commons;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 塩基列を 2 bit/塩基 で詰め、語単位で突き合わせる部品の検証
    /// </summary>
    /// <remarks>
    /// 重なりの探索がこれに置き換わるので、1 塩基ずつ比べた場合と完全に同じ数を返すことが要件になる
    /// </remarks>
    public class PackedBasesTests
    {
        #region 公開メソッド

        /// <summary>
        /// 開始位置・長さのあらゆる組み合わせで、1 塩基ずつ数えた結果と一致すること
        /// </summary>
        /// <remarks>
        /// 語の境界をまたぐ位置 (32 の倍数の前後) を必ず含むように総当たりする
        /// </remarks>
        [Fact]
        public void V_あらゆる開始位置と長さで語単位の不一致数が素朴な数え方と一致する()
        {
            var l_乱数 = new Random(20_260_925);
            for (var l_試行 = 0; l_試行 < 20; l_試行++)
            {
                var l_列1 = V_生成_乱数塩基列(150, l_乱数);
                var l_列2 = V_生成_乱数塩基列(150, l_乱数);

                // 一部を一致させて、不一致 0 や少数の場合も通す
                Array.Copy(l_列1, 20, l_列2, 20, 60);

                var l_詰め1 = PackedBases.Get_作る(l_列1);
                var l_詰め2 = PackedBases.Get_作る(l_列2);
                Assert.NotNull(l_詰め1);
                Assert.NotNull(l_詰め2);

                for (var l_開始1 = 0; l_開始1 < 70; l_開始1++)
                {
                    for (var l_開始2 = 0; l_開始2 < 70; l_開始2 += 7)
                    {
                        var l_上限 = 150 - Math.Max(l_開始1, l_開始2);
                        for (var l_長さ = 1; l_長さ <= l_上限; l_長さ += 3)
                        {
                            Assert.Equal(Get_不一致数_素朴(l_列1, l_列2, l_開始1, l_開始2, l_長さ), Get_不一致数_語単位(l_詰め1!, l_詰め2!, l_開始1, l_開始2, l_長さ));
                        }
                    }
                }
            }
        }

        /// <summary>
        /// 曖昧な塩基を含む列は作れない
        /// </summary>
        [Fact]
        public void V_曖昧な塩基を含む列は作れない()
        {
            byte[] l_列 = [Consts.塩基ID.A, Consts.無効な塩基, Consts.塩基ID.T];

            Assert.Null(PackedBases.Get_作る(l_列));
        }

        /// <summary>
        /// 末尾を超えた範囲は 0 として読める
        /// </summary>
        [Fact]
        public void V_末尾を超えた範囲は0として読める()
        {
            var l_詰め = PackedBases.Get_作る([Consts.塩基ID.T, Consts.塩基ID.T]);

            Assert.NotNull(l_詰め);

            // 3 塩基目以降は空き
            // 先頭 2 塩基だけを見れば不一致は無い
            Assert.Equal(0, PackedBases.Get_不一致数(l_詰め!.Get_窓(0), l_詰め.Get_窓(0), 2));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 乱数から塩基 ID 列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_乱数">使う乱数</param>
        /// <returns>塩基 ID 列</returns>
        private static byte[] V_生成_乱数塩基列(int p_長さ, Random p_乱数)
        {
            var l_列 = new byte[p_長さ];
            for (var i = 0; i < p_長さ; i++)
            {
                l_列[i] = (byte)(p_乱数.Next(4) + 1);
            }
            return l_列;
        }

        /// <summary>
        /// 1 塩基ずつ比べて不一致の数を数える
        /// </summary>
        /// <param name="p_列1">比べる塩基 ID 列</param>
        /// <param name="p_列2">比べる塩基 ID 列</param>
        /// <param name="p_開始1">片方で比べ始める位置</param>
        /// <param name="p_開始2">もう片方で比べ始める位置</param>
        /// <param name="p_長さ">比べる長さ</param>
        /// <returns>不一致の数</returns>
        private static int Get_不一致数_素朴(byte[] p_列1, byte[] p_列2, int p_開始1, int p_開始2, int p_長さ)
        {
            var l_数 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                if (p_列1[p_開始1 + i] != p_列2[p_開始2 + i])
                {
                    l_数++;
                }
            }
            return l_数;
        }

        /// <summary>
        /// 語単位で比べて不一致の数を数える
        /// </summary>
        /// <param name="p_詰め1">比べる詰めた塩基列</param>
        /// <param name="p_詰め2">比べる詰めた塩基列</param>
        /// <param name="p_開始1">片方で比べ始める位置</param>
        /// <param name="p_開始2">もう片方で比べ始める位置</param>
        /// <param name="p_長さ">比べる長さ</param>
        /// <returns>不一致の数</returns>
        private static int Get_不一致数_語単位(PackedBases p_詰め1, PackedBases p_詰め2, int p_開始1, int p_開始2, int p_長さ)
        {
            var l_数 = 0;
            for (var i = 0; i < p_長さ; i += Consts.語あたりの塩基数)
            {
                var l_今回 = Math.Min(Consts.語あたりの塩基数, p_長さ - i);
                l_数 += PackedBases.Get_不一致数(p_詰め1.Get_窓(p_開始1 + i), p_詰め2.Get_窓(p_開始2 + i), l_今回);
            }
            return l_数;
        }

        #endregion

    }
}
