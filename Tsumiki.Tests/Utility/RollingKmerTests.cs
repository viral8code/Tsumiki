using System.Diagnostics;
using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;
using Xunit.Abstractions;

namespace Tsumiki.Tests.Utility
{
    /// <summary>差分更新の正確性と合成データでの処理費用を検証する</summary>
    /// <param name="p_出力">テストログ</param>
    public class RollingKmerTests(ITestOutputHelper p_出力)
    {
        #region 公開メソッド

        /// <summary>境界長と曖昧塩基を含む全窓を従来の厳密キーと比較する</summary>
        /// <param name="p_長さ">窓の長さ</param>
        [Theory]
        [InlineData(1)]
        [InlineData(21)]
        [InlineData(31)]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(63)]
        [InlineData(64)]
        [InlineData(65)]
        [InlineData(73)]
        [InlineData(103)]
        [InlineData(127)]
        [InlineData(128)]
        public void V_全窓を従来キーと比較(int p_長さ)
        {
            var l_配列 = Get_合成配列(900) + "NRYACGT" + Get_合成配列(900);
            var l_窓 = new RollingKmer(p_長さ);
            for (var i = 0; i < l_配列.Length; i++)
            {
                var l_有効 = l_窓.Try追加(l_配列[i], out var l_キー);
                var l_期待 = i >= p_長さ - 1 && !l_配列.Substring(i - p_長さ + 1, p_長さ).Any(Util.Is曖昧塩基);
                Assert.Equal(l_期待, l_有効);
                if (l_期待)
                {
                    Assert.Equal(RepeatRMerVerifier.Get_正準値(l_配列.AsSpan(i - p_長さ + 1, p_長さ)), l_キー);
                }
            }
        }

        /// <summary>種の絞り込みが逆相補の完全一致を落とさない</summary>
        [Fact]
        public void V_種は両鎖で一致する()
        {
            var l_配列 = Get_合成配列(135);
            HashSet<ulong> l_種 = [(ulong)RepeatRMerVerifier.Get_正準値(l_配列.AsSpan(0, 31)).A_下位];
            Assert.True(LocalAssembler.Has種一致("NN" + l_配列, l_種, 31));
            Assert.True(LocalAssembler.Has種一致(Util.V_逆相補(l_配列) + "NN", l_種, 31));
            Assert.False(LocalAssembler.Has種一致(new string('N', 150), l_種, 31));
        }

        /// <summary>最小値の差分計算が端や空の窓でも単純走査と一致する</summary>
        [Fact]
        public void V_最小値列を単純走査と比較()
        {
            var l_乱数 = new Random(812);
            var l_値 = Enumerable.Range(0, 300).Select(_ => l_乱数.Next(-1, 101)).ToArray();
            foreach (var l_幅 in new[] { -1, 0, 1, 9, 43, 115, 400 })
            {
                var l_結果 = new int[320];
                KmerCarryOver.V_計算_最小値列(l_値, l_幅, l_結果, new int[l_値.Length]);
                for (var i = 0; i < l_結果.Length; i++)
                {
                    var l_期待 = l_値.Skip(i).Take(Math.Max(0, l_幅)).DefaultIfEmpty(0).Min();
                    Assert.Equal(l_期待, l_結果[i]);
                }
            }
        }

        /// <summary>合成配列の全窓で従来実装と差分更新の時間を記録する</summary>
        /// <param name="p_長さ">窓の長さ</param>
        [Theory]
        [InlineData(31)]
        [InlineData(73)]
        [InlineData(103)]
        public void V_合成データの速度を記録(int p_長さ)
        {
            var l_配列 = Get_合成配列(500_000);
            UInt128 l_従来合計 = 0;
            var l_時計 = Stopwatch.StartNew();
            for (var i = 0; i + p_長さ <= l_配列.Length; i++)
            {
                l_従来合計 ^= RepeatRMerVerifier.Get_正準値(l_配列.AsSpan(i, p_長さ)).A_下位;
            }
            var l_従来時間 = l_時計.Elapsed.TotalMilliseconds;
            l_時計.Restart();
            UInt128 l_更新合計 = 0;
            var l_窓 = new RollingKmer(p_長さ);
            foreach (var l_塩基 in l_配列)
            {
                if (l_窓.Try追加(l_塩基, out var l_キー))
                {
                    l_更新合計 ^= l_キー.A_下位;
                }
            }
            var l_更新時間 = l_時計.Elapsed.TotalMilliseconds;
            Assert.Equal(l_従来合計, l_更新合計);
            p_出力.WriteLine($"r={p_長さ}: baseline={l_従来時間:F1} ms, rolling={l_更新時間:F1} ms, ratio={l_従来時間 / l_更新時間:F2}");
        }

        /// <summary>無関係な長いリードを種で除外した場合の時間と確保量を記録する</summary>
        [Fact]
        public void V_局所候補の絞り込み速度を記録()
        {
            const int l_k長 = 135;
            var l_配列 = Get_合成配列(200_000);
            var l_アンカー = new string('A', l_k長);
            HashSet<KmerKey> l_索引 = [new KmerKey(l_アンカー.AsSpan()).Get_正規形()];
            var l_確保前 = GC.GetAllocatedBytesForCurrentThread();
            var l_時計 = Stopwatch.StartNew();
            var l_従来一致 = false;
            for (var i = 0; i + l_k長 <= l_配列.Length; i++)
            {
                l_従来一致 |= l_索引.Contains(new KmerKey(l_配列.AsSpan(i, l_k長)).Get_正規形());
            }
            var l_従来時間 = l_時計.Elapsed.TotalMilliseconds;
            var l_従来確保 = GC.GetAllocatedBytesForCurrentThread() - l_確保前;
            HashSet<ulong> l_種 = [0UL];
            l_確保前 = GC.GetAllocatedBytesForCurrentThread();
            l_時計.Restart();
            var l_候補一致 = LocalAssembler.Has種一致(l_配列, l_種, 31);
            var l_候補時間 = l_時計.Elapsed.TotalMilliseconds;
            var l_候補確保 = GC.GetAllocatedBytesForCurrentThread() - l_確保前;
            Assert.False(l_従来一致);
            Assert.False(l_候補一致);
            p_出力.WriteLine($"local negative: baseline={l_従来時間:F1} ms/{l_従来確保} B, seed={l_候補時間:F1} ms/{l_候補確保} B");
        }

        #endregion

        #region 内部メソッド

        /// <summary>機密データを使わない再現可能な塩基列を作る</summary>
        /// <param name="p_長さ">配列長</param>
        /// <returns>合成配列</returns>
        private static string Get_合成配列(int p_長さ)
        {
            var l_乱数 = new Random(912);
            return new string(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]).ToArray());
        }

        #endregion
    }
}
