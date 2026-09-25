using System.Buffers.Binary;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 予算を小さくして何度も吐き出させても、併合した結果が素朴に数えた結果と一致することの検証
    /// </summary>
    [Collection(中間データ置き場の集まり.名前)]
    public class CountingDB併合Tests : IDisposable
    {
        private readonly string _作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_countingdb_merge_" + Guid.NewGuid().ToString("N"));
        private readonly Parameters _元の設定 = ConfigurationManager.A_実行時引数;

        public CountingDB併合Tests()
        {
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        public void Dispose()
        {
            ConfigurationManager.A_実行時引数 = this._元の設定;
            Directory.Delete(this._作業ディレクトリ, true);
        }

        [Theory]
        [InlineData(21)]
        [InlineData(45)]
        [InlineData(99)]
        public void 何度吐き出しても素朴な数え上げと一致する(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1, A_メモリ予算 = "1M" };
            var l_乱数 = new Random(p_k長);
            var l_期待 = new Dictionary<(UInt128, UInt128), ulong>();
            var l_ビット数 = 2 * p_k長;
            using var l_数え = new CountingDB(this._作業ディレクトリ);
            for (var i = 0; i < 60_000; i++)
            {
                // 同じ値が何度も出るように種類を絞る
                var l_種 = (ulong)l_乱数.Next(20_000);
                var l_下位 = (UInt128)(l_種 * 0x9E37_79B9_7F4A_7C15UL) | ((UInt128)l_種 << 64);
                var l_上位 = l_ビット数 > 128 ? (UInt128)l_種 : 0;
                if (l_ビット数 < 128)
                {
                    l_下位 &= ((UInt128)1 << l_ビット数) - 1;
                }
                if (l_ビット数 > 128)
                {
                    l_上位 &= ((UInt128)1 << (l_ビット数 - 128)) - 1;
                }
                l_数え.V_登録_値((l_上位, l_下位));
                l_期待[(l_上位, l_下位)] = l_期待.GetValueOrDefault((l_上位, l_下位)) + 1;
            }

            var l_パス = l_数え.Get_統合ファイル();
            var l_パック長 = (p_k長 + 3) / 4;
            var l_中身 = File.ReadAllBytes(l_パス);
            var l_エントリ長 = l_パック長 + sizeof(ulong);
            Assert.Equal(0, l_中身.Length % l_エントリ長);
            var l_件数 = l_中身.Length / l_エントリ長;
            Assert.Equal(l_期待.Count, l_件数);

            var l_期待の並び = l_期待.OrderBy(static x => x.Key).ToList();
            var l_バイト列 = new byte[l_パック長];
            for (var i = 0; i < l_件数; i++)
            {
                CountingDB.V_変換_パック済みバイト列(l_期待の並び[i].Key, p_k長, l_バイト列);
                Assert.True(l_中身.AsSpan(i * l_エントリ長, l_パック長).SequenceEqual(l_バイト列), $"{i} 件目の k-mer が違う");
                Assert.Equal(l_期待の並び[i].Value, BinaryPrimitives.ReadUInt64LittleEndian(l_中身.AsSpan((i * l_エントリ長) + l_パック長)));
            }
        }
    }
}
