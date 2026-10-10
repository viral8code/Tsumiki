using Tsumiki.Commons;
using Tsumiki.Cores.Mapping;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Mapping;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core.Mapping
{
    /// <summary>
    /// 配置候補の選び方と並び (種の多い順、同数なら最初に現れた順、上位 16 件) が、OrderByDescending と Take による参照実装と一致することの検証
    /// </summary>
    public class ReadMapperCandidateOrderTests
    {
        #region 定数

        /// <summary>
        /// ReadMapper の候補上限 (非公開のため、参照実装でも同じ値を使う)
        /// </summary>
        private const int C_参照実装の候補上限 = 16;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// テストに必要な共有状態を初期化する
        /// </summary>
        public ReadMapperCandidateOrderTests()
        {
            ConfigurationManager.A_実行時引数 = new Parameters();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 反復を含む参照に対するリード (両鎖、誤り・挿入・欠失あり) で、候補の並びと内容が参照実装と一致する
        /// </summary>
        /// <param name="p_種">乱数の種</param>
        [Theory]
        [InlineData(1)]
        [InlineData(2)]
        [InlineData(3)]
        public void Get_配置候補群_候補の並びは参照実装と同じになる(int p_種)
        {
            var l_乱数 = new Random(p_種);
            var l_反復 = Get_ランダム配列(l_乱数, 100);
            var l_参照 = string.Concat(Enumerable.Range(0, 20).SelectMany(_ => new[] { Get_ランダム配列(l_乱数, 150), l_反復 }));
            var l_マッパー = new ReadMapper([l_参照]);
            var l_候補超え = 0;
            for (var i = 0; i < 300; i++)
            {
                var l_リード = Get_リード(l_乱数, l_参照);
                if (l_マッパー.Get_候補数(l_リード).Count > C_参照実装の候補上限)
                {
                    l_候補超え++;
                }

                var l_実際 = l_マッパー.Get_配置候補群(l_リード);
                var l_期待 = Get_参照候補群(l_マッパー, l_リード);
                Assert.Equal(l_期待.Count, l_実際.Count);
                for (var j = 0; j < l_期待.Count; j++)
                {
                    Assert.Equal(l_期待[j].A_配列番号, l_実際[j].A_配列番号);
                    Assert.Equal(l_期待[j].A_Is逆鎖, l_実際[j].A_Is逆鎖);
                    Assert.Equal(l_期待[j].A_スコア, l_実際[j].A_スコア);
                    Assert.Equal(l_期待[j].A_整列位置群, l_実際[j].A_整列位置群);
                }
            }

            Assert.True(l_候補超え > 0);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 参照実装: OrderByDescending と Take で候補を選び、非公開の整列で配置にする
        /// </summary>
        /// <param name="p_マッパー">配置する ReadMapper</param>
        /// <param name="p_リード">リード</param>
        /// <returns>変更前と同じ方法で作った配置候補</returns>
        private static List<リード配置> Get_参照候補群(ReadMapper p_マッパー, string p_リード)
        {
            var l_順 = p_マッパー.Get_候補数(p_リード).OrderByDescending(x => x.Value).Take(C_参照実装の候補上限).Select(x => x.Key).ToList();
            var l_逆相補 = Util.V_逆相補_曖昧塩基あり(p_リード);
            var l_結果 = new List<リード配置>();
            foreach (var l_候補 in l_順)
            {
                var l_照合 = l_候補.A_Is逆鎖 ? l_逆相補 : p_リード;
                l_結果.Add(p_マッパー.Get_整列(p_リード, l_照合, l_候補));
            }

            return l_結果;
        }

        /// <summary>
        /// 乱数で ACGT の配列を作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <param name="p_長さ">長さ</param>
        /// <returns>配列</returns>
        private static string Get_ランダム配列(Random p_乱数, int p_長さ)
        {
            return new string([..Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[p_乱数.Next(4)])]);
        }

        /// <summary>
        /// 参照の一部から、誤り・挿入・欠失を入れ、半分は逆相補にしたリードを作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <param name="p_参照">参照</param>
        /// <returns>リード</returns>
        private static string Get_リード(Random p_乱数, string p_参照)
        {
            var l_長さ = p_乱数.Next(80, 201);
            var l_開始 = p_乱数.Next(0, p_参照.Length - l_長さ);
            var l_元 = p_参照.Substring(l_開始, l_長さ);
            var l_組 = new List<char>();
            foreach (var l_文字 in l_元)
            {
                var l_確率 = p_乱数.NextDouble();
                if (l_確率 < 0.03D)
                {
                    l_組.Add(Consts.塩基文字[p_乱数.Next(4)]);
                }
                else if (l_確率 < 0.04D)
                {
                    continue;
                }
                else if (l_確率 < 0.05D)
                {
                    l_組.Add(Consts.塩基文字[p_乱数.Next(4)]);
                    l_組.Add(l_文字);
                }
                else
                {
                    l_組.Add(l_文字);
                }
            }

            var l_リード = new string([..l_組]);
            return p_乱数.Next(2) == 0 ? l_リード : Util.V_逆相補_曖昧塩基あり(l_リード);
        }

        #endregion
    }
}
