using System.Reflection;
using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Mapping;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Mapping;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core.Mapping
{
    /// <summary>
    /// 帯つき整列の最適化版 (Get_整列) が、変更前の参照実装 (Get_整列_参照) と全候補で同じ配置・スコア・信頼度・整列位置を返すことの検証
    /// </summary>
    public class ReadMapperAlignmentEquivalenceTests
    {
        #region 定数

        /// <summary>
        /// 1 種あたりのリードの本数
        /// </summary>
        private const int C_リード数 = 2000;

        /// <summary>
        /// 帯域幅 (端の判定に使う。ReadMapper の値と合わせる)
        /// </summary>
        private const int C_帯域幅 = 24;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// テストに必要な共有状態を初期化する
        /// </summary>
        public ReadMapperAlignmentEquivalenceTests()
        {
            ConfigurationManager.A_実行時引数 = new Parameters();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 乱数の種ごとに、反復と端を含む参照に対するリード (両鎖、誤り・挿入・欠失・N、端からはみ出すもの、短いもの) の全候補で、最適化版と参照実装が一致する
        /// </summary>
        /// <param name="p_種">乱数の種</param>
        [Theory]
        [InlineData(1)]
        [InlineData(2)]
        [InlineData(3)]
        [InlineData(4)]
        [InlineData(5)]
        public void Get_整列_全候補で参照実装と同じ結果になる(int p_種)
        {
            var l_乱数 = new Random(p_種);
            var l_参照群 = Get_参照群(l_乱数);
            var l_マッパー = new ReadMapper(l_参照群);
            var l_整列 = typeof(ReadMapper).GetMethod("Get_整列", BindingFlags.Instance | BindingFlags.NonPublic)
                ?? throw new InvalidOperationException("Get_整列 が見つからない");
            var l_整列_参照 = typeof(ReadMapper).GetMethod("Get_整列_参照", BindingFlags.Instance | BindingFlags.NonPublic)
                ?? throw new InvalidOperationException("Get_整列_参照 が見つからない");

            var l_比較数 = 0;
            var l_端の候補数 = 0;
            var l_得点が正の数 = 0;
            for (var n = 0; n < C_リード数; n++)
            {
                var l_リード = Get_リード(l_乱数, l_参照群);
                foreach (var l_候補 in l_マッパー.Get_候補数(l_リード).Keys)
                {
                    var l_照合リード = l_候補.A_Is逆鎖 ? Util.V_逆相補_曖昧塩基あり(l_リード) : l_リード;
                    var l_実際 = (リード配置)l_整列.Invoke(l_マッパー, new object[] { l_リード, l_照合リード, l_候補 })!;
                    var l_期待 = (リード配置)l_整列_参照.Invoke(l_マッパー, new object[] { l_リード, l_照合リード, l_候補 })!;

                    Assert.Equal(l_期待.A_配列番号, l_実際.A_配列番号);
                    Assert.Equal(l_期待.A_Is逆鎖, l_実際.A_Is逆鎖);
                    Assert.Equal(l_期待.A_スコア, l_実際.A_スコア);
                    Assert.Equal(l_期待.A_信頼度, l_実際.A_信頼度);
                    Assert.Equal(l_期待.A_整列位置群.Count, l_実際.A_整列位置群.Count);
                    for (var k = 0; k < l_期待.A_整列位置群.Count; k++)
                    {
                        Assert.Equal(l_期待.A_整列位置群[k], l_実際.A_整列位置群[k]);
                    }

                    var l_参照長 = l_参照群[l_候補.A_配列番号].Length;
                    if (l_候補.A_対角線 < C_帯域幅 || l_候補.A_対角線 + l_照合リード.Length + C_帯域幅 > l_参照長)
                    {
                        l_端の候補数++;
                    }

                    if (l_期待.A_スコア > 0)
                    {
                        l_得点が正の数++;
                    }

                    l_比較数++;
                }
            }

            // 候補が十分に比べられたことを確認する (端の候補と、得点が正の配置の両方を含む)
            Assert.True(l_比較数 >= C_リード数, $"比べた候補が少なすぎる: {l_比較数}");
            Assert.True(l_端の候補数 > 0, "端に近い候補が 1 つも比べられていない");
            Assert.True(l_得点が正の数 > 0, "得点が正の配置が 1 つも比べられていない");
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 長さ 300〜5,000 の参照配列を 3 本作る。各配列には反復 (同じ断片の複製) を入れ、先頭と末尾にも入れて端に近い候補を作る
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <returns>参照配列群</returns>
        private static List<string> Get_参照群(Random p_乱数)
        {
            List<string> l_参照群 = [];
            for (var k = 0; k < 3; k++)
            {
                var l_長さ = p_乱数.Next(300, 5001);
                var l_配列 = Get_ランダム配列(p_乱数, l_長さ).ToCharArray();
                var l_反復 = Get_ランダム配列(p_乱数, p_乱数.Next(60, 151));
                var l_開始位置群 = new List<int> { 0, l_長さ - l_反復.Length };
                var l_追加数 = p_乱数.Next(1, 4);
                for (var m = 0; m < l_追加数; m++)
                {
                    l_開始位置群.Add(p_乱数.Next(0, l_長さ - l_反復.Length + 1));
                }

                foreach (var l_位置 in l_開始位置群)
                {
                    l_反復.CopyTo(0, l_配列, l_位置, l_反復.Length);
                }

                l_参照群.Add(new string(l_配列));
            }

            return l_参照群;
        }

        /// <summary>
        /// 参照群のどれかから、端からはみ出す位置を含めて切り出し、誤り・挿入・欠失・N を入れ、半分は逆鎖にしたリードを作る
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <param name="p_参照群"></param>
        /// <returns>リード</returns>
        private static string Get_リード(Random p_乱数, IReadOnlyList<string> p_参照群)
        {
            // 5% は参照と無関係なリード (候補が少ない場合も含める)
            if (p_乱数.Next(20) == 0)
            {
                return Get_ランダム配列(p_乱数, p_乱数.Next(30, 251));
            }

            var l_参照 = p_参照群[p_乱数.Next(p_参照群.Count)];
            var l_長さ = p_乱数.Next(30, 251);
            // 端から最大 30 塩基はみ出す
            var l_開始 = p_乱数.Next(-30, l_参照.Length - l_長さ + 31);
            var l_誤り率 = p_乱数.NextDouble() * 0.08;
            var l_塩基 = new StringBuilder(l_長さ + 8);
            for (var p = l_開始; p < l_開始 + l_長さ; p++)
            {
                var l_元 = p >= 0 && p < l_参照.Length ? l_参照[p] : Get_塩基(p_乱数);
                var l_値 = p_乱数.NextDouble();
                if (l_値 < 0.002)
                {
                    l_塩基.Append('N');
                }
                else if (l_値 < 0.002 + l_誤り率)
                {
                    l_塩基.Append(Get_別の塩基(p_乱数, l_元));
                }
                else
                {
                    l_塩基.Append(l_元);
                }
            }

            // 挿入・欠失を 0〜3 個入れる
            var l_挿入欠失数 = p_乱数.Next(4);
            for (var k = 0; k < l_挿入欠失数; k++)
            {
                var l_位置 = p_乱数.Next(l_塩基.Length);
                if (p_乱数.Next(2) == 0)
                {
                    l_塩基.Insert(l_位置, Get_塩基(p_乱数));
                }
                else
                {
                    l_塩基.Remove(l_位置, 1);
                }
            }

            var l_リード = l_塩基.ToString();
            return p_乱数.Next(2) == 0 ? l_リード : Util.V_逆相補_曖昧塩基あり(l_リード);
        }

        /// <summary>
        /// ACGT のランダムな配列を作る
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <param name="p_長さ"></param>
        /// <returns>配列</returns>
        private static string Get_ランダム配列(Random p_乱数, int p_長さ)
        {
            var l_配列 = new char[p_長さ];
            for (var i = 0; i < p_長さ; i++)
            {
                l_配列[i] = Get_塩基(p_乱数);
            }

            return new string(l_配列);
        }

        /// <summary>
        /// ACGT のうち 1 塩基を返す
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <returns>塩基</returns>
        private static char Get_塩基(Random p_乱数)
        {
            return "ACGT"[p_乱数.Next(4)];
        }

        /// <summary>
        /// p_元 と異なる ACGT の塩基を返す
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <param name="p_元"></param>
        /// <returns>塩基</returns>
        private static char Get_別の塩基(Random p_乱数, char p_元)
        {
            char l_塩基;
            do
            {
                l_塩基 = Get_塩基(p_乱数);
            } while (l_塩基 == p_元);

            return l_塩基;
        }

        #endregion
    }
}
