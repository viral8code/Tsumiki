using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 配列の端を、リードの続きの合意で伸ばす処理 (ReadEndExtender) の検証
    /// </summary>
    public class ReadEndExtenderTests
    {
        #region 定数

        /// <summary>
        /// 項目 S
        /// </summary>
        private const string C_項目_S = "S";

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 150;

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 1000;

        /// <summary>
        /// 配列 S の右端の位置 (G[0..700])
        /// </summary>
        private const int C_右端 = 700;

        /// <summary>
        /// 配列 S の左端の位置 (G[300..1000])
        /// </summary>
        private const int C_左端 = 300;

        /// <summary>
        /// 伸ばさないときに比べる、配列の長さの下限 (伸ばした配列は元より長い)
        /// </summary>
        private const int C_元の長さ = 700;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 右端の続きがリードで合意されるとき、配列は G の先頭からの部分列で長くなり、伸ばした端は 1 か所
        /// </summary>
        [Fact]
        public void Get_伸ばした配列群_右端の続きがある_Gの先頭からの部分列で長くなる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)> { (C_項目_S, l_ゲノム[..C_右端]) };
            var (l_結果, l_伸ばした端の数, l_伸ばした塩基数) = Get_実行(l_配列群, Get_右端のリード群(l_ゲノム));
            Assert.Single(l_結果);
            Assert.Equal(C_項目_S, l_結果[0].A_ID);
            var l_配列 = l_結果[0].A_配列;
            Assert.True(l_配列.Length > C_元の長さ);
            Assert.Equal(l_ゲノム[..l_配列.Length], l_配列);
            Assert.True(l_伸ばした端の数 >= 1);
            Assert.Equal(1, l_伸ばした端の数);
            Assert.Equal(l_配列.Length - C_右端, l_伸ばした塩基数);
        }

        /// <summary>
        /// 端の末尾 10 塩基が読み違いでも、その 10 塩基は置き直され、配列は G の先頭からの部分列になる
        /// </summary>
        [Fact]
        public void Get_伸ばした配列群_端の読み違いあり_読み違いが消えGの先頭からの部分列になる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_読み違い = Get_末尾を別の塩基にした(l_ゲノム[..C_右端], 10);
            var l_配列群 = new List<(string A_ID, string A_配列)> { (C_項目_S, l_読み違い) };
            var (l_結果, l_伸ばした端の数, _) = Get_実行(l_配列群, Get_右端のリード群(l_ゲノム));
            Assert.Single(l_結果);
            var l_配列 = l_結果[0].A_配列;
            Assert.True(l_配列.Length > C_元の長さ);
            Assert.Equal(l_ゲノム[..l_配列.Length], l_配列);
            Assert.Equal(1, l_伸ばした端の数);
        }

        /// <summary>
        /// リードが 3 本だけなら、合意がとれず配列は変わらない
        /// </summary>
        [Fact]
        public void Get_伸ばした配列群_リードが3本_配列は変わらない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)> { (C_項目_S, l_ゲノム[..C_右端]) };
            var l_リード群 = Get_右端のリード群(l_ゲノム).Take(3).ToList();
            var (l_結果, l_伸ばした端の数, l_伸ばした塩基数) = Get_実行(l_配列群, l_リード群);
            Assert.Single(l_結果);
            Assert.Equal(l_ゲノム[..C_右端], l_結果[0].A_配列);
            Assert.Equal(0, l_伸ばした端の数);
            Assert.Equal(0, l_伸ばした塩基数);
        }

        /// <summary>
        /// 錨の先の続きが半分ずつ割れるとき、合意は 30 bp 以下で止まり、配列は変わらない
        /// </summary>
        [Fact]
        public void Get_伸ばした配列群_続きが割れる_合意が止まり配列は変わらない()
        {
            var l_ゲノム = Get_ゲノム();
            var l_別のゲノム = Get_続きが別のゲノム(l_ゲノム);
            var l_配列 = l_ゲノム[..C_右端];
            var l_配列群 = new List<(string A_ID, string A_配列)> { (C_項目_S, l_配列) };
            var l_リード群 = Get_交互のリード群(l_ゲノム, l_別のゲノム);
            var (l_結果, l_伸ばした端の数, l_伸ばした塩基数) = Get_実行(l_配列群, l_リード群);
            Assert.Single(l_結果);
            Assert.Equal(l_配列, l_結果[0].A_配列);
            Assert.Equal(0, l_伸ばした端の数);
            Assert.Equal(0, l_伸ばした塩基数);
        }

        /// <summary>
        /// 左端の続きがリードで合意されるとき、配列は G の末尾までの部分列で長くなり、伸ばした端は 1 か所
        /// </summary>
        [Fact]
        public void Get_伸ばした配列群_左端の続きがある_Gの末尾までの部分列で長くなる()
        {
            var l_ゲノム = Get_ゲノム();
            var l_配列群 = new List<(string A_ID, string A_配列)> { (C_項目_S, l_ゲノム[C_左端..]) };
            var (l_結果, l_伸ばした端の数, _) = Get_実行(l_配列群, Get_左端のリード群(l_ゲノム));
            Assert.Single(l_結果);
            var l_配列 = l_結果[0].A_配列;
            Assert.True(l_配列.Length > C_元の長さ);
            Assert.Equal(l_ゲノム[^l_配列.Length..], l_配列);
            Assert.Equal(1, l_伸ばした端の数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を固定した乱数で、一意な 1,000 塩基のゲノムを作る
        /// </summary>
        /// <returns>ゲノム</returns>
        private static string Get_ゲノム()
        {
            return Get_ランダムな塩基(new Random(20240701), C_ゲノム長);
        }

        /// <summary>
        /// 指定した乱数で、A/C/G/T からなる塩基列を作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <param name="p_長さ">長さ</param>
        /// <returns>塩基列</returns>
        private static string Get_ランダムな塩基(Random p_乱数, int p_長さ)
        {
            const string l_塩基 = Consts.塩基文字;
            var l_配列 = new char[p_長さ];
            for (var i = 0; i < p_長さ; i++)
            {
                l_配列[i] = l_塩基[p_乱数.Next(4)];
            }

            return new string(l_配列);
        }

        /// <summary>
        /// 末尾の指定した長さの塩基を、別の塩基に置き換える
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_長さ">置き換える末尾の長さ</param>
        /// <returns>置き換えた配列</returns>
        private static string Get_末尾を別の塩基にした(string p_配列, int p_長さ)
        {
            const string l_塩基 = Consts.塩基文字;
            var l_配列 = p_配列.ToCharArray();
            for (var i = l_配列.Length - p_長さ; i < l_配列.Length; i++)
            {
                l_配列[i] = l_塩基[(l_塩基.IndexOf(l_配列[i]) + 1) % l_塩基.Length];
            }

            return new string(l_配列);
        }

        /// <summary>
        /// 位置 670 の塩基だけを変え、以降は乱数で作る別のゲノムを作る (前半はゲノムと同じ)
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <returns>続きが割れた別のゲノム</returns>
        private static string Get_続きが別のゲノム(string p_ゲノム)
        {
            const string l_塩基 = Consts.塩基文字;
            var l_乱数 = new Random(20240702);
            var l_位置の塩基 = l_塩基[(l_塩基.IndexOf(p_ゲノム[670]) + 1) % l_塩基.Length];
            return p_ゲノム[..670] + l_位置の塩基 + Get_ランダムな塩基(l_乱数, p_ゲノム.Length - 671);
        }

        /// <summary>
        /// 右端の錨を含む、G[560..] から始まる長さ 150 のリードを 5 塩基おきに 20 本作る<br/>
        /// 奇数番目は逆相補にする
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <returns>リード群</returns>
        private static List<string> Get_右端のリード群(string p_ゲノム)
        {
            var l_リード群 = new List<string>();
            for (var i = 0; i < 20; i++)
            {
                var l_始め = 560 + 5 * i;
                var l_リード = p_ゲノム.Substring(l_始め, C_リード長);
                l_リード群.Add(i % 2 == 1 ? Util.V_逆相補_曖昧塩基あり(l_リード) : l_リード);
            }

            return l_リード群;
        }

        /// <summary>
        /// 左端の錨を含む、G[140..290] から始まる長さ 150 のリードを 5 塩基おきに 31 本作る<br/>
        /// 奇数番目は逆相補にする
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <returns>リード群</returns>
        private static List<string> Get_左端のリード群(string p_ゲノム)
        {
            var l_リード群 = new List<string>();
            for (var i = 0; i <= 30; i++)
            {
                var l_始め = 140 + 5 * i;
                var l_リード = p_ゲノム.Substring(l_始め, C_リード長);
                l_リード群.Add(i % 2 == 1 ? Util.V_逆相補_曖昧塩基あり(l_リード) : l_リード);
            }

            return l_リード群;
        }

        /// <summary>
        /// 右端の錨を含むリードを、偶数番目は G から、奇数番目は別のゲノム H から作る<br/>
        /// 開始位置は 560 から 5 塩基おき、逆相補は 4 本ごとに後ろの 2 本
        /// </summary>
        /// <param name="p_ゲノム">ゲノム G</param>
        /// <param name="p_別のゲノム">続きが割れた別のゲノム H</param>
        /// <returns>リード群</returns>
        private static List<string> Get_交互のリード群(string p_ゲノム, string p_別のゲノム)
        {
            var l_リード群 = new List<string>();
            for (var i = 0; i < 20; i++)
            {
                var l_元 = i % 2 == 0 ? p_ゲノム : p_別のゲノム;
                var l_始め = 560 + 5 * i;
                var l_リード = l_元.Substring(l_始め, C_リード長);
                l_リード群.Add(i % 4 >= 2 ? Util.V_逆相補_曖昧塩基あり(l_リード) : l_リード);
            }

            return l_リード群;
        }

        /// <summary>
        /// 伸ばす処理を、リード群の索引を作って実行する
        /// </summary>
        /// <param name="p_配列群">配列群</param>
        /// <param name="p_リード群">リード群</param>
        /// <returns>伸ばした配列群、伸ばした端の数、伸ばした塩基数</returns>
        private static (List<(string A_ID, string A_配列)> A_結果, int A_伸ばした端の数, int A_伸ばした塩基数) Get_実行(List<(string A_ID, string A_配列)> p_配列群, List<string> p_リード群)
        {
            var l_索引 = ReadMinimizerIndex.V_構築(() => p_リード群);
            var l_結果 = ReadEndExtender.Get_伸ばした配列群(p_配列群, l_索引, C_リード長, 0, out var l_伸ばした端の数, out var l_伸ばした塩基数);
            return (l_結果, l_伸ばした端の数, l_伸ばした塩基数);
        }

        #endregion
    }
}
