using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 全リードのペアで端どうしを結び、局所組み立てで間を埋めて配列を繋ぐ処理 (PairBridgeFiller) の検証
    /// </summary>
    public class PairBridgeFillerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki pairbridge tests
        /// </summary>
        private const string C_項目_tsumiki_pairbridge_tests = "tsumiki_pairbridge_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 6000;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 12345;

        /// <summary>
        /// 配列 A の長さ (A = G[0..2800])
        /// </summary>
        private const int C_A長 = 2800;

        /// <summary>
        /// 配列 B の開始 (B = G[2850..6000])
        /// </summary>
        private const int C_B開始 = 2850;

        /// <summary>
        /// 断片の長さ (挿入長)
        /// </summary>
        private const int C_挿入長 = 300;

        /// <summary>
        /// リードの長さ
        /// </summary>
        private const int C_リード長 = 100;

        /// <summary>
        /// 断片の開始の間隔 (密)
        /// </summary>
        private const int C_間隔_密 = 7;

        /// <summary>
        /// 断片の開始の間隔 (疎)
        /// </summary>
        private const int C_間隔_疎 = 200;

        #endregion

        #region フィールド

        /// <summary>
        /// テストごとの作業ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 作業ディレクトリを作り、実行時引数を既定に戻す
        /// </summary>
        public PairBridgeFillerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_pairbridge_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = 21,
                A_スレッド数 = 1
            };
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 作業ディレクトリを消す
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 基本: A と B の間 50 bp を、挿入長 300 のペアで埋めて 1 本に繋ぐ
        /// </summary>
        [Fact]
        public void 基本_間50bpをペアと局所組み立てで埋めて1本に繋ぐ()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_ライブラリ = V_書き出し_ペア("base", l_G, V_断片の開始群(0, 5700, C_間隔_密));

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Single(l_結果);
            Assert.Equal(l_G, l_結果[0].A_配列);
            Assert.Equal(C_挿入長, l_統計.A_挿入長);
        }

        /// <summary>
        /// 向き: B を逆相補にして渡しても 1 本になり、G か G の逆相補と一致する
        /// </summary>
        [Fact]
        public void 向き_Bを逆相補で渡しても1本になりGかその逆相補と一致する()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = Util.V_逆相補_曖昧塩基あり(l_G[C_B開始..]);
            var l_ライブラリ = V_書き出し_ペア("rc", l_G, V_断片の開始群(0, 5700, C_間隔_密));

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out _);

            Assert.Single(l_結果);
            var l_配列 = l_結果[0].A_配列;
            Assert.True(l_配列 == l_G || l_配列 == Util.V_逆相補_曖昧塩基あり(l_G));
        }

        /// <summary>
        /// 当て方の位置: A の順向き 2700..2800 は (頂点 0, 位置 2700)、d は 100
        /// </summary>
        [Fact]
        public void 当て方の位置_順向きの末尾リードは頂点0の位置2700に当たる()
        {
            var (l_頂点, l_表) = V_A単独の索引();
            var l_A = V_A();

            var l_当たり = PairBridgeFiller.Get_当たり(l_A.Substring(2700, C_リード長), l_表, l_頂点);

            V_当たりを確かめる(l_当たり, 0, 2700, l_頂点, 100);
        }

        /// <summary>
        /// 当て方の位置: A の 2700..2800 の逆相補は (頂点 1, 位置 0)、d は 2800
        /// </summary>
        [Fact]
        public void 当て方の位置_末尾の逆相補リードは頂点1の位置0に当たる()
        {
            var (l_頂点, l_表) = V_A単独の索引();
            var l_A = V_A();

            var l_当たり = PairBridgeFiller.Get_当たり(Util.V_逆相補_曖昧塩基あり(l_A.Substring(2700, C_リード長)), l_表, l_頂点);

            V_当たりを確かめる(l_当たり, 1, 0, l_頂点, 2800);
        }

        /// <summary>
        /// 当て方の位置: A の 0..100 の逆相補は (頂点 1, 位置 2700)、頂点 1 の右端まで d は 100
        /// </summary>
        [Fact]
        public void 当て方の位置_先頭の逆相補リードは頂点1の位置2700に当たる()
        {
            var (l_頂点, l_表) = V_A単独の索引();
            var l_A = V_A();

            var l_当たり = PairBridgeFiller.Get_当たり(Util.V_逆相補_曖昧塩基あり(l_A.Substring(0, C_リード長)), l_表, l_頂点);

            V_当たりを確かめる(l_当たり, 1, 2700, l_頂点, 100);
        }

        /// <summary>
        /// 反復の端は繋がない: A の末尾 1,000 bp を含む別の配列 C があると、A の右端の種が捨てられ、繋がらない
        /// </summary>
        [Fact]
        public void 反復の端は繋がない_末尾1000bpが別配列にもあれば繋がらない()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_C = V_生成_ゲノム(500, 7001) + l_A[(C_A長 - 1000)..] + V_生成_ゲノム(500, 7002);
            var l_ライブラリ = V_書き出し_ペア("repeat", l_G, V_断片の開始群(0, 5700, C_間隔_密));

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B), ("c", l_C)], l_ライブラリ, out var l_統計);

            Assert.Equal(3, l_結果.Count);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
        }

        /// <summary>
        /// 票が足りない: ペアの間隔が 200 bp なら A と B を結ぶペアが無く、元の 2 本のまま
        /// </summary>
        [Fact]
        public void 票が足りないと繋がらない_間隔200bpでは元の2本のまま()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_ライブラリ = V_書き出し_ペア("sparse", l_G, V_断片の開始群(0, 5700, C_間隔_疎));

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out _);

            Assert.Equal(2, l_結果.Count);
            Assert.Equal(l_A, l_結果[0].A_配列);
            Assert.Equal(l_B, l_結果[1].A_配列);
        }

        /// <summary>
        /// 埋まらないときは繋がない: 間の配列を含むリードを除くと、票はあっても埋まらず元の 2 本のまま
        /// </summary>
        [Fact]
        public void 埋まらないときは繋がない_間の配列を含むリードが無ければ元の2本のまま()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_開始群 = V_断片の開始群(0, 5700, C_間隔_密)
                .Where(s => !V_重なる(s, s + C_リード長) && !V_重なる(s + 200, s + 300))
                .ToList();
            var l_ライブラリ = V_書き出し_ペア("gap", l_G, l_開始群);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Equal(2, l_結果.Count);
            Assert.Equal(l_A, l_結果[0].A_配列);
            Assert.Equal(l_B, l_結果[1].A_配列);
            Assert.Equal(1, l_統計.A_候補の辺の数);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
        }

        /// <summary>
        /// スレッド数 1 と 4 で、出力が同じになる
        /// </summary>
        [Fact]
        public void スレッド数1と4で出力が同じになる()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_配列群 = new List<(string A_ID, string A_配列)> { ("a", l_G[..C_A長]), ("b", l_G[C_B開始..]) };
            var l_ライブラリ = V_書き出し_ペア("threads", l_G, V_断片の開始群(0, 5700, C_間隔_密));

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 1 };
            var l_一 = PairBridgeFiller.Get_繋いだ配列群(l_配列群, l_ライブラリ, out var l_統計一);
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 4 };
            var l_四 = PairBridgeFiller.Get_繋いだ配列群(l_配列群, l_ライブラリ, out var l_統計四);

            Assert.Equal(l_一, l_四);
            Assert.Equal(l_統計一, l_統計四);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 区間 [p_始, p_終) が間 [C_A長, C_B開始) と重なるか
        /// </summary>
        /// <param name="p_始">区間の始め</param>
        /// <param name="p_終">区間の終わり (含まない)</param>
        /// <returns>重なるなら true</returns>
        private static bool V_重なる(int p_始, int p_終)
        {
            return p_始 < C_B開始 && p_終 > C_A長;
        }

        /// <summary>
        /// 断片の開始の並びを作る
        /// </summary>
        /// <param name="p_始">最初の開始</param>
        /// <param name="p_終">最後の開始 (含む)</param>
        /// <param name="p_間隔">間隔</param>
        /// <returns>開始の並び</returns>
        private static List<int> V_断片の開始群(int p_始, int p_終, int p_間隔)
        {
            var l_一覧 = new List<int>();
            for (var s = p_始; s <= p_終; s += p_間隔)
            {
                l_一覧.Add(s);
            }

            return l_一覧;
        }

        /// <summary>
        /// 乱数の種を固定したゲノム (ACGT のみ) を作る
        /// </summary>
        /// <param name="p_長さ">長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>ゲノム</returns>
        private static string V_生成_ゲノム(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// A (= G[0..2800]) を作る
        /// </summary>
        /// <returns>配列 A</returns>
        private static string V_A()
        {
            return V_生成_ゲノム(C_ゲノム長, C_ゲノムの種)[..C_A長];
        }

        /// <summary>
        /// A だけを対象にした頂点の文字列と種の索引を作る
        /// </summary>
        /// <returns>頂点の文字列と種の索引</returns>
        private static (string[] 頂点, Dictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> 表) V_A単独の索引()
        {
            var l_配列群 = new List<(string A_ID, string A_配列)> { ("a", V_A()) };
            var l_頂点 = PairBridgeFiller.Get_頂点群(l_配列群, out var l_対象);
            var l_表 = PairBridgeFiller.Get_種の索引(l_頂点, l_対象);
            return (l_頂点, l_表);
        }

        /// <summary>
        /// 当たりの頂点と位置、および d (頂点の長さ − 位置) を確かめる
        /// </summary>
        /// <param name="p_当たり">当たり</param>
        /// <param name="p_頂点">期待する頂点</param>
        /// <param name="p_位置">期待する位置</param>
        /// <param name="p_頂点群">頂点の文字列</param>
        /// <param name="p_d">期待する d</param>
        private static void V_当たりを確かめる((int A_頂点, int A_位置)? p_当たり, int p_頂点, int p_位置, string[] p_頂点群, int p_d)
        {
            Assert.True(p_当たり.HasValue);
            Assert.Equal(p_頂点, p_当たり.Value.A_頂点);
            Assert.Equal(p_位置, p_当たり.Value.A_位置);
            Assert.Equal(p_d, p_頂点群[p_当たり.Value.A_頂点].Length - p_当たり.Value.A_位置);
        }

        /// <summary>
        /// 断片ごとに FASTQ のペア (R1 = 先頭 100 bp、R2 = 末尾 100 bp の逆相補) を書き出す
        /// </summary>
        /// <param name="p_接頭辞">ファイル名の接頭辞</param>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_開始群">断片の開始の並び</param>
        /// <returns>ライブラリ (順リードと逆リードのパス)</returns>
        private List<(string A_順リード, string A_逆リード)> V_書き出し_ペア(string p_接頭辞, string p_ゲノム, IReadOnlyList<int> p_開始群)
        {
            var l_パス1 = Path.Combine(this._作業ディレクトリ, p_接頭辞 + "_1.fq");
            var l_パス2 = Path.Combine(this._作業ディレクトリ, p_接頭辞 + "_2.fq");
            using var l_書き込み1 = new StreamWriter(l_パス1);
            using var l_書き込み2 = new StreamWriter(l_パス2);
            for (var i = 0; i < p_開始群.Count; i++)
            {
                var l_s = p_開始群[i];
                var l_R1 = p_ゲノム.Substring(l_s, C_リード長);
                var l_R2 = Util.V_逆相補_曖昧塩基あり(p_ゲノム.Substring(l_s + C_挿入長 - C_リード長, C_リード長));
                l_書き込み1.WriteLine($"@pair{i}/1\n{l_R1}\n+\n{new string('I', C_リード長)}");
                l_書き込み2.WriteLine($"@pair{i}/2\n{l_R2}\n+\n{new string('I', C_リード長)}");
            }

            return [(l_パス1, l_パス2)];
        }

        #endregion
    }
}
