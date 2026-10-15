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

        /// <summary>
        /// 間 30 bp の場面での配列 B の開始 (B = G[2830..6000])
        /// </summary>
        private const int C_B開始_間30bp = 2830;

        /// <summary>
        /// 間 30 bp の場面で、A と B を跨ぐ断片の長さ (間の見込みが負になる長さ)
        /// </summary>
        private const int C_跨ぐ断片長_間30bp = 370;

        /// <summary>
        /// 反復の端の場面で、B の先頭に複製する A の末尾の長さ
        /// </summary>
        private const int C_反復の長さ = 300;

        /// <summary>
        /// 反復の端の場面で、A と B を跨ぐ断片の長さ (B の先頭の複製のぶん、読みが B の一意な部分まで届く長さ)
        /// </summary>
        private const int C_跨ぐ断片長_反復 = 450;

        /// <summary>
        /// 重なりで畳む場面で、A の末尾と B の先頭が一致する長さ (本当の重なり)
        /// </summary>
        private const int C_重なりの長さ = 100;

        /// <summary>
        /// 多コピーの場面で足す、単一コピーの配列の本数 (基準の密度を作るのに要る頂点数を超える数)
        /// </summary>
        private const int C_多コピーの配列数 = 10;

        /// <summary>
        /// 多コピーの場面で足す、単一コピーの配列の長さ
        /// </summary>
        private const int C_多コピーの配列の長さ = 1000;

        /// <summary>
        /// 多コピーの場面で足す配列の乱数の種の始まり
        /// </summary>
        private const int C_多コピーの種 = 20000;

        /// <summary>
        /// 多コピーの場面で 2.5 倍にする B の先頭の範囲 (端の窓 C_窓長 と同じ 1000 bp)
        /// </summary>
        private const int C_B先頭の窓 = 1000;

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
        /// 印で繋ぐ: 間を跨ぐリードが無く埋まらなくても、A と B を結ぶペアが 10 組以上あれば、N 列で 1 本に繋ぐ。
        /// 両側の片の端 30 bp を削り、N 列は見込み 50 bp + 削った 60 bp = 110 個
        /// </summary>
        [Fact]
        public void 印で繋ぐ_埋まらなくても票が10組以上あればN列で1本に繋ぐ()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_断片群 = V_場面の断片群(C_B開始, V_断片の開始群(2650, 2698, 3), C_挿入長);
            var l_ライブラリ = V_書き出し_断片("mark", l_G, l_断片群);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Single(l_結果);
            Assert.Equal(l_A[..^30] + new string('N', 110) + l_B[30..], l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_候補の辺の数);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
            Assert.Equal(1, l_統計.A_印で繋いだ辺の数);
        }

        /// <summary>
        /// 印の下限未満: 埋まらず、A と B を結ぶペアが 10 組未満 (6 組) なら、印では繋がず元の 2 本のまま
        /// </summary>
        [Fact]
        public void 印で繋がない_埋まらず票が10組未満なら元の2本のまま()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_跨ぐ開始群 = V_断片の開始群(2650, 2698, 3).Where((s, i) => i % 3 == 0).ToList();
            var l_断片群 = V_場面の断片群(C_B開始, l_跨ぐ開始群, C_挿入長);
            var l_ライブラリ = V_書き出し_断片("fewmark", l_G, l_断片群);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Equal(2, l_結果.Count);
            Assert.Equal(l_A, l_結果[0].A_配列);
            Assert.Equal(l_B, l_結果[1].A_配列);
            Assert.Equal(1, l_統計.A_候補の辺の数);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
            Assert.Equal(0, l_統計.A_印で繋いだ辺の数);
        }

        /// <summary>
        /// 間の見込みが負: 間 30 bp で、跨ぐ断片の長さ 370 により間の見込みは −40 (旧下限 −30 より小さい)。
        /// 候補になり、票が 10 組以上なので、埋まれば埋めて、埋まらなければ印で 1 本に繋ぐ
        /// </summary>
        [Fact]
        public void 間の見込みが負でも候補になり票が10組以上なら1本に繋ぐ()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始_間30bp..];
            var l_跨ぐ開始群 = V_断片の開始群(2612, 2650, 2);
            var l_断片群 = V_場面の断片群(C_B開始_間30bp, l_跨ぐ開始群, C_跨ぐ断片長_間30bp);
            var l_ライブラリ = V_書き出し_断片("neg", l_G, l_断片群);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Equal(1, l_統計.A_候補の辺の数);
            Assert.Equal(C_挿入長, l_統計.A_挿入長);
            Assert.Equal(1, l_統計.A_埋めて繋いだ辺の数 + l_統計.A_印で繋いだ辺の数);
            Assert.Single(l_結果);
            if (l_統計.A_印で繋いだ辺の数 == 1)
            {
                Assert.Equal(l_A[..^30] + new string('N', 70) + l_B[30..], l_結果[0].A_配列);
            }
        }

        /// <summary>
        /// 反復の端は繋がない (共有の 25-mer): B の先頭に A の末尾 300 bp と同じ配列を置く。票があっても反復の端として除き、繋がない
        /// </summary>
        [Fact]
        public void 反復の端は繋がない_先頭に末尾の複製があれば除いて繋がない()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_A[(C_A長 - C_反復の長さ)..] + l_G[C_B開始..];
            var l_断片群 = V_場面の断片群(C_B開始, V_断片の開始群(2440, 2470, 5), C_跨ぐ断片長_反復);
            var l_ライブラリ = V_書き出し_断片("dupend", l_G, l_断片群);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Equal(2, l_結果.Count);
            Assert.Equal(l_A, l_結果[0].A_配列);
            Assert.Equal(l_B, l_結果[1].A_配列);
            Assert.Equal(1, l_統計.A_候補の辺の数);
            Assert.Equal(1, l_統計.A_反復の端で除いた数);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
            Assert.Equal(0, l_統計.A_重なりで畳んだ辺の数);
        }

        /// <summary>
        /// 重なりで畳む: A の末尾 100 bp と B の先頭 100 bp が同じ (本当の重なり)。間を読むリードが無く埋まらないので、重なりで畳んで A + B[100..] に繋ぐ
        /// </summary>
        [Fact]
        public void 重なりで畳む_末尾と先頭が100bp重なれば畳んで1本に繋ぐ()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[(C_A長 - C_重なりの長さ)..];
            var l_断片群 = V_場面の断片群(C_B開始, V_断片の開始群(2510, 2670, 10), C_挿入長);
            var l_ライブラリ = V_書き出し_断片("overlap", l_G, l_断片群);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群([("a", l_A), ("b", l_B)], l_ライブラリ, out var l_統計);

            Assert.Single(l_結果);
            Assert.Equal(l_A + l_B[C_重なりの長さ..], l_結果[0].A_配列);
            Assert.Equal(l_G, l_結果[0].A_配列);
            Assert.Equal(1, l_統計.A_重なりで畳んだ辺の数);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
            Assert.Equal(0, l_統計.A_印で繋いだ辺の数);
        }

        /// <summary>
        /// 多コピーの片はつながない: 埋まらないが票が多い場面で、単一コピーの配列を 10 本足し、B の先頭 1000 bp から出るリードを 2.5 倍にすると、
        /// B の先頭の端は多コピーの片と見なして繋がず、統計の多コピーの端で除いた数は 1
        /// </summary>
        [Fact]
        public void 多コピーの片はつながない_B先頭のリードを増やすと繋がず多コピーの端で除く()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_配列群 = V_多コピーの配列群(l_G[..C_A長], l_G[C_B開始..]);
            var l_ライブラリ = V_書き出し_多コピーの場面("multi_head", l_G, p_B先頭を増やす: true);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群(l_配列群, l_ライブラリ, out var l_統計);

            Assert.Equal(l_配列群.Count, l_結果.Count);
            Assert.Contains(l_配列群[0].A_配列, l_結果.Select(x => x.A_配列));
            Assert.Contains(l_配列群[1].A_配列, l_結果.Select(x => x.A_配列));
            Assert.DoesNotContain(l_結果, x => x.A_配列.Contains('N'));
            Assert.Equal(1, l_統計.A_多コピーの端で除いた数);
            Assert.Equal(0, l_統計.A_印で繋いだ辺の数);
            Assert.Equal(0, l_統計.A_埋めて繋いだ辺の数);
        }

        /// <summary>
        /// 多コピーの片の検査: 同じ場面で B の先頭のリードを増やさなければ、今どおり印の N 列で繋ぎ、多コピーの端で除く数は 0
        /// </summary>
        [Fact]
        public void 多コピーの片の検査_B先頭のリードを増やさなければ今どおり印で繋ぐ()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_A = l_G[..C_A長];
            var l_B = l_G[C_B開始..];
            var l_配列群 = V_多コピーの配列群(l_A, l_B);
            var l_ライブラリ = V_書き出し_多コピーの場面("multi_plain", l_G, p_B先頭を増やす: false);

            var l_結果 = PairBridgeFiller.Get_繋いだ配列群(l_配列群, l_ライブラリ, out var l_統計);

            Assert.Equal(l_配列群.Count - 1, l_結果.Count);
            Assert.Contains(l_A[..^30] + new string('N', 110) + l_B[30..], l_結果.Select(x => x.A_配列));
            Assert.Equal(1, l_統計.A_印で繋いだ辺の数);
            Assert.Equal(0, l_統計.A_多コピーの端で除いた数);
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
            return V_書き出し_断片(p_接頭辞, p_ゲノム, V_断片群(p_開始群, C_挿入長));
        }

        /// <summary>
        /// 断片ごとに FASTQ のペアを書き出す (断片の長さは断片ごと。R1 = 先頭 100 bp、R2 = 末尾 100 bp の逆相補)
        /// </summary>
        /// <param name="p_接頭辞">ファイル名の接頭辞</param>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_断片群">断片の開始と長さの並び</param>
        /// <returns>ライブラリ (順リードと逆リードのパス)</returns>
        private List<(string A_順リード, string A_逆リード)> V_書き出し_断片(string p_接頭辞, string p_ゲノム, IReadOnlyList<(int A_開始, int A_断片長)> p_断片群)
        {
            var l_パス1 = Path.Combine(this._作業ディレクトリ, p_接頭辞 + "_1.fq");
            var l_パス2 = Path.Combine(this._作業ディレクトリ, p_接頭辞 + "_2.fq");
            using var l_書き込み1 = new StreamWriter(l_パス1);
            using var l_書き込み2 = new StreamWriter(l_パス2);
            for (var i = 0; i < p_断片群.Count; i++)
            {
                var (l_s, l_長さ) = p_断片群[i];
                var l_R1 = p_ゲノム.Substring(l_s, C_リード長);
                var l_R2 = Util.V_逆相補_曖昧塩基あり(p_ゲノム.Substring(l_s + l_長さ - C_リード長, C_リード長));
                l_書き込み1.WriteLine($"@pair{i}/1\n{l_R1}\n+\n{new string('I', C_リード長)}");
                l_書き込み2.WriteLine($"@pair{i}/2\n{l_R2}\n+\n{new string('I', C_リード長)}");
            }

            return [(l_パス1, l_パス2)];
        }

        /// <summary>
        /// 開始の並びに同じ断片の長さを付ける
        /// </summary>
        /// <param name="p_開始群">断片の開始の並び</param>
        /// <param name="p_断片長">断片の長さ</param>
        /// <returns>断片の開始と長さの並び</returns>
        private static List<(int A_開始, int A_断片長)> V_断片群(IEnumerable<int> p_開始群, int p_断片長)
        {
            return p_開始群.Select(s => (s, p_断片長)).ToList();
        }

        /// <summary>
        /// 印・間の場面の断片の並びを作る: A と B の中だけの挿入長 300 の断片 (A と B を跨がないもの) と、A と B を跨ぐ断片
        /// </summary>
        /// <param name="p_B開始">配列 B の開始 (ゲノムの位置)</param>
        /// <param name="p_跨ぐ開始群">A と B を跨ぐ断片の開始の並び</param>
        /// <param name="p_跨ぐ断片長">A と B を跨ぐ断片の長さ</param>
        /// <returns>断片の開始と長さの並び</returns>
        private static List<(int A_開始, int A_断片長)> V_場面の断片群(int p_B開始, IEnumerable<int> p_跨ぐ開始群, int p_跨ぐ断片長)
        {
            var l_一覧 = V_断片群(V_断片の開始群(0, 5700, C_間隔_密).Where(s => s + C_挿入長 <= C_A長 || s >= p_B開始), C_挿入長);
            l_一覧.AddRange(V_断片群(p_跨ぐ開始群, p_跨ぐ断片長));
            return l_一覧;
        }

        /// <summary>
        /// 多コピーの場面の配列群を作る: A と B、それに単一コピーの配列を C_多コピーの配列数 本足す
        /// </summary>
        /// <param name="p_A">配列 A</param>
        /// <param name="p_B">配列 B</param>
        /// <returns>配列群 (ID と配列)</returns>
        private static List<(string A_ID, string A_配列)> V_多コピーの配列群(string p_A, string p_B)
        {
            var l_配列群 = new List<(string A_ID, string A_配列)> { ("a", p_A), ("b", p_B) };
            for (var k = 0; k < C_多コピーの配列数; k++)
            {
                l_配列群.Add(($"x{k}", V_生成_ゲノム(C_多コピーの配列の長さ, C_多コピーの種 + k)));
            }

            return l_配列群;
        }

        /// <summary>
        /// 多コピーの場面のリードを書き出す: 印・間の場面の断片に、必要なら B の先頭 1000 bp の断片を 2.5 倍にし、
        /// 足した単一コピーの配列それぞれの断片も書き出す
        /// </summary>
        /// <param name="p_接頭辞">ファイル名の接頭辞</param>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_B先頭を増やす">true なら B の先頭 1000 bp の断片を 2.5 倍にする</param>
        /// <returns>ライブラリ (順リードと逆リードのパス)</returns>
        private List<(string A_順リード, string A_逆リード)> V_書き出し_多コピーの場面(string p_接頭辞, string p_ゲノム, bool p_B先頭を増やす)
        {
            var l_断片群 = V_場面の断片群(C_B開始, V_断片の開始群(2650, 2698, 3), C_挿入長);
            if (p_B先頭を増やす)
            {
                var l_追加 = new List<(int A_開始, int A_断片長)>();
                var l_番号 = 0;
                foreach (var l_断片 in l_断片群)
                {
                    if (l_断片.A_開始 >= C_B開始 && l_断片.A_開始 < C_B開始 + C_B先頭の窓)
                    {
                        l_追加.Add(l_断片);
                        if (l_番号 % 2 == 0)
                        {
                            l_追加.Add(l_断片);
                        }

                        l_番号++;
                    }
                }

                l_断片群.AddRange(l_追加);
            }

            var l_ライブラリ = V_書き出し_断片(p_接頭辞, p_ゲノム, l_断片群);
            for (var k = 0; k < C_多コピーの配列数; k++)
            {
                var l_配列 = V_生成_ゲノム(C_多コピーの配列の長さ, C_多コピーの種 + k);
                l_ライブラリ.AddRange(V_書き出し_断片(p_接頭辞 + "_" + k, l_配列, V_断片群(V_断片の開始群(0, 700, C_間隔_密), C_挿入長)));
            }

            return l_ライブラリ;
        }

        #endregion
    }
}
