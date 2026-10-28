using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ロングリードを錨に当て、互いに一意な端の組 (橋) を決める処理 (LongReadBridger) の検証
    /// </summary>
    public class LongReadBridgerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki longread tests
        /// </summary>
        private const string C_項目_tsumiki_longread_tests = "tsumiki_longread_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 錨の長さ
        /// </summary>
        private const int C_錨長 = 3000;

        /// <summary>
        /// 配列に無い間の長さ
        /// </summary>
        private const int C_間の長さ = 2000;

        /// <summary>
        /// 反復の長さ
        /// </summary>
        private const int C_反復長 = 1200;

        /// <summary>
        /// 錨の端から、跨ぐリードが錨に入る長さ
        /// </summary>
        private const int C_錨に入る長さ = 1500;

        /// <summary>
        /// 端の反復の場面で、跨ぐリードが A に入る長さ (反復 1200 bp の手前の一意な部分を 1300 bp 通る)
        /// </summary>
        private const int C_反復の手前まで入る長さ = 2500;

        /// <summary>
        /// リードに入れる置換の割合
        /// </summary>
        private const double C_置換の割合 = 0.03D;

        /// <summary>
        /// 一様に覆うリードの長さ
        /// </summary>
        private const int C_覆うリード長 = 2000;

        /// <summary>
        /// 一様に覆うリードの始まりの間隔
        /// </summary>
        private const int C_覆うリードの間隔 = 250;

        /// <summary>
        /// 錨の位置の許す誤差
        /// </summary>
        private const int C_位置の誤差 = 30;

        /// <summary>
        /// 間の長さの許す誤差
        /// </summary>
        private const int C_間の誤差 = 100;

        /// <summary>
        /// 当て方の場面のリードの始まり (錨上)
        /// </summary>
        private const int C_当て方の始まり = 1000;

        /// <summary>
        /// 当て方の場面のリードの終わり (錨上)
        /// </summary>
        private const int C_当て方の終わり = 2500;

        /// <summary>
        /// 経路を探す k-mer の長さ
        /// </summary>
        private const int C_経路のk = 31;

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
        public LongReadBridgerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_longread_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
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
        /// 当て方: 錨の 1000〜2500 を 3% の置換で写したリードは順向きの当たり 1 つ、その逆相補は逆向きの当たり 1 つ
        /// </summary>
        [Fact]
        public void 当て方_写したリードは錨の同じ区間に当たり逆相補は逆向きになる()
        {
            var l_錨 = V_生成_配列(C_錨長, 1);
            var l_リード = V_置換する(l_錨[C_当て方の始まり..C_当て方の終わり], 101);
            var l_索引 = LongReadBridger.Get_索引([("a", l_錨)], [true]);

            var l_順 = LongReadBridger.Get_当たり群(l_リード, l_索引, [C_錨長]);
            var l_逆 = LongReadBridger.Get_当たり群(Util.V_逆相補_曖昧塩基あり(l_リード), l_索引, [C_錨長]);

            var l_順の当たり = Assert.Single(l_順);
            Assert.True(l_順の当たり.A_Is順向き);
            Assert.InRange(l_順の当たり.A_錨始, C_当て方の始まり - C_位置の誤差, C_当て方の始まり + C_位置の誤差);
            Assert.InRange(l_順の当たり.A_錨終, C_当て方の終わり - C_位置の誤差, C_当て方の終わり + C_位置の誤差);
            var l_逆の当たり = Assert.Single(l_逆);
            Assert.False(l_逆の当たり.A_Is順向き);
            Assert.InRange(l_逆の当たり.A_錨始, C_当て方の始まり - C_位置の誤差, C_当て方の始まり + C_位置の誤差);
            Assert.InRange(l_逆の当たり.A_錨終, C_当て方の終わり - C_位置の誤差, C_当て方の終わり + C_位置の誤差);
        }

        /// <summary>
        /// 二次の当たり: 錨 A と、A に 10% の置換を入れた似た錨 A' があるとき、A から写したリードの当たりは A だけを残す
        /// </summary>
        [Fact]
        public void 二次の当たり_似た錨への重なる当たりは捨てる()
        {
            var l_A = V_生成_配列(C_錨長, 41);
            var l_似たA = V_置換する(V_置換する(V_置換する(l_A, 42), 43), 44);
            var l_リード = V_置換する(l_A[C_当て方の始まり..C_当て方の終わり], 45);
            var l_索引 = LongReadBridger.Get_索引([("a", l_A), ("a2", l_似たA)], [true, true]);

            var l_当たり群 = LongReadBridger.Get_当たり群(l_リード, l_索引, [C_錨長, C_錨長]);

            var l_当たり = Assert.Single(l_当たり群);
            Assert.Equal(0, l_当たり.A_錨);
        }

        /// <summary>
        /// 曖昧な当たり: 同じ配列の錨が 2 本あるとき、そこから写したリードの当たりはどちらとも決められないので捨てる
        /// </summary>
        [Fact]
        public void 曖昧な当たり_同じ配列の錨への当たりは両方捨てる()
        {
            var l_A = V_生成_配列(C_錨長, 51);
            var l_リード = V_置換する(l_A[C_当て方の始まり..C_当て方の終わり], 52);
            var l_索引 = LongReadBridger.Get_索引([("a", l_A), ("a2", l_A)], [true, true]);

            var l_当たり群 = LongReadBridger.Get_当たり群(l_リード, l_索引, [C_錨長, C_錨長]);

            Assert.Empty(l_当たり群);
        }

        /// <summary>
        /// 基本の橋: A + R (配列に無い 2000 bp) + B を跨ぐ 3 本で、橋 A+ → B+ を決め、間は 2000±100
        /// </summary>
        [Fact]
        public void 基本_跨ぐ3本で橋を決め間の長さを見込む()
        {
            var (l_A, l_B, l_リード群) = V_A間Bの場面(3);

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", l_B)], l_リード群, out var l_統計);

            var l_橋 = Assert.Single(l_橋群);
            Assert.Equal((0, 2), (l_橋.A_u, l_橋.A_x));
            Assert.Equal(3, l_橋.A_票);
            Assert.InRange(l_橋.A_間の長さ, C_間の長さ - C_間の誤差, C_間の長さ + C_間の誤差);
            Assert.Equal(3, l_統計.A_リード数);
            Assert.Equal(1, l_統計.A_橋数);
        }

        /// <summary>
        /// 票が足りない: 跨ぐリードが 2 本なら橋を決めない
        /// </summary>
        [Fact]
        public void 票が足りない_跨ぐ2本では橋を決めない()
        {
            var (l_A, l_B, l_リード群) = V_A間Bの場面(2);

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", l_B)], l_リード群, out var l_統計);

            Assert.Empty(l_橋群);
            Assert.Equal(1, l_統計.A_候補数);
        }

        /// <summary>
        /// 向き: B を逆相補で渡すと、橋 A+ → B− を決める
        /// </summary>
        [Fact]
        public void 向き_Bを逆相補で渡すと逆向きの橋を決める()
        {
            var (l_A, l_B, l_リード群) = V_A間Bの場面(3);

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", Util.V_逆相補_曖昧塩基あり(l_B))], l_リード群, out _);

            var l_橋 = Assert.Single(l_橋群);
            Assert.Equal((0, 3), (l_橋.A_u, l_橋.A_x));
        }

        /// <summary>
        /// 多コピーの錨: A+IS+B と C+IS+D で、IS は多コピーとして飛ばし、橋 A→B と C→D だけを決める
        /// </summary>
        [Fact]
        public void 多コピーの錨_反復を飛ばして前後の錨を繋ぐ()
        {
            var l_A = V_生成_配列(C_錨長, 11);
            var l_B = V_生成_配列(C_錨長, 12);
            var l_C = V_生成_配列(C_錨長, 13);
            var l_D = V_生成_配列(C_錨長, 14);
            var l_IS = V_生成_配列(C_反復長, 15);
            var l_G1 = l_A + l_IS + l_B;
            var l_G2 = l_C + l_IS + l_D;
            var l_リード群 = new List<string>();
            l_リード群.AddRange(V_跨ぐリード群(l_G1, C_錨長 - C_錨に入る長さ, C_錨長 + C_反復長 + C_錨に入る長さ, 3, 200));
            l_リード群.AddRange(V_跨ぐリード群(l_G2, C_錨長 - C_錨に入る長さ, C_錨長 + C_反復長 + C_錨に入る長さ, 3, 300));

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", l_B), ("c", l_C), ("d", l_D), ("is", l_IS)], l_リード群, out var l_統計);

            Assert.Equal([(0, 2), (4, 6)], l_橋群.Select(x => (x.A_u, x.A_x)).Order());
            Assert.Equal(1, l_統計.A_多コピーの錨数);
            Assert.All(l_橋群, x => Assert.InRange(x.A_間の長さ, C_反復長 - C_間の誤差, C_反復長 + C_間の誤差));
        }

        /// <summary>
        /// 端に畳み込まれた反復: A の末尾 1200 bp が C と D の間にもある。A の反復の部分だけを通るリードは数えず、橋 A→B を決める
        /// </summary>
        [Fact]
        public void 端の反復_反復の部分だけの当たりは数えずAとBを繋ぐ()
        {
            var l_A = V_生成_配列(C_錨長, 21);
            var l_R1 = V_生成_配列(1000, 22);
            var l_B = V_生成_配列(C_錨長, 23);
            var l_C = V_生成_配列(C_錨長, 24);
            var l_D = V_生成_配列(C_錨長, 25);
            var l_G1 = l_A + l_R1 + l_B;
            var l_G2 = l_C + l_A[(C_錨長 - C_反復長)..] + l_D;
            var l_リード群 = new List<string>();
            l_リード群.AddRange(V_覆うリード群(l_G1, 400));
            l_リード群.AddRange(V_覆うリード群(l_G2, 500));
            l_リード群.AddRange(V_跨ぐリード群(l_G1, C_錨長 - C_反復の手前まで入る長さ, C_錨長 + 1000 + C_錨に入る長さ, 3, 600));
            l_リード群.AddRange(V_跨ぐリード群(l_G2, C_錨長 - C_錨に入る長さ, C_錨長 + C_反復長 + C_錨に入る長さ, 3, 700));

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", l_B), ("c", l_C), ("d", l_D)], l_リード群, out var l_統計);

            Assert.Contains(l_橋群, x => (x.A_u, x.A_x) == (0, 2));
            Assert.DoesNotContain(l_橋群, x => (x.A_u, x.A_x) != (0, 2) && (x.A_u / 2 == 0 || x.A_x / 2 == 0));
            Assert.Equal(0, l_統計.A_多コピーの錨数);
        }

        /// <summary>
        /// 端まで届かない当たり: A の先頭 1500 と B の先頭 1500 をつないだキメラのリード 3 本では、A の端まで 1500 残るので橋を決めない
        /// </summary>
        [Fact]
        public void 端まで届かない当たり_キメラのリードでは橋を決めない()
        {
            var l_A = V_生成_配列(C_錨長, 31);
            var l_B = V_生成_配列(C_錨長, 32);
            var l_キメラ = l_A[..C_錨に入る長さ] + l_B[..C_錨に入る長さ];

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", l_B)], V_跨ぐリード群(l_キメラ, 0, l_キメラ.Length, 3, 800), out var l_統計);

            Assert.Empty(l_橋群);
            Assert.Equal(0, l_統計.A_候補数);
        }

        /// <summary>
        /// 途中から入る当たり: C + (A の末尾 1200 bp) + D だけを跨ぐリードでは、A の途中から入るので A の端の橋を決めない (深さでは反復に見えない場面)
        /// </summary>
        [Fact]
        public void 途中から入る当たり_錨の途中から入ったリードは端の証拠にしない()
        {
            var l_A = V_生成_配列(C_錨長, 61);
            var l_B = V_生成_配列(C_錨長, 62);
            var l_C = V_生成_配列(C_錨長, 63);
            var l_D = V_生成_配列(C_錨長, 64);
            var l_G2 = l_C + l_A[(C_錨長 - C_反復長)..] + l_D;
            var l_リード群 = V_跨ぐリード群(l_G2, C_錨長 - C_錨に入る長さ, C_錨長 + C_反復長 + C_錨に入る長さ, 3, 900);

            var l_橋群 = LongReadBridger.Get_決めた橋([("a", l_A), ("b", l_B), ("c", l_C), ("d", l_D)], l_リード群, out _);

            Assert.DoesNotContain(l_橋群, x => x.A_u / 2 == 0 || x.A_x / 2 == 0);
        }

        /// <summary>
        /// 繋ぐ (重なり): A = G[0..3000]、B = G[2908..6000] を、橋 A+→B+ (間 −92) で畳み、1 本の G にする
        /// </summary>
        [Fact]
        public void 繋ぐ_重なりで畳んで1本にする()
        {
            var l_G = V_生成_配列(6000, 71);
            var l_A = l_G[..C_錨長];
            var l_B = l_G[(C_錨長 - 92)..];

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, -92)], null, null, null, C_経路のk, out var l_繋ぎ方群, out var l_統計);

            Assert.Equal([l_G], l_結果.Select(x => x.A_配列));
            Assert.Equal(92, Assert.Single(l_繋ぎ方群).A_重なり);
            Assert.Equal(1, l_統計.A_重なりで畳んだ数);
        }

        /// <summary>
        /// 繋がない (縦列反復の重なり): 23 bp の単位が 4 回並ぶ所で A と B が重なると、一致する重なりは周期を持つ見かけのものなので畳まない
        /// </summary>
        [Fact]
        public void 繋がない_縦列反復の中の周期を持つ重なりは使わない()
        {
            var (l_A, l_B, _) = V_縦列反復の場面(69, 23);

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, -46)], null, null, null, C_経路のk, out var l_繋ぎ方群, out _);

            Assert.Equal(2, l_結果.Count);
            Assert.Null(Assert.Single(l_繋ぎ方群).A_重なり);
        }

        /// <summary>
        /// 繋ぐ (縦列反復を経路で): A が単位 1 つと 18 bp で終わり、B が単位 3 つ目から始まる (本当の間は 5 bp)。見かけの重なり 41 bp は使わず、
        /// 縦列反復を含む材料の経路で単位 4 つの G に戻す
        /// </summary>
        [Fact]
        public void 繋ぐ_縦列反復は見かけの重なりでなく経路で埋める()
        {
            var (l_A, l_B, l_G) = V_縦列反復の場面(41, 46);
            var l_材料 = new List<string> { l_G[(C_錨長 - 300)..(C_錨長 + 92 + 300)] };

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, -41)], l_材料, null, null, C_経路のk, out var l_繋ぎ方群, out _);

            Assert.Equal([l_G], l_結果.Select(x => x.A_配列));
            Assert.Equal(5, Assert.Single(l_繋ぎ方群).A_埋める配列!.Length);
        }

        /// <summary>
        /// 繋ぐ (経路): G1 = A + IS + B、G2 = C + IS + D の組み立てグラフの材料から、橋 A+→B+ (間 1200) を IS の経路で埋めて 1 本の G1 にする
        /// </summary>
        [Fact]
        public void 繋ぐ_グラフの経路で埋めて1本にする()
        {
            var (l_A, l_B, l_IS, l_材料) = V_ISの経路の場面();

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, C_反復長)], l_材料, null, null, C_経路のk, out var l_繋ぎ方群, out var l_統計);

            Assert.Equal([l_A + l_IS + l_B], l_結果.Select(x => x.A_配列));
            Assert.Equal(C_反復長, Assert.Single(l_繋ぎ方群).A_埋める配列!.Length);
            Assert.Equal(1, l_統計.A_経路で埋めた数);
        }

        /// <summary>
        /// 繋ぐ (経路、逆向き): B を逆相補で渡し、橋 A+→B− で埋めると、G1 かその逆相補になる
        /// </summary>
        [Fact]
        public void 繋ぐ_逆相補のBでも経路で埋める()
        {
            var (l_A, l_B, l_IS, l_材料) = V_ISの経路の場面();

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", Util.V_逆相補_曖昧塩基あり(l_B))], [new ロングリードの橋(0, 3, 5, 0, C_反復長)], l_材料, null, null, C_経路のk, out _, out _);

            var l_配列 = Assert.Single(l_結果).A_配列;
            var l_G1 = l_A + l_IS + l_B;
            Assert.True(l_配列 == l_G1 || l_配列 == Util.V_逆相補_曖昧塩基あり(l_G1));
        }

        /// <summary>
        /// 繋がない: 材料に間の配列が無ければ unreachable、見込みの長さが外れれば (600、本当は 1200) 範囲外で繋がない
        /// </summary>
        [Fact]
        public void 繋がない_材料が無いか長さが外れれば元のまま()
        {
            var (l_A, l_B, _, l_材料) = V_ISの経路の場面();

            var l_材料なし = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, C_反復長)], [], null, null, C_経路のk, out var l_繋ぎ方群1, out var l_統計1);
            var l_長さ違い = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, 600)], l_材料, null, null, C_経路のk, out var l_繋ぎ方群2, out _);

            Assert.Equal([l_A, l_B], l_材料なし.Select(x => x.A_配列));
            Assert.Equal("unreachable", Assert.Single(l_繋ぎ方群1).A_理由);
            Assert.Equal(1, l_統計1.A_決まらなかった数);
            Assert.Equal([l_A, l_B], l_長さ違い.Select(x => x.A_配列));
            Assert.NotNull(Assert.Single(l_繋ぎ方群2).A_理由);
        }

        /// <summary>
        /// 道選び (反復の 2 つのコピー): R1 と、100 bp ごとに置換を入れた R2 が材料にあり、経路は一意でないが、A→B を跨ぐ 5 本の多数決で G1 = A + R1 + B に戻す
        /// </summary>
        [Fact]
        public void 道選び_反復のコピーをロングリードの多数決で選ぶ()
        {
            var (l_A, l_B, l_G1, l_材料) = V_2コピーの反復の場面();
            var l_区間群 = V_跨ぐ区間群(l_G1, l_A.Length, l_G1.Length - l_B.Length, 5, 1000);

            var l_L2 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, C_反復長)], l_材料, null, null, C_経路のk, out var l_L2の繋ぎ方, out _);
            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, C_反復長)], l_材料, [l_区間群], null, C_経路のk, out var l_繋ぎ方群, out var l_統計);

            Assert.Equal(2, l_L2.Count);
            Assert.Equal("ambiguous", Assert.Single(l_L2の繋ぎ方).A_理由);
            Assert.Equal([l_G1], l_結果.Select(x => x.A_配列));
            Assert.True(Assert.Single(l_繋ぎ方群).A_Isロングリードで選んだ);
            Assert.Equal(1, l_統計.A_ロングリードで選んだ数);
        }

        /// <summary>
        /// 道選び (リードが足りない): 跨ぐ区間が 2 本なら繋がない
        /// </summary>
        [Fact]
        public void 道選び_区間が2本なら繋がない()
        {
            var (l_A, l_B, l_G1, l_材料) = V_2コピーの反復の場面();
            var l_区間群 = V_跨ぐ区間群(l_G1, l_A.Length, l_G1.Length - l_B.Length, 2, 1000);

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, C_反復長)], l_材料, [l_区間群], null, C_経路のk, out var l_繋ぎ方群, out _);

            Assert.Equal(2, l_結果.Count);
            Assert.NotNull(Assert.Single(l_繋ぎ方群).A_理由);
        }

        /// <summary>
        /// 道選び (縦列反復): 単位 100 bp が 4 つ並ぶ所で、材料には単位の輪がある。経路は一意でないが、跨ぐ 5 本の多数決で単位 4 つの G に戻す (見込み 430 では単位 4 つと 5 つの経路が L2 の幅に入る)
        /// </summary>
        [Fact]
        public void 道選び_縦列反復の単位の数をロングリードで決める()
        {
            var l_単位 = V_生成_配列(100, 91);
            var l_A = V_生成_配列(C_錨長, 92);
            var l_B = V_生成_配列(C_錨長, 93);
            var l_G = l_A + string.Concat(Enumerable.Repeat(l_単位, 4)) + l_B;
            List<string> l_材料 = [l_A[^300..] + l_単位 + l_単位, l_単位 + l_単位 + l_B[..300]];
            var l_区間群 = V_跨ぐ区間群(l_G, l_A.Length, l_G.Length - l_B.Length, 5, 1100);

            var l_L2 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, 430)], l_材料, null, null, C_経路のk, out var l_L2の繋ぎ方, out _);
            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 5, 0, 430)], l_材料, [l_区間群], null, C_経路のk, out var l_繋ぎ方群, out _);

            Assert.Equal("ambiguous", Assert.Single(l_L2の繋ぎ方).A_理由);
            Assert.Equal([l_G], l_結果.Select(x => x.A_配列));
            Assert.True(Assert.Single(l_繋ぎ方群).A_Isロングリードで選んだ);
        }

        /// <summary>
        /// 直す: 材料には反復の多数派の配列だけがあり、ゲノムの反復には 1 か所の置換がある。跨ぐ 8 本と短いリードの確かめで、埋めた配列をゲノムどおりに直す
        /// </summary>
        [Fact]
        public void 直す_コピー固有の置換をロングリードと短いリードで直す()
        {
            var (l_A, l_B, l_G, l_材料) = V_コピー固有の置換の場面(true);
            var l_区間群 = V_跨ぐ区間群(l_G, l_A.Length, l_G.Length - l_B.Length, 8, 1200);

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 8, 0, C_反復長)], l_材料, [l_区間群], V_短いリードの確かめ(l_G), C_経路のk, out _, out var l_統計);

            Assert.Equal([l_G], l_結果.Select(x => x.A_配列));
            Assert.Equal(1, l_統計.A_直した塩基数);
        }

        /// <summary>
        /// 直さない: 短いリードの確かめが通らなければ、多数派の配列のまま埋める
        /// </summary>
        [Fact]
        public void 直さない_短いリードで確かめられなければ直さない()
        {
            var (l_A, l_B, l_G, l_材料) = V_コピー固有の置換の場面(true);
            var l_区間群 = V_跨ぐ区間群(l_G, l_A.Length, l_G.Length - l_B.Length, 8, 1200);

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 8, 0, C_反復長)], l_材料, [l_区間群], _ => false, C_経路のk, out _, out var l_統計);

            Assert.NotEqual(l_G, Assert.Single(l_結果).A_配列);
            Assert.Equal(l_G.Length, l_結果[0].A_配列.Length);
            Assert.Equal(0, l_統計.A_直した塩基数);
        }

        /// <summary>
        /// 直さない: 区間が埋めた配列と一致していれば何も直さない
        /// </summary>
        [Fact]
        public void 直さない_区間が経路と一致していれば直さない()
        {
            var (l_A, l_B, l_G, l_材料) = V_コピー固有の置換の場面(false);
            var l_区間群 = V_跨ぐ区間群(l_G, l_A.Length, l_G.Length - l_B.Length, 8, 1300);

            var l_結果 = LongReadBridger.Get_繋いだ配列群([("a", l_A), ("b", l_B)], [new ロングリードの橋(0, 2, 8, 0, C_反復長)], l_材料, [l_区間群], V_短いリードの確かめ(l_G), C_経路のk, out _, out var l_統計);

            Assert.Equal([l_G], l_結果.Select(x => x.A_配列));
            Assert.Equal(0, l_統計.A_直した塩基数);
        }

        /// <summary>
        /// 引数: -long はカンマ区切りで複数のパスを持ち、無いパスは例外
        /// </summary>
        [Fact]
        public void 引数_複数のパスを持ち無いパスは例外()
        {
            var l_リード = Path.Combine(this._作業ディレクトリ, "r.fq");
            var l_ロング1 = Path.Combine(this._作業ディレクトリ, "l1.fq");
            var l_ロング2 = Path.Combine(this._作業ディレクトリ, "l2.fq");
            foreach (var l_パス in new[] { l_リード, l_ロング1, l_ロング2 })
            {
                File.WriteAllText(l_パス, "@r\nACGT\n+\nIIII\n");
            }

            var l_引数 = ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_リード, Consts.引数キー.ロングリードのパス, $"{l_ロング1},{l_ロング2}"]);

            Assert.Equal([l_ロング1, l_ロング2], l_引数.A_ロングリードのパス群);
            _ = Assert.ThrowsAny<Exception>(() => ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_リード, Consts.引数キー.ロングリードのパス, Path.Combine(this._作業ディレクトリ, "none.fq")]));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 一意な配列 3000 bp + 23 bp の単位 4 つ + 一意な配列 3000 bp を作り、A を単位の中で終わらせ、B を単位の中から始める
        /// </summary>
        /// <param name="p_Aに入る単位の長さ">A に含める反復の長さ</param>
        /// <param name="p_Bの始まり">B が始まる反復の中の位置</param>
        /// <returns>A、B、ゲノム</returns>
        private static (string A_A, string A_B, string A_G) V_縦列反復の場面(int p_Aに入る単位の長さ, int p_Bの始まり)
        {
            var l_単位 = V_生成_配列(23, 72);
            var l_G = V_生成_配列(C_錨長, 73) + string.Concat(Enumerable.Repeat(l_単位, 4)) + V_生成_配列(C_錨長, 74);
            return (l_G[..(C_錨長 + p_Aに入る単位の長さ)], l_G[(C_錨長 + p_Bの始まり)..], l_G);
        }

        /// <summary>
        /// 反復 R1 (1200 bp) と、R1 に 100 bp ごとに置換を入れた R2 を作り、G1 = A + R1 + B、G2 = C + R2 + D の材料を作る
        /// </summary>
        /// <returns>A、B、G1、材料</returns>
        private static (string A_A, string A_B, string A_G1, List<string> A_材料) V_2コピーの反復の場面()
        {
            var l_A = V_生成_配列(C_錨長, 101);
            var l_B = V_生成_配列(C_錨長, 102);
            var l_C = V_生成_配列(C_錨長, 103);
            var l_D = V_生成_配列(C_錨長, 104);
            var l_R1 = V_生成_配列(C_反復長, 105);
            var l_文字 = l_R1.ToCharArray();
            for (var i = 50; i < l_文字.Length; i += 100)
            {
                l_文字[i] = Consts.塩基文字[(Consts.塩基文字.IndexOf(l_文字[i]) + 1) % 4];
            }

            var l_R2 = new string(l_文字);
            List<string> l_材料 = [l_A[^300..] + l_R1 + l_B[..300], l_C[^300..] + l_R2 + l_D[..300]];
            return (l_A, l_B, l_A + l_R1 + l_B, l_材料);
        }

        /// <summary>
        /// 反復 R (1200 bp) の多数派の配列だけを材料にし、ゲノム A + R' + B を作る (置換ありなら R' は R の 600 bp 目を置換したもの)
        /// </summary>
        /// <param name="p_Is置換あり">ゲノムの反復に置換を入れるか</param>
        /// <returns>A、B、ゲノム、材料</returns>
        private static (string A_A, string A_B, string A_G, List<string> A_材料) V_コピー固有の置換の場面(bool p_Is置換あり)
        {
            var l_A = V_生成_配列(C_錨長, 111);
            var l_B = V_生成_配列(C_錨長, 112);
            var l_R = V_生成_配列(C_反復長, 113);
            var l_文字 = l_R.ToCharArray();
            if (p_Is置換あり)
            {
                l_文字[600] = Consts.塩基文字[(Consts.塩基文字.IndexOf(l_文字[600]) + 1) % 4];
            }

            List<string> l_材料 = [l_A[^300..] + l_R + l_B[..300]];
            return (l_A, l_B, l_A + new string(l_文字) + l_B, l_材料);
        }

        /// <summary>
        /// ゲノムの 41-mer (両向き) に含まれるかで、短いリードの確かめを真似る
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <returns>確かめる関数</returns>
        private static Func<string, bool> V_短いリードの確かめ(string p_ゲノム)
        {
            var l_長さ = (2 * LongReadBridger.C_確かめるkmerの片側の長さ) + 1;
            var l_集合 = new HashSet<string>();
            for (var i = 0; i + l_長さ <= p_ゲノム.Length; i++)
            {
                _ = l_集合.Add(p_ゲノム.Substring(i, l_長さ));
            }

            return x => l_集合.Contains(x) || l_集合.Contains(Util.V_逆相補_曖昧塩基あり(x));
        }

        /// <summary>
        /// ゲノムの [A の終わり − 500, B の始まり + 500] を、種を変えた置換で写した区間 (u→x の向き) を作る
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_Aの終わり">A の終わりの位置</param>
        /// <param name="p_Bの始まり">B の始まりの位置</param>
        /// <param name="p_本数">本数</param>
        /// <param name="p_種">乱数の種の始まり</param>
        /// <returns>区間</returns>
        private static List<string> V_跨ぐ区間群(string p_ゲノム, int p_Aの終わり, int p_Bの始まり, int p_本数, int p_種)
        {
            return V_跨ぐリード群(p_ゲノム, p_Aの終わり - 500, p_Bの始まり + 500, p_本数, p_種);
        }

        /// <summary>
        /// G1 = A + IS + B、G2 = C + IS + D の組み立てグラフの材料 (unitig が k−1 重なる形) を作る
        /// </summary>
        /// <returns>A、B、IS、材料</returns>
        private static (string A_A, string A_B, string A_IS, List<string> A_材料) V_ISの経路の場面()
        {
            var l_A = V_生成_配列(C_錨長, 81);
            var l_B = V_生成_配列(C_錨長, 82);
            var l_C = V_生成_配列(C_錨長, 83);
            var l_D = V_生成_配列(C_錨長, 84);
            var l_IS = V_生成_配列(C_反復長, 85);
            var l_重なり = C_経路のk - 1;
            List<string> l_材料 =
            [
                l_A[^300..] + l_IS[..l_重なり],
                l_IS,
                l_IS[^l_重なり..] + l_B[..300],
                l_C[^300..] + l_IS[..l_重なり],
                l_IS[^l_重なり..] + l_D[..300],
            ];
            return (l_A, l_B, l_IS, l_材料);
        }

        /// <summary>
        /// A + R (配列に無い) + B のゲノムと、A の末尾 1500 から B の先頭 1500 までを跨ぐリードを作る
        /// </summary>
        /// <param name="p_本数">跨ぐリードの本数</param>
        /// <returns>A、B、リード</returns>
        private static (string A_A, string A_B, List<string> A_リード群) V_A間Bの場面(int p_本数)
        {
            var l_A = V_生成_配列(C_錨長, 2);
            var l_R = V_生成_配列(C_間の長さ, 3);
            var l_B = V_生成_配列(C_錨長, 4);
            var l_G = l_A + l_R + l_B;
            return (l_A, l_B, V_跨ぐリード群(l_G, C_錨長 - C_錨に入る長さ, C_錨長 + C_間の長さ + C_錨に入る長さ, p_本数, 100));
        }

        /// <summary>
        /// ゲノムの同じ区間を、種を変えた置換で写したリードを作る
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_始め">区間の始まり</param>
        /// <param name="p_終わり">区間の終わり</param>
        /// <param name="p_本数">本数</param>
        /// <param name="p_種">乱数の種の始まり</param>
        /// <returns>リード</returns>
        private static List<string> V_跨ぐリード群(string p_ゲノム, int p_始め, int p_終わり, int p_本数, int p_種)
        {
            return [.. Enumerable.Range(0, p_本数).Select(i => V_置換する(p_ゲノム[p_始め..p_終わり], p_種 + i))];
        }

        /// <summary>
        /// ゲノム全体を一様に覆うリード (長さ C_覆うリード長、間隔 C_覆うリードの間隔) を作る
        /// </summary>
        /// <param name="p_ゲノム">ゲノム</param>
        /// <param name="p_種">乱数の種の始まり</param>
        /// <returns>リード</returns>
        private static List<string> V_覆うリード群(string p_ゲノム, int p_種)
        {
            var l_結果 = new List<string>();
            for (var s = 0; s + C_覆うリード長 <= p_ゲノム.Length; s += C_覆うリードの間隔)
            {
                l_結果.Add(V_置換する(p_ゲノム.Substring(s, C_覆うリード長), p_種 + s));
            }

            return l_結果;
        }

        /// <summary>
        /// 乱数の配列を作る
        /// </summary>
        /// <param name="p_長さ">長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>配列</returns>
        private static string V_生成_配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return new string([.. Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
        }

        /// <summary>
        /// 配列に C_置換の割合 の置換を入れる
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>置換を入れた配列</returns>
        private static string V_置換する(string p_配列, int p_種)
        {
            var l_乱数 = new Random(p_種);
            var l_文字 = p_配列.ToCharArray();
            for (var i = 0; i < l_文字.Length; i++)
            {
                if (l_乱数.NextDouble() < C_置換の割合)
                {
                    l_文字[i] = Consts.塩基文字[(Consts.塩基文字.IndexOf(l_文字[i]) + 1 + l_乱数.Next(3)) % 4];
                }
            }

            return new string(l_文字);
        }

        #endregion
    }
}
