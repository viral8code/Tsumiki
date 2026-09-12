using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig 間の隣接を de Bruijn グラフから厳密に構築する UnitigGraph の検証
    /// </summary>
    /// <remarks>
    /// 旧実装 (リードマッピング由来の隣接候補 + 任意長オーバーラップ探索) は実データで平均 2.96 塩基という偶然の一致で unitig を接着していたため、「辺が張られる条件」そのものをここで固定する
    /// </remarks>
    public class UnitigGraphTests
    {
        #region 定数

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一方の unitig の末尾がもう一方の unitig の先頭へ伸びると辺が張られる
        /// </summary>
        [Fact]
        public void V_一方のunitigの末尾がもう一方の先頭へ伸びると辺が張られる()
        {
            // unitig A の末尾 k-1 塩基が unitig B の先頭 k-1 塩基と一致する構成
            // A の末尾 k-mer から 1 塩基伸ばすと、ちょうど B の先頭 k-mer になる
            const int l_k長 = 8;
            const string l_shared = "CGTTACA"; // k-1 = 7 塩基の重なり
            var l_先頭配列 = "GCTAAAGACAATTAC" + l_shared;      // 末尾が shared
            var l_中間配列 = l_shared + "GGATCCTTAGGCAAT";      // 先頭が shared

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            var l_bForward = ContigMaker.Get_頂点番号(2);

            Assert.Contains(l_bForward, l_グラフ.A_出辺[l_aForward]);

            // 逆鎖対称性: A→B があるなら B' →A' も存在しなければならない
            // これが崩れると順鎖側と逆鎖側で別々の経路が組まれ、同じ領域が
            // 2 通りに組み立てられてしまう
            Assert.Contains(l_aForward ^ 1, l_グラフ.A_出辺[l_bForward ^ 1]);

            // 入次数は双子の出次数で表せる
            Assert.Equal(l_グラフ.A_出辺[l_bForward ^ 1].Count, l_グラフ.Get_入次数(l_bForward));
        }

        /// <summary>
        /// k-1 で重ならない unitig 同士には辺が張られない
        /// </summary>
        [Fact]
        public void V_k_1で重ならないunitig同士には辺が張られない()
        {
            // 互いに無関係な 2 本
            // 偶然の短い一致があっても辺は張られてはならない
            // (旧実装はここで平均 3 塩基程度の一致による誤結合を作っていた)
            const int l_k長 = 8;
            const string l_先頭配列 = "GCTAAAGACAATTACATAA";
            const string l_中間配列 = "TTGACCTGAATCCGGTTCA";

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            for (var v = 2; v < l_グラフ.A_出辺.Count; v++)
            {
                Assert.Empty(l_グラフ.A_出辺[v]);
            }
        }

        /// <summary>
        /// 一致する k-mer が相手の先頭ではなく途中にある場合は、辺が張られない
        /// </summary>
        [Fact]
        public void V_一致するkmerが相手の先頭でなく途中にある場合は辺が張られない()
        {
            // B の「途中」に一致する k-mer があっても、k-1 オーバーラップでの
            // 連結はできないため辺を張ってはならない (Position != 0 を弾く条件)
            // A の末尾 k-mer を 1 塩基伸ばした k-mer が B の 3 塩基目から始まるよう構成する
            const int l_k長 = 8;
            const string l_junction = "ACGGATCA"; // A の末尾から伸ばして得られる k-mer
            var l_先頭配列 = "GCTAAAGACAATTAC" + l_junction[..(l_k長 - 1)]; // 末尾 k-1 が junction の先頭 k-1
            var l_中間配列 = "TT" + l_junction + "CCTTAGGCAAT";          // junction は B の位置 2 に現れる

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            var l_bForward = ContigMaker.Get_頂点番号(2);
            Assert.DoesNotContain(l_bForward, l_グラフ.A_出辺[l_aForward]);
        }

        /// <summary>
        /// 末尾が 2 本の異なる unitig へ伸びる場合は、両方の辺が記録される
        /// </summary>
        [Fact]
        public void V_末尾が2本の異なるunitigへ伸びる場合は両方の辺が記録される()
        {
            // 分岐: A の末尾から B へも C へも伸びられる構成
            // 辺は両方張られ、どちらを選ぶかはリード支持に委ねられる
            const int l_k長 = 8;
            const string l_shared = "CGTTACA";
            var l_先頭配列 = "GCTAAAGACAATTAC" + l_shared;
            var l_中間配列 = l_shared + "GGATCCTTAGGCAAT";
            var l_末尾配列 = l_shared + "TGATCCTTAGGCAAT"; // 分岐点の 1 塩基だけ B と異なる

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列, l_末尾配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_グラフ.A_出辺[l_aForward].Count);
            Assert.Contains(ContigMaker.Get_頂点番号(2), l_グラフ.A_出辺[l_aForward]);
            Assert.Contains(ContigMaker.Get_頂点番号(3), l_グラフ.A_出辺[l_aForward]);
        }

        /// <summary>
        /// 単純バブル (u から 2 本に分かれ、それぞれ 1 本の unitig を経て同じ w へ再合流する) で、リード支持の高い枝だけが経路として残ることを確認する
        /// </summary>
        /// <remarks>
        /// 結合の採用条件を相互一意にした結果、再合流点 w の入次数が 2 のままだと u から w へ至る経路が一切結合されなくなるため、この処理が無いとバブルのたびに contig が千切れる
        /// </remarks>
        [Fact]
        public void V_単純バブルはリード支持の高い枝だけを残し逆鎖も対称に除去する()
        {
            const int l_k長 = 8;
            // k=8 で全 unitig を通じて重複する正規化 k-mer が無いことを確認済みの構成
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_b1 = "TTACGCAAGGATCCTGCACGT"; // u の末尾 7 塩基 + 'A' で始まり、w の先頭 7 塩基で終わる
            const string l_b2 = "TTACGCACTTAGCATGCACGT"; // 分岐点の 1 塩基だけ b1 と異なる同長の枝
            const string l_終点 = "TGCACGTAAGGCTTACCA";

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b1V = ContigMaker.Get_頂点番号(2);
            var l_b2V = ContigMaker.Get_頂点番号(3);
            var l_wV = ContigMaker.Get_頂点番号(4);

            // 前提: バブル構造が実際に構築されている
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_wV));

            // b1 側にだけリード支持を与える
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, l_b1V)] = 40UL,
                [(l_uV, l_b2V)] = 3UL,
            };

            var l_popped = l_グラフ.V_除去_単純バブル(l_ユニティグ一覧, l_支持, l_k長);

            Assert.Equal(1, l_popped);
            Assert.Equal([l_b1V], l_グラフ.A_出辺[l_uV]);
            Assert.Equal([l_wV], l_グラフ.A_出辺[l_b1V]);
            Assert.Empty(l_グラフ.A_出辺[l_b2V]);

            // 逆鎖側も対称に取り除かれていること (片側だけ消すと順鎖と逆鎖で
            // 別々の経路が組まれてしまう)
            Assert.Equal(1, l_グラフ.Get_入次数(l_wV));
            Assert.DoesNotContain(l_b2V ^ 1, l_グラフ.A_出辺[l_wV ^ 1]);
            Assert.DoesNotContain(l_uV ^ 1, l_グラフ.A_出辺[l_b2V ^ 1]);
        }

        /// <summary>
        /// 長さが大きく異なる分岐は、同じ領域の別表現 (バブル) ではなく本物の分岐 (反復配列の出入口など) である可能性が高いため、支持の低い側であっても勝手に経路から外してはならない
        /// </summary>
        [Fact]
        public void V_長さが大きく異なる分岐は単純バブルとして除去しない()
        {
            const int l_k長 = 8;
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_b1 = "TTACGCAAGGATCCTGCACGT";
            // b2 は b1 より大幅に長い (長さ比が既定の閾値 1.5 を超える)
            const string l_b2 = "TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT";
            const string l_終点 = "TGCACGTAAGGCTTACCA";

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);

            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40UL,
                [(l_uV, ContigMaker.Get_頂点番号(3))] = 3UL,
            };

            var l_popped = l_グラフ.V_除去_単純バブル(l_ユニティグ一覧, l_支持, l_k長);

            Assert.Equal(0, l_popped);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
        }

        /// <summary>
        /// バブルの枝が単一 unitig とは限らない
        /// </summary>
        /// <remarks>
        /// 分岐の無い (排他的な) 2 本の unitig をまたぐ枝同士でも、1 本の経路として検出・比較できること
        /// </remarks>
        [Fact]
        public void V_unitig2本からなる鎖を1本の枝として扱う()
        {
            const int l_k長 = 8;
            const string l_uTail = "TTACGCA"; // u の末尾 7 塩基 (既存テストと共通)
            const string l_wHead = "TGCACGT"; // w の先頭 7 塩基 (既存テストと共通)
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_終点 = "TGCACGTAAGGCTTACCA";

            // 枝 1 ・枝 2 は完全に無関係な乱数配列にする (配列類似度の検証 (第 2 段階) は
            // 別のテストで確かめるため、ここでは第 2 段階を無効にして多段構造の
            // 検出・長さ帯の比較そのものだけを確かめる)
            var l_joint1 = V_生成_乱数配列(7, p_乱数種: 501);
            var l_joint2 = V_生成_乱数配列(7, p_乱数種: 502);
            var l_middleA1 = V_生成_乱数配列(10, p_乱数種: 511);
            var l_middleA2 = V_生成_乱数配列(10, p_乱数種: 521);
            var l_middleB1 = V_生成_乱数配列(10, p_乱数種: 512);
            var l_middleB2 = V_生成_乱数配列(10, p_乱数種: 522);

            var l_b1a = l_uTail + l_middleA1 + l_joint1;
            var l_b1b = l_joint1 + l_middleB1 + l_wHead;
            var l_b2a = l_uTail + l_middleA2 + l_joint2;
            var l_b2b = l_joint2 + l_middleB2 + l_wHead;

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1a, l_b1b, l_b2a, l_b2b, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b1aV = ContigMaker.Get_頂点番号(2);
            var l_b1bV = ContigMaker.Get_頂点番号(3);
            var l_b2aV = ContigMaker.Get_頂点番号(4);
            var l_b2bV = ContigMaker.Get_頂点番号(5);
            var l_wV = ContigMaker.Get_頂点番号(6);

            // 前提: それぞれの枝が 2 unitig の分岐無しの鎖になっている
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
            Assert.Equal([l_b1bV], l_グラフ.A_出辺[l_b1aV]);
            Assert.Equal([l_b2bV], l_グラフ.A_出辺[l_b2aV]);
            Assert.Equal(2, l_グラフ.Get_入次数(l_wV));

            // 枝 1 にだけリード支持を与える
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, l_b1aV)] = 40UL,
                [(l_uV, l_b2aV)] = 3UL,
            };

            // 配列類似度の検証 (第 2 段階) は別のテストで確かめる
            // ここでは多段構造の検出・長さ帯の比較だけを見たいので無効にする
            var l_popped = l_グラフ.V_除去_単純バブル(l_ユニティグ一覧, l_支持, l_k長, p_類似度の下限: 0.0D);

            Assert.Equal(1, l_popped);
            // 枝 1 (2 unitig とも) は生き残り、枝 2 (2 unitig とも) は取り除かれる
            Assert.Equal([l_b1aV], l_グラフ.A_出辺[l_uV]);
            Assert.Equal([l_wV], l_グラフ.A_出辺[l_b1bV]);
            Assert.Empty(l_グラフ.A_出辺[l_b2aV]);
            Assert.Empty(l_グラフ.A_出辺[l_b2bV]);
            Assert.Equal(1, l_グラフ.Get_入次数(l_wV));
        }

        /// <summary>
        /// careful_bubble: 除去された側の経路の配列を、引き継ぎ先へ集められること
        /// </summary>
        /// <remarks>
        /// 「この k では敗者と判断したが、次の k は自分の証拠で判断し直せる」ため、配列自体は捨てない
        /// </remarks>
        [Fact]
        public void V_引き継ぎ先を指定すると除去された側の配列を集められる()
        {
            const int l_k長 = 8;
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_b1 = "TTACGCAAGGATCCTGCACGT";
            const string l_b2 = "TTACGCACTTAGCATGCACGT";
            const string l_終点 = "TGCACGTAAGGCTTACCA";

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b2V = ContigMaker.Get_頂点番号(3);

            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40UL,
                [(l_uV, l_b2V)] = 3UL,
            };

            List<string> l_carryOver = [];
            var l_popped = l_グラフ.V_除去_単純バブル(l_ユニティグ一覧, l_支持, l_k長, l_carryOver);

            Assert.Equal(1, l_popped);
            Assert.Equal([l_b2], l_carryOver);
        }

        /// <summary>
        /// 長さが揃っていても配列がまるで違う (たまたま長さが一致しただけの別の反復など) 場合は、同じ領域の別表現とは言えないため触らない
        /// </summary>
        [Fact]
        public void V_長さは同じでも配列が無関係な分岐は除去しない()
        {
            const int l_k長 = 8;
            const string l_uTail = "TTACGCA";
            const string l_wHead = "TGCACGT";
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_終点 = "TGCACGTAAGGCTTACCA";

            // 長さは完全に一致するが、中身は無関係な乱数配列
            var l_b1 = l_uTail + V_生成_乱数配列(20, p_乱数種: 601) + l_wHead;
            var l_b2 = l_uTail + V_生成_乱数配列(20, p_乱数種: 602) + l_wHead;

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);

            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40UL,
                [(l_uV, ContigMaker.Get_頂点番号(3))] = 3UL,
            };

            var l_popped = l_グラフ.V_除去_単純バブル(l_ユニティグ一覧, l_支持, l_k長);

            Assert.Equal(0, l_popped);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
        }

        /// <summary>
        /// 自己ループの辺は作らない
        /// </summary>
        [Fact]
        public void V_自己ループの辺は作らない()
        {
            // 自己ループを辺として持つと walk が同じ unitig を無限に伸ばしうるため、
            // 構築段階で除外していることを確認する
            const int l_k長 = 8;
            const string l_反復単位 = "ACGGATCT";
            var l_先頭配列 = l_反復単位 + "GCTAAAGA" + l_反復単位[..(l_k長 - 1)];

            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, l_k長, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            Assert.DoesNotContain(l_aForward, l_グラフ.A_出辺[l_aForward]);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ContigMaker のコンストラクタと同じ規則で kmerDict を組み立てる
        /// </summary>
        /// <remarks>
        /// 添字 2 u が unitig u の順鎖、2 u+1 が逆鎖
        /// </remarks>
        /// <param name="p_kmer長">k-mer 長</param>
        /// <param name="p_ユニティグ群">登録する unitig 配列群</param>
        /// <returns>ユニティグ一覧と kmer 辞書</returns>
        private static (List<string> A_ユニティグ一覧, Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> A_kmer辞書) V_構築(int p_kmer長, params string[] p_ユニティグ群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmer長, A_スレッド数 = 1 };
            List<string> l_ユニティグ一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> l_kmer辞書 = [];

            var l_ID = 1;
            foreach (var l_配列 in p_ユニティグ群)
            {
                l_ユニティグ一覧.Add(l_配列);
                l_ユニティグ一覧.Add(Util.V_逆相補(l_配列));

                for (var i = p_kmer長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - p_kmer長;
                    var l_キー = new KmerKey(l_配列.AsSpan(l_開始位置, p_kmer長));
                    var l_逆相補キー = l_キー.Get_逆相補();
                    var l_逆鎖開始位置 = l_配列.Length - i;
                    V_登録_kmer(l_kmer辞書, l_キー, l_ID, l_開始位置);
                    V_登録_kmer(l_kmer辞書, l_逆相補キー, -l_ID, l_逆鎖開始位置);
                }
                l_ID++;
            }
            return (l_ユニティグ一覧, l_kmer辞書);
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_キー">登録する k-mer</param>
        /// <param name="p_ID">ユニティグ ID</param>
        /// <param name="p_位置">ユニティグ内の開始位置</param>
        private static void V_登録_kmer(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存値))
            {
                if (l_既存値.Item1 is 曖昧kmer番号 || l_既存値.Item1 == p_ID)
                {
                    return;
                }
                p_辞書[p_キー] = (曖昧kmer番号, 0);
                return;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_乱数種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_長さ, int p_乱数種)
        {
            var l_乱数生成器 = new Random(p_乱数種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数生成器.Next(4)]));
        }

        #endregion

    }
}
