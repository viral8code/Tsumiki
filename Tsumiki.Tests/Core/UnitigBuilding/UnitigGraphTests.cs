using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig 間の隣接を de Bruijn グラフから厳密に構築する UnitigGraph の検証
    /// </summary>
    public class UnitigGraphTests
    {
        #region 定数

        /// <summary>
        /// 塩基配列 CGTTACA
        /// </summary>
        private const string C_塩基配列_CGTTACA = "CGTTACA";

        /// <summary>
        /// 塩基配列 GCTAAAGACAATTAC
        /// </summary>
        private const string C_塩基配列_GCTAAAGACAATTAC = "GCTAAAGACAATTAC";

        /// <summary>
        /// 塩基配列 GGATCCTTAGGCAAT
        /// </summary>
        private const string C_塩基配列_GGATCCTTAGGCAAT = "GGATCCTTAGGCAAT";

        /// <summary>
        /// 塩基配列 GCTAAAGACAATTACATAA
        /// </summary>
        private const string C_塩基配列_GCTAAAGACAATTACATAA = "GCTAAAGACAATTACATAA";

        /// <summary>
        /// 塩基配列 TTGACCTGAATCCGGTTCA
        /// </summary>
        private const string C_塩基配列_TTGACCTGAATCCGGTTCA = "TTGACCTGAATCCGGTTCA";

        /// <summary>
        /// 塩基配列 ACGGATCA
        /// </summary>
        private const string C_塩基配列_ACGGATCA = "ACGGATCA";

        /// <summary>
        /// 塩基配列 TT
        /// </summary>
        private const string C_塩基配列_TT = "TT";

        /// <summary>
        /// 塩基配列 CCTTAGGCAAT
        /// </summary>
        private const string C_塩基配列_CCTTAGGCAAT = "CCTTAGGCAAT";

        /// <summary>
        /// 塩基配列 TGATCCTTAGGCAAT
        /// </summary>
        private const string C_塩基配列_TGATCCTTAGGCAAT = "TGATCCTTAGGCAAT";

        /// <summary>
        /// 塩基配列 GCTAAAGACAATTACGCA
        /// </summary>
        private const string C_塩基配列_GCTAAAGACAATTACGCA = "GCTAAAGACAATTACGCA";

        /// <summary>
        /// 塩基配列 TTACGCAAGGATCCTGCACGT
        /// </summary>
        private const string C_塩基配列_TTACGCAAGGATCCTGCACGT = "TTACGCAAGGATCCTGCACGT";

        /// <summary>
        /// 塩基配列 TTACGCACTTAGCATGCACGT
        /// </summary>
        private const string C_塩基配列_TTACGCACTTAGCATGCACGT = "TTACGCACTTAGCATGCACGT";

        /// <summary>
        /// 塩基配列 TGCACGTAAGGCTTACCA
        /// </summary>
        private const string C_塩基配列_TGCACGTAAGGCTTACCA = "TGCACGTAAGGCTTACCA";

        /// <summary>
        /// 塩基配列 TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT
        /// </summary>
        private const string C_塩基配列_TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT = "TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT";

        /// <summary>
        /// 塩基配列 TTACGCA
        /// </summary>
        private const string C_塩基配列_TTACGCA = "TTACGCA";

        /// <summary>
        /// 塩基配列 TGCACGT
        /// </summary>
        private const string C_塩基配列_TGCACGT = "TGCACGT";

        /// <summary>
        /// 塩基配列 ACGGATCT
        /// </summary>
        private const string C_塩基配列_ACGGATCT = "ACGGATCT";

        /// <summary>
        /// 塩基配列 GCTAAAGA
        /// </summary>
        private const string C_塩基配列_GCTAAAGA = "GCTAAAGA";

        /// <summary>
        /// 塩基配列 GCATTGAGTCCA
        /// </summary>
        private const string C_塩基配列_GCATTGAGTCCA = "GCATTGAGTCCA";

        /// <summary>
        /// 塩基配列 TGGA
        /// </summary>
        private const string C_塩基配列_TGGA = "TGGA";

        /// <summary>
        /// 塩基配列 CCGGTTGAG
        /// </summary>
        private const string C_塩基配列_CCGGTTGAG = "CCGGTTGAG";

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        /// <summary>
        /// 項目 C
        /// </summary>
        private const string C_項目_C = "C";

        /// <summary>
        /// 項目 G
        /// </summary>
        private const string C_項目_G = "G";

        /// <summary>
        /// 項目 T
        /// </summary>
        private const string C_項目_T = "T";

        /// <summary>
        /// 空白
        /// </summary>
        private const string C_空白 = " ";

        /// <summary>
        /// 区切り
        /// </summary>
        private const string C_区切り = ",";

        /// <summary>
        /// 塩基配列 GGA
        /// </summary>
        private const string C_塩基配列_GGA = "GGA";

        /// <summary>
        /// 塩基配列 CCA
        /// </summary>
        private const string C_塩基配列_CCA = "CCA";

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int C_曖昧kmer番号 = int.MinValue;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一方の unitig の末尾がもう一方の unitig の先頭へ伸びると辺が張られる
        /// </summary>
        [Fact]
        public void V_一方のunitigの末尾がもう一方の先頭へ伸びると辺が張られる()
        {
            const int l_k長 = 8;
            const string l_shared = C_塩基配列_CGTTACA;
            var l_先頭配列 = C_塩基配列_GCTAAAGACAATTAC + l_shared;
            var l_中間配列 = l_shared + C_塩基配列_GGATCCTTAGGCAAT;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_aForward = ContigMaker.Get_頂点番号(1);
            var l_bForward = ContigMaker.Get_頂点番号(2);
            Assert.Contains(l_bForward, l_グラフ.A_出辺[l_aForward]);
            Assert.Contains(l_aForward ^ 1, l_グラフ.A_出辺[l_bForward ^ 1]);
            Assert.Equal(l_グラフ.A_出辺[l_bForward ^ 1].Count, l_グラフ.Get_入次数(l_bForward));
        }

        /// <summary>
        /// k-1 で重ならない unitig 同士には辺が張られない
        /// </summary>
        [Fact]
        public void V_k_1で重ならないunitig同士には辺が張られない()
        {
            const int l_k長 = 8;
            const string l_先頭配列 = C_塩基配列_GCTAAAGACAATTACATAA;
            const string l_中間配列 = C_塩基配列_TTGACCTGAATCCGGTTCA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
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
            const int l_k長 = 8;
            const string l_junction = C_塩基配列_ACGGATCA;
            var l_先頭配列 = C_塩基配列_GCTAAAGACAATTAC + l_junction[..(l_k長 - 1)];
            var l_中間配列 = C_塩基配列_TT + l_junction + C_塩基配列_CCTTAGGCAAT;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
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
            const int l_k長 = 8;
            const string l_shared = C_塩基配列_CGTTACA;
            var l_先頭配列 = C_塩基配列_GCTAAAGACAATTAC + l_shared;
            var l_中間配列 = l_shared + C_塩基配列_GGATCCTTAGGCAAT;
            var l_末尾配列 = l_shared + C_塩基配列_TGATCCTTAGGCAAT;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_中間配列, l_末尾配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_aForward = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_グラフ.A_出辺[l_aForward].Count);
            Assert.Contains(ContigMaker.Get_頂点番号(2), l_グラフ.A_出辺[l_aForward]);
            Assert.Contains(ContigMaker.Get_頂点番号(3), l_グラフ.A_出辺[l_aForward]);
        }

        /// <summary>
        /// 単純バブル (u から 2 本に分かれ、それぞれ 1 本の unitig を経て同じ w へ再合流する) で、リード支持の高い枝だけが経路として残ることを確認する
        /// </summary>
        [Fact]
        public void V_単純バブルはリード支持の高い枝だけを残し逆鎖も対称に除去する()
        {
            const int l_k長 = 8;
            const string l_u = C_塩基配列_GCTAAAGACAATTACGCA;
            const string l_b1 = C_塩基配列_TTACGCAAGGATCCTGCACGT;
            const string l_b2 = C_塩基配列_TTACGCACTTAGCATGCACGT;
            const string l_終点 = C_塩基配列_TGCACGTAAGGCTTACCA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b1V = ContigMaker.Get_頂点番号(2);
            var l_b2V = ContigMaker.Get_頂点番号(3);
            var l_wV = ContigMaker.Get_頂点番号(4);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_wV));
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, l_b1V)] = 40UL,
                [(l_uV, l_b2V)] = 3UL,
            };
            var l_popped = l_グラフ.V_除去_単純バブル(l_unitig一覧, l_支持, l_k長);
            Assert.Equal(1, l_popped);
            Assert.Equal([l_b1V], l_グラフ.A_出辺[l_uV]);
            Assert.Equal([l_wV], l_グラフ.A_出辺[l_b1V]);
            Assert.Empty(l_グラフ.A_出辺[l_b2V]);
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
            const string l_u = C_塩基配列_GCTAAAGACAATTACGCA;
            const string l_b1 = C_塩基配列_TTACGCAAGGATCCTGCACGT;
            const string l_b2 = C_塩基配列_TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT;
            const string l_終点 = C_塩基配列_TGCACGTAAGGCTTACCA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40UL,
                [(l_uV, ContigMaker.Get_頂点番号(3))] = 3UL,
            };
            var l_popped = l_グラフ.V_除去_単純バブル(l_unitig一覧, l_支持, l_k長);
            Assert.Equal(0, l_popped);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
        }

        /// <summary>
        /// バブルの枝が単一 unitig とは限らない
        /// </summary>
        [Fact]
        public void V_unitig2本からなる鎖を1本の枝として扱う()
        {
            const int l_k長 = 8;
            const string l_uTail = C_塩基配列_TTACGCA;
            const string l_wHead = C_塩基配列_TGCACGT;
            const string l_u = C_塩基配列_GCTAAAGACAATTACGCA;
            const string l_終点 = C_塩基配列_TGCACGTAAGGCTTACCA;
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
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1a, l_b1b, l_b2a, l_b2b, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b1aV = ContigMaker.Get_頂点番号(2);
            var l_b1bV = ContigMaker.Get_頂点番号(3);
            var l_b2aV = ContigMaker.Get_頂点番号(4);
            var l_b2bV = ContigMaker.Get_頂点番号(5);
            var l_wV = ContigMaker.Get_頂点番号(6);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
            Assert.Equal([l_b1bV], l_グラフ.A_出辺[l_b1aV]);
            Assert.Equal([l_b2bV], l_グラフ.A_出辺[l_b2aV]);
            Assert.Equal(2, l_グラフ.Get_入次数(l_wV));
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, l_b1aV)] = 40UL,
                [(l_uV, l_b2aV)] = 3UL,
            };
            var l_popped = l_グラフ.V_除去_単純バブル(l_unitig一覧, l_支持, l_k長, p_類似度の下限: 0.0D);
            Assert.Equal(1, l_popped);
            Assert.Equal([l_b1aV], l_グラフ.A_出辺[l_uV]);
            Assert.Equal([l_wV], l_グラフ.A_出辺[l_b1bV]);
            Assert.Empty(l_グラフ.A_出辺[l_b2aV]);
            Assert.Empty(l_グラフ.A_出辺[l_b2bV]);
            Assert.Equal(1, l_グラフ.Get_入次数(l_wV));
        }

        /// <summary>
        /// careful_bubble: 除去された側の経路の配列を、引き継ぎ先へ集められること
        /// </summary>
        [Fact]
        public void V_引き継ぎ先を指定すると除去された側の配列を集められる()
        {
            const int l_k長 = 8;
            const string l_u = C_塩基配列_GCTAAAGACAATTACGCA;
            const string l_b1 = C_塩基配列_TTACGCAAGGATCCTGCACGT;
            const string l_b2 = C_塩基配列_TTACGCACTTAGCATGCACGT;
            const string l_終点 = C_塩基配列_TGCACGTAAGGCTTACCA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b2V = ContigMaker.Get_頂点番号(3);
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40UL,
                [(l_uV, l_b2V)] = 3UL,
            };
            List<string> l_carryOver = [];
            var l_popped = l_グラフ.V_除去_単純バブル(l_unitig一覧, l_支持, l_k長, l_carryOver);
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
            const string l_uTail = C_塩基配列_TTACGCA;
            const string l_wHead = C_塩基配列_TGCACGT;
            const string l_u = C_塩基配列_GCTAAAGACAATTACGCA;
            const string l_終点 = C_塩基配列_TGCACGTAAGGCTTACCA;
            var l_b1 = l_uTail + V_生成_乱数配列(20, p_乱数種: 601) + l_wHead;
            var l_b2 = l_uTail + V_生成_乱数配列(20, p_乱数種: 602) + l_wHead;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_u, l_b1, l_b2, l_終点);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
            Dictionary<(int, int), ulong> l_支持 = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40UL,
                [(l_uV, ContigMaker.Get_頂点番号(3))] = 3UL,
            };
            var l_popped = l_グラフ.V_除去_単純バブル(l_unitig一覧, l_支持, l_k長);
            Assert.Equal(0, l_popped);
            Assert.Equal(2, l_グラフ.A_出辺[l_uV].Count);
        }

        /// <summary>
        /// 自己ループの辺は作らない
        /// </summary>
        [Fact]
        public void V_自己ループの辺は作らない()
        {
            const int l_k長 = 8;
            const string l_反復単位 = C_塩基配列_ACGGATCT;
            var l_先頭配列 = l_反復単位 + C_塩基配列_GCTAAAGA + l_反復単位[..(l_k長 - 1)];
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_aForward = ContigMaker.Get_頂点番号(1);
            Assert.DoesNotContain(l_aForward, l_グラフ.A_出辺[l_aForward]);
        }

        /// <summary>
        /// 残る続きが一意な分岐から、行き止まりの短い枝を外す
        /// </summary>
        [Fact]
        public void V_行き止まりの短い枝を外して続きを一意にする()
        {
            const int l_k長 = 8;
            const string l_先頭配列 = C_塩基配列_GCTAAAGACAATTAC;
            var l_続き配列 = l_先頭配列[^(l_k長 - 1)..] + C_塩基配列_GCATTGAGTCCA;
            var l_枝配列 = l_先頭配列[^(l_k長 - 1)..] + C_塩基配列_TGGA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_続き配列, l_枝配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_先頭 = ContigMaker.Get_頂点番号(1);
            var l_続き = ContigMaker.Get_頂点番号(2);
            var l_枝 = ContigMaker.Get_頂点番号(3);
            Assert.Equal(2, l_グラフ.A_出辺[l_先頭].Count);
            var l_外した数 = l_グラフ.V_除去_行き止まり枝(l_unitig一覧, p_枝長の上限: 12);
            Assert.Equal(1, l_外した数);
            Assert.Equal([l_続き], l_グラフ.A_出辺[l_先頭]);
            Assert.Empty(l_グラフ.A_出辺[l_枝 ^ 1]);
        }

        /// <summary>
        /// 続きの側に別の入口がある分岐は、反復の別コピーへの辺かもしれないので触らない
        /// </summary>
        [Fact]
        public void V_続きに別の入口があれば行き止まりの枝を外さない()
        {
            const int l_k長 = 8;
            const string l_先頭配列 = C_塩基配列_GCTAAAGACAATTAC;
            var l_続き配列 = l_先頭配列[^(l_k長 - 1)..] + C_塩基配列_GCATTGAGTCCA;
            var l_枝配列 = l_先頭配列[^(l_k長 - 1)..] + C_塩基配列_TGGA;
            var l_別の入口 = C_塩基配列_CCGGTTGAG + l_先頭配列[^(l_k長 - 1)..];
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_続き配列, l_枝配列, l_別の入口);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_外した数 = l_グラフ.V_除去_行き止まり枝(l_unitig一覧, p_枝長の上限: 12);
            Assert.Equal(0, l_外した数);
            Assert.Equal(2, l_グラフ.A_出辺[ContigMaker.Get_頂点番号(1)].Count);
        }

        /// <summary>
        /// 反復に隣接する unitig にはペアが無くても、その先の一意な鎖を結ぶペアで反復を解く
        /// </summary>
        [Fact]
        public void V_隣接unitigにペアが無くても足場のペアで短い反復を解く()
        {
            const int l_k長 = 8;
            var l_反復 = V_生成_乱数配列(16, p_乱数種: 1001);
            var l_入1 = V_生成_乱数配列(11, p_乱数種: 1002) + C_項目_A + l_反復[..(l_k長 - 1)];
            var l_入2 = V_生成_乱数配列(11, p_乱数種: 1003) + C_項目_C + l_反復[..(l_k長 - 1)];
            var l_上流1 = V_生成_乱数配列(12, p_乱数種: 1004) + l_入1[..(l_k長 - 1)];
            var l_上流2 = V_生成_乱数配列(12, p_乱数種: 1005) + l_入2[..(l_k長 - 1)];
            var l_出1 = l_反復[^(l_k長 - 1)..] + C_項目_G + V_生成_乱数配列(11, p_乱数種: 1006);
            var l_出2 = l_反復[^(l_k長 - 1)..] + C_項目_T + V_生成_乱数配列(11, p_乱数種: 1007);
            var l_下流1 = l_出1[^(l_k長 - 1)..] + V_生成_乱数配列(12, p_乱数種: 1008);
            var l_下流2 = l_出2[^(l_k長 - 1)..] + V_生成_乱数配列(12, p_乱数種: 1009);
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_上流1, l_上流2, l_入1, l_入2, l_反復, l_出1, l_出2, l_下流1, l_下流2);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_反復頂点 = ContigMaker.Get_頂点番号(5);
            var l_隣接 = string.Join(C_空白, Enumerable.Range(2, l_グラフ.A_出辺.Count - 2).Select(v => $"{v}->[{string.Join(C_区切り, l_グラフ.A_出辺[v])}]"));
            Assert.True(l_グラフ.A_出辺[l_反復頂点].Count == 2 && l_グラフ.Get_入次数(l_反復頂点) == 2, l_隣接);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(ContigMaker.Get_頂点番号(1), ContigMaker.Get_頂点番号(8))] = 30UL,
                [(ContigMaker.Get_頂点番号(2), ContigMaker.Get_頂点番号(9))] = 30UL,
            };
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, [], l_ペア連結, p_反復長の上限: 100, p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.Equal(1, l_解決数);
            var l_入1から = Assert.Single(l_グラフ.A_出辺[ContigMaker.Get_頂点番号(3)]);
            var l_入2から = Assert.Single(l_グラフ.A_出辺[ContigMaker.Get_頂点番号(4)]);
            Assert.NotEqual(l_入1から, l_入2から);
            Assert.Equal([ContigMaker.Get_頂点番号(6)], l_グラフ.A_出辺[l_入1から]);
            Assert.Equal([ContigMaker.Get_頂点番号(7)], l_グラフ.A_出辺[l_入2から]);
        }

        /// <summary>
        /// 入口も出口も 2 本ある頂点は、コピー数が 1 と推定されていても通り抜けさせない
        /// </summary>
        [Fact]
        public void Is通過可能_入口と出口が両方分岐する頂点はコピー数1でも通さない()
        {
            const int l_k長 = 8;
            var l_反復 = V_生成_乱数配列(16, p_乱数種: 1101);
            var l_入1 = V_生成_乱数配列(11, p_乱数種: 1102) + C_項目_A + l_反復[..(l_k長 - 1)];
            var l_入2 = V_生成_乱数配列(11, p_乱数種: 1103) + C_項目_C + l_反復[..(l_k長 - 1)];
            var l_出1 = l_反復[^(l_k長 - 1)..] + C_項目_G + V_生成_乱数配列(11, p_乱数種: 1104);
            var l_出2 = l_反復[^(l_k長 - 1)..] + C_項目_T + V_生成_乱数配列(11, p_乱数種: 1105);
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_入1, l_入2, l_反復, l_出1, l_出2);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_反復頂点 = ContigMaker.Get_頂点番号(3);
            Assert.True(l_グラフ.A_出辺[l_反復頂点].Count == 2 && l_グラフ.Get_入次数(l_反復頂点) == 2);
            Dictionary<int, int> l_コピー数 = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            Assert.False(l_グラフ.Is通過可能(l_コピー数, l_反復頂点));
            Assert.False(l_グラフ.Is通過可能(l_コピー数, l_反復頂点 ^ 1));
            Assert.True(l_グラフ.Is通過可能(l_コピー数, ContigMaker.Get_頂点番号(1)));
        }

        /// <summary>
        /// 入口 2・出口 2 の片方がこの頂点へ戻る短い脇道だけなら通り抜けを許し、脇道が長ければ許さない
        /// </summary>
        /// <param name="p_脇道の中身"></param>
        /// <param name="p_期待"></param>
        [Theory]
        [InlineData("TC", true)]
        [InlineData("TGACGGTAC", false)]
        public void Is通過可能_短い脇道だけの分岐は通す(string p_脇道の中身, bool p_期待)
        {
            const int l_k長 = 8;
            var l_反復 = V_生成_乱数配列(16, p_乱数種: 1201);
            var l_入 = V_生成_乱数配列(11, p_乱数種: 1202) + C_項目_A + l_反復[..(l_k長 - 1)];
            var l_出 = l_反復[^(l_k長 - 1)..] + C_項目_G + V_生成_乱数配列(11, p_乱数種: 1203);
            var l_脇道 = l_反復[^(l_k長 - 1)..] + p_脇道の中身 + l_反復[..(l_k長 - 1)];
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_入, l_反復, l_出, l_脇道);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_反復頂点 = ContigMaker.Get_頂点番号(2);
            var l_隣接 = string.Join(C_空白, Enumerable.Range(2, l_グラフ.A_出辺.Count - 2).Select(v => $"{v}->[{string.Join(C_区切り, l_グラフ.A_出辺[v])}]"));
            Assert.True(l_グラフ.A_出辺[l_反復頂点].Count == 2 && l_グラフ.Get_入次数(l_反復頂点) == 2, l_隣接);
            Dictionary<int, int> l_コピー数 = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1
            };
            Assert.Equal(p_期待, l_グラフ.Is通過可能(l_コピー数, l_反復頂点, l_unitig一覧));
            Assert.Equal(p_期待, l_グラフ.Is通過可能(l_コピー数, l_反復頂点 ^ 1, l_unitig一覧));
            Assert.False(l_グラフ.Is通過可能(l_コピー数, l_反復頂点));
        }

        /// <summary>
        /// 入口も出口も 2 本ある反復の頂点では、行き止まりの枝は外すが walk に通り抜けさせない
        /// </summary>
        [Fact]
        public void V_反復の頂点で枝を外したら通り抜けさせない()
        {
            const int l_k長 = 8;
            var l_反復 = V_生成_乱数配列(16, p_乱数種: 1301);
            var l_入1 = V_生成_乱数配列(11, p_乱数種: 1302) + C_項目_A + l_反復[..(l_k長 - 1)];
            var l_入2 = C_塩基配列_GGA + C_項目_C + l_反復[..(l_k長 - 1)];
            var l_出1 = l_反復[^(l_k長 - 1)..] + C_項目_G + V_生成_乱数配列(11, p_乱数種: 1304);
            var l_出2 = l_反復[^(l_k長 - 1)..] + C_項目_T + C_塩基配列_CCA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_入1, l_入2, l_反復, l_出1, l_出2);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_反復頂点 = ContigMaker.Get_頂点番号(3);
            Assert.True(l_グラフ.A_出辺[l_反復頂点].Count == 2 && l_グラフ.Get_入次数(l_反復頂点) == 2);
            Assert.Equal(2, l_グラフ.V_除去_行き止まり枝(l_unitig一覧, p_枝長の上限: 12));
            var l_続き = Assert.Single(l_グラフ.A_出辺[l_反復頂点]);
            Assert.Equal(1, l_グラフ.Get_入次数(l_反復頂点));
            Dictionary<int, int> l_コピー数 = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            Assert.False(l_グラフ.Is通過可能(l_コピー数, l_反復頂点, l_unitig一覧));
            Assert.False(l_グラフ.Is通過可能(l_コピー数, l_反復頂点 ^ 1, l_unitig一覧));
            Assert.False(l_グラフ.Is構造上一意な辺(l_反復頂点, l_続き));
        }

        /// <summary>
        /// 入口が 2 本残る反復の頂点で出口側の枝だけを外したなら、頂点から続きへの結合は許すが、頂点へ入る結合は塞いだままにする
        /// </summary>
        [Fact]
        public void V_入口が2本残る反復の頂点は続きとだけ繋げる()
        {
            const int l_k長 = 8;
            var l_反復 = V_生成_乱数配列(16, p_乱数種: 1311);
            var l_入1 = V_生成_乱数配列(11, p_乱数種: 1312) + C_項目_A + l_反復[..(l_k長 - 1)];
            var l_入2 = V_生成_乱数配列(11, p_乱数種: 1313) + C_項目_C + l_反復[..(l_k長 - 1)];
            var l_出1 = l_反復[^(l_k長 - 1)..] + C_項目_G + V_生成_乱数配列(11, p_乱数種: 1314);
            var l_出2 = l_反復[^(l_k長 - 1)..] + C_項目_T + C_塩基配列_CCA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_入1, l_入2, l_反復, l_出1, l_出2);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            var l_反復頂点 = ContigMaker.Get_頂点番号(3);
            Assert.True(l_グラフ.A_出辺[l_反復頂点].Count == 2 && l_グラフ.Get_入次数(l_反復頂点) == 2);
            Assert.Equal(1, l_グラフ.V_除去_行き止まり枝(l_unitig一覧, p_枝長の上限: 12));
            var l_続き = Assert.Single(l_グラフ.A_出辺[l_反復頂点]);
            Assert.Equal(2, l_グラフ.Get_入次数(l_反復頂点));
            Assert.True(l_グラフ.Is構造上一意な辺(l_反復頂点, l_続き));
            Dictionary<int, int> l_コピー数 = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            Assert.False(l_グラフ.Is通過可能(l_コピー数, l_反復頂点 ^ 1, l_unitig一覧));
            foreach (var l_入口 in l_グラフ.A_出辺[l_反復頂点 ^ 1])
            {
                Assert.False(l_グラフ.Is構造上一意な辺(l_入口 ^ 1, l_反復頂点));
            }
        }

        /// <summary>
        /// 上限より長い行き止まりの枝は外さない
        /// </summary>
        [Fact]
        public void V_長い行き止まりの枝は外さない()
        {
            const int l_k長 = 8;
            const string l_先頭配列 = C_塩基配列_GCTAAAGACAATTAC;
            var l_続き配列 = l_先頭配列[^(l_k長 - 1)..] + C_塩基配列_GCATTGAGTCCA;
            var l_枝配列 = l_先頭配列[^(l_k長 - 1)..] + C_塩基配列_TGGA;
            var (l_unitig一覧, l_kmer辞書) = V_構築(l_k長, l_先頭配列, l_続き配列, l_枝配列);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, l_k長, C_曖昧kmer番号);
            Assert.Equal(0, l_グラフ.V_除去_行き止まり枝(l_unitig一覧, p_枝長の上限: 5));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ContigMaker のコンストラクタと同じ規則で kmerDict を組み立てる
        /// </summary>
        /// <param name="p_kmer長">k-mer 長</param>
        /// <param name="p_unitig群">登録する unitig 配列群</param>
        /// <returns>unitig 一覧と kmer 辞書</returns>
        private static (List<string> A_unitig一覧, Dictionary<KmerKey, (int A_unitigID, int A_位置)> A_kmer辞書) V_構築(int p_kmer長, params string[] p_unitig群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = p_kmer長,
                A_スレッド数 = 1
            };
            List<string> l_unitig一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_unitigID, int A_位置)> l_kmer辞書 = [];
            var l_ID = 1;
            foreach (var l_配列 in p_unitig群)
            {
                l_unitig一覧.Add(l_配列);
                l_unitig一覧.Add(Util.V_逆相補(l_配列));
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

            return (l_unitig一覧, l_kmer辞書);
        }

        /// <summary>
        /// k-mer を、それが載る unitig と開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_キー">登録する k-mer</param>
        /// <param name="p_ID">unitig ID</param>
        /// <param name="p_位置">unitig 内の開始位置</param>
        private static void V_登録_kmer(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存値))
            {
                if (l_既存値.Item1 is C_曖昧kmer番号 || l_既存値.Item1 == p_ID)
                {
                    return;
                }

                p_辞書[p_キー] = (C_曖昧kmer番号, 0);
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
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数生成器.Next(4)]));
        }

        #endregion
    }
}
