using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig 間の隣接を de Bruijn グラフから厳密に構築する UnitigGraph の検証
    /// </summary>
    /// <remarks>
    /// 旧実装 (リードマッピング由来の隣接候補 + 任意長オーバーラップ探索) は
    /// 実データで平均 2.96 塩基という偶然の一致で unitig を接着していたため、
    /// 「辺が張られる条件」そのものをここで固定する
    /// </remarks>
    public class UnitigGraphTests
    {
        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// ContigMaker のコンストラクタと同じ規則で kmerDict を組み立てる
        /// </summary>
        /// <remarks>
        /// 添字 2 u が unitig u の順鎖、2 u+1 が逆鎖
        /// </remarks>
        /// <param name="p_kmerLength">k-mer 長</param>
        /// <param name="p_unitigs">登録する unitig 配列群</param>
        /// <returns>ユニティグ一覧と kmer 辞書</returns>
        private static (List<string> A_ユニティグ一覧, Dictionary<KmerKey, (int UnitigId, int Position)> A_kmer辞書) V_構築(
            int p_kmerLength,
            params string[] p_unitigs)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmerLength, A_スレッド数 = 1 };
            List<string> l_unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> l_kmerDict = [];

            var l_id = 1;
            foreach (var l_seq in p_unitigs)
            {
                l_unitigList.Add(l_seq);
                l_unitigList.Add(Util.V_逆相補(l_seq));

                for (var i = p_kmerLength; i <= l_seq.Length; i++)
                {
                    var l_startPos = i - p_kmerLength;
                    var l_key = new KmerKey(l_seq.AsSpan(l_startPos, p_kmerLength));
                    var l_revKey = l_key.Get_逆相補();
                    var l_revStartPos = l_seq.Length - i;
                    V_登録_kmer(l_kmerDict, l_key, l_id, l_startPos);
                    V_登録_kmer(l_kmerDict, l_revKey, -l_id, l_revStartPos);
                }
                l_id++;
            }
            return (l_unitigList, l_kmerDict);
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_dict">登録先の辞書</param>
        /// <param name="p_key">登録する k-mer</param>
        /// <param name="p_id">ユニティグ ID</param>
        /// <param name="p_position">ユニティグ内の開始位置</param>
        private static void V_登録_kmer(Dictionary<KmerKey, (int, int)> p_dict, KmerKey p_key, int p_id, int p_position)
        {
            if (p_dict.TryGetValue(p_key, out var l_existing))
            {
                if (l_existing.Item1 is 曖昧kmer番号 || l_existing.Item1 == p_id)
                {
                    return;
                }
                p_dict[p_key] = (曖昧kmer番号, 0);
                return;
            }
            p_dict[p_key] = (p_id, p_position);
        }

        /// <summary>
        /// 一方の unitig の末尾がもう一方の unitig の先頭へ伸びると辺が張られる
        /// </summary>
        [Fact]
        public void 一方のunitigの末尾がもう一方の先頭へ伸びると辺が張られる()
        {
            // unitig A の末尾 k-1 塩基が unitig B の先頭 k-1 塩基と一致する構成
            // A の末尾 k-mer から 1 塩基伸ばすと、ちょうど B の先頭 k-mer になる
            const int k = 8;
            const string l_shared = "CGTTACA"; // k-1 = 7 塩基の重なり
            var l_a = "GCTAAAGACAATTAC" + l_shared;      // 末尾が shared
            var l_b = l_shared + "GGATCCTTAGGCAAT";      // 先頭が shared

            var (l_unitigList, l_kmerDict) = V_構築(k, l_a, l_b);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            var l_bForward = ContigMaker.Get_頂点番号(2);

            Assert.Contains(l_bForward, l_graph.A_出辺[l_aForward]);

            // 逆鎖対称性: A→B があるなら B' →A' も存在しなければならない
            // これが崩れると順鎖側と逆鎖側で別々の経路が組まれ、同じ領域が
            // 2 通りに組み立てられてしまう
            Assert.Contains(l_aForward ^ 1, l_graph.A_出辺[l_bForward ^ 1]);

            // 入次数は双子の出次数で表せる
            Assert.Equal(l_graph.A_出辺[l_bForward ^ 1].Count, l_graph.Get_入次数(l_bForward));
        }

        /// <summary>
        /// k-1 で重ならない unitig 同士には辺が張られない
        /// </summary>
        [Fact]
        public void k_1で重ならないunitig同士には辺が張られない()
        {
            // 互いに無関係な 2 本
            // 偶然の短い一致があっても辺は張られてはならない
            // (旧実装はここで平均 3 塩基程度の一致による誤結合を作っていた)
            const int k = 8;
            const string l_a = "GCTAAAGACAATTACATAA";
            const string l_b = "TTGACCTGAATCCGGTTCA";

            var (l_unitigList, l_kmerDict) = V_構築(k, l_a, l_b);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            for (var v = 2; v < l_graph.A_出辺.Count; v++)
            {
                Assert.Empty(l_graph.A_出辺[v]);
            }
        }

        /// <summary>
        /// 一致する k-mer が相手の先頭ではなく途中にある場合は、辺が張られない
        /// </summary>
        [Fact]
        public void 一致するkmerが相手の先頭でなく途中にある場合は辺が張られない()
        {
            // B の「途中」に一致する k-mer があっても、k-1 オーバーラップでの
            // 連結はできないため辺を張ってはならない (Position != 0 を弾く条件)
            // A の末尾 k-mer を 1 塩基伸ばした k-mer が B の 3 塩基目から始まるよう構成する
            const int k = 8;
            const string l_junction = "ACGGATCA"; // A の末尾から伸ばして得られる k-mer
            var l_a = "GCTAAAGACAATTAC" + l_junction[..(k - 1)]; // 末尾 k-1 が junction の先頭 k-1
            var l_b = "TT" + l_junction + "CCTTAGGCAAT";          // junction は B の位置 2 に現れる

            var (l_unitigList, l_kmerDict) = V_構築(k, l_a, l_b);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            var l_bForward = ContigMaker.Get_頂点番号(2);
            Assert.DoesNotContain(l_bForward, l_graph.A_出辺[l_aForward]);
        }

        /// <summary>
        /// 末尾が 2 本の異なる unitig へ伸びる場合は、両方の辺が記録される
        /// </summary>
        [Fact]
        public void 末尾が2本の異なるunitigへ伸びる場合は両方の辺が記録される()
        {
            // 分岐: A の末尾から B へも C へも伸びられる構成
            // 辺は両方張られ、どちらを選ぶかはリード支持に委ねられる
            const int k = 8;
            const string l_shared = "CGTTACA";
            var l_a = "GCTAAAGACAATTAC" + l_shared;
            var l_b = l_shared + "GGATCCTTAGGCAAT";
            var l_c = l_shared + "TGATCCTTAGGCAAT"; // 分岐点の 1 塩基だけ B と異なる

            var (l_unitigList, l_kmerDict) = V_構築(k, l_a, l_b, l_c);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_graph.A_出辺[l_aForward].Count);
            Assert.Contains(ContigMaker.Get_頂点番号(2), l_graph.A_出辺[l_aForward]);
            Assert.Contains(ContigMaker.Get_頂点番号(3), l_graph.A_出辺[l_aForward]);
        }

        /// <summary>
        /// 単純バブル (u から 2 本に分かれ、それぞれ 1 本の unitig を経て
        /// 同じ w へ再合流する) で、リード支持の高い枝だけが経路として
        /// 残ることを確認する
        /// </summary>
        /// <remarks>
        /// 結合の採用条件を相互一意にした結果、再合流点 w の入次数が
        /// 2 のままだと u から w へ至る経路が一切結合されなくなるため、
        /// この処理が無いとバブルのたびに contig が千切れる
        /// </remarks>
        [Fact]
        public void 単純バブルはリード支持の高い枝だけを残し逆鎖も対称に除去する()
        {
            const int k = 8;
            // k=8 で全 unitig を通じて重複する正規化 k-mer が無いことを確認済みの構成
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_b1 = "TTACGCAAGGATCCTGCACGT"; // u の末尾7塩基 + 'A' で始まり、w の先頭7塩基で終わる
            const string l_b2 = "TTACGCACTTAGCATGCACGT"; // 分岐点の1塩基だけ b1 と異なる同長の枝
            const string l_w = "TGCACGTAAGGCTTACCA";

            var (l_unitigList, l_kmerDict) = V_構築(k, l_u, l_b1, l_b2, l_w);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b1V = ContigMaker.Get_頂点番号(2);
            var l_b2V = ContigMaker.Get_頂点番号(3);
            var l_wV = ContigMaker.Get_頂点番号(4);

            // 前提: バブル構造が実際に構築されている
            Assert.Equal(2, l_graph.A_出辺[l_uV].Count);
            Assert.Equal(2, l_graph.Get_入次数(l_wV));

            // b1 側にだけリード支持を与える
            Dictionary<(int, int), ulong> l_support = new()
            {
                [(l_uV, l_b1V)] = 40,
                [(l_uV, l_b2V)] = 3,
            };

            var l_popped = l_graph.V_除去_単純バブル(l_unitigList, l_support, k);

            Assert.Equal(1, l_popped);
            Assert.Equal([l_b1V], l_graph.A_出辺[l_uV]);
            Assert.Equal([l_wV], l_graph.A_出辺[l_b1V]);
            Assert.Empty(l_graph.A_出辺[l_b2V]);

            // 逆鎖側も対称に取り除かれていること (片側だけ消すと順鎖と逆鎖で
            // 別々の経路が組まれてしまう)
            Assert.Equal(1, l_graph.Get_入次数(l_wV));
            Assert.DoesNotContain(l_b2V ^ 1, l_graph.A_出辺[l_wV ^ 1]);
            Assert.DoesNotContain(l_uV ^ 1, l_graph.A_出辺[l_b2V ^ 1]);
        }

        /// <summary>
        /// 長さが大きく異なる分岐は、同じ領域の別表現 (バブル) ではなく
        /// 本物の分岐 (反復配列の出入口など) である可能性が高いため、
        /// 支持の低い側であっても勝手に経路から外してはならない
        /// </summary>
        [Fact]
        public void 長さが大きく異なる分岐は単純バブルとして除去しない()
        {
            const int k = 8;
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_b1 = "TTACGCAAGGATCCTGCACGT";
            // b2 は b1 より大幅に長い (長さ比が既定の閾値 1.5 を超える)
            const string l_b2 = "TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT";
            const string l_w = "TGCACGTAAGGCTTACCA";

            var (l_unitigList, l_kmerDict) = V_構築(k, l_u, l_b1, l_b2, l_w);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_graph.A_出辺[l_uV].Count);

            Dictionary<(int, int), ulong> l_support = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40,
                [(l_uV, ContigMaker.Get_頂点番号(3))] = 3,
            };

            var l_popped = l_graph.V_除去_単純バブル(l_unitigList, l_support, k);

            Assert.Equal(0, l_popped);
            Assert.Equal(2, l_graph.A_出辺[l_uV].Count);
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_length">作る長さ</param>
        /// <param name="p_seed">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_length, int p_seed)
        {
            var l_rng = new Random(p_seed);
            return string.Concat(Enumerable.Range(0, p_length).Select(_ => "ACGT"[l_rng.Next(4)]));
        }

        /// <summary>
        /// バブルの枝が単一 unitig とは限らない
        /// </summary>
        /// <remarks>
        /// 分岐の無い (排他的な) 2 本の
        /// unitig をまたぐ枝同士でも、1 本の経路として検出・比較できること
        /// </remarks>
        [Fact]
        public void unitig2本からなる鎖を1本の枝として扱う()
        {
            const int k = 8;
            const string l_uTail = "TTACGCA"; // u の末尾7塩基(既存テストと共通)
            const string l_wHead = "TGCACGT"; // w の先頭7塩基(既存テストと共通)
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_w = "TGCACGTAAGGCTTACCA";

            // 枝 1・枝 2 は完全に無関係な乱数配列にする (配列類似度の検証 (第 2 段階) は
            // 別のテストで確かめるため、ここでは第 2 段階を無効にして多段構造の
            // 検出・長さ帯の比較そのものだけを確かめる)
            var l_joint1 = V_生成_乱数配列(7, p_seed: 501);
            var l_joint2 = V_生成_乱数配列(7, p_seed: 502);
            var l_middleA1 = V_生成_乱数配列(10, p_seed: 511);
            var l_middleA2 = V_生成_乱数配列(10, p_seed: 521);
            var l_middleB1 = V_生成_乱数配列(10, p_seed: 512);
            var l_middleB2 = V_生成_乱数配列(10, p_seed: 522);

            var l_b1a = l_uTail + l_middleA1 + l_joint1;
            var l_b1b = l_joint1 + l_middleB1 + l_wHead;
            var l_b2a = l_uTail + l_middleA2 + l_joint2;
            var l_b2b = l_joint2 + l_middleB2 + l_wHead;

            var (l_unitigList, l_kmerDict) = V_構築(k, l_u, l_b1a, l_b1b, l_b2a, l_b2b, l_w);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b1aV = ContigMaker.Get_頂点番号(2);
            var l_b1bV = ContigMaker.Get_頂点番号(3);
            var l_b2aV = ContigMaker.Get_頂点番号(4);
            var l_b2bV = ContigMaker.Get_頂点番号(5);
            var l_wV = ContigMaker.Get_頂点番号(6);

            // 前提: それぞれの枝が 2 unitig の分岐無しの鎖になっている
            Assert.Equal(2, l_graph.A_出辺[l_uV].Count);
            Assert.Equal([l_b1bV], l_graph.A_出辺[l_b1aV]);
            Assert.Equal([l_b2bV], l_graph.A_出辺[l_b2aV]);
            Assert.Equal(2, l_graph.Get_入次数(l_wV));

            // 枝 1 にだけリード支持を与える
            Dictionary<(int, int), ulong> l_support = new()
            {
                [(l_uV, l_b1aV)] = 40,
                [(l_uV, l_b2aV)] = 3,
            };

            // 配列類似度の検証 (第 2 段階) は別のテストで確かめる
            // ここでは多段構造の検出・長さ帯の比較だけを見たいので無効にする
            var l_popped = l_graph.V_除去_単純バブル(l_unitigList, l_support, k, p_類似度の下限: 0.0);

            Assert.Equal(1, l_popped);
            // 枝 1(2 unitig とも) は生き残り、枝 2(2 unitig とも) は取り除かれる
            Assert.Equal([l_b1aV], l_graph.A_出辺[l_uV]);
            Assert.Equal([l_wV], l_graph.A_出辺[l_b1bV]);
            Assert.Empty(l_graph.A_出辺[l_b2aV]);
            Assert.Empty(l_graph.A_出辺[l_b2bV]);
            Assert.Equal(1, l_graph.Get_入次数(l_wV));
        }

        /// <summary>
        /// careful_bubble: 除去された側の経路の配列を、引き継ぎ先へ集められること
        /// </summary>
        /// <remarks>
        /// 「この k では敗者と判断したが、次の k は自分の証拠で判断し直せる」ため、
        /// 配列自体は捨てない
        /// </remarks>
        [Fact]
        public void 引き継ぎ先を指定すると除去された側の配列を集められる()
        {
            const int k = 8;
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_b1 = "TTACGCAAGGATCCTGCACGT";
            const string l_b2 = "TTACGCACTTAGCATGCACGT";
            const string l_w = "TGCACGTAAGGCTTACCA";

            var (l_unitigList, l_kmerDict) = V_構築(k, l_u, l_b1, l_b2, l_w);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            var l_b2V = ContigMaker.Get_頂点番号(3);

            Dictionary<(int, int), ulong> l_support = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40,
                [(l_uV, l_b2V)] = 3,
            };

            List<string> l_carryOver = [];
            var l_popped = l_graph.V_除去_単純バブル(l_unitigList, l_support, k, l_carryOver);

            Assert.Equal(1, l_popped);
            Assert.Equal([l_b2], l_carryOver);
        }

        /// <summary>
        /// 長さが揃っていても配列がまるで違う (たまたま長さが一致しただけの
        /// 別の反復など) 場合は、同じ領域の別表現とは言えないため触らない
        /// </summary>
        [Fact]
        public void 長さは同じでも配列が無関係な分岐は除去しない()
        {
            const int k = 8;
            const string l_uTail = "TTACGCA";
            const string l_wHead = "TGCACGT";
            const string l_u = "GCTAAAGACAATTACGCA";
            const string l_w = "TGCACGTAAGGCTTACCA";

            // 長さは完全に一致するが、中身は無関係な乱数配列
            var l_b1 = l_uTail + V_生成_乱数配列(20, p_seed: 601) + l_wHead;
            var l_b2 = l_uTail + V_生成_乱数配列(20, p_seed: 602) + l_wHead;

            var (l_unitigList, l_kmerDict) = V_構築(k, l_u, l_b1, l_b2, l_w);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, l_graph.A_出辺[l_uV].Count);

            Dictionary<(int, int), ulong> l_support = new()
            {
                [(l_uV, ContigMaker.Get_頂点番号(2))] = 40,
                [(l_uV, ContigMaker.Get_頂点番号(3))] = 3,
            };

            var l_popped = l_graph.V_除去_単純バブル(l_unitigList, l_support, k);

            Assert.Equal(0, l_popped);
            Assert.Equal(2, l_graph.A_出辺[l_uV].Count);
        }

        /// <summary>
        /// 自己ループの辺は作らない
        /// </summary>
        [Fact]
        public void 自己ループの辺は作らない()
        {
            // 自己ループを辺として持つと walk が同じ unitig を無限に伸ばしうるため、
            // 構築段階で除外していることを確認する
            const int k = 8;
            const string l_repeatUnit = "ACGGATCT";
            var l_a = l_repeatUnit + "GCTAAAGA" + l_repeatUnit[..(k - 1)];

            var (l_unitigList, l_kmerDict) = V_構築(k, l_a);
            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k, 曖昧kmer番号);

            var l_aForward = ContigMaker.Get_頂点番号(1);
            Assert.DoesNotContain(l_aForward, l_graph.A_出辺[l_aForward]);
        }
    }
}
