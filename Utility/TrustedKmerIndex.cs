using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-merの出現回数カウント(CountingDBによる厳密な外部マージソート)と、
    /// カットオフを通過した「信頼できるk-mer」の厳密な集合を保持するクラス。
    ///
    /// 以前はカットオフ後の所属判定に多重ハッシュのビット配列
    /// (Bloom filter)を使っていたが、複数ハッシュのうち1つが
    /// 事実上塩基の並び順に依存しない合計値に退化しておりハッシュの
    /// 独立性が低いこと、そもそも近似判定である以上フォールスポジティブに
    /// よるグラフ構造の誤判定(誤った分岐点検出・誤った隣接判定)を
    /// 原理的に排除できないことが問題だった。
    ///
    /// 7Mbp程度のバクテリアゲノムであれば、カットオフを通過した信頼できる
    /// k-merの総数は現実的にせいぜい数千万件程度に収まり、厳密な集合として
    /// メモリに保持できる規模である。そのためカットオフ後は厳密な集合に
    /// 置き換え、近似判定を完全に排除した。
    /// </summary>
    internal class TrustedKmerIndex : IDisposable
    {
        private readonly string _一時ディレクトリ;

        // k-mer カウント用のシャード。k-mer 自身のハッシュで振り分けるため、
        // ある k-mer は必ず1つのシャードにしか載らない。
        private CountingDB[]? _カウンタ群;

        // シャードごとのロック。k-mer をワーカー単位ではなくハッシュ値で
        // 振り分けるようにしたため、複数スレッドが同じシャードへ書きうる。
        private object[]? _シャードロック;

        // カットオフ実行後に確定する、カットオフを通過したk-merの厳密な集合
        // (常に正規化された形で保持する)。値はそのk-merの出現回数
        // (カバレッジ)。GraphSimplifier の低カバレッジ端の除去に使う。
        // 以降のグラフ探索(開始点判定・次数計算・UnitigMakerの伸長判定)は、
        // 値を見ない所属判定のみで行う。
        //
        // k<=32 の場合は 2bit パックした ulong 1個で直接引く高速経路、
        // 33<=k<=64 の場合は UInt128 の高速経路を使う。どちらも値型なので
        // KmerKey(ulong[] を毎回ヒープ確保する)と違い割り当てが発生しない。
        // ErrorCorrector は1リードあたり数百〜数千回この判定を呼ぶため、
        // ヒープ確保のオーバーヘッドが実データ規模で致命的に効く。
        // k>64 の場合のみ、従来通り厳密だが低速な KmerKey 経路にフォールバックする。
        private Dictionary<KmerKey, ulong>? _信頼kmer_大;

        private Dictionary<ulong, ulong>? _信頼kmer_小;

        // 33 <= k <= 64 用。150bp リードで k=31 のままだと 31bp 以上の反復配列が
        // すべて潰れてしまい contig N50 が伸びないため、k を 63 前後まで上げられる
        // ことが品質上きわめて重要になる。
        private Dictionary<UInt128, ulong>? _信頼kmer_中;

        private static bool 小経路を使うか => ConfigurationManager.A_実行時引数.A_k長 <= 32;

        private static bool 中経路を使うか => ConfigurationManager.A_実行時引数.A_k長 is > 32 and <= 64;

        public TrustedKmerIndex(string p_一時ディレクトリ)
        {
            this._一時ディレクトリ = p_一時ディレクトリ;
            var l_シャード数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            this._カウンタ群 = new CountingDB[l_シャード数];
            this._シャードロック = new object[l_シャード数];
            for (var i = 0; i < l_シャード数; i++)
            {
                this._カウンタ群[i] = new CountingDB(p_一時ディレクトリ, l_シャード数);
                this._シャードロック[i] = new object();
            }
        }

        /// <summary>
        /// 曖昧塩基(N など、候補が複数ある位置)を含む k-mer を、
        /// ありうる塩基の組み合わせすべてに展開して登録する。
        /// p_塩基候補列[i] はその位置で取りうる塩基ID の一覧。
        ///
        /// 展開は塩基ID の空間で行い、1件ずつ通常の登録へ渡す。
        /// こうすることで正規化(順鎖・逆鎖の寄せ)とシャード振り分けが
        /// 通常経路とまったく同じ扱いになる。
        ///
        /// 以前は CountingDB 側でパック済みバイト列を組み立てており、
        /// (a) 正規化されない (b) 組み立て途中のバッファをそのまま辞書のキーと
        /// して格納しており、その後も書き換え続けるため既に登録済みのキーが
        /// 壊れる、という2つの問題があった。
        /// </summary>
        public void V_登録_曖昧塩基あり(Span<byte[]> p_塩基候補列, int p_ワーカー番号)
        {
            if (this._カウンタ群 is null)
            {
                return;
            }
            var l_kmer = new byte[p_塩基候補列.Length];
            this.V_登録_組み合わせ展開(p_塩基候補列, 0, l_kmer, p_ワーカー番号);
        }

        private void V_登録_組み合わせ展開(Span<byte[]> p_塩基候補列, int p_位置, byte[] p_kmer, int p_ワーカー番号)
        {
            if (p_位置 == p_塩基候補列.Length)
            {
                this.V_登録(p_kmer.AsSpan(), p_ワーカー番号);
                return;
            }
            foreach (var l_塩基ID in p_塩基候補列[p_位置])
            {
                p_kmer[p_位置] = l_塩基ID;
                this.V_登録_組み合わせ展開(p_塩基候補列, p_位置 + 1, p_kmer, p_ワーカー番号);
            }
        }

        /// <summary>
        /// k-mer を1件カウントする。
        ///
        /// 振り分けはワーカー番号ではなく k-mer 自身のハッシュ値で行う。
        /// ワーカーごとに別の辞書へ入れていた頃は、同じ k-mer が最大で
        /// スレッド数ぶんの辞書に重複して載り、メモリも書き出し量も
        /// そのぶん膨らんでいた(実データ 100x・16スレッドでピーク12.5GB)。
        /// ハッシュで振り分ければ、ある k-mer は必ず1つのシャードにしか載らない。
        ///
        /// さらに、順鎖・逆鎖のうち辞書順で小さいほう(正規化形)に寄せてから
        /// 数える。以前は両向きを別キーとして数え、カットオフ時に合算していたため、
        /// エントリ数・書き出し量ともに2倍になっていた。
        /// </summary>
        public void V_登録(Span<byte> p_kmer, int p_ワーカー番号)
        {
            if (this._カウンタ群 is not { } l_カウンタ群)
            {
                return;
            }

            var l_パック済み = Get_正規化パック(p_kmer);
            var l_シャード = (int)(Get_ハッシュ(l_パック済み) % (uint)l_カウンタ群.Length);
            lock (this._シャードロック![l_シャード])
            {
                l_カウンタ群[l_シャード].V_登録_パック済み(l_パック済み);
            }
        }

        /// <summary>
        /// k-mer を、順鎖・逆鎖のうち塩基列として辞書順で小さいほうの向きで
        /// 2bit パックした byte 配列にする。
        ///
        /// パック後のバイト列の辞書順は塩基列の辞書順と一致する(先頭塩基が
        /// 上位ビットに来るため)ので、外部マージソートの順序とも整合する。
        /// </summary>
        private static byte[] Get_正規化パック(ReadOnlySpan<byte> p_kmer)
        {
            var l_順鎖を使うか = Get_順鎖が正規形か(p_kmer);
            var l_パック済み = new byte[(p_kmer.Length + 3) / 4];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                // 逆鎖側を採用する場合は、末尾から相補塩基を取り出す。
                // 相補は A(1)<->T(4), C(2)<->G(3) なので 5 - x で得られる。
                var l_塩基ID = l_順鎖を使うか ? p_kmer[i] : (byte)(5 - p_kmer[p_kmer.Length - 1 - i]);
                l_パック済み[i >> 2] |= (byte)((l_塩基ID - 1) << ((3 - (i & 3)) << 1));
            }
            return l_パック済み;
        }

        /// <summary>
        /// 順鎖側がその逆相補以下(辞書順)かどうか。確保なしで判定する。
        /// </summary>
        private static bool Get_順鎖が正規形か(ReadOnlySpan<byte> p_kmer)
        {
            int i = 0, j = p_kmer.Length - 1;
            while (i <= j)
            {
                var l_順鎖 = p_kmer[i];
                var l_逆鎖 = (byte)(5 - p_kmer[j]);
                if (l_順鎖 != l_逆鎖)
                {
                    return l_順鎖 < l_逆鎖;
                }
                i++;
                j--;
            }
            // 回文(自身が逆相補と一致)。どちらでも同じなので順鎖扱い。
            return true;
        }

        /// <summary>パック済みキーの FNV-1a ハッシュ。シャードの振り分けに使う。</summary>
        private static uint Get_ハッシュ(byte[] p_パック済みkmer)
        {
            var l_ハッシュ = 2166136261u;
            foreach (var l_バイト in p_パック済みkmer)
            {
                l_ハッシュ ^= l_バイト;
                l_ハッシュ *= 16777619u;
            }
            return l_ハッシュ;
        }

        /// <summary>
        /// kmer(順鎖・逆鎖いずれの向きでもよい)がカットオフを通過した
        /// 信頼できるk-mer集合に含まれるかどうかを厳密に判定する。
        /// </summary>
        public bool Get_含まれるか(Span<byte> p_kmer)
        {
            if (小経路を使うか)
            {
                return this._信頼kmer_小!.ContainsKey(Get_正規形_小(p_kmer));
            }
            if (中経路を使うか)
            {
                return this._信頼kmer_中!.ContainsKey(Get_正規形_中(p_kmer));
            }
            return this._信頼kmer_大!.ContainsKey(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// kmerの出現回数(カバレッジ)を返す。信頼できるk-mer集合に
        /// 含まれない場合は0を返す。
        /// </summary>
        public ulong Get_カバレッジ(Span<byte> p_kmer)
        {
            if (小経路を使うか)
            {
                return this._信頼kmer_小!.GetValueOrDefault(Get_正規形_小(p_kmer), 0UL);
            }
            if (中経路を使うか)
            {
                return this._信頼kmer_中!.GetValueOrDefault(Get_正規形_中(p_kmer), 0UL);
            }
            return this._信頼kmer_大!.GetValueOrDefault(new KmerKey(p_kmer).Get_正規形(), 0UL);
        }

        /// <summary>
        /// kmer(塩基ID 1-4、長さ32以下)を2bit/塩基でulong1個にパックする。
        /// 先頭塩基が最上位側、末尾塩基が最下位側に来る(空きビットは下位側に残る)。
        /// </summary>
        private static ulong Get_パック_小(ReadOnlySpan<byte> p_kmer)
        {
            var l_値 = 0UL;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | ((ulong)l_塩基ID - 1);
            }
            return l_値;
        }

        /// <summary>
        /// Get_パック_小 でパックした値の逆相補を、ヒープ確保なしで直接計算する。
        /// 2bitコドンごとに相補を取り(A&lt;-&gt;T, C&lt;-&gt;G)、
        /// 下位から順に取り出しつつ上位へ積み直すことでコドン順序も反転させる。
        /// </summary>
        private static ulong Get_逆相補_小(ulong p_パック済み, int p_長さ)
        {
            var l_残り = p_パック済み;
            var l_結果 = 0UL;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_コドン = l_残り & 0x3UL;
                l_結果 = (l_結果 << 2) | (l_コドン ^ 0x3UL);
                l_残り >>= 2;
            }
            return l_結果;
        }

        private static ulong Get_正規形_小(ReadOnlySpan<byte> p_kmer)
        {
            var l_パック済み = Get_パック_小(p_kmer);
            var l_逆相補 = Get_逆相補_小(l_パック済み, p_kmer.Length);
            return Math.Min(l_パック済み, l_逆相補);
        }

        /// <summary>Get_パック_小 の逆変換。末尾塩基が最下位ビット側にあるため、末尾から復元する。</summary>
        private static byte[] Get_復元_小(ulong p_パック済み, int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((p_パック済み & 0x3UL) + 1);
                p_パック済み >>= 2;
            }
            return l_塩基列;
        }

        /// <summary>
        /// Get_パック_小 の 128bit 版(k は 64 以下)。ビット配置の規約は同じで、
        /// kmer の先頭塩基が最上位側、末尾塩基が最下位側に来る。
        /// </summary>
        private static UInt128 Get_パック_中(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_値 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | (UInt128)(l_塩基ID - 1);
            }
            return l_値;
        }

        /// <summary>Get_逆相補_小 の 128bit 版。</summary>
        private static UInt128 Get_逆相補_中(UInt128 p_パック済み, int p_長さ)
        {
            var l_残り = p_パック済み;
            UInt128 l_結果 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_コドン = l_残り & 3;
                l_結果 = (l_結果 << 2) | (l_コドン ^ 3);
                l_残り >>= 2;
            }
            return l_結果;
        }

        private static UInt128 Get_正規形_中(ReadOnlySpan<byte> p_kmer)
        {
            var l_パック済み = Get_パック_中(p_kmer);
            var l_逆相補 = Get_逆相補_中(l_パック済み, p_kmer.Length);
            return l_パック済み < l_逆相補 ? l_パック済み : l_逆相補;
        }

        /// <summary>Get_パック_中 の逆変換。</summary>
        private static byte[] Get_復元_中(UInt128 p_パック済み, int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((ulong)(p_パック済み & 3) + 1);
                p_パック済み >>= 2;
            }
            return l_塩基列;
        }

        /// <summary>
        /// カットオフを通過した信頼できるk-merを(正規化された、いずれかの向きの)
        /// byte配列として1件ずつ列挙する。GraphSimplifier の tip 除去等、
        /// 集合全体を舐めて再判定する処理から使う。
        /// </summary>
        public IEnumerable<byte[]> Get_信頼kmer一覧()
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            if (小経路を使うか)
            {
                foreach (var l_パック済み in this._信頼kmer_小!.Keys)
                {
                    yield return Get_復元_小(l_パック済み, l_k長);
                }
            }
            else if (中経路を使うか)
            {
                foreach (var l_パック済み in this._信頼kmer_中!.Keys)
                {
                    yield return Get_復元_中(l_パック済み, l_k長);
                }
            }
            else
            {
                foreach (var l_キー in this._信頼kmer_大!.Keys)
                {
                    yield return l_キー.Get_塩基列(l_k長);
                }
            }
        }

        /// <summary>
        /// kmerを信頼できるk-mer集合から除去する(順鎖・逆鎖どちらの向きで
        /// 渡してもよい)。GraphSimplifier が tip の構成k-merを取り除く際に使う。
        /// </summary>
        public void V_除去(ReadOnlySpan<byte> p_kmer)
        {
            if (小経路を使うか)
            {
                _ = this._信頼kmer_小!.Remove(Get_正規形_小(p_kmer));
            }
            else if (中経路を使うか)
            {
                _ = this._信頼kmer_中!.Remove(Get_正規形_中(p_kmer));
            }
            else
            {
                _ = this._信頼kmer_大!.Remove(new KmerKey(p_kmer).Get_正規形());
            }
        }

        /// <summary>
        /// 現在の信頼できるk-mer集合を1回走査し、unitigの開始点となる
        /// k-merをすべて再検出する。ファイルの読み直しではなく
        /// インメモリの集合をそのまま使うため、tip 除去で集合を縮小した後の
        /// 再構築にも安価に使える。
        ///
        /// Get_信頼kmer一覧 は各座位ごとに正規化された向きのk-merを1つだけ返すが、
        /// 開始点判定は向き依存(そのk-mer自身の入次数を見る)である。ある座位が
        /// 「順鎖では分岐点の直後」でも「逆鎖(=正規形側)では分岐点そのものでは
        /// ない」ことがありえるため、正規形側だけを調べると本来開始点であるべき
        /// 向きを見逃す(小さなバブルのテストケースで実際に見逃しを確認:
        /// 分岐解消後も低カバレッジ側の枝が除去されないままになっていた)。
        /// そのため、各座位について両方の向きを個別に判定する。
        /// </summary>
        public List<byte[]> Get_開始kmer一覧()
        {
            List<byte[]> l_開始kmer = [];
            foreach (var l_kmer in this.Get_信頼kmer一覧())
            {
                if (this.Get_開始kmerか(l_kmer))
                {
                    l_開始kmer.Add(l_kmer);
                }

                var l_逆相補 = Util.V_逆相補(l_kmer).ToArray();
                if (this.Get_開始kmerか(l_逆相補))
                {
                    l_開始kmer.Add(l_逆相補);
                }
            }
            return l_開始kmer;
        }

        public List<byte[]> V_カットオフ(ulong p_カットオフ)
        {
            // 各シャードの CountingDB をそれぞれ統合し、
            // 出来上がった複数のソート済みファイルをさらに1本にマージする。
            var l_統合済みファイル = new List<string>();
            foreach (var l_カウンタ in this._カウンタ群!)
            {
                l_統合済みファイル.Add(l_カウンタ.Get_統合ファイル());
                l_カウンタ.Dispose();
            }
            this._カウンタ群 = null;

            var l_ファイルパス = CountingDB.Get_統合ファイル_シャード間(this._一時ディレクトリ, l_統合済みファイル);

            var l_パック長 = (ConfigurationManager.A_実行時引数.A_k長 + 3) / 4;
            var l_小経路 = 小経路を使うか;
            var l_中経路 = 中経路を使うか;
            var l_信頼kmer_大 = l_小経路 || l_中経路 ? null : new Dictionary<KmerKey, ulong>();
            var l_信頼kmer_小 = l_小経路 ? new Dictionary<ulong, ulong>() : null;
            var l_信頼kmer_中 = l_中経路 ? new Dictionary<UInt128, ulong>() : null;
            using (var l_読み込み = new BinaryReader(File.Open(l_ファイルパス, FileMode.Open, FileAccess.Read)))
            {
                ulong l_採用数 = 0;
                ulong l_総種類数 = 0;
                // 出現回数 -> その回数を持つユニークk-merの種類数。
                // エラー由来の低頻度k-merと真のゲノム由来k-merを分ける「谷」を
                // 推定するために、カットオフ判定と同じこのループで集計する
                // (このファイルはこの後削除されるため、ここでしか見られない)。
                Dictionary<ulong, long> l_ヒストグラム = [];
                while (Util.Get_続きがあるか(l_読み込み))
                {
                    var l_パック済み = l_読み込み.ReadBytes(l_パック長);
                    var l_出現回数 = l_読み込み.ReadUInt64();
                    l_総種類数 += 1;
                    l_ヒストグラム[l_出現回数] = l_ヒストグラム.GetValueOrDefault(l_出現回数, 0L) + 1;
                    if (l_出現回数 >= p_カットオフ)
                    {
                        l_採用数 += 1;
                        List<byte> l_塩基列 = [];
                        foreach (var l_バイト in l_パック済み)
                        {
                            l_塩基列.AddRange(Util.V_変換_塩基列(l_バイト));
                        }
                        var l_kmer = CollectionsMarshal.AsSpan(l_塩基列)[..ConfigurationManager.A_実行時引数.A_k長];
                        // カウント段階で既に正規形へ寄せてあるため、同じ正規形が
                        // 複数エントリとして現れることはない。それでも加算で受けて
                        // おけば、将来カウント側の正規化をやめた場合でも壊れない。
                        if (l_小経路)
                        {
                            var l_正規形 = Get_正規形_小(l_kmer);
                            l_信頼kmer_小![l_正規形] = l_信頼kmer_小.GetValueOrDefault(l_正規形, 0UL) + l_出現回数;
                        }
                        else if (l_中経路)
                        {
                            var l_正規形 = Get_正規形_中(l_kmer);
                            l_信頼kmer_中![l_正規形] = l_信頼kmer_中.GetValueOrDefault(l_正規形, 0UL) + l_出現回数;
                        }
                        else
                        {
                            var l_正規形 = new KmerKey(l_kmer).Get_正規形();
                            l_信頼kmer_大![l_正規形] = l_信頼kmer_大.GetValueOrDefault(l_正規形, 0UL) + l_出現回数;
                        }
                    }
                }
                Console.WriteLine("kmer count: " + l_総種類数);
                Console.WriteLine("good kmer: " + l_採用数);
                Console.WriteLine($"[Info] k-mer count histogram (count:#distinct kmers): {KmerHistogram.Get_要約(l_ヒストグラム)}");
                var l_推奨カットオフ = KmerHistogram.Get_推奨カットオフ(l_ヒストグラム);
                if (l_推奨カットオフ is { } l_推奨値)
                {
                    var l_注記 = l_推奨値 == p_カットオフ ? " (matches the cutoff currently in effect)" : $" (currently using -kc {p_カットオフ})";
                    Console.WriteLine($"[Info] Suggested k-mer cutoff from histogram valley: {l_推奨値}{l_注記}");
                }
                else
                {
                    Console.WriteLine("[Info] Could not identify a clear histogram valley to suggest a k-mer cutoff (spectrum may not be bimodal at this coverage).");
                }
            }
            File.Delete(l_ファイルパス);
            this._信頼kmer_大 = l_信頼kmer_大;
            this._信頼kmer_小 = l_信頼kmer_小;
            this._信頼kmer_中 = l_信頼kmer_中;

            Console.WriteLine("Search First k-mer");
            // 以前はここで一度カットオフ通過k-merをファイルへ書き出し、
            // 読み直して開始点を判定していた。厳密な集合をインメモリで
            // 保持するようになったため、その集合を直接走査すれば同じ結果が
            // 得られ、ディスクI/Oを1往復省略できる。
            return this.Get_開始kmer一覧();
        }

        /// <summary>
        /// 与えられた k-mer が unitig の開始点かどうかを判定する。
        ///
        /// 入次数が0または2個以上(=分岐点そのもの)であれば当然開始点になる。
        /// 入次数がちょうど1の場合でも、その唯一の予測元自身が分岐点
        /// (出次数が1でない)であれば、UnitigMaker の前進walkは予測元の時点で
        /// 停止してしまいこのk-merへは到達しない(=このk-merは誰からも
        /// 「walkで訪れてもらえない」)ため、新たなunitigの開始点として
        /// 別途扱う必要がある。
        /// これを見落とすと、分岐点の直後から始まる配列がunitig化されず
        /// 丸ごと欠落する(小さなテストケースで実際に発生を確認した)。
        /// </summary>
        public bool Get_開始kmerか(Span<byte> p_kmer)
        {
            var l_入次数 = this.Get_入次数(p_kmer, out var l_唯一の予測元);
            if (l_入次数 != 1)
            {
                return true;
            }
            return this.Get_出次数(l_唯一の予測元!) != 1;
        }

        /// <summary>
        /// kmer への入次数(前方に接続しうる異なる1塩基拡張の数)を数える。
        ///
        /// UnitigMaker の前進伸長規則は「kmerの先頭1文字を落とし、
        /// 末尾に候補塩基cを付加する」(後続 = kmer[1..] + c)。
        /// この関係の逆(予測元)を解くと、予測元 P は
        /// 「P[1..] + (kmerの末尾文字) == kmer」を満たす必要があり、
        /// P[1..] = kmer[..^1](kmerの末尾を落としたもの)、
        /// P[0] = 任意の候補塩基c、すなわち P = c + kmer[..^1] となる。
        ///
        /// 以前の実装は候補 = c + kmer[1..](kmerの"先頭"を落としたもの)
        /// を試しており、c = kmer[0] のとき候補が kmer 自身と一致して
        /// しまう(=常に最低1回は自己ヒットする)退化バグがあった。これにより
        /// 真の入次数0(=配列の先頭)が絶対に検出できず、開始点判定が
        /// 意図通りに機能していなかった。
        /// </summary>
        public int Get_入次数(Span<byte> p_kmer)
        {
            return this.Get_入次数(p_kmer, out _);
        }

        /// <summary>
        /// 入次数計算の本体。入次数がちょうど1だった場合、その唯一の
        /// 予測元(kmer長のbyte配列)も同時に返す(開始点判定が使う)。
        /// </summary>
        private int Get_入次数(Span<byte> p_kmer, out byte[]? p_唯一の予測元)
        {
            var l_候補 = new byte[p_kmer.Length];
            p_kmer[..^1].CopyTo(l_候補.AsSpan(1));
            var l_件数 = 0;
            byte[]? l_一致 = null;
            for (byte i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                l_候補[0] = i;
                if (this.Get_含まれるか(l_候補))
                {
                    l_件数++;
                    l_一致 = l_件数 == 1 ? (byte[])l_候補.Clone() : null;
                }
            }
            p_唯一の予測元 = l_件数 == 1 ? l_一致 : null;
            return l_件数;
        }

        /// <summary>
        /// kmer からの出次数(後方に接続しうる異なる1塩基拡張の数)を数える。
        /// UnitigMaker の前進伸長規則(kmer[1..] + c)そのものを試す。
        /// GraphSimplifier がunitigの末尾端の次数(=tip判定)を見る際に使う。
        /// </summary>
        public int Get_出次数(Span<byte> p_kmer)
        {
            var l_候補 = new byte[p_kmer.Length];
            p_kmer[1..].CopyTo(l_候補.AsSpan(0, p_kmer.Length - 1));
            var l_件数 = 0;
            for (byte i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                l_候補[^1] = i;
                if (this.Get_含まれるか(l_候補))
                {
                    l_件数++;
                }
            }
            return l_件数;
        }

        public void Dispose()
        {
            if (this._カウンタ群 != null)
            {
                foreach (var l_カウンタ in this._カウンタ群)
                {
                    l_カウンタ.Dispose();
                }
            }
        }
    }
}
