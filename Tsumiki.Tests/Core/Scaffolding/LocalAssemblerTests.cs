using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Evidence;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// GapFiller が埋められなかった scaffold のギャップを、その両端に実際にマップされた局所リードだけで再アセンブリして埋める処理の検証
    /// </summary>
    /// <remarks>
    /// AssemblyMerger (-mg) と違い、他の k の「既に確定した結論」を持ち込むのではなく、生リードから新しく証拠を集める<br/>
    /// ここでは生の FASTQ を直接与えて、局所アセンブリだけでギャップが埋まる/埋まらないことを検証する (GapFiller 側は経由しない)
    /// </remarks>
    public class LocalAssemblerTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public LocalAssemblerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_localasm_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            // LocalAssembler は内部で TrustedKmerIndex/KmerKey を使うため、
            // 現在の実行時引数の k 長を、テストで使う k (21) に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 1 };
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 短い文脈で曖昧な反復を長いリード文脈で識別する
        /// </summary>
        [Fact]
        public void V_局所文脈_短い反復を識別()
        {
            var l_左 = V_生成_ランダム配列(300, p_シード: 91);
            var l_右 = V_生成_ランダム配列(300, p_シード: 92);
            var l_反復 = V_生成_ランダム配列(25, p_シード: 93);
            var l_間 = V_生成_ランダム配列(5, p_シード: 94);
            var l_充填 = l_反復 + l_間 + l_反復;
            var l_正解 = l_左 + l_充填 + l_右;
            var l_パス = this.V_書き出し_スキャフォールド("adaptive.fasta", l_左 + new string('N', l_充填.Length) + l_右);
            var l_リード = this.V_書き出し_リード("adaptive.fq", [l_正解], p_リード長: 150);
            var l_固定結果 = LocalAssembler.Get_固定kの局所結果(new 局所ギャップ(0, 300, l_充填.Length, l_左, l_右), Get_証拠(FastqReader.Get_生リード列(l_リード, null)), 21, out var l_固定判定);
            Assert.Null(l_固定結果);
            Assert.Equal(ギャップ充填判定.一意でない, l_固定判定);
            var l_結果 = LocalAssembler.V_充填_ギャップ(l_パス, l_リード, string.Empty, 21);
            Assert.Equal(1, l_結果.A_埋めたギャップ数);
            Assert.Equal(l_正解, Get_単一配列(l_パス));
        }

        /// <summary>高 k ではリード間が切れる低深度橋を低 k で回復する</summary>
        [Fact]
        public void V_非連結なら低kで真の橋を回復する()
        {
            const int l_基準k = 37;
            var l_左 = V_生成_ランダム配列(300, 501);
            var l_充填 = V_生成_ランダム配列(90, 502);
            var l_右 = V_生成_ランダム配列(300, 503);
            var l_正解 = l_左 + l_充填 + l_右;
            var l_局所リード = new List<string>();
            for (var i = l_左.Length - 35; i <= l_左.Length + l_充填.Length - 20; i += 5)
            {
                l_局所リード.Add(l_正解.Substring(i, 50));
            }

            var l_ギャップ = new 局所ギャップ(0, l_左.Length, l_充填.Length, l_左, l_右);
            var l_証拠 = Get_証拠(l_局所リード);
            Assert.Null(LocalAssembler.Get_固定kの局所結果(l_ギャップ, l_証拠, l_基準k, out var l_固定判定));
            Assert.Equal(ギャップ充填判定.到達不能, l_固定判定);

            var l_結果 = LocalAssembler.Get_適応kの局所結果(l_ギャップ, l_証拠, l_基準k, out var l_判定);
            Assert.Equal(ギャップ充填判定.充填済み, l_判定);
            Assert.Equal(l_充填, l_結果);
        }

        /// <summary>
        /// base-k の結果も、fallback 候補と同じ「複数リードの支持」ゲートを通すべきことを確かめる
        /// </summary>
        /// <remarks>
        /// 現状は Get_適応kの局所結果 が base-k の非 null 結果を Has経路支持 を経ずに即返すため、
        /// 橋渡し全体をまたぐリードが 1 本しかなくても採用してしまう (P1 で修正する非対称性)
        /// </remarks>
        [Fact]
        public void Get_適応kの局所結果_base_kでも複数リードの支持が無ければ採用しない()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, 701);
            var l_橋 = V_生成_ランダム配列(40, 702);
            var l_右 = V_生成_ランダム配列(300, 703);
            var l_正解 = l_左 + l_橋 + l_右;
            // 橋渡し全体をまたぐリードが 1 本だけ (fallback 候補に課す「2 本以上」の基準を満たさない)
            var l_単一リード = new List<string> { l_正解.Substring(l_左.Length - 30, 30 + l_橋.Length + 30) };

            var l_ギャップ = new 局所ギャップ(0, l_左.Length, l_橋.Length, l_左, l_右);
            var l_結果 = LocalAssembler.Get_適応kの局所結果(l_ギャップ, Get_証拠(l_単一リード), l_k長, out _);

            Assert.Null(l_結果);
        }

        /// <summary>
        /// 逆相補で与えても同じ結果になることを確かめる
        /// </summary>
        /// <remarks>
        /// k-mer の正規化は順鎖・逆鎖どちらでも同じキーに寄せる設計のはずで、
        /// 入力の鎖の向きだけで結果が変わってはいけない
        /// </remarks>
        [Fact]
        public void Get_適応kの局所結果_逆相補でも同じ結果になる()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, 901);
            var l_橋 = V_生成_ランダム配列(60, 902);
            var l_右 = V_生成_ランダム配列(300, 903);
            var l_正解 = l_左 + l_橋 + l_右;
            var l_リード = new List<string>();
            for (var i = 0; i + 100 <= l_正解.Length; i += 1)
            {
                l_リード.Add(l_正解.Substring(i, 100));
            }

            var l_ギャップ = new 局所ギャップ(0, l_左.Length, l_橋.Length, l_左, l_右);
            var l_順鎖結果 = LocalAssembler.Get_適応kの局所結果(l_ギャップ, Get_証拠(l_リード), l_k長, out var l_順鎖判定);

            var l_逆相補ギャップ = new 局所ギャップ(0, l_左.Length, l_橋.Length, Util.V_逆相補(l_右), Util.V_逆相補(l_左));
            var l_逆相補リード = l_リード.Select(Util.V_逆相補).ToList();
            var l_逆鎖結果 = LocalAssembler.Get_適応kの局所結果(l_逆相補ギャップ, Get_証拠(l_逆相補リード), l_k長, out var l_逆鎖判定);

            Assert.Equal(ギャップ充填判定.充填済み, l_順鎖判定);
            Assert.Equal(l_橋, l_順鎖結果);
            Assert.Equal(ギャップ充填判定.充填済み, l_逆鎖判定);
            Assert.Equal(Util.V_逆相補(l_橋), l_逆鎖結果);
        }

        /// <summary>
        /// 同じリードを重複させて渡しても、1本のときと同じ判定になることを確かめる
        /// </summary>
        /// <remarks>
        /// PCR 重複のような同一配列の繰り返しを、独立した追加の支持として数えてはいけない
        /// </remarks>
        [Fact]
        public void Get_適応kの局所結果_同一リードの重複は追加の支持にならない()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, 911);
            var l_橋 = V_生成_ランダム配列(40, 912);
            var l_右 = V_生成_ランダム配列(300, 913);
            var l_正解 = l_左 + l_橋 + l_右;
            // 橋渡し全体をまたぐ「同じ」リードを重複させただけで、独立な証拠は1本ぶんしかない
            var l_単一リード = l_正解.Substring(l_左.Length - 30, 30 + l_橋.Length + 30);
            var l_重複リード群 = Enumerable.Repeat(l_単一リード, 50).ToList();

            var l_ギャップ = new 局所ギャップ(0, l_左.Length, l_橋.Length, l_左, l_右);
            var l_結果 = LocalAssembler.Get_適応kの局所結果(l_ギャップ, Get_証拠(l_重複リード群), l_k長, out _);

            Assert.Null(l_結果);
        }

        /// <summary>
        /// mate1・mate2 が両方とも橋渡し全体をまたいでいても、同一分子由来なので独立な証拠は1件にしかならない
        /// </summary>
        /// <remarks>
        /// 短いインサートで両 mate が橋渡し領域を読み通した場合、配列自体は (向きが逆なので) 異なる文字列になる<br/>
        /// pair ID を見ずに文字列の違いだけで独立性を数えると、この 2 本を誤って「2 件の独立した支持」として受理してしまう
        /// </remarks>
        [Fact]
        public void V_mate1とmate2は同一分子として1件にしか数えない()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, 921);
            var l_橋 = V_生成_ランダム配列(40, 922);
            var l_右 = V_生成_ランダム配列(300, 923);
            var l_正解 = l_左 + l_橋 + l_右;
            var l_scaffold = this.V_書き出し_スキャフォールド("matepair_dup.fasta", l_左 + new string('N', l_橋.Length) + l_右);

            // 同じ1組のペアの mate1・mate2 が両方とも橋渡し全体をまたぐ (短いインサートで読み通した状況)
            // 内容は異なる (mate2 は逆相補・オフセット違い) が、同一分子由来なので独立な証拠は1件のみ
            var l_mate1 = new List<string> { l_正解.Substring(l_左.Length - 30, 30 + l_橋.Length + 30) };
            var l_mate2 = new List<string> { Util.V_逆相補(l_正解.Substring(l_左.Length - 25, 25 + l_橋.Length + 25)) };
            var (l_read1, l_read2) = this.V_書き出し_リードペア("matepair_dup", l_mate1, l_mate2);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_scaffold, l_read1, l_read2, l_k長);

            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Contains('N', Get_単一配列(l_scaffold));
        }

        /// <summary>
        /// アンカーに直接ヒットしない (橋の中間だけをカバーする) リードでも、
        /// アンカーヒットリードとの重なりを介した近傍thread拡張で回収され、ギャップを埋められることを確かめる
        /// </summary>
        /// <remarks>
        /// P3b: anchor直接hitのみから、近傍のread threadへと段階的に回収する
        /// </remarks>
        [Fact]
        public void V_近傍thread拡張でアンカーに届かないリードも拾って埋める()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, p_シード: 1001);
            var l_橋 = V_生成_ランダム配列(200, p_シード: 1002);
            var l_右 = V_生成_ランダム配列(300, p_シード: 1003);
            var l_正解 = l_左 + l_橋 + l_右;

            // 80bp のリードがアンカーと 21bp 以上重ならないと k-mer 索引にヒットしないため、
            // 1回目 (アンカー直接ヒット) で回収できる範囲は左右それぞれアンカー境界から最大 79bp までに限られる
            // (260-379 側 / 420-541 側)。ゾーンB (330-470) は橋の中間だけをカバーし、
            // どちらのアンカーにも届かないが、A・Cの1回目回収範囲とはそれぞれ21bpを大きく超えて重なる
            var l_ゾーンA = l_正解.Substring(260, 120);
            var l_ゾーンB = l_正解.Substring(330, 140);
            var l_ゾーンC = l_正解.Substring(420, 120);

            var l_スキャフォールド = this.V_書き出し_スキャフォールド("threadexp.fasta", l_左 + new string('N', l_橋.Length) + l_右);
            var l_リード = this.V_書き出し_リード("threadexp.fq", [l_ゾーンA, l_ゾーンB, l_ゾーンC], p_リード長: 80);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールド, l_リード, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_埋めたギャップ数);
            Assert.Equal(l_正解, Get_単一配列(l_スキャフォールド));
        }

        /// <summary>
        /// 対照実験: 橋の中間だけをカバーするゾーンB が無ければ、アンカーヒットだけでは k-mer に本当の欠落があり埋まらない
        /// </summary>
        /// <remarks>
        /// 上のテストが「拡張回収のおかげで埋まった」ことを示すための対照。
        /// これが無いと、A+C だけで偶然埋まっていた可能性を排除できない
        /// </remarks>
        [Fact]
        public void V_近傍thread拡張が無ければ中間の欠落で埋まらない()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, p_シード: 1001);
            var l_橋 = V_生成_ランダム配列(200, p_シード: 1002);
            var l_右 = V_生成_ランダム配列(300, p_シード: 1003);
            var l_正解 = l_左 + l_橋 + l_右;

            var l_ゾーンA = l_正解.Substring(260, 120);
            var l_ゾーンC = l_正解.Substring(420, 120);

            var l_スキャフォールド = this.V_書き出し_スキャフォールド("threadexp_noB.fasta", l_左 + new string('N', l_橋.Length) + l_右);
            var l_リード = this.V_書き出し_リード("threadexp_noB.fq", [l_ゾーンA, l_ゾーンC], p_リード長: 80);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールド, l_リード, string.Empty, l_k長);

            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールド));
        }

        /// <summary>片側だけがアンカーに当たるペアでは mate も回収して橋へ使う</summary>
        [Fact]
        public void V_アンカーに当たったリードのmateも回収する()
        {
            const int l_k長 = 21;
            var l_左 = V_生成_ランダム配列(300, 601);
            var l_充填 = V_生成_ランダム配列(100, 602);
            var l_右 = V_生成_ランダム配列(300, 603);
            var l_正解 = l_左 + l_充填 + l_右;
            var l_scaffold = this.V_書き出し_スキャフォールド("mate.fasta", l_左 + new string('N', l_充填.Length) + l_右);

            // 接合部の両端 (左アンカー末尾のk-mer から右アンカー先頭のk-mer まで) を
            // 複数の異なるリードで覆うだけの余裕と密度を持たせる (Has経路支持 の基準を満たすため)
            var l_アンカー側 = new List<string>();
            var l_mate側 = new List<string>();
            for (var i = l_左.Length - 40; i <= l_左.Length + l_充填.Length - 10; i += 3)
            {
                l_アンカー側.Add(l_左.Substring(100 + l_アンカー側.Count, 50));
                l_mate側.Add(l_正解.Substring(i, 50));
            }
            var (l_read1, l_read2) = this.V_書き出し_リードペア("mate", l_アンカー側, l_mate側);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_scaffold, l_read1, l_read2, l_k長);

            Assert.Equal(1, l_統計.A_埋めたギャップ数);
            Assert.Equal(l_正解, Get_単一配列(l_scaffold));
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// リードがギャップを実際に跨いでいるとき、真の配列を復元することを検証する
        /// </summary>
        [Fact]
        public void V_リードがギャップを実際に跨ぐときは真の配列を復元する()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 1);
            var l_充填配列 = V_生成_ランダム配列(60, p_シード: 2);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 3);
            var l_正解 = l_接頭配列 + l_充填配列 + l_接尾配列;

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold.fasta", l_接頭配列 + new string('N', l_充填配列.Length) + l_接尾配列);
            var l_リードパス = this.V_書き出し_リード("reads.fq", [l_正解], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(1, l_統計.A_埋めたギャップ数);
            Assert.Equal(l_充填配列.Length, l_統計.A_埋めた塩基数);
            Assert.Equal(l_正解, Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// どちらのアンカー付近にも局所リードが集まらないとき、ギャップを N のまま残すことを検証する
        /// </summary>
        [Fact]
        public void V_どちらのアンカー付近にも局所リードが無いときはギャップをNのまま残す()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 10);
            var l_充填配列 = V_生成_ランダム配列(60, p_シード: 11);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 12);

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_noreads.fasta", l_接頭配列 + new string('N', l_充填配列.Length) + l_接尾配列);
            // まったく無関係な配列からリードを取る
            var l_無関係配列 = V_生成_ランダム配列(500, p_シード: 999);
            var l_リードパス = this.V_書き出し_リード("reads_unrelated.fq", [l_無関係配列], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_局所リードが集まらなかった数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// アンカー付近のリードはあるが、両者を橋渡しする配列 (ギャップの中身) を読んだリードが無い場合、経路が繋がらないので埋められない
        /// </summary>
        [Fact]
        public void V_局所リードがギャップを橋渡ししないときは経路無しとしてNのまま残す()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 20);
            var l_充填配列 = V_生成_ランダム配列(60, p_シード: 21);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 22);

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_nopath.fasta", l_接頭配列 + new string('N', l_充填配列.Length) + l_接尾配列);
            // prefix と suffix それぞれの内部だけを読んだリード (橋渡しは無い)
            var l_リードパス = this.V_書き出し_リード("reads_nopath.fq", [l_接頭配列, l_接尾配列], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_到達できなかった数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// 橋渡しの配列が 2 通りとも読まれている場合、どちらが正しいか決められないので N のまま残す (誤った配列で埋めるより安全)
        /// </summary>
        [Fact]
        public void V_橋渡しの配列が2通り同程度に支持されるときは推測せずNのまま残す()
        {
            const int l_k長 = 21;
            var l_接頭配列 = V_生成_ランダム配列(300, p_シード: 30);
            var l_接尾配列 = V_生成_ランダム配列(300, p_シード: 31);
            var l_主充填配列 = V_生成_ランダム配列(60, p_シード: 32);
            var l_代替充填配列 = V_生成_ランダム配列(60, p_シード: 33);

            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_ambiguous.fasta", l_接頭配列 + new string('N', l_主充填配列.Length) + l_接尾配列);
            var l_リードパス = this.V_書き出し_リード("reads_ambiguous.fq", [l_接頭配列 + l_主充填配列 + l_接尾配列, l_接頭配列 + l_代替充填配列 + l_接尾配列], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(1, l_統計.A_対象ギャップ数);
            Assert.Equal(0, l_統計.A_埋めたギャップ数);
            Assert.Equal(1, l_統計.A_一意に定まらなかった数);
            Assert.Contains('N', Get_単一配列(l_スキャフォールドパス));
        }

        /// <summary>
        /// scaffold にギャップが無いときは、対象ギャップ数をゼロとして報告することを検証する
        /// </summary>
        [Fact]
        public void V_スキャフォールドにギャップが無いときは対象ギャップ数をゼロとして報告する()
        {
            const int l_k長 = 21;
            var l_正解 = V_生成_ランダム配列(300, p_シード: 40);
            var l_スキャフォールドパス = this.V_書き出し_スキャフォールド("scaffold_nogap.fasta", l_正解);
            var l_リードパス = this.V_書き出し_リード("reads_nogap.fq", [l_正解], p_リード長: 100);

            var l_統計 = LocalAssembler.V_充填_ギャップ(l_スキャフォールドパス, l_リードパス, string.Empty, l_k長);

            Assert.Equal(0, l_統計.A_対象ギャップ数);
            Assert.Equal(l_正解, Get_単一配列(l_スキャフォールドパス));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// pair 情報の無い単独読み取りとして、配列群を読取証拠のリストへ変換する
        /// </summary>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>読取証拠のリスト</returns>
        private static List<読取証拠> Get_証拠(IEnumerable<string> p_配列群)
        {
            return [.. p_配列群.Select(x => new 読取証拠(x, ""))];
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_ランダム配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// scaffold を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_スキャフォールド(string p_ファイル名, string p_配列)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new FastaWriter(l_パス);
            l_書き込み.V_書き込み("SCAFFOLD1", p_配列);
            return l_パス;
        }

        /// <summary>
        /// ゲノムを覆うリードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_ゲノム群">元になるゲノム</param>
        /// <param name="p_リード長">リード長</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_リード(string p_ファイル名, IEnumerable<string> p_ゲノム群, int p_リード長)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_連番 = 0;
            foreach (var l_ゲノム in p_ゲノム群)
            {
                for (var i = 0; i + p_リード長 <= l_ゲノム.Length; i++)
                {
                    l_書き込み.WriteLine($"@r{l_連番}");
                    l_書き込み.WriteLine(l_ゲノム.Substring(i, p_リード長));
                    l_書き込み.WriteLine("+");
                    l_書き込み.WriteLine(new string('I', p_リード長));
                    l_連番++;
                }
            }
            return l_パス;
        }

        private (string A_リード1, string A_リード2) V_書き出し_リードペア(string p_接頭辞, IReadOnlyList<string> p_リード1, IReadOnlyList<string> p_リード2)
        {
            Assert.Equal(p_リード1.Count, p_リード2.Count);
            var l_パス1 = Path.Combine(this._作業ディレクトリ, p_接頭辞 + ".1.fq");
            var l_パス2 = Path.Combine(this._作業ディレクトリ, p_接頭辞 + ".2.fq");
            using var l_書き込み1 = new StreamWriter(l_パス1);
            using var l_書き込み2 = new StreamWriter(l_パス2);
            for (var i = 0; i < p_リード1.Count; i++)
            {
                l_書き込み1.WriteLine($"@pair{i}/1\n{p_リード1[i]}\n+\n{new string('I', p_リード1[i].Length)}");
                l_書き込み2.WriteLine($"@pair{i}/2\n{p_リード2[i]}\n+\n{new string('I', p_リード2[i].Length)}");
            }
            return (l_パス1, l_パス2);
        }

        /// <summary>
        /// FASTA から 1 本だけの配列を読み込む
        /// </summary>
        /// <param name="p_パス">読み込むパス</param>
        /// <returns>読み込んだ配列</returns>
        private static string Get_単一配列(string p_パス)
        {
            using var l_読み込み = new FastaReader(p_パス);
            Assert.True(l_読み込み.Has続き());
            return l_読み込み.Get_次の配列().A_配列;
        }

        #endregion

    }
}
