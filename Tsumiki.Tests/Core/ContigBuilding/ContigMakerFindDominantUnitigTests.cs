using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ContigMaker.FindDominantUnitig が返す LastMatchEndOffset が、「read 内での位置」ではなく「unitig 内での正しい位置」を指すことを直接検証する
    /// </summary>
    /// <remarks>
    /// 以前は kmerDict が (unitigId のみ) しか保持しておらず、FindDominantUnitig は read 内の最終ヒット k-mer の終端位置 (read 基準) をそのまま unitig 内終端位置として誤用していた<br/>
    /// unitig が read より十分短い間は両者がたまたま近い値になり問題が表面化しなかったが、tip clipping 導入後に unitig が大幅に長くなり、インサートサイズ自動推定 (この値を使う) が明後日の値 (中央値 30 bp 等) を返すようになったことで発覚した
    /// </remarks>
    public class ContigMakerFindDominantUnitigTests : IDisposable
    {
        #region 定数

        // k=8 で内部に k-mer 重複のないことを確認済みの 100 bp 配列

        /// <summary>
        /// 検証に使う唯一のユニティグ
        /// </summary>
        private const string ユニティグ配列 = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAG";

        #endregion

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
        public ContigMakerFindDominantUnitigTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_contigmaker_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

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
        /// 順鎖一致では read の長さではなく unitig 内での終端位置が返る
        /// </summary>
        [Fact]
        public void V_代表ユニティグ_順鎖一致はunitig内座標の終端位置を返す()
        {
            var l_コンティグ構築 = this.Get_コンティグ構築_単一ユニティグ(p_k長: 8);

            // read = unitig の [40,70) 部分 (30 bp)
            // read 自身の長さ (30) ではなく、
            // unitig 内での終端位置 (70) が返るはず
            var l_read = ユニティグ配列.Substring(40, 30);

            var l_ヒット = l_コンティグ構築.Get_代表ユニティグ(l_read);

            Assert.Equal(1, l_ヒット.A_ユニティグID); // 正の値 = 順鎖でのヒット
            Assert.Equal(70, l_ヒット.A_最終一致終端位置);
            Assert.Equal(ユニティグ配列.Length, l_ヒット.A_ユニティグ長);
        }

        /// <summary>
        /// 逆相補一致では unitig を逆向きに見た座標系での終端位置が返る
        /// </summary>
        [Fact]
        public void V_代表ユニティグ_逆相補一致は逆向きunitig座標の終端位置を返す()
        {
            var l_コンティグ構築 = this.Get_コンティグ構築_単一ユニティグ(p_k長: 8);

            // 元の [40,70) を逆相補した read
            // unitig 全体を逆相補した向きで見ると、
            // 元の区間 [40,70) は [100-70, 100-40) = [30,60) に写る
            var l_read = Util.V_逆相補(ユニティグ配列.Substring(40, 30));

            var l_ヒット = l_コンティグ構築.Get_代表ユニティグ(l_read);

            Assert.Equal(-1, l_ヒット.A_ユニティグID); // 負の値 = 逆鎖でのヒット
            Assert.Equal(60, l_ヒット.A_最終一致終端位置);
        }

        /// <summary>
        /// unitig の末尾ちょうどで一致すると終端位置が unitig の全長になる
        /// </summary>
        [Fact]
        public void V_代表ユニティグ_unitig末尾での一致は終端位置がunitig全長になる()
        {
            var l_コンティグ構築 = this.Get_コンティグ構築_単一ユニティグ(p_k長: 8);

            var l_read = ユニティグ配列[^20..]; // unitig の末尾 20 bp

            var l_ヒット = l_コンティグ構築.Get_代表ユニティグ(l_read);

            Assert.Equal(1, l_ヒット.A_ユニティグID);
            Assert.Equal(ユニティグ配列.Length, l_ヒット.A_最終一致終端位置);
            Assert.Equal(0, l_ヒット.A_末尾までの残り長);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ユニティグ 1 本だけを持つコンティグ構築を組み立てる
        /// </summary>
        /// <param name="p_k長">k 長</param>
        /// <returns>組み立てたコンティグ構築</returns>
        private ContigMaker Get_コンティグ構築_単一ユニティグ(int p_k長)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };
            var l_ユニティグパス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            File.WriteAllText(l_ユニティグパス, $">1\n{ユニティグ配列}\n");
            return new ContigMaker(l_ユニティグパス);
        }

        #endregion

    }
}
