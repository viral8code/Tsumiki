using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Pipeline;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Reporting;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// multi-k の k ごとの保存と読み戻しを検証する
    /// </summary>
    public class MultiKCheckpointTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki multik checkpoint
        /// </summary>
        private const string C_項目_tsumiki_multik_checkpoint = "tsumiki_multik_checkpoint_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 塩基配列 ACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGT = "ACGTACGT";

        /// <summary>
        /// 塩基配列 TTTTGGGG
        /// </summary>
        private const string C_塩基配列_TTTTGGGG = "TTTTGGGG";

        /// <summary>
        /// 項目 signature 2
        /// </summary>
        private const string C_項目_signature_2 = "signature-2";

        /// <summary>
        /// 塩基配列 ACGTAC
        /// </summary>
        private const string C_塩基配列_ACGTAC = "ACGTAC";

        /// <summary>
        /// 塩基配列 GGGCCC
        /// </summary>
        private const string C_塩基配列_GGGCCC = "GGGCCC";

        /// <summary>
        /// 項目 U1
        /// </summary>
        private const string C_項目_U1 = "U1+";

        /// <summary>
        /// 項目 abcdef0123456789
        /// </summary>
        private const string C_項目_abcdef0123456789 = "abcdef0123456789";

        /// <summary>
        /// ファイル名 unitigs fasta
        /// </summary>
        private const string C_ファイル名_unitigs_fasta = "unitigs.fasta";

        /// <summary>
        /// ファイル名 contigs fasta
        /// </summary>
        private const string C_ファイル名_contigs_fasta = "contigs.fasta";

        /// <summary>
        /// 書式 a ACGT
        /// </summary>
        private const string C_書式_a_ACGT = ">a\nACGT\n";

        /// <summary>
        /// テストで使う署名
        /// </summary>
        private const string C_署名 = "signature-1";

        #endregion

        #region 内部変数

        /// <summary>
        /// テスト専用の作業パス
        /// </summary>
        private readonly string _作業パス = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_multik_checkpoint + Guid.NewGuid().ToString(C_GUID書式));

        #endregion

        #region コンストラクタ

        /// <summary>
        /// テスト用ディレクトリを用意する
        /// </summary>
        public MultiKCheckpointTests()
        {
            _ = Directory.CreateDirectory(this._作業パス);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// テスト用ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業パス))
            {
                Directory.Delete(this._作業パス, recursive: true);
            }
        }

        /// <summary>
        /// 保存した結果と引き継ぎが、そのまま読み戻せること
        /// </summary>
        [Fact]
        public void V_保存したkの結果と引き継ぎを読み戻す()
        {
            var l_結果 = this.Get_結果();
            List<引き継ぎ配列> l_引き継ぎ = [
                new(C_塩基配列_ACGTACGT, [5, 6, 7, 300_000], 5, A_Is確定経路: true, A_分岐の継ぎ目位置: [1, 3], A_未観測の連続範囲: [(0, 2)]),
                new(C_塩基配列_TTTTGGGG, [1, 1, 1, 1], 5),
            ];
            MultiKCheckpoint.V_保存_k(this._作業パス, C_署名, l_結果, l_引き継ぎ, 3);
            Assert.True(MultiKCheckpoint.Is成功_読込_k(this._作業パス, C_署名, out var l_読んだ結果, out var l_読んだ引き継ぎ, out var l_合成リード数));
            Assert.Equal(l_結果, l_読んだ結果);
            Assert.Equal(3, l_合成リード数);
            Assert.Equal(2, l_読んだ引き継ぎ.Count);
            Assert.Equal(C_塩基配列_ACGTACGT, l_読んだ引き継ぎ[0].A_配列);
            Assert.Equal([5, 6, 7, 300_000], l_読んだ引き継ぎ[0].A_カバレッジ);
            Assert.True(l_読んだ引き継ぎ[0].A_Is確定経路);
            Assert.Equal([1, 3], l_読んだ引き継ぎ[0].A_分岐の継ぎ目位置!);
            Assert.Equal([(0, 2)], l_読んだ引き継ぎ[0].A_未観測の連続範囲!);
            Assert.Null(l_読んだ引き継ぎ[1].A_分岐の継ぎ目位置);
            Assert.Null(l_読んだ引き継ぎ[1].A_未観測の連続範囲);
        }

        /// <summary>
        /// 入力や設定が変わって署名が一致しなければ使わないこと
        /// </summary>
        [Fact]
        public void V_署名が違えば読み戻さない()
        {
            MultiKCheckpoint.V_保存_k(this._作業パス, C_署名, this.Get_結果(), [], 0);
            Assert.False(MultiKCheckpoint.Is成功_読込_k(this._作業パス, C_項目_signature_2, out _, out _, out _));
        }

        /// <summary>
        /// 記録があっても、結果のファイルが無くなっていれば使わないこと
        /// </summary>
        [Fact]
        public void V_結果のファイルが無ければ読み戻さない()
        {
            var l_結果 = this.Get_結果();
            MultiKCheckpoint.V_保存_k(this._作業パス, C_署名, l_結果, [], 0);
            File.Delete(l_結果.A_scaffoldパス!);
            Assert.False(MultiKCheckpoint.Is成功_読込_k(this._作業パス, C_署名, out _, out _, out _));
        }

        /// <summary>
        /// 合成リードは本数まで一致したときだけ読み戻すこと
        /// </summary>
        [Fact]
        public void V_合成リードは本数が一致したときだけ読み戻す()
        {
            List<引き継ぎ配列> l_合成リード = [new(C_塩基配列_ACGTAC, [9, 9], 5), new(C_塩基配列_GGGCCC, [4, 4], 5)];
            MultiKCheckpoint.V_保存_合成リード(this._作業パス, C_署名, l_合成リード);
            Assert.True(MultiKCheckpoint.Is成功_読込_合成リード(this._作業パス, C_署名, 2, out var l_読んだ));
            Assert.Equal([C_塩基配列_ACGTAC, C_塩基配列_GGGCCC], l_読んだ.Select(x => x.A_配列));
            Assert.False(MultiKCheckpoint.Is成功_読込_合成リード(this._作業パス, C_署名, 3, out _));
        }

        /// <summary>
        /// 保存した曖昧箇所の記録を、再開時にその k の記録と履歴へ読み戻すこと
        /// </summary>
        [Fact]
        public void V_曖昧箇所の記録を読み戻す()
        {
            const int l_k長 = 997;
            AmbiguityRecorder.V_開始(l_k長);
            AmbiguityRecorder.V_記録(曖昧箇所の種別.僅差, C_項目_U1, 3.5D, 1.25D, 7L, C_項目_abcdef0123456789);
            var l_記録 = AmbiguityRecorder.Get_記録(l_k長);
            AmbiguityRecorder.V_保存(this._作業パス, l_k長);
            AmbiguityRecorder.V_開始(l_k長);
            AmbiguityRecorder.V_読込(this._作業パス, l_k長);
            Assert.Equal(l_記録, AmbiguityRecorder.Get_記録(l_k長));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ファイルを伴う k の結果を作る
        /// </summary>
        /// <returns>k の結果</returns>
        private アセンブリ実行結果 Get_結果()
        {
            var l_パス群 = new[]
            {
                C_ファイル名_unitigs_fasta,
                C_ファイル名_contigs_fasta,
                Consts.Scaffoldファイル名,
                Consts.GFAファイル名
            }.Select(x => Path.Combine(this._作業パス, x)).ToArray();
            foreach (var l_パス in l_パス群)
            {
                File.WriteAllText(l_パス, C_書式_a_ACGT);
            }

            return new アセンブリ実行結果(63, l_パス群[0], l_パス群[1], l_パス群[2], 11UL, 160.25D, l_パス群[3], new 整合性検査結果(1L, 2L, 3L, 4L, 5L, 6L), コピー数基準の出所.Weighted);
        }

        #endregion
    }
}
