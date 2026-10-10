using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// リードファイルからの代表リード長の抽出
    /// </summary>
    public class ReadLengthSnifferTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki readlength tests
        /// </summary>
        private const string C_項目_tsumiki_readlength_tests = "tsumiki_readlength_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 uniform fq
        /// </summary>
        private const string C_ファイル名_uniform_fq = "uniform.fq";

        /// <summary>
        /// ファイル名 trimmed fq
        /// </summary>
        private const string C_ファイル名_trimmed_fq = "trimmed.fq";

        /// <summary>
        /// ファイル名 limited fq
        /// </summary>
        private const string C_ファイル名_limited_fq = "limited.fq";

        /// <summary>
        /// ファイル名 empty fq
        /// </summary>
        private const string C_ファイル名_empty_fq = "empty.fq";

        /// <summary>
        /// ファイル名 pair 1 fq
        /// </summary>
        private const string C_ファイル名_pair_1_fq = "pair.1.fq";

        /// <summary>
        /// ファイル名 pair 2 fq
        /// </summary>
        private const string C_ファイル名_pair_2_fq = "pair.2.fq";

        /// <summary>
        /// ファイル名 single fq
        /// </summary>
        private const string C_ファイル名_single_fq = "single.fq";

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
        public ReadLengthSnifferTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_readlength_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// すべて同じ長さのリードでは、その長さを代表リード長として返すことを検証する
        /// </summary>
        [Fact]
        public void V_すべて同じ長さのリードではその長さを返す()
        {
            var l_パス = this.V_書き出し_FASTQ(C_ファイル名_uniform_fq, 150, 150, 150, 150);
            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(l_パス));
        }

        /// <summary>
        /// トリミング済みのデータではリード長がばらつく
        /// </summary>
        [Fact]
        public void V_トリミング済みデータでは平均や最大ではなく中央値を返す()
        {
            var l_パス = this.V_書き出し_FASTQ(C_ファイル名_trimmed_fq, 35, 40, 148, 150, 150, 150, 151);
            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(l_パス));
        }

        /// <summary>
        /// 標本上限に達したらそこで読み取りを打ち切ることを検証する
        /// </summary>
        [Fact]
        public void V_標本上限に達したら打ち切る()
        {
            var l_長さ一覧 = new[]
            {
                200,
                200
            }.Concat(Enumerable.Repeat(50, 100)).ToArray();
            var l_パス = this.V_書き出し_FASTQ(C_ファイル名_limited_fq, l_長さ一覧);
            Assert.Equal(200, ReadLengthSniffer.Get_代表リード長(l_パス, p_標本上限: 2));
        }

        /// <summary>
        /// 空のファイルでは null を返すことを検証する
        /// </summary>
        [Fact]
        public void V_空ファイルではnullを返す()
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_empty_fq);
            File.WriteAllText(l_パス, string.Empty);
            Assert.Null(ReadLengthSniffer.Get_代表リード長(l_パス));
        }

        /// <summary>
        /// k は「どちらのリードからも k-mer が取れる」必要があるため、ペアで長さが違う場合は短いほうに合わせる
        /// </summary>
        [Fact]
        public void V_ペアで長さが違う場合は短いほうを返す()
        {
            var l_先行パス = this.V_書き出し_FASTQ(C_ファイル名_pair_1_fq, 150, 150, 150);
            var l_後続パス = this.V_書き出し_FASTQ(C_ファイル名_pair_2_fq, 100, 100, 100);
            Assert.Equal(100, ReadLengthSniffer.Get_代表リード長(l_先行パス, l_後続パス));
        }

        /// <summary>
        /// シングルエンドでは、無い方のファイルを無視することを検証する
        /// </summary>
        [Fact]
        public void V_シングルエンドでは無い方のファイルを無視する()
        {
            var l_先行パス = this.V_書き出し_FASTQ(C_ファイル名_single_fq, 150, 150, 150);
            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(l_先行パス, null));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_リード長群">書き出すリードの長さ</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_FASTQ(string p_ファイル名, params int[] p_リード長群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_ID = 1;
            foreach (var l_長さ in p_リード長群)
            {
                l_書き込み.WriteLine($"@read{l_ID++}");
                l_書き込み.WriteLine(new string('A', l_長さ));
                l_書き込み.WriteLine('+');
                l_書き込み.WriteLine(new string('I', l_長さ));
            }

            return l_パス;
        }

        #endregion
    }
}
