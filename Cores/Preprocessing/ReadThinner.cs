using System.Globalization;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// 深すぎる入力のリードのペアを、目標の深度まで一様に間引く
    /// </summary>
    internal static class ReadThinner
    {
        #region 定数

        /// <summary>
        /// 推定ゲノムサイズを求めるときに使う一時ディレクトリの名前
        /// </summary>
        private const string C_推定の作業ディレクトリ名 = "thinning";

        /// <summary>
        /// 書き込みの途中経過を置く一時ファイルの接尾辞
        /// </summary>
        private const string C_一時ファイルの接尾辞 = ".tmp";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 目標の深度へ間引くときの採用割合を求める
        /// </summary>
        /// <param name="p_総塩基数">全ライブラリの全リードの塩基数の和</param>
        /// <param name="p_推定ゲノムサイズ">推定ゲノムサイズ (bp)</param>
        /// <param name="p_目標深度">目標の深度 (x)</param>
        /// <returns>採用割合、間引かないときは null</returns>
        public static double? Get_採用割合(long p_総塩基数, long p_推定ゲノムサイズ, int p_目標深度)
        {
            if (p_目標深度 <= 0 || p_推定ゲノムサイズ <= 0 || p_総塩基数 <= 0)
            {
                return null;
            }

            var l_深度 = (double)p_総塩基数 / p_推定ゲノムサイズ;
            if (l_深度 <= p_目標深度 * Consts.間引きを始める倍率)
            {
                return null;
            }

            return p_目標深度 / l_深度;
        }

        /// <summary>
        /// 0 始まりの番号のリードを採るかを決める (乱数を使わず、割合に従って一様に採る)
        /// </summary>
        /// <param name="p_番号">0 始まりの番号</param>
        /// <param name="p_割合">採用割合</param>
        /// <returns>採るなら真</returns>
        public static bool Is採用(long p_番号, double p_割合)
        {
            return (long)Math.Floor((p_番号 + 1) * p_割合) > (long)Math.Floor(p_番号 * p_割合);
        }

        /// <summary>
        /// 1 ライブラリのリードを間引いて書き出す
        /// </summary>
        /// <param name="p_順リード">入力の順リード</param>
        /// <param name="p_逆リード">入力の逆リード、単一リードなら null か空</param>
        /// <param name="p_出力1">出力の順リード</param>
        /// <param name="p_出力2">出力の逆リード、単一リードなら null</param>
        /// <param name="p_割合">採用割合</param>
        /// <param name="p_Is圧縮">真なら出力を gzip で書く</param>
        /// <returns>元の件数 (ペアなら組の数) と採用件数</returns>
        public static (long A_元の件数, long A_採用件数) V_間引き_ファイル(string p_順リード, string? p_逆リード, string p_出力1, string? p_出力2, double p_割合, bool p_Is圧縮 = false)
        {
            var l_Hasペア = !string.IsNullOrWhiteSpace(p_逆リード);
            var l_直接書く = 中間データ置き場.A_Is有効;
            var l_書込先1 = l_直接書く ? p_出力1 : p_出力1 + C_一時ファイルの接尾辞;
            var l_書込先2 = l_Hasペア && p_出力2 is not null ? (l_直接書く ? p_出力2 : p_出力2 + C_一時ファイルの接尾辞) : null;

            long l_元の件数 = 0;
            long l_採用件数 = 0;
            using (var l_読込1 = new FastqReader(p_順リード))
            using (var l_読込2 = l_Hasペア ? new FastqReader(p_逆リード!) : null)
            using (var l_書込1 = new FastqWriter(l_書込先1, p_Is圧縮))
            using (var l_書込2 = l_書込先2 is null ? null : new FastqWriter(l_書込先2, p_Is圧縮))
            {
                while (l_読込1.Has続き())
                {
                    var l_レコード1 = l_読込1.Get_次のレコード();
                    var l_採る = Is採用(l_元の件数, p_割合);
                    if (l_読込2 is null)
                    {
                        if (l_採る)
                        {
                            l_書込1.V_書き込み(l_レコード1.A_ID, l_レコード1.A_配列, l_レコード1.A_クオリティ);
                            l_採用件数++;
                        }

                        l_元の件数++;
                        continue;
                    }

                    if (!l_読込2.Has続き())
                    {
                        throw new InvalidDataException($"{p_逆リード}: 順リード {p_順リード} より先に終わっている。");
                    }

                    var l_レコード2 = l_読込2.Get_次のレコード();
                    if (l_採る)
                    {
                        l_書込1.V_書き込み(l_レコード1.A_ID, l_レコード1.A_配列, l_レコード1.A_クオリティ);
                        l_書込2!.V_書き込み(l_レコード2.A_ID, l_レコード2.A_配列, l_レコード2.A_クオリティ);
                        l_採用件数++;
                    }

                    l_元の件数++;
                }

                if (l_読込2 is not null && l_読込2.Has続き())
                {
                    throw new InvalidDataException($"{p_順リード}: 逆リード {p_逆リード} より先に終わっている。");
                }
            }

            if (!l_直接書く)
            {
                File.Move(l_書込先1, p_出力1, overwrite: true);
                if (l_書込先2 is not null)
                {
                    File.Move(l_書込先2, p_出力2!, overwrite: true);
                }
            }

            return (l_元の件数, l_採用件数);
        }

        /// <summary>
        /// 全ライブラリの全リードの塩基数の和を求める
        /// </summary>
        /// <param name="p_引数">入力の設定</param>
        /// <returns>総塩基数</returns>
        public static long Get_総塩基数(Parameters p_引数)
        {
            long l_総数 = 0;
            foreach (var (l_順リード, l_逆リード) in p_引数.A_ライブラリ群)
            {
                l_総数 += FastqReader.Get_生リード列(l_順リード, l_逆リード).Sum(x => (long)x.Length);
            }

            return l_総数;
        }

        /// <summary>
        /// 入力の k-mer の計数から推定ゲノムサイズを求める
        /// </summary>
        /// <param name="p_引数">入力の設定 (推定に使う k 長はエラー訂正 k 長)</param>
        /// <param name="p_一時ディレクトリ">作業用の置き場</param>
        /// <returns>推定ゲノムサイズ (bp)、推定できなければ null</returns>
        public static long? Get_推定ゲノムサイズ(Parameters p_引数, string p_一時ディレクトリ)
        {
            var l_元設定 = ConfigurationManager.A_実行時引数;
            var l_元モデル = ConfigurationManager.A_スペクトルモデル;
            var l_設定 = p_引数.Get_複製();
            l_設定.Set_推定k長(p_引数.A_エラー訂正k長);
            l_設定.A_Is曖昧塩基許容 = false;
            var l_作業パス = Path.Combine(p_一時ディレクトリ, C_推定の作業ディレクトリ名);
            try
            {
                ConfigurationManager.A_実行時引数 = l_設定;
                _ = Directory.CreateDirectory(l_作業パス);
                using var l_索引 = new TrustedKmerIndex(l_作業パス);
                KmerCounting.V_読込_リードペア(l_設定, l_索引, p_Is進行状況出力: false);
                return KmerHistogram.Get_解析結果(l_索引.Get_出現回数ヒストグラム())?.A_推定ゲノムサイズ;
            }
            finally
            {
                ConfigurationManager.A_実行時引数 = l_元設定;
                ConfigurationManager.A_スペクトルモデル = l_元モデル;
                if (Directory.Exists(l_作業パス))
                {
                    Directory.Delete(l_作業パス, recursive: true);
                }
            }
        }

        /// <summary>
        /// 間引きの採用割合を、1 行のテキストとして書き出す (文化に依らず "R" 書式)
        /// </summary>
        /// <param name="p_パス">書き出すパス</param>
        /// <param name="p_採用割合">採用割合</param>
        internal static void V_書き出し_採用割合(string p_パス, double p_採用割合)
        {
            File.WriteAllText(p_パス, p_採用割合.ToString("R", CultureInfo.InvariantCulture) + "\n");
        }

        /// <summary>
        /// V_書き出し_採用割合 で書いた採用割合を読む
        /// </summary>
        /// <param name="p_パス">読むパス</param>
        /// <returns>採用割合、ファイルが無い・読めないときは null</returns>
        internal static double? Get_読込_採用割合(string p_パス)
        {
            if (!File.Exists(p_パス))
            {
                return null;
            }

            var l_行 = File.ReadAllText(p_パス).Trim();
            return double.TryParse(l_行, NumberStyles.Float, CultureInfo.InvariantCulture, out var l_割合) ? l_割合 : null;
        }

        #endregion
    }
}
