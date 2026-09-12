using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 合成データに対してエラー訂正を実際に走らせ、その効き目を測る
    /// </summary>
    /// <remarks>
    /// tools/simulate_reads.py が出力した、正解の errors.tsv 付きの合成データを使う<br/>
    /// 注入したエラーのうち何割を正しく真の塩基へ戻せたか (recall) と、逆に正しかった塩基を誤って書き換えてしまった割合 (誤訂正率) を測る
    /// </remarks>
    /// <remarks>
    /// 合成データが存在しない場合はスキップする (通常の CI/dotnet test の対象外、手動で tools/simulate_reads.pyを実行した後に手動で実行する想定)
    /// </remarks>
    public class ErrorCorrectorGroundTruthValidation
    {
        #region 定数

        // Bash tool 経由 (Git Bash/MSYS) で python tools/simulate_reads.py --out-dir /tmp/tsumiki_synth
        // を実行した場合の実際の出力先 (MSYS が/tmp をこの Windows パスへ解決する)
        // .NET のファイル API は MSYS のパス変換を経由しないため、Windows 形式で直接指定する

        /// <summary>
        /// 合成データを置くディレクトリ
        /// </summary>
        private static readonly string _合成データディレクトリ = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "Temp", "tsumiki_synth");

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 正解データに対してエラー訂正の recall と誤訂正率を測定する
        /// </summary>
        [Fact]
        public void V_正解データに対して訂正精度を測定する()
        {
            var l_参照パス = Path.Combine(_合成データディレクトリ, "reference.fasta");
            var l_リード1パス = Path.Combine(_合成データディレクトリ, "reads.1.fq");
            var l_リード2パス = Path.Combine(_合成データディレクトリ, "reads.2.fq");
            var l_エラーパス = Path.Combine(_合成データディレクトリ, "errors.tsv");
            if (!File.Exists(l_参照パス) || !File.Exists(l_リード1パス) || !File.Exists(l_エラーパス))
            {
                return; // 合成データ未生成、tools/simulate_reads.py --out-dir /tmp/tsumiki_synth で生成してから実行する
            }

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 31, A_kmerカットオフ = 2UL, A_スレッド数 = 8 };

            var l_出力ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_ec_validation_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_出力ディレクトリ);
            var l_訂正済み1 = Path.Combine(l_出力ディレクトリ, "corrected.1.fq");
            var l_訂正済み2 = Path.Combine(l_出力ディレクトリ, "corrected.2.fq");

            try
            {
                ErrorCorrector.V_訂正_リードファイル(l_リード1パス, l_リード2パス, l_出力ディレクトリ, l_訂正済み1, l_訂正済み2);

                // read_id -> mate -> position -> true_base (注入されたエラーの正解)
                var l_正解エラー = new Dictionary<(string A_リードID, int A_ペア番号, int A_位置), char>();
                foreach (var l_行 in File.ReadLines(l_エラーパス).Skip(1))
                {
                    var l_列 = l_行.Split('\t');
                    var l_リードID = l_列[0];
                    var l_ペア番号 = int.Parse(l_列[1]);
                    var l_位置 = int.Parse(l_列[2]);
                    var l_正解塩基 = l_列[3][0];
                    l_正解エラー[(l_リードID, l_ペア番号, l_位置)] = l_正解塩基;
                }

                var l_元のリード = V_読み込み_リードID別(l_リード1パス, 1);
                foreach (var l_組 in V_読み込み_リードID別(l_リード2パス, 2))
                {
                    l_元のリード[l_組.Key] = l_組.Value;
                }

                var l_訂正数 = 0;
                var l_未訂正数 = 0;
                var l_誤訂正数 = 0; // 元々正しかった塩基を誤って書き換えてしまった数
                var l_変更位置総数 = 0;

                V_検証_ファイル(l_訂正済み1, 1, l_元のリード, l_正解エラー, ref l_訂正数, ref l_未訂正数, ref l_誤訂正数, ref l_変更位置総数);
                V_検証_ファイル(l_訂正済み2, 2, l_元のリード, l_正解エラー, ref l_訂正数, ref l_未訂正数, ref l_誤訂正数, ref l_変更位置総数);

                var l_注入エラー総数 = l_正解エラー.Count;
                var l_recall = l_注入エラー総数 == 0 ? 0.0D : (double)l_訂正数 / l_注入エラー総数;
                var l_誤訂正率 = l_変更位置総数 == 0 ? 0.0D : (double)l_誤訂正数 / l_変更位置総数;

                Console.WriteLine($"Injected errors: {l_注入エラー総数}");
                Console.WriteLine($"Fixed back to true base (recall): {l_訂正数} ({l_recall:P2})");
                Console.WriteLine($"Still wrong (not fixed, or fixed to a different wrong base): {l_未訂正数}");
                Console.WriteLine($"Total positions changed by corrector: {l_変更位置総数}");
                Console.WriteLine($"Of those, changed a previously-CORRECT base to something wrong (false corrections): {l_誤訂正数} ({l_誤訂正率:P2})");

                // 大まかな健全性チェック: recall は意味のある水準まで達し、
                // 誤訂正率は低く抑えられているべき
                Assert.True(l_recall > 0.5D, $"Expected recall > 50%, got {l_recall:P2}");
                Assert.True(l_誤訂正率 < 0.05D, $"Expected false-correction rate < 5%, got {l_誤訂正率:P2}");
            }
            finally
            {
                if (Directory.Exists(l_出力ディレクトリ))
                {
                    Directory.Delete(l_出力ディレクトリ, recursive: true);
                }
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// リードを ID から引ける形で読み込む
        /// </summary>
        /// <param name="p_パス">読み込むパス</param>
        /// <param name="p_ペア番号">ペアのどちら側か</param>
        /// <returns>ID から引けるリード</returns>
        private static Dictionary<(string, int), string> V_読み込み_リードID別(string p_パス, int p_ペア番号)
        {
            var l_結果 = new Dictionary<(string, int), string>();
            using var l_リーダー = new 簡易FASTQ読み込み(p_パス);
            while (l_リーダー.Get_続きがあるか())
            {
                var (l_生ID, l_配列) = l_リーダー.Get_次のリード();
                var l_ID = l_生ID.TrimStart('@').Split('/')[0];
                l_結果[(l_ID, p_ペア番号)] = l_配列;
            }
            return l_結果;
        }

        /// <summary>
        /// 訂正済みのリードが真の配列へ近づいているかを確かめる
        /// </summary>
        /// <param name="p_訂正済みパス">訂正済みリードのパス</param>
        /// <param name="p_ペア番号">ペアのどちら側か</param>
        /// <param name="p_元のリード">訂正前のリード</param>
        /// <param name="p_正解エラー">真の配列</param>
        /// <param name="p_訂正数">正しく訂正できた位置数の累計</param>
        /// <param name="p_未訂正数">訂正できなかった位置数の累計</param>
        /// <param name="p_誤訂正数">正しかった塩基を誤って書き換えた位置数の累計</param>
        /// <param name="p_変更位置総数">訂正で書き換わった位置数の累計</param>
        private static void V_検証_ファイル(string p_訂正済みパス, int p_ペア番号, Dictionary<(string, int), string> p_元のリード, Dictionary<(string A_リードID, int A_ペア番号, int A_位置), char> p_正解エラー, ref int p_訂正数, ref int p_未訂正数, ref int p_誤訂正数, ref int p_変更位置総数)
        {
            using var l_リーダー = new 簡易FASTQ読み込み(p_訂正済みパス);
            while (l_リーダー.Get_続きがあるか())
            {
                var (l_生ID, l_訂正済み配列) = l_リーダー.Get_次のリード();
                var l_リードID = l_生ID.TrimStart('@').Split('/')[0];
                if (!p_元のリード.TryGetValue((l_リードID, p_ペア番号), out var l_元の配列))
                {
                    continue;
                }

                for (var l_位置 = 0; l_位置 < l_訂正済み配列.Length && l_位置 < l_元の配列.Length; l_位置++)
                {
                    var l_エラーだったか = p_正解エラー.TryGetValue((l_リードID, p_ペア番号, l_位置), out var l_正解塩基);
                    var l_変更されたか = l_訂正済み配列[l_位置] != l_元の配列[l_位置];

                    if (l_変更されたか)
                    {
                        p_変更位置総数++;
                    }

                    if (l_エラーだったか)
                    {
                        if (l_訂正済み配列[l_位置] == l_正解塩基)
                        {
                            p_訂正数++;
                        }
                        else
                        {
                            p_未訂正数++;
                        }
                    }
                    else if (l_変更されたか)
                    {
                        // 元々エラーではなかった (=正しかった) 位置を書き換えてしまった
                        p_誤訂正数++;
                    }
                }
            }
        }

        #endregion

    }

    /// <summary>
    /// FASTQ を「id, 配列」の 2 行単位として読むだけの軽量リーダー (品質行は無視)
    /// </summary>
    /// <param name="p_パス"></param>
    internal sealed class 簡易FASTQ読み込み(string p_パス) : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 読み込み中のストリーム
        /// </summary>
        private readonly StreamReader _読み込み = new(p_パス);

        #endregion

        #region 公開メソッド

        /// <summary>
        /// まだ読めるリードがあるか
        /// </summary>
        /// <returns>続きがあれば true</returns>
        public bool Get_続きがあるか() => !this._読み込み.EndOfStream;

        /// <summary>
        /// 次のリードの ID と配列を返す
        /// </summary>
        /// <returns></returns>
        public (string A_ID, string A_配列) Get_次のリード()
        {
            var l_ID = this._読み込み.ReadLine()!;
            var l_配列 = this._読み込み.ReadLine()!;
            _ = this._読み込み.ReadLine(); // '+'
            _ = this._読み込み.ReadLine(); // quality
            return (l_ID, l_配列);
        }

        /// <summary>
        /// 読み込みストリームを解放する
        /// </summary>
        public void Dispose() => this._読み込み.Dispose();

        #endregion

    }

}
