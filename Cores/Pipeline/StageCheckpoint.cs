using System.Security.Cryptography;
using System.Text;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 工程の入力と出力の一致をハッシュで確認する
    /// </summary>
    internal static class StageCheckpoint
    {
        #region 公開メソッド

        /// <summary>
        /// 設定と入力内容から工程の識別値を作る
        /// </summary>
        /// <param name="p_引数">工程開始時の設定</param>
        /// <returns>工程の識別値</returns>
        public static string Get_入力署名(Parameters p_引数)
        {
            var l_設定 = p_引数.Get_複製();
            l_設定.A_Is再開 = false;

            // 入力が前の工程の出力なら、その工程の記録で識別する
            // 記録には入力の署名と出力の内容ハッシュが入っているので、入力そのものを読み直さずに済む
            // (前の工程の中間ファイルを消してあっても再開できる)
            var l_入力の識別 = string.Join("\n", l_設定.A_ライブラリ群.Select(x =>
                Get_保存済み記録(x.A_リード1) ?? (Get_ハッシュ(x.A_リード1) + "\n" + Get_ハッシュ(x.A_リード2))));
            var l_本文 = typeof(StageCheckpoint).Assembly.ManifestModule.ModuleVersionId + "\n" + l_設定 + "\n" + l_入力の識別;
            return Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(l_本文)));
        }

        /// <summary>
        /// その出力を作った工程の記録を返す
        /// </summary>
        /// <param name="p_出力">工程の主出力</param>
        /// <returns>記録の中身、無ければ null</returns>
        public static string? Get_保存済み記録(string? p_出力)
        {
            return !string.IsNullOrWhiteSpace(p_出力) && File.Exists(p_出力 + ".sha256") ? File.ReadAllText(p_出力 + ".sha256") : null;
        }

        /// <summary>
        /// ファイルの内容を識別する
        /// </summary>
        /// <param name="p_パス">入力パス、空なら入力なし</param>
        /// <returns>SHA-256 値または空文字列</returns>
        public static string Get_ハッシュ(string? p_パス)
        {
            if (string.IsNullOrWhiteSpace(p_パス))
            {
                return string.Empty;
            }
            using var l_入力 = File.OpenRead(p_パス);
            return Convert.ToHexString(SHA256.HashData(l_入力));
        }

        /// <summary>
        /// 完了記録が現在の入力と出力に一致するか調べる
        /// </summary>
        /// <remarks>
        /// 入力の識別は Get_入力署名 が済ませているので、ここで見るのは出力の実体と内容だけ
        /// </remarks>
        /// <param name="p_署名">入力署名</param>
        /// <param name="p_出力">工程の主出力</param>
        /// <param name="p_対出力">対になる出力</param>
        /// <returns>完全に一致すれば true</returns>
        public static bool Is再利用可能(string p_署名, string p_出力, string? p_対出力)
        {
            return File.Exists(p_出力 + ".sha256") && File.Exists(p_出力) && (p_対出力 is null || File.Exists(p_対出力)) && File.ReadAllText(p_出力 + ".sha256") == p_署名 + "\n" + Get_ハッシュ(p_出力) + "\n" + Get_ハッシュ(p_対出力);
        }

        /// <summary>
        /// 出力が完成した工程の検証情報を最後に保存する
        /// </summary>
        /// <param name="p_署名">入力署名</param>
        /// <param name="p_出力">工程の主出力</param>
        /// <param name="p_対出力">対になる出力</param>
        public static void V_保存(string p_署名, string p_出力, string? p_対出力)
        {
            var l_一時パス = p_出力 + ".sha256.tmp";
            File.WriteAllText(l_一時パス, p_署名 + "\n" + Get_ハッシュ(p_出力) + "\n" + Get_ハッシュ(p_対出力));
            File.Move(l_一時パス, p_出力 + ".sha256", overwrite: true);
        }

        #endregion
    }
}
