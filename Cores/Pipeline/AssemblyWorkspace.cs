namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 作業ディレクトリの中間成果物を管理する
    /// </summary>
    internal static partial class AssemblyWorkspace
    {
        #region 正規表現

        /// <summary>
        /// 中間生成物・成果物判定用正規表現 (ファイル)
        /// </summary>
        /// <returns></returns>
        [System.Text.RegularExpressions.GeneratedRegex(@"^(preprocessed|corrected)\.[12]\.fq(\.sha256(\.tmp)?)?$")]
        private static partial System.Text.RegularExpressions.Regex _成果物Regex_ファイル();

        /// <summary>
        /// 中間生成物・成果物判定用正規表現 (ディレクトリ)
        /// </summary>
        /// <returns></returns>
        [System.Text.RegularExpressions.GeneratedRegex(@"^(k|anchor)[0-9]+$")]
        private static partial System.Text.RegularExpressions.Regex _成果物Regex_ディレクトリ();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 作業ディレクトリから中間ファイルだけを消す
        /// </summary>
        /// <param name="p_作業ディレクトリ">中間ファイルを消す対象の作業ディレクトリ</param>
        /// <remarks>
        /// 最終成果物とログは残す<br/>
        /// 作業ディレクトリは利用者が受け取る成果物の置き場でもあるので、消してよいのは k ごとの途中経過や訂正済みリードのほうだけになる<br/>
        /// ログを残すのは、消す指定をした実行こそ後から確かめる手段がそれしかなくなるため
        /// </remarks>
        internal static void V_削除_中間ファイル(string p_作業ディレクトリ)
        {
            foreach (var l_ディレクトリ in Directory.EnumerateDirectories(p_作業ディレクトリ))
            {
                var l_名前 = Path.GetFileName(l_ディレクトリ);
                if ((_成果物Regex_ディレクトリ().IsMatch(l_名前) || l_名前 is "error_correction" or "validation") && !File.GetAttributes(l_ディレクトリ).HasFlag(FileAttributes.ReparsePoint))
                {
                    Directory.Delete(l_ディレクトリ, recursive: true);
                }
            }
            foreach (var l_ファイル in Directory.EnumerateFiles(p_作業ディレクトリ))
            {
                if (_成果物Regex_ファイル().IsMatch(Path.GetFileName(l_ファイル)) || Path.GetFileName(l_ファイル) is "polished.fasta" or "merged_scaffolds.fasta")
                {
                    File.Delete(l_ファイル);
                }
            }
        }

        #endregion
    }
}
