namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 作業ディレクトリの中間成果物を管理する
    /// </summary>
    internal static partial class AssemblyWorkspace
    {
        #region 定数

        /// <summary>
        /// エラー訂正の作業ディレクトリ名
        /// </summary>
        public const string C_エラー訂正ディレクトリ名 = "error_correction";

        /// <summary>
        /// 検査の作業ディレクトリ名
        /// </summary>
        public const string C_検査ディレクトリ名 = "validation";

        /// <summary>
        /// 統合した成果物の接頭辞
        /// </summary>
        public const string C_統合接頭辞 = "merged_";

        /// <summary>
        /// ポリッシュ済み配列のファイル名
        /// </summary>
        private const string C_ポリッシュ済みファイル名 = "polished.fasta";

        /// <summary>
        /// 統合した Scaffold のファイル名
        /// </summary>
        private const string C_統合Scaffoldファイル名 = C_統合接頭辞 + Tsumiki.Commons.Consts.Scaffoldファイル名;

        /// <summary>
        /// 統合した Contig のファイル名
        /// </summary>
        private const string C_統合Contigファイル名 = C_統合接頭辞 + AssemblyPipeline.C_Contigファイル名;

        /// <summary>
        /// 中間ファイルを判定する正規表現
        /// </summary>
        private const string C_成果物Regex_ファイル = @"^(preprocessed|corrected)\.[12]\.fq(\.sha256(\.tmp)?)?$";

        /// <summary>
        /// 中間ディレクトリを判定する正規表現
        /// </summary>
        private const string C_成果物Regex_ディレクトリ = @"^(k|anchor)[0-9]+$";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 作業ディレクトリから中間ファイルだけを消す
        /// </summary>
        /// <param name="p_作業ディレクトリ">中間ファイルを消す対象の作業ディレクトリ</param>
        public static void V_削除_中間ファイル(string p_作業ディレクトリ)
        {
            foreach (var l_ディレクトリ in Directory.EnumerateDirectories(p_作業ディレクトリ))
            {
                var l_名前 = Path.GetFileName(l_ディレクトリ);
                if ((Get_成果物Regex_ディレクトリ().IsMatch(l_名前) || l_名前 is C_エラー訂正ディレクトリ名 or C_検査ディレクトリ名) && !File.GetAttributes(l_ディレクトリ).HasFlag(FileAttributes.ReparsePoint))
                {
                    Directory.Delete(l_ディレクトリ, recursive: true);
                }
            }

            foreach (var l_ファイル in Directory.EnumerateFiles(p_作業ディレクトリ))
            {
                var l_名前 = Path.GetFileName(l_ファイル);
                if (Get_成果物Regex_ファイル().IsMatch(l_名前) || l_名前 is C_ポリッシュ済みファイル名 or C_統合Scaffoldファイル名 or C_統合Contigファイル名)
                {
                    File.Delete(l_ファイル);
                }
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 中間生成物・成果物判定用正規表現 (ファイル)
        /// </summary>
        /// <returns></returns>
        [System.Text.RegularExpressions.GeneratedRegex(C_成果物Regex_ファイル)]
        private static partial System.Text.RegularExpressions.Regex Get_成果物Regex_ファイル();

        /// <summary>
        /// 中間生成物・成果物判定用正規表現 (ディレクトリ)
        /// </summary>
        /// <returns></returns>
        [System.Text.RegularExpressions.GeneratedRegex(C_成果物Regex_ディレクトリ)]
        private static partial System.Text.RegularExpressions.Regex Get_成果物Regex_ディレクトリ();

        #endregion
    }
}
