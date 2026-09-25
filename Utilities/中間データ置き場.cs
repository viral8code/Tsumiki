using System.Collections.Concurrent;
using System.IO.Compression;
using Tsumiki.IO;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 中間データ (前処理・訂正済みリード、断片、k-mer 計数の一時ファイル) を読み書きする口
    /// </summary>
    internal static class 中間データ置き場
    {
        #region 定数

        /// <summary>
        /// 1 つの塊の大きさ
        /// </summary>
        internal const int C_塊の大きさ = 16 << 20;

        /// <summary>
        /// ファイルとして開くときのバッファの大きさ
        /// </summary>
        private const int C_ファイルのバッファ = 1 << 20;

        /// <summary>
        /// パスの大文字小文字を区別するかは OS で違う
        /// </summary>
        private static readonly StringComparer C_パスの比較 = OperatingSystem.IsWindows() ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;

        /// <summary>
        /// 名前の比較 (パスの比較と揃える)
        /// </summary>
        private static readonly StringComparison C_名前の比較 = OperatingSystem.IsWindows() ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal;

        #endregion

        #region 内部変数

        /// <summary>
        /// パスごとの圧縮済みの中身
        /// </summary>
        private static readonly ConcurrentDictionary<string, (List<byte[]> A_塊群, long A_元の長さ)> _置き場 = new(C_パスの比較);

        /// <summary>
        /// 取り込んだ入力のキー
        /// </summary>
        private static readonly ConcurrentDictionary<string, byte> _取り込み済み = new(C_パスの比較);

        #endregion

        #region プロパティ

        /// <summary>
        /// 中間データをメモリに置くか
        /// </summary>
        public static bool A_Is有効 { get; set; }

        /// <summary>
        /// いま置き場が抱えている圧縮済みの大きさ (バイト)
        /// </summary>
        public static long A_使用量 => _置き場.Values.SelectMany(static x => x.A_塊群).Sum(static x => (long)x.Length);

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 中間データを書き込むストリームを開く
        /// </summary>
        /// <param name="p_パス">書き込み先</param>
        /// <returns>閉じた時点で書き込みが確定するストリーム</returns>
        public static Stream Get_書込ストリーム(string p_パス)
        {
            塩基列控え.V_無効化(p_パス);
            if (!A_Is有効)
            {
                return new FileStream(p_パス, FileMode.Create, FileAccess.Write, FileShare.Read, C_ファイルのバッファ, FileOptions.SequentialScan);
            }

            var l_キー = Get_キー(p_パス);
            _ = _置き場.TryRemove(l_キー, out _);
            var l_塊 = new 塊書込ストリーム();
            return new 計数ストリーム(new DeflateStream(l_塊, CompressionLevel.Fastest), l_元の長さ => _置き場[l_キー] = (l_塊.A_塊群, l_元の長さ));
        }

        /// <summary>
        /// ディスク上のファイルを置き場へ写す
        /// </summary>
        /// <param name="p_パス">写すファイル</param>
        public static void V_取り込み(string p_パス)
        {
            if (!A_Is有効 || _置き場.ContainsKey(Get_キー(p_パス)))
            {
                return;
            }

            _取り込み済み[Get_キー(p_パス)] = 0;
            using var l_元 = new FileStream(p_パス, FileMode.Open, FileAccess.Read, FileShare.Read, C_ファイルのバッファ, FileOptions.SequentialScan);
            using var l_先 = Get_書込ストリーム(p_パス);
            l_元.CopyTo(l_先);
        }

        /// <summary>
        /// 読み込むストリームを開く
        /// </summary>
        /// <param name="p_パス">読み込み元</param>
        /// <returns>置き場にあればその中身、無ければファイル</returns>
        public static Stream Get_読込ストリーム(string p_パス)
        {
            return _置き場.TryGetValue(Get_キー(p_パス), out var l_項目)
                ? new 計数ストリーム(new DeflateStream(new 塊読込ストリーム(l_項目.A_塊群), CompressionMode.Decompress), l_項目.A_元の長さ)
                : new FileStream(p_パス, FileMode.Open, FileAccess.Read, FileShare.Read, C_ファイルのバッファ, FileOptions.SequentialScan);
        }

        /// <summary>
        /// 置き場かディスクのどちらかにあるか
        /// </summary>
        /// <param name="p_パス">調べるパス</param>
        /// <returns></returns>
        public static bool Is存在(string? p_パス)
        {
            return !string.IsNullOrWhiteSpace(p_パス) && (_置き場.ContainsKey(Get_キー(p_パス)) || File.Exists(p_パス));
        }

        /// <summary>
        /// 置き場かディスクから消す
        /// </summary>
        /// <param name="p_パス">消すパス</param>
        public static void V_削除(string p_パス)
        {
            塩基列控え.V_無効化(p_パス);
            var l_キー = Get_キー(p_パス);
            if (!_置き場.TryRemove(l_キー, out _) && !_取り込み済み.ContainsKey(l_キー) && File.Exists(p_パス))
            {
                File.Delete(p_パス);
            }
        }

        /// <summary>
        /// ディレクトリ直下で、名前が接頭辞で始まり接尾辞で終わるものを置き場とディスクの両方から列挙する
        /// </summary>
        /// <param name="p_ディレクトリ">探す場所</param>
        /// <param name="p_接頭辞">名前の先頭</param>
        /// <param name="p_接尾辞">名前の末尾</param>
        /// <returns></returns>
        public static IReadOnlyList<string> Get_一覧(string p_ディレクトリ, string p_接頭辞, string p_接尾辞)
        {
            var l_場所 = Get_キー(p_ディレクトリ);
            var l_置き場 = _置き場.Keys.Where(x => string.Equals(Path.GetDirectoryName(x), l_場所, C_名前の比較) && Is名前が一致(Path.GetFileName(x), p_接頭辞, p_接尾辞));
            var l_ディスク = Directory.Exists(p_ディレクトリ)
                ? Directory.EnumerateFiles(p_ディレクトリ).Where(x => Is名前が一致(Path.GetFileName(x), p_接頭辞, p_接尾辞))
                : [];
            return [.. l_置き場, .. l_ディスク];
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// 置き場を空にする
        /// </summary>
        public static void V_消去()
        {
            _置き場.Clear();
            _取り込み済み.Clear();
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 同じ場所を指す別表記を 1 つにまとめたキー
        /// </summary>
        /// <param name="p_パス"></param>
        /// <returns></returns>
        private static string Get_キー(string p_パス)
        {
            return Path.GetFullPath(p_パス);
        }

        /// <summary>
        /// 名前が接頭辞で始まり接尾辞で終わるか
        /// </summary>
        /// <param name="p_名前"></param>
        /// <param name="p_接頭辞"></param>
        /// <param name="p_接尾辞"></param>
        /// <returns></returns>
        private static bool Is名前が一致(string p_名前, string p_接頭辞, string p_接尾辞)
        {
            return p_名前.StartsWith(p_接頭辞, C_名前の比較) && p_名前.EndsWith(p_接尾辞, C_名前の比較);
        }

        #endregion

    }
}
