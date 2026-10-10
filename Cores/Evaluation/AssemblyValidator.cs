using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかの自己検査
    /// </summary>
    internal static class AssemblyValidator
    {
        #region 公開メソッド

        /// <summary>
        /// FASTA のアセンブリを、k-mer インデックスが保持する信頼できる k-mer 集合と突き合わせる
        /// </summary>
        /// <param name="p_FASTAパス">突き合わせる FASTA のパス</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_単一コピー基準値">その k-mer が何回現れてよいかをカバレッジから見積もるための基準値</param>
        /// <returns>自己検査の結果</returns>
        public static 整合性検査結果? Get_検査結果(string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長, double p_単一コピー基準値)
        {
            using var l_計測 = new StageTimer($"assembly-validation k={p_k長}");
            Dictionary<UInt128, int> l_観測 = [];
            var l_延べ数 = 0L;
            using (var l_読み込み = new FastaReader(p_FASTAパス))
            {
                while (l_読み込み.Has続き())
                {
                    var l_配列 = l_読み込み.Get_次の配列().A_配列;
                    for (var i = 0; i + p_k長 <= l_配列.Length; i++)
                    {
                        if (!KmerPacking.Is成功_正規化キー(l_配列, i, p_k長, out var l_正規形))
                        {
                            continue;
                        }

                        l_観測[l_正規形] = l_観測.GetValueOrDefault(l_正規形) + 1;
                        l_延べ数++;
                    }
                }
            }

            var (l_信頼kmer数, l_取りこぼし数, l_出しすぎ種類数, l_余分な延べ数) = p_kmerインデックス.Get_信頼kmer一覧().AsParallel().WithDegreeOfParallelism(Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数)).Aggregate(() => (0L, 0L, 0L, 0L), (l_途中, l_kmer) => Get_突き合わせた集計(l_途中, l_kmer, l_観測, p_kmerインデックス, p_単一コピー基準値), (l_左, l_右) => (l_左.Item1 + l_右.Item1, l_左.Item2 + l_右.Item2, l_左.Item3 + l_右.Item3, l_左.Item4 + l_右.Item4), l_合計 => l_合計);
            return new 整合性検査結果(l_信頼kmer数, l_延べ数, l_観測.Count, l_取りこぼし数, l_出しすぎ種類数, l_余分な延べ数);
        }

        /// <summary>
        /// 自己検査の結果をログへ出力する
        /// </summary>
        /// <param name="p_ラベル">どの成果物に対する検査かを表す名前</param>
        /// <param name="p_検査結果">自己検査の結果、調べられなかった場合は null</param>
        public static void V_出力_検査結果(string p_ラベル, 整合性検査結果? p_検査結果)
        {
            if (p_検査結果 is not { } l_結果)
            {
                Logger.V_出力(メッセージID.検査_対象外のk長, p_ラベル);
                return;
            }

            Logger.V_出力(メッセージID.検査_取りこぼし, p_ラベル, l_結果.A_信頼kmer数, l_結果.A_取りこぼし数, l_結果.A_取りこぼし率);
            Logger.V_出力(メッセージID.検査_出しすぎ, p_ラベル, l_結果.A_アセンブリ内の延べ数, l_結果.A_余分な延べ数, l_結果.A_出しすぎ率, l_結果.A_出しすぎkmer種類数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 信頼できる k-mer 1 つをアセンブリ内の出現数と突き合わせ、途中の集計 (信頼 k-mer 数・取りこぼし数・出しすぎ種類数・余分な延べ数) に足す
        /// </summary>
        /// <param name="p_途中">ここまでの集計</param>
        /// <param name="p_kmer">信頼できる k-mer</param>
        /// <param name="p_観測">アセンブリ内の正規化キーごとの出現数</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_単一コピー基準値">その k-mer が何回現れてよいかをカバレッジから見積もるための基準値</param>
        /// <returns>足した後の集計</returns>
        private static (long, long, long, long) Get_突き合わせた集計((long, long, long, long) p_途中, byte[] p_kmer, Dictionary<UInt128, int> p_観測, TrustedKmerIndex p_kmerインデックス, double p_単一コピー基準値)
        {
            var (l_信頼kmer数, l_取りこぼし数, l_出しすぎ種類数, l_余分な延べ数) = p_途中;
            l_信頼kmer数++;
            var l_出現数 = p_観測.GetValueOrDefault(KmerPacking.Get_正規化キー(p_kmer));
            if (l_出現数 == 0)
            {
                return (l_信頼kmer数, l_取りこぼし数 + 1L, l_出しすぎ種類数, l_余分な延べ数);
            }

            if (p_単一コピー基準値 <= 0D)
            {
                return (l_信頼kmer数, l_取りこぼし数, l_出しすぎ種類数, l_余分な延べ数);
            }

            var l_期待コピー数 = Math.Max(1, (int)Math.Round(p_kmerインデックス.Get_カバレッジ(p_kmer) / p_単一コピー基準値));
            return l_出現数 > l_期待コピー数 ? (l_信頼kmer数, l_取りこぼし数, l_出しすぎ種類数 + 1L, l_余分な延べ数 + l_出現数 - l_期待コピー数) : (l_信頼kmer数, l_取りこぼし数, l_出しすぎ種類数, l_余分な延べ数);
        }

        #endregion
    }
}
