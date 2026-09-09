using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Core.Preprocessing
{
    /// <summary>
    /// FASTQ ファイルを読み進めて TrustedKmerIndex へ k-mer を登録する処理
    /// (曖昧塩基を無視する既定経路)
    /// </summary>
    /// <remarks>
    /// 本パイプラインと ErrorCorrector の
    /// 事前カウントパスの両方から呼べるよう切り出したもの
    /// </remarks>
    internal static class KmerCounting
    {
        /// <summary>
        /// FASTQ を 1 本のスレッドで順に読み進めつつ、ワーカー群へ配って並列に登録する
        /// </summary>
        /// <remarks>
        /// 読み取りを 1 本に保つのはディスク I/O をシーケンシャルなままにするため
        /// </remarks>
        public static void V_読込_リードファイル(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            ulong l_総リード数 = 0UL;
            ulong l_ログ回数 = 0UL;
            var l_カウンタロック = new object();

            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 64,
                Get_リード列(p_ファイルパス),
                (l_リード, l_ワーカー番号) =>
                {
                    V_登録_1リード(l_リード, p_kmerインデックス, l_ワーカー番号);

                    var l_ログ出力するか = false;
                    ulong l_ログ値 = 0UL;
                    lock (l_カウンタロック)
                    {
                        l_総リード数++;
                        if (l_総リード数 % Consts.進捗ログ間隔 == 0)
                        {
                            l_ログ回数++;
                            l_ログ出力するか = true;
                            l_ログ値 = l_ログ回数 * Consts.進捗ログ間隔;
                        }
                    }
                    if (l_ログ出力するか)
                    {
                        Logger.V_出力(メッセージID.リード読込の進捗, l_ログ値);
                    }
                });

            Logger.V_出力(メッセージID.リード読込完了, l_総リード数, Path.GetFileName(p_ファイルパス));
        }

        /// <summary>
        /// リード 1(・指定があればリード 2) を、-ab の有無に応じた経路で
        /// TrustedKmerIndex へ読み込む
        /// </summary>
        /// <remarks>
        /// AssemblyPipeline と MultiKAssembler の
        /// どちらも (単一 k・複数 k の違いだけで) 同じ読み込み手順を必要とするため
        /// ここにまとめる
        /// </remarks>
        public static void V_読込_リードペア(
            Parameters p_引数, TrustedKmerIndex p_kmerインデックス, bool p_進行状況を出力するか = false)
        {
            var l_ペアエンドか = !string.IsNullOrWhiteSpace(p_引数.A_リード2のパス);
            if (p_進行状況を出力するか)
            {
                Logger.V_出力(l_ペアエンドか ? メッセージID.リード1の読込開始 : メッセージID.単一リードの読込開始);
            }
            V_読込_1ファイル(p_引数.A_リード1のパス, p_引数.A_曖昧塩基を許容するか, p_kmerインデックス);

            if (!l_ペアエンドか)
            {
                return;
            }

            if (p_進行状況を出力するか)
            {
                Logger.V_出力(メッセージID.リード2の読込開始);
            }
            V_読込_1ファイル(p_引数.A_リード2のパス, p_引数.A_曖昧塩基を許容するか, p_kmerインデックス);
        }

        /// <summary>
        /// 1 ファイルを読み込んで k-mer を数える
        /// </summary>
        /// <param name="p_パス">読み込むリードのパス</param>
        /// <param name="p_曖昧塩基を許容するか">曖昧塩基を展開して数えるか</param>
        /// <param name="p_kmerインデックス">数え上げ先</param>
        private static void V_読込_1ファイル(string p_パス, bool p_曖昧塩基を許容するか, TrustedKmerIndex p_kmerインデックス)
        {
            if (p_曖昧塩基を許容するか)
            {
                V_読込_リードファイル_曖昧塩基あり(p_パス, p_kmerインデックス);
            }
            else
            {
                V_読込_リードファイル(p_パス, p_kmerインデックス);
            }
        }

        /// <summary>
        /// 曖昧塩基を許容する経路
        /// </summary>
        /// <remarks>
        /// 呼ばれる頻度が低い想定のため未並列
        /// </remarks>
        public static void V_読込_リードファイル_曖昧塩基あり(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス)
        {
            ulong l_件数 = 0UL;
            ulong l_ログ回数 = 0UL;
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_Phredオフセット = ConfigurationManager.A_実行時引数.A_Phredオフセット;
            var l_クオリティカットオフ = ConfigurationManager.A_実行時引数.A_クオリティカットオフ;

            using var l_読み込み = new FastqReader(p_ファイルパス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_リード = l_読み込み.Get_次のリード();
                if (l_リード.A_塩基候補列!.Count < l_k長)
                {
                    continue;
                }
                var l_塩基候補 = CollectionsMarshal.AsSpan(l_リード.A_塩基候補列);
                for (var i = 0; i < l_リード.A_クオリティ.Length; i++)
                {
                    if (l_リード.A_クオリティ[i] - l_Phredオフセット - l_クオリティカットオフ < 0)
                    {
                        l_塩基候補[i] = [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.G, Consts.塩基ID.T];
                    }
                }
                p_kmerインデックス.V_登録_曖昧塩基あり(l_塩基候補[..l_k長], 0);
                for (var i = l_k長; i < l_リード.A_塩基候補列.Count; i++)
                {
                    p_kmerインデックス.V_登録_曖昧塩基あり(l_塩基候補.Slice(i - l_k長 + 1, l_k長), 0);
                }
                if (++l_件数 == Consts.進捗ログ間隔)
                {
                    Logger.V_出力(メッセージID.リード読込の進捗, ++l_ログ回数 * Consts.進捗ログ間隔);
                    l_件数 = 0;
                }
            }
            Logger.V_出力(メッセージID.リード読込完了, (l_ログ回数 * Consts.進捗ログ間隔) + l_件数, Path.GetFileName(p_ファイルパス));
        }

        /// <summary>
        /// FASTQ を順に読み進めてリードを返す
        /// </summary>
        private static IEnumerable<リードデータ> Get_リード列(string p_ファイルパス)
        {
            using var l_読み込み = new FastqReader(p_ファイルパス);
            while (l_読み込み.Get_続きがあるか())
            {
                yield return l_読み込み.Get_次のリード_軽量();
            }
        }

        /// <summary>
        /// 1 リード分の k-mer 抽出・品質フィルタリング・登録
        /// </summary>
        /// <remarks>
        /// 逆相補側を別途登録してはいけない<br/>
        /// TrustedKmerIndex.V_登録 が
        /// 正規形へ寄せて数えるため、二重計上になる
        /// </remarks>
        private static void V_登録_1リード(リードデータ p_リード, TrustedKmerIndex p_kmerインデックス, int p_ワーカー番号)
        {
            var l_塩基列 = p_リード.A_塩基列!;
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            if (l_塩基列.Length < l_k長)
            {
                return;
            }

            var l_Phredオフセット = ConfigurationManager.A_実行時引数.A_Phredオフセット;
            var l_クオリティカットオフ = ConfigurationManager.A_実行時引数.A_クオリティカットオフ;

            var l_低品質数 = 0;
            var l_クオリティ = p_リード.A_クオリティ.AsSpan();
            var l_塩基 = l_塩基列.AsSpan();

            for (var i = 0; i < l_k長; i++)
            {
                if (l_塩基[i] == Consts.無効な塩基 ||
                    l_クオリティ[i] - l_Phredオフセット - l_クオリティカットオフ < 0)
                {
                    l_低品質数++;
                }
            }
            if (l_低品質数 == 0)
            {
                p_kmerインデックス.V_登録(l_塩基[..l_k長], p_ワーカー番号);
            }

            for (var i = l_k長; i < l_塩基列.Length; i++)
            {
                if (l_塩基[i - l_k長] == Consts.無効な塩基 ||
                    l_クオリティ[i - l_k長] - l_Phredオフセット - l_クオリティカットオフ < 0)
                {
                    l_低品質数--;
                }
                if (l_塩基[i] == Consts.無効な塩基 ||
                    l_クオリティ[i] - l_Phredオフセット - l_クオリティカットオフ < 0)
                {
                    l_低品質数++;
                }
                if (l_低品質数 == 0)
                {
                    p_kmerインデックス.V_登録(l_塩基.Slice(i - l_k長 + 1, l_k長), p_ワーカー番号);
                }
            }
        }
    }
}
