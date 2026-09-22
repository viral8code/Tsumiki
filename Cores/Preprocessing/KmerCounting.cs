using System.Runtime.InteropServices;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// FASTQ ファイルを読み進めて TrustedKmerIndex へ k-mer を登録する処理 (曖昧塩基を無視する既定経路)
    /// </summary>
    /// <remarks>
    /// 本パイプラインと ErrorCorrector の事前カウントパスの両方から呼べるよう切り出したもの
    /// </remarks>
    internal static class KmerCounting
    {
        #region 公開メソッド

        /// <summary>
        /// FASTQ を 1 本のスレッドで順に読み進めつつ、ワーカー群へ配って並列に登録する
        /// </summary>
        /// <param name="p_ファイルパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <remarks>
        /// 読み取りを 1 本に保つのはディスク I/O をシーケンシャルなままにするため
        /// </remarks>
        public static void V_読込_リードファイル(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス, int p_Phredオフセット)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_総リード数 = 0UL;

            // 128 塩基までは正規形のパック値をリード上で転がして作り、ワーカーごとに束ねて渡す
            var l_束群 = ConfigurationManager.A_実行時引数.A_k長 <= TrustedKmerIndex.パック値のk上限
                ? Enumerable.Range(0, l_スレッド数).Select(_ => new KmerCountBatch(p_kmerインデックス)).ToArray()
                : null;

            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 64, Get_レコード列(p_ファイルパス), (l_レコード, l_ワーカー番号) =>
                {
                    if (l_束群 is not null)
                    {
                        V_登録_1リード_パック値(l_レコード.A_配列, l_レコード.A_クオリティ, l_束群[l_ワーカー番号], p_Phredオフセット);
                    }
                    else
                    {
                        V_登録_1リード(Util.V_変換_塩基列(l_レコード.A_配列), l_レコード.A_クオリティ, p_kmerインデックス, p_Phredオフセット);
                    }

                    var l_件数 = Interlocked.Increment(ref l_総リード数);
                    if (l_件数 % Consts.進捗ログ間隔 == 0UL)
                    {
                        Logger.V_出力(メッセージID.リード読込の進捗, l_件数);
                    }
                });

            foreach (var l_束 in l_束群 ?? [])
            {
                l_束.V_吐き出し();
            }

            Logger.V_出力(メッセージID.リード読込完了, l_総リード数, Path.GetFileName(p_ファイルパス));
        }

        /// <summary>
        /// リード 1 (・指定があればリード 2) を、-ab の有無に応じた経路で TrustedKmerIndex へ読み込む
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_Is進行状況出力"></param>
        /// <remarks>
        /// AssemblyPipeline と MultiKAssembler のどちらも (単一 k ・複数 k の違いだけで) 同じ読み込み手順を必要とするためここにまとめる
        /// </remarks>
        public static void V_読込_リードペア(Parameters p_引数, TrustedKmerIndex p_kmerインデックス, bool p_Is進行状況出力 = false)
        {
            using var l_計測 = new StageTimer($"kmer-count k={p_引数.A_k長}");
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (A_リード1, A_リード2) = p_引数.A_ライブラリ群[i];
                var l_Phred = p_引数.Get_Phredオフセット(i);
                var l_Isペアエンド = !string.IsNullOrWhiteSpace(A_リード2);
                if (p_Is進行状況出力)
                {
                    Logger.V_出力(l_Isペアエンド ? メッセージID.リード1の読込開始 : メッセージID.単一リードの読込開始);
                }
                V_読込_1ファイル(A_リード1, p_引数.A_Is曖昧塩基許容, p_kmerインデックス, l_Phred);

                if (!l_Isペアエンド)
                {
                    continue;
                }

                if (p_Is進行状況出力)
                {
                    Logger.V_出力(メッセージID.リード2の読込開始);
                }
                V_読込_1ファイル(A_リード2, p_引数.A_Is曖昧塩基許容, p_kmerインデックス, l_Phred);
            }
        }

        /// <summary>
        /// 曖昧塩基を許容する経路
        /// </summary>
        /// <param name="p_ファイルパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <remarks>
        /// 呼ばれる頻度が低い想定のため未並列
        /// </remarks>
        public static void V_読込_リードファイル_曖昧塩基あり(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス, int p_Phredオフセット)
        {
            var l_件数 = 0UL;
            var l_ログ回数 = 0UL;
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_Phredオフセット = p_Phredオフセット;
            var l_クオリティカットオフ = ConfigurationManager.A_実行時引数.A_クオリティカットオフ;

            using var l_読み込み = new FastqReader(p_ファイルパス);
            while (l_読み込み.Has続き())
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
                        l_塩基候補[i] = [];
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
                    l_件数 = 0UL;
                }
            }
            Logger.V_出力(メッセージID.リード読込完了, (l_ログ回数 * Consts.進捗ログ間隔) + l_件数, Path.GetFileName(p_ファイルパス));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 ファイルを読み込んで k-mer を数える
        /// </summary>
        /// <param name="p_パス">読み込むリードのパス</param>
        /// <param name="p_Is曖昧塩基許容">曖昧塩基を展開して数えるか</param>
        /// <param name="p_kmerインデックス">数え上げ先</param>
        private static void V_読込_1ファイル(string p_パス, bool p_Is曖昧塩基許容, TrustedKmerIndex p_kmerインデックス, int p_Phredオフセット)
        {
            if (p_Is曖昧塩基許容)
            {
                V_読込_リードファイル_曖昧塩基あり(p_パス, p_kmerインデックス, p_Phredオフセット);
            }
            else
            {
                V_読込_リードファイル(p_パス, p_kmerインデックス, p_Phredオフセット);
            }
        }

        /// <summary>
        /// FASTQ を順に読み進めて塩基列とクオリティを返す
        /// </summary>
        /// <param name="p_ファイルパス"></param>
        /// <returns></returns>
        private static IEnumerable<(string A_配列, string A_クオリティ)> Get_レコード列(string p_ファイルパス)
        {
            using var l_読み込み = new FastqReader(p_ファイルパス);
            while (l_読み込み.Has続き())
            {
                var (_, l_配列, l_クオリティ) = l_読み込み.Get_次のレコード();
                yield return (l_配列, l_クオリティ);
            }
        }

        /// <summary>
        /// 1 リード分の k-mer を、正規形のパック値を転がしながら束へ溜める (k &lt;= 128)
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_クオリティ"></param>
        /// <param name="p_束"></param>
        /// <remarks>
        /// 低品質の塩基は曖昧塩基と同じく窓を切る壁として渡し、V_登録_1リード と同じ k-mer だけを数える
        /// </remarks>
        private static void V_登録_1リード_パック値(string p_配列, string p_クオリティ, KmerCountBatch p_束, int p_Phredオフセット)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            if (p_配列.Length < l_k長)
            {
                return;
            }

            var l_品質下限 = p_Phredオフセット + ConfigurationManager.A_実行時引数.A_クオリティカットオフ;
            var l_窓 = new RollingKmer(l_k長);
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_塩基 = p_クオリティ[i] < l_品質下限 ? 'N' : p_配列[i];
                if (l_窓.Try追加(l_塩基, out var l_キー))
                {
                    p_束.V_追加(l_キー.A_上位, l_キー.A_下位);
                }
            }
        }

        /// <summary>
        /// 1 リード分の k-mer 抽出・品質フィルタリング・登録
        /// </summary>
        /// <param name="p_塩基列"></param>
        /// <param name="p_クオリティ"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <remarks>
        /// 逆相補側を別途登録してはいけない<br/>
        /// TrustedKmerIndex.V_登録 が正規形へ寄せて数えるため、二重計上になる
        /// </remarks>
        private static void V_登録_1リード(byte[] p_塩基列, string p_クオリティ, TrustedKmerIndex p_kmerインデックス, int p_Phredオフセット)
        {
            var l_塩基列 = p_塩基列;
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            if (l_塩基列.Length < l_k長)
            {
                return;
            }

            var l_Phredオフセット = p_Phredオフセット;
            var l_クオリティカットオフ = ConfigurationManager.A_実行時引数.A_クオリティカットオフ;

            var l_低品質数 = 0;
            var l_クオリティ = p_クオリティ.AsSpan();
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
                p_kmerインデックス.V_登録(l_塩基[..l_k長]);
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
                    p_kmerインデックス.V_登録(l_塩基.Slice(i - l_k長 + 1, l_k長));
                }
            }
        }

        #endregion
    }
}
