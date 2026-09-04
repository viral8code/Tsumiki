using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// FASTQ ファイルを読み進めて TrustedKmerIndex へ k-mer を登録する処理
    /// (曖昧塩基を無視する既定経路)。本パイプラインと ErrorCorrector の
    /// 事前カウントパスの両方から呼べるよう切り出したもの。
    /// </summary>
    internal static class KmerCounting
    {
        /// <summary>
        /// FASTQ を1本のスレッドで順に読み進めつつ(ディスクI/Oはシーケンシャルなまま)、
        /// 読み取ったリードを BlockingCollection 経由でワーカースレッド群に配る
        /// プロデューサー/コンシューマ方式。
        /// </summary>
        public static void V_読込_リードファイル(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            ulong l_総リード数 = 0;
            ulong l_ログ回数 = 0;
            var l_カウンタロック = new object();

            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 64,
                Get_リード列(p_ファイルパス),
                (l_リード, l_ワーカー番号) =>
                {
                    V_登録_1リード(l_リード, p_kmerインデックス, l_ワーカー番号);

                    var l_ログ出力するか = false;
                    ulong l_ログ値 = 0;
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
                        Console.WriteLine(l_ログ値 + " reads Loaded");
                    }
                });

            Console.WriteLine($"Loaded {l_総リード数} reads from {Path.GetFileName(p_ファイルパス)}");
        }

        /// <summary>
        /// 曖昧塩基を許容する経路。現状シングルスレッドのまま(ワーカー番号固定)で、
        /// 呼ばれる頻度が低い想定のため並列化の優先度を下げている。
        ///
        /// KmerCounting 側と同様、以前はここでもリードを逆相補にしてもう一度
        /// すべての k-mer を登録していた。登録側が正規化するようになった今は
        /// 完全に冗長で、すべてのカウントをちょうど2倍にしてしまうため削除した。
        /// </summary>
        public static void V_読込_リードファイル_曖昧塩基あり(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス)
        {
            ulong l_件数 = 0;
            ulong l_ログ回数 = 0;
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
                    Console.WriteLine((++l_ログ回数 * Consts.進捗ログ間隔) + " reads Loaded");
                    l_件数 = 0;
                }
            }
            Console.WriteLine($"Loaded {(l_ログ回数 * Consts.進捗ログ間隔) + l_件数} reads from {Path.GetFileName(p_ファイルパス)}");
        }
        /// <summary>
        /// FASTQ を順に読み進めてリードを返す。ReadPipeline へ渡す供給元として使う。
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
        /// 1リード分の k-mer 抽出・品質フィルタリング・登録を行う。
        ///
        /// かつてはこの後、リードを逆相補にしてもう一度すべての k-mer を
        /// 登録していた。カウント段階で正規化していなかったため、どちらの
        /// 向きから問い合わせても当たるようにするには両向きを入れる必要が
        /// あったからである。
        ///
        /// 現在は TrustedKmerIndex.V_登録 が順鎖・逆鎖のうち辞書順で小さいほう
        /// (正規化形)に寄せてから数えるため、この2周目は完全に冗長であり、
        /// すべてのカウントをちょうど2倍にしてしまう(実データのヒストグラムが
        /// 偶数のカウントしか持たない、という形で表面化した)。
        /// 削除したことで、登録回数・メモリ・ディスク書き込みがいずれも半減する。
        /// </summary>
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
            var l_クオリティ = p_リード.A_クオリティ.ToCharArray().AsSpan();
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
