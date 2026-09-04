using System.Collections.Concurrent;
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

            using var l_キュー = new BlockingCollection<リードデータ>(boundedCapacity: l_スレッド数 * 64);

            var l_ワーカー = new Task[l_スレッド数];
            for (var w = 0; w < l_スレッド数; w++)
            {
                var l_ワーカー番号 = w;
                l_ワーカー[w] = Task.Run(() =>
                {
                    foreach (var l_リード in l_キュー.GetConsumingEnumerable())
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
                    }
                });
            }

            using (var l_読み込み = new FastqReader(p_ファイルパス))
            {
                while (l_読み込み.Get_続きがあるか())
                {
                    l_キュー.Add(l_読み込み.Get_次のリード_軽量());
                }
            }
            l_キュー.CompleteAdding();

            Task.WaitAll(l_ワーカー);

            Console.WriteLine($"Loaded {l_総リード数} reads from {Path.GetFileName(p_ファイルパス)}");
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
