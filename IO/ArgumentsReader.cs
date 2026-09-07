using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.IO
{
    internal class ArgumentsReader
    {
        public static Parameters Get_実行時引数(string[] p_引数列)
        {
            var l_引数 = new Parameters();
            try
            {
                var l_位置 = 0;
                while (l_位置 < p_引数列.Length)
                {
                    var l_キー = p_引数列[l_位置++];
                    switch (l_キー)
                    {
                        case Consts.引数キー.リード1のパス:
                            l_引数.A_リード1のパス = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.リード2のパス:
                            l_引数.A_リード2のパス = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.k長:
                            l_引数.Set_k長一覧(
                                p_引数列[l_位置++].Split(',', StringSplitOptions.RemoveEmptyEntries
                                    | StringSplitOptions.TrimEntries).Select(int.Parse));
                            break;

                        case Consts.引数キー.kmerカットオフ:
                            l_引数.A_kmerカットオフ = ulong.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.Phredオフセット:
                            l_引数.A_Phredオフセット = int.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.クオリティカットオフ:
                            l_引数.A_クオリティカットオフ = int.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.メモリ予算:
                            l_引数.A_メモリ予算 = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.インサートサイズ:
                            l_引数.A_インサートサイズ = int.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.一時ディレクトリ:
                            l_引数.A_一時ディレクトリ = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.スレッド数:
                            l_引数.A_スレッド数 = int.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.ペア結合閾値:
                            l_引数.A_ペア結合閾値 = decimal.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.ペア支持数閾値:
                            l_引数.A_ペア支持数閾値 = ulong.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.ヘルプ:
                            l_引数.A_ヘルプモードか = true;
                            break;

                        case Consts.引数キー.曖昧塩基を許容:
                            l_引数.A_曖昧塩基を許容するか = true;
                            break;

                        case Consts.引数キー.エラー訂正:
                            l_引数.A_エラー訂正するか = true;
                            break;

                        case Consts.引数キー.前処理:
                            l_引数.A_前処理するか = true;
                            break;

                        case Consts.引数キー.マルチk:
                            l_引数.A_マルチkか = true;
                            break;

                        case Consts.引数キー.マージ:
                            l_引数.A_マージするか = true;
                            break;

                        case Consts.引数キー.引き継ぎなし:
                            l_引数.A_引き継ぐか = false;
                            break;

                        case Consts.引数キー.SuperRead:
                            l_引数.A_SuperReadを作るか = true;
                            break;

                        case Consts.引数キー.反復r_mer検証:
                            l_引数.A_反復をrMerで検証するか = true;
                            break;

                        case Consts.引数キー.局所アセンブリ:
                            l_引数.A_局所アセンブリするか = true;
                            break;

                        case Consts.引数キー.積極性モード:
                            V_適用_積極性モード(l_引数, p_引数列[l_位置++]);
                            break;

                        default:
                            Logger.V_出力_警告(Logger.Get_メソッド名(), new ArgumentException($"Unknown argment: {l_キー}"));
                            break;
                    }
                }
            }
            catch (Exception ex)
            {
                Logger.V_出力_エラー(Logger.Get_メソッド名(), ex);
                Environment.Exit(1);
            }

            // ヘルプ表示だけを求められている場合は、リードパスの必須チェックを行わない。
            // (以前は -h のみを指定してもここで「Please set read path」エラーになり
            //  ヘルプが表示できなかった)
            if (l_引数.A_ヘルプモードか)
            {
                return l_引数;
            }

            if (string.IsNullOrWhiteSpace(l_引数.A_リード1のパス))
            {
                l_引数.A_リード1のパス = l_引数.A_リード2のパス;
                l_引数.A_リード2のパス = string.Empty;
            }

            if (string.IsNullOrWhiteSpace(l_引数.A_リード1のパス))
            {
                Logger.V_出力_エラー(Logger.Get_メソッド名(), new ArgumentException("Please set read path"));
                Environment.Exit(0);
            }

            return l_引数;
        }

        /// <summary>
        /// -pu(優勢閾値)と-pc(支持数閾値)を、完全性/正確性のどちらに倒すかの
        /// 1軸で束ねて適用する。個別に -pu/-pc を後ろに書けばそちらで上書きできる
        /// (通常の CLI 引数と同じく、後に書いたものが勝つ)。
        /// </summary>
        private static void V_適用_積極性モード(Parameters p_引数, string p_モード名)
        {
            switch (p_モード名)
            {
                case Consts.積極性モード名.保守的:
                    p_引数.A_ペア結合閾値 = Consts.保守的モードのペア結合閾値;
                    p_引数.A_ペア支持数閾値 = Consts.保守的モードのペア支持数閾値;
                    break;

                case Consts.積極性モード名.標準:
                    p_引数.A_ペア結合閾値 = Consts.ペア結合閾値の既定値;
                    p_引数.A_ペア支持数閾値 = Consts.ペア支持数閾値の既定値;
                    break;

                case Consts.積極性モード名.積極的:
                    p_引数.A_ペア結合閾値 = Consts.積極的モードのペア結合閾値;
                    p_引数.A_ペア支持数閾値 = Consts.積極的モードのペア支持数閾値;
                    break;

                default:
                    throw new ArgumentException(
                        $"Unknown mode \"{p_モード名}\": expected one of " +
                        $"{Consts.積極性モード名.保守的}, {Consts.積極性モード名.標準}, {Consts.積極性モード名.積極的}");
            }
        }
    }
}
