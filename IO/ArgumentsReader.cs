using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.IO
{
    /// <summary>
    /// コマンドライン引数を読み取る
    /// </summary>
    internal class ArgumentsReader
    {
        #region 公開メソッド

        /// <summary>
        /// コマンドライン引数を実行時引数へ組み立てて返す
        /// </summary>
        /// <param name="p_引数列">コマンドライン引数</param>
        /// <returns>実行時引数</returns>
        public static Parameters Get_実行時引数(string[] p_引数列)
        {
            var l_引数 = new Parameters();
            l_引数.V_適用_既定の機能();
            HashSet<string> l_指定済み = [];
            try
            {
                var l_位置 = 0;
                while (l_位置 < p_引数列.Length)
                {
                    var l_キー = p_引数列[l_位置++];
                    _ = l_指定済み.Add(l_キー);
                    switch (l_キー)
                    {
                        case Consts.引数キー.リード1のパス:
                            l_引数.A_リード1のパス = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.リード2のパス:
                            l_引数.A_リード2のパス = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.シングルのパス:
                            l_引数.A_シングルのパス = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.k長:
                            l_引数.Set_k長一覧(p_引数列[l_位置++].Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries).Select(int.Parse));
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

                        case Consts.引数キー.品質トリム閾値:
                            l_引数.A_品質トリム閾値 = int.Parse(p_引数列[l_位置++]);
                            if (l_引数.A_品質トリム閾値 < 0)
                            {
                                throw new ArgumentException($"{Consts.引数キー.品質トリム閾値} must be 0 or greater (0 disables trimming)");
                            }
                            break;

                        case Consts.引数キー.メモリ予算:
                            l_引数.A_メモリ予算 = p_引数列[l_位置++];
                            break;

                        case Consts.引数キー.インサートサイズ:
                            l_引数.A_インサートサイズ = int.Parse(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.一時ディレクトリ削除:
                            l_引数.A_Is一時ディレクトリ削除 = true;
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
                            l_引数.A_Isヘルプモード = true;
                            break;

                        case Consts.引数キー.バージョン:
                            l_引数.A_Isバージョンモード = true;
                            break;

                        case Consts.引数キー.曖昧塩基を許容:
                            l_引数.A_Is曖昧塩基許容 = true;
                            break;

                        case Consts.引数キー.エラー訂正なし:
                            l_引数.A_Isエラー訂正 = false;
                            break;

                        case Consts.引数キー.前処理なし:
                            l_引数.A_Is前処理 = false;
                            break;

                        case Consts.引数キー.マルチkなし:
                            l_引数.A_Isマルチk = false;
                            break;

                        case Consts.引数キー.マージ:
                            l_引数.A_Isマージ = true;
                            break;

                        case Consts.引数キー.引き継ぎなし:
                            l_引数.A_Is引き継ぎ = false;
                            break;

                        case Consts.引数キー.SuperReadなし:
                            l_引数.A_IsSuperRead作成 = false;
                            break;

                        case Consts.引数キー.反復r_mer検証なし:
                            l_引数.A_Is反復rMer検証 = false;
                            break;

                        case Consts.引数キー.局所アセンブリなし:
                            l_引数.A_Is局所アセンブリ = false;
                            break;

                        case Consts.引数キー.言語:
                            l_引数.A_言語 = Get_言語(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.積極性モード:
                            V_適用_積極性モード(l_引数, p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.GFA出力なし:
                            l_引数.A_IsGFA出力 = false;
                            break;

                        case Consts.引数キー.ポリッシュなし:
                            l_引数.A_Isポリッシュ = false;
                            break;

                        case Consts.引数キー.環状閉鎖検証なし:
                            l_引数.A_Is環状閉鎖検証 = false;
                            break;

                        case Consts.引数キー.救済kmerなし:
                            l_引数.A_Is救済kmer使用 = false;
                            break;

                        case Consts.引数キー.コピー数基準:
                            l_引数.A_コピー数基準の出所 = Get_コピー数基準の出所(p_引数列[l_位置++]);
                            break;

                        case Consts.引数キー.低カバレッジ端トリミングなし:
                            l_引数.A_Is低カバレッジ端トリミング = false;
                            break;

                        case Consts.引数キー.再開:
                            l_引数.A_Is再開 = true;
                            break;

                        case Consts.引数キー.オンメモリ:
                            l_引数.A_Isオンメモリ = true;
                            break;

                        case Consts.引数キー.ログ水準:
                            l_引数.A_ログ水準 = Get_ログ水準(p_引数列[l_位置++]);
                            break;

                        default:
                            throw new ArgumentException($"Unknown argument: {l_キー}");
                    }
                }
            }
            catch (Exception l_例外)
            {
                Logger.V_出力_エラー(Logger.Get_メソッド名(), l_例外);
                throw;
            }

            if (l_引数.A_Isヘルプモード || l_引数.A_Isバージョンモード)
            {
                return l_引数;
            }

            if (string.IsNullOrWhiteSpace(l_引数.A_リード1のパス) && !string.IsNullOrWhiteSpace(l_引数.A_リード2のパス))
            {
                l_引数.A_リード1のパス = l_引数.A_リード2のパス;
                l_引数.A_リード2のパス = string.Empty;
            }

            if (l_引数.A_ライブラリ数 == 0)
            {
                Logger.V_出力_エラー(Logger.Get_メソッド名(), new ArgumentException("Please set read path"));
                throw new ArgumentException("Please set read path");
            }

            if (!l_引数.Isライブラリ数が一致)
            {
                var l_例外 = new ArgumentException("-1 and -2 must list the same number of comma-separated libraries");
                Logger.V_出力_エラー(Logger.Get_メソッド名(), l_例外);
                throw l_例外;
            }

            if (Get_相反する指定(l_指定済み, l_引数) is { } l_理由)
            {
                var l_例外 = new ArgumentException(l_理由);
                Logger.V_出力_エラー(Logger.Get_メソッド名(), l_例外);
                throw l_例外;
            }

            return l_引数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 同時に指定すると意味が通らない組を探す
        /// </summary>
        /// <param name="p_指定済み">コマンドラインに書かれたキー</param>
        /// <param name="p_引数">組み立て済みの実行時引数</param>
        /// <returns>見つかれば理由、無ければ null</returns>
        internal static string? Get_相反する指定(IReadOnlySet<string> p_指定済み, Parameters p_引数)
        {
            return p_指定済み.Contains(Consts.引数キー.オンメモリ) && p_指定済み.Contains(Consts.引数キー.メモリ予算)
                ? $"{Consts.引数キー.オンメモリ} and {Consts.引数キー.メモリ予算} cannot be used together: {Consts.引数キー.オンメモリ} keeps k-mer counting runs in memory instead of spilling them to disk under the {Consts.引数キー.メモリ予算} budget"
                : p_指定済み.Contains(Consts.引数キー.オンメモリ) && p_指定済み.Contains(Consts.引数キー.再開)
                ? $"{Consts.引数キー.オンメモリ} and {Consts.引数キー.再開} cannot be used together: {Consts.引数キー.オンメモリ} leaves no intermediate files to resume from"
                : p_指定済み.Contains(Consts.引数キー.マルチkなし) && p_引数.A_k長一覧.Count > 1
                ? $"{Consts.引数キー.マルチkなし} cannot be used with more than one value for {Consts.引数キー.k長}: a comma-separated {Consts.引数キー.k長} asks to try each value and keep the best"
                : p_指定済み.Contains(Consts.引数キー.インサートサイズ) && p_引数.A_ライブラリ群.Count(x => !string.IsNullOrWhiteSpace(x.A_リード2)) > 1
                ? $"{Consts.引数キー.インサートサイズ} cannot be used with more than one paired library: a single insert size would be applied to every library; omit it to estimate each library separately"
                : null;
        }

        /// <summary>
        /// -lang に渡された言語名を解釈する
        /// </summary>
        /// <param name="p_言語名"></param>
        /// <returns></returns>
        private static 言語 Get_言語(string p_言語名)
        {
            return p_言語名 switch
            {
                Consts.言語名.日本語 => 言語.日本語,
                Consts.言語名.英語 => 言語.英語,
                Consts.言語名.中国語 => 言語.中国語,
                _ => throw new ArgumentException($"Unknown language \"{p_言語名}\": expected one of {Consts.言語名.日本語}, {Consts.言語名.英語}, {Consts.言語名.中国語}"),
            };
        }

        /// <summary>
        /// -log に渡された水準名を解釈する
        /// </summary>
        /// <param name="p_水準名"></param>
        /// <returns></returns>
        private static ログ水準 Get_ログ水準(string p_水準名)
        {
            return p_水準名 switch
            {
                Consts.ログ水準名.最小 => ログ水準.最小,
                Consts.ログ水準名.標準 => ログ水準.標準,
                Consts.ログ水準名.詳細 => ログ水準.詳細,
                _ => throw new ArgumentException($"Unknown log level \"{p_水準名}\": expected one of {Consts.ログ水準名.最小}, {Consts.ログ水準名.標準}, {Consts.ログ水準名.詳細}"),
            };
        }

        /// <summary>
        /// -cnb に渡された単一コピー基準の出所を解釈する
        /// </summary>
        /// <param name="p_出所名"></param>
        /// <returns></returns>
        private static コピー数基準の出所 Get_コピー数基準の出所(string p_出所名)
        {
            return p_出所名 switch
            {
                Consts.コピー数基準の出所.スペクトラム => コピー数基準の出所.Spectrum,
                Consts.コピー数基準の出所.重みづけ => コピー数基準の出所.Weighted,
                _ => throw new ArgumentException($"Unknown copy-number baseline \"{p_出所名}\": expected {Consts.コピー数基準の出所.スペクトラム} or {Consts.コピー数基準の出所.重みづけ}"),
            };
        }

        /// <summary>
        /// -pu (優勢閾値) と -pc (支持数閾値) を、完全性/正確性のどちらに倒すかの 1 軸で束ねて適用する
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_モード名"></param>
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
                    throw new ArgumentException($"Unknown mode \"{p_モード名}\": expected one of {Consts.積極性モード名.保守的}, {Consts.積極性モード名.標準}, {Consts.積極性モード名.積極的}");
            }
        }

        #endregion
    }
}
