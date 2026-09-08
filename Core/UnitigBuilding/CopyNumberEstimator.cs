using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Tsumiki.Model.UnitigBuilding;
using Tsumiki.Utility;

namespace Tsumiki.Core.UnitigBuilding
{
    /// <summary>
    /// 各 unitig のコピー数をカバレッジから推定する。
    /// n 回現れる配列にはリードが n 倍集まるので、平均カバレッジ / 基準値 を丸める。
    ///
    /// 反復かどうかをグラフの形ではなく量的な根拠で判定できる点が要点。
    /// 入次数2・出次数2でも単一コピー(バブルの残骸)でありうるし、
    /// 次数1でも高カバレッジならタンデムリピートを1本に潰している疑いがある。
    /// 経路探索では「この unitig を何回まで使ってよいか」の予算にもなる。
    ///
    /// 基準値は長さ加重中央値。単純平均や単純中央値だと本数の多い短い断片に
    /// 引きずられ、ゲノムの大部分を占める単一コピー領域の水準から外れる。
    /// </summary>
    internal static class CopyNumberEstimator
    {
        /// <summary>
        /// これを下回るカバレッジ比の unitig は、コピー数を推定できるだけの
        /// 根拠が無いとみなして 1 として扱う(0 コピーにはしない)。
        /// </summary>
        private const double 多コピーとみなす比の下限 = 1.5;

        /// <summary>
        /// コピー数の上限。これを超える比が出た場合、rRNA オペロンのような
        /// 高コピー反復か、あるいはカバレッジ異常のどちらかで区別がつかない。
        /// 経路探索の予算としては大きすぎると探索が発散するため頭打ちにする。
        /// </summary>
        private const int コピー数の上限 = 12;

        /// <summary>
        /// 「染色体側の確定成分と繋がりが無い、独立した島」を単一の複製単位
        /// (プラスミド等)とみなすために要求する、島の合計長の下限。
        /// AssemblyStatsReporter の「比較可能」しきい値と同じ 500bp を使う。
        /// 短い島は偶然の孤立(flanking配列がtrimで消えた等)や短い反復配列との
        /// 区別がつきにくいため、この下限より短い場合は判定しない
        /// (Unicycler が短いセグメントの単一コピー確定に厳しい条件を課している
        /// のと同じ理由)。
        /// </summary>
        private const int 孤立複製単位とみなす最小合計長 = 500;

        /// <summary>
        /// unitig ID(1始まり) -> その unitig を構成する k-mer の平均カバレッジ、を計算する。
        /// </summary>
        public static Dictionary<int, double> Get_カバレッジ(
            TrustedKmerIndex p_kmerインデックス,
            IReadOnlyDictionary<int, string> p_ユニティグ配列,
            int p_k長)
        {
            Dictionary<int, double> l_カバレッジ = [];
            foreach (var (l_ID, l_配列) in p_ユニティグ配列)
            {
                if (l_配列.Length < p_k長)
                {
                    l_カバレッジ[l_ID] = 0;
                    continue;
                }

                var l_塩基列 = Util.V_変換_塩基列(l_配列);

                ulong l_合計 = 0;
                var l_件数 = 0;
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    l_合計 += p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長));
                    l_件数++;
                }
                l_カバレッジ[l_ID] = l_件数 == 0 ? 0 : (double)l_合計 / l_件数;
            }
            return l_カバレッジ;
        }

        /// <summary>
        /// カバレッジからコピー数を推定する。
        ///
        /// p_グラフ を渡すと、大域基準値との比だけでは「多コピー」に見える
        /// unitig を、排他的な鎖(分岐の無い一続きの隣接)で繋がった近傍の
        /// カバレッジと比較し直す(Get_修正_接続による単一コピー再判定 参照)。
        /// これにより、染色体全体とは異なるカバレッジ水準を持つプラスミドの
        /// 単一コピー領域を、反復と誤判定しにくくなる。渡さない場合は
        /// 従来どおり大域基準値との比だけで判定する。
        /// </summary>
        public static コピー数推定結果 Get_推定結果(
            IReadOnlyDictionary<int, double> p_カバレッジ,
            IReadOnlyDictionary<int, int> p_ユニティグ長,
            UnitigGraph? p_グラフ = null)
        {
            // k-mer スペクトルの2成分混合モデルが適合できていれば、その単一コピー平均を
            // 基準値に使う。カットオフと同じモデルから導くことで、unitig の長さ加重
            // 中央値という別のヒューリスティックとの食い違いを無くす
            // (KmerCutoffSelector.V_解決_kmerカットオフ 参照)。適合に失敗している場合は
            // 従来どおり unitig カバレッジの長さ加重中央値にフォールバックする。
            var l_モデル基準値 = ConfigurationManager.A_スペクトルモデル?.A_単一コピー平均;
            var l_基準値 = l_モデル基準値 is { } l_値 && l_値 > 0
                ? l_値
                : Get_長さ加重中央値(p_カバレッジ, p_ユニティグ長);

            Dictionary<int, int> l_コピー数 = [];
            foreach (var (l_ID, l_カバレッジ値) in p_カバレッジ)
            {
                if (l_基準値 <= 0)
                {
                    l_コピー数[l_ID] = 1;
                    continue;
                }

                var l_比 = l_カバレッジ値 / l_基準値;
                if (l_比 < 多コピーとみなす比の下限)
                {
                    // 単一コピー(あるいは低カバレッジで判断できない)。
                    // 0 にはしない: 実際に配列は存在しており、経路から
                    // 締め出してしまうと組み立てられなくなる。
                    l_コピー数[l_ID] = 1;
                    continue;
                }

                l_コピー数[l_ID] = Math.Clamp((int)Math.Round(l_比), 1, コピー数の上限);
            }

            if (p_グラフ is { } l_グラフ)
            {
                V_修正_孤立した複製単位を単一コピーとみなす(l_グラフ, p_カバレッジ, p_ユニティグ長, l_コピー数);
                V_修正_接続による単一コピー再判定(l_グラフ, p_カバレッジ, l_コピー数);
            }

            return new コピー数推定結果(l_基準値, p_カバレッジ, l_コピー数);
        }

        /// <summary>
        /// 大域基準値との比では多コピーに見える unitig 群のうち、染色体側の
        /// (大域基準値と一致する、確定済みの)成分とグラフ上まったく繋がりが
        /// 無い「島」を成すものを、独立した複製単位(高コピープラスミド等)の
        /// 単一コピー領域とみなす。
        ///
        /// 分散型の反復配列(rRNAオペロン等)は複数の異なるゲノム上の文脈を
        /// 前後に持つため、通常は染色体側の成分と繋がった分岐点になる
        /// (孤立した島にはならない)。したがって「染色体と繋がりが無い、
        /// 内部でカバレッジが一貫した島」という条件は、反復の誤判定を
        /// 招きにくい。
        /// </summary>
        private static void V_修正_孤立した複製単位を単一コピーとみなす(
            UnitigGraph p_グラフ,
            IReadOnlyDictionary<int, double> p_カバレッジ,
            IReadOnlyDictionary<int, int> p_ユニティグ長,
            Dictionary<int, int> p_コピー数)
        {
            var l_成分ID = Get_連結成分(p_グラフ, p_コピー数.Keys);

            HashSet<int> l_確定済み成分 = [.. p_コピー数
                .Where(x => x.Value <= 1)
                .Select(x => l_成分ID[x.Key])];

            var l_未確定の島一覧 = p_コピー数.Keys
                .GroupBy(id => l_成分ID[id])
                .Where(g => !l_確定済み成分.Contains(g.Key));

            foreach (var l_島 in l_未確定の島一覧)
            {
                var l_島の合計長 = l_島.Sum(id => (long)p_ユニティグ長.GetValueOrDefault(id, 0));
                if (l_島の合計長 < 孤立複製単位とみなす最小合計長)
                {
                    continue;
                }

                var l_島内カバレッジ = l_島
                    .Select(id => p_カバレッジ.GetValueOrDefault(id, 0.0))
                    .Where(x => x > 0)
                    .ToList();
                if (l_島内カバレッジ.Count == 0)
                {
                    continue;
                }

                var l_局所基準値 = StatsUtil.Get_中央値(l_島内カバレッジ);
                if (l_局所基準値 <= 0)
                {
                    continue;
                }

                var l_内部で一貫しているか = l_島.All(id =>
                {
                    var l_値 = p_カバレッジ.GetValueOrDefault(id, 0.0);
                    return l_値 <= 0 || l_値 / l_局所基準値 < 多コピーとみなす比の下限;
                });
                if (!l_内部で一貫しているか)
                {
                    // 島の内部でもカバレッジ水準がばらついている
                    // (=島の中に反復がある)。触らない。
                    continue;
                }

                foreach (var id in l_島)
                {
                    p_コピー数[id] = 1;
                }
            }
        }

        /// <summary>
        /// unitig をグラフ上の連結成分に分ける(向きは無視し、辺があれば
        /// 繋がっているとみなす)。辺 v→w があれば逆鎖対称性より w^1→v^1 も
        /// あるため、各 unitig の両頂点(順鎖・逆鎖)の出辺だけを辿れば
        /// 入ってくる辺も含めて全方向を辿ったことになる。
        /// </summary>
        private static Dictionary<int, int> Get_連結成分(UnitigGraph p_グラフ, IEnumerable<int> p_ユニティグID一覧)
        {
            Dictionary<int, int> l_成分ID = [];
            var l_次の成分ID = 0;

            foreach (var l_開始ID in p_ユニティグID一覧)
            {
                if (l_成分ID.ContainsKey(l_開始ID))
                {
                    continue;
                }

                var l_現在の成分ID = l_次の成分ID++;
                Queue<int> l_キュー = new([l_開始ID]);
                l_成分ID[l_開始ID] = l_現在の成分ID;

                while (l_キュー.Count > 0)
                {
                    var l_ID = l_キュー.Dequeue();
                    foreach (var l_頂点 in (ReadOnlySpan<int>)[2 * l_ID, (2 * l_ID) + 1])
                    {
                        foreach (var l_隣接頂点 in p_グラフ.A_出辺[l_頂点])
                        {
                            var l_隣接ID = l_隣接頂点 >> 1;
                            if (l_成分ID.TryAdd(l_隣接ID, l_現在の成分ID))
                            {
                                l_キュー.Enqueue(l_隣接ID);
                            }
                        }
                    }
                }
            }

            return l_成分ID;
        }

        /// <summary>
        /// 大域基準値との比では多コピーに見える unitig を、接続構造を使って
        /// 再判定する(Unicycler の copy depth propagation の簡略版)。
        ///
        /// 対象の unitig から、分岐の無い(出次数1かつ行き先の入次数も1という
        /// 意味で排他的な)辺だけを辿って両方向に伸ばせるだけ伸ばし、
        /// 到達できた unitig 群を「排他的成分」とする。この成分が2本以上から
        /// なり、かつ成分内でのカバレッジの中央値に対する対象の比が
        /// 多コピーとみなす比の下限 を下回るなら、大域基準値とは水準が
        /// 違うだけの単一コピー領域(高コピープラスミドの背骨など)と判断し、
        /// コピー数を1に修正する。
        ///
        /// 排他的な辺だけを辿るため、途中に本物の分岐(反復の入口・合流)が
        /// あれば成分はそこで止まる。したがって成分内のカバレッジが
        /// 実際に反復を含んでいれば、その反復自身は今回の対象にならない限り
        /// 誤って巻き込まれない。
        /// </summary>
        private static void V_修正_接続による単一コピー再判定(
            UnitigGraph p_グラフ, IReadOnlyDictionary<int, double> p_カバレッジ, Dictionary<int, int> p_コピー数)
        {
            foreach (var l_ID in p_コピー数.Keys.ToList())
            {
                if (p_コピー数[l_ID] <= 1)
                {
                    continue;
                }
                if (!p_カバレッジ.TryGetValue(l_ID, out var l_自身のカバレッジ) || l_自身のカバレッジ <= 0)
                {
                    continue;
                }

                var l_成分 = Get_排他的成分(p_グラフ, l_ID);
                if (l_成分.Count < 2)
                {
                    // 排他的に繋がる相手がいない。比較材料が無いので判定しない。
                    continue;
                }

                var l_成分内カバレッジ = l_成分
                    .Select(x => p_カバレッジ.GetValueOrDefault(x, 0.0))
                    .Where(x => x > 0)
                    .ToList();
                if (l_成分内カバレッジ.Count < 2)
                {
                    continue;
                }

                var l_局所基準値 = StatsUtil.Get_中央値(l_成分内カバレッジ);
                if (l_局所基準値 <= 0)
                {
                    continue;
                }

                if (l_自身のカバレッジ / l_局所基準値 < 多コピーとみなす比の下限)
                {
                    p_コピー数[l_ID] = 1;
                }
            }
        }

        /// <summary>
        /// 指定した unitig から、排他的な辺(出次数1かつ行き先の入次数も1)だけを
        /// 両方向へ辿って到達できる unitig ID の集合(自分自身を含む)を返す。
        /// 両方の頂点(順鎖・逆鎖)から辿ることで、鎖を両方向に伸ばす。
        /// </summary>
        private static HashSet<int> Get_排他的成分(UnitigGraph p_グラフ, int p_ユニティグID)
        {
            var l_開始1 = 2 * p_ユニティグID;
            var l_開始2 = l_開始1 ^ 1;

            HashSet<int> l_訪問済み頂点 = [l_開始1, l_開始2];
            HashSet<int> l_結果 = [p_ユニティグID];
            Queue<int> l_キュー = new([l_開始1, l_開始2]);

            while (l_キュー.Count > 0)
            {
                var l_現在 = l_キュー.Dequeue();
                var l_出辺 = p_グラフ.A_出辺[l_現在];
                if (l_出辺.Count != 1)
                {
                    continue;
                }

                var l_次 = l_出辺[0];
                if (p_グラフ.Get_入次数(l_次) != 1 || !l_訪問済み頂点.Add(l_次))
                {
                    continue;
                }

                _ = l_結果.Add(l_次 >> 1);
                if (l_訪問済み頂点.Add(l_次 ^ 1))
                {
                    l_キュー.Enqueue(l_次 ^ 1);
                }
                l_キュー.Enqueue(l_次);
            }

            return l_結果;
        }

        /// <summary>
        /// 長さで重み付けしたカバレッジの中央値。ゲノムの大部分を占める
        /// 単一コピー領域の水準を推定するために使う。
        /// </summary>
        private static double Get_長さ加重中央値(
            IReadOnlyDictionary<int, double> p_カバレッジ,
            IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            var l_組 = p_カバレッジ
                .Where(x => p_ユニティグ長.ContainsKey(x.Key) && x.Value > 0)
                .Select(x => ((long)p_ユニティグ長[x.Key], x.Value));
            return StatsUtil.Get_長さ加重中央値(l_組);
        }

        /// <summary>
        /// 推定結果の要約をコンソールへ出力する。
        /// 「単一コピーが何本・何bp、2コピー以上が何本・何bp」が分かると、
        /// 反復配列がアセンブリのどれだけを占めているかが把握できる。
        /// </summary>
        public static void V_出力_推定結果(コピー数推定結果 p_推定結果, IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            Logger.V_出力(メッセージID.単一コピー基準値, p_推定結果.A_単一コピー基準値);

            var l_コピー数別 = p_推定結果.A_コピー数
                .GroupBy(x => x.Value)
                .OrderBy(x => x.Key)
                .Select(x => (A_コピー数: x.Key,
                              A_本数: x.Count(),
                              A_塩基数: x.Sum(y => (long)p_ユニティグ長.GetValueOrDefault(y.Key, 0))))
                .ToList();

            var l_要約 = string.Join(", ", l_コピー数別.Select(x => $"x{x.A_コピー数}: {x.A_本数} unitig(s)/{x.A_塩基数:N0}bp"));
            Logger.V_出力(メッセージID.コピー数の要約, l_要約);

            var l_反復塩基数 = l_コピー数別.Where(x => x.A_コピー数 >= 2).Sum(x => x.A_塩基数);
            var l_総塩基数 = l_コピー数別.Sum(x => x.A_塩基数);
            if (l_総塩基数 > 0)
            {
                Logger.V_出力(メッセージID.反復配列の割合, l_反復塩基数, l_総塩基数, 100.0 * l_反復塩基数 / l_総塩基数);
            }
        }
    }
}
