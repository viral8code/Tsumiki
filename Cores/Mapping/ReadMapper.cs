using Tsumiki.Commons;
using Tsumiki.Models.Mapping;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Mapping
{
    /// <summary>
    /// 種ヒットのチェインと帯域制限整列でリードを配置する
    /// </summary>
    /// <remarks>
    /// 索引は全参照位置を持つが、各種について保持する反復由来の位置には上限を設ける<br/>
    /// 反復を無制限に展開すると、配列の大部分が低複雑度である入力で処理時間が発散するため
    /// </remarks>
    internal sealed class ReadMapper
    {
        #region 定数

        /// <summary>
        /// 種の長さ
        /// </summary>
        private const int 種長 = 21;

        /// <summary>
        /// 参照で種を取る間隔
        /// </summary>
        private const int 種間隔 = 4;

        /// <summary>
        /// 1 種で保持する位置の上限
        /// </summary>
        private const int 種ヒット上限 = 16;

        /// <summary>
        /// 整列で許容する対角線からのずれ
        /// </summary>
        private const int 帯域幅 = 24;

        /// <summary>
        /// 一致の得点
        /// </summary>
        private const int 一致得点 = 2;

        /// <summary>
        /// 不一致の罰点
        /// </summary>
        private const int 不一致罰点 = -3;

        /// <summary>
        /// ギャップを開く罰点
        /// </summary>
        private const int ギャップ開始罰点 = -5;

        /// <summary>
        /// ギャップを延長する罰点
        /// </summary>
        private const int ギャップ延長罰点 = -1;

        /// <summary>
        /// 配置に必要な最小スコア
        /// </summary>
        private const int 最小スコア = 30;

        #endregion

        #region 内部変数

        /// <summary>
        /// 配置先の配列群
        /// </summary>
        private readonly IReadOnlyList<string> _参照配列群;

        /// <summary>
        /// 種から参照上の位置を引く索引
        /// </summary>
        private readonly Dictionary<UInt128, List<種ヒット>> _種索引;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 参照配列群の種索引を構築する
        /// </summary>
        /// <param name="p_参照配列群"></param>
        public ReadMapper(IReadOnlyList<string> p_参照配列群)
        {
            this._参照配列群 = p_参照配列群;
            this._種索引 = [];
            this.V_構築_種索引();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// リードを最も尤もらしい参照位置へ配置する
        /// </summary>
        /// <param name="p_リード"></param>
        /// <returns>配置結果</returns>
        public リード配置 Get_配置(string p_リード)
        {
            if (p_リード.Length < 種長)
            {
                return リード配置.配置なし;
            }

            var l_候補数 = this.Get_候補数(p_リード);
            if (l_候補数.Count == 0)
            {
                return リード配置.配置なし;
            }

            List<リード配置> l_配置候補 = [];
            foreach (var l_候補 in l_候補数.OrderByDescending(x => x.Value).Take(種ヒット上限))
            {
                l_配置候補.Add(this.Get_整列(p_リード, l_候補.Key));
            }

            var l_最良 = l_配置候補.MaxBy(x => x.A_スコア);

            if (l_最良.A_スコア < 最小スコア)
            {
                return リード配置.配置なし;
            }

            var l_次善スコア = 0;
            foreach (var l_配置 in l_配置候補)
            {
                if (Is異なる配置(l_最良, l_配置))
                {
                    l_次善スコア = Math.Max(l_次善スコア, l_配置.A_スコア);
                }
            }
            var l_信頼度 = Math.Clamp((l_最良.A_スコア - Math.Max(0, l_次善スコア)) * 3, 0, 60);
            return l_最良 with { A_信頼度 = l_信頼度 };
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種索引を構築する
        /// </summary>
        private void V_構築_種索引()
        {
            for (var i = 0; i < this._参照配列群.Count; i++)
            {
                var l_配列 = this._参照配列群[i];
                for (var j = 0; j + 種長 <= l_配列.Length; j += 種間隔)
                {
                    if (!KmerPacking.TryGet_パック(l_配列, j, 種長, out var l_順鎖))
                    {
                        continue;
                    }

                    this.V_登録_種(l_順鎖, new 種ヒット(i, j, false));
                    this.V_登録_種(KmerPacking.Get_逆相補(l_順鎖, 種長), new 種ヒット(i, j, true));
                }
            }
        }

        /// <summary>
        /// 種を索引へ登録する
        /// </summary>
        /// <param name="p_種"></param>
        /// <param name="p_ヒット"></param>
        private void V_登録_種(UInt128 p_種, 種ヒット p_ヒット)
        {
            if (!this._種索引.TryGetValue(p_種, out var l_ヒット群))
            {
                l_ヒット群 = [];
                this._種索引[p_種] = l_ヒット群;
            }

            if (l_ヒット群.Count < 種ヒット上限)
            {
                l_ヒット群.Add(p_ヒット);
            }
        }

        /// <summary>
        /// 種の対角線ごとに候補を数える
        /// </summary>
        /// <param name="p_リード"></param>
        /// <returns>候補とそれを支持する種の数</returns>
        private Dictionary<(int A_配列番号, bool A_Is逆鎖, int A_対角線), int> Get_候補数(string p_リード)
        {
            Dictionary<(int A_配列番号, bool A_Is逆鎖, int A_対角線), int> l_候補数 = [];
            for (var i = 0; i + 種長 <= p_リード.Length; i++)
            {
                if (!KmerPacking.TryGet_パック(p_リード, i, 種長, out var l_種) || !this._種索引.TryGetValue(l_種, out var l_ヒット群))
                {
                    continue;
                }

                foreach (var l_ヒット in l_ヒット群)
                {
                    var l_リード位置 = l_ヒット.A_Is逆鎖 ? p_リード.Length - i - 種長 : i;
                    var l_候補 = (l_ヒット.A_配列番号, l_ヒット.A_Is逆鎖, l_ヒット.A_参照位置 - l_リード位置);
                    l_候補数[l_候補] = l_候補数.GetValueOrDefault(l_候補) + 1;
                }
            }
            return l_候補数;
        }

        /// <summary>
        /// 2 つの配置が別の場所にあるかを返す
        /// </summary>
        /// <param name="p_基準"></param>
        /// <param name="p_比較対象"></param>
        /// <returns>別の配置なら true</returns>
        private static bool Is異なる配置(リード配置 p_基準, リード配置 p_比較対象)
        {
            if (p_基準.A_配列番号 != p_比較対象.A_配列番号 || p_基準.A_Is逆鎖 != p_比較対象.A_Is逆鎖)
            {
                return true;
            }

            if (p_基準.A_整列位置群.Count == 0 || p_比較対象.A_整列位置群.Count == 0)
            {
                return true;
            }

            var l_基準の開始 = p_基準.A_整列位置群.Min(x => x.A_参照位置);
            var l_基準の終了 = p_基準.A_整列位置群.Max(x => x.A_参照位置);
            var l_比較対象の開始 = p_比較対象.A_整列位置群.Min(x => x.A_参照位置);
            var l_比較対象の終了 = p_比較対象.A_整列位置群.Max(x => x.A_参照位置);
            var l_重なり = Math.Min(l_基準の終了, l_比較対象の終了) - Math.Max(l_基準の開始, l_比較対象の開始) + 1;
            return l_重なり * 2 < Math.Min(p_基準.A_整列位置群.Count, p_比較対象.A_整列位置群.Count);
        }

        /// <summary>
        /// 候補の近傍で半大域整列を行う
        /// </summary>
        /// <param name="p_リード"></param>
        /// <param name="p_候補"></param>
        /// <returns>整列した配置</returns>
        private リード配置 Get_整列(string p_リード, (int A_配列番号, bool A_Is逆鎖, int A_対角線) p_候補)
        {
            var l_照合リード = p_候補.A_Is逆鎖 ? Util.V_逆相補_曖昧塩基あり(p_リード) : p_リード;
            var l_参照 = this._参照配列群[p_候補.A_配列番号];
            var l_開始 = Math.Max(0, p_候補.A_対角線 - 帯域幅);
            var l_終了 = Math.Min(l_参照.Length, p_候補.A_対角線 + l_照合リード.Length + 帯域幅);
            if (l_終了 - l_開始 < 種長)
            {
                return リード配置.配置なし;
            }

            var l_幅 = l_終了 - l_開始;
            var l_得点 = Enumerable.Repeat(int.MinValue / 4, (l_照合リード.Length + 1) * (l_幅 + 1)).ToArray();
            var l_挿入得点 = Enumerable.Repeat(int.MinValue / 4, l_得点.Length).ToArray();
            var l_削除得点 = Enumerable.Repeat(int.MinValue / 4, l_得点.Length).ToArray();
            var l_経路 = new byte[l_得点.Length];
            for (var j = 0; j <= l_幅; j++)
            {
                l_得点[j] = 0;
                l_削除得点[j] = 0;
            }

            for (var i = 1; i <= l_照合リード.Length; i++)
            {
                l_得点[i * (l_幅 + 1)] = ギャップ開始罰点 + (i - 1) * ギャップ延長罰点;
                l_挿入得点[i * (l_幅 + 1)] = l_得点[i * (l_幅 + 1)];
                l_経路[i * (l_幅 + 1)] = 1;
            }

            for (var i = 1; i <= l_照合リード.Length; i++)
            {
                var l_中心 = i + 帯域幅;
                var l_左 = Math.Max(1, l_中心 - 帯域幅 * 2);
                var l_右 = Math.Min(l_幅, l_中心 + 帯域幅 * 2);
                for (var j = l_左; j <= l_右; j++)
                {
                    var l_添字 = i * (l_幅 + 1) + j;
                    var l_左上添字 = (i - 1) * (l_幅 + 1) + j - 1;
                    var l_上添字 = (i - 1) * (l_幅 + 1) + j;
                    var l_左添字 = i * (l_幅 + 1) + j - 1;
                    var l_対角 = l_得点[l_左上添字] + (l_照合リード[i - 1] == l_参照[l_開始 + j - 1] ? 一致得点 : 不一致罰点);
                    l_挿入得点[l_添字] = Math.Max(l_得点[l_上添字] + ギャップ開始罰点, l_挿入得点[l_上添字] + ギャップ延長罰点);
                    l_削除得点[l_添字] = Math.Max(l_得点[l_左添字] + ギャップ開始罰点, l_削除得点[l_左添字] + ギャップ延長罰点);
                    l_得点[l_添字] = Math.Max(l_対角, Math.Max(l_挿入得点[l_添字], l_削除得点[l_添字]));
                    l_経路[l_添字] = l_得点[l_添字] == l_対角 ? (byte)0 : l_得点[l_添字] == l_挿入得点[l_添字] ? (byte)1 : (byte)2;
                }
            }

            var l_末尾 = 0;
            for (var j = 1; j <= l_幅; j++)
            {
                if (l_得点[l_照合リード.Length * (l_幅 + 1) + j] > l_得点[l_照合リード.Length * (l_幅 + 1) + l_末尾])
                {
                    l_末尾 = j;
                }
            }

            var l_最終スコア = l_得点[l_照合リード.Length * (l_幅 + 1) + l_末尾];
            List<整列位置> l_位置群 = [];
            for (var i = l_照合リード.Length; i > 0;)
            {
                var l_経路種別 = l_経路[i * (l_幅 + 1) + l_末尾];
                if (l_経路種別 == 0)
                {
                    var l_リード位置 = p_候補.A_Is逆鎖 ? p_リード.Length - i : i - 1;
                    l_位置群.Add(new 整列位置(l_リード位置, l_開始 + l_末尾 - 1));
                    i--;
                    l_末尾--;
                }
                else if (l_経路種別 == 1)
                {
                    i--;
                }
                else
                {
                    l_末尾--;
                }
            }

            l_位置群.Reverse();
            return new リード配置(p_候補.A_配列番号, p_候補.A_Is逆鎖, l_最終スコア, 0, l_位置群);
        }

        #endregion
    }
}
