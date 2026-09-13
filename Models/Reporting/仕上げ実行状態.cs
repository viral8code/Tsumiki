namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// v0.2 仕上げ経路の実行状態
    /// </summary>
    internal enum 仕上げ実行状態
    {
        /// <summary>
        /// 仕上げ経路を実行していない
        /// </summary>
        未実行,

        /// <summary>
        /// 仕上げ経路を完了した
        /// </summary>
        完了,

        /// <summary>
        /// 証拠不足のため draft を返した
        /// </summary>
        証拠不足,

        /// <summary>
        /// 資源上限のため draft を返した
        /// </summary>
        資源制限,

        /// <summary>
        /// 入力が対象モデルに適合しないため draft を返した
        /// </summary>
        モデル不適合,
    }

    /// <summary>
    /// v0.2 仕上げ経路の状態と理由
    /// </summary>
    /// <param name="A_状態"></param>
    /// <param name="A_理由コード"></param>
    internal readonly record struct 仕上げ実行結果(仕上げ実行状態 A_状態, string? A_理由コード)
    {
        #region 公開メソッド

        /// <summary>
        /// 状態と理由の組を検証する
        /// </summary>
        public void V_検証()
        {
            if (!Enum.IsDefined(this.A_状態))
            {
                throw new ArgumentOutOfRangeException(nameof(this.A_状態), "Unknown finish status");
            }

            if (this.A_状態 is 仕上げ実行状態.未実行 or 仕上げ実行状態.完了)
            {
                if (!string.IsNullOrWhiteSpace(this.A_理由コード))
                {
                    throw new ArgumentException("A completed or disabled finish result must not have a reason code", nameof(this.A_理由コード));
                }
                return;
            }

            if (string.IsNullOrWhiteSpace(this.A_理由コード))
            {
                throw new ArgumentException("A draft finish result must have a reason code", nameof(this.A_理由コード));
            }
        }

        #endregion
    }
}
