namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// v0.2 仕上げ経路の固定設定
    /// </summary>
    internal sealed class 仕上げ設定
    {
        #region 定数

        /// <summary>
        /// 対応するレポート schema
        /// </summary>
        public const int 対応schemaバージョン = 2;

        /// <summary>
        /// 対応する設定 revision
        /// </summary>
        public const int 対応設定revision = 1;

        #endregion

        #region プロパティ

        /// <summary>
        /// v0.2 仕上げ経路を使うか
        /// </summary>
        public bool A_Is有効 { get; init; } = false;

        /// <summary>
        /// レポート schema
        /// </summary>
        public int A_schemaバージョン { get; init; } = 対応schemaバージョン;

        /// <summary>
        /// 設定 revision
        /// </summary>
        public int A_設定revision { get; init; } = 対応設定revision;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 設定がこの実装で扱えることを検証する
        /// </summary>
        public void V_検証()
        {
            if (this.A_schemaバージョン != 対応schemaバージョン)
            {
                throw new NotSupportedException($"Unsupported finish report schema version: {this.A_schemaバージョン}");
            }

            if (this.A_設定revision != 対応設定revision)
            {
                throw new NotSupportedException($"Unsupported finish settings revision: {this.A_設定revision}");
            }
        }

        #endregion
    }
}
