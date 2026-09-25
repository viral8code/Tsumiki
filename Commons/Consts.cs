namespace Tsumiki.Commons
{
    /// <summary>
    /// 定数定義クラス
    /// </summary>
    internal static class Consts
    {
        #region 定数

        /// <summary>
        /// バージョン
        /// </summary>
        public const string バージョン = "0.4.0";

        /// <summary>
        /// 保守的モードのペア結合閾値
        /// </summary>
        public const decimal 保守的モードのペア結合閾値 = 0.9M;

        /// <summary>
        /// 保守的モードのペア支持数閾値
        /// </summary>
        public const ulong 保守的モードのペア支持数閾値 = 15UL;

        /// <summary>
        /// 積極的モードのペア結合閾値
        /// </summary>
        public const decimal 積極的モードのペア結合閾値 = 0.65M;

        /// <summary>
        /// 積極的モードのペア支持数閾値
        /// </summary>
        public const ulong 積極的モードのペア支持数閾値 = 5UL;

        /// <summary>
        /// 自動選択する k 長の上限
        /// </summary>
        public const int 自動k長の上限 = 63;

        /// <summary>
        /// マルチ k で試す k の個数
        /// </summary>
        public const int マルチkで試す個数 = 6;

        /// <summary>
        /// 試す価値があるとみなす k-mer カバレッジの下限
        /// </summary>
        public const double マルチkの最小kmerカバレッジ = 10D;

        /// <summary>
        /// k-mer カウント時にメモリ上へ保持するカウントの総量 (バイト)
        /// </summary>
        public const long メモリ予算の既定値 = 768L * 1_024L * 1_024L;

        /// <summary>
        /// Phred オフセットの既定値
        /// </summary>
        public const int Phredオフセットの既定値 = 33;

        /// <summary>
        /// クオリティカットオフの既定値
        /// </summary>
        public const int クオリティカットオフの既定値 = 1;

        /// <summary>
        /// 3' 末端の品質トリムの閾値の既定値
        /// </summary>
        public const int 品質トリム閾値の既定値 = 0;

        /// <summary>
        /// 一時ディレクトリの既定値
        /// </summary>
        public const string 一時ディレクトリの既定値 = "temp";

        /// <summary>
        /// 許容 Phred オフセット
        /// </summary>
        public static readonly int[] 許容Phredオフセット = [33, 64,];

        /// <summary>
        /// ペア結合閾値の既定値
        /// </summary>
        public const decimal ペア結合閾値の既定値 = 0.8M;

        /// <summary>
        /// ペア支持数閾値の既定値
        /// </summary>
        public const ulong ペア支持数閾値の既定値 = 10UL;

        /// <summary>
        /// 短い反復解決の拒否権で、経路の接合点が「実際にリードに読まれている」と認めるのに必要な、接合点を跨ぐ r-mer の最小本数
        /// </summary>
        public const int r_mer接合点支持の閾値の既定値 = 4;

        /// <summary>
        /// scaffold ファイル名
        /// </summary>
        public const string Scaffoldファイル名 = "scaffolds.fasta";

        /// <summary>
        /// GFA ファイル名
        /// </summary>
        public const string GFAファイル名 = "assembly.gfa";

        /// <summary>
        /// 環状に閉じた複製単位であることを示す、配列 ID 中の目印
        /// </summary>
        public const string 環状の目印 = "circular";

        /// <summary>
        /// 環状に閉じた経路を「複製単位が 1 周組み上がった」とみなす最小の長さ
        /// </summary>
        public const int 環状として数える最小長 = 1_000;

        /// <summary>
        /// レポートファイル名
        /// </summary>
        public const string レポートファイル名 = "assembly.report.json";

        /// <summary>
        /// 曖昧箇所ファイル名
        /// </summary>
        public const string 曖昧箇所ファイル名 = "assembly.ambiguous.tsv";

        /// <summary>
        /// 実行中に出した内容を全量残すファイル
        /// </summary>
        public const string ログファイル名 = "Tsumiki.log";

        /// <summary>
        /// 進捗ログ間隔
        /// </summary>
        public const ulong 進捗ログ間隔 = 100_000UL;

        /// <summary>
        /// 無効な塩基
        /// </summary>
        public const byte 無効な塩基 = 5;

        /// <summary>
        /// 1 個の 64 bit へ格納する塩基数
        /// </summary>
        public const int ワードあたりの塩基数 = 32;

        /// <summary>
        /// ギャップ充填 (GapFiller・LocalAssembler 共通) で、推定ギャップ長に対して許容する誤差 (塩基)
        /// </summary>
        public const int ギャップ充填の長さの余裕幅 = 30;

        /// <summary>
        /// ギャップ充填 (GapFiller・LocalAssembler 共通) の対象とするギャップ長の上限
        /// </summary>
        public const int ギャップ充填のギャップ長上限 = 500;

        #endregion

        #region 定数クラス

        /// <summary>
        /// コマンドライン引数のキー
        /// </summary>
        public static class 引数キー
        {
            #region 定数

            /// <summary>
            /// リード 1 のパス
            /// </summary>
            public const string リード1のパス = "-1";

            /// <summary>
            /// リード 2 のパス
            /// </summary>
            public const string リード2のパス = "-2";

            /// <summary>
            /// シングルエンドリードのパス
            /// </summary>
            public const string シングルのパス = "-s";

            /// <summary>
            /// k 長
            /// </summary>
            public const string k長 = "-k";

            /// <summary>
            /// kmer カットオフ
            /// </summary>
            public const string kmerカットオフ = "-kc";

            /// <summary>
            /// Phred オフセット
            /// </summary>
            public const string Phredオフセット = "-p";

            /// <summary>
            /// クオリティカットオフ
            /// </summary>
            public const string クオリティカットオフ = "-q";

            /// <summary>
            /// 3' 末端の品質トリムの閾値
            /// </summary>
            public const string 品質トリム閾値 = "-qt";

            /// <summary>
            /// インサートサイズ
            /// </summary>
            public const string インサートサイズ = "-i";

            /// <summary>
            /// 曖昧塩基を許容
            /// </summary>
            public const string 曖昧塩基を許容 = "-ab";

            /// <summary>
            /// ヘルプ
            /// </summary>
            public const string ヘルプ = "-h";

            /// <summary>
            /// 一時ディレクトリ
            /// </summary>
            public const string 一時ディレクトリ = "-t";

            /// <summary>
            /// 一時ディレクトリ削除
            /// </summary>
            public const string 一時ディレクトリ削除 = "-rt";

            /// <summary>
            /// 言語
            /// </summary>
            public const string 言語 = "-lang";

            /// <summary>
            /// バージョン
            /// </summary>
            public const string バージョン = "-v";

            /// <summary>
            /// スレッド数
            /// </summary>
            public const string スレッド数 = "-th";

            /// <summary>
            /// ペア結合閾値
            /// </summary>
            public const string ペア結合閾値 = "-pu";

            /// <summary>
            /// ペア支持数閾値
            /// </summary>
            public const string ペア支持数閾値 = "-pc";

            /// <summary>
            /// エラー訂正を行わない
            /// </summary>
            public const string エラー訂正なし = "-nec";

            /// <summary>
            /// 前処理を行わない
            /// </summary>
            public const string 前処理なし = "-npp";

            /// <summary>
            /// メモリ予算
            /// </summary>
            public const string メモリ予算 = "-mem";

            /// <summary>
            /// マルチ k を行わない
            /// </summary>
            public const string マルチkなし = "-nmk";

            /// <summary>
            /// マージ
            /// </summary>
            public const string マージ = "-mg";

            /// <summary>
            /// 引き継ぎなし
            /// </summary>
            public const string 引き継ぎなし = "-nc";

            /// <summary>
            /// 合成リードを作らない
            /// </summary>
            public const string SuperReadなし = "-nsr";

            /// <summary>
            /// 反復の r-mer 検証を行わない
            /// </summary>
            public const string 反復r_mer検証なし = "-nrv";

            /// <summary>
            /// 局所アセンブリを行わない
            /// </summary>
            public const string 局所アセンブリなし = "-nla";

            /// <summary>
            /// 積極性モード
            /// </summary>
            public const string 積極性モード = "-mode";

            /// <summary>
            /// GFA を出力しない
            /// </summary>
            public const string GFA出力なし = "-ngfa";

            /// <summary>
            /// ポリッシュを行わない
            /// </summary>
            public const string ポリッシュなし = "-npo";

            /// <summary>
            /// 環状閉鎖検証を行わない
            /// </summary>
            public const string 環状閉鎖検証なし = "-ncc";

            /// <summary>
            /// 低頻度 k-mer を救済しない
            /// </summary>
            public const string 救済kmerなし = "-nmy";

            /// <summary>
            /// コピー数基準
            /// </summary>
            public const string コピー数基準 = "-cnb";

            /// <summary>
            /// 低カバレッジ端トリミングなし
            /// </summary>
            public const string 低カバレッジ端トリミングなし = "-nt";

            /// <summary>
            /// 再開
            /// </summary>
            public const string 再開 = "-rs";

            /// <summary>
            /// 中間データをメモリに置く
            /// </summary>
            public const string オンメモリ = "-inmem";

            /// <summary>
            /// ログ水準
            /// </summary>
            public const string ログ水準 = "-log";

            #endregion
        }

        /// <summary>
        /// -lang が束ねる値
        /// </summary>
        public static class 言語名
        {
            #region 定数

            /// <summary>
            /// 英語
            /// </summary>
            public const string 英語 = "en";

            /// <summary>
            /// 日本語
            /// </summary>
            public const string 日本語 = "ja";

            /// <summary>
            /// 中国語
            /// </summary>
            public const string 中国語 = "zh";

            #endregion
        }

        /// <summary>
        /// -log に指定できる水準名
        /// </summary>
        public static class ログ水準名
        {
            #region 定数

            /// <summary>
            /// 最小
            /// </summary>
            public const string 最小 = "quiet";

            /// <summary>
            /// 標準
            /// </summary>
            public const string 標準 = "normal";

            /// <summary>
            /// 詳細
            /// </summary>
            public const string 詳細 = "verbose";

            #endregion
        }

        /// <summary>
        /// 文言の先頭に付ける目印
        /// </summary>
        public static class ログ目印
        {
            #region 定数

            /// <summary>
            /// 内部の判断過程
            /// </summary>
            public const string 詳細 = "[Debug]";

            /// <summary>
            /// 完全長の判定
            /// </summary>
            public const string 完全性 = "[Complete]";

            /// <summary>
            /// レポートの出力先
            /// </summary>
            public const string レポート = "[Report]";

            #endregion
        }

        /// <summary>
        /// -mode に指定できるモード名
        /// </summary>
        public static class 積極性モード名
        {
            #region 定数

            /// <summary>
            /// 保守的
            /// </summary>
            public const string 保守的 = "conservative";

            /// <summary>
            /// 標準
            /// </summary>
            public const string 標準 = "normal";

            /// <summary>
            /// 積極的
            /// </summary>
            public const string 積極的 = "bold";

            #endregion
        }

        /// <summary>
        /// -cnb に指定できる基準
        /// </summary>
        public static class コピー数基準の出所
        {
            #region 定数

            /// <summary>
            /// k-mer スペクトル混合モデル
            /// </summary>
            public const string スペクトラム = "spectrum";

            /// <summary>
            /// unitig カバレッジの長さ加重中央値
            /// </summary>
            public const string 重みづけ = "weighted";

            #endregion
        }

        /// <summary>
        /// 塩基の内部表現
        /// </summary>
        public static class 塩基ID
        {
            #region 定数

            /// <summary>
            /// アデニン
            /// </summary>
            public const byte A = 1;

            /// <summary>
            /// シトシン
            /// </summary>
            public const byte C = 2;

            /// <summary>
            /// グアニン
            /// </summary>
            public const byte G = 3;

            /// <summary>
            /// チミン
            /// </summary>
            public const byte T = 4;

            #endregion
        }

        #endregion
    }
}
