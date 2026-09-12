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
        public const string バージョン = "0.1";

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
        /// <remarks>
        /// k が 64 を超えると 2 bit パックが UInt128 に収まらず高速経路から外れるため、そこで頭を抑える (偶数を避けるので実際に選ばれる最大値は 63)
        /// </remarks>
        public const int 自動k長の上限 = 63;

        /// <summary>
        /// -mk で試す k の個数
        /// </summary>
        /// <remarks>
        /// 増やすほど実行時間が線形に伸びる<br/>
        /// ヘルプに実行時間の目安として出るため、ここに置いている
        /// </remarks>
        public const int マルチkで試す個数 = 6;

        /// <summary>
        /// 試す価値があるとみなす k-mer カバレッジの下限
        /// </summary>
        /// <remarks>
        /// k を上げると 1 リードから取れる k-mer が減るため、カバレッジの薄いデータで高い k を試しても時間を捨てるだけになる<br/>
        /// 最初の k の結果から各 k のカバレッジを予測し、これを下回るものは実行前に捨てる
        /// </remarks>
        public const double マルチkの最小kmerカバレッジ = 10.0D;

        /// <summary>
        /// k-mer カウント時にメモリ上へ保持するカウントの総量 (バイト)
        /// </summary>
        /// <remarks>
        /// 大きいほどフラッシュ回数が減って I/O が軽くなる代わりメモリを使う
        /// </remarks>
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
        /// 一時ディレクトリの既定値
        /// </summary>
        public const string 一時ディレクトリの既定値 = "temp";

        /// <summary>
        /// 許容 Phred オフセット
        /// </summary>
        public static readonly int[] 許容Phredオフセット = [33, 64];

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
        /// <remarks>
        /// 名前を付ける側 (ContigMaker/Scaffolder) と、それを根拠に数える側 (AssemblyScorer/CircularClosureVerifier/CompletenessValidator) が別々に文字列を持つと、片方だけ変えたときに黙って 0 件になる
        /// </remarks>
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
        /// <remarks>
        /// 画面をどれだけ静かにしても、また一時ディレクトリを消す指定があっても、これだけは残す
        /// </remarks>
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
        /// 1 個の 64 bit 語へ格納する塩基数
        /// </summary>
        public const int 語あたりの塩基数 = 32;

        /// <summary>
        /// ギャップ充填 (GapFiller ・ LocalAssembler 共通) で、推定ギャップ長に対して許容する誤差 (塩基)
        /// </summary>
        /// <remarks>
        /// インサートサイズ推定のばらつきがそのままギャップ長推定のばらつきになるため、ぴったりの長さだけを探すと現実にはまず当たらない
        /// </remarks>
        public const int ギャップ充填の長さの余裕幅 = 30;

        /// <summary>
        /// ギャップ充填 (GapFiller ・ LocalAssembler 共通) の対象とするギャップ長の上限
        /// </summary>
        /// <remarks>
        /// これより長いギャップは探索空間が広すぎるうえ、推定長の誤差も大きく一意に定まる見込みが薄いため対象外とする
        /// </remarks>
        public const int ギャップ充填のギャップ長上限 = 500;

        #endregion

        #region 定数クラス

        /// <summary>
        /// コマンドライン引数のキー
        /// </summary>
        /// <remarks>
        /// 定数をまとめるための入れ子であり、値を持つ型 (エンティティ) ではないため Model 配下へは展開しない
        /// </remarks>
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
            /// エラー訂正
            /// </summary>
            public const string エラー訂正 = "-ec";

            /// <summary>
            /// 前処理
            /// </summary>
            public const string 前処理 = "-pp";

            /// <summary>
            /// メモリ予算
            /// </summary>
            public const string メモリ予算 = "-mem";

            /// <summary>
            /// マルチ k
            /// </summary>
            public const string マルチk = "-mk";

            /// <summary>
            /// マージ
            /// </summary>
            public const string マージ = "-mg";

            /// <summary>
            /// 引き継ぎなし
            /// </summary>
            public const string 引き継ぎなし = "-nc";

            /// <summary>
            /// 合成リードを作る
            /// </summary>
            public const string SuperRead = "-sr";

            /// <summary>
            /// 反復の r-mer 検証
            /// </summary>
            public const string 反復r_mer検証 = "-rv";

            /// <summary>
            /// 局所アセンブリ
            /// </summary>
            public const string 局所アセンブリ = "-la";

            /// <summary>
            /// 積極性モード
            /// </summary>
            public const string 積極性モード = "-mode";

            /// <summary>
            /// GFA 出力
            /// </summary>
            public const string GFA出力 = "-gfa";

            /// <summary>
            /// ポリッシュ
            /// </summary>
            public const string ポリッシュ = "-po";

            /// <summary>
            /// 環状閉鎖検証
            /// </summary>
            public const string 環状閉鎖検証 = "-cc";

            /// <summary>
            /// 救済 kmer
            /// </summary>
            public const string 救済kmer = "-my";

            /// <summary>
            /// 再開
            /// </summary>
            public const string 再開 = "-rs";

            /// <summary>
            /// ログ水準
            /// </summary>
            public const string ログ水準 = "-log";

            #endregion
        }

        /// <summary>
        /// -mode が束ねる値<br/>
        /// </summary>
        /// <remarks>
        /// -lang に指定できる言語名<br/>
        /// 利用者が本当に決めたいのは「完全性と正確性のどちらに倒すか」の 1 軸であり、-pu (優勢閾値) と -pc (支持数閾値) を個別の数値として露出するより、Unicycler の --mode {conservative, normal, bold} のようなプリセットのほうが意図を素直に表せる<br/>
        /// normal は既定値そのもの<br/>
        /// </remarks>
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
        /// <remarks>
        /// 行の種類を表すと同時に、どの水準で画面に出すかの判定にも使う (<see cref="Logger"/>) <br/>
        /// 言語によらず同じ綴りにすること
        /// </remarks>
        public static class ログ目印
        {
            #region 定数

            /// <summary>
            /// 内部の判断過程
            /// </summary>
            /// <remarks>
            /// 既定では画面に出さない
            /// </remarks>
            public const string 詳細 = "[Debug]";

            /// <summary>
            /// 完全長の判定
            /// </summary>
            /// <remarks>
            /// 静かにしていても出す
            /// </remarks>
            public const string 完全性 = "[Complete]";

            /// <summary>
            /// レポートの出力先
            /// </summary>
            /// <remarks>
            /// 静かにしていても出す
            /// </remarks>
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
        /// 塩基の内部表現
        /// </summary>
        /// <remarks>
        /// A/C/G/T は塩基記号そのものなので英字のまま残す (日本語にするとかえって読みにくいため)
        /// </remarks>
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
