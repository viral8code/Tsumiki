namespace Tsumiki.Common
{
    /// <summary>
    /// 定数定義クラス
    /// </summary>
    internal static class Consts
    {
        /// <summary>
        /// バージョン
        /// </summary>
        public const string バージョン = "0.1";

        /// <summary>
        /// 本アセンブラの作者
        /// </summary>
        public static readonly List<string> 作者一覧 = [
            "viral",
            ];

        /// <summary>
        /// コマンドライン引数のキー<br/>
        /// 定数をまとめるための入れ子であり、値を持つ型 (エンティティ) ではないため Model 配下へは展開しない
        /// </summary>
        public static class 引数キー
        {
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
        }

        /// <summary>
        /// -mode が束ねる値<br/>
        /// </summary>
        /// <remarks>
        /// -lang に指定できる言語名<br/>
        /// 利用者が本当に決めたいのは「完全性と正確性の
        /// どちらに倒すか」の 1 軸であり、-pu (優勢閾値) と -pc (支持数閾値) を
        /// 個別の数値として露出するより、Unicycler の --mode {conservative, normal, bold} のような
        /// プリセットのほうが意図を素直に表せる<br/>
        /// normal は既定値そのもの<br/>
        /// </remarks>
        public static class 言語名
        {
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
        }

        /// <summary>
        /// -log に指定できる水準名
        /// </summary>
        public static class ログ水準名
        {
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
        }

        /// <summary>
        /// 文言の先頭に付ける目印<br/>
        /// 行の種類を表すと同時に、どの水準で画面に出すかの判定にも使う (<see cref="Logger"/>)<br/>
        /// 言語によらず同じ綴りにすること
        /// </summary>
        public static class ログ目印
        {
            /// <summary>
            /// 内部の判断過程<br/>
            /// 既定では画面に出さない
            /// </summary>
            public const string 詳細 = "[Debug]";

            /// <summary>
            /// 完全長の判定<br/>
            /// 静かにしていても出す
            /// </summary>
            public const string 完全性 = "[Complete]";

            /// <summary>
            /// レポートの出力先<br/>
            /// 静かにしていても出す
            /// </summary>
            public const string レポート = "[Report]";
        }

        /// <summary>
        /// -mode に指定できるモード名
        /// </summary>
        public static class 積極性モード名
        {
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
        }

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
        /// インサートサイズ未指定表示
        /// </summary>
        public const string インサートサイズ未指定表示 = "unspecified";

        /// <summary>
        /// -k も、リード長の標本抽出も当てにできなかった場合の最後の拠り所<br/>
        /// 通常は実際のリード長から自動選択されるため、この値は使われない
        /// </summary>
        public const int k長の既定値 = 31;

        /// <summary>
        /// 自動選択する k 長の、リード長に対する比<br/>
        /// 大きいほど反復を跨げるが、1 リードから取れる k-mer の本数が減ってカバレッジが痩せる
        /// </summary>
        public const double 自動k長のリード長比 = 0.6D;

        /// <summary>
        /// 自動選択する k 長の上限<br/>
        /// k が 64 を超えると 2 bit パックが
        /// UInt128 に収まらず高速経路から外れるため、
        /// そこで頭を抑える (偶数を避けるので実際に選ばれる最大値は 63)
        /// </summary>
        public const int 自動k長の上限 = 63;

        /// <summary>
        /// k 長の自動選択を行うために最低限必要なリード長<br/>
        /// これより短いリードでは k を十分に取れず、自動選択しても意味がない
        /// </summary>
        public const int 自動k長に必要な最小リード長 = 32;

        /// <summary>
        /// -mk で試す k の個数<br/>
        /// 増やすほど実行時間が線形に伸びる<br/>
        /// ヘルプに実行時間の目安として出るため、ここに置いている
        /// </summary>
        public const int マルチkで試す個数 = 6;

        /// <summary>
        /// 行き止まりの短い unitig を tip (エラー由来) とみなすカバレッジの上限<br/>
        /// グラフ全体のカバレッジ基準値に対する比<br/>
        /// </summary>
        public const double tipとみなすカバレッジ比 = 0.5D;

        /// <summary>
        /// 自動で試す k の上限 (リード長に対する比)<br/>
        /// 単一 k の自動選択より大きく取る<br/>
        /// カバレッジが十分あればリード長に近い k のほうが良いことがあり、
        /// その領域は実際に試さないと分からない
        /// </summary>
        public const double マルチk上限のリード長比 = 0.9D;

        /// <summary>
        /// 自動で試す k の下限
        /// </summary>
        public const int マルチkの下限 = 21;

        /// <summary>
        /// 候補を評価する物差し (アンカー) の k を、候補の最小 k からどれだけ
        /// 下げるか<br/>
        /// 候補と同じ k にすると、その候補だけが自分と同じ条件で
        /// 作った k-mer 集合で測られ、物差しとして中立でなくなる
        /// </summary>
        public const int アンカーk長の候補からの差 = 2;

        /// <summary>
        /// アンカー k の下限<br/>
        /// これ以上小さくすると、細菌のゲノム規模でも
        /// 無関係な座位が同じ k-mer を持ちはじめ、物差しとして荒くなる
        /// </summary>
        public const int アンカーk長の下限 = 11;

        /// <summary>
        /// 試す価値があるとみなす k-mer カバレッジの下限<br/>
        /// k を上げると 1 リードから取れる k-mer が減るため、カバレッジの薄い
        /// データで高い k を試しても時間を捨てるだけになる<br/>
        /// 最初の k の結果から
        /// 各 k のカバレッジを予測し、これを下回るものは実行前に捨てる
        /// </summary>
        public const double マルチkの最小kmerカバレッジ = 10.0D;

        /// <summary>
        /// -kc も、k-mer スペクトルの解析も当てにできなかった場合の最後の拠り所<br/>
        /// 通常は k-mer スペクトルから自動選択される
        /// </summary>
        public const int kmerカットオフの既定値 = 2;

        /// <summary>
        /// k-mer カウント時にメモリ上へ保持するカウントの総量 (バイト)<br/>
        /// 大きいほどフラッシュ回数が減って I/O が軽くなる代わりメモリを使う
        /// </summary>
        public const long メモリ予算の既定値 = 768L * 1024 * 1024;

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
        /// フラグメント長の経験分布を刻むビン幅
        /// </summary>
        public const int フラグメント長のビン幅 = 5;

        /// <summary>
        /// 同一のギャップ長から出たとみなす既知長のばらつき幅の下限<br/>
        /// ライブラリが極端に狭いときに窓が潰れないようにする
        /// </summary>
        public const int 既知長のばらつき幅の下限 = 25;

        /// <summary>
        /// スキャフォールド辺を認めるのに必要な、距離が揃っているペアの本数<br/>
        /// 反復解決が使う -pc とは数える対象が違うので別に持つ
        /// </summary>
        public const ulong スキャフォールド支持数の下限 = 3UL;

        /// <summary>
        /// 短い反復解決の拒否権 (-rv) で使う r-mer 長を、その k でのアセンブリの
        /// k 長にこれだけ足して決める (r = k + この値)<br/>
        /// head と repeat、repeat と tail は de Bruijn グラフの辺である以上、
        /// 必ず k-1 塩基を共有しており、その共有区間は repeat 自身の配列にも
        /// そのまま現れる<br/>
        /// r <= k だと、接合点を跨ぐと判定した窓もこの共有区間の
        /// 内側に収まってしまい、head だけ・repeat だけを読んだリードでも
        /// 真になる (=対応付けの正しさを何も検定できない)<br/>
        /// r を k より
        /// 確実に長く取ることで、共有区間の外側まで踏み込んだリードでなければ
        /// 真になり得ない窓だけを見られるようにする<br/>
        /// ulong に 2 bit パックする都合上 r は 32 が上限で、k + この値が
        /// それを超えるkでは検証自体をスキップする (高い k では反復自体が
        /// 少なく、他の判定で十分間に合っていることが多い)
        /// </summary>
        public const int rMer長のk超過分の既定値 = 10;

        /// <summary>
        /// r-mer のふるいに確保するビット数の上限 (2^31 bit = 256 MB)
        /// </summary>
        public const long rMerふるいのビット数上限 = 1L << 31;

        /// <summary>
        /// r-mer のふるいで 1 件あたりに立てるビットの数
        /// </summary>
        public const int rMerふるいのハッシュ数 = 4;

        /// <summary>
        /// r-mer 検証を行うために、1 本のリードから最低限取れてほしい窓の数
        /// </summary>
        /// <remarks>
        /// r がリード長に近づくと窓が数個しか取れず、r-mer のカバレッジがリードのカバレッジの数 % まで落ちる<br/>
        /// そうなると拒否権が正しい経路まで棄却しはじめるため、痩せすぎる手前で見送る
        /// </remarks>
        public const int rMer検証に必要な窓数 = 20;

        /// <summary>
        /// 短い反復解決の拒否権で、経路の接合点が「実際にリードに読まれている」と
        /// 認めるのに必要な、接合点を跨ぐ r-mer の最小本数
        /// </summary>
        public const int r_mer接合点支持の閾値の既定値 = 4;

        /// <summary>
        /// ユニティグファイル名
        /// </summary>
        public const string ユニティグファイル名 = "unitigs.fasta";

        /// <summary>
        /// コンティグファイル名
        /// </summary>
        public const string コンティグファイル名 = "contigs.fasta";

        /// <summary>
        /// スキャフォールドファイル名
        /// </summary>
        public const string スキャフォールドファイル名 = "scaffolds.fasta";

        /// <summary>
        /// 利用者が受け取る最終成果物
        /// </summary>
        /// <remarks>
        /// unitigs / contigs / scaffolds が採用した k の各段階の出力であるのに対し、
        /// こちらは統合、短い配列の除外、ポリッシュまで通したものになる<br/>
        /// どれを使えばよいかがファイル名だけで分かるようにする
        /// </remarks>
        public const string 最終アセンブリファイル名 = "assembly.fasta";

        /// <summary>
        /// GFA ファイル名
        /// </summary>
        public const string GFAファイル名 = "assembly.gfa";

        /// <summary>
        /// 環状に閉じた複製単位であることを示す、配列 ID 中の目印<br/>
        /// 名前を付ける側 (ContigMaker/Scaffolder) と、それを根拠に数える側
        /// (AssemblyScorer/CircularClosureVerifier/CompletenessValidator) が
        /// 別々に文字列を持つと、片方だけ変えたときに黙って 0 件になる
        /// </summary>
        public const string 環状の目印 = "circular";

        /// <summary>
        /// 環状に閉じた経路を「複製単位が 1 周組み上がった」とみなす最小の長さ<br/>
        /// de Bruijn グラフにはホモポリマーや短いタンデム反復に由来する
        /// 極小の閉路が多数ある (実データで 1 bp〜150 bp の閉路が k あたり
        /// 10 本前後現れた)<br/>
        /// これらを複製単位として数えると、環状本数が
        /// 候補選択の最優先キーである以上、k の選択がその雑音で決まってしまう<br/>
        /// 既知の自然プラスミドで最も小さいものが 1 kb 前後なので、そこで切る<br/>
        /// 下回る閉路も配列としては出力する<br/>
        /// 数えないだけ
        /// </summary>
        public const int 環状として数える最小長 = 1000;

        /// <summary>
        /// レポートファイル名
        /// </summary>
        public const string レポートファイル名 = "assembly.report.json";

        /// <summary>
        /// 曖昧箇所ファイル名
        /// </summary>
        public const string 曖昧箇所ファイル名 = "assembly.ambiguous.tsv";

        /// <summary>
        /// リードに裏付けの無い箇所の一覧
        /// </summary>
        public const string 支持のない箇所ファイル名 = "assembly.unsupported.tsv";

        /// <summary>
        /// リードの支持を問うときの r-mer 長
        /// </summary>
        /// <remarks>
        /// アセンブリの k とは別に短く取る<br/>
        /// k と同じにすると、その k で組んだ配列は定義上すべて支持されてしまい、
        /// 後段 (ギャップ充填、局所アセンブリ、ポリッシュ) が持ち込んだ配列だけを見逃す<br/>
        /// 31 なら偶然の一致がまず起きず、反復の内側でも支持を問える
        /// </remarks>
        public const int 支持検査のr長 = 31;

        /// <summary>
        /// 作業ディレクトリに置く最終成果物のファイル名
        /// </summary>
        /// <remarks>
        /// 中間ファイルの削除で消してはいけないものの一覧でもある
        /// </remarks>
        public static readonly string[] 最終成果物のファイル名 =
        [
            ユニティグファイル名,
            コンティグファイル名,
            スキャフォールドファイル名,
            最終アセンブリファイル名,
            GFAファイル名,
            レポートファイル名,
            曖昧箇所ファイル名,
            支持のない箇所ファイル名,
            ログファイル名,
        ];

        /// <summary>
        /// ポリッシュ結果の一時的な置き場<br/>
        /// 最後に最終成果物へ被せる
        /// </summary>
        public const string ポリッシュ済みファイル名 = "polished.fasta";

        /// <summary>
        /// 実行中に出した内容を全量残すファイル<br/>
        /// 画面をどれだけ静かにしても、
        /// また一時ディレクトリを消す指定があっても、これだけは残す
        /// </summary>
        public const string ログファイル名 = "Tsumiki.log";

        /// <summary>
        /// 塩基の内部表現<br/>
        /// A/C/G/T は塩基記号そのものなので英字のまま残す
        /// (日本語にするとかえって読みにくいため)
        /// </summary>
        public static class 塩基ID
        {
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
        }

        /// <summary>
        /// 進捗ログ間隔
        /// </summary>
        public const ulong 進捗ログ間隔 = 100_000;

        /// <summary>
        /// 無効な塩基
        /// </summary>
        public const byte 無効な塩基 = 5;

        /// <summary>
        /// ユニティグ数の上限
        /// </summary>
        public const int ユニティグ数の上限 = 100_000;

        /// <summary>
        /// インサートサイズが未指定の場合の自動推定に使う、単一 unitig へ両リードが
        /// マップされたペアの最小標本数<br/>
        /// これに満たない場合は推定を諦め、
        /// ペアエンド由来のスキャフォールディングをスキップする
        /// </summary>
        public const int インサートサイズ標本数の下限 = 30;

        /// <summary>
        /// ギャップ長が推定上 0 以下になった場合に最低限挿入する N の数
        /// </summary>
        public const int ギャップ長の下限 = 1;

        /// <summary>
        /// ギャップ充填 (GapFiller・LocalAssembler 共通) で、推定ギャップ長に
        /// 対して許容する誤差 (塩基)<br/>
        /// インサートサイズ推定のばらつきが
        /// そのままギャップ長推定のばらつきになるため、ぴったりの長さだけを
        /// 探すと現実にはまず当たらない
        /// </summary>
        public const int ギャップ充填の長さの余裕幅 = 30;

        /// <summary>
        /// ギャップ充填 (GapFiller・LocalAssembler 共通) の対象とするギャップ長の
        /// 上限<br/>
        /// これより長いギャップは探索空間が広すぎるうえ、推定長の誤差も
        /// 大きく一意に定まる見込みが薄いため対象外とする
        /// </summary>
        public const int ギャップ充填のギャップ長上限 = 500;

        /// <summary>
        /// ペアの重なりで合成リードを作るときに要求する最小の重なり長
        /// </summary>
        /// <remarks>
        /// 前処理の相互訂正より厳しくしているのは、こちらは重ねた結果を 1 本の配列として下流へ渡すため、
        /// 偶然の一致で繋ぐと存在しない接合をグラフへ持ち込んでしまうから
        /// </remarks>
        public const int ペア結合の最小重なり長 = 60;

        /// <summary>
        /// ペアの重なりで合成リードを作るときに許す不一致率
        /// </summary>
        public const double ペア結合の許容不一致率 = 0.05D;
    }
}
