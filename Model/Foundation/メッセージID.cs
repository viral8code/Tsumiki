namespace Tsumiki.Model.Foundation
{
    /// <summary>
    /// 画面に出す文言の識別子
    /// </summary>
    /// <remarks>
    /// 文言そのものは <see cref="Common.MessageCatalog"/> が
    /// 言語ごとに持つ<br/>
    /// 呼び出し側はこの ID だけを指す
    /// </remarks>
    internal enum メッセージID
    {
        /// <summary>
        /// de Bruijn グラフの要約情報
        /// </summary>
        デブルーイングラフの要約,

        /// <summary>
        /// 先読みで解決できた分岐の数
        /// </summary>
        先読みで解決した分岐数,

        /// <summary>
        /// GFA ファイルの出力が完了した
        /// </summary>
        GFA出力完了,

        /// <summary>
        /// 分岐選択に使った重みの内訳
        /// </summary>
        分岐選択の重み内訳,

        /// <summary>
        /// グラフ単純化が収束した
        /// </summary>
        単純化の収束,

        /// <summary>
        /// グラフ単純化を打ち切った
        /// </summary>
        単純化の打ち切り,

        /// <summary>
        /// 除去したバブルの数
        /// </summary>
        バブル除去数,

        /// <summary>
        /// 解決した反復の数
        /// </summary>
        反復解決数,

        /// <summary>
        /// 辺選択の内訳
        /// </summary>
        辺選択の内訳,

        /// <summary>
        /// 反復を通り抜けるとして棄却した結合の数
        /// </summary>
        反復通り抜けで棄却した結合数,

        /// <summary>
        /// 相互一意性検査を通って残った結合の数
        /// </summary>
        相互一意で残った結合数,

        /// <summary>
        /// 環状コンティグが見つかった
        /// </summary>
        環状コンティグあり,

        /// <summary>
        /// 環状コンティグが見つからなかった
        /// </summary>
        環状コンティグなし,

        /// <summary>
        /// 確定辺として採った標本の数
        /// </summary>
        確定辺標本数,

        /// <summary>
        /// 確定辺標本の中央値
        /// </summary>
        確定辺標本の中央値,

        /// <summary>
        /// 短すぎる unitig を除外した
        /// </summary>
        短すぎるユニティグの除外,

        /// <summary>
        /// 曖昧塩基を含む k-mer を登録した
        /// </summary>
        曖昧なkmer登録,

        /// <summary>
        /// ペアの隣接候補の数
        /// </summary>
        ペア隣接候補数,

        /// <summary>
        /// 同一 unitig 内のペアの向きの集計
        /// </summary>
        同一ユニティグのペア向き集計,

        /// <summary>
        /// 同一 unitig 内の断片長の分布
        /// </summary>
        同一ユニティグの断片長分布,

        /// <summary>
        /// 同一 unitig 内の断片長の中央値
        /// </summary>
        同一ユニティグの断片長中央値,

        /// <summary>
        /// ペアのリード ID が一致しない
        /// </summary>
        ペアリードIDの不一致,

        /// <summary>
        /// 統合できる接合点が無かった
        /// </summary>
        統合できる接合点なし,

        /// <summary>
        /// 統合した接合点の数
        /// </summary>
        統合した接合点数,

        /// <summary>
        /// 候補一覧の見出し行
        /// </summary>
        候補一覧の見出し,

        /// <summary>
        /// 候補一覧の明細行
        /// </summary>
        候補一覧の明細,

        /// <summary>
        /// 統計の対象となるファイルが無い
        /// </summary>
        統計_ファイルなし,

        /// <summary>
        /// 統計情報
        /// </summary>
        統計,

        /// <summary>
        /// 長さで絞り込んだ統計情報
        /// </summary>
        統計_長さで絞り込み,

        /// <summary>
        /// 検査対象外の k 長
        /// </summary>
        検査_対象外のk長,

        /// <summary>
        /// 取りこぼしの検査結果
        /// </summary>
        検査_取りこぼし,

        /// <summary>
        /// 出しすぎの検査結果
        /// </summary>
        検査_出しすぎ,

        /// <summary>
        /// k-mer インデックスの構築を開始した
        /// </summary>
        kmerインデックス構築開始,

        /// <summary>
        /// k-mer カットオフを適用した
        /// </summary>
        kmerカットオフ適用,

        /// <summary>
        /// 引き継ぎで追加した k-mer の数
        /// </summary>
        引き継ぎで追加したkmer数,

        /// <summary>
        /// tip 除去を開始した
        /// </summary>
        tip除去開始,

        /// <summary>
        /// unitig の構築を開始した
        /// </summary>
        ユニティグ構築開始,

        /// <summary>
        /// グラフが複雑すぎる
        /// </summary>
        グラフが複雑すぎる,

        /// <summary>
        /// リードのマッピングを開始した
        /// </summary>
        リードのマッピング開始,

        /// <summary>
        /// unitig の結合を開始した
        /// </summary>
        ユニティグ結合開始,

        /// <summary>
        /// r-mer 検証を見送った
        /// </summary>
        rMer検証の見送り,

        /// <summary>
        /// コンティグの構築が完了した
        /// </summary>
        コンティグ構築完了,

        /// <summary>
        /// スキャフォールディングを開始した
        /// </summary>
        スキャフォールディング開始,

        /// <summary>
        /// ギャップ充填を開始した
        /// </summary>
        ギャップ充填開始,

        /// <summary>
        /// k のカバレッジが薄すぎて省略した
        /// </summary>
        kが薄すぎて省略,

        /// <summary>
        /// k の処理開始を示す見出し
        /// </summary>
        kの開始見出し,

        /// <summary>
        /// その k ではアセンブリできなかった
        /// </summary>
        kでアセンブリできず,

        /// <summary>
        /// 単一の k だけが成功した
        /// </summary>
        単一のkのみ成功,

        /// <summary>
        /// アンカーとする k-mer 集合を構築した
        /// </summary>
        アンカーkmer集合の構築,

        /// <summary>
        /// アンカーの k-mer スペクトルが二峰性でない
        /// </summary>
        アンカースペクトルが二峰でない,

        /// <summary>
        /// 候補を評価できない
        /// </summary>
        候補を評価できない,

        /// <summary>
        /// 採用した k の値
        /// </summary>
        採用したk,

        /// <summary>
        /// 統合処理を開始した
        /// </summary>
        統合開始,

        /// <summary>
        /// 統合結果を評価できない
        /// </summary>
        統合結果を評価できない,

        /// <summary>
        /// 統合前の評価結果
        /// </summary>
        統合前の評価,

        /// <summary>
        /// 統合後の評価結果
        /// </summary>
        統合後の評価,

        /// <summary>
        /// 統合結果が骨格アセンブリに勝てなかった
        /// </summary>
        統合が骨格に勝てず,

        /// <summary>
        /// 統合結果を採用した
        /// </summary>
        統合結果を採用,

        /// <summary>
        /// エラー訂正用の k-mer スペクトルを構築した
        /// </summary>
        エラー訂正_スペクトル構築,

        /// <summary>
        /// エラー訂正を開始した
        /// </summary>
        エラー訂正_訂正開始,

        /// <summary>
        /// ファイルごとのエラー訂正統計
        /// </summary>
        エラー訂正_ファイル別統計,

        /// <summary>
        /// リードの読み込みが完了した
        /// </summary>
        リード読込完了,

        /// <summary>
        /// リード 2 の読み込みを開始した
        /// </summary>
        リード2の読込開始,

        /// <summary>
        /// 前処理の集計結果
        /// </summary>
        前処理統計,

        /// <summary>
        /// SuperRead の集計結果
        /// </summary>
        SuperRead統計,

        /// <summary>
        /// SuperRead 生成時の重なりに関する統計
        /// </summary>
        SuperRead統計_重なり,

        /// <summary>
        /// SuperRead 生成時の曖昧箇所に関する統計
        /// </summary>
        SuperRead統計_曖昧,

        /// <summary>
        /// ギャップ充填の対象が無い
        /// </summary>
        ギャップ充填_対象なし,

        /// <summary>
        /// ギャップ充填の集計結果
        /// </summary>
        ギャップ充填統計,

        /// <summary>
        /// 局所アセンブリの対象が無い
        /// </summary>
        局所アセンブリ_対象なし,

        /// <summary>
        /// 局所アセンブリの集計結果
        /// </summary>
        局所アセンブリ統計,

        /// <summary>
        /// インサートサイズが不明なためスキャフォールディングを省略した
        /// </summary>
        スキャフォールディング省略_インサートサイズ不明,

        /// <summary>
        /// インサートサイズを添えてスキャフォールディングを開始した
        /// </summary>
        スキャフォールディング開始_インサートサイズ,

        /// <summary>
        /// コンティグが無いためスキャフォールディングを省略した
        /// </summary>
        スキャフォールディング省略_コンティグなし,

        /// <summary>
        /// 配列内部を指したペア候補の数
        /// </summary>
        内部を指したペア候補,

        /// <summary>
        /// 未配置の配列を指したペア候補の数
        /// </summary>
        未配置を指したペア候補,

        /// <summary>
        /// 閾値適用後に残ったスキャフォールド辺
        /// </summary>
        閾値後のスキャフォールド辺,

        /// <summary>
        /// スキャフォールドの出力が完了した
        /// </summary>
        スキャフォールド出力完了,

        /// <summary>
        /// 同一 unitig 内のペアからのインサートサイズ推定
        /// </summary>
        インサートサイズ推定_同一ユニティグ,

        /// <summary>
        /// 確定辺からのインサートサイズ推定
        /// </summary>
        インサートサイズ推定_確定辺,

        /// <summary>
        /// インサートサイズ推定に必要な標本が不足している
        /// </summary>
        インサートサイズ推定_標本不足,

        /// <summary>
        /// インサートサイズ推定に使った全標本数
        /// </summary>
        インサートサイズ推定_全標本,

        /// <summary>
        /// 単一コピーとみなすカバレッジの基準値
        /// </summary>
        単一コピー基準値,

        /// <summary>
        /// コピー数推定の要約
        /// </summary>
        コピー数の要約,

        /// <summary>
        /// 反復配列とみなされた配列の割合
        /// </summary>
        反復配列の割合,

        /// <summary>
        /// グラフ単純化の反復回数
        /// </summary>
        グラフ単純化の反復,

        /// <summary>
        /// r-mer 検証によって棄却された件数
        /// </summary>
        rMer検証による棄却,

        /// <summary>
        /// 無視した例外の見出し
        /// </summary>
        例外を無視_見出し,

        /// <summary>
        /// 例外を無視したメソッド名
        /// </summary>
        例外を無視_メソッド,

        /// <summary>
        /// 処理停止の見出し
        /// </summary>
        停止_見出し,

        /// <summary>
        /// 処理を停止したメソッド名
        /// </summary>
        停止_メソッド,

        /// <summary>
        /// 出力時刻のタイムスタンプ
        /// </summary>
        タイムスタンプ,

        /// <summary>
        /// 明示指定した Phred オフセットが推定値と一致しない
        /// </summary>
        Phred_明示指定と不一致,

        /// <summary>
        /// Phred オフセットを自動判定した
        /// </summary>
        Phred_自動判定,

        /// <summary>
        /// Phred オフセット検査の警告
        /// </summary>
        Phred_検査の警告,

        /// <summary>
        /// 混合モデルから求めた k-mer カットオフ
        /// </summary>
        kmerカットオフ_混合モデル,

        /// <summary>
        /// k-mer スペクトルの谷が不明でカットオフを決められない
        /// </summary>
        kmerカットオフ_谷が不明,

        /// <summary>
        /// スペクトルの谷から求めた k-mer カットオフ
        /// </summary>
        kmerカットオフ_スペクトル,

        /// <summary>
        /// k-mer 出現回数のヒストグラム
        /// </summary>
        kmerヒストグラム,

        /// <summary>
        /// k-mer スペクトルの谷が不明
        /// </summary>
        スペクトルの谷が不明,

        /// <summary>
        /// k-mer スペクトルの谷とピークの位置
        /// </summary>
        スペクトルの谷とピーク,

        /// <summary>
        /// ゲノムサイズとカバレッジの推定結果
        /// </summary>
        ゲノムサイズとカバレッジの推定,

        /// <summary>
        /// リード長が不明で k を自動選択できない
        /// </summary>
        k自動選択_リード長不明,

        /// <summary>
        /// 自動選択した k がリード長に対して長すぎる
        /// </summary>
        k自動選択_kが長すぎる,

        /// <summary>
        /// リード長が短く k の自動選択に不向き
        /// </summary>
        k自動選択_リード長が短い,

        /// <summary>
        /// k を自動選択した結果
        /// </summary>
        k自動選択,

        /// <summary>
        /// 開始点となる k-mer の探索
        /// </summary>
        開始kmerの探索,

        /// <summary>
        /// 一時ディレクトリを削除せず残した
        /// </summary>
        一時ディレクトリを残した,

        /// <summary>
        /// 観測されたリード長
        /// </summary>
        リード長の観測値,

        /// <summary>
        /// 一時ディレクトリが既に存在する
        /// </summary>
        一時ディレクトリが既にある,

        /// <summary>
        /// 指定されたパスの確認結果
        /// </summary>
        パスの確認,

        /// <summary>
        /// ペアが無いため前処理を省略した
        /// </summary>
        前処理省略_ペアなし,

        /// <summary>
        /// 前処理を開始した
        /// </summary>
        前処理開始,

        /// <summary>
        /// エラー訂正を開始した
        /// </summary>
        エラー訂正開始,

        /// <summary>
        /// 開発中で未対応の機能
        /// </summary>
        開発中,

        /// <summary>
        /// コンティグの総延長
        /// </summary>
        コンティグ総延長,

        /// <summary>
        /// リードファイルのパス
        /// </summary>
        リードファイルのパス,

        /// <summary>
        /// 試す k の一覧
        /// </summary>
        試すk一覧,

        /// <summary>
        /// リード読み込みの進捗
        /// </summary>
        リード読込の進捗,

        /// <summary>
        /// リード 1 の読み込みを開始した
        /// </summary>
        リード1の読込開始,

        /// <summary>
        /// シングルエンドリードの読み込みを開始した
        /// </summary>
        単一リードの読込開始,

        /// <summary>
        /// スキャフォールド候補辺の数
        /// </summary>
        スキャフォールド候補辺数,

        /// <summary>
        /// 理想本数モデルが利用できた
        /// </summary>
        理想本数モデルあり,

        /// <summary>
        /// 理想本数モデルが利用できなかった
        /// </summary>
        理想本数モデルなし,

        /// <summary>
        /// リード 1 と 2 で推定 Phred オフセットが一致しない
        /// </summary>
        Phred_ファイル間で不一致,

        /// <summary>
        /// Phred オフセットを確定できなかった
        /// </summary>
        Phred_未確定,

        /// <summary>
        /// 出現した k-mer の種類数
        /// </summary>
        kmer種類数,

        /// <summary>
        /// 採用した k-mer の数
        /// </summary>
        採用kmer数,

        /// <summary>
        /// アセンブリを行えなかった
        /// </summary>
        アセンブリ不能,

        /// <summary>
        /// ツールの概要説明
        /// </summary>
        概要_説明,

        /// <summary>
        /// 作者情報
        /// </summary>
        概要_作者,

        /// <summary>
        /// バージョン情報
        /// </summary>
        概要_バージョン,

        /// <summary>
        /// 使い方の見出し
        /// </summary>
        ヘルプ_使い方,

        /// <summary>
        /// 入力に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_入力,

        /// <summary>
        /// k-mer と品質に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_kmerと品質,

        /// <summary>
        /// ペアエンドに関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_ペアエンド,

        /// <summary>
        /// 前処理に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_前処理,

        /// <summary>
        /// マルチ k に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_マルチk,

        /// <summary>
        /// 反復配列の安全策に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_反復の安全策,

        /// <summary>
        /// 出力とその他に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_出力とその他,

        /// <summary>
        /// -1 オプションのヘルプ文
        /// </summary>
        ヘルプ_リード1,

        /// <summary>
        /// -2 オプションのヘルプ文
        /// </summary>
        ヘルプ_リード2,

        /// <summary>
        /// 曖昧塩基オプションのヘルプ文
        /// </summary>
        ヘルプ_曖昧塩基,

        /// <summary>
        /// -k オプションのヘルプ文
        /// </summary>
        ヘルプ_k長,

        /// <summary>
        /// -kc オプションのヘルプ文
        /// </summary>
        ヘルプ_kmerカットオフ,

        /// <summary>
        /// -p オプションのヘルプ文
        /// </summary>
        ヘルプ_Phred,

        /// <summary>
        /// クオリティカットオフオプションのヘルプ文
        /// </summary>
        ヘルプ_クオリティカットオフ,

        /// <summary>
        /// メモリ予算オプションのヘルプ文
        /// </summary>
        ヘルプ_メモリ予算,

        /// <summary>
        /// インサートサイズオプションのヘルプ文
        /// </summary>
        ヘルプ_インサートサイズ,

        /// <summary>
        /// ペア結合閾値オプションのヘルプ文
        /// </summary>
        ヘルプ_ペア結合閾値,

        /// <summary>
        /// ペア支持数閾値オプションのヘルプ文
        /// </summary>
        ヘルプ_ペア支持数閾値,

        /// <summary>
        /// 積極性モードオプションのヘルプ文
        /// </summary>
        ヘルプ_積極性モード,

        /// <summary>
        /// 前処理オプションのヘルプ文
        /// </summary>
        ヘルプ_前処理,

        /// <summary>
        /// エラー訂正オプションのヘルプ文
        /// </summary>
        ヘルプ_エラー訂正,

        /// <summary>
        /// マルチ k オプションのヘルプ文
        /// </summary>
        ヘルプ_マルチk,

        /// <summary>
        /// 引き継ぎなしオプションのヘルプ文
        /// </summary>
        ヘルプ_引き継ぎなし,

        /// <summary>
        /// SuperRead オプションのヘルプ文
        /// </summary>
        ヘルプ_SuperRead,

        /// <summary>
        /// マージオプションのヘルプ文
        /// </summary>
        ヘルプ_マージ,

        /// <summary>
        /// 反復 r-mer 検証オプションのヘルプ文
        /// </summary>
        ヘルプ_反復rMer検証,

        /// <summary>
        /// 局所アセンブリオプションのヘルプ文
        /// </summary>
        ヘルプ_局所アセンブリ,

        /// <summary>
        /// GFA 出力オプションのヘルプ文
        /// </summary>
        ヘルプ_GFA出力,

        /// <summary>
        /// 一時ディレクトリオプションのヘルプ文
        /// </summary>
        ヘルプ_一時ディレクトリ,

        /// <summary>
        /// 一時ディレクトリ削除オプションのヘルプ文
        /// </summary>
        ヘルプ_一時ディレクトリ削除,

        /// <summary>
        /// スレッド数オプションのヘルプ文
        /// </summary>
        ヘルプ_スレッド数,

        /// <summary>
        /// 言語オプションのヘルプ文
        /// </summary>
        ヘルプ_言語,

        /// <summary>
        /// バージョンオプションのヘルプ文
        /// </summary>
        ヘルプ_バージョン,

        /// <summary>
        /// ヘルプオプション自体のヘルプ文
        /// </summary>
        ヘルプ_ヘルプ,

        /// <summary>
        /// 前段からの引き継ぎ統合を開始した
        /// </summary>
        引き継ぎの統合開始,

        /// <summary>
        /// 引き継ぎ統合の進捗
        /// </summary>
        引き継ぎの統合進捗,

        /// <summary>
        /// 引き継ぎの準備を開始した
        /// </summary>
        引き継ぎの準備開始,

        /// <summary>
        /// 引き継ぎの準備が完了した
        /// </summary>
        引き継ぎの準備完了,

        /// <summary>
        /// SuperRead の橋渡し処理を開始した
        /// </summary>
        SuperRead橋渡し開始,

        /// <summary>
        /// SuperRead の橋渡し処理の進捗
        /// </summary>
        SuperRead橋渡し進捗,

        /// <summary>
        /// 短い配列を除外した
        /// </summary>
        短い配列を除外,

        /// <summary>
        /// 最終成果物の出力先
        /// </summary>
        最終成果物,

        /// <summary>
        /// 中間ファイルを削除した
        /// </summary>
        中間ファイルを削除した,

        /// <summary>
        /// リードによる支持検査を開始した
        /// </summary>
        支持検査の開始,

        /// <summary>
        /// リードによる支持検査の結果
        /// </summary>
        支持検査の結果,

        /// <summary>
        /// 支持のない箇所を書き出した
        /// </summary>
        支持のない箇所を書き出した,

        /// <summary>
        /// リードの支持に関する検査項目
        /// </summary>
        検査項目_リードの支持,

        /// <summary>
        /// ポリッシュ処理を開始した
        /// </summary>
        ポリッシュ開始,

        /// <summary>
        /// ポリッシュ用の索引を構築した
        /// </summary>
        ポリッシュの索引構築,

        /// <summary>
        /// ポリッシュの種となる配列が無い
        /// </summary>
        ポリッシュの種が無い,

        /// <summary>
        /// ポリッシュのためのリードマッピングを開始した
        /// </summary>
        ポリッシュのマッピング開始,

        /// <summary>
        /// ポリッシュで求めた深度の中央値
        /// </summary>
        ポリッシュの深度中央値,

        /// <summary>
        /// ポリッシュを行えなかった
        /// </summary>
        ポリッシュを行えず,

        /// <summary>
        /// ポリッシュのマッピング結果
        /// </summary>
        ポリッシュのマッピング結果,

        /// <summary>
        /// ポリッシュによる訂正結果
        /// </summary>
        ポリッシュの訂正結果,

        /// <summary>
        /// ポリッシュで見つかった深度不足の箇所
        /// </summary>
        ポリッシュの深度不足,

        /// <summary>
        /// 環状の閉じ目の検証を開始した
        /// </summary>
        閉じ目の検証開始,

        /// <summary>
        /// 検証対象となる環状の配列が無い
        /// </summary>
        閉じ目_環状の配列が無い,

        /// <summary>
        /// 閉じ目をリードで裏付けられた
        /// </summary>
        閉じ目を裏付けた,

        /// <summary>
        /// 閉じ目をリードで裏付けられなかった
        /// </summary>
        閉じ目を裏付けられず,

        /// <summary>
        /// グラフ被覆に関する検査項目
        /// </summary>
        検査項目_グラフ被覆,

        /// <summary>
        /// コピー数整合に関する検査項目
        /// </summary>
        検査項目_コピー数整合,

        /// <summary>
        /// 深度の連続性に関する検査項目
        /// </summary>
        検査項目_深度の連続性,

        /// <summary>
        /// 未解決のギャップに関する検査項目
        /// </summary>
        検査項目_未解決のギャップ,

        /// <summary>
        /// 接合点の支持に関する検査項目
        /// </summary>
        検査項目_接合点の支持,

        /// <summary>
        /// 競合経路に関する検査項目
        /// </summary>
        検査項目_競合経路,

        /// <summary>
        /// 環状閉鎖に関する検査項目
        /// </summary>
        検査項目_環状閉鎖,

        /// <summary>
        /// 検査判定が合格であることを示す文言
        /// </summary>
        検査判定_合格,

        /// <summary>
        /// 検査判定が不合格であることを示す文言
        /// </summary>
        検査判定_不合格,

        /// <summary>
        /// 検査判定が判定不能であることを示す文言
        /// </summary>
        検査判定_判定不能,

        /// <summary>
        /// 完全性判定結果の見出し
        /// </summary>
        完全性の見出し,

        /// <summary>
        /// 完全性判定の検査結果を示す行
        /// </summary>
        完全性の検査行,

        /// <summary>
        /// 完全長と判定されたことを示す文言
        /// </summary>
        完全長と判定,

        /// <summary>
        /// 完全長に届かなかったことを示す文言
        /// </summary>
        完全長に届かず,

        /// <summary>
        /// 完全長に届かなかった理由の一覧
        /// </summary>
        完全性の未達理由,

        /// <summary>
        /// 救済 k-mer の処理を開始した
        /// </summary>
        救済kmerの開始,

        /// <summary>
        /// 救済処理の対象外となる k 長
        /// </summary>
        救済_対象外のk長,

        /// <summary>
        /// 救済した k-mer の数
        /// </summary>
        救済したkmer数,

        /// <summary>
        /// レポートファイルを書き出した
        /// </summary>
        レポートを書き出した,

        /// <summary>
        /// 曖昧箇所の記録を書き出した
        /// </summary>
        曖昧箇所を書き出した,

        /// <summary>
        /// 再開時に中間ファイルを再利用した
        /// </summary>
        再開_中間ファイルを再利用,

        /// <summary>
        /// 完全性の検証に関するヘルプ節の見出し
        /// </summary>
        ヘルプ節_完全性の検証,

        /// <summary>
        /// ポリッシュオプションのヘルプ文
        /// </summary>
        ヘルプ_ポリッシュ,

        /// <summary>
        /// 環状閉鎖検証オプションのヘルプ文
        /// </summary>
        ヘルプ_環状閉鎖検証,

        /// <summary>
        /// 救済 k-mer オプションのヘルプ文
        /// </summary>
        ヘルプ_救済kmer,

        /// <summary>
        /// 再開オプションのヘルプ文
        /// </summary>
        ヘルプ_再開,

        /// <summary>
        /// レポート出力オプションのヘルプ文
        /// </summary>
        ヘルプ_レポートの説明,

        /// <summary>
        /// 分岐の無い閉路を検出したことを示す文言
        /// </summary>
        分岐のない閉路,

        /// <summary>
        /// 孤立した環状の複製単位を検出したことを示す文言
        /// </summary>
        孤立した環状の複製単位,

        /// <summary>
        /// ログ水準オプションのヘルプ文
        /// </summary>
        ヘルプ_ログ水準,

        /// <summary>
        /// ログファイルの保存先
        /// </summary>
        ログの保存先,

        /// <summary>
        /// 短すぎる閉路を検出したことを示す文言
        /// </summary>
        短すぎる閉路,

        /// <summary>
        /// 合成リード (SuperRead) を再利用したことを示す文言
        /// </summary>
        合成リードを再利用,
    }
}
