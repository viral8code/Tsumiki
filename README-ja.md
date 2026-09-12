# Tsumiki

**Tsumiki** は、シングルエンドおよびペアエンドのショートリードに対応した C# 製の実験的ゲノムアセンブラです。信頼性の高い k-mer から de Bruijn グラフを構築し、リードの接続情報を用いてコンティグやスキャフォールドの構築、および各種検証レポートの出力を行います。

> [!CAUTION]
> **本プロジェクトは現在開発中（アクティブ開発段階）です**
> 内部で生成される完全性評価指標は、入力リードに対する整合性チェック結果であり、完成したゲノムやリファレンスゲノム比較評価の代用となるものではありません。また、実データにおけるアセンブリ精度や計算リソース要件の検証は途上にあります。

---

## 主な機能

- **k-mer 自動最適化**: リード長や頻度に基づく k-mer 長・カットオフの自動決定（手動指定も可）
- **メモリ効率の高い k-mer カウント**: 指定したメモリバジェット内で動作するディスクバック型カウント
- **グラフ構造の解析・最適化**: de Bruijn グラフの簡略化、Unitig 構築、コピー数を考慮したリピート処理、ペアエンド情報に基づくスキャフォールディング
- **前処理 & エラー訂正**: リードのオーバーラップに基づく前処理や、k-mer スペクトラムを用いたエラー訂正（オプション）
- **マルチ k 処理**: 複数 k-mer でのアセンブリ実行、配列の引き継ぎ（Carry-over）、最良候補の選定（オプション）
- **局所リファインメント & ポリッシング**: 有界コンテキスト長での局所ギャップ再アセンブリ、塩基置換のポリッシング、環状接合部の検証（オプション）
- **多様な出力フォーマット**: FASTA、GFA1（グラフ構造）、JSON（詳細レポート）、TSV（各種エビデンスサマリー）
- **多言語 CLI サポート**: 日本語・英語・中国語によるコマンドラインメッセージ出力

---

## ビルド環境・動作手順

本プロジェクトは **Windows 向け .NET 10** (`net10.0-windows`) をターゲットとしています。.NET 10 SDK を事前にインストールしてください。

### ビルドと実行

リポジトリ直下で以下のコマンドを実行します。

```powershell
# リリースビルドとヘルプの表示（日本語）
dotnet build Tsumiki.csproj -c Release
dotnet run --project Tsumiki.csproj -c Release --no-build -- -h -lang ja
```

### フレームワーク依存ビルドのパブリッシュ (Native AOT 無効)

```powershell
dotnet publish Tsumiki.csproj -c Release -p:PublishAot=false -o publish
dotnet publish/Tsumiki.dll -h -lang ja
```

*※デフォルトの設定では Native AOT が有効化されています。Native AOT パブリッシュには適切な C++ コンパイル環境が必要です。上記の例では `-p:PublishAot=false` で明示的に無効化しています。現在、Linux および macOS での動作検証・ドキュメント化は未対応です。*

---

## クイックスタート

### ペアエンドリードのアセンブリ

```powershell
dotnet publish/Tsumiki.dll -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t assembly-run -lang ja
```

### シングルエンドリードのアセンブリ

```powershell
dotnet publish/Tsumiki.dll -1 reads.fastq.gz -t single-end-run -lang ja
```

### マルチ k・各種オプションの有効化

```powershell
dotnet publish/Tsumiki.dll -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -k 31,51,63 -pp -ec -po -cc -gfa -t multi-k-run -th 8 -mem 2G -lang ja
```

* ※ `-k` にカンマ区切りで複数の数値を渡すと、自動的にマルチ k モードで実行されます。指定する k 値はリード長未満に設定してください。
* ※ 入力ファイルは FASTQ 形式（gzip 圧縮対応）に対応しています。ペアエンドの場合はファイル間でリードの並び順が一致している必要があります。
* ※ 曖昧な塩基（IUPAC コード）や低クオリティな塩基を含むウィンドウは標準で除外されます。`-ab` を指定することで、最大 64 パターンまでの塩基展開を許容できます。

---

## コマンドラインオプション一覧

`-h` オプションで全パラメータの詳細を確認できます。

| オプション | 説明 | デフォルト値 |
|---|---|---|
| `-1 <path>` | フォワード (R1) またはシングルエンド FASTQ パス | **必須** |
| `-2 <path>` | リバース (R2) FASTQ パス | なし |
| `-ab` | IUPAC 曖昧塩基の候補展開を許可 | Off |
| `-k <n[,n...]>` | k-mer 長（カンマ区切りで複数指定可） | リードから自動推定 (上限 63) |
| `-kc <n>` | 信頼できる k-mer の最小カウント閾値 | スペクトラムから自動推定 |
| `-p <33\|64>` | Phred クオリティスコアのオフセット指定 | 自動推定 (デフォルト 33) |
| `-q <n>` | 信頼できる最小塩基クオリティ (Phred Score) | `1` |
| `-mem <size>` | k-mer カウント用のメモリ上限 (例: `512M`, `2G`) | `768M` |
| `-th <n>` | 使用するスレッド数 | CPU の論理プロセッサ数 |
| `-i <n>` | 推定インサートサイズ (bp) | アライメント結果から自動推定 |
| `-pu <ratio>` | ペア接続の採択に必要な優位度比率 | `0.8` |
| `-pc <n>` | ショートリピート解決に必要な最小ペアリード数 | `10` |
| `-mode <conservative\|normal\|bold>` | ペア接続閾値のプリセット | `normal` |
| `-pp` | アダプター除去 & リード間オーバーラップ訂正 | Off |
| `-ec` | k-mer スペクトラムを用いたリードエラー訂正 | Off |
| `-mk` | 複数 k 候補の生成・評価を有効化 | Off |
| `-nc` | k 実行間での配列引き継ぎ（Carry-over）を無効化 | 有効 (Off 指定で無効) |
| `-sr` | オーバーラップペアからの合成配列引き継ぎ | Off |
| `-mg` | 異言語 / 別 k アセンブリ結果の統合・マージ | Off |
| `-rv` | リピート接合部での厳密な r-mer エビデンスを要求 | Off |
| `-la` | 局所ギャップアセンブリ (曖昧領域の長 k 再試行) | Off |
| `-my` | 低頻度 k-mer の救出処理 | Off |
| `-po` | 最終アセンブリの置換ポリッシング (修復) | Off |
| `-cc` | 環状接合部を跨ぐリードの整合性検証 | Off |
| `-gfa` | GFA1 形式での Unitig グラフ出力 | Off |
| `-t <path>` | 出力・一時ファイルの保存先ディレクトリス | `temp` |
| `-rs` | 以前の前処理・訂正結果をキャッシュとして再利用 | Off |
| `-rt` | 正常完了時におよび中間ファイルを自動削除 | Off |
| `-lang <ja\|en\|zh>` | 表示言語 (`ja`: 日本語, `en`: 英語, `zh`: 中国語) | `ja` |
| `-log <quiet\|normal\|verbose>` | ログ出力の詳細度 | `normal` |
| `-v` | バージョン情報の表示 | — |
| `-h` | ヘルプメッセージの表示 | — |

> **メモリ使用量に関する注意**  
> `-mem` は k-mer カウント処理のバジェット領域を指します。グラフ展開やマッピングインデックス構築、マルチ k 処理時には別途メモリが割り当てられます。特に `-rv` オプションなどは、高ノイズリード使用時にメモリ使用量が大幅に増加する可能性があります。

---

## 出力ディレクトリ構造

指定した出力ディレクトリ (`-t`) には以下のファイルが生成されます。

```text
assembly-run/
├── assembly.fasta           # 【最終出力】アセンブリ結果の FASTA ファイル
├── unitigs.fasta            # 構築された Unitig 配列
├── contigs.fasta            # スキャフォールディング前のコンティグ配列
├── scaffolds.fasta          # スキャフォールド配列 (生成時のみ)
├── assembly.report.json     # アセンブリ統計・整合性検証レポート
├── assembly.provenance.json # 実行環境、入力ファイルの SHA-256、設定ログ
├── assembly.ambiguous.tsv   # 曖昧領域のサマリー
├── assembly.unsupported.tsv # リードサポートが不十分な区間一覧
├── assembly.gfa             # GFA1 グラフデータ (-gfa 有効時)
├── Tsumiki.log              # 実行ログ
├── k31/                     # 各 k 値での中間計算結果
└── validation/              # 最終検証用インデックスデータ
```

> **注意**  
> 最終的な解析結果には **`assembly.fasta`** を使用してください。フィルタリングやポリッシング適用済みの確定結果が保持されます。

---

## 中断と再開 (Resume)

既存の出力ディレクトリが存在する場合、安全のため再実行は拒否されます。前処理やエラー訂正の中間成果物を再利用して再実行する場合は、`-rs` オプションを付与してください。

```powershell
dotnet publish/Tsumiki.dll -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -pp -ec -t assembly-run -rs -lang ja
```

※入力ファイルやパラメータ、ビルドハッシュに変更がない場合のみ中間処理が再利用されます。設定を変更した場合は新しい出力ディレクトリを指定することを推奨します。

---

## テストと開発

```powershell
dotnet test Tsumiki.Tests/Tsumiki.Tests.csproj --verbosity minimal -p:GenerateDocumentationFile=true -warnaserror
```

* リポジトリに同梱されている `read.1.fq`, `read.2.fq` はユニットテスト用のダミーファイルです。アセンブリ精度や速度のベンチマークには使用しないでください。
* ローカルでのベンチマークテストやデータ検証を行う場合は、Git 除外設定がされている `local-data/` ディレクトリ配下などを使用してください。

---

## リポジトリ構造

```text
Tsumiki/
├── Program.cs             # エントリーポイント
├── Cores/
│   ├── Pipeline/          # パイプライン制御、マルチk統合、検証処理
│   ├── Preprocessing/     # トリミング、エラー訂正、k-merカウント
│   ├── UnitigBuilding/    # グラフ構築、グラフ単純化、コピー数推定
│   ├── ContigBuilding/    # リードマッピング、コンティグ探索
│   ├── Scaffolding/       # スキャフォールディング、ギャップリファインメント
│   └── Evaluation/        # 候補評価、レポート生成・データ出力
├── IO/                    # 入出力・ファイル読み込み
├── Utilities/             # 各種ユーティリティ・インデックス構造
├── Models/                # データモデル定義
├── Commons/               # 共通設定・メッセージ定義
└── Tsumiki.Tests/         # ユニットテスト・回帰テスト
```

---

## ライセンス

[LICENSE.txt](LICENSE.txt) に従って MIT ライセンスで提供されています。
