// Tsumiki.Common.ConfigurationManager は静的な共有可変状態(Arguments.Kmer 等)を
// 保持しており、多くのテストがそれを一時的に書き換えて動作を検証する。
// xUnit の既定であるテストクラス間の並列実行を許すと、異なるクラスが
// 同時に Arguments を書き換え合ってフレーキーになりうるため、
// このアセンブリ全体でテスト並列化を無効化する。
[assembly: Xunit.CollectionBehavior(DisableTestParallelization = true)]
