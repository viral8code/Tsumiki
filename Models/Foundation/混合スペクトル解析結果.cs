namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// k-mer 出現回数ヒストグラムを「誤り成分 (裾の重い幾何分布) ＋ 真の k-mer 成分 (単一コピー平均 λ の整数倍に山を持つ、コピー数上限までの負の二項混合) 」の 2 成分混合モデルとして EM 推定した結果
    /// </summary>
    /// <param name="A_単一コピー平均"></param>
    /// <param name="A_カットオフ"></param>
    /// <param name="A_信頼下限"></param>
    /// <param name="A_誤り成分の混合比"></param>
    /// <param name="A_誤り成分の平均"></param>
    /// <param name="A_過分散">単一コピー成分の裾の重さ、大きいほどポアソン分布に近い</param>
    /// <param name="A_反復回数"></param>
    internal record 混合スペクトル解析結果(double A_単一コピー平均, ulong A_カットオフ, ulong A_信頼下限, double A_誤り成分の混合比, double A_誤り成分の平均, double A_過分散, int A_反復回数);
}
