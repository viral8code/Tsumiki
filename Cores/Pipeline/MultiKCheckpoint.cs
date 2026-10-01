using System.Text;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// multi-k の各 k の結果と次の k への引き継ぎを保存し、再開時に読み戻す
    /// </summary>
    internal static class MultiKCheckpoint
    {
        #region 定数

        /// <summary>
        /// k ごとの記録のファイル名
        /// </summary>
        private const string C_k記録ファイル名 = "resume.bin";

        /// <summary>
        /// 合成リードの記録のファイル名
        /// </summary>
        private const string C_合成リード記録ファイル名 = "superreads.bin";

        /// <summary>
        /// 書き込み途中の記録の拡張子
        /// </summary>
        private const string C_一時拡張子 = ".tmp";

        /// <summary>
        /// 記録の形式の版
        /// </summary>
        private const int C_形式の版 = 1;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// その k の記録を消す
        /// </summary>
        /// <param name="p_k作業ディレクトリ"></param>
        public static void V_削除_k(string p_k作業ディレクトリ)
        {
            var l_パス = Path.Combine(p_k作業ディレクトリ, C_k記録ファイル名);
            if (File.Exists(l_パス))
            {
                File.Delete(l_パス);
            }
        }

        /// <summary>
        /// その k の結果と、合成リードを除いた次の k への引き継ぎを保存する
        /// </summary>
        /// <param name="p_k作業ディレクトリ"></param>
        /// <param name="p_署名">入力・設定・k の候補一覧から作った識別値</param>
        /// <param name="p_結果"></param>
        /// <param name="p_引き継ぎ">合成リードを除いた引き継ぎ</param>
        /// <param name="p_合成リード数">続けて渡す合成リードの本数</param>
        public static void V_保存_k(string p_k作業ディレクトリ, string p_署名, アセンブリ実行結果 p_結果, IReadOnlyList<引き継ぎ配列> p_引き継ぎ, int p_合成リード数)
        {
            V_書き込み(Path.Combine(p_k作業ディレクトリ, C_k記録ファイル名), p_署名, l_書き込み =>
            {
                V_書き込み_結果(l_書き込み, p_結果);
                l_書き込み.Write(p_合成リード数);
                V_書き込み_引き継ぎ群(l_書き込み, p_引き継ぎ);
            });
        }

        /// <summary>
        /// その k の記録を読み戻す
        /// </summary>
        /// <param name="p_k作業ディレクトリ"></param>
        /// <param name="p_署名"></param>
        /// <param name="p_結果"></param>
        /// <param name="p_引き継ぎ">合成リードを除いた引き継ぎ</param>
        /// <param name="p_合成リード数"></param>
        /// <returns>署名が一致し、結果のファイルが揃っていれば true</returns>
        public static bool Try読込_k(string p_k作業ディレクトリ, string p_署名, out アセンブリ実行結果? p_結果, out List<引き継ぎ配列> p_引き継ぎ, out int p_合成リード数)
        {
            p_結果 = null;
            p_引き継ぎ = [];
            p_合成リード数 = 0;
            アセンブリ実行結果? l_結果 = null;
            List<引き継ぎ配列> l_引き継ぎ = [];
            var l_合成リード数 = 0;
            if (!Try読み込み(Path.Combine(p_k作業ディレクトリ, C_k記録ファイル名), p_署名, l_読み込み =>
            {
                l_結果 = Get_結果(l_読み込み);
                l_合成リード数 = l_読み込み.ReadInt32();
                l_引き継ぎ = Get_引き継ぎ群(l_読み込み);
            }))
            {
                return false;
            }

            if (l_結果 is null || !File.Exists(l_結果.A_unitigパス) || !File.Exists(l_結果.A_contigパス) || (l_結果.A_scaffoldパス is { } l_scaffold && !File.Exists(l_scaffold)))
            {
                return false;
            }

            p_結果 = l_結果;
            p_引き継ぎ = l_引き継ぎ;
            p_合成リード数 = l_合成リード数;
            return true;
        }

        /// <summary>
        /// 最初の k で作った合成リードを保存する
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_署名"></param>
        /// <param name="p_合成リード"></param>
        public static void V_保存_合成リード(string p_一時ディレクトリ, string p_署名, IReadOnlyList<引き継ぎ配列> p_合成リード)
        {
            V_書き込み(Path.Combine(p_一時ディレクトリ, C_合成リード記録ファイル名), p_署名, l_書き込み => V_書き込み_引き継ぎ群(l_書き込み, p_合成リード));
        }

        /// <summary>
        /// 保存した合成リードを読み戻す
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_署名"></param>
        /// <param name="p_本数">k の記録にある本数</param>
        /// <param name="p_合成リード"></param>
        /// <returns>署名と本数が一致すれば true</returns>
        public static bool Try読込_合成リード(string p_一時ディレクトリ, string p_署名, int p_本数, out List<引き継ぎ配列> p_合成リード)
        {
            List<引き継ぎ配列> l_合成リード = [];
            var l_読めたか = Try読み込み(Path.Combine(p_一時ディレクトリ, C_合成リード記録ファイル名), p_署名, l_読み込み => l_合成リード = Get_引き継ぎ群(l_読み込み));
            p_合成リード = l_読めたか && l_合成リード.Count == p_本数 ? l_合成リード : [];
            return l_読めたか && l_合成リード.Count == p_本数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 版と署名に続けて本文を書き、書き終えてから記録の名前に移す
        /// </summary>
        /// <param name="p_パス"></param>
        /// <param name="p_署名"></param>
        /// <param name="p_本文"></param>
        private static void V_書き込み(string p_パス, string p_署名, Action<BinaryWriter> p_本文)
        {
            var l_一時パス = p_パス + C_一時拡張子;
            using (var l_書き込み = new BinaryWriter(new BufferedStream(File.Create(l_一時パス), 1 << 20), Encoding.UTF8))
            {
                l_書き込み.Write(C_形式の版);
                l_書き込み.Write(p_署名);
                p_本文(l_書き込み);
            }

            File.Move(l_一時パス, p_パス, overwrite: true);
        }

        /// <summary>
        /// 版と署名が一致するときだけ本文を読む
        /// </summary>
        /// <param name="p_パス"></param>
        /// <param name="p_署名"></param>
        /// <param name="p_本文"></param>
        /// <returns>読めたら true</returns>
        private static bool Try読み込み(string p_パス, string p_署名, Action<BinaryReader> p_本文)
        {
            if (!File.Exists(p_パス))
            {
                return false;
            }

            try
            {
                using var l_読み込み = new BinaryReader(new BufferedStream(File.OpenRead(p_パス), 1 << 20), Encoding.UTF8);
                if (l_読み込み.ReadInt32() != C_形式の版 || l_読み込み.ReadString() != p_署名)
                {
                    return false;
                }

                p_本文(l_読み込み);
                return true;
            }
            catch (EndOfStreamException)
            {
                return false;
            }
        }

        /// <summary>
        /// k の結果を書く
        /// </summary>
        /// <param name="p_書き込み"></param>
        /// <param name="p_結果"></param>
        private static void V_書き込み_結果(BinaryWriter p_書き込み, アセンブリ実行結果 p_結果)
        {
            p_書き込み.Write(p_結果.A_k長);
            p_書き込み.Write(p_結果.A_unitigパス);
            p_書き込み.Write(p_結果.A_contigパス);
            V_書き込み_省略可(p_書き込み, p_結果.A_scaffoldパス);
            p_書き込み.Write(p_結果.A_kmerカットオフ);
            p_書き込み.Write(p_結果.A_単一コピー基準値);
            V_書き込み_省略可(p_書き込み, p_結果.A_GFAパス);
            p_書き込み.Write(p_結果.A_整合性検査 is not null);
            if (p_結果.A_整合性検査 is { } l_検査)
            {
                p_書き込み.Write(l_検査.A_信頼kmer数);
                p_書き込み.Write(l_検査.A_アセンブリ内の延べ数);
                p_書き込み.Write(l_検査.A_アセンブリ内の種類数);
                p_書き込み.Write(l_検査.A_取りこぼし数);
                p_書き込み.Write(l_検査.A_出しすぎkmer種類数);
                p_書き込み.Write(l_検査.A_余分な延べ数);
            }

            p_書き込み.Write((int)p_結果.A_実際のコピー数基準);
        }

        /// <summary>
        /// k の結果を読む
        /// </summary>
        /// <param name="p_読み込み"></param>
        /// <returns></returns>
        private static アセンブリ実行結果 Get_結果(BinaryReader p_読み込み)
        {
            var l_k長 = p_読み込み.ReadInt32();
            var l_unitigパス = p_読み込み.ReadString();
            var l_contigパス = p_読み込み.ReadString();
            var l_scaffoldパス = Get_省略可(p_読み込み);
            var l_カットオフ = p_読み込み.ReadUInt64();
            var l_基準値 = p_読み込み.ReadDouble();
            var l_GFAパス = Get_省略可(p_読み込み);
            整合性検査結果? l_検査 = p_読み込み.ReadBoolean()
                ? new 整合性検査結果(p_読み込み.ReadInt64(), p_読み込み.ReadInt64(), p_読み込み.ReadInt64(), p_読み込み.ReadInt64(), p_読み込み.ReadInt64(), p_読み込み.ReadInt64())
                : null;
            var l_コピー数基準 = (コピー数基準の出所)p_読み込み.ReadInt32();
            return new アセンブリ実行結果(l_k長, l_unitigパス, l_contigパス, l_scaffoldパス, l_カットオフ, l_基準値, l_GFAパス, l_検査, l_コピー数基準);
        }

        /// <summary>
        /// 引き継ぎ配列の一覧を書く
        /// </summary>
        /// <param name="p_書き込み"></param>
        /// <param name="p_引き継ぎ群"></param>
        private static void V_書き込み_引き継ぎ群(BinaryWriter p_書き込み, IReadOnlyList<引き継ぎ配列> p_引き継ぎ群)
        {
            p_書き込み.Write(p_引き継ぎ群.Count);
            foreach (var l_配列 in p_引き継ぎ群)
            {
                p_書き込み.Write(l_配列.A_配列);
                p_書き込み.Write7BitEncodedInt(l_配列.A_k長);
                p_書き込み.Write(l_配列.A_Is確定経路);
                p_書き込み.Write7BitEncodedInt(l_配列.A_カバレッジ.Length);
                foreach (var l_値 in l_配列.A_カバレッジ)
                {
                    p_書き込み.Write7BitEncodedInt(l_値);
                }

                p_書き込み.Write(l_配列.A_分岐の継ぎ目位置 is not null);
                if (l_配列.A_分岐の継ぎ目位置 is { } l_継ぎ目)
                {
                    p_書き込み.Write7BitEncodedInt(l_継ぎ目.Count);
                    foreach (var l_位置 in l_継ぎ目)
                    {
                        p_書き込み.Write7BitEncodedInt(l_位置);
                    }
                }

                p_書き込み.Write(l_配列.A_未観測の連続範囲 is not null);
                if (l_配列.A_未観測の連続範囲 is { } l_範囲群)
                {
                    p_書き込み.Write7BitEncodedInt(l_範囲群.Count);
                    foreach (var (A_開始, A_終了) in l_範囲群)
                    {
                        p_書き込み.Write7BitEncodedInt(A_開始);
                        p_書き込み.Write7BitEncodedInt(A_終了);
                    }
                }
            }
        }

        /// <summary>
        /// 引き継ぎ配列の一覧を読む
        /// </summary>
        /// <param name="p_読み込み"></param>
        /// <returns></returns>
        private static List<引き継ぎ配列> Get_引き継ぎ群(BinaryReader p_読み込み)
        {
            var l_本数 = p_読み込み.ReadInt32();
            var l_結果 = new List<引き継ぎ配列>(l_本数);
            for (var i = 0; i < l_本数; i++)
            {
                var l_配列 = p_読み込み.ReadString();
                var l_k長 = p_読み込み.Read7BitEncodedInt();
                var l_Is確定経路 = p_読み込み.ReadBoolean();
                var l_カバレッジ = new int[p_読み込み.Read7BitEncodedInt()];
                for (var j = 0; j < l_カバレッジ.Length; j++)
                {
                    l_カバレッジ[j] = p_読み込み.Read7BitEncodedInt();
                }

                List<int>? l_継ぎ目 = null;
                if (p_読み込み.ReadBoolean())
                {
                    l_継ぎ目 = new List<int>(p_読み込み.Read7BitEncodedInt());
                    for (var j = l_継ぎ目.Capacity; j > 0; j--)
                    {
                        l_継ぎ目.Add(p_読み込み.Read7BitEncodedInt());
                    }
                }

                List<(int A_開始, int A_終了)>? l_範囲群 = null;
                if (p_読み込み.ReadBoolean())
                {
                    l_範囲群 = new List<(int, int)>(p_読み込み.Read7BitEncodedInt());
                    for (var j = l_範囲群.Capacity; j > 0; j--)
                    {
                        l_範囲群.Add((p_読み込み.Read7BitEncodedInt(), p_読み込み.Read7BitEncodedInt()));
                    }
                }

                l_結果.Add(new 引き継ぎ配列(l_配列, l_カバレッジ, l_k長, l_Is確定経路, l_継ぎ目, l_範囲群));
            }

            return l_結果;
        }

        /// <summary>
        /// null になりうる文字列を書く
        /// </summary>
        /// <param name="p_書き込み"></param>
        /// <param name="p_値"></param>
        private static void V_書き込み_省略可(BinaryWriter p_書き込み, string? p_値)
        {
            p_書き込み.Write(p_値 is not null);
            if (p_値 is not null)
            {
                p_書き込み.Write(p_値);
            }
        }

        /// <summary>
        /// null になりうる文字列を読む
        /// </summary>
        /// <param name="p_読み込み"></param>
        /// <returns></returns>
        private static string? Get_省略可(BinaryReader p_読み込み)
        {
            return p_読み込み.ReadBoolean() ? p_読み込み.ReadString() : null;
        }

        #endregion
    }
}
