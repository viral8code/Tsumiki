#!/usr/bin/env python3
"""
合成データを作り、アセンブルし、真の配列と突き合わせるところまでを一括で行う
回帰チェック。誤アセンブリがあれば終了コード 1 を返す。

N50 や総延長だけを見ていると、誤って繋いだときにむしろ数字が良くなるため
品質の劣化に気付けない。実際、反復配列を1回しか使えないまま通り抜けて
中間を飛ばした contig が出力され、N50 が 99,974 から 199,945 へ「改善」して
いたことがある。真の配列と突き合わせて初めて分かる種類の不具合なので、
変更のたびにこれを回す。

使い方:
    python tools/ground_truth_check.py --exe path/to/Tsumiki.exe [--work-dir DIR]

反復を含む条件と含まない条件の両方を試す。反復ありの条件は、上記の
「反復を通り抜ける誤連結」を再現する最小のケースになっている。
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (名前, シミュレータへの追加引数, アセンブラの -k)
CASES = [
    ("反復なし", ["--repeat-count", "0"], 63),
    ("2コピー反復150bp", ["--repeat-count", "2", "--repeat-length", "150"], 63),
    ("3コピー反復150bp", ["--repeat-count", "3", "--repeat-length", "150"], 63),
]


def run(cmd, cwd=None):
    # Windows の日本語ロケールではコンソール出力が CP932 になりうる。
    # UTF-8 決め打ちで読むとデコードに失敗して検査自体が落ちるため、
    # 置換文字で読み飛ばす(検査に必要なのは終了コードと標準出力の要約だけ)。
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                            encoding='utf-8', errors='replace')
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(f"コマンドが失敗した: {' '.join(str(c) for c in cmd)}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True, help="Tsumiki の実行ファイル")
    parser.add_argument("--work-dir", default="ground_truth_work")
    parser.add_argument("--genome-length", type=int, default=300_000)
    parser.add_argument("--coverage", type=float, default=60.0)
    parser.add_argument("--error-rate", type=float, default=0.005)
    args = parser.parse_args()

    work = Path(args.work_dir).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    failures = []
    for name, extra, kmer in CASES:
        print("=" * 70)
        print(f"ケース: {name}")
        print("=" * 70)

        case_dir = work / name.replace(" ", "_")
        reads_dir = case_dir / "reads"
        asm_dir = case_dir / "asm"
        reads_dir.mkdir(parents=True)
        asm_dir.mkdir(parents=True)

        run([sys.executable, str(HERE / "simulate_reads.py"),
             "--out-dir", str(reads_dir),
             "--genome-length", str(args.genome_length),
             "--coverage", str(args.coverage),
             "--error-rate", str(args.error_rate),
             *extra])

        run([args.exe,
             "-1", str(reads_dir / "reads.1.fq"),
             "-2", str(reads_dir / "reads.2.fq"),
             "-k", str(kmer), "-kc", "3", "-th", "4",
             "-t", "tmp"], cwd=str(asm_dir))

        assembly = asm_dir / "scaffolds.fasta"
        if not assembly.exists():
            assembly = asm_dir / "contigs.fasta"

        result = subprocess.run(
            [sys.executable, str(HERE / "validate_assembly.py"),
             "--reference", str(reads_dir / "reference.fasta"),
             "--assembly", str(assembly)],
            capture_output=True, text=True, encoding='utf-8', errors='replace')
        print(result.stdout)
        if result.returncode != 0:
            failures.append(name)

    print("=" * 70)
    if failures:
        print(f"誤アセンブリが見つかったケース: {', '.join(failures)}")
        return 1
    print("すべてのケースで誤アセンブリなし。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
