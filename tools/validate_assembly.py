#!/usr/bin/env python3
"""
真の配列が分かっている合成データに対して、アセンブリ結果に誤アセンブリが
無いかを検査する。

なぜ必要か:
N50 や総延長は「どれだけ長く繋がったか」しか見ておらず、誤って繋いだ場合
むしろ良い数字になる。実際、反復配列 R が2回現れる A-R-B-R-C というゲノムで
「A-R-C」という中間を飛ばした contig が出力され、N50 は 99,974 から 199,945 へ
「改善」していた。個々の隣接(A→R と R→C)はどちらも本物なので、辺を検査する
だけでは見抜けない。真の配列と突き合わせて初めて分かる。

判定:
  完全一致    ... contig 全体が真の配列(順鎖または逆相補)にそのまま現れる
  キメラ      ... 途中までは一致するが、その先が別の場所に由来する = 誤った連結
  不一致      ... 一致する部分がほとんど無い(エラーが残っている等)

使い方:
    python tools/validate_assembly.py --reference ref.fasta --assembly scaffolds.fasta
終了コード 0 = 誤アセンブリなし、1 = 誤アセンブリあり。
"""

import argparse
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def load_fasta(path):
    names, seqs, buf, name = [], [], [], None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    names.append(name)
                    seqs.append("".join(buf))
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        names.append(name)
        seqs.append("".join(buf))
    return names, seqs


def longest_matching_prefix(contig, reference):
    """contig の先頭から、reference 内にそのまま現れる最長の長さを二分探索で求める。"""
    lo, hi = 0, len(contig)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if contig[:mid] in reference:
            lo = mid
        else:
            hi = mid - 1
    return lo


def classify(contig, reference):
    """1本の contig を判定する。戻り値は (種別, 詳細文字列)。"""
    for oriented, label in ((contig, "順鎖"), (revcomp(contig), "逆相補")):
        if oriented in reference:
            start = reference.index(oriented)
            return "ok", f"ref[{start}:{start + len(oriented)}] {label} と完全一致"

    # 完全一致しない。どちらの向きでより長く一致するかを見て、破断点を特定する。
    best_len, best_seq, best_label = 0, None, None
    for oriented, label in ((contig, "順鎖"), (revcomp(contig), "逆相補")):
        n = longest_matching_prefix(oriented, reference)
        if n > best_len:
            best_len, best_seq, best_label = n, oriented, label

    if best_len < 100:
        return "mismatch", f"真の配列とほとんど一致しない(最長一致 {best_len}bp)"

    head_start = reference.index(best_seq[:best_len])
    remainder = best_seq[best_len:]
    # 破断点の先が別の場所に由来していれば、それは誤った連結(キメラ)。
    probe = remainder[:200] if len(remainder) >= 200 else remainder
    if probe and probe in reference:
        tail_start = reference.index(probe)
        return "chimera", (
            f"{best_label}: 先頭 {best_len}bp は ref[{head_start}:{head_start + best_len}] に一致するが、"
            f"その先は ref[{tail_start}:] 由来。つまり離れた2箇所を誤って連結している"
        )
    return "chimera", (
        f"{best_label}: 先頭 {best_len}bp までしか一致せず、残り {len(remainder)}bp は"
        " 真の配列上に見つからない"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--min-length", type=int, default=0,
                        help="この長さ未満の配列は検査対象から外す(既定: 0 = すべて)")
    args = parser.parse_args()

    _, refs = load_fasta(args.reference)
    if len(refs) != 1:
        print(f"注意: リファレンスが {len(refs)} 本ある。連結して扱う。")
    reference = "".join(refs)

    names, seqs = load_fasta(args.assembly)
    pairs = [(n, s) for n, s in zip(names, seqs) if len(s) >= args.min_length]

    counts = {"ok": 0, "chimera": 0, "mismatch": 0}
    covered = 0
    print(f"リファレンス {len(reference):,}bp / アセンブリ {len(pairs)} 本 "
          f"{sum(len(s) for _, s in pairs):,}bp\n")

    for name, seq in sorted(pairs, key=lambda p: -len(p[1])):
        kind, detail = classify(seq, reference)
        counts[kind] += 1
        if kind == "ok":
            covered += len(seq)
        mark = {"ok": "  ", "chimera": "★ ", "mismatch": "▲ "}[kind]
        print(f"{mark}{name:<16} {len(seq):>9,}bp  {detail}")

    print()
    print(f"完全一致: {counts['ok']} 本   キメラ(誤連結): {counts['chimera']} 本   "
          f"不一致: {counts['mismatch']} 本")
    print(f"完全一致した配列が覆うリファレンス: {covered:,}bp "
          f"({100.0 * covered / len(reference):.2f}%)")

    bad = counts["chimera"] + counts["mismatch"]
    if bad:
        print(f"\n誤アセンブリが {bad} 本ある。")
        return 1
    print("\n誤アセンブリなし。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
