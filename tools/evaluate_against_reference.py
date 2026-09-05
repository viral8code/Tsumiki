#!/usr/bin/env python3
"""
実データのアセンブリを、既知のリファレンスゲノムと突き合わせて評価する。

validate_assembly.py との違い:
あちらは合成データ向けで、contig がリファレンスに「そのまま部分文字列として
現れるか」を見る完全一致判定だった。実データではリファレンスとの間に
SNP・indel・株差が必ずあるため、その判定では全ての contig が「不一致」になり
何も分からない。

ここでは QUAST が行うのと同じ考え方で、k-mer をアンカーにした共線性で評価する:
  1. リファレンスの k-mer のうち、ゲノム中に1回しか現れないものを位置つきで索引する
     (反復配列由来の k-mer は、どの位置に対応するか決められないので使わない)
  2. contig の k-mer を順に引き、当たったアンカーを (contig上の位置, リファレンス上の位置)
     の対として並べる
  3. その対の列が単調に進んでいる限り、同じ「整列ブロック」とみなす。
     参照先の配列が変わる、向きが変わる、あるいは進み方が急に飛ぶ箇所を
     誤アセンブリ(breakpoint)として切る

報告する指標:
  breakpoints    ... 誤アセンブリの箇所数。これが本命の精度指標
  NGA50          ... 誤アセンブリで切ったあとの整列ブロックで計算した N50。
                     GAGE / GAGE-B が "corrected N50" と呼ぶものに相当する。
                     素の N50 は誤って繋ぐほど良くなるため、単体では精度を表さない
  genome fraction... リファレンスのうちアセンブリが覆えた割合
  duplication    ... 整列ブロックの総延長 / 覆えたリファレンス長。
                     1 を大きく超えるなら同じ領域を重複して出力している

使い方:
    python tools/evaluate_against_reference.py --reference ref.fna --assembly scaffolds.fasta
"""

import argparse
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# ヒアドキュメント越しに書き込むとエスケープが二重に解釈されるため、
# 制御文字はリテラルではなく定数で持つ。
TAB = chr(9)
NEWLINE = chr(10)

# アンカーに使う k-mer 長。短すぎると偶然一致が増え、長すぎると SNP に当たって
# アンカーが取れなくなる。31 なら 4^31 通りあり細菌ゲノムで偶然一致はまず起きず、
# かつ SNP 密度が 1/1000 程度なら 3% 程度のアンカーしか失わない。
ANCHOR_K = 31

# 反復配列由来として索引から除外する印。
REPEATED = -1


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
                name = line[1:].split()[0] if len(line) > 1 else ""
                buf = []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        names.append(name)
        seqs.append("".join(buf))
    return names, seqs


CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def iter_canonical_kmers(seq, k):
    """(位置, 正規形の整数表現, 正鎖ならTrue) を順に返す。

    N を含む窓は飛ばす。逆相補と小さいほうを正規形とすることで、
    どちらの向きで書かれていても同じ鍵になる。
    """
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    fwd = rev = 0
    valid = 0
    for i, base in enumerate(seq):
        code = CODE.get(base)
        if code is None:
            valid = 0
            fwd = rev = 0
            continue
        fwd = ((fwd << 2) | code) & mask
        rev = (rev >> 2) | ((3 - code) << shift)
        valid += 1
        if valid >= k:
            if fwd <= rev:
                yield i - k + 1, fwd, True
            else:
                yield i - k + 1, rev, False


def build_reference_index(seqs, k):
    """リファレンスの一意な k-mer を {正規形: (配列番号, 位置, 正鎖か)} で返す。"""
    index = {}
    for seq_id, seq in enumerate(seqs):
        for pos, key, forward in iter_canonical_kmers(seq, k):
            if key in index:
                index[key] = REPEATED
            else:
                index[key] = (seq_id, pos, forward)
    return {key: value for key, value in index.items() if value != REPEATED}


def split_into_blocks(anchors, tolerance):
    """アンカー列を、共線性が保たれている区間ごとに切り分ける。

    anchors: (contig上の位置, 配列番号, リファレンス上の位置, 向きが一致するか)
    リファレンス上の進み方が contig 上の進み方と食い違う量が tolerance を
    超えたら、そこを誤アセンブリとみなして切る。
    """
    blocks = []
    current = []
    for anchor in anchors:
        if not current:
            current = [anchor]
            continue
        prev = current[-1]
        same_place = anchor[1] == prev[1] and anchor[3] == prev[3]
        if same_place:
            contig_step = anchor[0] - prev[0]
            ref_step = (anchor[2] - prev[2]) if anchor[3] else (prev[2] - anchor[2])
            if ref_step > 0 and abs(ref_step - contig_step) <= tolerance:
                current.append(anchor)
                continue
        blocks.append(current)
        current = [anchor]
    if current:
        blocks.append(current)
    return blocks


def n50(lengths, total):
    """長さの列に対する N50。total を分母に使うと NGA50 になる。"""
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= total:
            return length
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--min-contig", type=int, default=500,
                        help="この長さ未満の contig は無視する(abyss-fac 等と揃えるため)")
    parser.add_argument("--min-block", type=int, default=200,
                        help="この長さ未満の整列ブロックは雑音とみなして数えない")
    parser.add_argument("--tolerance", type=int, default=1000,
                        help="共線性のずれをこの範囲まで許す(indel を誤アセンブリと呼ばないため)")
    parser.add_argument("--label", default="")
    parser.add_argument("--breakpoints", default=None,
                        help="誤アセンブリの位置を TSV で書き出す先(contig名, contig上の位置)")
    args = parser.parse_args()

    ref_names, ref_seqs = load_fasta(args.reference)
    ref_total = sum(len(s) for s in ref_seqs)
    index = build_reference_index(ref_seqs, ANCHOR_K)

    asm_names, asm_seqs = load_fasta(args.assembly)
    kept = [(n, s) for n, s in zip(asm_names, asm_seqs) if len(s) >= args.min_contig]

    covered = [bytearray(len(s)) for s in ref_seqs]
    block_lengths = []
    breakpoint_count = 0
    unaligned = 0
    contigs_with_breakpoints = 0

    breakpoints = []
    for name, seq in kept:
        anchors = []
        for pos, key, forward in iter_canonical_kmers(seq, ANCHOR_K):
            hit = index.get(key)
            if hit is not None:
                seq_id, ref_pos, ref_forward = hit
                anchors.append((pos, seq_id, ref_pos, forward == ref_forward))
        if not anchors:
            unaligned += 1
            continue

        blocks = [b for b in split_into_blocks(anchors, args.tolerance)
                  if b[-1][0] - b[0][0] + ANCHOR_K >= args.min_block]
        if not blocks:
            unaligned += 1
            continue

        if len(blocks) > 1:
            breakpoint_count += len(blocks) - 1
            contigs_with_breakpoints += 1
            # 切れ目は、直前ブロックの末尾と次ブロックの先頭の中間に置く。
            for 前, 後 in zip(blocks, blocks[1:]):
                breakpoints.append((name, (前[-1][0] + 後[0][0]) // 2))

        for block in blocks:
            block_lengths.append(block[-1][0] - block[0][0] + ANCHOR_K)
            seq_id = block[0][1]
            low = min(a[2] for a in block)
            high = max(a[2] for a in block) + ANCHOR_K
            for i in range(low, min(high, len(covered[seq_id]))):
                covered[seq_id][i] = 1

    covered_bp = sum(sum(c) for c in covered)
    aligned_total = sum(block_lengths)
    asm_total = sum(len(s) for _, s in kept)

    label = f"[{args.label}] " if args.label else ""
    print(f"{label}reference: {args.reference}")
    print(f"{label}  reference sequences : {len(ref_seqs)}, total {ref_total:,} bp")
    print(f"{label}  unique {ANCHOR_K}-mer anchors : {len(index):,}")
    print(f"{label}assembly: {args.assembly}")
    print(f"{label}  contigs (>= {args.min_contig} bp) : {len(kept)}, total {asm_total:,} bp")
    print(f"{label}  N50 (raw)          : {n50([len(s) for _, s in kept], asm_total):,}")
    print(f"{label}  misassemblies      : {breakpoint_count} breakpoint(s) in {contigs_with_breakpoints} contig(s)")
    print(f"{label}  NGA50 (corrected)  : {n50(block_lengths, ref_total):,}")
    print(f"{label}  aligned blocks     : {len(block_lengths)}, total {aligned_total:,} bp")
    print(f"{label}  genome fraction    : {covered_bp / ref_total * 100:.2f}% ({covered_bp:,} / {ref_total:,} bp)")
    if covered_bp:
        print(f"{label}  duplication ratio  : {aligned_total / covered_bp:.3f}")
    print(f"{label}  contigs with no anchor : {unaligned}")

    if args.breakpoints:
        with open(args.breakpoints, "w") as f:
            f.write("contig" + TAB + "position" + NEWLINE)
            for name, pos in breakpoints:
                f.write(f"{name}{TAB}{pos}{NEWLINE}")
        print(f"{label}  breakpoints written to : {args.breakpoints}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
