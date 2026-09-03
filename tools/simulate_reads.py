#!/usr/bin/env python3
"""
Local, self-contained synthetic ground-truth generator for developing/testing
Tsumiki's assembly pipeline (error correction, graph simplification, etc.)
without needing any external download or real (confidential) data.

Generates a random reference genome and simulates paired-end short reads from
it with a configurable per-base substitution error rate and Normal-distributed
insert size. Also writes a sidecar TSV recording every injected error
(read id, 0-based position within the read, true base, mutated base) so a
correction step's output can be checked directly against ground truth
(precision/recall), not just eyeballed via N50.

Usage:
    python tools/simulate_reads.py --out-dir <dir> [options]

Outputs (in --out-dir):
    reference.fasta   the true genome
    reads.1.fq         forward reads
    reads.2.fq         reverse reads
    errors.tsv          read_id \t mate(1|2) \t position(0-based) \t true_base \t mutated_base
"""

import argparse
import random
from pathlib import Path

BASES = "ACGT"
TAB = chr(9)
NEWLINE = chr(10)
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def random_genome(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def insert_repeats(genome: str, count: int, repeat_length: int, rng: random.Random):
    """
    ゲノム中に、同一配列の反復を count 箇所に埋め込む。

    反復配列はアセンブリが途切れる主要因であり、de Bruijn グラフ上では
    1個の頂点に潰れて入次数・出次数が複数になる。リード自体は反復の内側から
    読まれてもどのコピー由来か区別できないため、これを解くにはコピー全体を
    跨いだペアエンドの情報が要る。repeat resolution が正しく働くかを
    真値つきで確かめるために、その状況を意図的に作る。

    埋め込み位置は互いに十分離し(反復同士が隣接して構造が入れ子にならないよう)、
    元のゲノム長は変えない(同じ長さの区間を置き換える)。
    戻り値は (反復入りゲノム, [埋め込み開始位置...], 反復配列)。
    """
    if count <= 0 or repeat_length <= 0:
        return genome, [], ""

    unit = "".join(rng.choice(BASES) for _ in range(repeat_length))
    bases = list(genome)

    # 位置は等間隔に散らし、周囲に十分な固有配列が残るようにする。
    span = len(genome) // (count + 1)
    if span <= repeat_length * 3:
        raise SystemExit("genome is too short for the requested number/length of repeats")

    positions = []
    for i in range(count):
        start = span * (i + 1)
        bases[start:start + repeat_length] = list(unit)
        positions.append(start)

    return "".join(bases), positions, unit


def mutate_base(true_base: str, rng: random.Random) -> str:
    choices = [b for b in BASES if b != true_base]
    return rng.choice(choices)


def random_quality(length: int, rng: random.Random, min_q: int, max_q: int) -> str:
    # Phred33。読みの端でクオリティが落ちる典型的な傾向を軽く模す。
    quals = []
    for i in range(length):
        edge_penalty = 4 if i < 5 or i >= length - 5 else 0
        q = max(2, rng.randint(min_q, max_q) - edge_penalty)
        quals.append(chr(q + 33))
    return "".join(quals)


def apply_errors(seq: str, error_rate: float, rng: random.Random):
    """seq に置換エラーを注入する。戻り値は (変異後の配列, [(position, true_base, mutated_base), ...])。"""
    bases = list(seq)
    injected = []
    for i, b in enumerate(bases):
        if rng.random() < error_rate:
            mutated = mutate_base(b, rng)
            injected.append((i, b, mutated))
            bases[i] = mutated
    return "".join(bases), injected


def write_fasta(path: Path, seq_id: str, seq: str):
    with path.open("w") as f:
        f.write(f">{seq_id}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", required=True, help="output directory (created if missing)")
    parser.add_argument("--genome-length", type=int, default=300_000, help="synthetic genome length in bp (default: 300000)")
    parser.add_argument("--coverage", type=float, default=50.0, help="target sequencing depth (default: 50x)")
    parser.add_argument("--read-length", type=int, default=150, help="read length in bp (default: 150)")
    parser.add_argument("--insert-size", type=int, default=350, help="mean fragment insert size (default: 350)")
    parser.add_argument("--insert-sd", type=float, default=30.0, help="insert size standard deviation (default: 30)")
    parser.add_argument("--error-rate", type=float, default=0.01, help="per-base substitution error rate (default: 0.01)")
    parser.add_argument("--min-qual", type=int, default=30, help="min simulated Phred quality (default: 30)")
    parser.add_argument("--max-qual", type=int, default=40, help="max simulated Phred quality (default: 40)")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for reproducibility (default: 42)")
    parser.add_argument("--repeat-count", type=int, default=0,
                        help="number of identical repeat copies to embed in the genome (default: 0 = no repeats)")
    parser.add_argument("--repeat-length", type=int, default=200,
                        help="length of each embedded repeat copy in bp (default: 200)")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genome = random_genome(args.genome_length, rng)
    genome, repeat_positions, repeat_unit = insert_repeats(
        genome, args.repeat_count, args.repeat_length, rng)
    if repeat_positions:
        print(f"embedded {len(repeat_positions)} identical repeat(s) of {len(repeat_unit)}bp "
              f"at positions {repeat_positions}")
        rows = ["start" + TAB + "length"]
        rows += [str(pos) + TAB + str(len(repeat_unit)) for pos in repeat_positions]
        (out_dir / "repeats.tsv").write_text(NEWLINE.join(rows) + NEWLINE)
    write_fasta(out_dir / "reference.fasta", "synthetic_reference", genome)

    num_pairs = int((args.genome_length * args.coverage) / (2 * args.read_length))

    with (out_dir / "reads.1.fq").open("w") as f1, \
         (out_dir / "reads.2.fq").open("w") as f2, \
         (out_dir / "errors.tsv").open("w") as ferr:
        ferr.write("read_id\tmate\tposition\ttrue_base\tmutated_base\n")

        written = 0
        attempts = 0
        max_attempts = num_pairs * 20
        while written < num_pairs and attempts < max_attempts:
            attempts += 1
            frag_len = max(args.read_length * 2, int(rng.gauss(args.insert_size, args.insert_sd)))
            if frag_len > args.genome_length:
                continue
            start = rng.randint(0, args.genome_length - frag_len)
            fragment = genome[start:start + frag_len]

            read1_true = fragment[:args.read_length]
            read2_true = revcomp(fragment[-args.read_length:])

            read_id = f"SIM{written}"

            read1_seq, err1 = apply_errors(read1_true, args.error_rate, rng)
            read2_seq, err2 = apply_errors(read2_true, args.error_rate, rng)

            qual1 = random_quality(args.read_length, rng, args.min_qual, args.max_qual)
            qual2 = random_quality(args.read_length, rng, args.min_qual, args.max_qual)

            f1.write(f"@{read_id}/1\n{read1_seq}\n+\n{qual1}\n")
            f2.write(f"@{read_id}/2\n{read2_seq}\n+\n{qual2}\n")

            for pos, true_b, mut_b in err1:
                ferr.write(f"{read_id}\t1\t{pos}\t{true_b}\t{mut_b}\n")
            for pos, true_b, mut_b in err2:
                ferr.write(f"{read_id}\t2\t{pos}\t{true_b}\t{mut_b}\n")

            written += 1

    total_errors = sum(1 for _ in (out_dir / "errors.tsv").open()) - 1
    print(f"Wrote {written} read pairs ({written * 2 * args.read_length} bp, ~{args.coverage:.1f}x of a {args.genome_length}bp genome)")
    print(f"Injected {total_errors} substitution errors (target rate {args.error_rate:.4f})")
    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
