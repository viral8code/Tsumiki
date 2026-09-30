# Tsumiki

**Tsumiki** is a de novo genome assembler for short reads (single-end and paired-end).
It builds a de Bruijn graph from trusted k-mers and uses multiple k values and read-pair information to produce contigs and scaffolds.

> [!CAUTION]
> Tsumiki is experimental software under active development. Validation is limited to the twelve bacterial short-read datasets (HiSeq and MiSeq) in [Benchmark](#benchmark);
> behaviour outside those conditions (eukaryotes, metagenomes, reads longer than 250 bp, and so on) is unknown.
> If you use the results for analysis, check correctness against a reference genome with a tool such as QUAST.

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Benchmark](#benchmark)
- [Input](#input)
- [Options](#options)
- [Output](#output)
- [Re-running](#re-running)
- [Tips](#tips)
- [License](#license)

## Requirements

- Windows (x64)
- Linux (x64): tested on Ubuntu 24.04. The executables on the Releases page require glibc 2.34 or later (Ubuntu 22.04, Debian 12, RHEL 9, or newer)
- The [.NET 10 SDK](https://dotnet.microsoft.com/download) is required to build
- Linux (arm64) and macOS are not yet tested

## Installation

Build from source.

```bash
git clone https://github.com/viral8code/Tsumiki.git
```

```bash
cd Tsumiki
```

To build a self-contained executable (Native AOT), you need the "Desktop development with C++" workload of Visual Studio.

```bash
dotnet publish Tsumiki.csproj -c Release -o publish
```

Without a C++ toolchain, disable Native AOT. `Tsumiki.exe` then requires the .NET 10 runtime.

```bash
dotnet publish Tsumiki.csproj -c Release -p:PublishAot=false -o publish
```

The installation is complete if the help message is shown.

```bash
publish/Tsumiki.exe -h -lang en
```

## Quick start

Paired-end reads:

```bash
Tsumiki.exe -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t out -lang en
```

Single-end reads:

```bash
Tsumiki.exe -1 reads.fastq.gz -t out -lang en
```

Specifying threads, memory, and k values:

```bash
Tsumiki.exe -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -k 31,51,71 -th 16 -mem 4G -t out -lang en
```

The final assembly is written to `out/assembly.fasta`.

By default, all major stages (preprocessing, error correction, multi-k, local gap assembly, polishing, and so on) are enabled, so specifying the input and output directory is usually enough.

## Benchmark

Results for the eight HiSeq datasets and four MiSeq datasets of [GAGE-B](https://ccb.jhu.edu/gage_b/) (raw reads, no trimming), assembled with default settings. Evaluated with [QUAST](https://quast.sourceforge.net/) 5.3.0 (`--min-contig 500`) on 16 threads with `-mem 8G`.

| Dataset | Platform | GC% | Depth | Reference level | Contigs | Total length | N50 | NA50 | Misassemblies | Local misassemblies | Genome fraction |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| A. hydrophila SSU | HiSeq | 61.5 | 256x | Scaffold (2) | 42 | 4,846,968 | 272,482 | 272,482 | 1 | 5 | 98.43% |
| B. cereus VD118 | HiSeq | 35.3 | 242x | Scaffold (10) | 215 | 5,633,188 | 70,173 | 68,612 | 10 | 9 | 98.40% |
| B. fragilis HMW 615 | HiSeq | 43.5 | 266x | Scaffold (14) | 142 | 5,277,391 | 116,434 | 107,957 | 15 | 4 | 97.96% |
| M. abscessus 6G-0125-R | HiSeq | 64.1 | 105x | Contig (5) | 60 | 5,127,199 | 165,957 | 165,957 | 0 | 0 | 99.68% |
| R. sphaeroides 2.4.1 | HiSeq | 68.8 | 224x | **Complete** (7) | 131 | 4,591,173 | 128,436 | 128,436 | 1 | 1 | 98.66% |
| S. aureus M0927 | HiSeq | 32.8 | 301x | Scaffold (12) | 62 | 2,829,974 | 122,540 | 122,540 | 2 | 1 | 98.36% |
| V. cholerae CP1032(5) | HiSeq | 47.5 | 94x | Contig (17) | 107 | 3,917,494 | 97,957 | 97,251 | 3 | 1 | 98.34% |
| X. axonopodis UA323 | HiSeq | 65.1 | 332x | Contig (151) | 142 | 4,906,110 | 108,063 | 61,310 | 57 | 2 | 99.83% |
| B. cereus ATCC 10987 | MiSeq | 35.5 | 89x | **Complete** (2) | 70 | 5,320,399 | 157,085 | 157,082 | 3 | 2 | 97.79% |
| M. abscessus 6G-0125-R | MiSeq | 64.1 | 100x | Contig (5) | 95 | 5,140,051 | 104,862 | 104,862 | 1 | 0 | 99.83% |
| R. sphaeroides 2.4.1 | MiSeq | 68.8 | 101x | **Complete** (7) | 112 | 4,554,974 | 113,007 | 113,007 | 2 | 1 | 98.50% |
| V. cholerae CP1032(5) | MiSeq | 47.5 | 101x | Contig (17) | 125 | 3,937,853 | 106,307 | 106,307 | 6 | 1 | 98.79% |

HiSeq reads are 100-101 bp and MiSeq reads are 250-251 bp. The B. cereus HiSeq and MiSeq datasets are different strains.

### Mixing HiSeq and MiSeq

For the three species with HiSeq and MiSeq data from the same strain, both sets of reads were also given to a single run (same conditions as above).

| Dataset | Depth | Reference level | Contigs | Total length | N50 | NA50 | Misassemblies | Local misassemblies | Genome fraction |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| M. abscessus 6G-0125-R | 205x | Contig (5) | 48 | 5,136,173 | 249,022 | 249,022 | 1 | 0 | 99.82% |
| R. sphaeroides 2.4.1 | 325x | **Complete** (7) | 117 | 4,582,568 | 142,745 | 142,745 | 2 | 1 | 99.17% |
| V. cholerae CP1032(5) | 195x | Contig (17) | 81 | 3,926,794 | 246,185 | 246,185 | 5 | 1 | 98.60% |

For all three species, NA50 is larger than with either platform alone (at least 1.5 times for M. abscessus and V. cholerae). The misassembly count falls between those of HiSeq alone and MiSeq alone.

### How to read these numbers

**Reference level** describes how completely the reference genome itself is assembled. `Complete` is a finished genome, `Scaffold` is supercontigs with N-filled gaps, and `Contig` is fragments with no adjacency information; the number in parentheses is the number of sequences in the reference.

**The more fragmented the reference, the more misassemblies are counted.** A break in the reference is reported as a disagreement even where the sample is genuinely contiguous. Across these twelve datasets the correlation between the number of reference sequences and the misassembly count is **+0.97**. The 57 misassemblies for X. axonopodis largely reflect its reference being split into 151 sequences, while the same assembly has a genome fraction of 99.83%, the highest level of the twelve. Misassembly counts are therefore not comparable across datasets.

Finished reference genomes are available for R. sphaeroides and B. cereus ATCC 10987 (MiSeq); there the misassembly counts are 1 (HiSeq), 2 (MiSeq), and 3 respectively.

**Misassembly counts also include structural differences between the sequenced strain and the reference.** For R. sphaeroides, the same junction joining positions around 36.5 kb and 64.5 kb of the reference is counted once in each of the HiSeq, MiSeq, and mixed assemblies. The 40 bases spanning this junction occur verbatim in 176 raw reads (HiSeq), while the reference contains nothing similar. This is most likely a genuine difference in the sequenced strain rather than an assembly error.

**Depth does not drive accuracy.** Depth spans 89x to 332x, a 3.7-fold range, yet its correlation with NA50 is **-0.06**.

### Reproducing

```bash
Tsumiki.exe -1 <reads_1.fastq> -2 <reads_2.fastq> -mg -th 16 -mem 8G -t <output> -lang en
```

To mix HiSeq and MiSeq, list the files comma-separated for `-1` and `-2`.

```bash
Tsumiki.exe -1 <hiseq_1.fastq>,<miseq_1.fastq> -2 <hiseq_2.fastq>,<miseq_2.fastq> -mg -th 16 -mem 8G -t <output> -lang en
```

```bash
quast.py -o <output> -r <reference.fna> --min-contig 500 <output>/assembly.fasta
```

## Input

- FASTQ (gzip-compressed `.gz` files are supported)
- For paired-end data, pass R1 and R2 as separate files with reads in the same order (interleaved files are not supported)
- The quality offset (Phred+33 / +64) is detected automatically
- By default, regions containing ambiguous bases such as `N` are excluded from k-mer counting

## Options

Run `Tsumiki.exe -h` for the full list.

### Basic

| Option | Description | Default |
|---|---|---|
| `-1 <path>` | R1 FASTQ (or the reads for single-end data). Comma-separate to give several libraries | **Required** |
| `-2 <path>` | R2 FASTQ. Comma-separate in the same order as `-1` | None |
| `-s <path>` | Single-end FASTQ. Comma-separate to give several; they are added to the `-1`/`-2` libraries | None |
| `-t <path>` | Output directory | `temp` |
| `-th <int>` | Number of threads | Logical processor count |
| `-mem <size>` | Memory limit for k-mer counting (e.g. `512M`, `4G`) | 1/8 of available memory, between `768M` and `4G` |
| `-lang <ja\|en\|zh>` | Message language | `ja` |
| `-log <quiet\|normal\|verbose>` | Console verbosity | `normal` |
| `-v` | Show version | |
| `-h` | Show help | |

### k-mer and quality

| Option | Description | Default |
|---|---|---|
| `-k <int[,int...]>` | k-mer length. With a comma-separated list, assembles with each k and keeps the best result | Selected from read length |
| `-kc <int>` | Minimum count for a trusted k-mer | Selected from the k-mer spectrum |
| `-q <int>` | Minimum trusted base quality | `1` |
| `-qt <int>` | During preprocessing, trim the 3' ends of paired reads where quality falls below this value (same method as BWA `-q`). `0` disables it. Not applied to single-end reads (`-s`) | `0` |
| `-p <33\|64>` | Phred offset | Auto-detected |
| `-ab` | Expand ambiguous (IUPAC) bases instead of skipping them | Off |

### Paired-end and scaffolding

| Option | Description | Default |
|---|---|---|
| `-mode <conservative\|normal\|bold>` | How aggressively contigs are joined. `conservative` reduces misjoins; `bold` produces longer sequences | `normal` |
| `-i <int>` | Insert size (bp) | Estimated from mapped pairs |
| `-pu <decimal>` | Dominance ratio required to accept a connection (overrides `-mode` if placed after it) | `0.8` |
| `-pc <int>` | Read pairs required to resolve a short repeat (overrides `-mode` if placed after it) | `10` |

Values set by each `-mode` preset:

| `-mode` | `-pu` | `-pc` |
|---|---|---|
| `conservative` | 0.9 | 15 |
| `normal` | 0.8 | 10 |
| `bold` | 0.65 | 5 |

### Enabling and disabling stages

All of the following stages are enabled by default. To disable one, use the corresponding option starting with `-n`.

| Option | Stage it disables |
|---|---|
| `-npp` | Preprocessing (adapter read-through trimming and overlap-based correction) |
| `-nec` | k-mer spectrum error correction |
| `-nmk` | Multi-k assembly (assemble with a single k) |
| `-nc` | Sequence carry-over between k values |
| `-nsr` | Carry-over of synthetic reads built from overlapping pairs |
| `-nrv` | Read support required before duplicating repeats (r-mer verification) |
| `-nla` | Local assembly of unclosed gaps |
| `-nmy` | Rescue of low-count k-mers flanked by trusted k-mers |
| `-nt` | Trimming of low-coverage graph ends |
| `-npo` | Polishing (substitution error correction) |
| `-ncc` | Verification of circular junctions |
| `-ngfa` | Graph output in GFA1 format |

```bash
Tsumiki.exe -npo -ngfa -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t out
```

Stages that are off by default and other settings:

| Option | Description | Default |
|---|---|---|
| `-mg` | Fill unresolved junctions with sequence from other k values (may increase misassemblies) | Off |
| `-cnb <spectrum\|weighted>` | How the baseline depth for copy-number estimation is chosen | `weighted` |

### Output directory and intermediate data

| Option | Description | Default |
|---|---|---|
| `-rs` | Reuse preprocessed and corrected reads in an existing output directory | Off |
| `-rt` | Delete intermediate files after a successful run (final results and log are kept) | Off |
| `-inmem` | Keep intermediate reads and k-mer counting runs in memory instead of on disk. Greatly reduces disk reads and writes at the cost of more memory | Off |

### Options that cannot be combined

If any of these combinations is given, Tsumiki prints the reason and exits without starting.

| Combination | Reason |
|---|---|
| `-inmem` and `-mem` | `-inmem` keeps k-mer counting runs in memory, whereas `-mem` assumes they spill to disk |
| `-inmem` and `-rs` | `-inmem` leaves no intermediate files on disk to resume from |
| `-nmk` and `-k` with several values | A comma-separated `-k` means "try each value and keep the best" |
| `-i` and two or more paired-end libraries | A single value would be applied to every library, breaking the assumptions for libraries with different insert sizes |

## Output

Results are written to the directory given by `-t`.

| File | Contents |
|---|---|
| `assembly.fasta` | **Final assembly. Use this file in most cases** |
| `scaffolds.fasta` | Scaffolds (before polishing) |
| `contigs.fasta` | Contigs (before scaffolding) |
| `unitigs.fasta` | Unitigs |
| `assembly.agp` | Layout of `assembly.fasta` (AGP v2.1). Gaps of estimated length are type `N`; links supported by read pairs whose overlap could not be confirmed are written as 100 `N`s and marked as type `U` (unknown length) |
| `assembly.gfa` | Assembly graph (GFA1), viewable with tools such as [Bandage](https://rrwick.github.io/Bandage/) |
| `assembly.report.json` | Sequence statistics and verification results |
| `assembly.ambiguous.tsv` | Locations that remained ambiguous |
| `assembly.unsupported.tsv` | Intervals not supported by reads |
| `assembly.provenance.json` | Record of settings, input file hashes, and so on |
| `Tsumiki.log` | Full execution log (saved regardless of `-log`) |

Intermediate results for each k (e.g. `k31/`) and verification files are also created. Use `-rt` to remove them automatically.

## Re-running

If the output directory already exists, Tsumiki stops to avoid mixing results.
To re-run while reusing the previous preprocessing and error-correction results, run with the same input and options plus `-rs`.

```bash
Tsumiki.exe -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t out -rs
```

If the input files or settings have changed, the affected stages are redone automatically.
Only preprocessing and error correction are reused; the assembly itself always runs from the beginning.

## Tips

- **Memory**: `-mem` limits k-mer counting only. Graph construction and other stages use additional memory, so actual usage will be higher.
- **Run time**: Multi-k repeats the assembly for several k values and takes longer. To save time, specify a single k, e.g. `-k 55`.
- **Choosing k**: If you specify `-k`, use values smaller than the read length.
- **Checking correctness**: A large N50 or a sequence marked as circular does not by itself mean the assembly is correct.

## License

MIT License. See [LICENSE.txt](LICENSE.txt) for details.
