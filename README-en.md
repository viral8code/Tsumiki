# Tsumiki

**Tsumiki** is a de novo genome assembler for short reads (single-end and paired-end).
It builds a de Bruijn graph from trusted k-mers and uses multiple k values and read-pair information to produce contigs and scaffolds.

> [!CAUTION]
> Tsumiki is experimental software under active development. Its accuracy and resource requirements on real data have not yet been thoroughly validated.
> If you use the results for analysis, check correctness against a reference genome with a tool such as QUAST.

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input](#input)
- [Options](#options)
- [Output](#output)
- [Re-running](#re-running)
- [Tips](#tips)
- [License](#license)

## Requirements

- Windows (x64)
- The [.NET 10 SDK](https://dotnet.microsoft.com/download) is required to build
- Linux and macOS are not yet tested

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
| `-1 <path>` | R1 FASTQ (or the reads for single-end data) | **Required** |
| `-2 <path>` | R2 FASTQ | None |
| `-t <path>` | Output directory | `temp` |
| `-th <int>` | Number of threads | Logical processor count |
| `-mem <size>` | Memory limit for k-mer counting (e.g. `512M`, `4G`) | `768M` |
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

With the default profile (`-profile standard`), all of the following are enabled:

- Preprocessing (adapter read-through trimming and overlap-based correction)
- k-mer spectrum error correction
- Multi-k assembly with sequence carry-over between k values
- r-mer verification of repeats
- Local assembly of unclosed gaps
- Rescue of low-count k-mers
- Polishing (substitution error correction)
- Verification of circular junctions
- GFA output

To enable only some of them, start from `-profile legacy`, which disables all of them, and add what you need. Put `-profile` before other options.

```bash
Tsumiki.exe -profile legacy -ec -po -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t out
```

| Option | Description | standard |
|---|---|---|
| `-profile <standard\|legacy>` | Set the combination of stages at once | `standard` |
| `-pp` | Preprocessing | On |
| `-ec` | Error correction | On |
| `-mk` | Multi-k assembly | On |
| `-sr` | Carry synthetic reads built from overlapping pairs to the next k | On |
| `-rv` | Require read support before duplicating repeats | On |
| `-la` | Local assembly of unclosed gaps | On |
| `-my` | Rescue low-count k-mers flanked by trusted k-mers | On |
| `-po` | Polishing | On |
| `-cc` | Verify circular junctions | On |
| `-gfa` | Write the graph in GFA1 format | On |
| `-nc` | Disable sequence carry-over between k values | (carry over) |
| `-nt` | Disable trimming of low-coverage graph ends | (trim) |
| `-mg` | Fill unresolved junctions with sequence from other k values (may increase misassemblies) | Off |
| `-cnb <spectrum\|weighted>` | How the baseline depth for copy-number estimation is chosen | `weighted` |

### Output directory

| Option | Description | Default |
|---|---|---|
| `-rs` | Reuse preprocessed and corrected reads in an existing output directory | Off |
| `-rt` | Delete intermediate files after a successful run (final results and log are kept) | Off |

## Output

Results are written to the directory given by `-t`.

| File | Contents |
|---|---|
| `assembly.fasta` | **Final assembly. Use this file in most cases** |
| `scaffolds.fasta` | Scaffolds (before polishing) |
| `contigs.fasta` | Contigs (before scaffolding) |
| `unitigs.fasta` | Unitigs |
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
