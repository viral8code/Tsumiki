# Tsumiki

Tsumiki is an experimental genome assembler for single-end and paired-end short reads. It builds a de Bruijn graph from trusted k-mers, resolves supported connections, and produces contigs, scaffolds, and an evidence report. The implementation is written in C#.

Tsumiki is under active development. Its internal completeness labels describe checks against the input library; they do not establish a finished genome or replace reference-based assessment. Comparative accuracy and resource requirements on real datasets have not yet been established.

## Features

- Automatic k-mer length and abundance cutoff selection, with explicit overrides.
- Disk-backed k-mer counting with a configurable counting memory budget.
- Graph simplification, unitig construction, copy-number-aware repeat handling, and paired-end scaffolding.
- Optional overlap-based read preprocessing and k-mer-based error correction.
- Optional multi-k assembly, sequence carry-over, and candidate selection.
- Optional local gap assembly with bounded context-length refinement, substitution polishing, and circular-junction verification.
- FASTA output, optional GFA1 graph output, JSON reports, and tab-separated evidence summaries.
- Japanese, English, and Chinese command-line messages.

## Build

The current project targets **.NET 10 for Windows** (`net10.0-windows`). Install the .NET 10 SDK, then run from the repository root:

```powershell
dotnet build Tsumiki.csproj -c Release
dotnet run --project Tsumiki.csproj -c Release --no-build -- -h -lang en
```

To publish a framework-dependent build without Native AOT:

```powershell
dotnet publish Tsumiki.csproj -c Release -p:PublishAot=false -o publish
dotnet publish/Tsumiki.dll -h -lang en
```

The project enables Native AOT by default. Native publishing requires the matching native compilation toolchain. The commands above explicitly disable it for the published build. Linux and macOS builds are not currently documented or validated.

## Quick start

Paired-end assembly:

```powershell
dotnet publish/Tsumiki.dll -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t assembly-run -lang en
```

Single-end assembly:

```powershell
dotnet publish/Tsumiki.dll -1 reads.fastq.gz -t single-end-run -lang en
```

Explicit k values, preprocessing, correction, and final checks:

```powershell
dotnet publish/Tsumiki.dll -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -k 31,51,63 -pp -ec -po -cc -gfa -t multi-k-run -th 8 -mem 2G -lang en
```

This last example enables optional processing; it is not a validated best-performing preset. Choose k values below the usable read length. A comma-separated `-k` list automatically enables multi-k execution.

Input files must be FASTQ, optionally gzip-compressed. Supply mates in corresponding record order in separate files. By default, windows containing ambiguous or low-quality bases are skipped during counting. `-ab` expands IUPAC ambiguity into candidate bases, up to 64 combinations per window; low-quality windows and larger ambiguity expansions are skipped. Expanded alternatives are hypotheses, not independent observations. Keep the original reads outside the output directory.

## Command-line options

`-h` prints the complete option list. Options take separate arguments; boolean switches enable a feature by their presence.

| Option | Description | Default |
|---|---|---|
| `-1 <path>` | Forward or single-end FASTQ | Required |
| `-2 <path>` | Reverse FASTQ | None |
| `-ab` | Allow ambiguous input bases | Off |
| `-k <n[,n...]>` | k-mer length or ascending candidate set | Inferred from reads, automatic initial k capped at 63 |
| `-kc <n>` | Minimum trusted k-mer count | Inferred from the spectrum |
| `-p <33|64>` | Explicit Phred encoding offset | Inferred, with 33 as fallback |
| `-q <n>` | Minimum trusted base quality | 1 |
| `-mem <size>` | Counting memory budget, e.g. `512M`, `2G` | 768 MB |
| `-th <n>` | Worker thread count | Logical processor count |
| `-i <n>` | Expected insert size | Estimated from mapped pairs |
| `-pu <ratio>` | Dominance required for paired connections | 0.8 |
| `-pc <n>` | Minimum pair support for short-repeat resolution | 10 |
| `-mode <conservative|normal|bold>` | Preset for pair thresholds | `normal` |
| `-pp` | Adapter read-through trimming and overlap correction | Off |
| `-ec` | k-mer-spectrum read correction | Off |
| `-mk` | Generate and evaluate multiple k candidates | Off |
| `-nc` | Disable sequence carry-over between k runs | Carry-over enabled |
| `-sr` | Carry synthetic sequences from overlapping pairs | Off |
| `-mg` | Merge sequence from alternative k assemblies | Off |
| `-rv` | Require exact r-mer evidence at both repeat junctions | Off |
| `-la` | Local gap assembly, retrying ambiguous regions at longer k | Off |
| `-my` | Rescue some low-count k-mers between trusted anchors | Off |
| `-po` | Polish substitutions in the final assembly | Off |
| `-cc` | Check reads spanning circular junctions | Off |
| `-gfa` | Write a GFA1 unitig graph | Off |
| `-t <path>` | Output and intermediate directory | `temp` |
| `-rs` | Reuse verified preprocessing/correction outputs | Off |
| `-rt` | Remove recognized intermediates after a successful run | Off |
| `-lang <ja|en|zh>` | Message language | `ja` |
| `-log <quiet|normal|verbose>` | Console verbosity | `normal` |
| `-v` | Show version | — |
| `-h` | Show help | — |

`-mem` controls k-mer counting, not the process's total memory. Graphs, mapping indices, and optional algorithms allocate additional memory. In particular, `-rv` retains exact r-mer keys and may consume substantial memory on noisy reads; it no longer uses a Bloom filter as positive evidence. Multi-k runs and final validation require additional passes over the reads. Merging and aggressive pair thresholds can increase misassemblies; assess correctness as well as contiguity.

## Output

All results are placed under `-t`:

```text
assembly-run/
  assembly.fasta
  unitigs.fasta
  contigs.fasta
  scaffolds.fasta              # when produced
  assembly.report.json
  assembly.provenance.json
  assembly.ambiguous.tsv
  assembly.unsupported.tsv     # when support checking returns a result
  assembly.gfa                 # with -gfa
  Tsumiki.log
  k31/                        # intermediate results per attempted k
  validation/                 # final validation index
```

**Use `assembly.fasta` as the final sequence output.** It contains the selected assembly after short-sequence filtering and any requested polishing. Root-level contig/scaffold files preserve earlier stages and may differ from this final file. The GFA describes the unitig graph, not every subsequent scaffold or polishing edit.

`assembly.report.json` summarizes sequence statistics, unresolved gaps, ambiguity, and available consistency checks. A missing check means unknown, not passed. `assembly.unsupported.tsv` identifies intervals without the required exact read-window support; this is not a reference alignment or a basewise accuracy estimate.

Final read support, polishing, and circular-junction checks use the original uncorrected input reads. Final k-mer consistency is recomputed after sequence edits. These reads come from the same library used for assembly and are **not independent holdout evidence**. The provenance file records their paths and SHA-256 hashes, the final FASTA hash, build identifier, and input settings.

With `-po`, depth is measured again on the polished sequence. With `-rv`, supported windows must occur exactly in the original reads at both junctions; the count is a number of windows, not a number of independent molecules.

## Resuming a run

An existing output directory is rejected unless `-rs` is specified. Re-run with the original input paths and processing options:

```powershell
dotnet publish/Tsumiki.dll -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -pp -ec -t assembly-run -rs -lang en
```

Preprocessing and correction outputs are reused only when their completion records match the build, settings, input contents, and output contents. Older outputs without these records are regenerated. This is stage reuse, not a checkpoint of graph traversal; assembly stages run again. Prefer a new output directory for a different experiment.

The program returns exit code 0 on successful completion or help/version display and 1 when an execution error is caught. A run without an assembly result is an error.

## Algorithm overview

1. Inspect read length and quality encoding, then select an initial k.
2. Optionally preprocess mate overlaps and correct reads against a trusted k-mer spectrum.
3. Count k-mers, select an abundance threshold, simplify the de Bruijn graph, and extract unitigs.
4. Estimate copy number and map read evidence to construct supported contig connections. Resolve some ambiguous branches by bounded look-ahead.
5. Use paired-read links for scaffolding and optionally reassemble unresolved gaps.
6. For multi-k execution, evaluate candidates using a shared anchor k-mer set and select an assembly; optionally attempt merging.
7. Filter short final sequences, optionally polish and verify circular closures, then check the final sequence and write reports.

Copy-number estimates, candidate scores, and bounded search are heuristics. Short reads cannot uniquely resolve every repeat. A circular header or a large N50 alone does not demonstrate correct reconstruction.

For `-la`, a region that remains ambiguous at the assembly k is also examined at k+10 and k+20 where reads are long enough. Conflicting paths remain unresolved. An accepted refinement must have window support from at least two distinct read sequences. Distinct sequences are not equivalent to independent molecules, and this experimental refinement has been tested on small known-answer cases rather than validated genome-wide.

## Development and validation

```powershell
dotnet test Tsumiki.Tests/Tsumiki.Tests.csproj --verbosity minimal -p:GenerateDocumentationFile=true -warnaserror
```

The repository's bundled `read.1.fq` and `read.2.fq` are artificially generated fixtures. They are not suitable evidence for real-data assembly quality or performance. Unit tests with known small examples establish specific invariants, not whole-genome accuracy.

Local evaluation data belong in the ignored `local-data/` directory. For GAGE-B comparisons, record the organism, library, preprocessing, read selection, reference accession and strain, command, build, elapsed time, and peak memory. Confirm that the sample and reference strains match before interpreting discrepancies as assembler errors. Report reference-aligned correctness and genome recovery alongside contiguity.

### Source layout

| Directory / file | Responsibility |
|---|---|
| `Program.cs` | Process entry point |
| `Cores/Pipeline/` | Application orchestration, read preparation, single-/multi-k execution, final validation, checkpoints |
| `Cores/Preprocessing/` | Trimming, correction, k-mer counting, carry-over |
| `Cores/UnitigBuilding/` | Graph construction, simplification, copy-number estimation, look-ahead |
| `Cores/ContigBuilding/` | Mapping evidence and contig traversal |
| `Cores/Scaffolding/` | Scaffold construction and gap assembly |
| `Cores/Evaluation/`, `Cores/Output/` | Selection, validation, and reports |
| `IO/`, `Utilities/`, `Models/`, `Commons/` | Input/output, indexing, models, configuration, messages |
| `Tsumiki.Tests/` | Unit and regression tests |

## License

See [LICENSE.txt](LICENSE.txt) for the MIT license text.
