# Parity & performance pipeline

This is the **integration / real-data / performance layer** that sits on top of
the per-tool `*_test.go` parity suites. Where those tests assert parity for one
tool on small crafted inputs, this layer runs every tool over **real data**
against the **vendored upstream binaries**, measuring both byte-for-byte parity
and performance, and keeps a battery of **crafted-input correctness** suites
(differential fuzzing, edge cases, container round-trip, conformance corpora).

> The earlier **synthetic** combinatorics harness (the `PIPELINE_SCALE`-tiered
> `parity-pipeline` / `full-validation` drivers and the `parity-bench`
> micro-benchmarks under `pipeline/bench`) has been **retired** and replaced by
> the real-data `realbench` / `realparity` harnesses below.

## Layout

### Real-data parity + performance

- **`realbench`** (`pipeline/realbench`, CLI `pipeline/cmd/realbench`) — the
  full-port real-data benchmark **and** parity harness. It drives our ports
  against the upstream binaries over the **whole ported surface** (every
  subcommand of samtools, bcftools, the ~41 bed\* tools, seqtk, fastp, sickle,
  skewer, prinseq, vcftools, mosdepth, bgzip, tabix, htsfile) on a real GIAB
  HG002 / GRCh38 dataset, at three tiers (`chr20` | `exome` | `wgs`). Each cell
  reports **parity** (byte-exact after provenance-stripping; BGZF/BAM/CRAM
  decoded first, never compared as raw framing) as PASS / DIFF / SKIP / ERROR,
  plus **performance** (wall / CPU / peak RSS for each side and the
  ours/upstream ratios). Its **Seqera / Nextflow entrypoint** is `test/nextflow/` —
  see [`../test/nextflow/README.md`](../test/nextflow/README.md) for how the real inputs
  are staged and region-subset (with the **upstream** tools, so the per-tier
  inputs never come from the code under test) and how to run it on the Seqera
  Platform / AWS Batch or locally under Docker.
- **`realparity`** (`pipeline/cmd/realparity`) — the companion real-data
  samtools/bcftools differential parity + performance runner over whole-genome,
  **multi-contig** inputs (a GIAB-class BAM/CRAM/VCF + indexed reference). Pure
  our-vs-upstream differential testing: no truth set is needed. Every input is
  optional; a cell whose input is absent SKIPs, so it exits cleanly (all-SKIP)
  on a machine with no real data.
- **`giab`** (`pipeline/giab`, CLI `pipeline/cmd/giab-concordance`) — the GIAB
  real-data biological-**concordance** harness (against the GIAB benchmark VCF +
  high-confidence BED). See
  [`../docs/GIAB_CONCORDANCE.md`](../docs/GIAB_CONCORDANCE.md).

### Crafted-input correctness suites (unchanged)

- **`difffuzz`** (`pipeline/difffuzz`, CLI `pipeline/cmd/diff-fuzz`) — the
  differential fuzzer that cross-checks our stdout/stderr/exit against upstream
  on randomised inputs. See
  [`../docs/DIFFERENTIAL_FUZZING.md`](../docs/DIFFERENTIAL_FUZZING.md).
- **`edgecases`** — a battery of crafted edge-case parity tests (CRAM
  references, calmd MD/NM, bcftools norm, etc.).
- **`roundtrip`** — container encode→decode (and the bidirectional
  ours-writes/upstream-reads **and** upstream-writes/ours-reads **interop**)
  round-trip checks for BGZF, BAM, CRAM, VCF↔BCF, FASTQ + `.bai`/`.csi`/`.tbi`.
- **`conformance`** — the htslib `test/` + htscodecs conformance corpora run
  through our binaries.

### Shared machinery

- **`fixtures`** — deterministic, cross-consistent fixture writers
  (FASTA/BAM/CRAM/VCF/BED/FASTQ/GFF) used by the crafted-input suites.
- **`matrix`** — the declarative `Entry`/`Registry` model (`matrix.Register`,
  `ExpandSpec`/`Combos` flag-sweep) for parity cells.
- **`runner`** — the comparison core: `StripProvenance` (drops SAM `@PG`/`@CO`
  and VCF version/command/date headers so tool-version stamps don't cause false
  divergence), `CompareByteExact`, gzip-aware output-file comparison, and the
  report writers. Both `realbench` and `realparity` reuse it, so "parity" means
  the same thing everywhere in the repo.
- **`stats`** — Wilson + Clopper-Pearson binomial confidence intervals
  (stdlib-only) for the manuscript parity-rate statistics.

## Comparison semantics

Parity is **byte-exact equality after provenance-stripping**. Binary outputs
(BGZF BAM/CRAM bytes) are **not** compared as raw framing — our klauspost
deflate backend frames blocks differently from htslib though both decode
identically — so the harness compares the **decoded** stream (e.g. `samtools
view` to SAM, or a `gunzip` of a `.gz` output) instead. A handful of genuinely
floating-point or unseeded-RNG outputs are compared structurally within a
numeric tolerance. `runner.StripProvenance` / `runner.CompareByteExact`
implement this and are shared by every layer above.

## Relationship to the per-tool parity suites

The `tools/<tool>/.../*_test.go` suites remain the **authoritative** parity
checks (crafted edge cases, upstream test corpora, error-path stderr). This
pipeline is the layer **on top**: `realbench`/`realparity` exercise real-sized
inputs and performance against upstream, while the crafted-input suites
(`difffuzz`/`edgecases`/`roundtrip`/`conformance`) catch
interaction/round-trip/conformance regressions that small unit fixtures don't.
