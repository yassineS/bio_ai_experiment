# Parity & performance pipeline

This is the **integration / combinatorics / performance layer** that sits on
top of the per-tool `*_test.go` parity suites. Where those tests assert parity
for one tool on small crafted inputs, this pipeline:

1. generates **real-sized, valid, cross-consistent** fixtures (FASTA, BAM,
   CRAM, VCF, BED) at a configurable scale,
2. runs a **combinatorics matrix** of command/flag combinations through **our**
   tool binary and the **vendored upstream** binary,
3. **compares** their output (byte-exact, similarity, or directory contents),
   capturing wall-clock timing on every run, and
4. emits a machine-readable `report.json` and a human `report.md`.

It is the keystone framework; the per-tool matrices and the full bench set are
populated incrementally (see [Extending](#extending)).

## Quick start

```bash
# smoke (CI-sized, runs in seconds) — proves the loop end to end
go run ./pipeline/cmd/parity-pipeline -scale=smoke

# default scale (small, ~5 MB), filtered to two tools, verbose
go run ./pipeline/cmd/parity-pipeline -scale=small -tools=samtools,bcftools -v

# force fixture regeneration
go run ./pipeline/cmd/parity-pipeline -scale=medium -update-fixtures

# performance benchmarks (separate harness; see Benchmarks below)
PIPELINE_SCALE=medium go test -bench=. ./pipeline/bench
```

Exit status is **non-zero** if any entry **DIVERGE**s or **ERROR**s. `SIMILAR`
matches (heuristic paths within tolerance) and `SKIP`s do **not** fail the run.

### Driver flags

| flag | meaning |
|------|---------|
| `-scale` | `smoke` \| `small` \| `medium` \| `large` (default `$PIPELINE_SCALE` or `small`) |
| `-tools` | comma-separated tool filter (default: all registered) |
| `-bench` | print the bench-harness invocation hint |
| `-update-fixtures` | regenerate fixtures even if a valid cache exists |
| `-out` | report output dir (default `pipeline/.fixtures/<scale>/report`) |
| `-seed` | RNG seed for fixture generation (default 1) |
| `-v` | verbose logging |

## Upstream binaries

The pipeline locates the **vendored upstream binaries** under `reference_code/`:

| key | path |
|-----|------|
| `samtools` | `reference_code/samtools/samtools` |
| `bcftools` | `reference_code/bcftools/bcftools` |
| `bgzip` / `tabix` | `reference_code/htslib/{bgzip,tabix}` |
| `bedtools` | `reference_code/bedtools/bin/bedtools` |

These are expected to already be **built** (exactly like the existing per-tool
live parity tests assume). In an **isolated git worktree** the submodules may be
unpopulated; symlink the binaries from the main checkout — they are **not**
committed (they live under submodule paths):

```bash
ln -sf /path/to/main/reference_code/samtools/samtools reference_code/samtools/samtools
ln -sf /path/to/main/reference_code/bcftools/bcftools reference_code/bcftools/bcftools
ln -sf /path/to/main/reference_code/htslib/bgzip       reference_code/htslib/bgzip
ln -sf /path/to/main/reference_code/htslib/tabix       reference_code/htslib/tabix
mkdir -p reference_code/bedtools/bin
ln -sf /path/to/main/reference_code/bedtools/bin/bedtools reference_code/bedtools/bin/bedtools
```

A missing binary fails with an actionable hint (the build/submodule command),
never a silent skip.

## Fixtures

Fixtures are **deterministic** (seeded `math/rand`, default seed 1) and
**cross-consistent** (BAM/CRAM/VCF/BED all reference the same contigs). The raw
text (FASTA sequence, SAM records, VCF lines, BED intervals) is generated in Go;
the **vendored upstream tools** then turn it into the valid binary/indexed
formats:

| fixture | built by |
|---------|----------|
| `ref.fa` + `.fai` | generated text, indexed by `samtools faidx` |
| `reads.bam` + `.bai`/`.csi` | SAM text → `samtools sort` → `samtools index` |
| `reads.cram` + `.crai` | `samtools view -C -T ref` |
| `variants.vcf.gz` + `.tbi` | VCF text → `bgzip` → `tabix -p vcf` |
| `variants.vcf` | the same VCF text, uncompressed |
| `intervals.bed`, `intervals12.bed`, `genome.txt` | generated text over the same coordinate space |

They are cached under **`pipeline/.fixtures/<scale>/`** (gitignored) with a
`manifest.json`; the generator regenerates only when the cache is missing or
stale (manifest version / scale / seed change, or a file went missing).

### Scale tiers

`PIPELINE_SCALE` (or `-scale`) selects the tier; approximate on-disk footprint
of the whole set:

| tier | reference | reads | variants | intervals | ≈ size | use |
|------|-----------|-------|----------|-----------|--------|-----|
| `smoke` | 2 × 20 kb | 2 000 | 400 | 500 | a few hundred KB | CI |
| `small` | 4 × 250 kb | 40 000 | 8 000 | 6 000 | ~5 MB | default |
| `medium` | 8 × 2 Mb | 300 000 | 60 000 | 40 000 | ~50 MB | benchmarks |
| `large` | 16 × 12 Mb | 2.5 M | 400 000 | 250 000 | ~500 MB | heavy benchmarks |

## The matrix

`pipeline/matrix` defines the declarative model. An `Entry` is one runnable
parity case:

```go
type Entry struct {
    Tool, Subcommand, UpstreamTool string
    UsesSubcommand                 bool
    Name                           string
    Args                           []string   // with {bam} {cram} {vcf} {vcf_plain} {bed} {bed12} {fasta} {genome} placeholders
    Input                          InputKind
    Compare                        CompareMode // ByteExact | Similarity | DirContents
    Heavy                          bool
    Skip                           string
}
```

Tools register entries into a global `Registry` from `init()` (`matrix.Register`).
The driver consumes `matrix.Default()`.

### Combinatorics expander (curated, NOT power set)

`ExpandSpec.Expand()` turns a subcommand's flag set + value choices into:

1. a **baseline** entry (no extra flags),
2. **one entry per single flag value** (every flag exercised in isolation),
3. exactly the **curated multi-flag `Combos`** the author lists.

We **deliberately do not** generate the `2^N` power set: for tools like
`samtools view` (dozens of flags) it is intractable, mostly meaningless or
mutually exclusive, and unreadable. Size stays **linear-plus-curated**
(`N + |Combos|`). Adding a specific interaction is one line — append a `Combo`.
For the rare *bounded* case where a small full cross-product really is
meaningful, `matrix.CrossProduct(flags...)` builds it explicitly (and you pass
the result as `Combos`), keeping the no-power-set policy visible at the call
site.

## Comparison modes

Implemented in `pipeline/runner/compare.go`:

- **`ByteExact`** (default) — provenance-stripped stdout equality. The strip
  helper drops SAM `@PG`/`@CO` lines and VCF `##...Command`/`##source`/
  `##fileDate`/version headers so tool-version stamps don't cause false
  divergence; data lines and structural headers are preserved.
- **`Similarity`** — for documented heuristic / unseeded-RNG / libm-float paths:
  structural comparison with numeric fields within a relative epsilon (`1e-6`),
  recording the **max deviation** observed. Use this where byte-exact is not
  expected but structural + numeric agreement is.
- **`DirContents`** — placeholder for multi-file outputs; currently falls back
  to byte-exact of stdout. (Follow-on: compare the set + stripped contents of an
  output directory.)

Binary outputs (BGZF BAM/CRAM bytes) are **not** compared byte-exact — our
klauspost deflate backend frames blocks differently from htslib though both
decode identically. Compare the **decoded** stream (e.g. `samtools view` to SAM)
for parity, and use the bench harness for the encode path's timing.

## Report

Both written under `-out` (default `pipeline/.fixtures/<scale>/report/`):

- **`report.json`** — full machine-readable results + summary tally.
- **`report.md`** — per-tool tables of combos run, PASS/SIMILAR/DIVERGE/SKIP/
  ERROR counts, ours/upstream timing, and the `ours/upstream` **ratio** for
  `Heavy` entries.

## Benchmarks

`pipeline/bench` is a `go test -bench` harness that times **our** binary against
the **upstream** binary on a `medium`/`large` fixture and reports the ratio:

```bash
PIPELINE_SCALE=medium go test -bench=. ./pipeline/bench
PIPELINE_SCALE=large  go test -bench=. -benchtime=3x ./pipeline/bench
```

Each benchmark reports custom metrics `ours_ms/op`, `upstream_ms/op`, and
`ratio_ours/up`. Three representative benches ship as the pattern
(`samtools view` BAM→BAM, `bcftools view` filter, `bedtools intersect`).

## Extending

**Add a matrix entry** — in `pipeline/matrix`, either append to a smoke spec or
add a new file with an `init()` that calls `matrix.Register(...)`. Prefer
`ExpandSpec` for flag sweeps; add `Combos` for interactions. Use `{...}`
placeholders for fixture paths.

**Add a benchmark** — in `pipeline/bench`, add a `BenchmarkXxx` that calls
`benchPair(b, ourTool, upstreamKey, ourArgs, upstreamArgs)`; reuse `benchManifest`
for the fixture.

**Add a fixture** — extend the writers in `pipeline/fixtures` and record the new
file in the manifest; bump `manifestVersion` so caches invalidate.

## Relationship to the per-tool parity suites

The `tools/<tool>/.../*_test.go` suites remain the **authoritative** parity
checks (crafted edge cases, upstream test corpora, error-path stderr). This
pipeline is the layer **on top**: it exercises real-sized inputs, the
combinatorics of flags, and performance — catching scale/interaction/perf
regressions that small unit fixtures don't.
