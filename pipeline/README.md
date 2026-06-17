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
| `seqtk` | `reference_code/seqtk/seqtk` |
| `sickle` | `reference_code/sickle/sickle` |
| `skewer` | `reference_code/skewer/skewer` |
| `fastp` | `reference_code/fastp/fastp` |
| `vcftools` | `reference_code/vcftools/src/cpp/vcftools` |
| `prinseq` | `reference_code/prinseq/prinseq-lite.pl` (run via `perl`) |
| `mosdepth` | `$MOSDEPTH_BIN` or the temp-dir release cache (linux/amd64 only) |

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
ln -sf /path/to/main/reference_code/seqtk/seqtk           reference_code/seqtk/seqtk
ln -sf /path/to/main/reference_code/sickle/sickle         reference_code/sickle/sickle
ln -sf /path/to/main/reference_code/skewer/skewer         reference_code/skewer/skewer
ln -sf /path/to/main/reference_code/fastp/fastp           reference_code/fastp/fastp
mkdir -p reference_code/vcftools/src/cpp reference_code/prinseq
ln -sf /path/to/main/reference_code/vcftools/src/cpp/vcftools reference_code/vcftools/src/cpp/vcftools
ln -sf /path/to/main/reference_code/prinseq/prinseq-lite.pl   reference_code/prinseq/prinseq-lite.pl
```

`sickle`/`skewer` build from source with `make` in their submodule; `prinseq`
is a Perl script (the runner invokes a `*.pl` upstream path through `perl`).
**mosdepth** ships only as a linux/amd64 GitHub release asset (a Nim project, not
built from source): the runner resolves it from `$MOSDEPTH_BIN` or the temp-dir
cache the per-tool mosdepth parity test populates; on other platforms the
mosdepth matrix entries Skip with a clear reason.

A missing binary fails with an actionable hint (the build/submodule command),
never a silent skip.

## Fixtures

Fixtures are **deterministic** (seeded `math/rand`, default seed 1) and
**cross-consistent** (BAM/CRAM/VCF/BED all reference the same contigs). The raw
text (FASTA sequence, SAM records, VCF lines, BED intervals) is generated in Go;
the **vendored upstream tools** then turn it into the valid binary/indexed
formats:

| fixture | placeholder | built by |
|---------|-------------|----------|
| `ref.fa` + `.fai` | `{fasta}` | generated text, indexed by `samtools faidx` |
| `reads.bam` + `.bai`/`.csi` | `{bam}` | SAM text → `samtools sort` → `samtools index` |
| `reads.cram` + `.crai` | `{cram}` | `samtools view -C -T ref` |
| `variants.vcf.gz` + `.tbi` | `{vcf}` | VCF text → `bgzip` → `tabix -p vcf` |
| `variants.vcf` | `{vcf_plain}` | the same VCF text, uncompressed |
| `variants.multi.vcf(.gz)` + `.tbi` | `{vcf_multi_plain}` | multi-sample VCF → `bgzip` → `tabix` (vcftools relatedness/het/LD) |
| `intervals.bed`, `intervals12.bed`, `genome.txt` | `{bed}` `{bed12}` `{genome}` | generated text over the same coordinate space |
| `annotations.gff3` | `{gff}` | gene/mRNA/exon/CDS rows over the same contigs (bed\* / bcftools csq) |
| `reads.fastq` (+ `.gz`) | `{fastq}` `{fastq_gz}` | seeded single-end reads with adapters + low-quality/N tails |
| `reads_R1.fastq`, `reads_R2.fastq` | `{fastq1}` `{fastq2}` | matched-name paired-end reads for PE QC tools |

The **FASTQ** reads have realistic per-read length jitter (±10% of the tier's
mean), high base qualities in the body, the Illumina TruSeq 3' adapter on ~25%
of reads, and a degraded low-quality / N tail on ~30% — so quality-trim,
adapter-trim, and N-handling flags all have real work to do. Both a plain and a
gzip variant are written for the single-end file.

They are cached under **`pipeline/.fixtures/<scale>/`** (gitignored) with a
`manifest.json`; the generator regenerates only when the cache is missing or
stale (manifest version / scale / seed change, or a file went missing). The
`{out}` placeholder is special: it is **not** a fixture but a per-invocation
output prefix the runner assigns to entries that declare `OutputFiles` (see
"Comparison modes").

### Scale tiers

`PIPELINE_SCALE` (or `-scale`) selects the tier; approximate on-disk footprint
of the whole set:

| tier | reference | reads | variants | intervals | FASTQ reads | genes | ≈ size | use |
|------|-----------|-------|----------|-----------|-------------|-------|--------|-----|
| `smoke` | 2 × 20 kb | 2 000 | 400 | 500 | 2 000 | 60 | < 3 MB | CI |
| `small` | 4 × 250 kb | 40 000 | 8 000 | 6 000 | 40 000 | 800 | ~30 MB | default |
| `medium` | 8 × 2 Mb | 300 000 | 60 000 | 40 000 | 300 000 | 6 000 | ~150 MB | benchmarks |
| `large` | 16 × 12 Mb | 2.5 M | 400 000 | 250 000 | 2 M | 40 000 | ~1 GB | heavy benchmarks |

(The FASTQ reads, GFF genes, and multi-sample VCF are generated at every tier
alongside the original BAM/CRAM/VCF/BED set; the multi-sample VCF reuses the
per-tier variant count.)

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

### Tool families currently covered

| family | file | tools | how compared |
|--------|------|-------|--------------|
| htslib core | `smoke.go` | `samtools view`, `bcftools view`/`query` | byte-exact stdout |
| bedtools | `smoke.go` + `bedtools.go` | all 41 bed\* tools (each maps to one `bedtools <sub>`) | byte-exact stdout (text tools) / byte-exact **output files** (`bedsplit`) |
| QC / format | `qc.go` | `seqtk` (23 subcommands), `prinseq`, `sickle`, `skewer`, `fastp` | seqtk byte-exact stdout; prinseq/sickle byte-exact **output files**; skewer byte-exact stdout (per-side args); fastp documented Skips + one Similarity |
| vcftools | `vcftools.go` | `vcftools` (freq/counts/depth/pi/TsTv/het/relatedness/recode/…) | byte-exact **output files** (the `<prefix>.<ext>` the mode writes) |
| mosdepth | `mosdepth.go` | `mosdepth` (per-base, `--by`, thresholds, fast-mode, mapq) | byte-exact **output files** (decompressed `.bed.gz` + summary/dist) |

QC, vcftools, and mosdepth consume the FASTQ / GFF / multi-sample-VCF fixtures
and exercise the output-file comparison path. A handful of entries are
**documented `Skip`s** rather than failures, each recording a *real* divergence
the matrix surfaced (so it is visible without breaking the run):

- **sickle** — our CLI hardcodes window-size 10; upstream uses dynamic
  `int(0.1·len)`. Comparison entries pass `-w 0` (the upstream-faithful
  setting) and are byte-exact; one baseline entry is Skipped to record the bug.
- **vcftools** — `--geno-r2`/`--hap-r2`/`--012` abort with a glibc
  buffer-overflow and `--LROH` segfaults in this upstream build (real upstream
  crashes); `--hardy` differs only by libc `-nan` vs Go `NaN`; `--recode`
  reorders INFO keys (ours alphabetical, upstream source-order) — a real port
  gap owned by the vcftools agent.
- **mosdepth** — `--by` regions depths are byte-exact, but our summary omits
  upstream's per-region rows and we don't emit `region.dist.txt` (gaps owned by
  the mosdepth agent).
- **fastp** — default filtering, the `cut_tail` window boundary, and the
  JSON/HTML report stamps make generic byte/similarity comparison impractical;
  the per-tool suite owns fastp's byte-exact algorithm validation.

The **bedtools** matrix (`bedtools.go`) registers curated combinatorics for all
41 bed\* tools over the shared `{bed}`/`{bed12}`/`{genome}`/`{fasta}`/`{bam}`
fixtures (baseline + per-flag + curated Combos, never the power set). Most are
byte-exact stdout; `bedsplit` uses the output-file path. The RNG tools
(`bedsample`/`bedrandom`/`bedshuffle`) are byte-exact under a fixed `-seed` —
our port reproduces upstream's MT19937 stream exactly. The matrix runs
PASS/SKIP-only; its documented `Skip`s record *real* divergences (each names a
concrete root cause + owner) flagged for follow-up by the bedtools agent:

- **bednuc** — PANICS on every invocation: the `seq` long flag is registered
  twice (`cliflag.BoolVar` then `fs.BoolVar`), so `flag.Var` aborts with "flag
  redefined: seq". The whole tool is unrunnable.
- **interval-sort tie-break** — equal-`(chrom,start)` records are ordered by end
  ascending where upstream uses a stable sort preserving input order. Surfaces
  in `bedsort` (default + `-sizeD`/score keys), `bedmap`/`bedmerge`
  `-o collapse|distinct`, and `bedwindow` B-hit order. Order-independent ops
  (mean/sum/min/max/count/…) are byte-exact and run; `-sizeA`/`-chrThenSizeA`
  sorts are byte-exact.
- **bedwindow (join)** — `-w` default is 0 vs upstream 1000, BED12 `-b` records
  are truncated to 6 columns, plus the tie-break order above; the `-v`/`-c`
  (A-only) outputs are byte-exact.
- **bedfisher** — under-counts overlaps (13134 vs the true 14356 that both
  `bedtools intersect` ports agree on), skewing the contingency table/p-values.
- **bed12tobed6** — drops the BED12 score (emits 0) on each split record.
- **bedexpand** — a trailing comma in the expanded column yields an extra empty
  expansion row.
- **bedmakewindows** — the default `-i none` is rejected by our own parser; `-i
  src` emits an empty name column. `-i winnum`/`-i srcwinnum` are byte-exact.
- **bedsummary** — different output table (missing chrom_length/frac columns)
  and the CLI does not accept upstream's required `-g`.
- **bedannotate** — prepends a `# <file>` header upstream omits and orders
  records differently.
- **bedsubtract** — no reciprocal `-r` flag (upstream has it); other flags pass.
- **bedsplit** — the default `size` heuristic bin-packs differently; `-a simple`
  is byte-exact.
- **bedtobam** — raw BGZF stdout (block framing differs; decoded records are
  identical, per project policy binary BGZF is never byte-compared).
- **bedtag** — different model (upstream tags a BAM/writes BAM; ours is
  BED-in/BED-out), not comparable.
- **bedoverlap / bedunionbedg / bedpairtobed / bedpairtopair** — need a
  pre-joined stream or BEDPE/BedGraph input the fixture corpus does not generate;
  all match out-of-band on crafted inputs and are owned by the per-tool suites.

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
- **Output-file comparison** (the `OutputFiles` mechanism) — for tools that
  write to an output **prefix** instead of stdout (vcftools, mosdepth, and the
  trimming QC tools driven via `-o`). An entry lists the output suffixes it
  produces (e.g. `[".frq"]`, `[".regions.bed.gz", ".mosdepth.summary.txt"]`);
  the runner gives **each side its own temp directory**, resolves the `{out}`
  placeholder to `<tmpdir>/out`, runs both tools, then compares each named
  output file between the two directories. Files ending in `.gz` are
  **decompressed** before comparison (so BGZF block-framing differences are
  irrelevant), and provenance is stripped exactly as for stdout. Each file is
  compared with the entry's `Compare` mode (`ByteExact` or `Similarity`); a
  file present on only one side is a divergence. `CompareOutputFiles` and
  `readMaybeGzip` in `compare.go` implement this; `Compare: DirContents` is the
  conventional label for these entries.
- **Per-side args** (`OurArgs` / `UpstreamArgs`) — for tools whose CLI shape
  genuinely differs from upstream's (our `skewer`/`fastp`/`prinseq`/`sickle`
  are subcommand- or `-i/-o`-based while upstream `skewer` is flat with
  positionals, `prinseq` is the `prinseq-lite.pl` Perl script, etc.). When set,
  these replace the shared `Args` for the respective side; `{placeholder}`
  substitution (including `{out}`) still applies, as does the per-side
  subcommand prepend. A `*.pl` upstream binary path is invoked through `perl`.

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
