# Full validation flow (`pipeline/cmd/full-validation`)

A single orchestrated pass that runs the **entire** validation matrix and the
performance sweep, writes a consolidated report per scale tier, and exits
non-zero if anything diverges — so it can gate a release. It ties together three
stages:

1. **Parity matrix** — every command × its flag-combo cells (bcftools, samtools,
   bedtools, vcftools, the QC tools, htslib utils) is run against the vendored
   **upstream** binary and compared **byte-for-byte** (provenance stripped) or,
   for the few genuinely floating-point outputs, within a numeric **tolerance**
   (similarity mode). Reuses `pipeline/runner` + `pipeline/matrix`.
2. **Round-trip validation** — for each container format the suite implements,
   encode→decode (or decode→re-encode) a fixture and require the payload to come
   back intact, cross-checked against upstream where the format is
   reference-driven. Covers **BGZF, BAM, CRAM, VCF↔BCF, FASTQ** (`pipeline/roundtrip`).
3. **Performance / scalability bench** — wall-clock, CPU (user+sys), and **peak
   RSS** for our binary vs upstream on every load-bearing cell, at the chosen
   tier(s). Reuses `pipeline/bench`.

The orchestrator prints a per-scale **verdict** (PASS only when parity has zero
DIVERGE/ERROR and every round-trip passes) plus the orchestrator's own peak
heap, and writes `report.{json,md}` (parity), `roundtrip.md`, and
`bench.{json,md}` under `pipeline/.fixtures/<scale>/full-validation/`.

## Running it

```bash
# Quick local sanity across the small tiers (CI-friendly, seconds–minutes):
go run ./pipeline/cmd/full-validation -scales=smoke,small -reps=2

# One tool / tier while iterating:
go run ./pipeline/cmd/full-validation -scales=medium -tools=bcftools

# The release / manuscript gate — the LARGE tier (run on a fat node, see below):
go run ./pipeline/cmd/full-validation -scales=large -reps=5

# Skip the bench (parity + round-trip only):
go run ./pipeline/cmd/full-validation -scales=medium -skip-bench
```

Flags: `-scales` (smoke|small|medium|large, comma-separated), `-tools` (parity
filter), `-reps` (bench repetitions), `-seed`, `-update-fixtures`, `-skip-bench`,
`-out`.

## Running the LARGE tier (fat node)

The large tier is **not** runnable on a small-disk / CI box. Requirements:

- **~30 GB+ scratch.** `fixtures.Generate` materialises the full large manifest
  (192 Mb reference, 2.5 M reads, 400 k variants, 250 k intervals) **including
  the `mpileup`/`call` truth VCFs**, which alone reach ~19 GB at large scale.
  (Earlier large-tier attempts on the CI box aborted in `bcftools mpileup` at
  100 % disk — see `docs/METRICS.md` §6.)
- **Upstream binaries built.** The submodules under `reference_code/` must be
  compiled (the harness builds them on demand via `git submodule update --init
  && make`, which needs a C/C++ toolchain). Our binaries are built into
  `pipeline/.fixtures/bin/` automatically.
- Expect a long run: `bcftools call` and `samtools mpileup` at large are the
  slow cells (minutes each); `-reps=5` multiplies the bench.

Suggested invocation on a node with scratch:

```bash
export TMPDIR=/path/to/big/scratch          # fixtures + temp files land here
go run ./pipeline/cmd/full-validation -scales=medium,large -reps=5 \
  -out=/path/to/results
# inspect /path/to/results-style dirs, or the default
# pipeline/.fixtures/{medium,large}/full-validation/
```

## Interpreting the output

- **Parity** `report.md` lists every cell with PASS / SIMILAR / DIVERGE / SKIP /
  ERROR. The gate is **DIVERGE == 0 && ERROR == 0** (SIMILAR is an accepted
  numeric-tolerance match; the one standing SIMILAR is `bcftools call` QUAL — a
  glibc libm last-ULP difference, documented in `docs/METRICS.md` §5).
- **Round-trip** `roundtrip.md` lists each format check PASS / FAIL / SKIP.
- **Bench** `bench.md` is the wall/CPU/RSS table (ratios = ours/upstream).

## Known finding surfaced by this flow

The round-trip stage currently reports **one FAIL: `bam-via-cram` (CRAM)**. It is
a **pre-existing** CRAM *encoder* edge case (reproduces on the pre-`#428` base),
not a regression — and it is exactly the kind of bug single-operation parity
missed: the medium tier happens to encode byte-identically to upstream, so
`view BAM→CRAM` passes there, but on the smoke fixture our `view -C` produces a
CRAM whose SEQ cannot be reconstructed (our decode → `NNNN`, upstream decode →
`====`). Isolation (2×2 encode/decode matrix) shows our **decoder is correct**
(it decodes upstream's CRAM byte-for-byte) — the fault is in our CRAM **encoder**
on small inputs. Tracked for a dedicated fix; the round-trip check is written to
compare our round-trip against upstream's (the correct oracle for a
reference-driven format), so it will pass automatically once the encoder is
fixed.
