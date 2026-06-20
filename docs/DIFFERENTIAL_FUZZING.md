# Differential fuzzing harness

The differential fuzzer (`pipeline/difffuzz` + `pipeline/cmd/diff-fuzz`)
stress-tests our ported tools against their vendored upstream originals on
**adversarial inputs**. It is manuscript experiment C2/P1: where the parity
matrix (`pipeline/cmd/parity-pipeline`) checks agreement on *valid* fixtures,
the fuzzer checks agreement on *mutated, malformed, and random* inputs — the
inputs a real user eventually feeds a drop-in replacement by accident.

For each tool/subcommand **target** it:

1. **Generates fuzzed inputs** three ways:
   - **(a) mutation** of a valid seed fixture — bit/byte flips, truncation,
     duplicated/reordered records, and boundary numeric values
     (`0`, `-1`, `2^31`, `2^63-1`, …);
   - **(b) structured random** generation of the target's format (VCF, SAM,
     BED, FASTA, FASTQ) — syntactically plausible documents that parse far
     enough to reach tool *logic*, so divergences here are semantic;
   - **(c) raw random** bytes — to exercise the never-crash / clean-error
     parity path (both binaries should reject the same garbage the same way,
     and neither should segfault).

   All randomness flows from a single `-seed`, so a run is reproducible.

2. **Runs both binaries** — ours (built from `tools/<tool>/cmd/<tool>`) and the
   vendored upstream (resolved from `reference_code/` via
   `pipeline/internal/upstream`) — with identical arguments on each input.

3. **Diffs stdout, stderr, AND exit code**, classifying every divergence:
   `stdout-differs` / `stderr-differs` / `exitcode-differs` / `one-crashed` /
   `both-crashed`. Error-handling parity matters for drop-in behavior, so the
   harness deliberately diffs stderr and the exit code, not just stdout.

4. **Minimizes** a divergence-triggering input by delta-debugging (line-level
   then byte-level `ddmin`, plus a trailing-truncation squeeze) so the report
   carries a small reproducer instead of a 200 KB blob.

5. **Optionally captures Go coverage** of our binary across the run, so the
   report can state what fraction of our code the fuzzing exercised.

6. **Writes `difffuzz.{json,md}`** — per-target input counts, divergence counts
   by class, the minimized reproducers (text and base64), and coverage.

## Provenance / tolerance handling

A "divergence" must mean a *real* behavioral difference, not a benign version
stamp. The classifier therefore normalizes stdout and stderr with the **exact
same** `StripProvenance` the parity harness uses (exported from
`pipeline/runner`): it drops `@PG`/`@CO` SAM headers, `##<tool>_*Command=` /
`##source=` / `##fileDate=` VCF headers, the stats-report version/command/CWD
banner lines, and timing lines — then compares. So an `@PG` line that differs
only by the tool version is **not** flagged; a genuine data line that differs
**is**.

## Running it locally (now)

Everything runs offline against the vendored binaries already under
`reference_code/`. A quick smoke run takes seconds:

```bash
# 3 targets, 40 iterations each — runs in a few seconds
go run ./pipeline/cmd/diff-fuzz -quick

# verbose: log each divergence + the minimization result
go run ./pipeline/cmd/diff-fuzz -quick -v
```

Reports land in `pipeline/.fixtures/<scale>/difffuzz/difffuzz.{json,md}`
(the `.fixtures/` tree is gitignored).

If a target's upstream binary is missing, that target is **SKIPped** (not a
failure) with an actionable reason; the rest of the run proceeds.

### Flags

| Flag | Default | Meaning |
|---|---|---|
| `-quick` | off | few iterations on a 3-target subset (seconds) |
| `-iters N` | 200 (40 in `-quick`) | fuzzed inputs per target |
| `-seed N` | 1 | RNG seed (reproducible) |
| `-targets a,b` | all | comma-separated target-name filter |
| `-coverage` | off | capture Go statement coverage of our binaries |
| `-timeout D` | 10s | per-invocation timeout (a hang → `one-crashed`) |
| `-scale S` | small | seed-fixture scale: `smoke`/`small`/`medium` |
| `-minimize-steps N` | 500 | max delta-debugging steps per reproducer |
| `-out DIR` | `.fixtures/<scale>/difffuzz` | report output dir |
| `-v` | off | verbose progress logging |

Exit status is `1` iff a **strict** parity divergence is found
(`stdout`/`stderr`/`exitcode`/`one-crashed`); `both-crashed` (both reject the
garbage) and SKIPs do **not** fail the run, since both binaries agreeing to
reject an input is not a drop-in mismatch.

### Default targets

`bcftools view`, `bcftools query -f …`, `samtools view`, `samtools flagstat`,
`bedtools merge`, `bedtools intersect` (self-intersection), and `bgzip -d`.
Add or edit them in `pipeline/difffuzz/targets.go`. A target declares its tool,
upstream key, subcommand, argument template (`{in}` → a temp file holding the
fuzzed input; absent → input is piped on stdin), input format, and the seed
fixture key it mutates.

## Coverage capture

With `-coverage`, the harness builds each tool with `go build -cover`, runs it
with `GOCOVERDIR` pointed at a per-run directory so every invocation appends
coverage counters, then renders the percentage via
`go tool covdata percent`. The figure appears in both reports. This needs the
Go ≥ 1.20 binary-coverage toolchain; if it is unavailable the run degrades
gracefully (coverage reported as "not captured", fuzzing unaffected).

To obtain coverage manually for a one-off:

```bash
go build -cover -o /tmp/bcftools ./tools/bcftools/cmd/bcftools
mkdir /tmp/cov
GOCOVERDIR=/tmp/cov /tmp/bcftools view some.vcf >/dev/null
go tool covdata percent -i=/tmp/cov
```

## Scaling it on an external machine

The quick mode is for CI-speed confidence; real fuzzing wants orders of
magnitude more iterations and wall time. On a bigger box:

```bash
# Long campaign: every default target, 100k inputs each, coverage on.
go run ./pipeline/cmd/diff-fuzz -iters=100000 -coverage -scale=medium -v \
  -out=/data/difffuzz-run-$(date +%F)

# Parallelize across seeds (independent reproducible shards) — N shells/nodes:
for s in $(seq 1 32); do
  go run ./pipeline/cmd/diff-fuzz -iters=50000 -seed="$s" \
    -out="/data/difffuzz/shard-$s" &
done
wait
```

Notes for a large external run:

- **Seeds shard cleanly.** Each `-seed` yields an independent, reproducible
  input stream; run one seed per core/node and union the JSON reports. A
  reproducer's `input_b64` field replays the exact bytes that triggered a bug.
- **Pre-build the binaries once.** Our tool binaries are cached under
  `pipeline/.fixtures/bin`; the first target build is amortized across the run.
- **Bound each run with `-timeout`.** Adversarial input can hang a parser; a
  timeout is classified as `one-crashed` (a finding) rather than wedging the
  campaign.
- **Bigger seed fixtures find deeper bugs.** `-scale=medium` (or `large`) seeds
  the mutation strategy with larger valid files, so byte/record mutations reach
  more code paths. Regenerate with the parity pipeline's
  `-update-fixtures` if needed.
- **Coverage is the stopping signal.** Track the reported coverage %; when more
  iterations stop raising it, the fuzzer has saturated the reachable surface
  for the current target/seed-fixture set — time to add targets or richer
  structured generators rather than more iterations.
- **Triage by class.** `stdout-differs` on structured input is the highest-value
  bucket (a semantic data bug); `exitcode-differs` / `stderr-differs` on raw
  input are usually error-message wording differences — real, but lower
  priority for drop-in parity. The Markdown report groups reproducers by class.

## What it has already found

A `-quick` run surfaces real divergences, e.g.:

- **`bcftools view` large-float formatting**: a QUAL of `4294967296` prints as
  `4294967296` (ours) vs `4.29497e+09` (upstream htslib) — htslib switches to
  scientific notation for large floats and our serializer does not.
- **`samtools flagstat` "mate mapped to a different chr"**: on structured SAM
  records with an empty/odd `RNEXT`, our count is `0` where upstream counts the
  mate as mapped elsewhere — a semantic difference in the flag interpretation.
- **`bedtools merge` input validation**: ours accepts (exit 0) several
  malformed inputs upstream rejects (inconsistent field counts; unsorted input
  under the implicit sorted-merge assumption) — validation/robustness gaps.

These are exactly the parity gaps the harness exists to expose; file them
against `docs/PARITY_ROADMAP.md`. (The remaining `stderr-differs` /
`exitcode-differs` on raw-random bytes are mostly error-message wording, which
the manuscript treats as benign for drop-in behavior.)
