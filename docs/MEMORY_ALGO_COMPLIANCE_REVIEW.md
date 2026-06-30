# Memory & algorithmic compliance review (2026-06-30)

A parallel, per-tool/subcommand review of whether the Go port's **memory
complexity** and **math/algorithm** match the vendored upstream. Triggered by the
memory figure showing most tools at ~2× and several at >5× RSS.

## Headline answer

**Most of the high ratios are benign, not algorithmic errors.** The flagged
tools (sickle 11.5×, seqtk seq 5.9×, bed* 4–6×) were measured on **1–3 MB
synthetic inputs**, where Go's fixed runtime floor (~8–15 MB: heap arenas, GC
headroom, goroutine stacks, buffered I/O) dwarfs upstream's tiny C footprint. A
2–5× ratio there is the floor, not over-allocation — the same Go-vs-C runtime
floor characterised for the CRAM `view` residual (live heap leaner than
upstream's whole RSS). A finding only counts as **real** when our extra memory
**grows with input size** faster than upstream's algorithm requires.

By that test, the review surfaced **one genuine memory over-allocation**, plus a
correctness bug and a process lesson:

## Confirmed real findings

| # | Area | Severity | What | Status |
|---|------|----------|------|--------|
| 1 | CRAM `--threads` parallel decode (`pkg/htsgo/cram/parallel.go`) | **medium** | `containerJobBuffer=32` across 3 channels → up to **96 fully-decoded containers in flight**, each holding all its slices' `[]*sam.Record`. For WGS CRAM that is ~10–100× the single-threaded working set. htslib (`cram_flush_container_mt`) keeps only ~2–4 containers/thread and recycles immediately. **Scales with input.** | Open — fix candidate |
| 2 | `samtools consensus` bayesian (`consensus_bayesian.go`) | **high (crash)** | `81f87cc` removed `qual>100`/`qual2>100` clamps; a non-conformant BAM (base qual >100; valid SAM is Phred 0–93) panics with index-out-of-range on the `[101]` prob tables. | **Fixed** (`b6efe89`, clamps restored; byte-inert for valid SAM) |

### Finding 1 — CRAM parallel decode buffering (fix candidate)

The single-threaded CRAM `view -b -T` path is now ≤2× wall-neutral (the prior
work). The **`-@`/threaded** path, however, buffers up to `3 × containerJobBuffer
= 96` decoded containers. Recommended fix: cap `containerJobBuffer` near
`threads+2` and add back-pressure between `collect()` and `fillNextSliceParallel`
so consumed containers are freed before more decode. Must stay byte-identical
(records emitted in strict file order) and `-race` clean.

## Process lesson (caught, not shipped)

The review fan-out used agents with write access, and — before being reverted —
some agents **edited parity-critical code** mid-review. The adversarial synthesis
then caught that one of those edits was a **real high-severity regression**:

- **`vcftools --hap-r2` D-statistic** was rewritten from `D = pAB - pA*pB`
  (float64, the upstream `rel_x11 - p1*q1` form) to an integer-numerator form,
  justified by an **x87 80-bit extended-precision** argument. An exhaustive
  comparison (635,371 count configs, n=2..60) showed the **original matched
  upstream 100%; the new form only 29.78%**. The x87 premise is false — both
  oracle platforms (vendored aarch64 Linux, CI x86-64 SSE2) use IEEE-754 double.
  This was **reverted** (never committed).

Anti-patterns worth enforcing repo-wide:

1. **No FP-parity claim should cite "x87/80-bit extended precision" without
   verifying the oracle arch actually uses it.** Both oracles here are strict
   IEEE-754 double; x87 has not been default on either arch for ~20 years.
2. **Byte-exact statistic changes must be gated by a live-binary parity test,
   not a self-asserting unit test** that pins what the author *believes* upstream
   does (`TestComputeHapR2_IntegerZeroD` pinned the regressed value and would
   have passed CI green).

## Coverage caveat — this review is PARTIAL

The run was contaminated: because review agents dirtied the working tree, many
later agents reviewed *that uncommitted diff* instead of cleanly assessing their
assigned unit's memory profile (the confirmed-finding tool/area labels are
cross-wired, e.g. `bedjaccard` → `consensus_bayesian.go`). So the only
**clean, verified** memory result is Finding 1 (CRAM parallel path); the broad
per-tool memory sweep did **not** reliably complete.

Known load-all/per-position patterns from `docs/METRICS.md` (`bcftools isec`
×30, `samtools stats` ×26, `samtools depth` per-position arrays, `bcftools norm`
×17) were **not** cleanly re-assessed here and remain candidates.

## Recommendation

Re-run a **clean, read-only** review (Explore agents, no write access, on a clean
tree) to comprehensively answer "where does our memory genuinely scale worse than
upstream?" — focused on the htslib-core load-all paths above. Then fix the
confirmed real ones (starting with the CRAM parallel buffer) tool-by-tool with
byte-exact + `-race` gates.
