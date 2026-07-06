# Heavy-tier validation runbook (local testing box)

This is the **hand-off runbook** for the manuscript experiments that cannot run
in the cloud/CI sandbox because they need a fat node, large scratch disk, bulk
data download, or LLM-API budget. It is written to be executed by an **agent
(Claude Code) running on the local testing box**, but a human can follow it too.

The cheap, in-sandbox Tier-A pieces are already done and committed under
`test/manuscript/results/` (conformance, differential fuzzing, parity CIs,
flag-compat, the transpiler counterfactual). What remains here is the
**resource-heavy** set: real-data (GIAB-file) parity + interop + perf, the
large-tier parity+performance run,
and the K-run reproducibility study. Two further items (a timed human-port
anchor and a second independent bug-corpus labeler) need *people*, not a box,
and are out of scope for this runbook — see "Not in scope" at the end.

## 0. Prerequisites (the box)

| Resource | Minimum | Comfortable |
|---|---|---|
| CPU cores | 8 | 16–32 |
| RAM | 16 GB | 32–64 GB |
| Scratch disk (`$TMPDIR`) | 60 GB | 1 TB (for GIAB-full WGS) |
| Toolchain | Go (version from `go.mod`), `git`, a C/C++ toolchain (`gcc`/`g++`, `make`, `autoconf`, `automake`, `libtool`), `zlib`/`bzip2`/`liblzma`/`libcurl` dev headers | — (biological-concordance tooling — `hap.py`/RTG `vcfeval` — is **out of scope**; not used) |
| Network | outbound HTTPS to GitHub + the NIST GIAB FTP/HTTPS mirror | — |

> **Important — do NOT install `libdeflate-dev`.** htslib built against
> libdeflate emits different BGZF bytes and breaks the byte-exact parity
> goldens. Use the zlib deflate backend (the default).

## 1. Repo setup

```bash
git clone https://github.com/yassineS/bio_ai_experiment.git
cd bio_ai_experiment
# Base off main once PRs #436/#437 are merged; otherwise off the results branch:
git checkout main   # or: git checkout claude/manuscript-tier-a
git checkout -b claude/heavy-tier-results

# Build the upstream binaries once, serially (the harnesses also build on
# demand, but pre-building avoids the concurrent-make race):
git submodule update --init --recursive reference_code/htslib reference_code/bcftools \
  reference_code/samtools reference_code/bedtools
(cd reference_code/htslib && autoreconf -i && ./configure && make -j"$(nproc)")
(cd reference_code/bcftools && make -j"$(nproc)")
(cd reference_code/samtools && make -j"$(nproc)")
(cd reference_code/bedtools && make -j"$(nproc)")

go build ./...    # sanity
```

All result artifacts must be **copied out** of the git-ignored
`pipeline/.fixtures/` tree into a tracked `test/manuscript/results/<experiment>/`
directory and committed (see §6).

## 2. H1 — Large-tier parity + performance run

Runs the full parity matrix, the round-trip suite, and the performance/RSS
bench at the **large** tier. See `docs/FULL_VALIDATION.md` and
`pipeline/bench/README.md`.

```bash
export TMPDIR=/path/to/big/scratch        # ≥30 GB free; fixtures land here
# Full gate at medium+large with 5 bench reps:
go run ./pipeline/cmd/full-validation -scales=medium,large -reps=5 \
  -out="$PWD/test/manuscript/results/large_tier"
```

- **Footprint:** the large manifest (192 Mb ref, 2.5 M reads, 400 k variants,
  250 k intervals) plus the `mpileup`/`call` truth VCFs reaches **~19–30 GB**;
  earlier CI attempts aborted at 100 % disk (`docs/METRICS.md` §6).
- **Slow cells:** `bcftools call`, `samtools mpileup` at large (minutes each);
  `-reps=5` multiplies the bench.
- **Pass gate:** parity `report.md` must show `DIVERGE == 0 && ERROR == 0`
  (SIMILAR is acceptable for the documented floating-point cells).

### H1a — Performance statistics upgrade (required for the manuscript)

`parity-bench` currently reduces wall/CPU to the **min** over reps (and RSS to
the max). The manuscript (claim C3) needs **median + IQR + a CI / Hodges-Lehmann
estimate** of the ours/upstream ratio. Do this:

1. Enhance `pipeline/bench/cmd/parity-bench` (and `pipeline/bench`) to record the
   **raw per-rep samples** (wall, CPU, RSS for each side) in `bench.json`, not
   just the reduced min/max. Keep it stdlib-only.
2. Add a stats helper (reuse/extend `pipeline/stats`) to compute, per cell:
   median, IQR, and a bootstrap (or Hodges-Lehmann) CI of the ratio
   `t_ours / t_upstream`. `pipeline/stats` already has Wilson/Clopper-Pearson;
   add a `MedianIQR` + ratio-CI function with tests.
3. Re-run H1 with `-reps>=10` for stable quantiles and write
   `test/manuscript/results/large_tier/performance.md`: per-cell median ratio +
   IQR + CI, cold vs warm if you capture it, and the **slow cells reported
   plainly** (mpileup, call, isec — do not hide regressions).

### H1b — Hardware spec (required)

Record the box so the numbers are anchored:

```bash
{ echo "## Hardware"; date -u;
  lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Socket|MHz';
  grep MemTotal /proc/meminfo; uname -rsm;
  go version; gcc --version | head -1; } \
  > test/manuscript/results/large_tier/hardware.md
```

## 3. H2 — Real-data parity + round-trip interop + performance (multi-contig)

The GIAB biological-concordance experiment (truth set / `hap.py` / `vcfeval`
variant-calling accuracy) is **out of scope** for this project and is **not**
part of the manuscript — not deferred, not planned. The project validates
byte-exact parity against the upstream oracles, not variant-calling concordance.
Instead the GIAB files are used as **real, whole-genome,
multi-contig inputs** for pure our-vs-upstream **differential parity** and
**performance**, plus **bidirectional round-trip interop**. No truth set is
needed — the upstream binary is the oracle.

> **Multi-contig is required.** Run on the **full multi-chromosome** BAM/VCF, not
> a single-contig subset: cross-contig behaviour (BAI/CSI multi-ref bins, RNEXT
> `=` vs mate-on-another-contig, coordinate sort across contigs, per-contig
> `idxstats`) is exactly what a one-chromosome run misses.

### H2a — Real-data parity + perf (`realparity`)

`pipeline/cmd/realparity` runs the ports against the upstream binaries on a real
reference + BAM/CRAM/VCF and reports, per cell, byte-exact-after-provenance
parity (the repo's exact definition: `runner.StripProvenance` +
`CompareByteExact`) **and** wall/CPU/peak-RSS with ours/upstream ratios. Cells
span samtools (`view`/`flagstat`/`idxstats`/`stats`/`depth`/`quickcheck`/BAM/
`sort`/CRAM) and bcftools (`view`/`norm -f`/`stats`/`query`); cross-contig cells
are flagged `†`.

```bash
go run ./pipeline/cmd/realparity \
  -ref /data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -bam /data/HG002/HG002.GRCh38.bam \
  -vcf /data/HG002/HG002.GRCh38.vcf.gz \
  -reps 5 -out "$PWD/test/manuscript/results/realdata/HG002_GRCh38" -v
```

- Use the **whole-genome** BAM/VCF (do not pass `-region` for the full run; use a
  region only for a quick smoke).
- **Pass gate:** every cell `PASS` (zero `DIVERGE`); the command exits non-zero on
  a genuine divergence. A cell `SKIP`s if its input is absent.
- Commit `report.{json,md}` under `test/manuscript/results/realdata/<sample>/`.

### H2b — Bidirectional round-trip interop + the full flag matrix (`full-validation`)

The full-flag parity matrix (all 53 CLIs) **and** the bidirectional container
interop (`ours-writes/upstream-reads` *and* `upstream-writes/ours-reads` for
BGZF/BAM/CRAM/VCF.gz/BCF/FASTQ + `.bai`/`.csi`/`.tbi` region queries, on
multi-contig fixtures) already run inside `full-validation` (§2 H1). Run it at
the large tier to exercise them at scale; the interop checks use the harness's
own multi-contig fixtures and need no external data. Nothing extra to configure
here beyond H1 — just confirm `report.md` shows `0 DIVERGE/ERROR` and
`roundtrip.md` shows every interop check `PASS`.

## 4. H3 — K-run process-reproducibility (needs LLM-API budget)

Quantifies the nondeterminism of the *porting* process (claim C7 / `04`). This
re-runs the **porting agent**, not the test harness, so it consumes API tokens.

Protocol:

1. Pick a small, representative sample of tools/subcommands (e.g. `sickle`, one
   `bcftools` subcommand, one `bed*` tool) — 3–4 units.
2. For each, run the porting agent **K = 5–10** times from the same clean
   starting state and the same task prompt, each in an isolated worktree.
3. Record per run: first-pass-compile (Y/N), first-pass-parity (Y/N), number of
   agent turns, number of human interventions, wall-clock, tokens/$ (from the
   API usage). Reach-parity = does it pass the same `*Upstream*` parity tests.
4. Report `J/K` reached parity per unit and the intervention distribution in
   `test/manuscript/results/reproducibility/krun.md`. This is the honest
   treatment of agent nondeterminism the red-team asked for.

> Capture the **agent cost** here too (tokens/$, interventions, retries per
> tool) — it is also the data claim C5 needs ("never report agent-cost alone;
> pair it with the supervision cost").

## 5. Sanity check before the long runs

```bash
go run ./pipeline/cmd/full-validation -scales=smoke,small -reps=2   # minutes
go run ./pipeline/cmd/realparity                                    # SKIPs w/o data, exit 0
```

Both should exit 0. If `full-validation` smoke/small is not clean, fix that
before spending hours on large.

## 6. Commit-back protocol

1. Copy every report out of the git-ignored `pipeline/.fixtures/` into the
   tracked results tree, e.g.:

   ```bash
   mkdir -p test/manuscript/results/large_tier
   cp pipeline/.fixtures/large/full-validation/*.{md,json} test/manuscript/results/large_tier/
   ```

2. Keep raw multi-GB inputs OUT of git (only the `.md`/`.json` reports + the
   hardware/spec files belong in the repo).
3. Lint and commit:

   ```bash
   npx --yes markdownlint-cli2@0.13.0 "test/manuscript/results/**/*.md"
   git add test/manuscript/results pipeline/   # pipeline/ only if you changed bench/stats code
   git commit -m "manuscript/results: large-tier perf + real-data parity + K-run"
   git push -u origin claude/heavy-tier-results
   ```

4. Open a PR against `main` (draft) so the numbers flow into the manuscript.

## Not in scope for this box (need people)

- **Timed human-port anchor** (claim C5 ★★): 1–2 expert devs port a scoped slice
  (`sickle` + a few `bcftools` subcommands) to the same parity bar, time-tracked.
- **Second independent labeler** for the bug-corpus κ (claim C7): re-label
  `test/manuscript/bug_corpus.md` independently and compute inter-rater κ.

Both are calibration/validation by humans, not compute — schedule separately.
