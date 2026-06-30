# Local-agent mission prompt

Paste the block below into Claude Code running **on the local testing box** (the
machine that has the GIAB files and the scratch disk). Fill in the four data
paths at the top first. The agent's authoritative instructions live in
[`HEAVY_TIER_RUNBOOK.md`](HEAVY_TIER_RUNBOOK.md); this prompt drives it.

---

```text
You are Claude Code running on a local testing box that has large scratch disk,
network access, and a set of GIAB data files already on disk. Your job is to run
the manuscript's HEAVY-TIER validation for this repository and commit the result
artifacts back. Work autonomously; only stop to ask me if a decision is genuinely
ambiguous or a run fails for a real (non-transient) reason.

DATA ON THIS BOX

This box has 250G free total, use that carefully. 

For fasta reference pull from here: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/. test samtools indexing.

For testing fastq operations use one untrimmed fastq pair from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/ to test adaptor trimming. use fastp before and afer and compare the trimmed version from the repo.

For BAM operations tests use one of bam files from the same repo pointed above. Test samtools merge etc. with another bam file from the same location. after merging using reference-based cram compression.

For vcf files, use https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/ files.

(If any path is wrong or missing, tell me before proceeding. Do NOT subset to a
single chromosome — the whole point is MULTI-CONTIG coverage.)

SOURCE OF TRUTH: read test/manuscript/results/HEAVY_TIER_RUNBOOK.md in full and
follow it. Key rules from it: do NOT install libdeflate-dev (it changes BGZF
bytes and breaks parity); pre-build the upstream binaries serially; copy reports
out of the git-ignored pipeline/.fixtures/ into tracked test/manuscript/results/;
keep raw multi-GB inputs OUT of git.

DO NOT run the GIAB biological-concordance experiment (no truth set / hap.py /
vcfeval). We use the GIAB files purely as real, multi-contig inputs for
our-vs-upstream differential parity, round-trip interop, and performance — the
upstream binary is the oracle.

STEPS (commit after each; push; one draft PR for the whole branch):

0. Setup. From a fresh clone, `git checkout main` (or the latest results branch)
   then `git checkout -b claude/heavy-tier-results`. Build the upstream binaries
   per the runbook (htslib/bcftools/samtools/bedtools, no libdeflate). Run the
   sanity pass: `go run ./pipeline/cmd/full-validation -scales=smoke,small -reps=2`
   and `go run ./pipeline/cmd/realparity` (no args → SKIPs, exit 0). Proceed only
   if both are clean.

1. REAL-DATA parity + performance (multi-contig). Run:
     export TMPDIR=$SCRATCH
     go run ./pipeline/cmd/realparity -ref $REF -bam $BAM -vcf $VCF -reps 5 \
       -out "$PWD/test/manuscript/results/realdata/HG002_GRCh38" -v
   Use the WHOLE-GENOME files (no -region). Gate: every cell PASS, zero DIVERGE
   (the command exits non-zero on a real divergence). Commit report.{json,md}.
   If a cell DIVERGEs, capture the minimized diff and tell me — do not "fix" it
   by relaxing the comparison.

2. LARGE-TIER parity + round-trip interop + performance (full flag matrix, all
   tools, bidirectional interop on multi-contig fixtures, perf). Run:
     go run ./pipeline/cmd/full-validation -scales=medium,large -reps=5 \
       -out "$PWD/test/manuscript/results/large_tier"
   Gate: parity report.md shows 0 DIVERGE and 0 ERROR; roundtrip.md shows every
   interop check PASS. The large tier needs ~30 GB scratch and the slow cells
   (bcftools call, samtools mpileup) take minutes. If it aborts at 100% disk in
   mpileup (a known failure mode), report the stage and stop. Copy the
   report.{json,md}, roundtrip.md, bench.{json,md} into the results dir.

   2a. PERF STATS UPGRADE (required, small code change). parity-bench reduces
   timings to the min over reps; the manuscript needs median + IQR + a ratio CI.
   Enhance pipeline/bench (+ cmd/parity-bench) to record the RAW per-rep wall/CPU/
   RSS samples in bench.json, and add a stdlib-only MedianIQR + bootstrap (or
   Hodges-Lehmann) ratio-CI helper to pipeline/stats (which already has the
   binomial CIs). Re-run with -reps>=10 and write
   test/manuscript/results/large_tier/performance.md: per-cell median ours/upstream
   ratio + IQR + CI, with the slow cells reported PLAINLY (do not hide
   regressions). gofmt/vet/build/test must stay green, incl. `go test -short`.

   2b. HARDWARE SPEC (required). Write test/manuscript/results/large_tier/
   hardware.md with CPU model/cores, RAM, kernel, Go and gcc versions (the
   runbook has the exact snippet).

3. K-RUN REPRODUCIBILITY — ONLY if you have LLM-API budget for it (it re-runs the
   PORTING agent, not the test harness). Follow runbook §4: pick 3-4 tool units,
   K=5-10 runs each from a clean state, record first-pass-compile/parity, turns,
   interventions, tokens/$; write test/manuscript/results/reproducibility/krun.md.
   If you do not have budget, SKIP this and say so — do not fake it.

COMMIT-BACK: per experiment, copy reports out of pipeline/.fixtures/ into the
tracked test/manuscript/results/ tree; run
`npx --yes markdownlint-cli2@0.13.0 "test/manuscript/results/**/*.md"` (0 errors);
`gofmt -l`, `go vet ./...`, `go build ./...`, and `go test -short ./...` must be
clean for any code you touched; commit with a clear message; push to
origin/claude/heavy-tier-results. When done, open a DRAFT PR against main titled
"Heavy-tier validation results (real-data parity + interop + perf)".

HONESTY MANDATE: report slow cells and regressions plainly; if a run fails for a
real reason, report the exact stage and stop rather than guessing or relaxing a
gate. The deliverable is honest numbers, not green-looking ones.
```

---

Notes for the human before pasting:

- Fill `REF`/`BAM`/`VCF`/`SCRATCH` with the actual paths on the box. The BAM/VCF
  must be **indexed** and span **all chromosomes**.
- For long unattended runs, launch inside `tmux`/`screen`; grant the agent
  broader autonomy (accept-edits / skip-permissions on a sandboxed box) so it
  does not block on each command.
- If you have several GIAB samples, run step 1 once per sample (e.g.
  `realdata/HG002_GRCh38`, `realdata/HG005_GRCh38`).
