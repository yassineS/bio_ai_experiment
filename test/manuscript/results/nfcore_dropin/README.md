# nf-core drop-in demo — `samtools/flagstat`

**Our Go `samtools` is a drop-in for the nf-core `samtools/flagstat` module:
the module ran unchanged, and its output is byte-identical to upstream.**

## What this shows (C4 usability)

The **real, unmodified** nf-core `samtools/flagstat` module is vendored verbatim
under [`modules/nf-core/samtools/flagstat/main.nf`](modules/nf-core/samtools/flagstat/main.nf)
(fetched from `nf-core/modules` `master`; `main.nf`, `meta.yml` and
`environment.yml` are the upstream files, untouched). A minimal DSL2 wrapper
([`main.nf`](main.nf)) instantiates that process on a small 40,000-read BAM. The
module was run **twice on the same input**, changing nothing but which
`samtools` is first on `PATH`:

1. upstream htslib `samtools 1.22.1` (`bin/upstream/samtools`)
2. our Go re-implementation (`bin/ours/samtools`)

Both runs completed successfully and produced **byte-identical** `.flagstat`
output. The nf-core module file is never edited between the two runs — only
`PATH` changes — which is exactly the drop-in claim.

## How to reproduce

```bash
make ours upstream                       # build both binaries into bin/
bash test/manuscript/results/nfcore_dropin/run_dropin.sh
```

Requirements: `nextflow` on `PATH` (tested with Nextflow 26.04.4, DSL2).
Containers/conda are disabled in [`nextflow.config`](nextflow.config) so the
local executor resolves `samtools` from `PATH`.

## Evidence (captured run logs & outputs)

- [`evidence/upstream.stdout`](evidence/upstream.stdout) — Nextflow run with
  upstream samtools; shows `[SUCCESS] completed=1 failed=0` and which binary ran.
- [`evidence/ours.stdout`](evidence/ours.stdout) — same, with our samtools.
- [`evidence/upstream.flagstat`](evidence/upstream.flagstat),
  [`evidence/ours.flagstat`](evidence/ours.flagstat) — the module's emitted
  output files; `diff` reports **no differences**.

Captured run (abridged):

```
=== [upstream] samtools = .../bin/upstream/samtools ===
=== [upstream] samtools version line: samtools 1.22.1-19-ge406d9e ===
[PROCESS 37/9530e8] SAMTOOLS_FLAGSTAT (test)
[SUCCESS] completed=1 failed=0 cached=0

=== [ours] samtools = .../bin/ours/samtools ===
=== [ours] samtools version line: 0.1.0 ===
[PROCESS 6c/5dc383] SAMTOOLS_FLAGSTAT (test)
[SUCCESS] completed=1 failed=0 cached=0

=== DIFF upstream vs ours (.flagstat) ===
RESULT: BYTE-IDENTICAL — our samtools is a drop-in for nf-core samtools/flagstat
```

The only observable difference between the two runs is the version string the
module reports (`samtools version | sed ...`): upstream emits its htslib
version, ours emits the port's `0.1.0`. The scientific output — the `.flagstat`
report the module exists to produce — is identical.

## Files

```
nfcore_dropin/
├── README.md                                  # this file
├── main.nf                                     # thin DSL2 wrapper (the only added Nextflow)
├── nextflow.config                             # local executor, no container/conda
├── run_dropin.sh                               # runs the module twice, diffs output
├── data/test.bam (+ .bai)                      # 40k-read fixture (from pipeline/.fixtures/small)
├── modules/nf-core/samtools/flagstat/          # the REAL nf-core module, vendored verbatim
│   ├── main.nf
│   ├── meta.yml
│   └── environment.yml
└── evidence/                                   # captured run logs + output .flagstat files
    ├── upstream.stdout / upstream.flagstat
    └── ours.stdout / ours.flagstat
```
