# Project metrics

Quantitative summary of the `bio_ai_experiment` port, for the manuscript. All
numbers are reproducible from the repository (the command used is given beside
each); regenerate with the snippets below after material changes.

> Honesty note: this file states scope, wins, **and** the current
> shortcomings (the performance hotspots and the one documented divergence). It
> is intended to be cited as-is, so nothing here is rounded up.

## 1. Scope (what "the port" is)

The project does **not** attempt to port every bioinformatics tool. Its scope is
a deliberately bounded set: the **htslib ecosystem core**, the **bedtools
interval-analysis surface**, and a **QC / format tool set**. Within that scope,
subcommand coverage is complete.

| Layer | Tools | Subcommands / binaries |
|---|---|---|
| htslib core | samtools, bcftools, bgzip, tabix, htsfile | samtools 24 subcommands, bcftools 24 subcommands |
| bedtools surface | ~37 `bed*` tools (bedintersect, bedmap, …) | one drop-in binary each |
| QC / format | seqtk, prinseq, sickle, skewer, fastp, vcftools, mosdepth | per-tool CLIs |

- **13 in-scope tool families**, shipped as **53 drop-in POSIX CLI binaries**
  (`git ls-files 'tools/*/cmd/*/main.go' | … | sort -u | wc -l`).
- The repository vendors **53 upstream submodules** under `reference_code/`, but
  ~40 of those (bwa, STAR, salmon, …) are **reference material from the phase-1
  tool survey** (`analysis/`, the top-~200 ranking), *not* in-scope-but-unported
  tools. They are vendored for provenance, not slated for porting (see
  `CLAUDE.md`: "no longer taking on new tools").

So the correct manuscript claim is **"complete reimplementation of a bounded,
named tool set,"** not "a subset of an open-ended porting effort."

## 2. Lines of code

| | LOC | how |
|---|---|---|
| Upstream source of the ported tools (C/C++/Perl: htslib, htscodecs, samtools, bcftools, bedtools, seqtk, sickle, skewer, fastp, prinseq, vcftools) | **~441,300** | `find … -name '*.c' -o … \| wc -l` |
| Our Go — non-test | **218,116** | `git ls-files '*.go' \| grep -v _test \| xargs wc -l` |
| Our Go — tests | 142,575 | `git ls-files '*_test.go' \| xargs wc -l` |
| Our Go — total | 360,691 | |

> Caveat: the upstream figure is the *full* source of those tools (including
> paths we do not exercise and build scaffolding; htslib alone is ~120k LOC of a
> shared library we only partially reimplement in `pkg/htsgo`). It is a
> scope-honest "size of the corpus," **not** a 1:1 reimplementation ratio.

### Retained code by phase (git-blame attribution on the final tree)

The project ran in two phases (see §3). Attributing every line in the *current*
tree to the phase whose commit last touched it:

| | non-test Go LOC | share |
|---|---|---|
| Phase 1 (retained) | 11,291 | **5.2 %** |
| Phase 2 (retained) | 206,825 | **94.8 %** |

Phase 2 authored ~95 % of the surviving production code — it was effectively a
rewrite, not an incremental extension of phase 1.
(`git blame --line-porcelain` per file, bucketed by committer-year; boundary
2026-01-01.)

## 3. Process: phases, commits, PRs

The two phases are unambiguous in history. Phase 1 used `copilot/*` branches
(Oct–Dec 2025); phase 2 used `claude/*` branches (May–Jun 2026) and opened with
the `CLAUDE.md` guidance file.

| metric | Phase 1 | Phase 2 | total |
|---|---|---|---|
| Commits | 123 | 709 | **832** |
| Merged PRs | 21 (#2–22) | 335 (#23–413) | **356** |
| Bug-fix commits (subject matches fix/bug/correct/…) | 5 | 88 | 93 |

> "Phase-1 bugs fixed in phase-2" is **not** cleanly countable: with 95 % of
> phase-1 code rewritten, those defects were *replaced*, not patched. The first
> phase-2 commits ("implement real stats, reject the rest, fix overstated docs")
> indicate phase 1 carried stubbed / overstated implementations that phase 2
> discarded. We report the rewrite fraction (above) rather than a fabricated
> bug-fix count.

### Upstream interaction

- **Upstream bugs identified & documented:** **21** (`docs/UPSTREAM_BUGS.md`),
  a mix of fix-on-port (our output is corrected) and reproduced-for-parity.
- **Upstream PRs filed:** 0 filed, **1 staged ready-to-file**
  (`docs/upstream-reports/bcftools-csq-prnup-recycled-pos.md`) — this
  environment's GitHub scope cannot open issues on `samtools/bcftools`.

## 4. Documentation

| metric | value | how |
|---|---|---|
| Markdown files (excl. vendored) | 108 | `git ls-files '*.md' \| grep -v reference_code \| wc -l` |
| Words of documentation | ~195,400 | `… \| xargs wc -w` |
| Per-tool `README.md` | 51 | `git ls-files 'tools/*/README.md' \| wc -l` |
| `docs/` reference documents | 18 | |
| Doc lines authored in phase 2 (retained) | 20,879 of 34,001 (**61 %**) | blame attribution |

Phase 2 also created the status/parity documentation system from scratch:
`PROJECT_STATUS.md`, `docs/PARITY_ROADMAP.md`, `docs/UPSTREAM_BUGS.md`, the docs
map, and the per-tool READMEs.

## 5. Correctness (parity)

Differential testing harness: `pipeline/cmd/parity-pipeline` runs OUR binary and
the vendored UPSTREAM binary on shared, seeded fixtures and compares output
byte-for-byte (compressed BAM/BGZF compared after decode, since framing
legitimately differs).

Small tier, **400 invocations across 53 tools**:

| status | count | meaning |
|---|---|---|
| PASS (byte-exact / decoded byte-exact) | **398** | identical output |
| SIMILAR (within tolerance) | 1 | `bcftools call` QUAL — a glibc-version-dependent libm last-ULP rounding (every other field byte-exact; values agree to ~6 sig figs); accepted at 2e-5 |
| SKIP (documented upstream bug) | 1 | `bcftools csq` `@pos` marker — an upstream recycled-pointer bug; **our output is the correct one** (`docs/UPSTREAM_BUGS.md`) |
| DIVERGE | **0** | — |

Backed by **~4,700 Go test functions** (unit + parity).

> Rigor caveat for the paper: this is empirical byte-equivalence over a curated
> matrix of invocations on **seeded synthetic** fixtures, on a single platform
> (linux/amd64) and a single upstream version each. It is not a formal proof,
> nor (yet) validated on real public data at production scale. Subcommand
> coverage is complete; exhaustive *flag-combination* coverage is not formally
> measured.

## 6. Performance & scalability

Measured by `pipeline/bench` (`cmd/parity-bench`): wall-clock, CPU (user+sys)
and peak RSS for each side, from `wait4` `rusage`. Ratio = ours / upstream
(**< 1.0 = we use less**).

> The `smoke`/`small` tiers are dominated by Go process startup + GC, so their
> ratios overstate our cost. The **medium** tier (16 Mb reference, 300 k reads,
> 60 k variants) is the first tier where steady-state behaviour shows; **large**
> is pending (to be run once the hotspots below are addressed).

### Medium tier — selected cells (wall ratio, ours/upstream)

Wins / parity:

| operation | wall× |
|---|---|
| sickle se (FASTQ trim) | **0.51** |
| samtools view BAM→BAM | **0.82** |
| samtools sort | **0.88** |
| bed sort | 1.02 |
| bedtools coverage / genomecov | 1.29 / 1.30 |

Current hotspots (honest — these are open optimization targets, not wins):

| operation | wall× | note |
|---|---|---|
| samtools view BAM→CRAM | **71.5** | CRAM **encoder** — the standout regression (decode is fine, ×1.40) |
| bedtools intersect (self / pair) | **18.1 / 16.8** | the flagship interval op |
| bcftools isec | 5.1 | + RSS ×17 |
| samtools mpileup | 2.9 | RSS **×45** — a memory hotspot |
| samtools stats / bcftools query | 4.2 / 3.9 | per-record overhead on light scans |

The suite turns the prior hand-wavy "faster than the originals" into a precise,
per-operation characterization: genuine wins on heavy re-encode / sort / FASTQ
paths, and concrete, named bottlenecks (CRAM encode, interval intersect,
mpileup memory) under active work before the large-tier headline run.
