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
> 60 k variants) is the first tier where steady-state behaviour shows.

### Medium tier — all 22 cells (measured 2026-06-19, `-reps 5`)

`wall×` = ours/upstream wall-clock (**< 1.0 = we are faster**); `RSS×` flagged
where memory is the open cost. Full per-axis data: `pipeline/.fixtures/medium/bench/bench.{md,json}`.

Wins / parity (wall× ≤ 1.05):

| operation | wall× | note |
|---|---|---|
| sickle se (FASTQ trim) | **0.42** | |
| samtools view BAM→CRAM | **0.72** | CRAM **encoder** — see "transformed" below |
| samtools view BAM→BAM | **0.72** | |
| samtools sort | **0.74** | RSS ×3.0 |
| bedtools intersect (pair / self) | **0.74 / 0.82** | now faster than upstream (was ×1.42 / 1.57); output-path allocs −78% |
| bedtools coverage | **0.80** | now faster than upstream (was ×1.23); allocs −53%, bytes −89% |
| bcftools view | **0.98** | |
| bed sort | **1.05** | |

Modest overhead (1.05 < wall× < 2) — at or near the pure-Go inflate/parse floor:

| operation | wall× | note |
|---|---|---|
| bcftools stats | 1.14 | was ×2.32 — Variant reuse + scratch buffers cut allocs −93% |
| samtools depth | 1.16 | RSS ×100 (per-position arrays) |
| bedtools genomecov | ~1.3 | histogram default is bound by the 16M-int depth array; per-record allocs now ~constant, the **per-base `-d` mode is 7.3× faster** (48M→28 allocs) |
| seqtk seq | 1.31 | tiny op (~170 ms); per-record FASTQ alloc removed (was ×3.12) |
| samtools view CRAM→BAM | 1.34 | CRAM decode |
| samtools stats | 1.41 | RSS ×26; per-record `ReadInto` reuse (−28% allocs) on the single-threaded path; threaded default unchanged |
| bcftools query | 1.61 | parse-bound (htsgo `Scanner.Text` floor); `%POS` per-record alloc removed (allocs/record halved) |
| samtools flagstat | 1.66 | |
| bcftools norm | 1.72 | RSS ×17 |

Remaining hotspots (wall× ≥ 2 — open optimization targets, not wins):

| operation | wall× | note |
|---|---|---|
| samtools mpileup | 3.23 | RSS now **×8.7** (was ×204); CPU ×4.8 (was ×5.2) — see "transformed" |
| bcftools call | 2.17 | **alloc/parse-bound, not libm** (was ×3.04; profiled ~35% alloc/GC, ~25% maps, ~14% float parse, **~2% libm**); VCF reuse + in-place INFO/FORMAT + single-pass list parse cut allocs 32M→18M (−44%), bytes −69%, CPU ×3.4→×2.3 |
| bed merge | 3.48 | tiny op (43 ms) — startup-dominated; allocs cut ~103k→23k |
| bcftools isec | ~4.0 | RSS **×108 → ×28** — per-contig streaming (forward cursors, gated on a sorted/contiguous pre-scan; falls back to load-all otherwise); peak 538→120 MB on the medium fixtures; wall ~neutral |

**Transformed this optimization cycle** (medium tier, before → after):

| operation | was | now |
|---|---|---|
| samtools view BAM→CRAM (encode) | ×71.5 | **×0.72** |
| bedtools intersect (self / pair) | ×18.1 / 16.8 | **×0.82 / 0.74** (now faster than upstream) |
| bedtools coverage | ×1.23 | **×0.80** (now faster than upstream; bytes −89%) |
| bcftools stats | ×2.32 | **×1.14** (allocs −93%) |
| samtools mpileup — **RSS** | ×204 | **×8.7** (CPU ×5.2 → ×4.8) |
| samtools stats | ×4.2 | **×1.41** |
| bcftools query | ×3.9 | **×1.66** |
| bcftools call | ×3.04 | **×2.17** (alloc cut, not faster math) |
| seqtk seq | ×3.12 | **×1.31** |
| bcftools isec — **RSS** | ×108 | **×28** (peak 538→120 MB; per-contig streaming) |
| bed merge | ×4.08 | ×3.48 |
| samtools view BAM→BAM / sort / sickle se | ×0.82 / 0.88 / 0.51 | ×0.72 / 0.74 / 0.42 |

The headline turnaround is CRAM **encode**, which went from a 71× regression to
**faster than upstream** (×0.72) without cgo. The second is `mpileup` **memory**:
it built a per-position event matrix for the whole contig at once (peak RSS
×204); it now streams a single coordinate-sorted input tile by tile, dropping
peak RSS to ×8.7 and trimming CPU (×5.2 → ×4.8), with wall roughly flat (×2.8 →
×3.2 — the tiled walk trades a little scheduling/GC overhead for the memory
collapse, and output stays byte-exact). The remaining wall-time gaps sit in
two honest buckets: (1) tiny ops (`bed merge`, and the now-much-better
`seqtk seq`) whose sub-200 ms runtime is dominated by Go startup/GC, not
throughput; and (2) genuinely heavier paths — `bcftools call` was previously
labelled "libm-bound", but profiling disproved that (only ~2% libm; ~35%
alloc/GC, ~25% maps, ~14% float parse). Its lever was allocation, not faster
math: reusing the input Variant (`ReadInto`), rewriting INFO/FORMAT in place
instead of copying per record, and single-pass comma-list parsing cut
allocations 32M→18M (−44%) and took it from ×3.04 to **×2.17** — all
output-neutral. `isec`/`depth` carry **memory** (RSS), not wall-time, as their
primary remaining cost. CPU-bound
scan paths (stats, query, flagstat, depth) are now at the pure-Go inflate/parse
floor: closing them further would require cgo into libdeflate, which the
project deliberately
forgoes to keep a single static, memory-safe binary (see `CLAUDE.md`).

### Large tier — disk-bound in this environment

A `large`-tier run (192 Mb reference, 2.5 M reads, 400 k variants) was
attempted but is **not reproducible on a small-disk node**: `fixtures.Generate`
materialises the full manifest for the tier unconditionally — including the
`mpileup` and `call` truth VCFs — which comes to **~19 GB** at large scale and
exhausted the container's root filesystem during fixture generation (the run
aborted in `bcftools mpileup` with a write error at 100 % disk). The tier needs
a fat node with ~30 GB+ scratch (consistent with the `bench/README.md` "run on
a fat node / HPC" guidance). Medium is the largest tier that both fits here and
shows steady-state behaviour, so it is the headline tier for these numbers.
