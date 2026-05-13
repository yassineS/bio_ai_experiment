# Session handoff — 2026-05-13

Brief notes for whoever picks this up next (most likely me in another Claude
Code session). The session that produced these notes was approaching its usage
limit, so the to-do list and the rationale matter more than fresh tooling work.

## What landed in this batch of sessions

Every PR below has been **merged to `main`**. See `git log` for the squash
commits if you need the exact diffs.

### Audit + corrections (May 12)

- PR #23 — added `CLAUDE.md` (project conventions for AI assistants).
- PR #24 — vcftools: fixed `--site-pi` formula; implemented `--window-pi`,
  `--TsTv-by-count`, `--depth`; **rejected** the unimplemented options
  with a clear error instead of silently accepting them; corrected the
  overstated `tools/PORTING_STATUS.md` (real coverage numbers, dropped
  "feature parity" claims).
- PR #25 — vcftools `--TsTv-by-qual`.
- PR #26 — seqtk `subseq` (upstream-compatible: name list or BED region).
- PR #27 — skewer: fixed **two real bugs** while pushing coverage 46→90%:
  `TrimSingleEnd --auto-detect` was a no-op (consumed input and "rewound"
  by re-wrapping the now-exhausted reader); `detectAdapters` was
  non-deterministic (Go map iteration).
- PR #28 — first `PORTING_STATUS.md` refresh.

### Implementation push (May 13)

- PR #29 — fastp sliding-window quality trimming (`--cut_front/--cut_tail/--cut_right` + `--cut_window_size/--cut_mean_quality`).
- PR #30 — bedmerge `bedtools merge`-style `-c`/`-o` column aggregation (sum/min/max/mean/median/count/count_distinct/distinct/collapse/first/last/mode/antimode).
- PR #31 — prinseq tests 55→90%.
- PR #32 — seqtk `mergepe` + `cutN`.
- PR #33 — sickle Phred-encoding auto-detect (default `-t auto`, peek-based — no consume-and-rewind bug, the lesson from PR #27).
- PR #34 — vcftools Weir & Cockerham 1984 Fst (`--weir-fst-pop`, `--fst-window-size/--fst-window-step`); the worked-example test matches the textbook to 6 decimals.
- PR #35 — PORTING_STATUS refresh.
- PR #36 — formalised **POSIX-compliant CLIs** as a parity requirement (alongside 100% feature parity) in README + docs/CLI_CONVENTIONS.md.
- PR #37 — skewer coverage 46→99.7% (then 100% via PR #40).
- PR #38 — bedmerge coverage 79→**100%**.
- PR #39 — prinseq coverage 90→99.9% + fixed a real bug: `graph.go`'s three `generate*SVG` helpers were dropping every `fmt.Fprintf` error, masking truncated/broken-pipe writes.
- PR #40 — skewer `ToJSON` drop unreachable error return (99.7→**100%**).
- PR #41 — bedtools: `bedsort`, `bedslop`, `bedcomplement`.
- PR #42 — bedtools: `bedgenomecov`, `bedjaccard`.
- PR #43 — bedtools: `bedsubtract`, `bedflank`, `bedclosest`.

## Current state of `tools/`

| Tool | Coverage | Status |
|------|---------:|--------|
| **bedmerge** | 100.0% | bedtools-merge parity-ish, with `-c`/`-o` |
| **skewer** | 100.0% | adapter trimming, auto-detect works |
| **prinseq** | 99.9% | only `report.go:15-17` defensive line uncovered |
| **bedjaccard** | 96.3% | new |
| **bedslop** | 95.2% | new |
| **bedcomplement** | 94.6% | new |
| **bedgenomecov** | 94.0% | new |
| **bedsubtract** | 93.7% | new |
| **bedclosest** | 92.8% | new |
| **bedflank** | 92.2% | new |
| **bedsort** | 91.6% | new |
| **sickle** | 81.9% | Phred auto-detect on by default |
| **bedintersect** | 75.1% | (untouched this round; could use a coverage push) |
| **fastp** | 66.6% | + sliding-window cut |
| **seqtk** | 66.2% | 8 subcommands now |
| **vcftools** | 58.3% | ~50/147 options; the largest remaining gap is LD |

For the cross-tool view see [`tools/PORTING_STATUS.md`](../tools/PORTING_STATUS.md).

## Outstanding work — explicit user asks

These are the things the user said to do "next" (i.e. in the next session):

1. **Finish what we started.** "Finish all the tools we started working on."
   That means closing the gap to **100% feature parity AND POSIX-compliant CLI**
   (per `docs/CLI_CONVENTIONS.md`) for each of the 8 bed* tools, plus
   `bedintersect`, `vcftools`, `fastp`, `seqtk`, `sickle`. The largest pieces are:
   - **vcftools**: LD analysis (`--geno-r2`, `--hap-r2`, all the `--ld-window-*`
     options). The other smaller misses are listed in
     `tools/vcftools/ROADMAP.md` / `FEATURE_COMPARISON.md`.
   - **fastp**: HTML / JSON reports, automatic adapter detection.
   - **seqtk**: `mutfa`, `randbase`, `hpc`.
   - **sickle**: validated byte-for-byte parity check against the C original.
   - All bed* tools: validated parity against upstream `bedtools <subcommand>`
     using its test suite (haven't done this yet — coverage numbers don't
     guarantee output identity).

2. **bgzip/gzip support, when needed.** The user proposed:
   - `compress/gzip` (stdlib) is already used via
     `pkg/bioformats/iohelper`. Reuse that.
   - `github.com/biogo/hts/bgzf` for proper **BGZF** support (block-gzip used
     by tabix, BAM, BCF, gzipped-and-indexed VCF). **There is no stdlib BGZF
     equivalent.** Adding `biogo/hts/bgzf` would be this repo's **first
     third-party Go dependency** — CLAUDE.md currently says "avoid adding
     external dependencies without a strong reason". BGZF support for
     vcftools (`.vcf.gz` is BGZF, not vanilla gzip) is the strong reason.
     Plan: add it to `pkg/bioformats/iohelper` so it's the one place that
     touches the third-party dep, gate it on the magic-byte sniff
     (`1f 8b 08 04 ...` plus the BGZF EOF block), and keep regular gzip on
     stdlib. Then wire it through vcftools (and anywhere else that needs
     to read a `.gz` that might actually be BGZF, e.g. tabix-indexed VCF).
   - This is a Go 1.21-compatible dep — verify before adding (`go get`).

3. **Re-evaluate the tool ranking.** The user's exact framing:
   > "What we want from this project is deliver value to people, not
   > re-implement tools no one uses anymore."

   The current ranking lives in
   `analysis/top_50_packages_for_improvement.md` and `analysis/top_50_packages.json`.
   Some entries are clearly dated (PRINSEQ-lite, sickle, skewer — all
   ~2012-era tools that have been supplanted by **fastp** for QC and by
   modern wrappers like **Trim Galore** for adapter trimming). Concrete
   suggestion for the re-rank:
   - Pull recent **bioconda** download stats (`bioconda-recipes` repo or
     anaconda.org JSON API) for the last 12 months, not the lifetime total.
   - Cross-check with **PyPI** / **CRAN** / **Bioconductor** download stats
     where the upstream lives there.
   - Look at recent (post-2022) Nature Methods / NAR Bioinformatics tool
     papers to spot what's actually being adopted (e.g. `minimap2`,
     `mafft`, `mmseqs2`, `STAR`, `salmon`, `mosdepth`, `bcftools`, `samtools`,
     `tabix`, `bgzip`, modern callers like `deepvariant`, `clair3`).
   - Be ruthless about **demoting**: PRINSEQ-lite is functionally dead;
     sickle is unmaintained since 2014; skewer hasn't been touched since 2016.
     Keep what *we already did* (sunk cost is fine, the work is done) but
     don't pour more effort into 100% parity for tools no one runs.

## Smaller follow-ups noted but deferred

- **prinseq quirk** (`tools/prinseq/pkg/prinseq/prinseq.go`): the
  `trimQualityLeft` / `trimQualityRight` helpers assume Phred+33 regardless
  of the `QualType` option. Flagged but not fixed because the "correct"
  behaviour couldn't be justified against the README in the session that
  found it. Should be checked against upstream PRINSEQ-lite source and
  either fixed or documented.
- **bedintersect coverage push** (75% → ≥90%) — same pattern as bedmerge/skewer.
- **`cmd/` entry points have 0% coverage** across the whole repo. We could
  factor each tool's main into a thin `run(args, stdin, stdout, stderr) int`
  function and test that. Not urgent — the CLI is exercised end-to-end via
  per-tool integration tests already.
- **Validated parity** for every tool against its upstream's own test suite,
  not just "the tests we wrote pass". This is the only way to actually claim
  "feature parity"; coverage numbers alone don't do it.

## Conventions reminder (cheap to forget)

- Single Go module at the root, **no per-tool `go.mod`**.
- No third-party deps yet. Adding `biogo/hts/bgzf` (see above) is the first
  one and should be a deliberate, documented decision in CLAUDE.md when it
  happens.
- POSIX-compliant CLIs are now a **hard requirement** for "done" status (see
  `docs/CLI_CONVENTIONS.md`). New tools use `pkg/cliflag` for short + long
  flags from day one.
- One PR per tool (or per coherent feature). Keep `tools/PORTING_STATUS.md`
  in a separate doc-only PR after a wave of merges, not bundled in.
- Agents run inside the main checkout (despite the worktree-isolation API);
  brief them to scope `git add tools/<their tool>` rather than `git add -A`,
  and pick non-overlapping directories so the parallel waves don't stomp
  each other.

## Where I left off

`main` at the time of writing this file:

```
5b0a184 bedtools: add bedsubtract, bedflank, bedclosest (#43)
58d5909 bedtools: add bedgenomecov, bedjaccard (#42)
57dc852 bedtools: add bedsort, bedslop, bedcomplement (#41)
de66eea skewer: drop unreachable error return in ToJSON (99.7 -> 100%) (#40)
2d8a00d prinseq: 90% -> ~100% coverage + fix silent SVG-write-error bug (#39)
```

No PRs open. Working tree clean.
