# htsgo roadmap

A proposal to consolidate this repo's scattered format/index/compression
code into a single htslib-equivalent Go library, **htsgo**, and to treat
that library as a first-class deliverable alongside the CLI tools.

**Status: structural migration complete (PR-A–I landed; CRAM and
hfile remain as separate multi-PR follow-ups governed by
`docs/CRAM_DESIGN.md` and a future `HFILE_DESIGN.md`).** This doc
captures the design
decisions so the migration PRs that follow are mechanical, not
re-litigated.

## Why

We are already, accidentally, half-way to a Go htslib. Pieces of it live
in `pkg/bioformats/` (FASTA, FASTQ, BED, VCF, SAM, BCF, GFF, iohelper),
and other pieces — the ones a tool wrote first and never extracted —
live inside individual tool packages:

| Capability                | Currently lives in                                          | Used by (today)                      |
|---------------------------|-------------------------------------------------------------|--------------------------------------|
| BGZF read/write           | `tools/bgzip/pkg/bgzip/`                                    | bgzip, tabix, samtools, bcftools     |
| BGZF (.gzi) index         | `tools/bgzip/pkg/bgzip/index.go`                            | bgzip, samtools faidx, bcftools      |
| BAM read/write            | `tools/samtools/pkg/samtools/` (view.go, sort.go, …)        | samtools                             |
| BAM (.bai) index          | `tools/samtools/pkg/samtools/bai.go`                        | samtools                             |
| Tabix / CSI index         | `tools/tabix/pkg/tabix/` (binning, csi, voffset)            | tabix                                |
| BCF read/write (records)  | `pkg/bioformats/bcf/` (incl. `writer.go`)                   | bcftools `view` (BCF output already works) |
| BCF region/index          | `tools/bcftools/pkg/bcftools/region_bcf.go, index.go`       | bcftools (separate from `pkg/bioformats/bcf`) |
| Virtual offset arithmetic | Duplicated in `tools/tabix/pkg/tabix/voffset.go` and others | tabix, samtools, bcftools            |

Symptoms:

- A `pkg/bioformats/bcf` parser **and writer** already exist
  (`pkg/bioformats/bcf/writer.go` exposes `NewWriter`,
  `NewWriterFromVCFHeader`, `Write(*vcf.Variant)`, etc.), used
  today by `tools/bcftools/pkg/bcftools/view.go` to emit BCF.
  Meanwhile `tools/bcftools` has its own BCF region/index code
  separate from the parser. Either the two share the same on-disk
  record layout (so there's silent duplication) or they don't (so
  behaviour diverges across tools). Both are bad.
- BAM is read inside `tools/samtools/` only — when we port a second
  BAM-aware tool (mosdepth already wants it; future picard/gatk
  re-implementations will), the choice is "copy" or "extract", and
  "copy" keeps winning.
- `pkg/bioformats/README.md` lists SAM/BAM, BGZF, and index support as
  "Future Enhancements" while all of them already exist scattered
  across the tools/ tree. The README is wrong because the surface is
  fragmented.
- The four remaining vcftools flags, the not-yet-started CRAM port,
  faidx, and any future cloud I/O all hit the same missing-foundations
  wall.

## What htsgo is

An in-tree Go library, sited at `pkg/htsgo/`, that owns every
bioinformatics format, index, and compression primitive used by more
than one tool in this repo. Mirrors htslib's scope, not its C API —
the public Go API is idiomatic Go.

What htsgo is **not**:

- Not a CGO wrapper around libhts. The whole point is pure Go.
- Not a separate Go module. Stays inside this repo's single module
  (`github.com/yassineS/bio_ai_experiment`) per CLAUDE.md.
- Not a public/external library at v1. Internal consumers (our tools)
  are the only contract that matters until the surface is stable.

## Treating htsgo as a tool

Per the project convention each tool has a `README.md`, a row in
`tools/PORTING_STATUS.md`, parity notes, and (where applicable) CLI
binaries. htsgo gets the same treatment:

- `pkg/htsgo/README.md` — replaces the current
  `pkg/bioformats/README.md`, lists the full surface, links to
  sub-package docs, names a maintainer agent role.
- A `htsgo` row in `tools/PORTING_STATUS.md` tracking the parity matrix
  (formats supported × {decode, encode, index, region-query}) against
  htslib version pinned in `reference_code/htslib`.
- The CLI binaries htslib itself ships are ported as actual tools that
  thin-wrap htsgo:
  - `tools/bgzip/` (exists) — moves its `pkg/bgzip/` contents into
    `pkg/htsgo/bgzf/`, keeps the CLI.
  - `tools/tabix/` (exists) — same move for `pkg/tabix/`.
  - `tools/htsfile/` (new, small) — htslib's format-sniffer CLI.
  - `tools/faidx/` (new) — wraps `pkg/htsgo/faidx/` once that lands.
    (Note: `samtools faidx` stays as a samtools sub-command for parity,
    but it calls the same package.)

This means htsgo has both a library-tracking surface (the README +
PORTING_STATUS row) and CLI-tracking surfaces (the wrapper tools), and
no part of it is invisible to the porting status spreadsheet.

## Target surface (proposed package layout)

```
pkg/htsgo/
├── README.md                 # consolidated; supersedes pkg/bioformats/README.md
├── iohelper/                 # gzip + stdin/stdout (already exists; moved as-is)
├── bgzf/                     # BGZF reader/writer + .gzi index (from tools/bgzip/)
├── hfile/                    # virtual file: local, gzip, bgzip, http(s), s3, gcs
├── fasta/                    # (moved from pkg/bioformats/fasta)
├── fastq/                    # (moved)
├── bed/                      # (moved)
├── gff/                      # (moved)
├── faidx/                    # FASTA random-access index (.fai). NEW.
├── sam/                      # text SAM (moved from pkg/bioformats/sam)
├── bam/                      # binary SAM + .bai index. NEW; extracted from
│                              #   tools/samtools/ (folds .bai into bam/ to
│                              #   avoid the reader-writes-index cyclic-import
│                              #   trap; htslib does the same)
├── cram/                     # see docs/CRAM_DESIGN.md; .crai sits as cram/index.go
│   └── codec/                # rANS 4x8, rANS 4x16, LZMA — sole 3rd-party-dep zone
├── vcf/                      # text VCF (moved)
├── bcf/                      # binary VCF (moved from pkg/bioformats/bcf); adds WRITE path
├── tabix/                    # tabix/CSI index (moved from tools/tabix/)
├── region/                   # parser for "chr1:100-200" syntax, shared by all
└── internal/                 # binary readers, kstring-like helpers, etc.
```

`pkg/bioformats/` becomes a deprecation shim — re-exports for one
release, then deleted — so tool imports can migrate in a single
mechanical sweep PR.

## Gap analysis

What's missing today, in priority order:

### P0 — unblocks shipped tools

1. **Wire vcftools BCF flags onto the existing `bcf.Writer`.** The
   `pkg/bioformats/bcf` package **already exposes a working
   writer** (`NewWriter`, `NewWriterFromVCFHeader`, `Write(*vcf.Variant)`,
   `WriteRecord`, `Flush`) used by `tools/bcftools/pkg/bcftools/view.go`
   to emit BCF. What's missing is wiring the four vcftools BCF
   flags (`--bcf`, `--diff-bcf`, `--recode-bcf`, `--contigs`) onto
   that writer — see `tools/vcftools/README.md:240-241`. This is a
   tool-level integration job, not a library-level write-path port.
   **~3-5 days**, decoupled from the A-F extraction chain — can land
   independently before, during, or after the moves.
2. **BAM extracted from samtools** into `pkg/htsgo/bam`. No new
   functionality, just lifts the existing code to a shared home so
   mosdepth and future tools stop copying it. ~3 days incl. tests.
3. **BGZF extracted from `tools/bgzip`** into `pkg/htsgo/bgzf`, same
   reasoning. Concurrent with the BAM extraction. ~2 days.
4. **Tabix/CSI extracted** into `pkg/htsgo/tabix`. ~2 days.

### P1 — unblocks the next wave

5. **faidx** (`pkg/htsgo/faidx`). Enables `samtools faidx`, reference
   FASTA random-access for variant callers, CRAM reference resolution.
   ~1 week.
6. **region/ package** — one parser for `chr:start-end` syntax, used by
   every region-aware sub-command across samtools/bcftools/tabix. ~2
   days; mostly extracting + de-duplicating what exists.

### P2 — opens the format frontier

7. **CRAM read** (v3.0 first, then v2.1 + v3.1). See
   `docs/CRAM_DESIGN.md`. ~4-6 weeks. This is the one place
   third-party deps are sanctioned (rANS, LZMA), confined to
   `pkg/htsgo/cram/codec/`.
8. **CRAM write** (v3.0). ~3-4 weeks on top of read.

### P3 — opens cloud workflows

9. **hfile virtual filesystem** with HTTP(S) ranged reads and (later)
   S3/GCS. Only matters once a real user needs cloud I/O; ranking it
   last avoids speculative complexity. ~2-3 weeks for HTTPS+ranged;
   cloud backends each ~1 week.

## Migration plan

Phased so each PR is reviewable in isolation and `main` is always
green. The green-main mechanism for every extraction PR (C/D/E) is
**temporary re-export shims at the tool-package path**: lift the
code into `pkg/htsgo/<x>/`, then leave the old `tools/<tool>/pkg/<x>/`
exposing the same identifiers as type aliases + function re-exports.
Consumers (tests, other tools) keep working unchanged until PR-I
flips imports in one mechanical sweep. The pattern mirrors
`pkg/bioformats/`'s own shim role during PRs A/B.

- ~~**PR-A: skeleton.** Create empty `pkg/htsgo/`, move `iohelper` and
  the four uncontroversial format packages (fasta, fastq, bed, gff)
  from `pkg/bioformats/`. Add re-export shims so nothing breaks. Update
  `pkg/htsgo/README.md`. Also update the path reference in
  `docs/CRAM_DESIGN.md` (`pkg/bioformats/cram/codec/` →
  `pkg/htsgo/cram/codec/`) and the stale "Future Enhancements" line
  in `pkg/bioformats/README.md`. No new functionality.~~
  **Landed.** All five leaf packages moved with their tests; shims
  at the old `pkg/bioformats/{iohelper,fasta,fastq,bed,gff}/` paths
  re-export through type aliases + function variables so the ~220
  in-tree importers keep working unchanged. `docs/CRAM_DESIGN.md`
  paths updated; `pkg/bioformats/README.md` rewritten as a
  deprecation pointer; `pkg/htsgo/README.md` is the new canonical
  inventory.
- ~~**PR-B: SAM/VCF/BCF move.**~~ **Landed.** The three remaining
  `pkg/bioformats/` packages relocated to `pkg/htsgo/`. After this
  PR, `pkg/bioformats/` contains only re-export shims (the
  `bcf.Writer` API is unchanged; the in-tree `bcf` import inside
  `pkg/htsgo/bcf` was rewired to import `pkg/htsgo/vcf` directly
  rather than going back through the bioformats shim).
- ~~**PR-C: extract BGZF.**~~ **Landed.** Five source files moved
  from `tools/bgzip/pkg/bgzip/` to `pkg/htsgo/bgzf/`; the package
  name flipped from `bgzip` to `bgzf` to match htslib's terminology
  and the htsgo target tree. The old `tools/bgzip/pkg/bgzip/`
  directory now holds only a re-export shim (kept under the legacy
  `bgzip` package name for backwards compatibility). The
  `tools/bgzip/cmd/bgzip/` CLI binary still builds via the shim.
- ~~**PR-D: extract BAM + BAI.**~~ **Partially landed.** The BAI
  format primitives (`bai.go` / `bai_test.go`) and the BAM →
  BAI builder (`BuildBAI`) moved from `tools/samtools/pkg/samtools/`
  to `pkg/htsgo/bam/`. The BAM reader/writer themselves already
  migrated as part of PR-A/B (they live in `pkg/htsgo/sam/`
  alongside SAM since htslib keeps them together too); a future
  PR can re-split BAM into `pkg/htsgo/bam/` if the SAM/BAM
  decoupling becomes more important than the shared-record-types
  ergonomics. Shim at `tools/samtools/pkg/samtools/bai_shim.go`
  re-exports the moved BAI surface; the in-tree `samtools index`
  orchestration (the `Index` / `IndexFile` functions) keeps
  living in `tools/samtools/` because it's CLI-level, not
  format-level.
- ~~**PR-E: extract Tabix/CSI.**~~ **Landed.** Ten source files moved
  from `tools/tabix/pkg/tabix/` to `pkg/htsgo/tabix/` (Tabix/TBI
  + CSI format primitives + the binning helpers `Reg2bin` /
  `LinearTile` shared by BAI). Shim at the old path keeps the
  tabix CLI building. As a bonus this PR fixed the
  htsgo-imports-tools antipattern: `pkg/htsgo/sam/bam_writer.go`
  and `pkg/htsgo/bam/bai.go` previously imported the binning
  helpers from `tools/tabix/pkg/tabix` (library reaching into a
  tool package); both now import `pkg/htsgo/tabix` directly. The
  package's `pkg/htsgo/bgzf` import was also tightened (no longer
  routed through the bgzip shim).
- ~~**PR-F: extract region parser**~~ **Partially landed.** The
  canonical region parser (the most complete of the four
  scattered variants — samtools's) moved from
  `tools/samtools/pkg/samtools/region.go` to
  `pkg/htsgo/region/`. Surface: `Region`, `ParseRegion`,
  `ResolvedRegion`, `ResolveRegions`, `Region.OverlapsRef`. The
  BAI-specific `UnionChunks` aggregator moved alongside its
  `BAIIndex` / `BAIChunk` dependencies into `pkg/htsgo/bam/`
  (kept out of the region package so it stays format-agnostic).
  A `region_shim.go` in samtools re-exports both surfaces. The
  other variants (`tools/bcftools/pkg/bcftools/{view.go,gtcheck.go}`,
  `tools/tabix/cmd/tabix/main.go`) are file-local one-offs with
  slightly different semantics (0/1-based, comma-stripping,
  half-open vs inclusive); migrating them needs per-call-site
  semantic equivalence checking and is deferred to a follow-up.
- **PR-G: wire vcftools BCF flags.** Plug `--bcf` / `--diff-bcf` /
  `--recode-bcf` / `--contigs` onto the existing `bcf.Writer`. The
  scope is vcftools integration, not a `pkg/htsgo/bcf` write-path
  port (that path is already in `pkg/bioformats/bcf/writer.go`). **This
  PR is independent of A-F** and can land any time — it doesn't
  depend on the extractions and the extractions don't depend on it.
- **PR-H: faidx + region-iterator API + `tools/htsfile`** — partial
  landing. The **`tools/htsfile`** CLI shipped (a pure-Go htsfile
  re-implementation that sniffs SAM/BAM/CRAM/VCF/BCF/FASTA/FASTQ/
  BED/GFF wrapped in plain/gzip/BGZF, with the same one-line
  summary form upstream emits). Implementation in
  `tools/htsfile/pkg/htsfile/` (kept inside the tool tree because
  the heuristics are CLI-shaped). The other two halves —
  `pkg/htsgo/fasta` polish and a unified region-iterator API —
  are still TBD; they ride on no current consumer and can
  follow on demand.
- ~~**PR-I: drop `pkg/bioformats/` and tool-package shims**~~
  **Landed.** All ~220 importers swept from
  `pkg/bioformats/<x>` to `pkg/htsgo/<x>` in a single mechanical
  commit; the 10 shim files at the old paths plus
  `pkg/bioformats/README.md` removed; the now-empty
  `pkg/bioformats/`, `tools/bgzip/pkg/bgzip/`, and
  `tools/tabix/pkg/tabix/` directories all gone. bgzip imports
  routed through `pkg/htsgo/bgzf` with a `bgzip` package alias
  so call sites kept their existing `bgzip.X` qualifier.
  The two in-samtools sub-shims
  (`tools/samtools/pkg/samtools/bai_shim.go` and
  `region_shim.go`) were initially deferred, then removed in
  the PR-I follow-up: the samtools subcommands now qualify
  `bam.BAIIndex` / `region.ParseRegion` directly. The feared
  `Region`-type-vs-field collision turned out not to bite —
  the bare `Region` type is never referenced unqualified in
  the subcommand code (only the `.Region` field access on
  `ResolvedRegion`, which is untouched).

CRAM (PR-J onwards) and hfile (PR-K onwards) are their own multi-PR
projects governed by `docs/CRAM_DESIGN.md` and a future `HFILE_DESIGN.md`.

### Testing strategy during the migration

Tests **move with the code**, never get duplicated. For every
extraction PR (C/D/E/F):

1. Move `*_test.go` files alongside the code into `pkg/htsgo/<x>/`.
2. The shim package at the old path gets **no tests** (it's pure
   re-export; if it broke its tests would too).
3. Parity tests (`*_parity_test.go`) that depend on
   `reference_code/<tool>` submodules must be moved together with
   the submodule init line in their setup; the submodule itself
   stays at `reference_code/<tool>/` (CLAUDE.md convention).
4. Every PR's checklist requires `go test -race ./...`,
   `go vet ./...`, `gofmt -l .` (all empty/clean) before request
   for review. `markdownlint **/*.md` for the docs PRs.
5. Add `reference_code/biogo-hts` and `reference_code/htscodecs`
   as git submodules in PR-A so future contributors can do delta
   reviews against pinned upstream SHAs when porting chunks.

## Inspiration / source recycling

- **biogo/hts** (BSD-3, github.com/biogo/hts) covers BGZF, BAM, SAM,
  and tabix in pure Go. We read it as a reference and may port chunks
  verbatim **with proper attribution in package-level comments and a
  `NOTICE` file at repo root**. We're not vendoring it — we want
  control over the API and the option to evolve it. CRAM (their
  issue #54) is unimplemented upstream; we are not blocked on them.
- **htslib** (MIT, github.com/samtools/htslib) is the behavioural spec.
  Parity tests in tool packages (`*_parity_test.go`) already compare
  against the upstream binary; the same pattern extends to htsgo's
  format-level tests where useful.
- **htscodecs** (BSD-3 / MIT dual-licensed, github.com/samtools/htscodecs)
  is the C reference for the CRAM codec layer. Port-from, not
  link-against. Either licence is compatible with this repo.

License hygiene: every file derived from these gets a header
attributing the source and noting "ported to Go, adapted from <repo>
commit <sha>". The `NOTICE` file enumerates upstream attributions in
one place.

## Dependency policy

Unchanged from CLAUDE.md and `docs/CRAM_DESIGN.md`:

1. Go stdlib first.
2. In-tree pure-Go implementation second.
3. Third-party Go dep only inside `pkg/htsgo/cram/codec/` for rANS and
   LZMA. Any other proposed dep needs its own ADR.

## Non-goals

- **htslib C-API compatibility.** Go callers get a Go API.
- **Plug-in architecture.** htslib's plug-in mechanism (hfile, codec)
  is fine for a C library shipped as `.so`. We compile static binaries;
  registries with runtime registration are the wrong shape.
- **Multi-threading inside readers.** Concurrency is a tool-level
  decision (sort, mpileup). Readers stay single-goroutine and
  reentrant-by-cloning.
- **Zero-copy claims we can't back up.** Where we copy, we copy. Where
  we don't, we benchmark.

## Open questions

- Should `pkg/bioformats/` survive as a permanent re-export shim for
  external users (we have no announced external users today), or be
  deleted in PR-I? Default: delete.
- Where does `pkg/cliflag/` live? It's not htsgo-scoped. Stays at
  `pkg/cliflag/`.
- Naming: `pkg/htsgo/` vs `pkg/hts/`. Default: `pkg/htsgo/` —
  unambiguous, googleable, no clash with any well-known Go package.

## Success criteria

The roadmap is "done" when:

- `pkg/bioformats/` no longer exists (or is a one-line shim slated for
  removal).
- No tool package contains BGZF, BAM, BAI, Tabix, CSI, or region-parser
  code that another tool also needs.
- vcftools shows 100% flag parity (the four BCF binary I/O flags
  wired onto the existing `bcf.Writer` — see PR-G).
- A new `tools/htsfile/` CLI exists and round-trips through every
  format htsgo supports.
- `pkg/htsgo/README.md` is a true inventory, not a roadmap of vapour.
- CRAM v3.0 read works on the htslib test corpus; CRAM is the only
  sub-package with a third-party dep.
