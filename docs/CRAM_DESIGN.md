# CRAM design notes (work in progress)

CRAM is the htslib-stack format we have **not** yet ported. This doc
captures the up-front decisions so when we start the implementation we're
not re-litigating them. **Status: in execution — see
[`CRAM_ROADMAP.md`](CRAM_ROADMAP.md) for the live plan.**

> **Note (post-C1):** the "Dependency policy" and "Open questions"
> sections below pre-date the first implementation PR. The rANS
> port-vs-dep question was resolved **in favour of an in-tree pure-Go
> port** — `pkg/htsgo/cram/codec` ships rANS 4x8 with no third-party
> dep, proven byte-exact against the htscodecs corpus. The only
> sanctioned CRAM dep is now `ulikunitz/xz` for LZMA. The
> authoritative, up-to-date decisions live in
> [`CRAM_ROADMAP.md`](CRAM_ROADMAP.md) §1.2 and §5.

## Why CRAM is held up

- **Reference-dependent.** Reads are stored as differences from a
  FASTA; you cannot decode without the exact reference (verified via
  embedded MD5 + optional `m5:` URL / `REF_PATH` / `REF_CACHE`
  resolution).
- **Multiple custom compression codecs**, none in Go stdlib:
  - rANS 4x8 (CRAM v3.0)
  - rANS 4x16 (CRAM v3.1)
  - per-block: gzip ✅, bzip2 ✅, LZMA ❌, RAW ✅
- **Multiple version dialects.** v2.1 / v3.0 / v3.1 all in use, v4
  emerging. Real-world data has all of them.
- **Tree container format** (file → container → slice → block) with
  ~150 pages of spec.
- **Entropy-coder zoo**: external / byte_array_stop / byte_array_len /
  huffman / beta / subexp / gamma / golomb-rice. Each data series
  picks one.
- **Lossy modes**: quality-binning schemes change the embedded MD5,
  meaning even *reading* lossy CRAM is conditional on knowing the
  binning.
- **`.crai` index format** — different from `.bai`/`.csi`.

Conservative effort: **read-only port ≈ 4-6 weeks** focused work for
one engineer, **read+write ≈ 8-12 weeks**. So this is its own
multi-PR project, not a single slice.

## Dependency policy

**Owner has OK'd third-party dependencies for the CRAM codec layer.**
The "no third-party deps" rule (CLAUDE.md) is therefore relaxed for:

- rANS 4x8 / 4x16 implementations (no Go-stdlib equivalent; ~1,500
  lines of careful bit-twiddling to port from htslib's C SIMD-ish
  reference).
- LZMA decompression (no Go-stdlib equivalent; the `xz` package is
  the de-facto Go choice).

The dep MUST be confined to a single sub-package (proposed:
`pkg/htsgo/cram/codec/`) so the rest of the repo can still
honestly claim "stdlib only" for non-CRAM workflows.

Preference order remains:

1. Go stdlib.
2. In-tree Go implementation (the bgzip / tabix / sam / bcf packages
   all went this route).
3. Third-party Go dep (allowed for CRAM codec primitives only).

Candidate libraries to evaluate when the time comes:

- **rANS**: there is no widely-used pure-Go rANS library. Options:
  - Port htslib's `htscodecs` C reference to Go (the cleanest
    long-term answer).
  - Port `JKKBenchmarks/rans` (small, MIT, used by jankratochvil's
    Go fastcompress experiments).
  - Cgo binding to `htscodecs` (last resort — drags libc into builds).
- **LZMA**: `github.com/ulikunitz/xz` (BSD-2, mature, MIT-compatible).
- **bzip2**: stdlib `compress/bzip2` is read-only; for write use
  `github.com/dsnet/compress` if needed (BSD-3).

## Version-support matrix (proposed v1)

| Version | Decode | Encode | Notes |
|---------|:------:|:------:|-------|
| 2.1     | ✅     | ❌     | Legacy; archives still have it. Decode-only is enough. |
| 3.0     | ✅     | ✅     | Most-deployed format. Full read+write. |
| 3.1     | ✅     | ✅     | Adds rANS 4x16; full read+write. |
| 4.0     | ❌     | ❌     | Spec not finalised. Defer. |

Lossy modes: **decode** all of them (we have to, to read existing
files). **Encode**: lossless only in v1; lossy quality binning is its
own follow-up.

## Reference resolution

Match samtools's `REF_PATH` / `REF_CACHE` semantics:

- Honour `REF_PATH` (colon-separated list of URL templates with
  `%s` for the MD5).
- Honour `REF_CACHE` (filesystem cache root for downloaded refs).
- Prefer a passed-in `--reference FASTA`; fall back to env-resolved
  paths; fall back to the canonical EBI URL
  (`https://www.ebi.ac.uk/ena/cram/md5/%s`) only when the user
  explicitly opts in (`--use-cram-online-ref`).
- MD5 verification on every read; surface mismatches as a clear error,
  not a silent decode failure.

## Milestone breakdown (proposed)

1. **`pkg/htsgo/cram/codec/`** — rANS 4x8 + 4x16 + LZMA wrappers,
   with byte-for-byte test fixtures from htslib's `htscodecs`. Allowed
   to use third-party deps. ~2 weeks.
2. **`pkg/htsgo/cram` reader (v3.0)** — file/container/slice/block
   parser, data-series decoders, MD5+REF_PATH plumbing, `.crai` index
   read. Pure-Go on top of (1). ~3 weeks.
3. **CLI plumbing** — wire the new reader through
   `pkg/htsgo/iohelper` (auto-detect CRAM by magic bytes) so
   `samtools view`, `samtools depth`, `samtools fastq`,
   `samtools mpileup`, etc. all transparently accept CRAM input. ~1
   week.
4. **`pkg/htsgo/cram` reader (v3.1)** — additional rANS 4x16
   wiring + edge cases. ~1 week.
5. **`pkg/htsgo/cram` writer (v3.0)** — encode side. ~3 weeks.
6. **`pkg/htsgo/cram` writer (v3.1)**. ~1 week.
7. **`samtools view --output-fmt cram`** + `samtools index` for
   `.crai`. ~1 week.
8. **Lossy quality binning** (encode side). ~2 weeks. Optional.

**Total**: ~14 weeks engineering effort assuming no major unknowns.

## When to start

After:

- bcftools subcommand tail closure (`annotate`/`merge`/`isec`/`csq`)
- samtools subcommand tail closure (`markdup`/`merge`/`idxstats`/`stats`)

Both are smaller, parallelisable across agents, and individually
shippable. Better return per PR-day until they're done. Then we set
aside a dedicated CRAM block.

## Open questions for the start of the project

1. Confirm rANS port-vs-dep choice. (Owner has OK'd dep.)
2. Pick a specific rANS Go library or commit to porting `htscodecs`.
3. Decide on LZMA library (`ulikunitz/xz` is the obvious pick).
4. Decide whether v2.1 decode is required for v1 ship or can be
   deferred.
5. Decide on the EBI online-ref fallback default (off in v1?).
