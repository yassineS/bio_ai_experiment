# CRAM execution roadmap

Companion to [`CRAM_DESIGN.md`](CRAM_DESIGN.md) (the *why* and the
up-front decisions). This doc is the *how*: the shopping plan, the
PR-by-PR execution sequence, and the testing + compliance strategy.

**Status: in execution.** Updated as milestones land.

---

## 1. Shopping plan — what we need to acquire

### 1.1 Reference material (read-only, for porting + compliance)

| Item | Purpose | How |
|------|---------|-----|
| `samtools/htscodecs` | The C reference for rANS 4x8 / 4x16 + the canonical test vectors. | git submodule under `reference_code/htscodecs` |
| `samtools/htslib` | CRAM container reference (`cram/` subtree) + the spec-conformance CRAM test files. | git submodule under `reference_code/htslib` |
| CRAM spec PDF (v3.0, v3.1) | The normative format definition. | Linked, not vendored (hsformats.github.io). |

If the sandbox network policy blocks submodule init, the port
proceeds from spec knowledge and the codec PRs rely on
round-trip + property tests; a `t.Skip`-gated compliance-vector
suite is wired so it activates the moment the submodule is
present. This is a documented, tracked gap — never a silent one.

### 1.2 Third-party Go dependencies

Per `CRAM_DESIGN.md` the no-deps rule is relaxed *only* inside
`pkg/htsgo/cram/codec/`. Revised picks after a fresh assessment:

| Codec | Decision | Rationale |
|-------|----------|-----------|
| rANS 4x8 | **In-tree pure Go** (preference #2, not #3) | ~300-400 LOC decoder + similar encoder. A textbook static-rANS coder — well within the codebase's in-tree appetite (bgzf/tabix/bcf all went this route). No dep. |
| rANS 4x16 | **In-tree pure Go** | Same family; the 4x16 normalisation + the v3.1 "order bits" (RLE / bit-packing / striping transforms) add ~400 LOC. Still tractable in-tree. |
| LZMA (decode) | **`github.com/ulikunitz/xz`** (BSD-2) | Genuinely hard to port; LZMA is a rare optional per-block CRAM codec. This is the *one* sanctioned dep. Confined to `codec/lzma.go`. |
| gzip / bzip2 / RAW | **stdlib** (`compress/gzip`, `compress/bzip2`) | bzip2 stdlib is decode-only — fine, CRAM rarely bzip2-encodes and v1 write uses gzip/rANS. |

**Net effect:** the only new `go.mod` entry for the whole CRAM
project is `ulikunitz/xz`, and it is reachable only from
`pkg/htsgo/cram/codec/lzma.go`. Everything else stays
stdlib + gonum.

LZMA support is **deferred past the first reader PR** — a CRAM
block compressed with LZMA returns a clear `unsupported codec`
error until `codec/lzma.go` lands. Real-world CRAM almost never
uses it.

---

## 2. Execution roadmap — the PR sequence

Each row is one reviewable PR. "Gate" = what must be green to merge.

| PR | Scope | Gate |
|----|-------|------|
| **C1** | `pkg/htsgo/cram/codec/` — rANS 4x8 (order-0 + order-1), decode **and** encode. | Round-trip + property tests; compliance vectors if htscodecs submodule present. |
| **C2** | rANS 4x16 **order-0** decode + encode, plus the framing format byte and the X_CAT store-uncompressed fallback. Order-1 and the transform bits (PACK / RLE / STRIPE / X32) are rejected with a clear error. | Same as C1 — byte-exact compliance vs the `r4x16/*.0` vectors landed. |
| **C2.1** | rANS 4x16 **order-1** context model. | Byte-exact vs the `r4x16/*.1` vectors. |
| **C2.2** | The v3.1 transform bits — X_PACK (bit-packing), X_RLE, X_STRIPE (4-way), X_32 unrolling. | Byte-exact vs the remaining `r4x16/*.{4,5,64,…}` vectors. |
| **C3** | `pkg/htsgo/cram` container parser — file def, container, compression header, slice header, block. No data-series decode yet; just the tree walk + per-block decompress dispatch. | Parse + walk every container in a real v3.0 CRAM without error. |
| **C4** | CRAM v3.0 **read** — data-series decoders (the entropy-coder zoo: external / byte_array_stop / byte_array_len / huffman / beta / subexp / gamma / golomb-rice), record reconstruction, the reference-diff decode. | Decode a v3.0 CRAM to SAM records matching `samtools view` output. |
| **C5** | Reference resolution — `--reference`, `REF_PATH`, `REF_CACHE`, MD5 verify; `.crai` index read. | MD5-mismatch surfaced as a clear error; `.crai` region query works. |
| **C6** | CLI plumbing — `iohelper` CRAM magic-byte autodetect so `samtools view/depth/fastq/mpileup` accept CRAM input transparently. | The samtools subcommands round-trip a CRAM fixture. |
| **C7** | CRAM v3.1 read — wire rANS 4x16 + 3.1 edge cases. | v3.1 fixtures decode. |
| **C8** | CRAM v3.0 **write** — encoder, container/slice assembly, codec selection. | Written CRAM re-reads byte-identically through our own reader; `samtools view` reads it. |
| **C9** | `samtools view --output-fmt cram` + `samtools index` `.crai` write. | End-to-end SAM→CRAM→SAM via the CLI. |
| **C10** | CRAM v3.1 write. | v3.1 round-trip. |
| **C11** *(optional)* | Lossy quality binning (encode). | Documented, opt-in. |
| **C-LZMA** *(slots in before C4 if a fixture needs it)* | `codec/lzma.go` via `ulikunitz/xz`. | LZMA-compressed block decodes. |

v2.1 decode is **deferred** — see Open Questions resolution below.
It slots in as a C4.1 only if real archive data forces it.

---

## 3. Testing strategy

Three layers, every PR:

1. **Round-trip / property tests** (always). Encode→decode is the
   identity; decode→encode→decode is stable; random inputs across
   the size/symbol-distribution space don't panic and reproduce.
   These run with no external fixtures.
2. **Compliance vectors** (when the htscodecs / htslib submodules
   are present). Decode the upstream test corpus and assert
   byte-for-byte. Gated behind a fixture-presence check that
   `t.Skip`s cleanly when the submodule isn't initialised — so CI
   without submodules still passes and CI with them gets the
   stronger guarantee.
3. **Interop tests** (from C4 onward). Our decoder vs. real
   `samtools view` output on the same CRAM; our encoder's output
   fed back through real `samtools`. Gated behind the
   `/tmp`-built upstream `samtools` binary, same pattern the
   vcftools parity suite already uses.

Coverage target: ≥85% statements on `pkg/htsgo/cram/...`, matching
the rest of htsgo.

---

## 4. Compliance strategy

CRAM is a published spec (CRAM v3.0 = SAMtools spec; v3.1 likewise).
"Compliant" for this project means:

- **Decode**: every container/slice/block/data-series shape the
  spec permits is handled, or rejected with a precise error citing
  the spec section — never a silent wrong answer.
- **Encode**: output is re-readable by upstream `samtools` *and*
  by our own reader, and the embedded reference MD5s verify.
- **Version dialects**: v3.0 and v3.1 are first-class; v2.1 decode
  is best-effort/deferred; v4 is explicitly out of scope.
- **Reference integrity**: MD5 mismatch is always a hard error.
  We never decode against the wrong reference silently.

Every intentional deviation from upstream behaviour gets an entry
in `docs/UPSTREAM_BUGS.md` (fix-on-port) or a documented
limitation here.

---

## 5. Open-question resolutions (from CRAM_DESIGN.md §"Open questions")

1. **rANS port-vs-dep** → **port, in-tree pure Go.** The dep was
   pre-approved but pure-Go is cleaner and keeps the
   stdlib-only claim for all non-LZMA CRAM. (Revises CRAM_DESIGN.md.)
2. **Specific rANS library** → **N/A**, we port (see 1).
3. **LZMA library** → **`ulikunitz/xz`**, confined to
   `codec/lzma.go`, landing in its own PR (C-LZMA) only when a
   fixture needs it.
4. **v2.1 decode required for v1?** → **Deferred.** v3.0 + v3.1
   cover essentially all CRAM produced this decade. v2.1 is a
   fast-follow if archive data demands it.
5. **EBI online-ref fallback default** → **off in v1.** Network
   fetches only happen with an explicit `--use-cram-online-ref`.

---

## 6. Progress log

- **C1** — landed (#160). rANS 4x8 order-0 + order-1, decode and
  encode, byte-exact against the `r4x8` compliance vectors.
- **C2** — landed (#161). rANS 4x16 order-0 decode + encode in
  `pkg/htsgo/cram/codec/rans4x16.go`, byte-exact against the
  `r4x16/*.0` vectors (q4, q8, qvar, q40+dir) for both decode and
  encode. Includes the framing format byte, the big-endian varint
  size field, the X_CAT store-uncompressed fallback, and a decoder
  fuzz target.
- **C2.1** — landed. rANS 4x16 order-1 context model in
  `pkg/htsgo/cram/codec/rans4x16_o1.go`, byte-exact against the
  `r4x16/*.1` vectors for both decode and encode. Includes the
  10/12-bit table-precision auto-tune (`rans_compute_shift`, with the
  `fast_log` bit-trick replicated exactly), the order-1 frequency
  table with run-length-encoded zero gaps, and the optional
  rANS-O0-compressed frequency header. The PACK/RLE/STRIPE/X32
  transforms remain rejected with a clear error — split out to C2.2.
