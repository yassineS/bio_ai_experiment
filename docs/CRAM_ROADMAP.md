# CRAM execution roadmap

Companion to [`CRAM_DESIGN.md`](CRAM_DESIGN.md) (the *why* and the
up-front decisions). This doc is the *how*: the shopping plan, the
PR-by-PR execution sequence, and the testing + compliance strategy.

**Status: read/write + all block codecs + lossy quality binning
complete (C1–C11, C-LZMA, C-Arith, C-FQZComp, C-NameTok landed); the
C-EmbedRef closure pass landed embedded-reference decode, the `cF`
internal-tag strip, the htslib RG/PG aux-order rule, and a
feature-ordering fix that unblocked v2.1 decode.**
Pure-Go CRAM v3.0/v3.1 read and write, reference resolution, `.crai`
index read/write, and samtools-CLI integration are all merged, every
CRAM block compression method (0–8 — raw, gzip, bzip2, LZMA, rANS
4x8/4x16, the arith_dynamic range coder, fqzcomp and the name
tokeniser) is implemented, and opt-in lossy quality-score binning on
the encode side is landed. The CRAM roadmap is complete.

**Documented decode remainder** (small, precise, all behind clear
errors — never a silent wrong answer):

- **v2.1 slice-header record counter — CLOSED.** htslib reads a v2
  slice header's record-counter field as a 32-bit varint (ITF-8) and a
  v3 one as 64-bit (LTF-8) (`cram/cram_decode.c`,
  `cram_decode_slice_header`). The reader now threads the container's
  CRAM major version through `ParseDataContainer` →
  `parseSliceHeader(p, major)` so a v2 slice reads the counter as
  ITF-8 and a v3+ slice as LTF-8, matching htslib exactly. The
  `Container` type gained a `Major` field (populated by `Reader.Next`;
  a hand-built zero-value container is treated as v3+, preserving the
  historical LTF-8 default), and the `.crai` builder threads the file
  definition's major into its own slice-header parse. The two
  encodings coincide for every value < 2^28, so realistic v2.1 files
  already decoded byte-exactly; the fix additionally handles a
  record-counter ≥ 2^28 (a v2.1 file with ≥ ~268 M reads in the
  slices preceding the one being read), where the encodings diverge.
  Validation: the live-samtools v2.1 parity round-trip
  (`TestEmbeddedReferenceParity/v2.1`, `TestV21RecordCounterParity`)
  plus a focused unit test (`TestParseSliceHeaderRecordCounterWidth`)
  proving the ITF-8 vs LTF-8 branch decodes the counter and keeps
  every trailing field aligned for v2 vs v3 at the 2^28 boundary.
- **Network REF_PATH / EBI URL fetch.** An unresolvable reference is a
  clear error naming the missing MD5; the sandbox does no network
  reference fetch (CRAM_DESIGN §"Reference resolution"). Embedded
  reference, an explicit `--reference` FASTA, and the local REF_CACHE
  directory all work.
- **X_EXT (bzip2) *encode*** for the arith / name-tokeniser codecs —
  **DEFERRED (documented), not a correctness gap.** Decode works via
  stdlib `compress/bzip2`. Encode is unsupported because Go has no
  standard-library bzip2 *encoder* and none is sanctioned (CLAUDE.md
  permits exactly one CRAM third-party dep, `ulikunitz/xz`, for LZMA
  decode only). A correct in-tree bzip2 encoder is a large port — the
  full bzip2 pipeline (BWT with suffix-array sorting, MTF, the two RLE
  stages, and multi-table Huffman with selector coding) is roughly
  1.5–2.5 kLOC of carefully-tested code — for a *rare optional* codec
  that this writer never emits: `chooseBlockCompression` only ever
  selects raw / gzip / rANS-4x16, and the encoders never auto-select
  the X_EXT order bit. X_EXT encode is reached only if a caller
  explicitly requests it, and then returns a clear error
  (`arith: X_EXT (bzip2) encode is unsupported …`) — never silent
  wrong output (verified: the encoders auto-select X_CAT/raw, not
  X_EXT). Implementing it is not in scope; it would warrant its own
  conversation and dependency review if a real workload ever needs it.
- **CRAM v4.0** is out of scope (spec not finalised).

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
| **C2.2** | The v3.1 transform bits — X_PACK (bit-packing), X_RLE, X_STRIPE. | Byte-exact vs the `r4x16/*.{8,9,64,65,128,129,192,193}` vectors. |
| **C2.3** | X_32 — the 32-way unrolled rANS core (a distinct on-wire format from the 4-way coder). | Byte-exact vs the `r4x16/*.{4,5}` vectors. |
| **C3** | `pkg/htsgo/cram` container parser — file def, container, compression header, slice header, block. No data-series decode yet; just the tree walk + per-block decompress dispatch. | Parse + walk every container in a real v3.0 CRAM without error. |
| **C4a** | CRAM v3.0 **read**, part 1 — the data-series encoding layer: the encoding zoo (external / byte_array_stop / byte_array_len / huffman / beta / subexp / gamma / golomb / golomb-rice), the CORE-block bit reader, and the compression-header / slice-header parsers. | Decode every data series of a real v3.0 CRAM's slices without error. |
| **C4b** | CRAM v3.0 **read**, part 2 — record reconstruction from the data series, the reference-diff decode, SAM record emission. | Decode a v3.0 CRAM to SAM records matching `samtools view` output. |
| **C5** | Reference resolution — `--reference`, `REF_PATH`, `REF_CACHE`, MD5 verify; `.crai` index read. | MD5-mismatch surfaced as a clear error; `.crai` region query works. |
| **C6** | CLI plumbing — `iohelper` CRAM magic-byte autodetect so `samtools view/depth/fastq/mpileup` accept CRAM input transparently. | The samtools subcommands round-trip a CRAM fixture. |
| **C7** | CRAM v3.1 read — wire rANS 4x16 + 3.1 edge cases. | v3.1 fixtures decode. |
| **C8** | CRAM v3.0 **write** — encoder, container/slice assembly, codec selection. | Written CRAM re-reads byte-identically through our own reader; `samtools view` reads it. |
| **C9** | `samtools view --output-fmt cram` + `samtools index` `.crai` write. | End-to-end SAM→CRAM→SAM via the CLI. |
| **C10** | CRAM v3.1 write. | v3.1 round-trip. |
| **C11** | Lossy quality binning (encode). | Documented, opt-in. **Landed.** |
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
- **C2.2** — landed. The v3.1 transform layer (X_PACK bit-packing,
  X_RLE run-length, X_STRIPE N-way transpose) in
  `pkg/htsgo/cram/codec/rans4x16_transform.go`, byte-exact against the
  `r4x16/*.{8,9,64,65,128,129,192,193}` vectors for both decode and
  encode. Includes the X_NOSZ size-threading that STRIPE sub-streams
  need and the per-stripe brute-force method search. X_32 (the 32-way
  unrolled core) remains rejected with a clear error — split out to
  C2.3. This slice also fixed a latent C2.1 bug: the order-1
  per-context total must be the exact row sum, with the order-0
  alphabet tracked separately, to stay byte-exact under the transforms.
- **C2.3** — landed. The rANS Nx16 32-way coder in
  `pkg/htsgo/cram/codec/rans4x16_32.go` (the scalar `NX==32` path of
  htscodecs' rANS_static32x16pr.c), byte-exact against the
  `r4x16/*.{4,5}` vectors for both decode and encode. The X_32 format
  bit now selects the 32-state core; it composes with PACK/RLE/STRIPE.
  The order-1 32-way encoder reuses the shared `encodeFreq1RANS4x16`
  with `nway=32`. This completes the rANS 4x16 / Nx16 codec family.
- **C3** — landed. The CRAM container parser in `pkg/htsgo/cram/`
  (`filedef.go`, `container.go`, `block.go`, `itf8.go`, `reader.go`):
  the file-definition / container / block tree walk with ITF-8 and
  LTF-8 integer decoding and per-block decompression dispatch (raw,
  gzip, bzip2, rANS 4x8, rANS 4x16). Every container and block CRC32
  is validated during the parse, which doubles as the compliance
  oracle. Verified against the samtools CRAM v3.0/3.1 fixtures
  (256 blocks, all CRCs valid); truncated input errors cleanly. The
  alignment data-series decode is a later slice. Note: the container
  `length` field is a fixed int32, not ITF-8 as some spec prose
  implies — htslib writes it fixed-width so it can be back-patched.
- **C4a** — landed. The CRAM v3.0 data-series encoding layer in
  `pkg/htsgo/cram/` (`bitreader.go`, `encoding.go`, `huffman.go`,
  `codes.go`, `decode.go`, `compheader.go`, `sliceheader.go`,
  `series.go`, `slice.go`): the MSB-first CORE-block bit reader, the
  full encoding zoo (NULL, EXTERNAL, GOLOMB, HUFFMAN, BYTE_ARRAY_LEN,
  BYTE_ARRAY_STOP, BETA, SUBEXP, GOLOMB_RICE, GAMMA), and the
  compression-header (preservation / data-series / tag maps) and
  slice-header parsers. `ParseDataContainer` exposes a container's
  compression header and per-slice CORE + external block set;
  `Slice.DrainSeries` / `DrainTag` decode a self-delimiting series in
  full. Verified against the samtools v3.0 fixtures
  (`test_input_1_a.cram`, `7.quickcheck.cram30.ok.cram`): every
  EXTERNAL / BYTE_ARRAY_STOP / BYTE_ARRAY_LEN data series and tag of
  every slice drains byte-exactly, and the BETA-encoded `AP` series
  decodes from the CORE stream. Two judgement calls: (1) CORE-bitstream
  series (HUFFMAN/BETA/…) share one interleaved stream and cannot be
  isolated per-series without C4b's record traversal, so C4a verifies
  them structurally and drains only the self-delimiting EXTERNAL family;
  (2) a BYTE_ARRAY_LEN whose length and values sub-encodings name the
  same content id stores length and value bytes interleaved (not
  lengths-up-front) — both layouts are handled. Fuzz targets on the
  compression- and slice-header parsers; ≥90% coverage.
- **C4b** — landed. CRAM v3.0 record reconstruction in `pkg/htsgo/cram/`
  (`record.go`, `feature.go`, `reconstruct.go`, `tag.go`,
  `iterator.go`): the per-record CORE/data-series traversal, the 12
  read-feature codes, feature-list → SEQ/QUAL/CIGAR reconstruction,
  auxiliary-tag decode, downstream-mate resolution, and the
  `RecordReader` API (`OpenRecords`, `Read`/`ReadAll`, `WriteSAM`,
  yielding `sam.Record`). Verified against `test_input_1_a.cram`: 14 of
  its 15 records decode byte-identically to `test_input_1_a.sam`
  (QNAME/FLAG/RNAME/POS/MAPQ/CIGAR/mate fields/SEQ/QUAL/tags, including
  complex CIGARs and `B`-array tags); the 15th is the unmapped read
  `u1`, which CRAM stores without MAPQ or CIGAR — the decoded `0`/`*`
  is spec-correct and `test_input_1_a.sam` is the pre-encode input.
  Reference-backed CRAM (`7.quickcheck.cram30.ok`) decodes structurally
  with `NeedsReference()` set; full external-reference resolution is C5.
  Fuzz target on the record reader; ≥86% coverage on the new code.
- **C5** — landed. Reference resolution + `.crai` index read in
  `pkg/htsgo/cram/` (`reference.go`, `refresolve.go`, `subst.go`,
  `crai.go`). `RecordReader.SetReference` / `SetReferenceFASTA` attach an
  indexed FASTA (via `pkg/htsgo/fasta` faidx random access, memoising one
  contig); `SetRefCache` / `UseRefCacheFromEnv` attach the htslib
  REF_CACHE local cache, keyed on the contig's `@SQ` `M5` tag with the
  `%2s/%2s/%s` path layout. The record decoder now reconstructs mapped
  reads from the reference — match runs copy the reference span,
  substitution features resolve through the SM substitution matrix
  (`reconstruct.go` gained a reference cursor alongside the read cursor)
  — and verifies each slice header's reference MD5 (md5 of the
  upper-cased, base-only span); a mismatch is a hard error. The network
  REF_PATH URL fetch is out of scope: an unresolvable reference is a
  clear error naming the missing MD5. `ReadCRAI` / `OpenCRAI` parse the
  gzip-compressed `.crai` TSV; `CRAIIndex.Query` / `QueryRegion` return
  the slice entries overlapping a region. Verified against
  `7.quickcheck.cram30.ok.cram` + `dat/mpileup.ref.fa` (contig `17`): the
  reference-backed decode drops the 'N' count from 56532 to ~0, every
  base is valid, and pure-match reads equal the reference span; the
  REF_CACHE path round-trips via the `@SQ M5` key; the real
  `mpileup/ce#5b.cram.crai` parses. Fuzz target on the `.crai` parser;
  ≥85% coverage on the new code. Two judgement calls: (1) the REF_CACHE
  key is the `@SQ` `M5` tag (the whole-sequence digest), not the slice
  header MD5 (the slice-span digest) — matching htslib; the slice MD5
  still verifies the span after the whole sequence is fetched. (2) An
  unmapped (`-1`) or multi-reference (`-2`) slice resolves to a nil span
  and the per-record reference-derived bases fall back to the 'N' fill —
  a single-reference slice is the fully-resolved path.
- **C6** — landed. CRAM input wired into the samtools CLI. The samtools
  subcommands all open alignment input through one seam, `sam.NewReader`;
  the CRAM router cannot live in `pkg/htsgo/sam` (that would cycle, as
  `cram` imports `sam`), so it lives in `pkg/htsgo/alnio`, which may
  import both. `iohelper.DetectFormat` sniffs the `CRAM`/BAM/SAM magic;
  `alnio.NewReader` returns a `sam.Reader` — and `*cram.RecordReader`
  already satisfies that interface directly (matching `Header()` /
  `Read()`), so no adapter is needed. `view`, `depth`, `fastq` and
  `mpileup` now auto-detect and accept CRAM with no flag; `view -T`
  threads a reference FASTA into the CRAM decoder. Verified against
  `test_input_1_a.cram`. No import cycle (`go list -deps` confirmed).
- **C7** — landed. CRAM v3.1 read. v3.1 shares the v3.0 container and
  record format and differs only in the available block codecs; the
  rANS 4x16 family (C2–C2.3) was already wired into the block
  decompression dispatch, so the C3–C5 reader decodes v3.1 with no new
  production code. C7 locks that in: `pkg/htsgo/cram/v31_test.go`
  verifies the v3.1 fixture `cram_size/mpileup.1.cram` is recognised as
  CRAM 3.1, decodes to its full record set, genuinely exercises the
  rANS 4x16 block codec (every block decompresses to its declared
  size). The arith_dynamic range coder (method 6), the fqzcomp
  quality-score codec (method 7) and the name tokeniser (method 8)
  have since all landed in `pkg/htsgo/cram/codec/`. The name
  tokeniser (`nametok.go` / `nametok_encode.go`) decodes byte-exact
  against the full htscodecs tok3 vector set (110 vectors, 11 corpora
  x 10 levels); its encoder round-trips at every level. Encode is not
  byte-exact against the committed tok3 vectors — those were produced
  by an older htscodecs whose rANS 4x16 encoder did not apply the PACK
  transform in some token-block paths, so the current sub-codec emits
  different (often smaller) streams; the htscodecs test harness itself
  only ever decodes those vectors, never re-encodes to match them. All
  CRAM block compression methods (0-8) are now supported.
- **C8** — landed. The CRAM v3.0 writer (`writer.go`, `writeencode.go`,
  `writeheader.go`, `writefeature.go`, `writetag.go`): a `RecordWriter`
  inverting the reader — `sam.Record`s → data series → blocks → slices
  → containers, with the file definition, SAM-header container, ITF-8
  encoders and CRC32 assembly. Deliberately simple and correct over
  compact: reference-free (preservation `RR=false`, a mapped read's
  bases stored literally as `b` read-features per match run), an
  all-EXTERNAL codec set (no CORE bitstream / Huffman), gzip-or-raw
  blocks, one slice per container capped at 10000 records, every
  record detached so the mate fields round-trip without the
  downstream-mate optimisation. Oracle: `test_input_1_a.cram` → our
  reader → our writer → our reader yields equal records; round-trip
  tests cover mapped/unmapped reads, every CIGAR op, all aux types
  including `B`-arrays, mate pairs (same- and cross-ref), multi-
  reference slices, multi-container files and empty input. Unencodable
  shapes (CIGAR/SEQ length mismatch, unknown reference, bad aux type)
  are rejected with a clear error. `FuzzRecordWriter`; ≥89% coverage
  on the new code. One judgement call: tag values use `BYTE_ARRAY_LEN`
  (length-prefixed), not `BYTE_ARRAY_STOP` — a fixed-width binary tag
  value can contain any byte, so a stop delimiter is ambiguous.
- **C9** — landed. CRAM output from the CLI and `.crai` index write.
  `samtools view -C` / `--output-fmt cram` now writes CRAM: a
  `cramWriter` adapter in `pkg/htsgo/alnio` bridges `cram.RecordWriter`
  (header at construction) to the `sam.Writer` interface (separate
  `WriteHeader`), symmetric with C6's reader adapter — `alnio` imports
  both `sam` and `cram`, `sam` imports neither, no cycle. `samtools
  index` detects a CRAM input via `iohelper.DetectFormat` and writes a
  `.crai` instead of a `.bai`: `pkg/htsgo/cram/craiwrite.go`
  (`WriteCRAI`, the inverse of the C5 `ReadCRAI`) and `craibuild.go`
  (`BuildCRAI` — an offset-tracking CRAM walk emitting one entry per
  slice). Oracle: a SAM→CRAM→SAM round-trip through the `view` code
  path recovers the records, and a built `.crai` parses back and its
  overlap query returns the expected slices. `FuzzCRAIWriteRoundTrip`.
- **C10** — landed. CRAM v3.1 write. A `Version` type (`VersionV30`
  default, `VersionV31`) selects the output dialect; new constructors
  `NewRecordWriterVersion` / `CreateCRAMVersion` / `WriteCRAMVersion`
  carry it, and the existing v3.0 constructors delegate unchanged. The
  package-level mutable `writerVersion` var is gone — version is a
  per-`RecordWriter` field, so concurrent writers of different
  versions are safe. `chooseBlockCompression` is now per-version: v3.0
  keeps the raw/gzip candidate set; v3.1 additionally tries rANS 4x16
  (`codec.RANS4x16Encode`, block method 5) and keeps the smallest, raw
  always in the running so output never exceeds raw. Method 5 is never
  emitted for v3.0. Verified: a v3.1-written file has file definition
  3.1, genuinely contains rANS-4x16 blocks, and round-trips through
  the reader; a v3.0 file stays 3.0 with no method-5 block.
  `FuzzRecordWriterV31`. This completes the C1–C10 CRAM roadmap:
  pure-Go CRAM v3.0/v3.1 read and write, CLI-integrated.
- **C11** — landed. Opt-in lossy quality-score binning on the encode
  side. `pkg/htsgo/cram/qualbin.go` adds a `QualityBinning` enum —
  `BinningNone` (default), `BinningIllumina8`, `BinningIllumina4`,
  `BinningIllumina2` — each carrying a dense `[256]byte` lookup table.
  The tables are the standard Illumina recalibration schemes from the
  "Reducing Whole-Genome Data Storage Footprint" technical note: the
  canonical 8-level table (0-2→0, 3-9→6, 10-19→15, 20-24→22, 25-29→27,
  30-34→33, 35-39→37, 40+→40), the coarser 4-level table (0-9→0,
  10-19→15, 20-29→25, 30+→37) and the 2-level NovaSeq-style table
  (0-14→6, 15+→37). `QualityBinning.BinQuality` maps a quality slice
  through the table, returning a fresh slice — the caller's
  `*sam.Record` is never mutated — and the SAM no-quality sentinel
  `0xff` is passed through unchanged. The writer option mirrors the C10
  `Version` pattern: a `WriterOptions{Version, Binning}` struct and a
  `NewRecordWriterOpts` constructor (the other constructors delegate to
  it, zero-value unchanged), so the default writer is byte-for-byte
  identical to before — verified by a dedicated default-bytes-unchanged
  test. When a real scheme is set the writer maps each record's QUAL
  through the table just before it reaches the `b.qs` series, and
  appends a `@CO` provenance line to a *copy* of the embedded SAM
  header noting the lossy transform. CLI: the
  `samtools view --output-fmt-option qbin=8|4|2|none` flag (upstream's
  KEY=VALUE form) threads through `ViewOptions.CRAMQualityBinning` →
  `alnio.NewCRAMWriterOpts`; it is a no-op for SAM/BAM output. Oracle:
  round-trip with `BinningIllumina8`
  decodes QUAL equal to the binned input, `BinningNone` is exactly
  lossless, and a CLI `view -C` with/without the flag bins or leaves
  qualities verbatim. This completes the CRAM roadmap.
- **C-LZMA** — landed. LZMA block decompression (`codec/lzma.go`,
  `LZMADecode`/`LZMAEncode`). CRAM method-3 blocks are a complete `.xz`
  container stream (the `\xFD7zXZ\x00` magic) — htslib's CRAM writer
  uses liblzma's `lzma_easy_buffer_encode` — read via
  `github.com/ulikunitz/xz`, the one sanctioned CRAM third-party dep,
  imported only from `codec/lzma.go`. `block.go`'s method-3 dispatch
  now decodes instead of erroring. No LZMA CRAM fixture exists in the
  reference corpus, so an `.xz` round-trip is the oracle; a 1 GiB
  decompression ceiling and a fuzz target guard malformed input.
- **C-Arith** — landed. The htscodecs `arith_dynamic` adaptive range
  coder (CRAM block compression method 6) in `codec/arith.go` +
  `codec/arith_transform.go`: the carry-less byte-wise range coder,
  the order-0/order-1 adaptive frequency models, the RLE cores, and
  the X_PACK/X_RLE/X_STRIPE/X_CAT/X_EXT transform layer (which reuses
  the rANS 4x16 transform/varint/pack helpers — only the entropy core
  differs). `block.go` method-6 dispatch now decodes. Byte-exact
  against all 32 `arith/q*` + `u32` vectors for decode and 31/32 for
  encode — the one gap is `u32.4` (X_EXT, a bzip2 payload): decode
  works via stdlib `compress/bzip2`, but Go has no bzip2 encoder and
  none is in the sanctioned dependency set, so X_EXT *encode* returns
  a clear error. `FuzzArithDecode` (1M+ execs clean).
- **C-EmbedRef** — landed (audit + closure pass). Closed four real
  decode-parity gaps surfaced by comparing our reader against live
  `samtools view` of `samtools` `view -C` output (the SAM
  `dat/test_input_1_a.sam`, which samtools embeds-reference because no
  external reference is available for its M5-less `@SQ` lines):
  1. **Embedded reference.** A slice written with samtools' `embed_ref`
     carries its reference span as a block (the slice header's
     `EmbeddedRefID`). The decoder never consumed it, so every
     reference-match base reconstructed as `N`. `Slice.EmbeddedReference`
     decompresses the block and `RecordReader.resolveSliceReference`
     now prefers it over any external source — self-contained, no
     external FASTA/REF_CACHE needed, and (like htslib) trusted verbatim
     with no MD5 cross-check. This is the single highest-value fix: it
     turned `NNN…` sequences into byte-exact ones.
  2. **`cF` internal tag.** htslib writes a single-byte `cF` ("CRAM
     flags", type `C`) tag into the slice tag dictionary and *strips*
     it on read (`cram_decode.c`, "Remove cF tag"). The reader was
     surfacing it as a spurious `cF:i:` SAM aux. It is now drained from
     the data series but never emitted.
  3. **RG/PG aux order.** A record's `RG` comes from a dedicated CRAM
     data series, not the tag dictionary; htslib emits the dictionary
     tags first and appends `RG` last (unless `RG` is itself in the
     dictionary, in which case it is emitted in place). The reader had
     a heuristic that inserted `RG` *before* `PG`; it now matches
     htslib exactly.
  4. **Feature ordering.** The read-feature → CIGAR reconstruction only
     flushed the pending reference-match run before *read-consuming*
     features, so a deletion or reference skip preceded by a match run
     mis-ordered (e.g. `4M1D` → `1D` ahead of the `4M`). htslib emits
     the match gap before *every* feature; the reconstruction now does
     too.
  Together these unblocked **v2.1 decode** for the realistic case: the
  same fixture round-trips byte-exactly through `samtools view -C
  --output-fmt-option version=2.1`, closing the long-standing "v2.1
  decode deferred" item for everything but the > 2^28-read
  record-counter edge (since closed — see C-V21 below). Validation:
  live-upstream parity tests (`parity_test.go`, the
  `upstreamSamtoolsCram` `sync.Once` builder, `t.Fatalf` never
  `t.Skip`) assert v3.0 and v2.1 embedded-reference decode
  byte-for-byte against `samtools view`, plus a strong unit regression
  for the feature-ordering CIGAR fix and the `cF`/RG-order rules.
- **C-V21** — landed. The CRAM v2.1 slice-header record-counter edge.
  htslib reads the slice record-counter as ITF-8 (32-bit) for CRAM
  major version 2 and LTF-8 (64-bit) for v3+; the reader previously
  always read LTF-8 (the two coincide below 2^28). `parseSliceHeader`
  now takes the container's major version and selects ITF-8 for v2 /
  LTF-8 for v3+, matching `cram/cram_decode.c`. The version is threaded
  cleanly: `Container` gained a `Major` field (set by `Reader.Next`),
  `ParseDataContainer` reads it (a hand-built zero-value container
  defaults to v3+ for backward compatibility), and the `.crai` builder
  threads the file definition's major into its slice-header parse —
  important because every field after the counter (block count, content
  ids, MD5) shifts when the counter width is wrong. Validation:
  `TestParseSliceHeaderRecordCounterWidth` proves the ITF-8/LTF-8
  branch decodes the counter and keeps the trailing fields aligned at
  the 2^28 boundary for v2 vs v3; the live-samtools v2.1 round-trips
  (`TestEmbeddedReferenceParity/v2.1`, `TestV21RecordCounterParity`,
  the latter via a uniquely-named `upstreamSamtoolsCramV21` helper)
  assert byte-for-byte parity against `samtools view` and a monotonic,
  non-negative per-slice record counter. A ≥ 2^28-record fixture is
  impractical to build (it needs ~268 M reads), so the unit test covers
  the divergence directly and the live round-trip covers the realistic
  v2.1 file. This closes the last documented v2.1-decode item.
  Separately, an honest re-assessment of **X_EXT (bzip2) encode**
  confirms it stays deferred: a correct in-tree bzip2 encoder is a
  ~1.5–2.5 kLOC port (BWT + MTF + RLE + multi-table Huffman) for a rare
  optional codec the writer never emits, and no bzip2 encoder is
  sanctioned. The X_EXT encode path errors cleanly and is never
  auto-selected — no silent wrong output.
