# Parity roadmap

**Goal:** **1:1 feature parity** with the upstream tool for every Go port in
this repo. This file is the authoritative gap list per tool.

**Scope (current focus):** the project is **not taking on new tools**. The
full focus is to drive the tools already ported or started to complete parity.
Every item below is about an existing tool; do not add new ports.

**bcftools native-plugin tail — FULLY CLOSED (2026-06-16 wave).** Every
remaining "not supported by the native framework" rejection across all 41
in-tree plugins is gone, each closed and validated against the real upstream
`bcftools` 1.23.1 binary (CLI-to-CLI oracle), with a binary-free `TestUnit*`
layer added so the pure logic is verifiable on a clean/offline checkout:
host-level `-r/-R/-t/-T` region/target selection **and** the
`--regions-overlap`/`--targets-overlap` pos/record/variant modes; `-W/--write-index`;
curly-brace `{t1,t2,t3}` multi-threshold filter expansion; the full `split-vep`
machinery (gene-list/severity/columns-types/all transcript selectors); `setGT`
binomial/random/read-depth (byte-exact via an in-tree `drand48` port); `prune`
LD/r2/RD thresholding + annotation + soft-filter + rand + keep-sites; the full
`fill-tags` surface (`-S` pops, custom `TAG=func(EXPR)`, `-l`, every built-in
tag); the stats-plugin extras (`-o`, indel-stats `-p` PED, trio-stats `-a`,
contrast `-f`); the format/output tails (ad-bias, remove-overlaps, tag2tag
localized-allele, guess-ploidy `-g`, bin file-input); `gvcfz` and `frameshifts`;
trio-dnm3 float DMM/ALM/DNG models (+ a tolerance-aware proximity harness),
streaming targets and `-o`; gtisec arbitrary ploidy; fixref `--use-id`; and
`vrfs` (the pileup/VAF-profile plugin, byte-exact — `mpileup2 LEGACY_MODE`
realignment is a stub upstream). Several upstream bugs were fixed-on-port and
documented in `docs/UPSTREAM_BUGS.md`.

**Recently closed (2026-06-17 samtools bugfix wave).** Five
parity-pipeline-found `samtools` bugs were fixed and each byte-validated
against the live upstream binary (live-oracle parity tests in
`tools/samtools/pkg/samtools/parity_bugfixes_live_oracle_test.go`, plus
binary-free `TestUnit*` for every helper):

- **`depth -q`/`-Q` were swapped.** Upstream `-q`/`--min-BQ` is the per-base
  quality floor and `-Q`/`--min-MQ` is the read MAPQ floor (bam2depth.c
  `min_qual`/`min_mqual`); the CLI mapping is now correct (legacy
  `--min-baseq`/`--min-mapq` retained as aliases), and the per-base count is
  byte-identical to upstream including with `-a`/`-aa`/`-r`/`-b`.
- **`mpileup -a` omitted positions.** Upstream `-a` emits every position
  (1..LN) of every read-bearing contig (leading/interior/trailing fill to
  `sam_hdr_tid2len`), not just the covered extent; `-aa` adds read-free
  contigs. Both now match (the old "covered extent" clamp was the bug — this
  was a genuine bug, not intended behaviour).
- **`sort` coordinate tie-break.** The key is `(refID, pos, reverse-strand)`
  (bam_sort.c `bam1_cmp_core`), with equal records preserving input order via
  the stable sort — not the previous QNAME tie-break.
- **`view -s` subsample.** Now uses upstream's deterministic per-read-name
  hash `Wang(X31(qname) ^ seed) & 0xffffff < frac` with the glibc
  `srand()`/`rand()` seed transform, so the kept set is identical and mates
  stay paired (was a per-record Go RNG draw).
- **`addreplacerg -r`.** Default mode is now `overwrite_all` (matching
  upstream), the `-r` path prunes other `@RG` header lines and replaces an
  existing ID under `-w`, and every record's `RG:Z:` tag is set accordingly.

**Recently closed (2026-06-14 parity wave).** The remaining tractable feature
gaps were closed and parity-validated against the vendored upstream binaries:

- **bedtools family — input formats & block-awareness.** A shared
  `pkg/htsgo/alnbed` (SAM/BAM→BED12, CIGAR blocks→BED12 blocks) gives **BAM/SAM
  input** to `bedgenomecov` (`-ibam`, `-pc`, `-fs`), `bedjaccard`,
  `bedcoverage`, `bedspacing`, `bedgroupby`; **`--split`** to
  `bedjaccard`/`bedgenomecov`/`bedcoverage` — including **`bedcoverage --split`
  over a *blocked query* (`-a`)**: a BED12 record or a spliced (`N`-CIGAR) BAM
  alignment is split into its sub-blocks, overlap is counted only within those
  blocks (introns excluded) while the reported length-of-A and per-base depth
  vector still span the full `[start,end)`, and a B feature straddling an intron
  is counted once per overlapped query block — byte-validated vs upstream 2.31.1
  across default/`--counts`/`--depth`/`--hist`/`--mean` and `-s`/`-S`
  (`TestUpstreamParity_SplitBlockedQuery{BED12,BAM}`, plus binary-free
  `split_unit_test.go`). Two further `bedcoverage` divergences were closed
  (2026-06): (1) under `--split` the overlap-fraction thresholds
  `-f`/`-F`/`-r`/`-e` are now correctly **ignored** — upstream keeps the
  always-populated `BlockMgr` overlapSet, not the fraction-filtered resultSet,
  so any B overlapping any A block is counted whatever the fractions (the
  non-`--split` path is unchanged); (2) BED12 `blockSizes`/`blockStarts` are
  echoed **verbatim** (trailing comma preserved iff present) via new
  `bed.Record.RawBlockSizes`/`RawBlockStarts` retained by the BED reader, with
  BAM-derived records falling back to upstream's trailing-comma form. The `-e`
  (either A OR B) flag was also added. Oracle tests
  `TestUpstreamParity_{FractionUnderSplitIgnored,NonSplitFractionUnchanged,VerbatimBED12BlockEcho}`;
  binary-free `frac_echo_unit_test.go`. **GFF/VCF input** to `bedmerge`
  and **GFF** to `bedmap`. `bedclosest` gained multi-database
  `-b … -names/-filenames/-mdb` and `-s`/`-S`/`-N`. Flag-mapping bugs fixed:
  `bedsubtract -N`, `bedmerge -S`/`--delim`, `bedjaccard -S`, `bedcomplement -L`.
- **bedtools family — fisher/subtract/summary/annotate parity (2026-06-17).**
  Four pipeline-found divergences closed and byte-validated vs upstream 2.31.1
  (each tool gained a `live_parity_test.go` that runs the vendored binary and
  `t.Fatalf`s if it is absent, plus binary-free `TestUnit*`): (1) **bedfisher**
  was under-counting overlaps — the overlap counter binary-searched on `ChromEnd`
  over a start-sorted B slice (invalid: `ChromEnd` is not monotonic there, so a
  long early-starting B was skipped), which skewed n11/n12/n21/n22 and the
  p-values; it now scans every B with `ChromStart < A.End` and tests the end
  per record. Also `-m` now merges **both** inputs (A and B), matching upstream's
  `FileRecordMergeMgr`. (2) **bedsubtract** gained the reciprocal `-r` flag
  (require `-f` of both A and B). (3) **bedsummary** was rewritten to the exact
  upstream 10-column output (chrom_length / chrom_frac_genome / frac_all_ivls /
  frac_all_bp columns, genome-file row ordering, fixed-9 precision, the
  per-data-row trailing tab, the `-1` default row, literal `1.0` on the `all`
  row) and now requires `-g`. (4) **bedannotate** no longer prints a spurious
  header without `-names` and now reproduces upstream's per-chromosome/per-UCSC-bin
  record ordering (the `getBin` function was ported); `-names` is variadic and
  the header pads the `#` with `bedType-1` tabs.
- **bcftools.** `view -s` recomputes INFO/AC/AN (`-I` to suppress); `view -v/-V`
  type selectors; `query` position tokens (`%POS0/%END/%END0/%FIRST_ALT/%IS_TS`);
  three **BCF writer** encoding fixes (missing-value sentinels, GT-missing,
  Flag) — `bcftools` now reads our BCF byte-equivalently.
- **samtools.** mpileup text-path BAQ (`-B`/`-E`); `fastq -T '*'` all-tags and
  QNAME-based pairing (lone mates → `-s`). (depth `-a`/`-b`/`-q`/`-Q`,
  mpileup `-a`, sort tie-break, view `-s`, and addreplacerg `-r` parity were
  closed in the 2026-06-17 bugfix wave above.)

**Nothing is out of scope.** There are no "non-goals." Every remaining `t.Skip`
is a gap to close — including CRAM input/codecs, the live-binary/fixture gates
(vendor the fixtures or build the submodule binaries), perl-backed prinseq
parity, and every documented edge case. The target is 100% parity with bug
fixes, docs, and tests for all of it.

The project's stated goal is to make these tools faster, better tested, and
better documented than their originals — which requires that we actually
implement the same features. "Working subset" is a milestone, not the
destination.

> The companion file [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md) tracks bugs
> we identify in the originals as we go. We don't carry those over; we
> either fix on port or note and skip with a clear pointer.

The companion file [`../analysis/tool_ranking_2026.md`](../analysis/tool_ranking_2026.md)
ranks **next** tools to port. It is **not** a "deprioritise existing tools"
filter — existing tools all get carried to 1:1.

> **Salvage of PR #219 (this branch).** The high-value bcftools/samtools
> parity from the stale, conflicting PR #219 was lifted onto current
> `main` without that PR's ~320k lines of committed test-trace fixtures.
> **Salvaged and validated** (live-oracle byte-parity against the
> vendored upstream binaries, no committed goldens): bcftools `call`
> modes (`-m` / `-c` / `--gvcf` / `-C alleles` / `-G` / `--ploidy
> GRCh37|GRCh38` / `--ploidy-file`), bcftools `mpileup --indels-cns`
> (with the pure-Go `pkg/htsgo/edlib` port), and samtools `phase` (full
> `phase.c` upstream-schema emit + `-l`/`-e` site lists). **Deferred to
> a follow-up** (not in this branch): the new `pkg/bamtobed` +
> `tools/bamtobed` + `tools/bedgenomecov` tools, and the broad PR-#219
> live-oracle suites for the many *other* subcommands (csq, convert,
> roh, gtcheck, view-types, …) whose PR-#219 implementations diverge
> from main's current, separately-validated versions — those tests were
> trimmed to the salvaged surface to avoid regressing main.

## Definition of "1:1"

A tool is **1:1** when:

1. Every subcommand the upstream binary exposes is present in the Go port.
2. Every documented flag/option is recognised (either implemented or
   gracefully rejected with a clear error pointing at this roadmap).
3. For every input shape upstream accepts, our port produces the same
   **logical result**, modulo:
   - documented intentional deviations (recorded in the tool's README +
     this roadmap),
   - upstream bugs we chose to fix (recorded in [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md)),
   - **RNG / stochastic divergence** — see the policy below.
4. The validated-parity test suite (runs the upstream test corpus through
   our port) passes for every supported case, with explicit `t.Skip()` for
   each documented exception.

### RNG / stochastic-output policy (2026-05-15)

For subcommands that use randomness (`bedshuffle`, `bedsample`, `seqtk sample`,
`seqtk randbase`, `samtools view -s subsample`, etc.) the parity bar is:

- **Reproducibility within our tool**: same seed + same input → same output
  byte-for-byte, every time, across Go versions and platforms.
- **Logical equivalence with upstream**: our output must satisfy the same
  invariants upstream's output does (correct sampling fraction, no
  duplicates without replacement, strand filters honoured, etc.).
- **Byte-identical with upstream where upstream is seed-reproducible.** This
  is now achieved by porting upstream's exact RNG in pure Go (a few dozen
  LOC per generator), not by deviating to `math/rand`:
  - `seqtk sample`/`randbase` — glibc `drand48` / krand MT19937-64 port.
  - `bedsample` — `std::mt19937_64` port.
  - `bedshuffle` — `std::mt19937_64` port (mt19937.go) plus the exact
    genome-file-order projection and per-mode draw/retry order from
    `shuffleBed.cpp`; `bedshuffle -seed N` matches `bedtools shuffle -seed N`
    byte-for-byte (default/`-incl`/`-excl`/`-chrom`/`-chromFirst`/
    `-allowBeyondChromEnd`), asserted against the upstream binary.
  - `vcftools --max-indv` — glibc `rand()` + libstdc++ `std::random_shuffle`
    port (glibc_rand.go). Upstream seeds with `srand(time(NULL))` and exposes
    **no** seed flag, so a plain upstream run is genuinely non-reproducible;
    this port adds `--max-indv-seed N` and reproduces the exact selection
    upstream would emit for that seed (verified against a C harness that was
    in turn checked against the real upstream binary at matching epoch
    seconds).

Where upstream's only entropy source is an unseeded wall-clock value with no
user-facing seed (as in `vcftools --max-indv`), byte-matching a *given*
upstream run is impossible by construction; the port's seeded path matches the
algorithm byte-for-byte and is documented as such.

The parity-test infrastructure asserts byte parity directly (the `t.Skip(
"RNG byte-parity ...")` markers are being removed as each generator is ported).
Logical-invariant assertions (e.g. "every shuffled interval has the same
length as the input; every shuffled interval is on a chrom in the genome")
are kept as cheap additional checks.

---

## Current status (2026-06-11, post the ~70-PR wave)

A skimmable per-tool completion table lives in the top-level
[`PROJECT_STATUS.md`](../PROJECT_STATUS.md), which also carries the
**definitive remaining-gap list** and the **non-goals** list. Quick state:

- **Done (1:1):** `seqtk`, `sickle`, `skewer`, `fastp` (100%; poly-G/poly-X
  and sliding-window trimming byte-exact, multi-thread `--split` byte-exact
  vs upstream per thread count, JSON report sub-fields — curves, kmer_count,
  q40, sequencing, PE insert_size — byte/structurally-exact, SE adapter
  auto-detect similarity-bounded and observed byte-identical), `htsfile`.
- **Near done (small tails):** `prinseq` (~95%), `vcftools` (146/146
  flags, ~98%, output-column polish only), `bgzip`/`tabix` (~92%),
  `mosdepth` (~98%; CRAM input now landed), `bedtools` (37 bed* tools,
  no missing subcommands, ~95%).
- **Near done (variant-calling now landed):** `samtools` (~97%) and
  `bcftools` (~96%). The former "boulders" are closed: full multi-allelic
  `call`, mpileup SNP slices 1–4, **both** the legacy `bam2bcf_indel` and
  the `--indels-cns` (edlib) indel callers, the csq haplotype engine
  through **slice 4**, and `convert`'s GEN/HAP/TSV/gVCF modes are all
  implemented and live-oracle validated.

### Performance & memory scalability follow-ups (2026-06-25, from the real-data + large-tier perf)

The H2a real-data GIAB run and the large-tier bench (all **byte-exact** vs
upstream) surfaced memory-scalability gaps — correctness is fine, but a few
cells hold O(file) state where upstream streams. These would OOM at WGS scale:

- **`bcftools norm` buffers the whole VCF — FIXED** (commit `4e53d2f`). Was
  `readVCFAll` + global `sortVariants` → ~9.1 GiB on the HG002 VCF (436×). Now
  streams through a bounded `normWindow` reorder buffer (port of upstream
  `vcfnorm.c`'s `rbuf`): peak RSS **18 MB**, byte-exact across every mode on
  synthetic + the real HG002 VCF (4.0 M rows), and unsorted-input behaviour now
  matches upstream (the earlier global-sort `TestNormOrdersByHeaderContigNot…`
  was superseded — upstream streams cross-contig, it does not header-reorder).
- **`samtools depth` ~27× slower — FIXED** (commit `9766251`). It scanned the
  whole BAM even for `-r chr20`; now an indexed BGZF seek (`.csi`/`.bai`) +
  a lean depth-only decode (`ReadDepthInto`). `depth -a -r chr20` on the real
  BAM: **~72 s → ~10 s** (upstream ~8 s), byte-identical, 26 MB.
- **`samtools view` region→SAM** — fast path landed (BAM→SAM direct serialize
  from raw bytes, ~2.6×); the ~11× RSS on `sam_view_bam2cram` remains a
  secondary optimisation target (lower priority; correct today).
- **`samtools view` indexed region over-read — FIXED** (commit `3a5f489`). A
  chunk-bounded scan opened a fresh `bgzip.Reader` at a non-zero block offset
  but compared its (then stream-relative) virtual offset against the absolute
  BAI chunk end, so the bound fired ~startBlock too late and the scan over-read
  into later chunks, **emitting records many times** (a real `view chr20:1-1M
  chr20:30M-30.1M` over a GIAB BAM: ~26k vs ~20k). Root-caused to BGZF:
  `bgzf.NewReaderAt(r, baseCoff)` now makes a mid-file reader report absolute
  virtual offsets; used at all four chunk-bounded scan sites (view fast+decode,
  merge, depth). Byte-validated vs upstream.
- **`samtools view` multi-region semantics — FIXED** (commit `c1c32d1`). Default
  `view reg1 reg2` now walks each region in command-line order, emitting a
  record once per overlapping region (was a single deduplicated union scan, i.e.
  the `-M` semantics always); `-M` collapses to the dedup/coordinate-ordered
  union. Forward/reversed/overlapping queries byte-identical vs upstream in both
  modes. CRAM now honours `-M` too — it walks the regions merged into
  non-overlapping coordinate-sorted intervals (`mergeResolvedRegions`),
  byte-validated vs upstream for default + `-M`, overlapping + reversed.
- **`bcftools isec` memory — FIXED** (commits `58efa79`, `25186a6`). Buffered a
  whole contig (2.7 GB peak on a 3-contig 2.4 M-record pair); now a streaming
  k-way position-window merge with a record+byte-bounded batch — peak RSS
  **2.7 GB → 128 MB** (and 39× → ~6.5× upstream on the 24-sample large fixture,
  492 → ~80 MB). Byte-identical, and a colocated-`sites.txt` ordering parity bug
  was fixed alongside (commit `a669525`).
- **CRAM compression ~1.87× larger — FIXED** (NEW). The v3.0 block-compression
  chooser only tried raw/gzip; it never offered **rANS 4x8**, the entropy coder
  upstream leans on for CRAM (gzip loses badly to rANS on quality scores and base
  calls). `chooseBlockCompression` now offers the version-appropriate rANS codec
  (4x8 for v3.0, 4x16 for v3.1) at order 0 (trusted) and order 1 (round-trip
  verified, since order-1 doesn't round-trip every degenerate block) and keeps
  the smallest. On the synthetic medium fixture our `view -C` output went from
  **1.87× → 1.019×** upstream (within 2%), decodes byte-identically to upstream's
  own CRAM, and encode stayed *faster* (0.15 s vs 0.25 s). The deflate formats
  (VCF.gz/BCF/BGZF) remain ~6% larger (klauspost vs libdeflate — a deliberate
  speed/ratio trade, see CLAUDE.md); BAM is already smaller than upstream.
- **`bcftools isec` "15× wall" was an artifact — FIXED** (commit `ff46358`).
  `isec -p` writes its projection VCFs + `sites.txt` to disk; the large-tier
  bench ran with `TMPDIR` on the slow Docker bind mount, and our plain-VCF writer
  used the default 4 KiB `bufio` while `sites.txt` was written unbuffered (one
  `fmt.Fprintf` per site) — a per-`write()` syscall storm that the bind mount's
  latency dominated (medium fixture 3.9 s on the mount vs 0.5 s on tmpfs).
  256 KiB buffers on both paths: **3.9 s → 0.46 s** on the slow mount,
  byte-identical. On a fast disk isec is a steady **~2.1–2.5×** across tiers
  (rising gently with sample count), not 15×; the figures use the fast-disk
  numbers.
- **`bcftools isec` multi-sample FORMAT over-decode — FIXED** (G6). isec only
  reads CHROM/POS/REF/ALT/ID and re-emits records unchanged, yet it parsed every
  sample's FORMAT into per-sample maps and rebuilt them on write. The vcf reader
  gained a `KeepRawSamples` mode that keeps the FORMAT+sample columns verbatim
  (`Variant.RawTail`) and re-emits them on write — byte-identical to the map
  round-trip for a well-formed record (`serialize(parse(tail)) == tail`),
  validated by the streaming-vs-in-memory differential test and byte-exact vs
  upstream incl. PL. isec wall **2.0× → 1.33× upstream** (0.56 s → 0.37 s on the
  16-sample fixture).
- **`samtools view` CRAM→BAM decode memory — PARTIALLY FIXED (task #30).**
  Decoding a reference CRAM to BAM peaked at **5.5–5.8× upstream** (~110 MB vs
  ~20 MB on chr20). Landed (byte-identical, real GIAB): a per-record
  read-feature scratch buffer (removes 71 % of decode allocation churn — a CPU
  win), a per-slice streaming iterator (bounds multi-slice-container CRAMs;
  htslib-aligned), the reference window trimmed 8→2 MiB, and a soft
  `debug.SetMemoryLimit` scoped to CRAM `view` at the **measured knee of
  64 MiB**. Peak RSS **5.5–5.8× → 3.83–3.86×** (worst-of-five), **no wall
  regression** (1.44×, at baseline). **≤2× was shown infeasible without a
  deeper change:** a `GOMEMLIMIT`/`GOGC` sweep proved RSS floors at ~70 MB
  regardless of GC pressure — the genuine per-slice working set (a slice's
  decompressed data-series blocks + reconstructed records) — while upstream's C
  decoder is an exceptionally lean ~20 MB.
  **REMAINING (the path to ≤2×) — investigated in task #40, re-scoped:** the
  original "stream the per-slice series blocks" hypothesis was DISPROVEN by
  profiling — the decompressed CORE+external series blocks are only ~2.12 MB per
  slice and are freed each slice, not the holder. The ~70–75 MB peak is the OS
  arena high-water from materialising a whole slice of ~10k **fat Go
  `*sam.Record` (~600 B each: `Seq` as a string, `Qual`/`Cigar` slices, `[]sam.Aux`)**
  versus upstream decoding the same slice into packed `bam1_t` (~50 B). It is
  GC-immune (lowering GOGC does not move the peak — it is a transient burst, not
  heap-goal headroom), and batching the fat materialisation is structurally
  blocked (the CRAM data series are consumed in record order in a single pass,
  `resolveMates` needs the whole slice, and the retained `Seq`/`Aux` escape to
  the consumer so they cannot be pooled). **Reaching ≤2× therefore requires a
  leaner `sam.Record` representation (pack `SEQ` as 4-bit nibbles and `Aux` as
  bytes, à la `bam1_t`)** — a large, cross-cutting change to the shared `sam`
  type used by every tool, tracked as its own follow-up. Task #40 shipped a
  byte-identical decode-churn cleanup (an in-place MD/NM aux splice replacing a
  per-record reallocation, ~38 % fewer decode allocations, marginally faster) but
  peak RSS stays ~3.6×.
  **Task #43 — PARTIAL CLOSE, leaner-record direction TESTED via aux-packing
  (merged `885baec`, byte-exact, real GIAB).** The first slice of the leaner
  representation shipped: a `RawAux` passthrough on the
  `view -b -T <ref> <cram>` path — the CRAM decoder builds each record's aux as a
  raw on-disk BAM aux byte block and the BAM writer emits it verbatim (the eager
  `[]sam.Aux` path is untouched for every other consumer; gated, with a defensive
  lazy materialise). Decoded-record md5 is ours == upstream on both CRAMs and the
  roundtrip is clean, zero new test failures, #42's packed-sort path still green.
  Worst-of-five peak RSS: up_chr20.cram **3.83× → 3.55×** (77.3 → 71.7 MB,
  upstream 20.2), up_small.cram **3.85× → 3.69×** (72.8 → 69.7 MB, upstream 18.9);
  wall **1.54× → 1.49×**. It did NOT reach ≤2×: an A/B proved the ~72 MB peak is
  the GC-immune transient per-record decode-allocation BURST (the #40 root cause),
  not the resident aux fat — packing aux trims only ~6 MB.
  **REMAINING (the path to ≤2×, new follow-up):** the full `bam1_t`-style
  single-allocation packed record (name + cigar + seq + qual + aux in one `data`
  block) is needed to cut the per-slice decode-allocation burst. This is the
  large cross-cutting change (the #40 diagnosis's rejected Option C; ~441
  `.Seq`/`.Qual`/`.Aux` field sites across the tree), tracked as the named
  follow-up.
- **`samtools view -C` CRAM encode wall — FIXED (task #33).** Serial per-contig
  `@SQ M5` (reference MD5) hashing dominated header-write (≈80 % of CPU; whole
  reference genome hashed). Now parallelised across a worker pool, each worker on
  its own independent `fasta.RandomAccess` handle (no shared seek state); the M5
  math and emitted bytes are unchanged (`@SQ M5` byte-identical to upstream and
  to the pre-change build). Apples-to-apples encode wall (the M/MT `embed_ref`
  bail artifact neutralised) went **1.43× slower → 0.37–0.46× (faster than
  upstream)**. **Two separate, out-of-scope gaps surfaced:**
  (1) **CRAM size** — ours-written external-reference CRAM is ~12.8 % larger than
  upstream on real GIAB (e.g. 6.23 MB vs 5.52 MB). This is `chooseBlockCompression`
  codec selection (we offer raw / gzip-L7 / rANS-4x8; upstream's default profile
  and its rANS order-1 selection differ), NOT M5 — a size-parity follow-up,
  separate from the wall fix. (2) **`embed_ref` auto-fallback** — upstream, on a
  reference whose first contig name mismatches the FASTA (GIAB `M` vs hs37d5
  `MT`), bails the whole M5 loop and switches to an embedded-reference CRAM
  (smaller, no external `M5`); ours keeps hashing the matching contigs and writes
  a reference-based CRAM. Replicating upstream's auto-fallback would *change the
  emitted file* (embedded ref, no `M5`), so it is a deliberate behaviour-alignment
  decision, not a perf fix — tracked here, not folded into #33.
- **`samtools mpileup` near-indel depth — LARGELY FIXED (task #36).** Deletion
  `*` / ref-skip `>`/`<` placeholders bypassed the `-Q`/`--min-BQ` filter, so
  depth was over-counted at indel flanks (region 20: 542 positions `-B`, 3200
  BAQ-on, ours always higher). Now the filter applies to every pileup entry (the
  placeholders already carry the post-gap base's quality, matching upstream's
  `bam_get_qual[p->qpos]`): **−B 542 → 7 positions, BAQ-on 3200 → 2**, the
  `20:126156` case byte-exact, region 20:30–31 Mb + `-Q 0` 0-diff, no new test
  failures, RSS unchanged. **The residual 7/2 positions** were a separate,
  pre-existing **overlap-removal streaming-order** divergence: `applyOverlapRemoval`
  (`mpileup_overlap.go`) de-weighted mate-pair overlaps EAGERLY over the whole
  contig before accumulation, whereas htslib's `bam_plp` applies each pair's
  `tweak_overlap_quality` incrementally as the later mate is pushed (at its ref
  start), and the tweak both sums and zeroes qualities — so a deletion `*`
  borrowing the post-gap base's quality needs an original-vs-tweaked value that
  varies per column (a single `visibleFrom` threshold cannot express it).
  Follow-ups #41 and #44 both confirmed (with trace-level evidence) that the
  buffered accumulate-then-render model could not reproduce htslib's
  borrowed-quality timing — htslib's semantics (`sam.c tweak_overlap_quality`,
  `p->qpos = s->y`) are an in-place mutate whose per-column visibility is coupled
  to `bam_plp`'s streaming push/emit order, so whether a `*` borrows the pre- or
  post-tweak quality is DATA-DEPENDENT on read interleaving (trace-proven:
  20:17617199 resolves pre-tweak; 20:60709531 and 20:36570298 post-tweak). No
  per-pair / per-record rule expresses it; #44's trial heuristic closed 4 of 7
  windows but regressed 3, so it was reverted. The only faithful path was a true
  streaming `bam_plp` push-then-emit engine.
  **CLOSED by task #46 (merged a04f952).** The streaming engine now exists
  (`mpileup_cursor.go`): it pushes each read at its ref start, fires the read's
  overlap tweak at the push point, advances a column ring buffer, emits each
  column the instant a read past it is pushed, and reads the deletion `*`
  borrowed quality LIVE via a lazy `{rec, qpos}` reference — exactly htslib's
  push-vs-cursor timing. It is gated behind an overlaps-active flag, so the
  buffered / `-x` / consensus paths are byte-unchanged. MEASURED on real GIAB
  (bioval), ours vs upstream — **byte-exact GENOME-WIDE, not just contig 20**:
  - contig 20 `-B` **7 → 0** (md5 == upstream over 9,711,557 lines), BAQ-on
    **2 → 0**;
  - contig 21 `-B` **4 → 0**; contig 22 `-B` **10 → 0**, BAQ-on **8 → 0**;
    ours == upstream on every contig tested, and unchanged vs the old code where
    there is no overlap residual (chr1/2/X) — zero output regression.
  - Latent `-O` placeholder bug fixed along the way (del/refskip now print
    `qpos+1`, not `0`, matching upstream).
  - PERF: faster than the OLD buffered code — `mpileup -B -r 20` wall
    **8.6 s → 6.33 s** (0.74× the buffered path, ~2.2× upstream), peak RSS flat
    (122 vs 127 MB); consensus md5 unchanged; `go test` 0 new failures.
  The streaming-`bam_plp`-engine follow-up (formerly the open #44/#41 item) is
  **DONE**.
- **`bcftools mpileup -r` BCF path OOM — FIXED (task #47, merged 8ff7ba5 /
  2517dc5).** Discovered during task #46: the BCF/genotype-likelihood `-r` path
  did NOT seek the index — on a 10 kb region it scanned the whole file and OOMed
  (~11.4 GB, exit -9), the SAME class of bug as the consensus `-r` OOM fixed in
  #45 and the text-pileup `mpileup -r` index-seek landed in 295a98a / c3785be.
  `bcftools.MpileupFile` now (1) seeks the BAI/CSI index for `-r` regions (a
  bcftools-local region reader on the shared `pkg/htsgo/bam.UnionChunks`, the same
  fast path as samtools `mpileup -r` 295a98a / `consensus -r` #45) and (2)
  materialises the per-contig events array one 1M-position window at a time. The
  change is confined to `tools/bcftools/` (zero samtools/htsgo churn). MEASURED on
  real GIAB (bioval): `mpileup -r 20:30000000-30010000` went **11.4 GB / OOM (exit
  -9) → ~204 MB typical** (worst-of-5 892 MB transient GC), exit 0 (upstream
  bcftools 106 MB), **byte-exact** — BCF/VCF body indexed == linear == upstream
  for plain / `-a AD,DP,SP` / indel, and the per-window materialisation is 0-diff
  indexed-vs-linear across all 1M-window boundary cases. New
  `TestMpileupIndexedRegionMatchesLinear` + the mpileup parity tests green;
  `go test` 0 new failures. Note `samtools mpileup -g` itself was removed upstream
  ("use bcftools mpileup"), so the BCF path is only reached via `bcftools mpileup`.
  **Remaining:** whole-contig `mpileup -r 20` is now bounded (OOM → ~6.4 GB, exit
  0, byte-exact slice) but the whole-contig bucket+BAQ+snapshot phase is still
  resident — that is the SAME per-read-buffer streaming residual class flagged for
  samtools `mpileup -r 20` / `consensus -r 20` below (needs the streaming
  `bam_plp` materialisation), tracked there.
- **`samtools consensus` — discovered residuals (recorded during #44 baseline,
  out of #44 scope).** Several separate, independent gaps surfaced while
  establishing the mpileup baseline and are tracked here so they are not lost
  (the `-r` OOM and the `-f pileup` deletion + insertion-pad running-min
  quality gaps are now fixed; the bayesian `-f fasta` residual remains open):
  - **`consensus -r <region>` OOM — FIXED (task #45, merged 60c81b2).** The
    `-r` path used to ignore the `.bai`/`.csi` index and scan the whole file,
    OOMing (~11.5 GB SIGKILL even for a 1 kb window) on the full GIAB BAM —
    the SAME class of bug already fixed for `mpileup -r` in commit 295a98a.
    `ConsensusFile` now takes the indexed fast path (reuses
    `openBAMRegionReader`/`bam.UnionChunks`, the same machinery `mpileup -r`
    got), so a `-r` region seeks the index chunks instead of draining the
    whole 2.9 GB BAM; it falls back to the unchanged linear scan for
    no-index/CRAM/SAM/unsorted/stdin/`-a`/`-aa`. Two output-invariant engine
    memory optimisations (64 KiB pileupEvent tiling + post-`nmInit` field
    release) keep whole-contig `-r 20` from OOMing. MEASURED on real GIAB
    (bioval), worst-of-5: `consensus -r 20:30000000-30050000 -f fasta --mode
    simple` went **11518 MB / OOM → 23–24.7 MB** (upstream 15.9 MB; ~1.55×,
    was ~725×), md5 byte-exact == upstream (`cmp` IDENTICAL); ours-indexed ==
    ours-linear across 6 regions × 3 modes and a 0-diff whole-chr20 bayesian
    run; new `consensus_indexed_test.go`. **Remaining:** whole-contig
    `consensus -r 20` is now bounded (OOM → ~2.1–2.8 GB, exit 0, byte-exact)
    but still ~20× upstream's 137 MB — that is the SEPARATE per-read Go-record
    buffer, the same streaming residual class flagged for `mpileup -r` under
    295a98a, and is out of scope here.
  - **`consensus -f pileup --mode simple` — 4-position deletion-quality gap —
    FIXED (task #48, merged b3d520d / f6ec27e).** The 4 positions
    (20:30038166, 20:30073764, 20:30075378, 20:30193753) rendered a deletion
    `*` placeholder's per-read quality as the single post-gap base quality,
    whereas upstream (`consensus_pileup.c:195-202`) carries a running minimum
    `MIN(pre-gap base qual, post-gap base qual)`. Fixed by carrying a
    deletion-only `delPileupQual` field rendered ONLY by
    `writeConsensusPileupRow`, so the shared `pileupEvent.qual` used by mpileup
    is untouched (the #46 mpileup byte-exact path is preserved). MEASURED on
    real GIAB (bioval) vs upstream: the 4 documented positions 2-diff → **0**,
    the broad `20:30000000-30200000 -f pileup --mode simple` window
    **0-diff**, and the CDEL running-min rule closed genome-wide. NO-REGRESS:
    `-f fasta --mode simple` md5 `b057be94` unchanged, bayesian `73519a08`
    unchanged, `mpileup -B` vs upstream 0-diff, 0 new test failures.
  - **`consensus -f pileup --mode simple` — insertion-pad running-min —
    FIXED (task #49, merged 84a5809 / 8557aa1).** Upstream's insertion-pad
    running-minimum (`consensus_pileup.c:182-191`, the `p->nth < nth` branch
    for reads whose insertion is shorter than the nth inserted column) is now
    matched: the insertion column pad `*` placeholder quality carries a
    per-read running `MIN(p->qual, b_qual[seq_offset+1])` (0 past read end),
    render-scoped to `callConsensusSimpleInsertions` (no struct field, no
    change to the #48 `delPileupQual` path, the fasta/bayesian call, or
    mpileup). MEASURED on real GIAB (bioval) vs upstream: the residual at
    `20:31032760` (insertion columns nth=1/2/3) 6-diff → **0**, the broad
    `20:30000000-33000000 -f pileup --mode simple` window **0-diff**, the #48
    CDEL control stays 0-diff. NO-REGRESS: `-f fasta --mode simple` md5
    `b057be94` unchanged, bayesian `73519a08` unchanged, `mpileup -B` 0-diff,
    0 new test failures. Together with #48 this closes the `-f pileup --mode
    simple` deletion + insertion running-min quality parity.
  - **`consensus` default Bayesian `-f fasta` — ~1198 structural off-by-one
    diffs**, a separate, larger gap (not the overlap tweak, not the pileup
    running-minimum rules).
- **`samtools sort` wall — FIXED (task #39 then #42, now ≈ parity).** `sort -@4`
  was ~3x upstream at matched memory. Profiling showed `sort.SliceStable`/heap/
  comparator are negligible (<4%); the costs were excessive spilling (57 vs 4
  shards), ~22% GC from per-record allocation churn, and a `-@4` threading
  regression (each spill spawned a fresh parallel-BGZF pool while the sort never
  parallelised). Fixed by a single-threaded spill writer (parallel writer
  reserved for final output) + per-record encode/decode allocation pooling
  (byte-identical): full `-@4` **3.06x → 2.36x**, subset `-@4` **12.25x →
  4.39x**, peak RSS **0.78x** upstream, decoded records byte-identical.
  The interim attempt to scale the spill byte-budget by thread count was
  reverted — a buffered Go `*sam.Record`'s real heap footprint is several times
  its packed size, so it blew peak RSS to 2.01x for no wall gain; the byte budget
  is not a safe RSS lever.
  **CLOSED by task #42 — packed spill.** Profiling showed the decode→re-encode
  round-trip was ~33% of CPU (+ the 22% GC it drove) while `sort.SliceStable` was
  only ~0.5% — so the *parallel in-memory sort* idea was dropped as worthless and
  the **packed-spill format** alone delivered the win. Sort now copies on-disk BAM
  record bytes VERBATIM through buffer → spill → merge → output and decodes only
  the sort key (`sam.BAMReader.ReadRaw` / `sam.BAMWriter.WriteRaw` /
  `FindRawAuxTag`, a reusable raw-record fast path), so the per-record footprint
  ≈ upstream's `bam1_t`. Measured on real GIAB: full `-@4` **2.54× → 1.08×**
  (near parity), subset now *faster* than upstream (`-@1` 0.48×, `-@4` 0.52×,
  chr20 fits RAM and no longer spills), peak RSS **0.53×** upstream, spill shards
  57 → 18 (chr20 2 → 0). Byte-identical across coordinate / `-n` / `-N` / `-t`
  and the `-O sam` fallback — including the BAM `bin` field (0 diff over 39.2 M
  records; the sort path never mutates a record, so the raw passthrough is exact).
  SAM/CRAM input falls back to the decode path via `alnio.ErrNoRawRead`.
- The large-tier heavy cells (`sam_mpileup`, `bcf_call`, `bcf_isec`) were
  mis-reported as OOM; that was **disk exhaustion** (a ~17 GB temp output on a
  5 GB-free overlay), not RAM — all three are bounded and run with a big-disk
  `TMPDIR`. The parity *matrix* harness still buffers both outputs in RAM to
  byte-diff them and wants the stream-comparing approach `realparity` uses.

Genuinely-remaining real gaps (the deliverable — see PROJECT_STATUS.md for
the canonical version with effort sizing). Cloud I/O is **complete**
(streaming + every indexed region path, incl. `.crai` CRAM, live-validated
against S3/GCS/HTTPS). This now includes **AWS Signature Version 2**
(`HTS_S3_V2`, HMAC-SHA1; previously a deliberate non-closure that simply
errored — now implemented and pinned to the canonical AWS worked example)
and full **S3-compatible interop** — custom endpoints (`HTS_S3_HOST`),
`HTS_S3_ADDRESS_STYLE` path/virtual, and the GCS S3-interop XML+HMAC path
(`s3://` against `storage.googleapis.com`, V2 and V4) — so MinIO/Ceph/Wasabi/
Backblaze/GCS all work. **Azure Blob Storage** is now supported too, as a
Go-port extension beyond htslib (which has no native Azure backend):
`az://account/container/blob` and recognised `*.blob.core.windows.net` HTTPS
URLs, with SAS (`AZURE_STORAGE_SAS_TOKEN`), Shared Key (`AZURE_STORAGE_ACCOUNT`
and `AZURE_STORAGE_KEY`, hand-rolled HMAC-SHA256 signed per ranged GET so the
signature covers the `Range` header, pinned to an Azurite-dev-key known-answer
vector), Azure-AD bearer (`AZURE_STORAGE_TOKEN`) and anonymous auth, plus an
`AZURE_STORAGE_BLOB_ENDPOINT` override (see `pkg/htsgo/hfile/README.md`). CRAM `MD`/`NM`
regeneration and the network REF_PATH
fetch are **done**; `gtcheck` `-i/-e` filter expressions and the bedtools
KeyListOps ops are **done**. The former "non-goals" are now all implemented:

- **CRAM v4.0 decode** — uint7-varint + 64-bit positions, version-threaded,
  byte-for-byte parity vs an upstream-written v4.0 CRAM (the v4 transform
  codecs XPACK/XRLE/XDELTA are parsed but not yet decoded — no fixture
  exercises them);
- `gtcheck -c/--cluster` clustering (own design — upstream is an error stub);
- `bcftools convert` PLINK exporters, `bcftools som`, and `samtools tview`
  (text/HTML **and** the interactive `-d C` viewer) are all implemented. The
  `-d C` viewer is pure-Go (Linux raw-mode termios, no ncurses); on non-Linux
  platforms it reports that interactive mode requires Linux.

Recently closed (so no longer in this list): bcftools `concat --ligate`
phased ligation; mendelian2 `sites_not_diploid` (non-diploid records are
now counted and skipped); vcftools BCF-binary I/O
(`--bcf`/`--recode-bcf`/`--diff-bcf`/`--contigs`, roundtrip-tested — the
former htsgo-BCF-writer blocker has landed); the filter engine's
bare-INFO-tag resolution.

Recently closed: samtools `consensus` pileup `-a` placeholder rows (incl.
ref-skip columns, per-nth gap-fill duplication, INT_MIN quality, and
deletion ref-skip boundaries — 120-subtest live oracle); bgzip's `--test`
integrity check (long-only, since `-t` binds to `--threads`).

Implemented (this wave): bcftools `convert` PLINK exporters. Upstream
`vcfconvert.c` leaves the `--plink`/`--tped`/`--bin` option block
**commented out** (lines ~1697–1699) with no implementation, so there is
no upstream binary to diff against. The Go port implements them to the
PLINK1 file-format spec: `-p`/`--plink` (`.ped`+`.map`), `--tped`
(`.tped`+`.tfam`), and `--plink-bed` (binary `.bed`+`.bim`+`.fam`,
SNP-major, `A1=ALT`/`A2=REF`). Multi-allelic and no-ALT records are
skipped with a warning (PLINK is biallelic). Byte-exact `.bed` tests plus
text round-trip tests live in `convert_plink_test.go`; see the per-tool
README for the conventions.

Recently closed (this wave, treat as merged):

- **bcftools**: full multi-allelic `-m` `call`; `--ploidy GRCh37/38` +
  `--ploidy-file`; `convert` `--gvcf2vcf` (and its upstream prefix
  abbreviation `--gvcf`)/`--tsv2vcf`/`--gensample(2vcf)`/
  `--hapsample(2vcf)`/`--haplegendsample(2vcf)`; `annotate`; `consensus`
  chain + iupac; `mendelian2` rule engine; `gtcheck` (PL/cluster minus the
  dendrogram, plus `-O t|z[N]` text/BGZF output and `--n-matches` /
  `--distinctive-sites`); `filter -M/--mask-file`; `cnv --AF-file`; `roh -O z`;
  `query %INFO/<tag>` + `%SAMPLE`; csq **slice 4** (FORMAT/TBCSQ,
  `--unify-chr-names`, `--dump-gff`, `-O b|u|z`); mpileup legacy
  `bam2bcf_indel` **and** `--indels-cns`; BCF `-O u|b` output throughout;
  `-@` output-compression threading.
  > NOTE: a few CLI flag-help strings still read "accepted; v1 not
  > implemented" / "accepted, ignored" (e.g. `filter -M`,
  > `mpileup --indels-cns`); those strings are **stale** — the features are
  > fully wired (verified in code). The help text is a cosmetic follow-up.
- **samtools**: `coverage -A` ASCII histogram; `markdup -d` (optical-dup),
  `-s`/`-S` (stats / supplementary marking); `calmd -C` (cap MAPQ, gated on
  `>10`), `-e` (`=` for match), `-u` (uncompressed BAM); `consensus`
  per-position indel calling + `--het-only`; `phase` `-l/-e` site lists +
  `-b` per-haplotype BAM split with `-A`/`-F` (byte-exact `dump_aln`
  routing via an in-tree glibc-`drand48` port; live-oracle validated);
  mpileup `-g/-u` BCF/genotype-likelihood emit (delegates to the bcftools
  engine); CRAM r/w + X_EXT bzip2 encode; `-@` view/sort/markdup threading.
- **mosdepth**: `-d/--d4`, `-a/--fragment-mode`, `-q/--quantize`,
  `-t/--threads`, and `-m/--use-median` (all byte-identical to the
  upstream v0.3.14 binary; threads produces identical output for any
  count); **CRAM input** (`-f/--fasta`/`--reference` + `REF_CACHE`,
  auto-detected and routed through `pkg/htsgo/alnio`; per-base / regions /
  summary output byte-identical to the equivalent BAM and to upstream
  `mosdepth` on the CRAM). No deferred flags remain.
- **bedtools**: BAM input (`bedintersect`, `bedmulticov`); VCF/GFF input
  (`bedintersect`, `bedmultiinter`); `bedclosest` direction flags
  (`-D`/`-id`/`-iu`/`-fu`/`-fd`/`-t`); `bedintersect -c`; `bedsample` RNG
  byte-parity (std::mt19937_64 port).
- **vcftools**: 0 unsupported flags; `--freq2`/`--counts2` schema complete.
- **CRAM**: X_EXT bzip2 *encode* (in-tree pure-Go); v2.1 decode + v3.0/3.1
  read+write validated byte-for-byte against live `samtools`. Lossy read
  names (`lossy_names=1`) now decode: a detached record's real name is read
  from the mate block, and dropped duplicate names are reconstructed as
  `<prefix>:<n>` exactly as htslib's `cram_to_bam`
  (`TestLossyNamesReadParity`). The writer's `=`/`X` CIGAR handling
  (`--eqx` aligner output, folded to per-base features like M, matching
  htslib) is validated against live `samtools` (`TestEqXWriteParity`); the
  `B` (CIGAR back-step) op is rejected to match htslib's own CRAM encoder.

## Per-tool gap list

Numbers reflect state at 2026-06-09. The per-tool prose below is kept in
sync as gaps close; where a sub-note still reads "accepted but rejected"
for one of the #220–#225 items above, the summary in this header and in
`PROJECT_STATUS.md` is authoritative.

### `seqtk`

**Status:** 24 of 24 upstream subcommands. **1:1 PARITY ACHIEVED.**

Dispatch-table audit at `reference_code/seqtk/seqtk.c::main()`
lines 2099-2122 lists exactly these 24 `stk_*` entry points (`hrun`
and `hpc` are SEPARATE dispatch entries, not aliases — see
`stk_hpc` at seqtk.c:1692 vs `stk_hrun` at seqtk.c:1174): `comp`,
`fqchk`, `hety`, `gc`, `subseq`, `mutfa`, `mergefa`, `mergepe`,
`dropse`, `randbase`, `cutN`, `gap`, `listhet`, `famask`, `trimfq`,
`hrun`, `sample`, `seq`, `kfreq`, `rename`, `split`, `hpc`, `size`,
`telo`. All 24 are now implemented in
`tools/seqtk/cmd/seqtk/main.go`.

Earlier roadmap iterations listed `hpc-bg` as missing — that was a
misreading of the dispatch table. There is no `hpc-bg` subcommand
upstream (verified empirically: `seqtk hpc-bg` is rejected with
`[main] unrecognized command 'hpc-bg'. Abort!`). The genuine missing
entry was `hrun` (homopolymer-RUN finder, BED4 output); the name
collision with `hpc` confused the original audit.

Added this iteration: `listhet`, `hrun`.

- `listhet` — full upstream surface implemented (no flags, positional
  `<in.fa>` only; the Go cmd adds `-o/--output` as the project-wide
  file-redirect convenience). The 1:1-ported algorithm walks every
  FASTA record byte-by-byte and emits a TSV row `name\tpos1based\tbyte`
  for every byte whose `bitcnt_table[seq_nt16_table[b]] == 2` —
  i.e. the 2-base IUPAC codes R, Y, S, W, K, M and their lowercase
  counterparts. The byte is emitted in its original case. Byte-parity
  verified against `reference_code/seqtk` v1.5 on `ambig.fa`,
  `hety_basic.fa`, `hety_lowercase.fa`, and `small.fa` (the latter
  pins the no-output path for fixtures with zero hets).
- `hrun` — full upstream surface implemented (no flags, positional
  `<in.fa> [minLen]`; the Go cmd exposes the minLen knob as
  `-l/--min-len` AND still accepts the upstream positional form for
  compatibility). One BED4 row (`chrom\tstart\tend\tbase`, 0-based
  half-open) is emitted per maximal byte-identical run of length
  `>= minLen` (default 7). Two upstream quirks are mirrored
  byte-for-byte for parity:
  1. Comparison is BYTE-EXACT (no case-fold, no IUPAC fold) so
     `AAaa` is two runs of length 2.
  2. The "open trailing run" flush at `seqtk.c:1200` lives OUTSIDE
     the `kseq_read` loop, so it fires AT MOST ONCE per input,
     using the last record's name and the run state left over from
     that record. If the last record is empty, upstream reads its
     NUL-terminator (`ks->seq.s[0]`) which sets `l = 1`, silently
     swallowing the would-be flush for any `minLen >= 2`. Our port
     reproduces both behaviours; see the trace in
     `tools/seqtk/pkg/seqtk/hrun.go` and the
     `TestHrun_TrailingFlushAcrossRecords` test.

  Byte-parity verified against `reference_code/seqtk` v1.5 on
  `nruns.fa` (default, `-l 3`, `-l 2`) and `hety_basic.fa` (`-l 4`).

Previously added: `kfreq`, `telo`. Both are byte-for-byte parity
ports against `reference_code/seqtk` v1.5 (verified by piping the
hand-built fixtures under `tools/seqtk/testdata/parity/` through both
the upstream binary and the Go port and diffing). Previous iterations
added `fqchk`, `hety`, `famask`, `mergefa`, `rename`, `split`, `size`.

- `kfreq` — full upstream surface implemented: no flags, positional
  `<kmer> <in.fa>`. Per-record TSV row with `name`, length, strand
  (`+` if forward neighbour count strictly exceeds reverse, else
  `-`), neighbour-count, exact-count. The neighbour set is every
  k-mer at Hamming distance ≤ 1 from the target (including the
  target itself). Lowercase ACGT bytes count via `seq_nt6_table`;
  the rolling encoder is reset on any non-ACGT byte. Non-ACGT bytes
  in the k-mer return a typed error (upstream `assert()`s and
  aborts). Byte-parity verified against `reference_code/seqtk` v1.5
  on `kfreq_small.fa` (`AC`, `ACGT`, `AAAA`), `kfreq_edge.fa` (`AC`),
  and `kfreq_mixed.fa` (`AA`, `ACGT`, `CCGG`, `CCCTAA`).
- `telo` — full upstream flag surface implemented:
  `-m STR` (motif, default `CCCTAA`), `-p INT` (penalty, default 1;
  negative values are silently flipped), `-d INT` (max-drop, default
  2000), `-s INT` (min-score, default 300), `-P` (print scoring
  profile instead of BED). The 5' scan walks left-to-right looking
  for motif rotations on the forward strand; the 3' scan walks
  right-to-left querying the same hash set with reverse-complement
  bases, so flipping `-m` to the motif's reverse complement swaps
  which end the BED rows describe. BED rows go to stdout; the
  `<sum_telo>\t<sum_input>` summary line goes to stderr (matching
  upstream byte-for-byte — `Telo` takes separate `stdout` and
  `stderr` writers for this reason). Byte-parity verified against
  `reference_code/seqtk` v1.5 on `telo_basic.fa` (default,
  `-m TTAGGG`, `-P -s 0`), `telo_complex.fa` (default, `-s 100`,
  `-p 2 -d 500`), and `telo_edge.fa` (default).

- `fqchk` — full upstream surface implemented: `-q INT` (default 20,
  use `-q 0` to dump the per-quality distribution). Output covers the
  preamble line, the `POS … avgQ errQ [%low %high | %Qk…]` header, the
  `ALL` aggregate row, and the per-position rows — all matched
  byte-for-byte on `fqchk_mixed.fq` (default, `-q 0`, `-q 30`) and
  `small.fq` (default, `-q 0`). PHRED+33 is hard-coded upstream and
  thus here too.
- `hety` — full upstream surface implemented: `-w INT` (default 50000),
  `-t INT` (default 5), `-m`. Per-position classification (the
  `bitcnt > 2 ? 0 : == 2 ? 2 : 1` map) is reproduced verbatim from
  `seqtk.c:614-619`, so 2-base IUPAC codes (R, Y, S, W, K, M) count as
  heterozygous while 3-/4-base IUPAC codes (B, D, H, V, N, X) do not.
  Byte-parity verified against `reference_code/seqtk` v1.5 on
  `hety_basic.fa` (default, `-w 30`, `-w 30 -t 3`, `-w 30 -m`) and
  `hety_lowercase.fa` (`-w 6 -t 1 -m`). Note: the i==l terminator and
  the upstream rollback at `seqtk.c:604-606` are required to make the
  partial last window emit the same byte count as upstream — both are
  pinned by the parity fixtures.

The previous list also mentioned `seqshuf`, `gcdc`, and `cnregion` —
these are NOT upstream subcommands per the dispatch-table audit.
`cnregion` was dropped in PR #112; `seqshuf` and `gcdc` are likewise
not real entries and should be ignored.

Project-extension policy (PR #113): the previous roadmap entry
implying `pair` was a missing upstream subcommand was wrong —
upstream seqtk v1.5 has no `pair` subcommand at all (verified
against `reference_code/seqtk/seqtk.c::main()` dispatch, which
registers `comp`, `fqchk`, `hety`, `gc`, `subseq`, `mutfa`,
`mergefa`, `mergepe`, `dropse`, `randbase`, `cutN`, `gap`,
`listhet`, `famask`, `trimfq`, `hrun`/`hpc`, `sample`, `seq`,
`kfreq`, `rename`, `split`, `telo`, `size`). Per the 1:1 parity
mandate (and the `cnregion` precedent from PR #112), project-
extension subcommands are NOT shipped under existing tool names;
the project-extension `pair` introduced in the PR #113 first
commit was dropped before merge.

Bugfixes landed since the 2026-05-14 audit:

- `comp` — fixed `seq_nt16_table['U']`/`['u']` mapping. Was 8 (T);
  upstream `reference_code/seqtk/seqtk.c` line 189 defines
  `seq_nt16_table[85] = 15` and `[117] = 15` (i.e. U is treated as
  the 4-base ambiguity N, not as T). The bug caused U bases to count
  in the `#T` column instead of `#4`; the regression test
  `TestParity_Seqtk_Comp_UBaseFasta` now pins this against an
  upstream-generated fixture.

Note: `cnregion` was listed here before but is **not** an upstream seqtk
subcommand (verified against `reference_code/seqtk/seqtk.c` v1.5: the only
`stk_*` functions registered in `main()` are `seq`, `comp`, `sample`,
`subseq`, `mergefa`, `mutfa`, `mergepe`, `randbase`, `hety`, `gc`, `fqchk`,
`hrun`/`hpc`, `listhet`, `famask`, `trimfq`, `hpc-bg`/`hpc`, `seq`, `cutN`,
`gap`, and `kfreq` — no `cnregion`). Dropped from the gap list.

Option-tail gaps (per existing subcommand):

- `comp` — missing `-r REGION` to restrict to a BED region.
- `seq` — missing `-A` (force ASCII output), `-C` (mask sequence with N), `-M FILE` (mask regions), the `-T int` trim option.
- `sample` — full upstream surface ported byte-for-byte (`-s SEED`, `-2`
  two-pass, fraction and fixed-number modes; see the option-tail closure
  below). The byte-parity test is no longer skipped.
- `trimfq` — full upstream surface ported byte-for-byte (Mott default,
  `-q`/`-l`/`-b`/`-e`/`-L`; see the option-tail closure below). The
  byte-parity test is no longer skipped. (There is no `-B` flag upstream.)
- `subseq` — missing the regex-name mode.
- `mutfa` — missing the inverse `--inverse` mode.
- `gap` — full upstream surface implemented (`-l` only). Note: upstream's
  "gap" is any non-ACGT byte (via `seq_nt6_table`), not just literal N,
  so IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) all count.
  We match that byte-for-byte against `reference_code/seqtk` v1.5.
- `gc` — full upstream surface implemented (`-w`, `-f`, `-l`, `-x`). The
  algorithm is the upstream X-dropoff scan, not a sliding window. Output
  is BED4 (`chrom\tstart\tend\thits`). Byte-parity verified against
  `reference_code/seqtk` v1.5 on `gc_small.fa` fixtures.
- `rename` — no upstream flags (positional `<in.fq> [prefix]` only).
  Reproduces upstream's `cpy_kstr` early-return bug at
  `reference_code/seqtk/seqtk.c:1210`: a record without a comment
  inherits the previous record's comment until something with a
  comment overwrites it. This "sticky-comment" quirk is required to
  match upstream byte-for-byte and is covered by a dedicated regression
  test (`TestRename_StickyComment_Quirk`) plus the
  `TestParity_Seqtk_Rename` fixture comparisons.
- `split` — full upstream surface implemented (`-n INT`, `-l INT`,
  positional `<prefix> <in.fa>`). Output files are uncompressed and
  named `<prefix>.<5-digit 1-based>.fa` (literal `.fa` suffix even for
  FASTQ input), matching upstream byte-for-byte across `small.fa`,
  `small.fq`, and the `-l` line-wrap variants.
- `size` — no upstream flags. Emits a single
  `<num_records>\t<total_bases>\n` line, matching upstream
  byte-for-byte across `small.fa`, `small.fq`, `empty.fa`, `nruns.fa`.
- `famask` — no upstream flags (`getopt("")` at
  `reference_code/seqtk/seqtk.c:878`). Output is FASTA wrapped at 60
  columns; the three mask rules (`X` = keep, `x` = lowercase,
  anything else = overwrite) are reproduced 1:1 and verified
  byte-for-byte against `reference_code/seqtk` v1.5 on
  `famask_simple_*` and `famask_*` fixtures.
- `mergefa` — full upstream flag surface implemented:
  `-q INT`, `-i`, `-m`, `-r`, `-h` (`getopt("himrq:")` at
  `reference_code/seqtk/seqtk.c:774`). FASTA-mode outputs and the
  stderr `(same,diff,hom-het,het-hom,het-het)` counter line are
  byte-for-byte against upstream on `mergefa_a.fa`/`mergefa_b.fa`
  for the default, `-i`, `-m`, `-h`, `-q 20` (FASTQ), and the
  60-col wrap fixture. The `-r` path uses Go's `math/rand` (seeded
  with upstream's constant 11 by default); RNG byte-parity is
  explicitly NOT a goal per the policy in the section above.

**Option-tail closure (this iteration).** The remaining per-subcommand
flag gaps to upstream are now closed, each verified live against a
freshly built `reference_code/seqtk` binary (no goldens) by the
`upstreamSeqtkOpts` `sync.Once` builder in
`tools/seqtk/pkg/seqtk/opts_parity_test.go` (byte-for-byte, `t.Fatalf`
never `t.Skip`):

- `seq` — re-ported to upstream's full `stk_seq` flag set
  (`getopt("N12q:l:Q:aACrn:s:f:M:L:cVUX:SF:xR")`, seqtk.c:1392) in
  `tools/seqtk/pkg/seqtk/seq.go`. Now implemented and parity-tested:
  `-A`/`-a` (force FASTA, drop quality), `-C` (drop header comment),
  `-r` (reverse complement), `-l INT` (residues per line),
  `-L INT` (drop sequences shorter than INT), `-q INT`/`-X INT`/`-n CHAR`
  (quality masking to lowercase or CHAR), `-Q INT` (quality shift),
  `-U` (uppercase), `-V` (shift quality by `(-Q)-33`), `-N` (drop
  sequences with ambiguous bases), `-M FILE` (mask BED/name-list regions
  via a 1:1 port of `stk_reg_read`+`stk_mask`), and `-c` (mask the
  complement region). **No `-T` flag exists upstream** (the task's `-T`
  was a phantom — `seq`'s getopt has no `T`); skipped. `-1`/`-2`
  (odd/even read selection), `-s`/`-f` (seq-level fraction sampling) and
  `-S`/`-x`/`-F`/`-R` remain unported and are noted as a residual minor
  gap (not requested in this unit).
- `comp -r FILE` — restrict composition to BED/name-list regions
  (`CompWithRegions` in `comp.go`); each region emits a row whose lead
  columns are `name\tbeg\tend` rather than `name\tlen`, matching
  `stk_comp` (seqtk.c:512-514). Byte-parity tested.
- `sample -2` and `-s SEED` — `stk_sample` (seqtk.c:1228) two-pass mode
  and seeded reservoir sampler are ported byte-for-byte, including a 1:1
  port of upstream's krand MT19937-64 RNG (`sample.go`). This also
  closes the previously-skipped `sample`/RNG parity gap: `SampleN`,
  `SampleFraction`, and the `-2` two-pass path now match upstream
  exactly. `TestParity_Seqtk_Sample_UpstreamByteParity` now asserts
  byte-for-byte (`t.Fatalf`, no `t.Skip`) across fraction, fixed-number,
  two-pass, and FASTA fixtures. (The legacy `Sample`/every-Nth helper is
  retained only for back-compat callers and its structural-invariant test.)
- `randbase` — `RandbaseUpstream` (`mutations.go`) ports `stk_randbase`
  (seqtk.c:525) byte-for-byte, including a 1:1 reimplementation of glibc's
  `drand48` (the 48-bit LCG with the default seed-0 state upstream relies
  on) and the output layout (2-base codes drawn via `m = drand48() < 0.5`,
  3/4-base codes and N passed through, comment dropped, sequence wrapped
  at 60 columns). The CLI default path (no `-s`) uses it; `-s INT` is a
  seedable extension. `TestParity_Seqtk_Randbase_UpstreamByteParity` now
  asserts byte-for-byte (no `t.Skip`), closing the previously-skipped
  `randbase`/RNG parity gap.
- `trimfq -L INT` (retain at most INT bp from the 5'-end) plus the full
  Mott-algorithm trimming and `-q FLOAT`/`-l INT`/`-b INT`/`-e INT`
  paths are ported byte-for-byte (`TrimFQ` in `trimfq.go`, a 1:1 port of
  `stk_trimfq`, seqtk.c:361). This closes the previously-skipped
  `trimfq` parity gap; `TestParity_Seqtk_Trimfq_UpstreamByteParity` now
  asserts byte-for-byte (no `t.Skip`) across the Mott default, the
  window fallback, and the `-b/-e`/`-L` fixed-offset paths. **There is no
  `-B` flag upstream** (trimfq's
  getopt is `"l:q:b:e:L:"`; the task's `-B` was a phantom, the real
  option is lowercase `-b` = trim-from-left); skipped.
- `subseq` regex/name-pattern mode — **not an upstream feature.**
  Upstream `stk_subseq` (`getopt("tl:s")`, seqtk.c:678) does exact
  sequence-name matching against a BED or name-list file only; there is
  no regex mode. The existing `Subseq` already does exact name-list /
  BED matching. Nothing to add; documented as a non-gap.
- `mutfa --inverse` — **not an upstream feature.** Upstream `stk_mutfa`
  (seqtk.c:913) takes no flags at all (positional `<in.fa> <in.snp>`);
  there is no inverse-mask option. Documented as a non-gap.

**Validation:** no upstream-test-suite run yet. The option-tail flags
above are validated live against the built upstream binary (Tier:
live-upstream byte parity, `upstreamSeqtkOpts` builder).

### `prinseq-lite`

**Status:** 1:1 parity. The Go port now covers every behavioural
flag of `prinseq-lite.pl` 0.20.4 that the agreed scope called for,
including `--graph_data` (PR `claude/prinseq-graph-data-land`).
Three subcommands ship: `stats`, `filter`, `graph_data`.

Implemented flags (with the upstream Perl line numbers consulted for
each — `reference_code/prinseq/prinseq-lite.pl`, 0.20.4):

- `--out_format` (1=FASTA, 2=FASTA+QUAL, 3=FASTQ, 4=FASTQ+FASTA,
  5=FASTQ+FASTA+QUAL) — POD spec at lines 242-247; CLI parsing at
  lines 769-789; per-mode write branches at lines 1302-1348,
  3703-3714, 3737-3757.
- `--seq_id_mappings <file>` — POD at lines 293-295; coupling check
  with `--seq_id` at lines 945-948; file open at lines 1350-1358;
  per-record write at lines 3640-3648.
- `--ns_max_p <float>` — POD at lines 344-346; filter check at
  lines 3465-3470 (strict `>` against `(N_count * 100 / length)`).
- `--noniupac` — POD at lines 352-354; filter check at lines
  3478-3481 (`uc($seq) =~ /[^ACGTN]/o`, i.e. case-insensitive).
- `--phred64` — POD at lines 230-232; gate at lines 760-764. This
  is an INPUT encoding toggle (the original roadmap entry called it
  `--phreds`; upstream has no such flag — the only Phred-related
  option is `-phred64`). Implemented as a CLI alias for
  `--qual-type illumina`, so the existing Phred+64 decoder is
  reused. The QUAL output (`--out_format 2/5`) honours the chosen
  encoding when converting to decimal phred scores.

Bundled as a CLI alias only (`--seq_id <prefix>`): documented at
lines 470-472 of the upstream POD, implemented at line 3648. Needed
to make `--seq_id_mappings` useful, since upstream rejects
`--seq_id_mappings` without `--seq_id`.

Behavioural divergences vs upstream that are documented rather than
matched (PRINSEQ-lite has no formal test suite, so byte-for-byte
parity is not enforced):

- Multi-stream outputs (`--out_format 2/4/5`) use the value of
  `--output` directly as the prefix, appending literal `.fasta` /
  `.qual` to it. Upstream derives the prefix from `-out_good`, with
  randomised `_prinseq_good_XXXX` suffixes when no prefix is given;
  we require an explicit `--output` prefix and refuse to stream
  multiple files to stdout. The semantic restriction matches upstream
  (lines 801-802), the on-disk filename layout differs.
- The QUAL output uses the upstream `convertQualArrayToString`
  layout (two-character space-padded decimal, single-space separated,
  wrapped every `LINE_WIDTH=60` values; lines 45 and 2531-2546). The
  upstream `-line_width` knob is now exposed as `--line_width`; it
  overrides the default 60-column wrap for both FASTA and QUAL output
  (0 = no wrap).
- The `--seq_id` rename now **preserves** any trailing whitespace/comment
  from the original FASTA/FASTQ description, matching upstream
  `prinseq-lite.pl:3685-3704` (`$sid.($header ? ' '.$header : '')`). The
  earlier divergence (comment dropped) is **resolved**; see
  `docs/UPSTREAM_BUGS.md > prinseq`.

Graph-data (PR `claude/prinseq-graph-data-land`):

- `--graph_data [FILE]` (POD lines 393-401; emission block at
  lines 2050-2287). Ported in full, including the full stat
  collectors (`getSeqStats`, `getQualStats`, `generateStatsType`,
  `dinucOdds`, `checkForDupl`, `getTagFrequency`, `getBinVal`).
  Exposed as the `prinseq graph_data` subcommand. When `--graph_data`
  is given without a value, the Go port falls back to upstream's
  `<input>__.gd` default (lines 984-987).
- `--graph_stats CODES` (lines 994-1015): the upstream stat
  selector CSV (`ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn`). Supported and
  enforced — unknown codes return an error as upstream does.
- `--qual_noscale` (lines 989-993). **Done.** Wired on the
  `prinseq graph_data` subcommand (`cmd/prinseq/main.go`,
  `--qual_noscale`) through `GraphDataOptions.QualNoscale`. Upstream's
  `$scale` defaults to 1; setting `--qual_noscale` flips it to 0,
  which (a) suppresses the relative (100-bin) `quals` table during
  collection (`AddQual`, graphdata.go:525) while leaving the absolute
  per-position table (emitted as `qualsbin`) untouched, and (b) writes
  `"scale":0` instead of `"scale":1` into the emitted `.gd` JSON
  (graphdata.go:1571-1576). It is a graph-data-only knob — the
  `stats`/`filter` paths do no fixed-bin quality scaling, so the flag
  has no effect there (matching upstream, where `$scale` is referenced
  only by `getQualStats` and the graph-data JSON emitter). Validated
  byte-for-byte (semantic-tree-equal) against the live upstream Perl
  oracle for both scaled and `--qual_noscale` runs in
  `tools/prinseq/pkg/prinseq/qual_noscale_test.go`
  (`runUpstreamPrinseqQualNoscale`).
- The byte-level emit deviates from upstream in one
  way: **map keys are emitted in lexicographic order**, where
  upstream uses Perl-hash iteration order (Perl >= 5.18 randomises
  this every interpreter start). Documented as an intentional
  improvement in `docs/UPSTREAM_BUGS.md > prinseq` and validated
  via a JSON-normalised semantic diff against the upstream-shipped
  `example1.gd`. See `tools/PARITY_VALIDATION.md > prinseq parity
  validation` for the test layout.

Quality-trim window/step/rule + range filters (PR
`claude/festive-planck-n9o2lm-prinseq-trim-range`):

- `--trim_qual_window <int>` / `--trim_qual_step <int>` /
  `--trim_qual_type <min|max|mean|sum>` / `--trim_qual_rule <lt|gt|et>`
  — the full sliding-window quality-trimming machinery
  (`prinseq-lite.pl:3215-3287`). The window starts at the read end,
  advances by the step each time the rule passes, shrinks to the
  remaining bases near the far end, and aggregates the window score by
  type before comparing it against `--trim_qual_left/right`. Implemented
  in `trimQualityWindow` / `scanQualTrim`; defaults match upstream
  (window=1, step=1, type=min, rule=lt). The existing per-base
  `--trim_qual_left/right` flags now route through this code as the
  window=1 case.
- `--trim_to_len <int>` — hard-trim each read to at most N bases from
  the 5' end, applied after all other trimming and before the
  length/GC filters (`prinseq-lite.pl:3382-3385`). The `length > N`
  guard means equal-or-shorter reads are untouched.
- `--range_len <ranges>` / `--range_gc <ranges>` — comma-separated
  `min-max` range filters on the trimmed length and the integer GC
  percentage (`prinseq-lite.pl:3403, 3458`, `checkRange` at 2548).
  GC% is truncated toward zero (`sprintf("%d", gc*100/len)`) before the
  comparison, and the upstream AND-across-ranges semantics of
  `checkRange` are reproduced exactly (disjoint ranges therefore reject
  everything — mirrored deliberately).

These also surfaced upstream's `zero_length` filter
(`prinseq-lite.pl:3389-3392`): a read trimmed to length 0 is now dropped
rather than emitted as an empty record.

**Validation:** byte-for-byte against the live Perl oracle. The tests in
`tools/prinseq/pkg/prinseq/trim_range_parity_test.go` run
`perl reference_code/prinseq/prinseq-lite.pl` and the Go port on the same
fixtures and compare the `-out_good` output byte-for-byte (these flags use
only core Perl modules and produce deterministic output). Strong unit tests
live in `tools/prinseq/pkg/prinseq/trim_range_test.go` and CLI-wiring tests
in `tools/prinseq/cmd/prinseq/trim_range_cli_test.go`.

Sequence/header transforms + misc knobs (PR
`claude/festive-planck-n9o2lm-prinseq-transforms-misc`):

- `--seq_case <upper|lower>` (lines 912-918, 3664-3671) — force the
  emitted sequence case.
- `--dna_rna <dna|rna>` (lines 920-927, 3672-3679) — convert T<->U,
  preserving case.
- `--rm_header` (lines 133, 3651-3653) — drop the original trailing
  header comment, leaving only the identifier.
- `--no_qual_header` (lines 153, 792-793, 3686) — emit a bare `+`
  line in FASTQ output. Validated to require FASTQ output.
- `--line_width <int>` (lines 132, 934-936, 3699-3714) — wrap
  FASTA/QUAL output at N chars (0 = no wrap); feeds the shared
  `QualLineWidth`.
- `--seq_num <int>` (lines 110, 3147-3150) — keep only the first N
  records that pass all other filters.
- `--exact_only` (lines 154, 833-849, 3604-3628) — restrict
  duplicate detection to exact (and revcomp) duplicates. Our derep is
  exact-only by construction, so this is accepted and validated to
  require `--derep` 1/4/5.
- `--params <file>` / `--custom_params <string>` (lines 553-560,
  1069-1085, 2353-2369, 3484-3503) — parameters file and custom
  dinucleotide/repeat complexity rules. Rules are `<bases> <count>`
  (`%` count = percentage); evaluated against the upper-cased sequence.

Closed gaps:

- The PNG report generation flow (`prinseq-graphs.pl`) is
  **implemented** as the `prinseq graph_png` subcommand (matching the
  upstream `-i <gd> -png_all [-html_all] -o <prefix>` surface). It
  reads a `--graph_data` `.gd` file and renders the full upstream
  graph set — length (`_ld`), poly-A/T tail (`_td5`/`_td3`), N% (`_ns`),
  GC (`_gc`), DUST/entropy complexity (`_cd`/`_ce`), dinucleotide-odds
  PCA (`_pm`/`_pv`) + odds-ratio (`_or`), base-quality boxplots
  (`_qd`/`_qd2`), per-read mean-quality bar (`_qd3`), and the three
  duplication-level stacks (`_df`/`_dl`/`_dm`) — plus the optional
  HTML index, all in pure Go (stdlib `image`/`image/png` + a
  hand-rolled 5x7 bitmap font; gonum drives the PCA eigendecomposition
  per its sanctioned linalg scope). **PNG byte-identity is N/A** — the
  renderer is stdlib raster, not Perl Cairo/GD, so pixels differ by
  design; the asserted parity is the graph *set* (filenames) and the
  plotted *data series*, both unit-tested in `pnggraphs_test.go`. (The
  upstream `prinseq-graphs.pl` could not be run as a live oracle here
  because its `Cairo`/`Statistics::PCA` CPAN deps are not installable;
  the graph set was therefore derived from the Perl `generateGraphs`
  source.)

**Validation:** the transform/misc knobs are validated **live** against the
upstream Perl `prinseq-lite.pl` — `TestParityTransforms` in
`tools/prinseq/pkg/prinseq/transforms_test.go` runs the real script via a
uniquely-named `runUpstreamPrinseqXform` oracle and compares the good-output
file byte-for-byte for 17 flag/format combinations (skipping only when perl or
the submodule is unavailable). Strong unit tests cover the transform helpers
and the dedup/seq_num/line_width paths; CLI integration tests in
`tools/prinseq/cmd/prinseq/transforms_cli_test.go` cover `--params` precedence,
`--custom_params`, and the validation error paths. Earlier flags remain covered
by `missing_flags_test.go` (`--ns_max_p` boundaries, `--out_format`,
`--noniupac`, `--seq_id_mappings`, `--phred64`).

### `sickle`

**Status:** 2 subcommands (`se`, `pe`). **1:1 parity validated.**

A 15-case parity corpus at `tools/sickle/testdata/parity/` covers basic
SE, `-n` (truncate-N), `-x` (no-5'-trim), illumina (Phred+64), empty
input, all-low-quality, the `-q`/`-l` boundary, gzip input, the
short-read filter, PE with singletons, synced PE pass/fail, **strict
(`-q 30 -l 5`)**, **lax (`-q 0 -l 0`)**, and **strict PE singletons
(`-q 30 -l 10`)**. Each `case*.expected.fq` was generated by running
upstream sickle v1.33 with the documented flags and is asserted
byte-for-byte by `TestParity_Sickle_*` in
`tools/sickle/pkg/sickle/parity_test.go`. See
the `tools/sickle` README parity notes
for the per-case description and the audit's "discrepancies found and
fixed" log.

The **default sliding window** now matches upstream exactly: upstream
sickle has no `-w` flag and always uses a dynamic per-read window of
`int(0.1 * read_length)` (falling back to the full read length when that
rounds to 0). A previous bug defaulted our `-w` to a hardcoded `10`, which
trimmed ~1% of reads one window short on reads ≠ 100 bp. The default is now
`0` (dynamic); a positive `-w` is a Go-port extension that pins a fixed
window. This is validated live against the upstream binary on
varying-length reads (SE+PE) by `TestParityDynamicWindow*` in
`tools/sickle/cmd/sickle/dynamic_window_parity_test.go`, and the
length→window rule is pinned binary-free by `TestUnitResolveWindowSize`.

Outstanding items (Go-port extensions, not parity gaps):

- Auto-detect heuristic (PR #33) is a Go-port extension with no upstream
  equivalent; exercised by `encoding_test.go`.
- The fixed `-w N` window (used only when `N > 0`) is a Go-port extension;
  upstream has no `-w` flag.
- Gzip *output* level was not part of the audit (parity is asserted on
  the trimmed FASTQ records, not on the gzip container bytes).

### `skewer`

**Status:** 2 subcommands (`se`, `pe`). **1:1 parity validated** — all
14 parity cases byte-match upstream skewer 0.2.2.

A 14-case parity corpus at `tools/skewer/testdata/parity/` covers SE 3'
trim, SE 5' (`-m head`), `-m any`, min-overlap, qual+adapter, length
filter, empty input, adapter-at-end, gzip input, off-by-one boundary,
**no-adapter pass-through**, **long reads (>40 bp) with embedded
adapter**, **PE matrix-mode pass-through**, and **error-tolerant matcher
rejection** (1-mismatch at Q40). All fourteen byte-match upstream
skewer 0.2.2 (`-r 0.1` defaults).

The two previously-skipped cases were closed by porting the relevant
algorithms verbatim from `reference_code/skewer/src/matrix.cpp`:

- **case04 — PE matrix mode (`-m pe`).** Ported
  `cMatrix::findAdapterWithPE` and `CalcRevCompScore`
  (matrix.cpp:487-522, 726-851) into
  `tools/skewer/pkg/skewer/skewer.go` as
  `detectPairedTrim` / `calcRevCompScore`. Matrix mode is exposed via
  the `PEMatrixMode` option (set to `true` by the PE CLI to mirror
  upstream's default `-m pe`).
- **case05 — error-tolerant matcher.** Ported the quality-weighted
  penalty model from `cAdapter::align` and the `cMatrix::penalty[]`
  ramp (matrix.cpp:138-141, 297-435, 547-556) into
  `findAdapterWithQual` / `mismatchPenalty`. Each mismatch costs a
  quality-derived penalty (Q40+ → MAX_PENALTY=4.477), and the match is
  rejected when the cumulative penalty exceeds
  `dPenaltyPerErr * compareLen + 0.001`.

See the `tools/skewer` README parity notes
for the per-case description.

### `fastp`

**Status:** Single `fastp` command with sliding-window cut, auto adapter
detection, HTML+JSON reports, duplication evaluation, UMI processing,
**overlap-based PE base correction (`--correction`), overrepresentation
analysis (`-p`/`-P`), output splitting (`-s`/`-S`/`-d`), the overlap-driven
merge writer (`-m`/`--merge`/`--merged_out`/`--include_unmerged`),
`--adapter_fasta`, `--poly_x_min_len`, and the full set of `--disable_*`
toggles (`--disable_adapter_trimming` `-A`, `--disable_quality_filtering`
`-Q`, `--disable_length_filtering` `-L`, `--disable_trim_poly_g` `-G`)**.

Done in the `claude/festive-planck-n9o2lm-fastp-tail` PR:

- **Base-correction in PE overlap** (`-c`/`--correction`): verbatim port of
  upstream `OverlapAnalysis::analyze` (overlapanalysis.cpp) +
  `BaseCorrector::correctByOverlapAnalysis` (basecorrector.cpp) in
  `tools/fastp/pkg/fastp/overlap.go`. Byte-for-byte parity with upstream on
  corrected reads, plus matching `corrected_reads`/`corrected_bases` JSON
  fields. Knobs `--overlap_len_require` / `--overlap_diff_limit` /
  `--overlap_diff_percent_limit` wired through.
- **Overrepresented sequence analysis** (`-p`/`--overrepresentation_analysis`,
  `-P`/`--overrepresentation_sampling`): verbatim two-phase port of
  `Evaluator::computeOverRepSeq` + `Stats::statRead`/`overRepPassed` in
  `tools/fastp/pkg/fastp/overrepresentation.go`. Emits
  `overrepresented_sequences` under `read{1,2}_before_filtering` in JSON,
  matching upstream sequence set + counts exactly on the test fixture.
- **Splitting output** (`-s`/`--split`, `-S`/`--split_by_lines`,
  `-d`/`--split_prefix_digits`): `tools/fastp/pkg/fastp/split.go`. File
  naming and per-file boundaries match upstream byte-for-byte in BOTH
  single- and multi-thread mode (`-w 1` and `-w 4`), including the
  PACK_SIZE=256 rollover quantization. The port reproduces upstream's exact
  pack/thread assignment: the reader hands pack `i` to thread `i % thread`
  and each thread owns a strided, disjoint set of split files (start index
  = threadId, stride = thread), rolling its current file at a pack boundary
  once its accumulated count reaches the per-file size (input-read count
  capped at `N-1` for `--split N`; uncapped passed-read count for
  `--split_by_lines`). Because the pack→thread counter is deterministic and
  each thread consumes its packs in FIFO order, the per-file contents are
  fully determined by the thread count — so we match upstream for any `-w`.
  IMPORTANT: upstream is **not** thread-count-invariant (the task's earlier
  "preserves input order regardless of thread count" premise was incorrect):
  for the same input, `-w 1` and `-w 4` send different reads to a given
  numbered file (e.g. `-S 4000` over 6000 reads yields 6 files at `-w 1`, 7
  at `-w 2` with a trailing empty, 8 at `-w 4`). The contract is byte-parity
  with upstream **per thread count**, validated in
  `parity_split_mt_test.go` against the upstream binary for `-w 2/3/4`
  (SE byNumber, SE byLines, SE byLines+filter, PE) plus the `-w 1` cases in
  `parity_tail_test.go`. Binary-free `TestUnitSplit*` pin the assignment,
  the zero-padded numbering, and the (intended) thread-dependence.

Done in the `claude/festive-planck-n9o2lm-fastp-merge` PR (merge writer +
adapter-fasta + poly-X knob + disable-adapter flag):

- **Merge mode** (`-m`/`--merge`, `--merged_out`, `--include_unmerged`):
  verbatim port of upstream `OverlapAnalysis::merge`
  (overlapanalysis.cpp:152-183) plus the peprocessor merge dispatch
  (peprocessor.cpp:488-523) in `tools/fastp/pkg/fastp/merge.go`
  (`mergeOverlappedPair`) and `fastp.go` (`processMergePair`). The
  per-read pipeline is split into `trimRecord` + `filterRecord` so the
  merge overlap analysis interposes between trimming and filtering, as
  upstream does. Merge auto-enables base correction
  (options.cpp:115-117), reproduced here. Merged read names carry the
  upstream `merged_<len1>_<len2>` suffix. The legacy `--merge-overlap`
  heuristic is retained as a separate flag for back-compat.
- **`--adapter_fasta FILE`**: verbatim port of
  `AdapterTrimmer::trimByMultiSequences` / `trimBySequence`
  (adaptertrimmer.cpp:47-170) + `Matcher::matchWithOneInsertion`
  (matcher.cpp:10-54) in `merge.go`, with the FASTA loader mirroring
  `Options::loadFastaAdapters` (options.cpp:52-79) — adapters < 6 bp
  skipped, ordered by sorted contig name (upstream `std::map` order).
- **`--poly_x_min_len`**: the naive consecutive-base poly-X counter was
  replaced by a verbatim port of `PolyX::trimPolyX` (polyx.cpp:49-116)
  with its own length knob (`PolyXMinLen`), independent of poly-G's.
- **`--disable_adapter_trimming` (`-A`)**: explicit flag wired to gate the
  entire adapter block (matching upstream `adapter.enabled =
  !disable_adapter_trimming`).

Done in the cut_tail / disable-flags parity follow-up:

- **`--disable_quality_filtering` (`-Q`) / `--disable_length_filtering`
  (`-L`) / `--disable_trim_poly_g` (`-G`)**: the remaining upstream
  `--disable_*` toggles, wired in `cmd/fastp/main.go` and gated in
  `fastp.go::filterRecord` / `trimRecord`. Quality filtering gates the
  low-quality-percent check **and the N-base limit** (upstream keeps the N
  check inside the quality block, `filter.cpp:43-50`); length filtering gates
  both the too-short (`length_required`) and too-long (`length_limit`) checks
  (`filter.cpp:52-57`); `-G` force-disables poly-G even when `--trim_poly_g`
  was requested. Validated byte-exact in parity Case 20 (alone + combined) and
  the binary-free `TestUnitDisableQualityFiltering` /
  `TestUnitDisableLengthFiltering` / `TestUnitDisableTrimPolyG`.
- **cut_tail window-boundary at scale** — closed. The parity pipeline found
  `--cut_tail` diverging from upstream on ~1% of reads by ~1bp; the small
  per-tool fixtures never triggered it. Root cause: the Go port applied the
  sliding-window cut **after** adapter/poly-G trimming, whereas upstream runs
  `Filter::trimAndCut` **first** (`seprocessor.cpp:235`, before the poly-G and
  adapter blocks). Cutting an already-shortened read shifts the window math
  (the `s+w<l-tail` bound and the `t = t-w+1` rewind on the passing window),
  so any read where adapter/poly trimming fired drifted by ~1bp. Fixed by
  reordering `trimRecord` to upstream's `cut → poly-G → adapter → poly-X`
  order. The `slidingWindowCut` algorithm itself was already byte-exact
  (confirmed by Case 18, cut_tail with adapter disabled). Validated
  byte-exact at scale on the new 5000-read varying-length/low-Q-tail fixture
  `se_cuttail_scale.fq` in parity Case 17 (with adapter) and Case 18 (without);
  cut_front/cut_right confirmed unregressed at scale in Case 19; the boundary
  rule pinned binary-free in `TestUnitCutTailBoundary`.

Both formerly-remaining residuals are now **closed** (multi-thread split
file-boundary distribution; unemitted JSON sub-fields):

- **Multi-thread `--split` file-boundary distribution** — closed. See the
  splitting bullet above: byte-parity with upstream for both `-w 1` and `-w 4`
  (and `-w 2/3`), for `--split N` and `--split_by_lines`, SE and PE.
- **JSON report sub-fields** — closed. The per-read blocks now emit, for BOTH
  the before- and after-filtering streams, the full upstream
  `Stats::reportJson` shape: `q40_bases`; `quality_curves` with per-base
  `A`/`T`/`C`/`G` + `mean` (previously only `mean`); `content_curves` with
  `A`/`T`/`C`/`G`/`N`/`GC` in upstream order (previously a 5-key subset in the
  wrong order, before-only); the 1024-entry `kmer_count` 5-mer histogram
  (previously absent); and the real per-read `q20_bases`/`q30_bases` totals
  (previously hard-coded to 0, before-only). The top-level report also now
  emits `summary.sequencing` (the deterministic "single/paired end (N
  cycles[ + N cycles])" descriptor) and, for paired-end, the `insert_size`
  block (`peak`, `unknown`, full 512-bin `histogram`) reproduced from
  upstream's overlap-analysis insert-size binning
  (`PairEndProcessor::statInsertSize`). All of these are validated
  byte/structurally-exact against the upstream binary in
  `parity_json_fields_test.go` (the curves are deterministic but compared
  within a 1e-4 absolute tolerance because upstream emits them through C++
  6-significant-digit float formatting; integer counts and the histogram are
  compared exactly), with binary-free `TestUnit*` for the curve/kmer/insert
  builders.

  **Intentionally excluded (non-reproducible):** `summary.fastp_version` (the
  upstream version string vs our `tool.version`), the top-level `command`
  string, and our `tool.time` wall-clock field — these are not deterministic
  across runs/builds and are excluded from the parity comparison, consistent
  with the existing JSON parity policy.

Schema note: our per-read blocks additionally carry a `length_distribution`
object that upstream's JSON omits (upstream surfaces that data only in the HTML
report). It is a documented Go extension, `omitempty`, and does not affect the
upstream-field comparison. The merge-mode `merged_and_filtered` block rename is
emitted (the after-filtering block is renamed and `read2_after_filtering`
dropped when `-m`/`--merge` is set).

#### Validated-parity audit

21-case test corpus at `tools/fastp/pkg/fastp/parity_test.go` against
upstream fastp 1.0.1. **21 PASS, 0 SKIP** (Cases 17-20 added by the
cut_tail / disable-flags follow-up: cut_tail-at-scale with and without
adapter, cut_front/cut_right unregressed at scale, and the new `--disable_*`
toggles; deterministic transforms validated byte-exact, the heuristic SE
adapter-detect validated by a documented similarity bound). See
the `tools/fastp` README parity notes
for the case list.

Bugs in the Go port surfaced + fixed inline by the initial audit:

- **UMI tag format** was unconditionally `":UMI_<umi>"`. Upstream uses
  `":<umi>"` (no prefix) or `":<prefix>_<umi>"` (with prefix). Fixed.
- **Low-complexity definition** was "unique 2-mers / total 2-mers".
  Upstream uses "fraction of adjacent positions where seq[i] !=
  seq[i+1]". Fixed.
- **`low_complexity_reads` JSON counter** was missing. Added.

Bugs fixed inline by the `claude/fastp-algorithmic-fixes` follow-up PR:

- **PolyG mismatch tolerance**: upstream's `trimPolyG` tolerates 1
  mismatch per 8 bases scanned (capped at 5 total) and anchors on the
  last-G position (`reference_code/fastp/src/polyx.cpp:16-42`). The Go
  port now runs a verbatim port (`trimPolyG` in
  `tools/fastp/pkg/fastp/fastp.go`). `TestParity_Fastp_Case12_SEPolyG`
  is no longer skipped.
- **Sliding-window boundary** (`cut_front` / `cut_tail` / `cut_right`):
  `slidingWindowCut` is now a verbatim port of upstream's
  `filter.cpp:83-222`. Specifically: (a) cut_right walks the high-Q
  prefix inside the offending bad window, (b) cut_front and cut_tail
  rewind to the START of the qualifying window (`s+w-1` for front, `t-w+1`
  for tail) and then skip N's at the boundary, (c) the loop bound stays
  strictly `s + w < l` so the trailing w bases are never scanned.
  `TestParity_Fastp_Case13_SECutRight` and
  `TestParity_Fastp_Case14_SECutFrontTail` are no longer skipped. (The
  cut_tail / disable-flags follow-up later corrected the step **order** so
  the cut runs before adapter/poly trimming, fixing a ~1bp at-scale drift;
  see below.)

Bugs fixed inline by the `claude/fastp-adapter-autodetect` follow-up PR:

- **SE adapter auto-detect**: upstream's
  `Evaluator::evalAdapterAndReadNum` (`evaluator.cpp:295-526`),
  `Evaluator::checkKnownAdapters` (`evaluator.cpp:207-293`),
  `Evaluator::getAdapterWithSeed` (`evaluator.cpp:472-526`), and
  `NucleotideTree` (`nucleotidetree.cpp`) are now ported verbatim into
  `tools/fastp/pkg/fastp/adapter_autodetect.go` and
  `tools/fastp/pkg/fastp/known_adapters.go`. The Go port now reproduces
  upstream's behavior byte-for-byte, including the 10000-record gate at
  `evaluator.cpp:344` (below that threshold the evaluator returns ""
  and no adapter trimming is applied — same as upstream's "No adapter
  detected for read1" path). `TestParity_Fastp_Case15_SEAutoDetect` is
  no longer skipped.

  Follow-up (`claude/fastp-adapter-similarity`): the SE auto-detect path
  is genuinely heuristic / sampling-dependent, so it is now validated by
  a documented SIMILARITY BOUND against the upstream binary on a fixture
  that actually fires detection — `se_detect.fq` (12000 reads, above the
  gate). `TestParity_Fastp_Case16_SEAutoDetectFires` asserts: detected
  adapter equals upstream's (prefix within 3bp); per-read trimmed-length
  agreement >= 99% with no read off by > 2bp and >= 99.9% base identity;
  and adapter-trimmed reads/bases + passed-filter reads within 1%
  relative tolerance. Building that case surfaced a real bug: the single
  configured/detected 3' adapter was trimmed with a plain
  `strings.Index` substring search, which cannot match a partial adapter
  prefix at the 3' read end (the common read-through case) and tolerates
  no mismatch. It now routes through the verbatim
  `AdapterTrimmer::trimBySequence` (`adaptertrimmer.cpp:71-170`,
  `matchReq=4`, as `seprocessor.cpp:245` does), which matches partial 3'
  prefixes (`cmplen = min(rlen-pos, alen)`), tolerates 1 mismatch per 8
  overlapping bases, and applies the A-tailing negative start. With that
  fix the observed Case 16 agreement is in fact byte-identical
  (lenAgreement 1.0, maxLenDelta 0, baseIdentity 1.0); the heuristic
  CONTRACT remains the similarity bound. Binary-free `TestUnitTrimPolyG`,
  `TestUnitTrimPolyX`, `TestUnitSlidingWindowCut`,
  `TestUnitDetectAdapterSE` and `TestUnitFastqSimilarityHelper` pin all
  four helpers with the reference_code submodule unpopulated.

Features added by the `claude/festive-planck-n9o2lm-fastp-tail` PR
(`--correction` + overlap analysis, overrepresentation analysis, output
splitting), with live-upstream parity tests in
`tools/fastp/pkg/fastp/parity_tail_test.go` (uniquely-named
`upstreamFastp` builder, `sync.Once`, `t.Fatalf` when the binary is
available):

- `TestParity_Fastp_Correction`: corrected R1/R2 bytes and
  `corrected_reads`/`corrected_bases` JSON match upstream exactly.
- `TestParity_Fastp_Overrepresentation`: `overrepresented_sequences` set
  and counts match upstream exactly.
- `TestParity_Fastp_SplitByLines` / `TestParity_Fastp_SplitByNumber`:
  all numbered split files byte-for-byte identical to upstream `-w 1`.

Features added by the `claude/festive-planck-n9o2lm-fastp-merge` PR
(merge writer, `--adapter_fasta`, `--poly_x_min_len`,
`--disable_adapter_trimming`), with live-upstream parity tests in
`tools/fastp/pkg/fastp/parity_merge_test.go` (uniquely-named
`ensureUpstreamFastpMerge` builder, `sync.Once`, `t.Fatalf` — never
`t.Skip`):

- `TestParity_Fastp_Merge_Basic`: merged FASTQ bytes match upstream
  `--merge --merged_out` exactly on the overlapping `corr_*` fixtures.
- `TestParity_Fastp_Merge_IncludeUnmerged`: `--include_unmerged` routes
  surviving mates into the merge stream byte-for-byte on non-overlapping
  `pe_*` fixtures.
- `TestParity_Fastp_AdapterFasta`: `--adapter_fasta` output matches
  upstream byte-for-byte.
- `TestParity_Fastp_DisableAdapterTrimming`: `-A` output matches upstream.
- `TestParity_Fastp_PolyXMinLen`: `--poly_x_min_len` at 5/10/15 matches
  upstream `PolyX::trimPolyX` byte-for-byte.

**Validation:** **16-case core parity suite (16 passing, 0 `t.Skip`) plus
4 tail-feature live-parity tests plus 5 merge/adapter/poly-X live-parity
tests. 1:1 parity achieved** for `--correction`, overrepresentation
analysis, single-thread splitting, the overlap merge writer,
`--adapter_fasta`, `--poly_x_min_len`, and `--disable_adapter_trimming`.
The only documented residuals are the multi-thread split distribution and
the `merged_and_filtered` JSON block (merged FASTQ output bytes are
byte-identical).

### bedtools (all subcommands ported)

**Status:** ~95%. All upstream bedtools subcommands are ported across the
37 `bed*` tool dirs (no missing subcommands). BAM input (`bedintersect`,
`bedmulticov`, incl. CRAM on `bedmulticov`), VCF/GFF input (`bedintersect`,
`bedmultiinter`), `bedclosest` direction flags
(`-D`/`-id`/`-iu`/`-fu`/`-fd`/`-t`), `bedintersect -c`, and `bedsample`
RNG byte-parity have all landed. Remaining work is scattered per-subcommand
option-tail polish (and CRAM input on the BED-only tools). 141 passing parity
tests
against the upstream test suite (across PR #55 + Phase-3 wave 1 + wave
2 simple + wave 2 algo) + 17 new cases from wave 3 (PR #87) + 6 cases
from the reldist/fisher full-parity wave (PR #90) + 6 cases from the
nuc/annotate wave (PR #91) + 7 cases from the multicov/multiinter
wave (PR #92) + 6 cases from the igv/links wave (PR #93) + 10 cases
from the BEDPE pair-ops wave (PR #94, five each for `bedpairtobed`
and `bedpairtopair`, derived from the upstream source as neither
subcommand ships a `test/<name>/` subdir) + **9 newly passing**
parity cases from the column-ops + discrepancies wave (this PR —
`jaccard.t02/t03/t05/t06/t10/t11`, `merge.t15`, `map.t11`,
`map.t13`).
Phase-3 wave 1 (PR #78) added `bedgroupby`/`bed12tobed6`/`bedmakewindows`;
wave 2 simple (PR #80) added `bedexpand`/`bedgetfasta`/`bedsample`/`bedspacing`;
wave 2 algo added `bedcoverage`/`bedmap`/`bedshuffle`; wave 3 tail
(PR #87) added `bedcluster`/`bedsplit`/`bedsummary`/`bedtag`/`bedwindow`;
the reldist/fisher wave (PR #90) added `bedreldist`/`bedfisher`; the
nuc/annotate wave (PR #91) added `bednuc`/`bedannotate`; the
multicov/multiinter wave (PR #92) closed the last two of the six
originally-planned algorithmic subcommands; the igv/links wave
(PR #93) landed the two pure-format converters; **the BEDPE pair-ops
wave (this PR) closes the last two missing subcommands — bedtools now
has no missing subcommands**.

Missing subcommands: none. All algorithmic, BEDPE, and converter
subcommands are ported. Remaining bedtools work is option-tail polish.

Parity-pipeline bug fixes (this wave): five byte-parity defects the parity
pipeline surfaced were fixed and live-validated against the upstream
binary:

- **`bednuc`** — a duplicate `--seq` flag registration made the CLI panic
  (`flag redefined: seq`) on every invocation; the redundant
  registrations were removed. `bednuc` now runs and matches upstream
  `bedtools nuc` (default, `-seq`, `-s`).
- **`bed12tobed6`** — the score column was hardcoded to `0`; it now
  carries the parent record's score onto each emitted BED6 block, and
  `-n` numbering reverses for every non-`+` strand (not just `-`),
  matching `bed12ToBed6.cpp`. See UPSTREAM_BUGS.md
  (`bedtools-bed12tobed6-n-strand`).
- **`bedexpand`** — a trailing comma in an expanded column no longer
  emits a spurious empty final row; comma tokenization now matches C++
  `getline` (leading/interior empties preserved, single terminating
  empty dropped).
- **`bedmakewindows`** — the `-i` default (`ID_NONE`/BED3) was wrong
  (it produced a window-number column) and `-i src` over a genome file
  emitted an empty name; the default now yields BED3 and genome
  intervals annotate with the chromosome name. See UPSTREAM_BUGS.md
  (`bedtools-makewindows-i-none`).
- **`bedsplit`** — the default `-a size` heuristic re-sorted each
  output bin back into input order; it now preserves the
  size-descending insertion order, so the emitted shard files are
  byte-identical to upstream `bedtools split -a size`.

Each fix ships a live-binary parity test (`live_parity_test.go`) that
`t.Fatalf`s when the upstream binary is absent, plus binary-free
`TestUnit*` for changed helpers.

`bedintersect` behavioral-flags wave: the previously
recognized-but-unimplemented join/overlap output flags now run to
byte-for-byte parity with the live upstream binary:

- **`-loj`** (left outer join): every A record is echoed, with the
  overlapping B appended, or a null-B placeholder when A has no overlap.
  The null shape is derived from the B file's detected record type
  (`. -1 -1` for BED3, plus the right trailing `.`/`-1` columns for
  BED4/5/6/12/bedGraph/BED+ — mirroring `RecordOutputMgr::null` and the
  per-type `printNull` methods).
- **`-wo`** (write overlap): A, B and the overlapping base count per
  overlap. A-with-no-hits is omitted.
- **`-wao`** (write all + overlap): like `-wo` but A-with-no-hits is also
  emitted, with a null B and an overlap of `0`.
- **`-wa -wb`** combined: A and B columns side-by-side per overlap.
- **`-split`**: BED12 records are split into their blocks; a B is a hit
  only if a block of A overlaps a block of B, and the `-wo`/`-wao`
  overlap count is the summed per-block overlap. The `-f`/`-F`/`-r`
  fraction tests under `-split` use the combined non-redundant overlap
  over the A and B block sums (once per A, dropping all hits on failure),
  matching `BlockMgr::findBlockedOverlaps`; the non-split path keeps the
  per-record fraction test from `Record::sameChromIntersects`.

These modes echo the original A and B input columns verbatim and in the
original B-file order (not sorted), via a raw line-preserving parser in
`tools/bedintersect/pkg/bedintersect/join.go`. Zero-length records
(`start==end`) are expanded to `[p-1,p+1]` for detection (non-split and
split paths) with the upstream overlap-count corrections. The fix also
corrected the `-s` same-strand filter across all intersect paths: `.`,
`*` and a missing strand column are UNKNOWN and can never satisfy `-s`
(previously they were treated as wildcards, over-reporting hits).

Parity is enforced by `cmd/bedintersect/behavioral_parity_test.go`,
which builds the real upstream `bedtools` (uniquely-named
`upstreamBedtoolsBehavioral` builder, `sync.Once`) and asserts
byte-for-byte equality (`t.Fatalf`, never `t.Skip`) across the BED3–BED12
null shapes, `-split` block math, zero-length intervals, B-order
preservation, and the UNKNOWN-strand cases.

bedtools-gaps closure wave (this PR): the items previously listed as
`bedintersect`/`bedclosest` remainder are now implemented and pinned
byte-for-byte against the live upstream `bedtools` binary
(`upstreamBedtoolsGaps`, `sync.Once`, `t.Fatalf` — never `t.Skip`):

- **BAM / VCF / GFF inputs on `-a`/`-b`.** `tools/bedintersect/pkg/bedintersect/input.go`
  autodetects the format (BAM by BGZF/`BAM\x01` magic; otherwise BED/VCF/GFF
  by upstream `BedFile::parseLine` precedence) and converts each record to a
  common 0-based half-open `inRecord` that echoes its original columns
  verbatim. BAM alignments render as BED12 with CIGAR-N block splitting
  (`-bed`), VCF spans are `POS-1 .. POS-1+len(REF)` (the line echoes
  unclipped), GFF is `start-1 .. end`. Clipped-BAM output keeps the original
  thickStart/thickEnd/block columns, matching upstream `BamRecord::print`
  exactly (the blocks describe the original span, not the clip — verified
  against upstream `-split`/`-wo`). Parity:
  `cmd/bedintersect/gaps_parity_test.go` (`TestGapsParity_BAMInput`,
  `TestGapsParity_VCFGFFInput`, including the `-s` strandless case).
- **`-c` count column.** The raw column-preserving path echoes every original
  column before the trailing count (the old typed `bed.Writer` dropped a
  score of `0` and the strand column). Pinned by
  `TestGapsParity_CountColumnDrop`.
- **`bedclosest` directional flags (`-iu`/`-id`/`-fu`/`-fd`).** Implemented
  with the `-D`-orientation sweep semantics from upstream `CloseSweep`
  (`classifyStream`). The "requires `-D`" guard matches upstream exactly:
  any `-D` value (including `-D ref`) satisfies it — upstream only errors when
  `-D` is entirely absent (verified against the live binary; `_haveStrandedDistMode`
  is set for `-D ref` too). Pinned by `TestGapsParity_ClosestDirectional`.
- **`bedclosest -D a/-D b` signed-distance sign fix.** The non-directional
  `signedDistance` previously flipped the `-D b` sign on a reverse-strand B;
  upstream (`_bDist && _dbForward`) flips on a FORWARD-strand B. `signedDistance`
  now reuses `classifyStream` so the directional and non-directional paths
  agree, matching upstream across every geometry/strand/mode
  (`TestGapsParity_ClosestSignedDistance`, 36 combinations).
- **`-sortedtree`/`--tree` B index.** `opts.UseTree` now actually selects an
  augmented interval tree over `inRecord` (`treeFinder`/`inIntervalTree`)
  instead of being a silent no-op in an unreachable typed branch; the tree
  path re-sorts candidates into B-file order so output stays byte-identical
  to the linear scan and to upstream (`TestGapsParity_UseTree`). The dead
  typed `findOverlaps`/`Overlap`/`splitBlockOverlap`/`typedBlocks` helpers and
  the `bed.Record` `IntervalTree` alias they used were removed.

Residual `bedintersect` notes (still as before):

- **Ragged B files.** Upstream hard-errors when a B file's records have
  inconsistent field counts; the null-shape classification here keys off
  the first B record. Genuine BED is uniform, so this only differs on
  malformed input (ours is lenient where upstream aborts). The new
  BED/VCF/GFF reader does validate the per-line field count against the
  locked format and errors like upstream's type checker.

CRAM-input + option-tail wave (this PR): `bedmulticov` now accepts CRAM
inputs anywhere it accepts BAM; `bedmultiinter` now autodetects VCF/GFF
inputs and takes `-names` as a space-separated variadic list; `bedsample`
now reproduces upstream's seeded sampler byte-for-byte via an in-tree
`std::mt19937_64` port. All three closures are backed by live-upstream
parity tests (build the real `bedtools` binary, compare byte-for-byte,
`t.Fatalf` never `t.Skip`). See the per-subcommand notes below.

Vendored-fixture + CRAM/`-fullHeader` wave (this PR): a batch of parity
cases that were `t.Skip`'d only for missing fixtures or unwired input
paths are now real byte-for-byte assertions:

- `bedgenomecov` — the BAM/SAM/CRAM `-ibam` parity cases (upstream
  `genomecov.t1..t18`) are closed. The upstream SAM fixtures and the
  htsutil-built BAMs (`y.bam`/`empty.bam`/`merged.bam`), plus the empty
  CRAM and its `test_ref.fa` reference, are vendored under
  `tools/bedgenomecov/testdata/parity/aln/`, and each case asserts against
  a golden generated from the upstream `bedtools genomecov` binary.
  **CRAM input is now wired**: `RunBAM` routes through
  `pkg/htsgo/alnio.NewReaderWithReference`, which auto-detects SAM/BAM/CRAM;
  a new `-T/--reference` CLI flag (and `Options.CRAMReference`) supplies the
  CRAM decode FASTA, and `REF_CACHE`/`REF_PATH` are honoured. Closing these
  cases surfaced and fixed **two real parity bugs**: (1) per-base `-dz`
  output was 1-based but upstream emits **0-based** positions
  (`offset = _eachBaseZeroBased ? 0 : 1`); (2) BAM coverage did not break
  the alignment on a CIGAR `D` (deletion) op, whereas upstream's
  `getBamBlocks(..., breakOnDeletionOps=!_ignoreD, …)` always breaks on D
  (and on N only under `-split`). Both are now implemented and asserted.
- `bedreldist` — `reldist.t01..t03` (the large refseq/aluY/gerp fixtures)
  are closed. The upstream `.bed.gz` inputs are vendored under
  `tools/bedreldist/testdata/parity/` (kept gzip-compressed; decompressed
  at test time) and each case asserts byte-for-byte against an upstream
  `bedtools reldist` golden.
- `bedsplit` — `split.01..03` now assert against the upstream-script
  goldens using the vendored `randData.bed`; the size-mode (LPT) case
  matches upstream exactly and is a hard assertion (no longer skipped on a
  tie-breaking caveat).
- `bednuc` — `-fullHeader` is now a real assertion, not a skip. Empirically,
  the htslib shipped with bedtools builds the `.fai` on the first
  whitespace token even with `-fullHeader`, so a full multi-token header in
  the BED chrom column resolves to nothing and is skipped with a
  "size (0 bp)" warning, while a first-token chrom resolves exactly as in
  the default mode. This was cross-checked against both `bedtools nuc
  -fullHeader` and `bedtools getfasta -fullHeader` (the latter's
  `getfasta.t06` passes only because it uses a first-token BED chrom). The
  port reproduces this observable behaviour verbatim (no aliasing), and the
  not-found warning text now matches upstream's "Feature (...) beyond the
  length of ... size (0 bp).  Skipping." The dead best-effort
  full-header alias map was removed.
- `bedsample` — the two remaining `t.Skip`s (`sample.t01` "No input file
  given", `sample.new.t02` "Unrecognized parameter") are **intentionally**
  CLI-only and retained as documented cases: input defaulting/validation and
  unknown-flag rejection live in `cmd/bedsample/main.go` (the library always
  takes an `io.Reader`), so there is no library behaviour to assert. The
  skip messages now state this explicitly. The genome-fixture cases were
  already covered: `mainFile.bed` is vendored and the seeded sampler has
  live byte-for-byte upstream parity.

Resolved in the column-ops + discrepancies wave:

- `bedjaccard` now pre-merges A and B before computing intersection /
  union, matching upstream's `setUseMergedIntervals(true)`. Previously
  skipped parity cases `jaccard.t02 / t03 / t05 / t06 / t10 / t11` now
  pass byte-for-byte.
- `bedmerge -s` now matches upstream's per-strand merge semantics: `.`
  (unknown) strand records are dropped under `-s` and `+` / `-` groups
  merge independently before the (chrom, start) re-merge. Previously
  skipped `merge.t15` now passes.
- The full column-op vocabulary from upstream (`stdev`, `sstdev`,
  `absmin`, `absmax`, `cat`, `cat_uniq`) is wired into the shared
  `bedmerge.ApplyOp`, unblocking `bedgroupby` (`TestGroup_StdevSstdev`,
  `TestGroup_AdditionalOps`) and `bedmap` (`map.t11`, `map.t13`).

Option-tail gaps on the wave-2 additions:

- `bedgetfasta` — `-fullHeader` is now implemented: contigs are indexed
  by the full FASTA header line (whitespace included) via
  `pkg/htsgo/fasta.BuildIndexFullHeader` /
  `OpenRandomAccessFullHeader`, and `bedgetfasta -fullHeader` flows the
  flag through to the index build. Upstream `getfasta.t06` (the
  `-fullHeader` two-line case) and `t07` (the no-`-fullHeader` warning
  case) now both pass byte-for-byte. BGZF FASTA input is also wired
  through (this PR): `pkg/htsgo/fasta` now sniffs the BGZF magic
  in `OpenRandomAccess` / `OpenRandomAccessFullHeader` and routes to a
  new `OpenRandomAccessBGZF` that fully decompresses the payload
  in-memory and reuses the existing FAI index path. The `.gzi` sidecar
  (when present) is parsed for early validation via a stdlib-only
  little-endian reader in `pkg/htsgo/fasta/bgzf.go`; a samtools-
  compatible `.fa.gz.fai` is honoured when present, otherwise the
  index is rebuilt from the decompressed payload. Upstream
  `getfasta.t18` (BGZF FASTA + `-split` BED12) now passes byte-for-byte
  using the upstream `t.fa.gz` fixture. **Partial-decompression seek via
  `.gzi` is now implemented** (htsgo-gzi-bcf PR): when both the
  samtools-style `.fa.gz.fai` (uncompressed-stream offsets) and the
  `.fa.gz.gzi` block index are present, `OpenRandomAccessBGZF` serves
  `Fetch` through a `bgzf.SeekReader`-backed `io.ReaderAt` that inflates
  only the blocks overlapping each request — no whole-file decompress.
  The in-memory decompress remains as the fallback when no `.gzi` is
  present. The seek path is validated against upstream `bgzip -i`-produced
  `.gz`+`.gzi` sidecars in `pkg/htsgo/fasta/bgzf_test.go`.
- `bedsort` — `-header` is now implemented (this PR): leading
  `#`-prefixed comment, `track`, and `browser` directive lines are
  buffered and emitted verbatim ahead of the sorted body. Upstream
  `sort.t09` now passes byte-for-byte.
- `bedsort` — **tie-break fixed to match upstream input order on equal
  `(chrom, start)`** (this wave). Upstream `sortBed`
  (`loadBedFileIntoMapNoBin` → `sortByStart`) sorts each chromosome by
  `chromStart` alone, so equal-`(chrom, start)` records keep input order;
  it never uses `chromEnd` as a tie-break. The size-descending and
  `-chrThenScore{A,D}` comparators carry no secondary key and so likewise
  preserve that input order on key ties. The previous port broke ties on
  `chromEnd` ascending, diverging from upstream for the default, `-sizeD`,
  `-chrThenSizeD`, `-chrThenScoreA`, and `-chrThenScoreD` modes whenever the
  input was not already end-ordered. `Sort` now mirrors upstream's two-stage
  arrangement (stable `(chrom, start)`-only pass, then a stable mode-key pass)
  and is asserted byte-for-byte against the live `bedtools sort` binary across
  all seven modes plus `-faidx`, on an input rich in equal-key records.
- `bedwindow` — **three parity bugs fixed** (this wave):
  (1) the window is now added to **A**, not B (upstream `AddWindow` operates on
  the A feature and queries the B database), correcting the asymmetric `-l`/`-r`
  and strand `-sw` direction; (2) the **default window is 1000 bp** (the port
  defaulted to 0); (3) per-A **B-hit order now follows upstream's UCSC bin
  traversal** — B is binned by its original coordinates and hits are emitted
  finest-level-first, bin-number ascending, then B-file order (the same
  `binorder` logic bedintersect uses), instead of B-start order; and (4) records
  are kept as raw text so **BED12 (and wide) B records round-trip verbatim**
  rather than being truncated to 6 columns. `-sm`/`-Sm`/`-u`/`-c`/`-v` are
  matched; the `-c`/`-v` A-only paths are unchanged. New parity tests assert
  byte-for-byte equality against the live `bedtools window` binary for the bin
  hit-order, default-window, BED12-B, and strand cases.
- `bedmerge` — **order-sensitive column-op tie-break fixed** (this wave). The
  internal pre-sort now keys on `(chrom, start)` with input order preserved on
  ties (chromEnd is no longer a tie-break), matching the `bedtools sort`-ed
  stream upstream `merge` consumes — so `-o collapse|distinct` emit equal-start
  groups' values in input order. The merge itself was reimplemented as a faithful
  port of upstream's `FileRecordMergeMgr` state machine, including the per-strand
  `StrandQueue` priority queue: under `-s`, deferred opposite-strand records are
  pulled back out in `(chrom, start, end)` order, reproducing upstream's `-s`
  collapse/distinct ordering byte-for-byte (the previous +/- re-merge approach
  diverged on both group order and within-group value order). Order-independent
  ops (sum/mean/min/max/count) were already correct and are unchanged. New
  live-binary parity tests cover the equal-key collapse/distinct ordering across
  default, `-s`, `-S`, and `-d` paths.
- `bedmap` — **order-sensitive op tie-break fixed** (this wave). For each A
  interval the overlapping B values are now emitted in upstream's stream order —
  `(chrom, start)` with B-file order preserved on ties (chromEnd is no longer a
  tie-break) — so `-o collapse|distinct` match `bedtools map` byte-for-byte. Both
  the B-load sort and the per-A match re-sort were corrected; each B record
  carries its load-order index so the tree-query candidates can be restored to
  input order on equal starts. Order-independent ops (sum/mean/min/max/count)
  were already correct and are unchanged. New live-binary parity tests cover the
  equal-key collapse/distinct ordering (incl. `-s`/`-S` strand subsets) with
  sum/count guards against regressing the matching paths.
- `bedsample` — **byte-for-byte parity with upstream's seeded sampler is
  now achieved** (this wave). The reservoir replacement uses an in-tree,
  stdlib-only Go port of `std::mt19937_64` (`mt19937.go`) — the exact
  64-bit Mersenne Twister upstream's default (non-`USE_RAND`) build uses —
  together with upstream's exact rejection-sampling bound
  (`rand_range`: `max = mt.max() - (mt.max() % limit); do n = mt(); while
  (n >= max); return n % limit`). The fill phase draws no RNG and BED
  output is emitted in reservoir-slot order (upstream only re-sorts when
  the output type is BAM), mirroring `sampleFile.cpp`. A live-upstream
  test (`upstream_parity_test.go`) asserts byte-for-byte equality across
  several `(N, seed)` pairs plus the `-header` and `N == total` (fill-only)
  cases. Caveat: parity holds against a **default** upstream build; a
  `USE_RAND=1` build (glibc `rand()`) is platform-dependent and out of
  scope.
- `bedmulticov` — <a id="bedmulticov-bam"></a>BAM input is wired through
  via `pkg/htsgo/sam.NewBAMReader`; primary alignments contribute
  one interval each over their reference span, and `-q` MAPQ filter +
  `-D` per-A-interval depth cap are honoured. Upstream `multicov.t1`
  through `t4` and `t10` pass byte-for-byte.
  **`-split`** block-aware coverage on BAM CIGAR `N` ops is now
  implemented (this PR): the BAM index pass walks each alignment's
  CIGAR and emits one block per contiguous reference-consuming op-run
  (M/=/X, with D extending the current block — matching upstream's
  `breakOnDeletionOps=false`), skipping any `N`-op gap. Each alignment
  is counted at most once per A interval. When combined with
  `-f`, the threshold is applied to `total_block_overlap /
  sum_of_BAM_block_lengths` using strict `>` — a quirk of bedtools 2.x
  preserved here for byte-for-byte parity (mirrors
  `multiBamCov.cpp::FindBlockedOverlaps`). Upstream `multicov.t5`
  through `t9` now pass.
  **CRAM** input is now supported (this wave): `.cram` inputs are routed
  through `pkg/htsgo/alnio.NewReaderWithReference` (which dispatches CRAM
  to `pkg/htsgo/cram` and BAM to `pkg/htsgo/sam`), so every `-bams` /
  `-files` path that accepts BAM also accepts CRAM. The same `-q` MAPQ
  filter, `-D` depth cap, `-s`/`-S` strand filter, and `-split` block math
  apply unchanged (multicov reads only POS/CIGAR/FLAG/MAPQ, none of which
  need base reconstruction, so CRAM decodes correctly even without a
  reference; `-T`/`--reference` and `REF_CACHE` are honoured when present).
  A live-upstream parity test (`upstream_parity_test.go`) builds upstream
  `bedtools`, generates a `.crai` for an htslib-produced CRAM fixture via
  the in-tree `cram.CreateCRAI` (proving index interop), and asserts
  byte-for-byte parity for the default, `-q`, and `-s` cases.
- `bedmultiinter` — **VCF/GFF input is now implemented** (this wave). The
  input format is autodetected per file from the first non-header data
  line, mirroring upstream `BedFile::parseLine`'s precedence (BED when
  cols 2,3 are integers; VCF when col 2 is an integer and there are ≥8
  cols; GFF when there are exactly 8 or 9 cols with integer cols 4,5),
  with a `##fileformat=VCF` header forcing VCF. Coordinates are converted
  to 0-based half-open spans exactly as upstream does (VCF: `POS-1 ..
  POS-1+len(REF)`; GFF: `start-1 .. end`). A live-upstream parity test
  covers a mixed VCF+GFF+BED 3-way intersection and the upstream issue311
  single-record VCF/GFF fixtures, byte-for-byte. Input is still assumed
  sorted; out-of-order records within a single file are tolerated only
  because each file is re-sorted and merged before the sweep.
  CLI: `-names` is now space-separated variadic (`-names A B C`), matching
  upstream's `multiIntersectBedMain.cpp` argument loop (previously the
  port took a single comma-separated token, which broke drop-in
  compatibility).

Column-op closure: the shared `bedmerge.ApplyOp` (used by `bedmerge`,
`bedgroupby`, `bedmap`, and `bedcoverage`) now supports the full
upstream KeyListOps vocabulary: `stdev`, `sstdev`, `absmin`, `absmax`,
`cat`, `cat_uniq` (in addition to the previously-implemented sum,
min/max, mean, median, mode/antimode, count, count_distinct, distinct,
collapse, first, last). Done; no remaining gaps.

### `vcftools`

**Status:** **146 of 146 unique upstream long flags (100%)** after
long-tail wave 23. A complete `in_str ==` enumeration of
`parameters.cpp` finds 146 distinct upstream long flags; the port
registers and implements all of them. The four BCF-binary I/O flags
that earlier prose listed as blocked — `--bcf`, `--diff-bcf`,
`--recode-bcf`, `--contigs` — are all closed (waves 21–23): the
in-tree `pkg/htsgo/bcf` reader/writer (built on the in-tree BGZF
codec) supplies the binary I/O, so no external infrastructure is
needed. The PCA trio (`--pca`/`--pca-no-norm`/`--pca-snp-loadings`)
is fully implemented (wave 19) on top of gonum's symmetric
eigendecomposition — the former "LAPACK blocker" no longer applies.
The remaining work is per-output column-set polish (see the "Other"
list below), not flag-count gaps. Earlier header prose claiming
"142 of 146, blocked on 4 BCF flags" predated waves 19–23 and was
stale; this is the corrected count.

Closed in wave 1:

- **Inter-chromosomal LD**: `--interchrom-geno-r2`, `--interchrom-hap-r2` ✅
- **Chi-square LD**: `--geno-chisq` ✅
- **Relatedness**: `--relatedness` (Yang 2010), `--relatedness2` (KING-robust) ✅
- **Runs of homozygosity**: `--LROH` (+ `--LROH-min-variants`) ✅
- **Phased blocks**: `--phased-blocks` ✅
- **FILTER tag include/exclude**: `--remove-filtered`, `--keep-filtered` ✅
- **INFO selection in recode**: `--keep-INFO TAG`, `--remove-INFO TAG` ✅
- **INFO extraction**: `--get-INFO TAG[,TAG]` → `.INFO` ✅
- **INFO key order in recode**: the recoded `.recode.vcf` INFO column now
  preserves the SOURCE order of INFO keys (e.g. `DP;AF`, not the formerly
  emitted alphabetical `AF;DP`). Upstream prints raw `INFO_str` verbatim for
  `--recode-INFO-all` (`vcf_entry.cpp:311`) and walks the parsed INFO vector
  in source order for `--recode-INFO TAG` (`get_INFO`, `entry_getters.cpp:182`).
  Note the two paths render a bare flag differently: `--recode-INFO-all` keeps
  it bare; `--recode-INFO TAG` emits `KEY=1` because `set_INFO`
  (`vcf_entry_setters.cpp:253`) materialises a flag's value as `"1"`. Pinned
  byte-for-byte vs live upstream by `TestVcftools_RecodeINFOOrderUpstreamParity`
  (fixture `info_order.vcf`, non-alphabetical source order + flag keys) and the
  binary-free `TestUnitFilterRecodeInfoOrder`. ✅

Closed in wave 2:

- **LDhat output formats**: `--ldhat`, `--ldhat-geno` (paired
  `<prefix>.ldhat.sites` / `<prefix>.ldhat.locs`, byte-for-byte vs
  upstream) ✅
- **Phased-site filter**: `--phased` (composes with `--ldhat` per
  upstream's `phased_only` invariant) ✅

Closed in wave 3:

- **LDhelmet output format**: `--ldhelmet` (paired
  `<prefix>.ldhelmet.snps` / `<prefix>.ldhelmet.pos`, byte-for-byte vs
  upstream; implies `--phased` + `--remove-indels` per
  parameters.cpp:275, requires `--chr` per parameters.cpp:717) ✅
- **IMPUTE reference-panel output**: `--IMPUTE` (case-sensitive; emits
  `<prefix>.impute.legend` / `<prefix>.impute.hap` /
  `<prefix>.impute.hap.indv`, byte-for-byte vs upstream; implies
  `--phased`, biallelic-only, no missing data per parameters.cpp:255) ✅

Closed in wave 4:

- **`--diff-indv-map FILE`** — two-column whitespace-separated file that
  renames file-2 sample IDs before matching against file-1. Loader
  mirrors upstream `variant_file_diff.cpp:11-34`; mapping is applied when
  forming `commonPairs` and when classifying `INDV FILES` in
  `.diff.indv_in_files`. ✅
- **`--diff-discordance-matrix`** — emits
  `<prefix>.diff.discordance_matrix` with the 5x5 layout from upstream
  `variant_file_diff.cpp:944-1198`: header row of file-1 genotype labels,
  four data rows of file-2 genotype labels, biallelic + matching-ALT only,
  diploid only, byte-for-byte parity vs upstream. ✅

Closed in wave 5:

- **`--diff-switch-error`** — emits `<prefix>.diff.switch` (per-event
  log with `CHROM POS_START POS_END INDV` columns) and
  `<prefix>.diff.indv.switch` (per-individual rate with
  `INDV N_COMMON_PHASED_HET N_SWITCH SWITCH` columns), byte-for-byte vs
  upstream. Ported from `variant_file_diff.cpp:1207-1507`
  (`output_switch_error`). Plugged into the existing diff runner so it
  composes with `--diff-indv-map` (sample-ID renaming) without
  re-implementing the load-file-2 path. ✅
- **`--mendel <PED>`** — emits `<prefix>.mendel` listing Mendelian
  inconsistencies across trios defined in a four-column PED file
  (`family child father mother`; first line always skipped). Ported from
  `variant_file_output.cpp:5332-5470`
  (`output_mendel_inconsistencies`). The PED column ordering follows
  upstream's `ss >> family >> child >> father >> mother;` parse;
  `family_ids` for each trio is `<child>_<father>_<mother>`. ✅

The original brief proposed `--phase` (output format) as the second
target. After re-checking upstream `parameters.cpp`, no standalone
`--phase` flag exists — only `--phased` (already ported in wave 2 as a
site filter). `--mendel` was substituted as the next clear long-tail
target, per the brief's substitution clause.

Closed in wave 6 (this PR):

- **`--non-ref-af FLOAT`** — minimum non-reference allele frequency:
  every ALT's count/non-missing-chr ratio must be ≥ threshold. Ported
  from upstream `parameters.cpp:303` + `entry_filters.cpp:770-824`. We
  preserve the upstream quirk that `min_non_ref_af > 0` also drops
  monomorphic (no-ALT) sites via the `N_failed == N_alleles-1`
  fallback on line 814. ✅
- **`--non-ref-ac INT`** — minimum non-reference allele count: every
  ALT's per-site count must be ≥ threshold. Ported from upstream
  `parameters.cpp:302` + `entry_filters.cpp:869-920`. The
  monomorphic-fallback on line 912 is gated on `min_non_ref_ac_any`
  (NOT plain `min_non_ref_ac`), so `--non-ref-ac` alone deliberately
  does NOT drop monomorphic sites — verified against upstream and
  pinned in `TestParity_NonRefAF_DropsMonomorphic`. ✅

The brief originally requested `--FILTER-PASS-summary` and
`--remove-INFO-all`. Neither flag is registered in upstream
`parameters.cpp`: the real flag is `--FILTER-summary` (already
implemented since wave 0), and `--remove-INFO-all` simply does not
exist (the upstream registrations are `--keep-INFO-all` and
`--recode-INFO-all`, both already implemented). The brief's
substitution clause permits picking from the `--non-ref-af*` family,
so `--non-ref-af` and `--non-ref-ac` were chosen as the two new
flags. The remaining `--non-ref-af-any` / `--non-ref-ac-any` (and the
`--max-*` upper-bound counterparts) are closed in wave 7 below.

Closed in wave 7 (this PR):

- **`--max-non-ref-af FLOAT`** — upper-bound counterpart of
  `--non-ref-af`: drops the site if ANY ALT's freq > threshold (per-ALT
  immediate fail, entry_filters.cpp:807). Also drops monomorphic sites
  via the line-814 fallback (gate keyed on plain thresholds). ✅
- **`--max-non-ref-ac INT`** — upper-bound counterpart of
  `--non-ref-ac`: drops the site if ANY ALT's count > threshold
  (entry_filters.cpp:905). Like plain `--non-ref-ac`, does NOT drop
  monomorphic sites (the line-912 fallback is keyed on `_any`). ✅
- **`--non-ref-af-any FLOAT`** + **`--max-non-ref-af-any FLOAT`** —
  N_failed-counter variants registered at `parameters.cpp:304-305` and
  `:289-290`. Wired for command-line parity, but observably **NO-OPS**
  when used alone because upstream `entry_filters.cpp:814` gates the
  fallback on the PLAIN thresholds (`min_non_ref_af > 0` /
  `max_non_ref_af < 1.0`), not on the `_any` thresholds. The flags
  only have an effect when paired with their plain counterpart, in
  which case the fallback fires when EVERY ALT fails the `_any`
  threshold. We mirror this verbatim; pinned by
  `TestParity_NonRefAFAny_NoOp` and `TestParity_NonRefAF_03_Any_06`. ✅
- **`--non-ref-ac-any INT`** + **`--max-non-ref-ac-any INT`** — counter
  variants of the AC family. Unlike AF, the AC fallback at
  `entry_filters.cpp:912` IS gated on the `_any` thresholds, so these
  flags are functional standalone. Site dropped when N_failed equals
  N_alleles-1 (every ALT failed). Monomorphic sites (N_alleles=1)
  satisfy 0==0 and are dropped — counter to plain `--non-ref-ac` which
  keeps them. Pinned by `TestParity_NonRefACAny_2`,
  `TestParity_NonRefACAny_1_Chr20`, `TestParity_MaxNonRefACAny_2_Chr20`. ✅

Refactor: the wave-6 per-ALT early-return was lifted to an N_failed
accumulator pass (matching upstream's structure literally) so the AF
and AC `_any` fallbacks can decide post-loop. The plain flags still
short-circuit on the first failing ALT.

Closed in wave 8 (this PR):

- **`--hwe FLOAT`** — minimum per-site exact-test (Wigginton/Cao/Abecasis
  2005) Hardy-Weinberg p-value filter. Ported from upstream
  `parameters.cpp:254` (which also forces `max_alleles = 2` when the
  flag is supplied) + `entry_filters.cpp:922-946` + `entry.cpp:18-101`
  (the exact-test). The SNPHWE port is a line-for-line port of
  upstream's integer-arithmetic midpoint/walk algorithm so that
  filtering decisions are byte-identical. The CLI adapter mirrors the
  upstream `--hwe → max_alleles=2` coupling. Pinned by
  `TestParity_HWE_005_sample` (3-sample fixture, exercises the
  max_alleles=2 coupling), `TestParity_HWE_005_fixture` (20-sample
  fixture engineered so site `1:200` with counts (10,0,10) fails the
  exact test at p≈1.34e-6 — verified via upstream `--hardy`), and unit
  tests on `snpHWE` boundary cases. ✅
- **`--max-missing-count INT`** — maximum number of missing
  *chromosomes* (haploid alleles, NOT samples) tolerated per site.
  Ported from upstream `parameters.cpp:286` +
  `entry_filters.cpp:918`. The comparator is strict `>` so a site
  with exactly `INT` missing alleles is KEPT, but `INT+1` drops.
  Pinned with two parity tests bracketing the boundary case
  (`TestParity_MaxMissingCount_1` drops a 2-missing site,
  `TestParity_MaxMissingCount_2` keeps it). The Params struct grows
  a `MaxMissingCountSet` boolean so the CLI can distinguish
  "user passed 0" (drop any site with any missing call) from
  "user omitted the flag" (no filter); the CLI registers the flag via
  `flag.Func` to record both. ✅

Closed in wave 9 (this PR):

- **`--kept-sites`** — emits `<prefix>.kept.sites`, a two-column
  `CHROM\tPOS` TSV listing every site that survived all filters in
  input order. Ported from upstream `parameters.cpp:268` +
  `variant_file_output.cpp:4285-4326` (`output_kept_sites`). Upstream
  re-parses the input file specifically for this output and re-runs
  `entry::apply_filters` (`entry_filters.cpp:23`) on every entry; we
  piggy-back on the existing filter pipeline in `Run` instead, which is
  equivalent because the same filter gates apply to both code paths.
  Header is `CHROM\tPOS` (tab-separated, LF terminator), matching
  upstream's `out << "CHROM\\t" << "POS" << endl;`. Pinned by
  `TestParity_KeptSites_NoFilter` (all-sites-pass case),
  `TestParity_KeptSites_HWE` (HWE filter), and
  `TestParity_KeptSites_PosFilter` (chr+from-bp+to-bp filter). ✅
- **`--removed-sites`** — counterpart of `--kept-sites`: emits
  `<prefix>.removed.sites` with the same column layout, listing sites
  *dropped* by any filter. Ported from upstream `parameters.cpp:330` +
  `variant_file_output.cpp:4328-4373` (`output_removed_sites`). Pinned
  by `TestParity_RemovedSites_HWE` and `TestParity_RemovedSites_PosFilter`.
  Plus `TestKeptRemoved_Disjoint_And_Complete` (port-only invariant —
  upstream forbids both flags in one run via `num_outputs > 1`; we
  deliberately do not replicate that constraint per the CLAUDE.md
  "don't replicate upstream bugs" rule, and instead verify that the
  combined invocation partitions the input perfectly) and
  `TestKeptRemoved_Disabled_NoFiles` (neither flag → no file leaked). ✅

Implementation note: the trace writer lives in
`tools/vcftools/pkg/vcftools/sitetrace.go`. Each `continue` in the main
filter loop of `vcftools.go:Run` now calls
`siteTracker.recordRemoved(chrom, pos)` immediately before bailing out,
and the successful path calls `siteTracker.recordKept(chrom, pos)` just
before `keptSites++`. Both methods are no-ops when the corresponding
flag is not set (cheap nil-check on the bufio.Writer).

Closed in wave 10 (this PR):

- **`--remove-filtered-geno-all`** — sets GT to `./.` for every kept
  genotype whose FORMAT FT field is not "PASS" or ".". Ported from
  upstream `parameters.cpp:323` + `vcf_entry.cpp:580-608`
  (`filter_genotypes_by_filter_status` with `remove_all=true`). Only
  the GT slot is rewritten; other FORMAT fields (FT/DP/GQ/...) pass
  through unchanged, matching upstream's recode emission at
  `vcf_entry.cpp:320-368`. Sites with no FT FORMAT column are
  left untouched (mirrors upstream's early-return at
  `entry_filters.cpp:94-108`). Pinned by
  `TestParity_RemoveFilteredGenoAll` (byte-for-byte vs upstream on the
  new `ft_geno.vcf` fixture) and `TestRemoveFilteredGeno_NoFT_NoOp`
  (port-only invariant against `sample.vcf`, which has no FT column). ✅
- **`--remove-filtered-geno NAME`** (repeatable) — drops a genotype
  whose FT lists any of the named flags. Ported from
  `parameters.cpp:324` + `vcf_entry.cpp:601-605`. FT is parsed as a
  `;`-separated list per upstream's `vcf_entry_setters.cpp:188-212`
  (entries equal to "" or "." are dropped from the list). Pinned by
  `TestParity_RemoveFilteredGenoQ10` (single-flag) and
  `TestParity_RemoveFilteredGenoMulti` (two-flag invocation — the
  set behaviour in upstream's `geno_filter_flags_to_exclude`). ✅
- **`--max-indv N`** — caps the number of kept individuals at N.
  Ported from upstream `parameters.cpp:292` +
  `variant_file_filters.cpp:105-147` (`filter_individuals_randomly`).
  Upstream draws the random subset with glibc `srand()/rand()` driving
  libstdc++'s `std::random_shuffle`, but seeds with `srand(time(NULL))`
  and exposes **no** seed flag, so a plain upstream run's chosen subset
  changes every wall-clock second and is **not reproducible** (verified
  by running the upstream binary in consecutive seconds). This port
  therefore ports the exact glibc `rand()` generator (`glibc_rand.go`)
  and the exact `random_shuffle` swap sequence, and adds a
  **`--max-indv-seed N`** flag. With a seed the kept subset is
  **byte-for-byte identical** to what upstream would emit if seeded with
  the same value; without one it keeps the first N samples in header
  order (no reproducible upstream target exists). Byte-exact parity is
  pinned by `TestParity_MaxIndv_Seeded` (fixtures generated by a C
  harness that replays `filter_individuals_randomly` verbatim and was
  itself confirmed against the live upstream binary at matching epoch
  seconds), plus `TestGlibcRand*` for the generator and
  `TestVcftools_MaxIndvUpstream` for the live count/subset contract. ✅
- **`--keep-INFO-all`** — upstream-deprecated synonym for
  `--recode-INFO-all`. Both `parameters.cpp:267` and `:318` write
  to the same `recode_all_INFO` parameter bit; the CLI ORs them
  together so either flag (or both) produces identical output.
  Pinned by `TestKeepINFOAll_Synonym`. ✅
- **`--version`** — prints `VCFtools (0.1.18)` and exits, matching
  upstream `parameters.cpp:648-652` byte-for-byte. Hard-coded
  version string tracks the upstream submodule at port time;
  bump on rebase.  ✅

Implementation notes (wave 10):

- The FT-based filter lives inside `filterGenotypes` in
  `tools/vcftools/pkg/vcftools/vcftools.go` (next to the existing
  `--minDP/--maxDP/--minGQ` path). Two small helpers (`parseSampleFT`,
  `shouldDropByFT`) keep the hot loop branch-light and mirror the
  upstream getters/parsers documented above. The `sampleWithMissingGT`
  helper extracted from the duplicated DP/GQ paths is also reused by
  the FT path.
- `--max-indv` is wired through `buildSampleFilter`, which now returns
  a non-nil keep set when `MaxIndvSet` is true even with no
  identity-based filter. The cap iterates `header.Samples` in order
  so the truncation is deterministic.

Closed in wave 11 (this PR):

- **`--mask FILE`** — FASTA-style positional mask. The mask file has
  `>CHROM` headers followed by lines of digit characters (one per
  reference base, 1-based). A site at (CHROM, POS) is kept when its
  mask digit is `<= --mask-min` (default 0). Ported from upstream
  `parameters.cpp:280` + `entry_filters.cpp:674-752`
  (`filter_sites_by_mask`). Pinned by `TestParity_Mask_Default`,
  `TestParity_Mask_Min5`, `TestParity_Mask_Partial`. ✅
- **`--invert-mask FILE`** — same loader as `--mask`, but the keep/drop
  decision is flipped. Ported from upstream `parameters.cpp:262`.
  Pinned by `TestParity_InvertMask_Min5`. ✅
- **`--mask-min INT`** — maximum kept mask digit value, 0-9 (upstream
  errors when `> 9` at `parameters.cpp:720`; we additionally reject
  negatives at load time because they silently drop every site
  upstream — clearer to fail fast). Default 0. ✅

Implementation notes (wave 11):

- The mask reader is **forward-only**, mirroring upstream's stateful
  `ifstream` walk (entry_filters.cpp:680-688 keeps `mask_chr`,
  `mask_line`, and `mask_pos` as static state across calls). The Go
  port loads the file once into a `[]maskChromosome` and maintains a
  `(chromIdx, slabIdx)` cursor. The cursor never moves backwards, so a
  VCF presenting chr2 before chr1 against a mask listing chr1 then
  chr2 will lose the chr1 sites — this matches upstream behaviour
  (`TestMaskFilter_OutOfOrderVCFDrops` pins it).
- Header tokenisation matches upstream's
  `line.substr(1, line.find_first_of(" \t")-1)` (split on first
  whitespace after `>`; comments are discarded).
- Mutually-exclusive flag behaviour: `--mask` and `--invert-mask`
  share the same `mask_file` slot upstream (last one wins via
  parameters.cpp:262 vs :280). Go's `flag.String` does not preserve
  last-set ordering, so we apply the asymmetric rule "if
  `--invert-mask` is non-empty, override and set invert=true". This is
  observable only when both flags are supplied; documented in
  `main.go`.

Closed in wave 12 (this PR):

- **`--positions-overlap FILE`** — keep a record when ANY base in
  `[POS, POS+len(REF)-1]` matches a (CHROM, POS) entry in the file.
  Same two-column whitespace-separated format as `--positions`; `#`
  comments and blank lines tolerated. Ported from upstream
  `parameters.cpp:315` + `entry_filters.cpp:408-531`
  (`filter_sites_by_overlap_positions`, keep-branch). For 1-base REF
  records the behaviour reduces to plain `--positions`; the divergence
  appears on indels / MNPs, which is the entire reason upstream ships
  the overlap variant. Sites on chromosomes not named in the file are
  dropped (matches `entry_filters.cpp:515-516`). ✅
- **`--exclude-positions-overlap FILE`** — drop a record when ANY base
  in `[POS, POS+len(REF)-1]` matches a (CHROM, POS) entry. Ported from
  `parameters.cpp:221` + `entry_filters.cpp:533-547`. Inverse of
  `--positions-overlap`; sites on chromosomes not named pass through
  unchanged (matches the `chr_to_idx.find != end` guard at line 535). ✅

Implementation notes (wave 12):

- Upstream reuses the same `keep_positions`/`exclude_positions` state
  for both the plain `--positions` family (entry_filters.cpp:279-406)
  and the overlap family (lines 408-548). Consequence: combining
  `--positions` with `--positions-overlap` upstream silently degrades
  to overlap semantics for whichever file populated the set first
  (the second loader sees a non-empty set and skips). We keep the
  four filters as independent fields on `Params` so each behaves
  exactly as documented when used solo, and when combined both gates
  apply (a site must pass include AND not be excluded across the two
  flag pairs). Pinned by
  `TestPositionsOverlap_VsPlain_DivergesOnMultiBaseRef`.
- The upstream loop is half-open: `for ui=POS; ui<POS+REF.size(); ui++`.
  We mirror this with `for p := v.Pos; p < v.Pos+refLen; p++`. With a
  defensive guard `refLen=max(len(v.Ref),1)` we ensure at least the
  POS itself is tested for malformed VCFs with empty REF; valid VCFs
  always have `len(REF) >= 1` and the guard never triggers.
- File format reuses the existing `loadPositions` parser (same as
  `--positions`); separate `positionSet` instances are kept in `Run`
  for the four flags so the apply order is deterministic
  (`includePos` → `excludePos` → `includePosOverlap` →
  `excludePosOverlap` → BED → variant-type → frequency / quality /
  HWE / etc.).

Closed in wave 13 (this PR):

- **`--derived`** — when combined with `--freq` / `--counts`, reorder
  the allele columns so that the ancestral allele (INFO/AA,
  case-insensitive) appears first; drop sites where AA is missing,
  `.`, `?`, or does not match REF/ALT. Mirrors upstream
  `parameters.cpp:201` + `variant_file_output.cpp:67-159`
  (`output_frequency`, the `derived` branch). Implementation lives in
  `addFrequencyStat` (new `derivedSwap` flag on `siteFreqStat`) and
  the existing `outputFrequency` reorder loop. Multi-allelic sites
  are already dropped by our biallelic-only `--freq` restriction, so
  `--derived` only affects the biallelic subset (matches the subset
  this port emits at all under `--freq`/`--counts`). ✅
- **`--extract-FORMAT-info NAME`** — extract a per-genotype FORMAT
  field across all kept samples into a tab-separated
  `<prefix>.<NAME>.FORMAT` file. Header is `CHROM\tPOS\t<sample>...`;
  one data row per site whose FORMAT column lists NAME (sites lacking
  NAME in FORMAT are skipped entirely, matching upstream's
  `FORMAT_id_exists` gate). Samples whose colon-separated value
  vector is too short to reach NAME's index emit `.` (matches
  `vcf_entry.cpp:618` + the early `break` at line 637). Ported from
  upstream `parameters.cpp:222` +
  `variant_file_format_convert.cpp:1204-1263`
  (`output_FORMAT_information`). Single-valued upstream (the last
  value wins on the CLI). ✅

Implementation notes (wave 13):

- `--derived` is a *modifier*, not an output: it only takes effect
  when paired with `--freq` or `--counts`. Pinned by
  `TestDerived_NoFreqIsNoOp`. Upstream's
  `parameters.cpp:201` only flips a boolean; the reorder logic lives
  inside `output_frequency` (the boolean is also consumed by
  `output_indv_burden` and `output_indv_freq_burden`; both burden
  flags now honour `--derived` after wave 14 — see the wave-14
  section below).
- AA uppercasing — upstream calls
  `std::transform(AA.begin(), AA.end(), AA.begin(), ::toupper)` on
  line 78 / 439 / 564 before comparing against `e->get_allele(ui)`.
  We replicate this with `strings.ToUpper` on both AA and REF/ALT so
  `AA=a, REF=A` matches as expected. The case-insensitive match is
  pinned by the 1:400 site in `derived_fixture.vcf`.
- AA sentinels — upstream's `if ((AA == "?") || (AA == "."))` check
  appears at lines 79 / 440 / 565. We mirror it explicitly (empty,
  ".", "?"). Pinned by sites 1:500 / 2:100 in the same fixture.
- `--extract-FORMAT-info` shares its FORMAT-tag presence helper
  (`formatContains`) with the BEAGLE writer (declared once in
  `beagle.go`).
- The Go VCF parser only populates `sample.Data[key]` for keys whose
  colon-token slot exists in the per-sample string; absent slots
  therefore read back as missing-from-map. We treat that as upstream's
  "value vector too short → '.'" case (vcf_entry.cpp:618-637). Pinned
  by sites 1:100/S3, 1:400/S2 in `extract_format_fixture.vcf`.

Closed in wave 14 (this PR):

- **`--indv-burden`** — per-individual diploid-burden counts emitted to
  `<prefix>.iburden`. Header is
  `INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS`; with `--derived` the
  `_REF` / `_ALT` columns rename to `_ANC` / `_DER`. Non-diploid sites
  are skipped (upstream's `if (e->is_diploid() == false) continue;` at
  `variant_file_output.cpp:429-433`). With `--derived` the site's
  INFO/AA tag picks the ancestral-allele index; sites where AA is
  missing, `.`, `?`, or does not match any REF/ALT are skipped. Ported
  from upstream `parameters.cpp:257` +
  `variant_file_output.cpp:378-498` (`output_indv_burden`). ✅
- **`--indv-freq-burden`** — per-individual × per-allele-count matrix
  written to `<prefix>.ifreqburden`. For each kept diploid site,
  computes the per-allele count vector across kept individuals and
  for each kept individual increments the burden cell at column
  `allele_counts[geno_allele]` for each non-ref (or non-ancestral with
  `--derived`) allele the individual carries. Mirrors upstream
  `parameters.cpp:258` + `variant_file_output.cpp:501-627`
  (`output_indv_freq_burden` with `double_count_hom_alt=0`). ✅
- **`--indv-freq-burden2`** — same as `--indv-freq-burden` but with
  `double_count_hom_alt=1`, so a hom-alt genotype contributes 1 (not
  2) to the corresponding allele-count bin. Mirrors `vcftools.cpp:64`
  - the same `output_indv_freq_burden` routine. ✅

Implementation notes (wave 14):

- All three flags share `burden.go`: `indvBurdenRunner` for
  `--indv-burden` and `indvFreqBurdenRunner` for the two
  freq-burden variants (the latter takes a `doubleCountHomAlt`
  toggle for `--indv-freq-burden2`).
- **Upstream label-index bug preserved.**
  `output_indv_freq_burden` writes
  `out << meta_data.indv[indv_count];` at line 621 — that should be
  `meta_data.indv[ui]` (the original-index). With `--remove-indv` (or
  any sample-filter) dropping a non-trailing sample, the labels in
  `.ifreqburden` shift relative to the burden values (e.g. excluding
  `S2` from `[S1,S2,S3,S4]` yields labels `S1,S2,S3` for kept
  individuals `S1,S3,S4`). We mirror this byte-for-byte;
  pinned by `TestParity_IndvFreqBurden_LabelBug`. The
  `output_indv_burden` function at lines 488-497 correctly uses
  `meta_data.indv[ui]` and is not affected.
- The diploid-only check is per-site: any haploid genotype in a kept
  individual disqualifies the entire site. Upstream emits a one-off
  warning at `variant_file_output.cpp:431`; we do not re-emit the
  warning byte-for-byte but skip the site identically. Pinned by
  `TestIndvBurden_SkipsNonDiploid`.
- The `ancestralAlleleIndex` helper centralises the
  `variant_file_output.cpp:437-462` AA-resolution logic (uppercase
  match, missing-sentinel handling) so the same predicate is used by
  both burden runners and stays consistent with how
  `addFrequencyStat` already implements the `--derived` filter (see
  wave 13).
- The `diploidAlleles` helper splits a `a/b` / `a|b` GT into
  `(a1, a2)`, treating `.` as -1, and is shared by both runners. The
  existing `parseGTForLDhat` parser was close but its haploid branch
  silently coerces a missing single-allele GT into a `(−1, −2,
  true)` triple that the burden routines need to reject as
  non-diploid; `diploidAlleles` returns `ok=false` for any GT without
  a separator so the caller can do the diploid-skip up front.

**PCA family**: ✅ resolved in wave 19. `--pca`, `--pca-no-norm`,
and `--pca-snp-loadings INT` (upstream `parameters.cpp:308-310`,
`variant_file_output.cpp:4871-5249`) all land via gonum's symmetric
eigensolver (`gonum.org/v1/gonum/mat`'s `SymEigen`). The owner
sanctioned `gonum` as the second third-party-dep zone after the
CRAM codec carveout (CLAUDE.md). Parity goldens were generated by
rebuilding upstream vcftools with `--enable-pca` against system
`liblapack-dev` / `libblas-dev`; the resulting binary lives at
`/tmp/vcftools_lapack_install/bin/vcftools` in the dev sandbox.
Eigenvector signs are LAPACK-implementation-dependent (both
LAPACK and gonum), so parity tests use a per-column sign-tolerant
comparison; eigenvalues themselves are sign-invariant and compared
with a tight numerical tolerance. Wave 19 also fixes a latent
upstream bug — `output_PCA` reads past the end of the per-individual
M[i] vectors when any kept individual has a missing genotype, see
`docs/UPSTREAM_BUGS.md` for the writeup. Pinned by
`TestParity_PCA_Basic`, `TestParity_PCA_NoNorm`,
`TestParity_PCA_SNPLoadings`, and the in-pkg algebraic tests in
`pca_test.go`.

Closed in wave 15 (this PR):

- **`--hapcount BED`** — per-BED-bin haplotype-count summaries written
  to `<prefix>.hapcount`. Columns: `#CHROM BIN_START BIN_END N_SNP
  N_UNIQ_HAPS N_GROUPS {MULTIPLICITY:FREQ}...`. Implies `--phased`
  (upstream parameters.cpp:248 sets `phased_only=true`). Diploid-only
  per-site (upstream variant_file_output.cpp:1350-1354 skip otherwise).
  Bins must be non-overlapping (upstream errors otherwise at lines
  1208-1216). Ported from variant_file_output.cpp:1169-1401
  (`output_haplotype_count`) with three upstream bugs FIXED on port
  per CLAUDE.md and the wave-14 precedent (PR #138). See
  `docs/UPSTREAM_BUGS.md#fix-on-port-resolved` for the writeup:
    1. prev_bin_idx shift on within-chromosome bin transitions
       (lines 1314-1315) — old bin's counts were silently overwritten
       with the new bin's values.
    2. End-of-stream read-after-free (lines 1370-1400) — last
       chromosome's rows were silently dropped (or zeroed) on a
       glibc-built upstream binary.
    3. BED first-line unconditional skip (line 1183) — header-less
       BEDs silently lost one bin per invocation.
  Pinned by `TestHapcount_CorrectBinTransitions`,
  `TestHapcount_EndOfStreamFlush`, and `TestHapcount_BEDFirstLineWithData`
  in `tools/vcftools/pkg/vcftools/hapcount_test.go`. The
  `.expected.hapcount` fixture is hand-traced, NOT generated from
  the upstream binary (whose output is wrong). ✅
- **`--temp DIR`** — upstream parameters.cpp:341 stores DIR as the
  base path for `mkstemp` spill files used by the LD and
  format-convert paths (variant_file_output.cpp:1441,
  variant_file_format_convert.cpp:28/402/627/810/994). This port does
  not spill to disk for any of those paths so the flag is accepted
  for CLI parity but has no observable effect; `Run` logs to stderr
  that the value was parsed-but-unused. Pinned by
  `TestParams_TempDirAccepted`. ✅
- **`--gzdiff FILE`** — upstream parameters.cpp:237 sets
  `diff_file = FILE; diff_file_compressed = true;`, and
  vcf_file.cpp:21 then switches to the gzip reader. This port's
  `iohelper.OpenReader` already auto-sniffs gzip from the magic
  bytes, so `--gzdiff` is wired as a plain alias for `--diff`
  (last-set wins, matching upstream's shared `diff_file` slot
  semantics parameters.cpp:209 vs :237). Pinned by
  `TestGzdiffAliasesDiff`. ✅

Closed in wave 16 (this PR):

- **`--recode-INFO TAG`** — upstream-canonical name for the repeatable
  recode-INFO-column selector (parameters.cpp:319 →
  `recode_INFO_to_keep.insert(...)`). The port already implemented
  this semantic under the `--keep-INFO TAG` flag name; wave 16 adds
  the canonical spelling as a synonym that funnels into the same
  `keepINFOParts` slice in `tools/vcftools/cmd/vcftools/main.go`,
  matching the existing `--keep-INFO-all` ↔ `--recode-INFO-all`
  pattern. Pinned by `TestCLI_RecodeINFOAlias`,
  `TestCLI_RecodeINFOAlias_Repeatable`, and
  `TestCLI_RecodeINFOAlias_MixedWithKeepINFO` in
  `tools/vcftools/cmd/vcftools/aliases_cli_test.go`. ✅
- **`-c`** — upstream's short alias for `--stdout`
  (parameters.cpp:194). Wired with `flag.BoolVar` pointing at the
  same `useStdout` bool, so either spelling toggles streaming
  output. Pinned by `TestCLI_ShortStdoutFlag`. ✅

Wave 16 also enumerated the FULL upstream long-flag table for the
first time (146 distinct flags). Two flags were initially flagged as
gaps by an earlier wave's regex but actually exist in the port under
non-obvious wiring: `--help` (registered via `flag.BoolVar(help, "help"
…)`, not `flag.Bool("help", …)`). The remaining-gap set below is the
definitive list.

Note on the **`--keep-INFO` semantic gap**: ✅ resolved in wave 17.
Upstream's `--keep-INFO` (parameters.cpp:266 →
`site_INFO_flags_to_keep` → `entry_filters.cpp:1033-1063`) is a SITE
FILTER. Pre-wave-17 the Go port had this flag wired to the
recode-INFO-column selector semantic. Wave 17 separates the two:
`--keep-INFO TAG` now drives `Params.KeepINFO` (site filter that
errors on non-Flag tags, OR-composing across multiple tags), while
`--recode-INFO TAG` drives the new `Params.RecodeINFO` field
(recode-column selector). See `docs/UPSTREAM_BUGS.md` Fix-on-port
section for the migration note. The sibling `--remove-INFO`
divergence was resolved in wave 18 (see below).

**`--remove-INFO` semantic gap**: ✅ resolved in wave 18.
Upstream's `--remove-INFO` (parameters.cpp:328 →
`site_INFO_flags_to_remove` → `entry_filters.cpp:1068-1086`) is a
SITE FILTER that drops sites where the named Flag IS present
(OR-veto across multiple tags). Pre-wave-18 the Go port had this
flag wired as a recode-column stripper — a port-only invention with
no upstream equivalent. Wave 18 repoints `--remove-INFO TAG` at
`passRemoveINFOSite`, the polarity-inverted complement of wave 17's
`passKeepINFOSite`. Header validation (Type=Flag check) is shared
with the keep path via `validateFlagTypeINFO`. Composition with
`--keep-INFO` follows upstream's keep-then-remove ordering:
`TestRun_KeepAndRemoveINFO_Compose` and
`TestParity_KeepAndRemoveINFO_Compose` pin this. The dead
recode-column stripper code in `filterRecodeInfo` was deleted (no
CLI flag drives it now). See `docs/UPSTREAM_BUGS.md` Fix-on-port
section for the full migration note.

Flag history (definitive enumeration vs.
`reference_code/vcftools/src/cpp/parameters.cpp` — wave 16):

As of wave 16 the diff between upstream's `in_str == "--…"` table
(146 unique long flags) and the port's registered flags was **five
flags** (the PCA-deferred trio counted as registered). All five are
now closed — waves 19–23, recorded below for history:

- ~~**`--bcf` FILE** — BCF binary input (parameters.cpp:173).~~
  **Closed (wave 22).** Adapts `pkg/htsgo/bcf.Reader` (built on
  the in-tree BGZF reader and on the wave-21 BCF dictionary fixes)
  into the new `variantSource` interface so `Run` iterates BCF
  records through the same filter pipeline used by `--vcf`. Pinned
  by `TestRun_BCFInput_Roundtrip` (write-then-read symmetry) and
  `TestRun_BCFInput_ComposesWithFilters` (BCF input composes with
  `--chr` / `--from-bp`).
- ~~**`--diff-bcf` FILE** — second-file BCF input for `--diff-*`
  family (parameters.cpp:210).~~ **Closed (wave 23).**
  `loadDiffBCF` mirrors `loadDiffVCF` but routes through the
  wave-22 `bcfVariantSource` (BGZF + `bcf.Reader` + ToVariant).
  The shared `loadDiffFromSource` body drives both loaders so
  the (CHR,POS)-keyed `diffData` build is identical between VCF
  and BCF second files. CLI mutual-exclusion with `--diff` /
  `--gzdiff` mirrors upstream's last-set-wins slot semantics.
  Pinned by 5 tests including all five `--diff-*` outputs
  composed against a BCF second file.
- ~~**`--recode-bcf`** — emit BCF instead of VCF (parameters.cpp:317).~~
  **Closed (wave 21).** Layered the existing `pkg/htsgo/bcf.Writer`
  on top of `pkg/bgzip.NewWriter` and wired it parallel to the
  existing `--recode` text-VCF path; both flags may be combined.
  This wave also fixed three latent bugs in the shared BCF writer
  uncovered by interop with upstream's reader: (1) missing
  IDs were encoded as type-0 instead of zero-length-typed-char;
  (2) the unified INFO+FILTER+FORMAT dictionary numbering wasn't
  surfaced — entries now carry their IDX and the text header is
  emitted with the `,IDX=N` annotations htslib uses; (3) FORMAT
  field descriptors were encoded with the total flat length as
  `size` instead of the per-sample dimension. Pinned by
  `TestRun_RecodeBCF_Roundtrip` and
  `TestRun_RecodeBCF_HeaderHasIDXAnnotations`, plus interop
  testing with upstream vcftools 0.1.18 (decoded our `.recode.bcf`
  through `vcftools --bcf` → byte-identical VCF round-trip).
- ~~**`--contigs` FILE** — supplemental `##contig=` lines for BCF
  header construction (parameters.cpp:197 →
  `variant_file.cpp:45-69`).~~ **Closed (wave 22).**
  `augmentHeaderContigs` prepends `##contig=<ID=...>` MetaInfo
  lines to the parsed header when the source lacks contig
  declarations of its own (matching upstream's
  `has_contigs == false` gate). Accepts both bare contig names
  and full `##contig=<...>` forms. Pinned by
  `TestRun_ContigsFile_AddsContigLines`,
  `TestRun_ContigsFile_NoOpWhenHeaderAlreadyHasContigs`,
  `TestAugmentHeaderContigs_AcceptsMetaInfoForm`.
After wave 23 vcftools reaches **146/146 long flags** — the
complete `parameters.cpp` surface is exercised. The
remaining work is per-output column-set polish (see the
"Other" list below) and the multi-output `num_outputs > 1`
check upstream uses, which we deliberately don't replicate
per the CLAUDE.md "don't replicate upstream constraints" rule.

Other (per-output column-set gaps, not flag-count gaps):

- ~~**Per-individual output**: the per-individual `.imiss` row layout
  still has fields we don't emit (we have `--missing-indv`).~~
  **Closed.** `--missing-indv` now emits the upstream 5-column
  layout `INDV N_DATA N_GENOTYPES_FILTERED N_MISS F_MISS` with
  upstream-exact values: N_GENOTYPES_FILTERED counts genotypes dropped
  by a genotype-level filter (`--minDP/--maxDP/--minGQ/--remove-filtered-geno*`)
  and excludes them from N_DATA/N_MISS, matching
  `output_indv_missingness` (variant_file_output.cpp:776-845). Missing
  calls follow upstream's first-allele-only rule (`alleles.first == -1`),
  so `./1` counts as missing but `0/.` does not. F_MISS uses libstdc++
  `%g` formatting (six significant digits; `-nan` for a 0/0 ratio). The
  underlying `filterGenotypes` DP/GQ paths were also brought into line
  with upstream's `DP_idx/GQ_idx != -1` site-FORMAT gate and its
  missing-value-as-(-1) semantics. Pinned by the live-binary
  `TestVcftools_MissingIndvUpstreamParity` (no-filter, `--minGQ`,
  `--minDP` cases), plus `TestMissingIndv_GenoFilteredColumn` and
  `TestGenotypeIsMissing` unit tests.
- ~~**Diff family**: per-site / per-indv discordance outputs
  (`--diff-site-discordance`, `--diff-indv-discordance`) still emit a
  simpler column set than upstream's richer `.diff.sites` /
  `.diff.indv` schemas — see `variant_file_diff.cpp:635` for the gap.~~
  **Closed (wave 20).** `.diff.sites` now emits the upstream 7-column
  layout (`CHROM POS FILES MATCHING_ALLELES N_COMMON_CALLED N_DISCORD
  DISCORDANCE`), including file-1-only and file-2-only zero rows;
  `.diff.indv` now emits the 4-column layout (`INDV N_COMMON_CALLED
  N_DISCORD DISCORDANCE`) over the union of file-1 and effective
  file-2 samples in alphabetical order. Discordance values format
  via `%.6g` with `-nan` for 0/0, matching libstdc++'s default
  ostream output. Pinned by five parity tests against upstream
  goldens (`TestParity_DiffSiteDiscordance_{NoMap,WithMap,AltMismatch}`,
  `TestParity_DiffIndvDiscordance_{NoMap,WithMap}`). Residual
  deviations from upstream:
  1. Row ordering within `.diff.sites` follows file-1 streamed
     order with file-2-only sites appended in sorted-chrom-then-pos
     order, rather than upstream's strict merge sort — observable
     only when the two files have non-overlapping positions
     interleaved by chromosome.
  2. **REF-mismatch shared sites**: upstream
     `variant_file_diff.cpp:787-790` SKIPS B-sites where REFs
     differ (with a one-off `"Non-matching REF"` warning), treating
     the site as if it weren't shared. The Go port emits the row
     with `MATCHING_ALLELES=0` and accumulates discordance over
     it. Tracked as a separate follow-up.
  3. **REF=N/`.`/empty normalisation**: upstream replaces
     `REF1` with `REF2` (and vice versa) when one side is `N`,
     `.`, or empty before the alleles-match check
     (`variant_file_diff.cpp:780-783`). The Go port does a
     verbatim string compare. Same follow-up as (2).
  4. **`-nan` literal portability**: the port hardcodes
     `"-nan"` for division-by-zero discordance to match
     libstdc++'s default ostream output on glibc. A future
     golden regenerated on musl / macOS libc / MSVC could
     print `nan` or `NaN` instead. Re-running the goldens on
     a non-glibc system would require updating the literal.
- ~~**`--freq` / `--counts` float formatting**: the `.frq` allele-frequency
  column printed `%.6f` (`0.500000`), but upstream's `output_frequency`
  (`variant_file_output.cpp:131`) writes each freq straight to a default
  C++ ostream — `defaultfloat` with precision 6, i.e. six *significant*
  digits with trailing zeros stripped (`0.5`, `0.0833333`, `1`, `0`).~~
  **Closed.** A `formatFreq` helper (`statistics.go`) uses Go's
  `strconv.FormatFloat(v, 'g', 6, 64)`, which reproduces the C++ default
  ostream output byte-for-byte (verified over a 12.5M-ratio sweep against a
  live `ostringstream`, including the sub-`1e-4` scientific-notation
  threshold and round-half-to-even at the 6th digit). `--counts` was already
  integer-formatted and is unaffected. Pinned by the live-binary
  `TestVcftools_FreqUpstreamParity` (`--freq` + `--counts` byte-for-byte on
  the all-biallelic `freq_fmt_fixture.vcf`); the `freq.expected.frq` and
  `derived.expected.frq` goldens were regenerated to the upstream format.
  The same fix was applied to `outputFrequency2` (`--freq2`).
- ~~**`--freq2` / `--counts2` output schema**~~ **Closed.** Previously the Go
  port wrote `<prefix>.frq2` / `.frq.count2` with a PLINK-style
  `CHROM POS N_CHR REF_FREQ ALT_FREQ` layout. It now mirrors upstream's
  `suppress_allele_output` branch (parameters.cpp:198-225,
  variant_file_output.cpp:42-156): `--freq2` / `--counts2` write the SAME
  `.frq` / `.frq.count` files as `--freq` / `--counts` with the allele labels
  stripped — header `{FREQ}` / `{COUNT}` and bare tab-separated values — and
  the single `suppress_allele_output` toggle (set by either "2" flag) applies
  to all frequency/count output. The dedicated `outputFrequency2` /
  `outputCounts2` writers were removed in favour of a `suppressAlleles`
  parameter on `outputFrequency`. Pinned by `--freq2` / `--counts2` cases in
  the live-binary `TestVcftools_FreqUpstreamParity`.
- **Other**: small-format columns gaps tracked in
  `tools/PORTING_STATUS.md`.

Note: the brief mentioned `--haploid` as a possible wave-2 target. After
checking the upstream source (`reference_code/vcftools/src/cpp/`) there is
no `--haploid` flag — the closest thing is `--phased` (parameters.cpp:311
and entry_filters.cpp:989-1010), which we ported instead.

**Validation:** wave 1 adds header byte-for-byte parity tests for the new
output files; wave 2 ships full byte-for-byte parity tests for both
`.ldhat.sites` and `.ldhat.locs`; wave 3 adds byte-for-byte parity tests
for `.ldhelmet.snps` / `.ldhelmet.pos` and the IMPUTE bundle
(`.impute.legend` / `.impute.hap` / `.impute.hap.indv`); wave 4 adds
byte-for-byte parity tests for `.diff.discordance_matrix` (with and
without `--diff-indv-map`) and the mapped `.diff.indv_in_files` output
against upstream goldens (under `tools/vcftools/testdata/parity/`,
fixtures `diff_f1.vcf` / `diff_f2.vcf` / `diff_indv_map.txt`); wave 5
adds byte-for-byte parity tests for `.diff.switch` /
`.diff.indv.switch` (fixtures `switch_f1.vcf` / `switch_f2.vcf`) and
`.mendel` (fixtures `mendel.vcf` / `mendel.ped`); wave 6 adds
byte-for-byte parity tests for `--non-ref-af 0.3` / `--non-ref-af 0.5`
on `sample.vcf` and `--chr X --non-ref-ac 2|3` on the same fixture
(see `non_ref_af_*.expected.recode.vcf` and
`non_ref_ac_*_chrX.expected.recode.vcf`), plus a regression test that
pins the upstream `_any`-fallback asymmetry between the two flags;
wave 7 adds byte-for-byte parity tests for `--non-ref-ac-any 2`
(full VCF), `--chr 20 --non-ref-ac-any 1`,
`--chr 20 --max-non-ref-af 0.3`, `--chr 19 --max-non-ref-ac 2`,
`--chr 20 --max-non-ref-ac-any 2`, and
`--non-ref-af 0.3 --non-ref-af-any 0.6` (the only meaningful AF-any
usage), plus a `TestParity_NonRefAFAny_NoOp` regression that pins
upstream's documented "AF -any flags are no-ops alone" quirk by
asserting the port produces baseline output unchanged.
Full upstream-test-suite run still pending.

Upstream build note for golden generation: vcftools'
`variant_file_format_convert.cpp` LDhat writers allocate a stack array of
exactly `temp_dir.size()` bytes and then `strcpy` a longer string into
it (`new_tmp = temp_dir + "/vcftools.XXXXXX";` plus
`char tmpname[new_tmp.size()]; strcpy(tmpname, new_tmp.c_str());` at
lines 627-629 and again at 813-815, 658-665, 680-688, 841-845). Modern
glibc fortified strcpy aborts the run. To regenerate the goldens locally:

```
cd reference_code/vcftools/src/cpp
make clean
make CXXFLAGS='-O0 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0'
```

Filed as a parity-only workaround; we don't replicate the bug.

### `bgzip`

**Status:** 1 / 1 command, most flags.

Done:

- **`.gzi` random-access seek** (`bgzip -b N -s M`). `pkg/htsgo/bgzf`
  gains `SeekReader` (mirrors htslib `bgzf_useek`): given an
  `io.ReaderAt` over a BGZF stream and a `.gzi` block index, it
  binary-searches the index for the block owning a requested
  uncompressed offset, inflates only the overlapping blocks, and reads a
  region without decompressing the whole file. A sparse index is
  tolerated (it walks the compressed stream forward for unindexed
  blocks). Our existing `WriteGZI`/`ReadGZI` already matched htslib's
  on-disk `.gzi` format; the seek path is the missing consumer.
  **Validated live** against the upstream `bgzip` binary
  (`pkg/htsgo/bgzf/seek_test.go`): our `SeekReader` reproduces
  `bgzip -b/-s` region extraction byte-for-byte, our `WriteGZI` output is
  byte-identical to `bgzip -r`, and upstream `bgzip -b/-s` extracts
  correctly using a `.gzi` we wrote.

- **Multi-threaded compression** (`-@ / -t / --threads N`) — DONE.
  BGZF is block-parallel by construction, so `pkg/htsgo/bgzf.MultiWriter`
  deflates blocks across N worker goroutines and a single ordering
  collector writes the framed blocks back in their original sequence,
  yielding a valid BGZF stream (canonical EOF block, decodable blocks)
  for any thread count. Single- and multi-threaded output are
  byte-identical at equal block boundaries, and both decompress to the
  same plaintext. The `.gzi` index stays correct because block sizes and
  order are preserved.
- **Output-rename / `-o, --output FILE`** — DONE. Reading stdin without
  `-c` writes BGZF to stdout; `-o FILE` names the output explicitly (the
  upstream convention for naming stdin output). `-o -` is normalised to
  `-c` (stdout), matching htslib.

**Validation:** live parity against htslib's `bgzip`
(`TestBgzip_ThreadsUpstreamParity`, on-demand submodule build): our
`-@N` output decompresses to byte-identical plaintext via both upstream
`bgzip -d` and our own reader; upstream `bgzip -@N` output decompresses
via our reader; single- vs multi-thread outputs decompress-equal; and
the stream is structurally valid. Note: multi-threaded compressed bytes
are not guaranteed byte-identical to upstream (block boundaries differ),
so parity is asserted on the recovered plaintext plus structural
invariants. The `.gzi` seek path is validated byte-for-byte against the
live upstream `bgzip` binary (see above). Round-trips through `tabix`
continue to work. No full upstream-test suite run yet.

### `tabix`

**Status:** 1 / 1 command, most flags.

Done:

- **`--reheader FILE`** (`-r`) — replaces the leading meta-char header of a
  bgzipped file with the contents of FILE and re-emits a valid bgzipped
  stream to stdout (`pkg/htsgo/tabix.Reheader`). A trailing newline is
  appended to the replacement header when absent so the first data line is
  never merged. Note: upstream htslib's own `tabix -r` is **broken at the
  vendored commit** — it segfaults without `--threads N` and corrupts the
  data body for typical single-/multi-block inputs — so it cannot serve as a
  byte-for-byte producer oracle. Parity is instead validated by using
  upstream tabix as a *consumer*: it re-indexes and queries the Go output
  correctly and sees the replacement header verbatim
  (`TestTabix_ReheaderUpstreamParity`).
- **`--targets` strictness** (`-T`) — now a true overlap post-filter
  (`pkg/htsgo/tabix.Targets`), distinct from `-R`'s index-jump: only records
  overlapping the target intervals are emitted. Coordinate conventions match
  htslib regidx (BED filename → 0-based half-open; otherwise 1-based
  inclusive; chromosome-only lines select the whole chromosome). When `-T`
  is given without an explicit region the whole file is streamed, mirroring
  upstream's `.` region.

**Validation:** live upstream-binary parity tests
(`TestTabix_TargetsStrictUpstreamParity` byte-for-byte across tab/BED,
boundary, adjacent, chrom-only, and region-combined cases;
`TestTabix_ReheaderUpstreamParity` via the consumer round-trip above) plus
strong in-package unit tests. No full upstream-test-suite run yet.

### `samtools`

**Status:** 24 of ~25 subcommands (~96%). `view`, `sort`, `index`, `depth`,
`fastq`, `flagstat`, **`mpileup`** (wave-1 + tail wiring), PR #88's
wave-1 tail (`merge`, `coverage`, `idxstats`, `cat`, `reheader`,
`addreplacerg`, `fixmate`, `dict`, `split`, `quickcheck`), the
heavy-hitter pair `markdup` + `stats`, the calmd/import pair
(**`calmd`** + **`import`**), the niche pair landed in the
phase/targetcut PR (**`phase`** + **`targetcut`**), and
**`consensus`** (simple- and bayesian-mode FASTA/FASTQ/pileup; the
Gap5 posterior caller and the NM-halo MAPQ adjustment are byte-faithful
to upstream's default `MODE_RECALL`; `--het-only` and `--ignore-overlaps`
landed in the #220–#225 wave).

**mpileup `-g/-u` BCF/genotype-likelihood output is now DONE** (it delegates
to the ported bcftools mpileup engine — slices 1–4 plus both indel callers).
The single genuine remaining samtools gap is **`consensus` pileup `-a`
placeholder rows** (the `N\t0\t*\t*` deletion-only / zero-coverage emission;
the per-position calls themselves are already correct). Everything else is
either done or the cross-cutting multi-threading (`-@`) input-decode deferral.

Missing subcommands (in rough priority order):

- **`tview`** — alignment viewer, **all three display modes
  implemented**. The non-interactive **text (`-d T`) and HTML (`-d H`)**
  modes are verified byte-for-byte against the vendored upstream binary
  (`tools/samtools/pkg/samtools/tview.go`, `TestTviewLiveParity`):
  the ruler / reference / consensus (bam2bcf `bcf_call_glfgen` +
  errmod) / packed read rows (the `bam_lpileup.c` level-pool greedy
  packing), insertion-column expansion, the `-d`/`-p`/`-w`/`-s`/`-T`/`-i`
  flags, and the strand/match/mismatch/deletion/refskip characters all
  match. The interactive **`-d C` viewer** (and the bare default on a
  TTY) is a **pure-Go raw-mode loop, no ncurses**
  (`tview_interactive.go` + `tview_tty_linux.go`): it puts the terminal
  in cbreak/no-echo via the `TCGETS`/`TCSETS` termios ioctls, sizes the
  window via `TIOCGWINSZ`, reuses the same frame renderer, and handles
  the upstream `bam_tview_curses.c` key bindings (arrow/`hjkl` scroll,
  paging, `0`/Home, `g` goto, `.`/`i`/`r` toggles, `m`/`b`/`n`/`N`
  colour modes, `?` help, `q`/Esc quit). The key→action and
  action→state logic is unit-tested without a TTY; the termios core is
  the small untestable part. Piped `-d C` exits with a message pointing
  at `-d T` / `-d H`; non-Linux `-d C` reports it requires Linux.
- **`view` flag-tail**: `-X`/`--customized-index` (explicit index-file
  argument after `<in.bam>`) is implemented — the index kind (.bai or
  .csi) is auto-detected from the file's magic. `-L bed` landed as a
  linear-scan BED-region filter; `-M`/`--use-multi-region-iterator` is
  accepted but treated as a no-op since we always run the full
  intersection. `-d/-D` (tag-value filter) and `-N` (qname file) landed
  in the view-d-D-N PR.
- **`mpileup` text BAQ (`-B`/`-E`) — DONE.** The text-pileup path now
  applies BAQ realignment by default whenever a reference (`-f`) is
  supplied, matching upstream `bam_plcmd.c:442`
  (`sam_prob_realn(b, ref, ref_len, (MPLP_REDO_BAQ) ? 7 : 3)`):
  `applyTextMpileupBAQ` in `mpileup.go` fetches each contig once and runs
  `baq.SamProbRealn` in apply+extend mode on every bucketed read, lowering
  the per-base qualities in place before they feed the quality column and
  the `-Q` depth filter. `-B/--no-BAQ` disables it; `-E/--redo-BAQ` adds
  `baq.FlagRedo` (flag 7), recomputing BAQ and ignoring any pre-existing
  `BQ` tag. Previously `-B` was a silent no-op and `-E` was rejected.
  Byte-for-byte parity with `samtools mpileup` is confirmed for the
  default, `-B`, and `-E` modes (`TestMpileup_BAQ_LowersQualities` plus a
  live cross-check). (The `bcftools mpileup` genotype-likelihood path
  applied BAQ already — slice 3 below.)
- **`mpileup` BCF / genotype-likelihood output (`-g/-u`) — DONE.** `-aa`
  zero-fill of empty contigs is implemented (see
  `TestMpileup_AA_ZeroFillTableDriven` and the live byte-for-byte
  `TestParity_Mpileup_T12_AllPositionsZeroFill`, which diffs our output
  against the upstream `samtools mpileup -aa` binary on a small two-contig
  fixture, with and without a reference); the text-pileup path is complete.
  `-g` (BGZF-compressed BCF) and `-u` (uncompressed BCF) now emit the
  per-site genotype-likelihood records (`FORMAT/PL`, the `<*>` unseen
  allele, `INFO/DP/I16/QS/MQ0F`, `FORMAT/AD`).

  Upstream note: modern `samtools mpileup` (the vendored
  `reference_code/samtools` is 1.23.1) **removed** BCF/VCF output entirely
  (the `-g`/`-u` short options are no longer in the getopt string;
  upstream prints "using samtools mpileup to generate BCF or VCF files has
  been removed … please use bcftools mpileup instead"). `bcftools mpileup`
  is upstream's sanctioned replacement and is itself a thin driver over the
  same htslib bam2bcf genotype-likelihood pipeline. Rather than re-port
  bam2bcf a second time, `samtools mpileup -g/-u` **delegates** to the
  already-ported, golden-validated bcftools mpileup engine
  (`tools/bcftools/pkg/bcftools`: `errmod` MAQ model →
  `bcf_call_glfgen`/`combine`/`2bcf` → `pkg/htsgo/bcf` emit). The wiring
  lives in `tools/samtools/pkg/samtools/mpileup_bcf.go` (`MpileupBCF` /
  `MpileupBCFOptions`) and the cmd `-g/-u` branch in
  `cmd/samtools/main.go`. `-Q`/`-d` use the genotype-likelihood-caller
  defaults (min-BQ 1, max-depth 250) unless explicitly set; `-l` maps to
  bcftools `--targets-file`.

  Validation: `TestMpileupBCF_*` (in `mpileup_bcf_test.go`) build a
  FASTA+SAM fixture, run `MpileupBCF`, decode the emitted BCF back to VCF
  text through the in-tree BCF reader (`bcftools.ViewFile`), and assert the
  PL/`<*>`/`AD` records; `-g` and `-u` decode to identical records. The
  cmd-level `TestRunMpileupBCFWiring` drives the real CLI args. A live
  upstream `samtools mpileup -g` comparison is not possible (the path was
  removed upstream); the delegated bcftools mpileup engine is itself
  validated byte-for-byte against `bcftools mpileup` goldens
  (`mpileup.11.out`, `mpileup.12.out`, …) — see the bcftools section below.

  Remainder: the bcftools mpileup port's deferred indel caller
  (`bam2bcf_indel.c`) is inherited here — the SNP genotype-likelihood path
  is complete, full indel-row calling is the shared deferred item tracked
  in the bcftools section (slice 4 / `bam2bcf_indel`).

Plus:

- **CRAM** read/write throughout — DONE. Landed across PRs #162–#180
  (the rANS 4x8/4x16 codecs are in-tree pure Go; `ulikunitz/xz` is the
  only sanctioned third-party dep, confined to the LZMA block codec).
  A later audit + closure pass (C-EmbedRef) added embedded-reference
  decode, stripped the internal `cF` tag, matched htslib's RG/PG aux
  ordering, fixed read-feature → CIGAR ordering for deletions, and
  thereby unblocked **v2.1 decode** for the realistic case — all proven
  byte-for-byte against live `samtools view`. The v2.1 slice-header
  record-counter ITF-8/LTF-8 edge (files with ≥ 2^28 reads before the
  read slice) is now also **closed** (C-V21): `parseSliceHeader` takes
  the container's CRAM major version and reads the counter as ITF-8 for
  v2 / LTF-8 for v3+, matching htslib `cram_decode_slice_header`,
  validated by a 2^28-boundary unit test plus the live-samtools v2.1
  round-trip. **X_EXT bzip2 *encode* is now DONE**: an in-tree
  pure-Go bzip2 encoder (`pkg/htsgo/cram/codec/bzip2_encode.go`, stdlib
  only — Go's `compress/bzip2` is decode-only) backs both the X_EXT
  external codec and a new method-2 (bzip2) block-codec candidate; its
  output decodes byte-for-byte under Go `compress/bzip2`, the system
  `bzip2 -d` (libbz2), and upstream samtools (verified by
  `bzip2_parity_test.go`). Remaining gaps are behind clear errors: the
  network REF_PATH/EBI fetch (an unresolvable reference is a clear MD5
  error) and CRAM v4.0 (spec not final). See `docs/CRAM_DESIGN.md` and
  `docs/CRAM_ROADMAP.md`.
- **`.csi` index** — DONE (PR #189); `samtools index` emits both `.bai`
  and `.csi`, and readers auto-detect index kind from file magic.
- **Multi-threading (`-@`) — DONE for the BAM-writing subcommands
  (`view`, `sort`, `markdup`).** These now drive a parallel BGZF
  compressor when `-@/--threads N` (N > 1) is given, via the new
  `sam.NewBAMWriterThreads` back end built on `bgzf.MultiWriter` (the
  same block-parallel writer used by `bgzip`). BGZF is block-parallel by
  construction — every block is an independent gzip member — so the
  decoded BAM body is byte-identical regardless of thread count. `sort`
  also compresses its temporary external-merge shards in parallel.
  Validated by `tools/samtools/pkg/samtools/threads_test.go`: for
  `view`, `sort`, and `markdup` the decompressed BAM body of `-@ {2..8}`
  is byte-for-byte equal to `-@ 1`, and the decoded record set matches a
  live upstream `samtools` binary built from the vendored submodule
  (`TestThreads_View_UpstreamParity`, `TestThreads_Sort_UpstreamParity`,
  `t.Fatalf` never `t.Skip`). A benchmark (`BenchmarkThreads_ViewBAM`)
  shows the parallel path's throughput gain. Compressed *bytes* may
  differ from upstream's (block boundaries depend on buffering timing);
  the contract is decode-equality, exactly as for the `bgzip` MultiWriter.

  **Parallel *input* BGZF decode** now lands too: `bgzf.MultiReader`
  inflates BGZF blocks across `-@` worker goroutines and reorders them so
  the decoded byte stream is byte-identical for any thread count. It is
  wired into `view`, `flagstat`, `idxstats` (no-index scan), `stats`, and
  `depth` (and mosdepth `-t`), via `alnio.NewReaderThreaded` /
  `OpenReaderThreaded`. Isolated decode throughput roughly doubles at
  `-@ 2` and ~2.3x at `-@ 4` on a 4-core box (see
  `BenchmarkMultiReader`); record-bound subcommands gain less when parsing,
  not inflate, dominates.

  **Remaining single-threaded under `-@`** (the flag is accepted but
  currently a no-op for these): `markdup` pass-1 scan, CRAM *input* slice
  decode (CRAM uses its own container framing, not BGZF blocks — a parallel
  CRAM slice reader is future work) and CRAM encode, and the other
  non-BAM-writing subcommands not listed above (`mpileup`, `fastq`,
  `index`, `merge`, `cat`, `fixmate`, `reheader`, `addreplacerg`, `split`,
  `calmd`, `consensus`, `coverage`, `phase`, `targetcut`). The dominant
  IO-bound wins — parallel BGZF compression of BAM output and parallel BGZF
  decompression of BAM input — are covered.

**Genuine remaining samtools gaps** (everything else is done):

- **`mpileup` MAQ genotype-likelihood model — slices 1-4 DONE; only
  indel calling remains deferred.** The `mpileup` SNP MAQ-model port is
  complete: it was sliced into four parts:
  - **Slice 1 (DONE).** The MAQ error model (`errmod.c`) is ported to
    pure Go in the shared `pkg/htsgo/errmod` package (`errmod.Init` /
    `errmod.Cal`). Both `bcftools mpileup` and `samtools targetcut`
    consume this single implementation.
  - **Slice 2 (DONE).** The per-site genotype-likelihood pipeline —
    `bcf_call_glfgen` / `bcf_call_combine` / `bcf_call2bcf` from
    `bam2bcf.c` — is ported in
    `tools/bcftools/pkg/bcftools/bam2bcf.go`. `mpileup` now emits one
    BCF/VCF record per covered reference position with real MAQ PLs,
    the `<*>` "unseen" allele, the multi-allelic PL grid, and
    INFO/DP/I16/QS/MQ0F. `-O b` (BCF output) works through
    `pkg/htsgo/bcf`; the old `-O u/b` hard-rejection is gone. The
    upstream mpileup defaults are now applied (`min-BQ=1`,
    `max-BQ=60`, `delta-BQ=30`). The `delta_baseQ` neighbour-quality
    cap is implemented.
  - **Slice 3 (DONE).** BAQ realignment is wired into the pileup.
    `applyMpileupBAQ` in `mpileup.go` ports `mpileup.c`'s `mplp_realn`:
    every covered column is gated by the `MPLP_REALN_PARTIAL`
    heuristic, and each selected read is run once through
    `pkg/htsgo/baq.SamProbRealn` in apply+extend mode (`flag 3`,
    matching `sam_prob_realn(b, ref, ref_len, (flag & MPLP_REDO_BAQ) ?
    7 : 3)`) before its bases enter the pileup. `-B/--no-BAQ` disables
    BAQ; `-E/--redo-BAQ` adds `baq.FlagRedo` (`flag 7`). The default is
    PARTIAL realignment — upstream sets `MPLP_REALN | MPLP_REALN_PARTIAL`
    (mpileup.c:1389), so the per-column skip heuristic and the per-read
    spanning check apply. `-D/--full-BAQ` clears `MPLP_REALN_PARTIAL`
    (mpileup.c:1567), forcing full BAQ: every read on the chromosome is
    realigned. The `MPLP_REALN_PARTIAL` per-column skip heuristic depends
    on the per-column `p->indel` term, which needs indel detection — see
    slice 4 — so it is approximated from the `PLP_HAS_INDEL` CIGAR scan;
    exact for indel-free inputs. BAQ runs
    before `accumulateMpileupBases`, so the `delta_baseQ` cap and the
    `min_baseQ` filter inside `bcfCallGlfgen` see the BAQ-adjusted
    qualities — matching `mpileup.c`'s ordering.
  - **Slice 4 (DONE).** The per-site bias annotations are ported in
    `bam2bcf.go`: `calcVDB` (VDB, with the htslib `kf_erfc` rational
    approximation ported as `kfErfc` and the single-precision `float`
    arithmetic upstream relies on), `calcSegBias` (SGB) and
    `calcMWUBiasZ` (the standard-deviation-normalised Mann-Whitney U
    z-score that yields RPBZ / MQBZ / BQBZ / MQSBZ / SCBZ), plus the
    MQ0F fraction. The bias tallies (`ref_pos`/`alt_pos`,
    `ref_mq`/`alt_mq`, `ref_bq`/`alt_bq`, `fwd_mqs`/`rev_mqs`,
    `ref_scl`/`alt_scl`) accumulate per-sample in `bcfCallret` from the
    I16 path and from `getPosition` (the soft-clip-aware read-position
    port), then combine in `bcfCallCombine`. INFO floats — QS, I16 and
    every bias tag — are now rounded through `float32` and rendered with
    a 6-significant-digit `%g` (`formatFloat32G`), matching upstream's C
    `float` storage byte-for-byte. `MPLP_SMART_OVERLAPS` is ported too:
    `applySmartOverlaps` pairs the two mates of each proper pair and
    `tweakOverlapQuality` (a faithful port of htslib's
    `tweak_overlap_quality`, with the `cigar_iref2iseq` iterator and the
    Wang/X31 read-name hashes) merges the overlapping-span base
    qualities before BAQ. Byte-for-byte parity is verified by
    `TestMpileupSNPGoldens`: the full `mpileup/mpileup.11.out` golden
    (4001 covered positions, 87 SNP ALT records, two overlapping mate
    pairs) and the three-sample multi-BAM `mpileup/mpileup.1.out` golden
    match exactly — header and every SNP data record, all bias INFO tags
    included.

  **Remaining deferred work — indel calling only.** The one upstream
  `mpileup` path still unported is the indel caller (`bam2bcf_indel.c` /
  `bam2bcf_edlib.c`): indel candidate detection, the indel genotype
  likelihoods and the INDEL/IDV/IMF INFO tags. The single INDEL record
  of `mpileup.11.out` (17:302 `TA`) is consequently the only golden line
  not reproduced; `TestMpileupSNPGoldens` aligns records by `CHROM:POS`
  and skips it, and `TestMpileupGoldensDeferred` catalogues every other
  deferred golden with its precise reason (FORMAT tags beyond PL, `--ff`
  FLAG filtering, `-s/-S/-G` sample/read-group selection, IUPAC REF
  bases, indel/SCR fixtures).

  **Indel-caller sub-slicing (DONE — legacy caller at parity).** The
  legacy (non-`--indels-cns`) indel caller is fully ported: indel
  candidate generation, per-sample indel genotype likelihoods, and the
  indel-row BCF emission (`ALT`, `INFO/INDEL`, `INFO/IDV`, `INFO/IMF`,
  `FORMAT/PL`/`AD`). The work was broken into five sub-slices, all
  complete. A live-upstream parity harness now backs this:
  `mpileup_indel_parity_test.go` builds the upstream `bcftools` binary
  from the `reference_code/bcftools` submodule (once, via `sync.Once` in
  `upstreamBcftoolsMpileupIndel`) and diffs its `-Ov` output against the
  Go port — `TestMpileupIndelParity_Insertion_Live` asserts the full
  insertion INDEL record byte-for-byte on `indel-AD.2` (G→GTAAA… at
  11:75), and `TestMpileupIndelParity_Deletion_Live` asserts the
  candidate-generation outcome plus `INFO/IDV,IMF,DP` field-for-field on
  `indel-AD.1` (which carries three deletions and one insertion). Both
  `t.Fatalf` on any deterministic mismatch. The only remaining unported
  path is `--indels-cns` (the edlib consensus realigner, a separate
  algorithm — see the explicit note in sub-slice 4e.7 and below).

  **`--indels-cns` (DONE — salvaged from PR #219).** Upstream's
  consensus indel caller (`bam2bcf_edlib.c` / `bam2bcf_iaux.c`, reached
  via `mpileup --indels-cns`) realigns reads against an edlib-built
  consensus haplotype rather than the legacy
  probaln-against-each-candidate path. It is now ported in-tree:
  `bam2bcf_indelcns.go` implements the consensus indel caller and the
  glocal alignment scoring on top of a pure-Go edlib port
  (`pkg/htsgo/edlib/`, Myers bit-vector + NW). The `--indels-cns` flag
  is wired through `subcmds_mpileup.go` / `mpileup.go`. Live-oracle
  parity is asserted byte-for-byte against the upstream binary in
  `TestLiveMpileupIndelsCNS` on the `indel-AD.1` fixture (header +
  every record identical modulo provenance). The `--indels-2.0` /
  `--no-indels-cns` variants remain accepted-and-(mostly)-ignored. The
  legacy-caller homopolymer residuals catalogued below (single-ULP
  probaln drift on a handful of deep reads) are independent of
  `--indels-cns`; they have now been **closed** by matching htslib's
  `float`-width `g_qual2prob` table in `pkg/htsgo/baq` (see cluster 3
  in slice 4e below for the root cause and validation).

  - **4a + 4b (DONE).** Pileup data model + STR finder + indel
    candidate-type discovery helpers. `pileupBase` now carries the
    htslib `bam_pileup1_t.indel` field, an `aux` scratch word, and an
    optional `*sam.Record` back-pointer set only on indel-bearing
    columns; `accumulateMpileupBases` populates them by peeking at the
    next consuming CIGAR op. `bam2bcf_indel.go` adds `bcfCallauxIndel`
    (the indel-specific subset of upstream `bcf_callaux_t`) plus the
    static helpers `est_seqQ`, `est_indelreg`, `bcf_cgp_l_run`,
    `bcf_cgp_find_types`, `tpos2qpos`, and `get_pos`. The STR finder
    (upstream `str_finder.c` / `find_STR` / `find_STR64`) is ported as
    a reusable in-tree package `pkg/htsgo/strfinder` and unit-tested
    over hand-traced inputs (homopolymers, dinucleotide repeats,
    padding, lower-case filter). All CLI knobs (`--open-prob`,
    `--ext-prob`, `--tandem-qual`, `--min-ireads`, `--gap-frac`,
    `--indel-bias`, `--indel-size`) reach `bcfCallauxIndel` via
    `newBcfCallauxIndel`, but the value is not yet driven into
    emission.
  - **4c (DONE).** Per-sample consensus + reference-sample
    construction (`bcf_cgp_ref_sample`, `bcf_cgp_calc_cons`) and the
    Probaln-based alignment scoring core (`bcf_cgp_align_score`) are
    ported in `tools/bcftools/pkg/bcftools/bam2bcf_indel_align.go`.
    `bcfCgpRefSample` samples the per-sample 4-bit IUPAC reference and
    masks positions where ALT ≥30% with N (15); `bcfCgpCalcCons`
    majority-rule-builds the per-type insertion consensus (zeroing
    types whose consensus contains an N); `bcfCgpAlignScore` runs
    `baq.ProbalnGlocal` (BW = |typeLen|+3; PacBio CCS params for
    >1000-bp reads), applies the indel-bias clamp to the
    length-normalised score, and folds in the STR-finder fudge
    (`strfinder.FindSTR`). The returned word reproduces upstream's
    `(sc<<8) | min(255, l)` bit-pattern byte-for-byte, with the low
    byte rewritten as `min(255, (score&0xff)*0.8 + iscore*2)` after
    the STR fudge — verified by `TestBcfCgpAlignScore_BitPattern`.
  - **4d (DONE).** `bcf_call_gap_prep` + `bcf_cgp_compute_indelQ` are
    ported as `bcfCallGapPrep` and `bcfCgpComputeIndelQ` in the same
    file. `bcfCallGapPrep` orchestrates the full pipeline: cheap-reject
    on a clean column → `bcfCgpFindTypes` → `bcfCgpRefSample` →
    `bcfCgpCalcCons` → per-(read,type) `bcfCgpAlignScore`, scoring into
    a flat `N*nTypes` matrix. `bcfCgpComputeIndelQ` then folds those
    scores into per-read `p.aux` words (chosen-type<<16 | seqQ<<8 |
    indelQ), populates `bca.IndelTypes` / `bca.Inscns` / `bca.MaxIns`
    by sumq-sorting the candidate types (REF always at slot 0), and
    returns n_alt. Unit-tested at the bit-pattern level
    (`TestBcfCgpComputeIndelQ_BitPattern`) with a hand-derived
    expectation; the orchestrator has a clean-site -1 reject test and
    an indel-site smoke test.
  - **4e (DONE).** Indel-aware branches of `bcf_call_glfgen` /
    `bcf_call_combine` / `bcf_call2bcf` and emission of the
    INDEL/IDV/IMF INFO tags; golden and live-upstream tests land here.
    Broken into cross-cutting sub-slices since the targeted goldens
    require more than just the indel branch of glfgen+2bcf:
    - **4e.1 (DONE).** FORMAT/AD support and `-a/--annotate` parsing.
      `parseFormatFlag` (`mpileup.go`) is the Go port of
      `parse_format_flag` (mpileup.c:1141) with bit constants
      `B2BFmt*/B2BInfo*` matching `bam2bcf.h:46-75`.
      `validateMpileupOptions` seeds `opts.FmtFlag` from
      `DefaultMpileupFmtFlag` (mpileup.c:1399) and layers user tokens on
      top — `-AD` clears the bit, `AD`/`FORMAT/AD` sets it. `bcfCall2bcf`
      emits FORMAT/AD,ADF,ADR per-sample columns when the matching bit
      is set, with the per-allele depths drawn from the already-reordered
      `call.adf/adr` arrays. The corresponding `##FORMAT=<ID=AD,...>`
      header lines are emitted only when their bits are on.
      `TestMpileupSNPGoldens/multi-bam-region-format-AD` is the
      byte-for-byte golden check against `mpileup/mpileup.12.out`.
    - **4e.2 (DONE).** Indel-branch `bcf_call_glfgen`
      (`bam2bcf.c:300-460`) is now `bcfCallGlfgenIndel` /
      `bcfCallGlfgenCore` in `bam2bcf.go`: it consumes per-read `p.aux`
      words populated by `bcfCallGapPrep`, runs `errmodCal` on the indel
      bases, and populates `bcfCallret` with indel-specific QS / AD /
      I16 / bias tallies (using `is_diff = b ? 1 : 0`, mirroring
      `bam2bcf.c:350`). `bcfCallGapPrep` was extended to populate
      `bca.IrefPos / IaltPos / IrefMq / IaltMq / IrefScl / IaltScl`
      (`bam2bcf_indel.c:826-848`) from the per-read loop at t==0; the
      indel-flavored `getPos` helper supplies the read-position and
      soft-clip-length bins. `pileupBase.rec` is now populated for
      every read in the pile (not just indel-bearing columns) so the
      indel iref/ialt accumulation has the cigar available for all
      reads.
    - **4e.3 (DONE).** Indel-branch `bcf_call_combine` + `bcf_call2bcf`
      (`bam2bcf.c:1165-1198`, `bam2bcf.c:1211-1234`) are
      `bcfCallCombineIndel` (in `bam2bcf.go`) and `bcfCall2bcfIndel`
      (in `mpileup.go`): REF/ALT alleles built from
      `bca.Inscns`/`IndelTypes`/`IndelReg`, `INFO/INDEL` flag plus
      `IDV`/`IMF` emitted before `DP`/`I16`/`QS`, then the same bias
      subset as the SNP path (`VDB`/`SGB`/`RPBZ`/`MQBZ`/`MQSBZ`/`BQBZ`/
      `SCBZ`/`MQ0F`) followed by `FORMAT/PL`. The upstream "leaked"
      bias semantics (BQBZ/MQSBZ retain the last has-alt SNP's value
      since `bcf_callaux_clean` does not reset `call->mwu_*`) is modelled
      explicitly by a `biasLeak` struct threaded through the per-site
      driver. The driver in `emitChromMpileup` runs `bcfCallGapPrep`
      after the SNP emission and, when it returns ≥0, runs a second
      `bcfCallGlfgenIndel` + `bcfCallCombineIndel` + `bcfCall2bcfIndel`
      pass (matching `mpileup.c:589-613`). The full
      `mpileup/mpileup.11.out` golden — including the 17:302 INDEL row
      (`T → TA`, `BQBZ=-1.34164` inherited from the prior SNP combine
      at 17:237) — now matches byte-for-byte.
    - **4e.4 (DONE).** `--ambig-reads` (incAD / incAD0) ADF/ADR
      compensation (`bam2bcf.c:540-561`), `--skip-{all,any}-{set,unset}`
      BAM-flag filters (`mpileup.c:208-211`), and the three latent items
      the 4e.2+4e.3 review flagged:
        1. `p.is_del` modelling — `accumulateMpileupBases` now emits a
           pileupBase event for every reference column inside a read's
           `D` op (`isDel=true`, `b=0`); `bcfCallGlfgenCore` lets these
           reads through in the indel branch (matching upstream
           `bam2bcf.c:307`if (p->is_del && !is_indel) continue`).
        2. `p.is_refskip` modelling — same shape for CREF_SKIP (`N`)
           ops; `isRefskip` reads are dropped in both branches
           (`bam2bcf.c:301`).
        3. Cross-contig `biasLeak` reset — the `biasLeak` instance now
           lives on the per-run driver in `writeMpileupVCF` and threads
           through every `emitChromMpileup` call so the BQBZ / MQSBZ
           scalars persist across contigs, mirroring upstream's `conf->bc`
           lifetime. The leak is pre-initialised to `(0, ok=true)` so the
           very first indel record (before any has-alt SNP combine) sees
           BQBZ=MQSBZ=0 — matching upstream's C `bcf_call_t bc = {0}`
           default-initialisation.
      `--ambig-reads` is parsed via `parseAmbigReads` into
      `AmbigReadsMode` on `MpileupOptions`; the indel-branch glfgen
      stashes low-quality REF-looking reads in `adrRefMissed[]` /
      `adfRefMissed[]` per upstream and applies the `incAD` (proportional)
      or `incAD0` (claim as REF) compensation. The `--skip-*` strings go
      through `parseBAMFlagString` (Go port of htslib's `bam_str2flag`)
      into the `RflagSkip{Any,All}{Set,Unset}` masks consumed by
      `mpileupKeepRecord`. Byte-for-byte golden: `mpileup-filter.2.out`
      now matches (both the `--skip-all-unset READ1` and `--skip-any-unset
      READ1` forms; tested in `TestMpileupFilterGolden`).
    - **4e.5 (DONE).** INFO/FMT/SCR (soft-clipped reads counter). The
      SCR signal is a per-record `hasSoftClip` bit (true iff the read
      has any CIGAR S op, mirroring upstream's `PLP_HAS_SOFT_CLIP` set
      in `pileup_constructor`, `mpileup.c:317-323`). The bit is stamped
      on every `pileupBase` produced by `accumulateMpileupBases` so the
      SNP-branch `bcfCallGlfgenCore` can tally it pre-refskip into
      `bcfCallret.scr` (matching `bam2bcf.c:300`). `bcfCallCombine`
      folds the per-sample counts into `bcfCall.scrTotal` /
      `bcfCall.scr[]`, and `bcfCall2bcf` emits INFO/SCR (before I16) /
      FORMAT/SCR (after AD/ADF/ADR) when their bits are set.
      `parseFormatFlag` now also accepts the `FMT/` prefix in addition
      to `FORMAT/` (`SET_FMT_FLAG`, `mpileup.c:1120-1122`). Byte-for-
      byte golden `mpileup-SCR.out` (test.pl:1069) now matches; tested
      in `TestMpileupSCRGolden`.
    - **4e.6 (DONE).** INFO/NMBZ (per-read NM bias) — the Mann-Whitney
      U z-score over the per-read NM tag, split REF vs ALT. The port
      adds `getAuxNm` (Go counterpart of `get_aux_nm`, `bam2bcf.c:96`)
      which reads the BAM `NM:i:` tag, treats each indel CIGAR op as a
      single event (subtracting `len-1` for ops with `len>1`), counts
      soft-clip lengths as mismatches, then subtracts 1 for REF reads
      and 2 for ALT reads (the MNP-aware adjustment) and clamps to
      `[0, b2bNNm-1]`. Per-sample `refNm[32] / altNm[32]` histograms on
      `bcfCallret` are filled by both the SNP and indel branches of
      `bcfCallGlfgenCore` when any of the `B2BInfoNMBZ` / `B2BFmtNMBZ`
      / `B2BInfoNM` bits is set. `bcfCallCombine` sums them across
      samples and runs `calcMWUBiasZ`; `bcfCallCombineIndel` also
      folds in the matching SNP-pass tallies (mirroring upstream's
      shared `bca->ref_nm/alt_nm` accumulator). `bcfCall2bcf` and
      `bcfCall2bcfIndel` emit `INFO/NMBZ` between MQSBZ and SCBZ when
      `B2BInfoNMBZ` is set, and the header line is inserted in the
      same place. Byte-for-byte golden `annot-NMBZ.1.1.out` (test.pl
      line 1074) now matches; tested in `TestMpileupNMBZGolden`.
      `annot-NMBZ.2.1.out` byte-matches once the depth-cap port is in
      place (4e.8 below); `.3.1.out`'s SNP row byte-matches including
      `NMBZ=7.74597` while the indel row's QS / NMBZ / PL[0] residual
      tracks with the `bcfCall2bcfIndel` SCR-on-indel-rows polish item.
      FORMAT/NMBZ (per-sample) is not in this slice; only INFO/NMBZ.
    - **4e.7 (DONE).** Ports the legacy REF-rescue heuristic that
      `bcf_call_glfgen` applies in its indel branch
      (`bam2bcf.c:338-348`, originally Heng Li's e4e161068 fix for
      htslib issue #1446). At deeply-covered homopolymer / tandem-
      repeat sites `bcf_cgp_compute_indelQ` emits `indelQ = 0` for
      most REF-leaning reads — so without the rescue every REF read
      fails the indel-branch min-baseQ gate and lands in
      `ADR/ADF_ref_missed`, breaking I16 / QS / AD parity. The
      heuristic: when a read has no CIGAR indel (`p.indel == 0`) and
      either `q < _n/2` or `_n > 20`, reclassify it as REF (`b = 0`),
      promote `q` to the read's raw base quality at qpos, and rebuild
      `seqQ` as `(3*seqQ + 2*q)/8`; once `_n > 20`, cap that seqQ at
      40. `baseQ` for I16 stays at the pre-heuristic value in
      `p.aux>>8&0xff` (matching `bam2bcf.c:420`); the local seqQ is
      used only to cap `q` and gate min-baseQ. Two correctness fixes
      landed at the same time: (a) `seqQ = q = (p.aux & 0xff)` on the
      indel branch (upstream `bam2bcf.c:315` initialises both from
      the indelQ bits, not the saved seqQ bits), and (b) the SNP
      branch now stores `seqQ = b2bSeqQ` explicitly so the shared
      `if (q > seqQ) q = seqQ` step at `bam2bcf.c:459` is the same
      literal port. Goldens: `indel-AD.{2,3,4}.out` are byte-for-byte
      (tested in `TestMpileupIndelADGolden`, covering `--ambig-reads`
      default / `incAD` / `incAD0`); the `annot-NMBZ.3.1.out` indel
      row's I16 also matches byte-for-byte. The `--indels-cns` (edlib)
      path is a separate algorithm and stays deferred. Residual
      divergences live on `indel-AD.1.out` (~20 SNP rows with small
      I16 base-quality drifts plus four homopolymer-column indel rows
      with chosen-type off-by-1 assignments — DP ≤ 125 throughout so
      the depth cap is not in play) and the trailing QS/NMBZ/PL[0]
      columns of `annot-NMBZ.3.1.out`'s indel row.
    - **4e.8 (DONE).** Ports htslib's per-alignment-start depth cap.
      Upstream's `bam_plp_push` (reference\_code/htslib/sam.c:6090)
      drops a new read when `iter->pos == b->core.pos` and the
      pileup queue already holds `maxcnt` active reads, while our
      previous port truncated each per-column pile to `MaxDepth`
      reads instead. The new `applyMpileupDepthCap` walks the per-
      sample coordinate-sorted record stream, maintains a min-heap
      of in-flight end positions, and drops reads using exactly the
      htslib predicate (including the per-alignment-start trigger
      and the "mp->cnt includes one sentinel node" off-by-one). The
      cap runs BEFORE `applySmartOverlaps` so a capped-out read
      cannot leak a tweaked base quality into the surviving mate —
      upstream's `overlap_push` runs inside `bam_plp_push` after the
      cap test, so dropped reads never reach the overlap-quality
      merger. Byte-for-byte golden `annot-NMBZ.2.1.out` at chr6:75
      (raw coverage 449, capped to DP=283) now matches; tested in
      `TestMpileupDepthCapGolden`. Remaining mpileup residuals are
      the `--indels-cns` edlib path (separately deferred) and the
      indel-row QS/NMBZ/PL[0] columns at homopolymer columns on
      `indel-AD.1.out` and `annot-NMBZ.3.1.out`, both tracked under
      the `bcfCall2bcfIndel` SCR-on-indel-rows polish.
    - **4e.5 indel-row SCR (DONE).** `bcfCall2bcfIndel` now emits
      `INFO/SCR` (before I16) and `FORMAT/SCR` (after AD/ADF/ADR)
      when the corresponding `B2BInfoSCR` / `B2BFmtSCR` bits are
      set, mirroring the SNP-row code path. `bcfCallCombineIndel`
      copies the per-sample SCR tally from the SNP-pass
      `bcfCallret.scr` arrays (the indel branch of `bcfCallGlfgen`
      does not tally SCR — the SCR accumulator is gated on
      `!isIndel`, matching upstream's `bam2bcf.c:300`). Both the
      SNP and indel rows at the same column therefore report the
      same SCR counts. Regression: `TestMpileupSCROnIndelRow` (uses
      the `indel-AD.2.fa` / `indel-AD.2.bam` fixture, which has a
      homopolymer-anchored indel call and soft-clipped reads at
      `11:75`).
    - **Residual: `indel-AD.1.out` and `annot-NMBZ.3.1.out` indel-row
      drifts at homopolymer columns (DOCUMENTED, root cause traced).**
      A column-by-column diff against the goldens identifies three
      independent clusters:

      1. **Trailing N-REF rows past the FASTA end** (`indel-AD.1.out`
         positions `000000F:687-688`, two records with `REF=N`,
         `DP=1`, all I16 fields zero). Root cause: the 000000F
         contig is 686 bp in the FASTA, but the BAM contains a read
         whose CIGAR (`6M1D117M5D28M`) ends at reference position
         688 — two bases past the FASTA boundary. Upstream's
         pileup engine does not bound itself to the FASTA length:
         it walks reads' CIGARs and emits a column for every
         covered reference position, using `N` for the REF when
         the position is past the FASTA. Our port allocates
         `events[i]` of length `refLen` (the FASTA length) and the
         per-site loop in `emitChromMpileup` terminates at
         `pos0 < refLen` (mpileup.go:1080, 1095). Fix would
         require extending the events array to the maximum
         read-end across inputs, padding the REF with `N` for
         positions past `refLen`, and adjusting the regions/targets
         intersection to allow positions past the FASTA. Scoped as
         a follow-up; the two affected rows carry no biological
         signal (DP=1 with zero base-quality contribution).

      2. **SNP-row I16 base-quality-sum micro-drifts** (`indel-AD.1.out`
         ~12 columns near 446-624, all with I16 slots 4-5 off by a
         single read's base quality, e.g. 1204/45840 vs upstream's
         1205/45921 = one missing BQ=9 contribution). The reads
         involved cover both the BAQ-adjusted homopolymer at
         `000000F:537-540` and the columns straddling the FASTA
         boundary at 686. These are post-BAQ quality drifts on
         reads whose tail bases extend past the FASTA end —
         upstream's BAQ adjustment for those tail bases differs
         from ours, propagating a one-quality-unit shift back into
         the SNP-row I16 sums at columns where the affected reads
         contribute REF bases. Same root cause as cluster (1):
         FASTA-boundary handling.

      3. **Indel-row chosen-type off-by-one at homopolymer columns**
         (`indel-AD.1.out` at `000000F:537/538/658`, plus
         `annot-NMBZ.3.1.out` at `chr16:75`). At these columns the
         indel-row I16 fields agree on REF/ALT classification (so
         the `isDiff = b ? 1 : 0` split is identical), but the
         **per-allele** breakdown shifts because a handful of reads
         are assigned to a different non-REF indel type by
         `bcfCgpComputeIndelQ`. Concretely at `000000F:537`: ours
         classifies one extra read as type `-1` (the deletion),
         shifting I16 slot 3 (alt-rev count) from 27 to 28, the
         alt BQ-sum by exactly one BQ=40 contribution, the alt
         MQ-sum by MQ=60, and the alt min-dist sum by 25. This
         propagates into QS (per-type qsum), AD (per-allele
         counts), and (transitively) PL[2]. At `chr16:75`
         (`annot-NMBZ.3.1.out`) the I16 byte-matches because the
         swaps are between two ALT types (both `b!=0`, so isDiff
         stays 1), but QS shifts (1.45884 vs 1.43466 for type 1),
         NMBZ flips sign (+0.437589 vs -0.886523 — the indel-pass
         refNm/altNm split changes), and PL[0] for sample 1 falls
         from 255 to 226. Root cause: `bcfCgpComputeIndelQ` and
         `bcfCgpAlignScore` use the same encoding `score<<6 | t`
         and the same ascending-insertion-sort tie-break as
         upstream, and the orchestrator's `tbeg`/`tend`/`rStart`/
         `rEnd` / ref2-slice positioning matches upstream
         byte-for-byte (verified by `bcfCgpAlignScore`'s
         existing tests). The remaining divergence was in the
         `probaln_glocal` score itself for reads whose query
         window straddles a homopolymer run — at most-frequent /
         long homopolymers, two indel types could produce
         `probaln_glocal` returns that differ by a single Phred
         unit, and a one-unit drift in the underlying HMM
         likelihoods flipped the tie-break.

         **CLOSED (float-width alignment).** Root cause was a
         `float` vs `double` width mismatch, not FP-operation
         order. htslib declares its per-base error-probability
         table as a C `float` array (`static float
         g_qual2prob[256]`, populated by `pow(10, -i/10.)` stored
         into a `float`) and the per-query `qual[]` buffer is also
         `float`; the forward/backward recurrences only promote
         those values back to `double` when they enter the
         (double-typed) HMM transition expressions. The Go port in
         `pkg/htsgo/baq` held the same table at full `float64`
         precision, so each entry differed from htslib's
         float-rounded value by ~1e-8..1e-11. That tiny difference
         propagated through the scaled forward DP and, at long
         homopolymer columns on deep reads, flipped the integer
         `(int)(Pr1 + .499)` Phred score by one unit. Rounding the
         `qual2prob` table through `float32` (so it equals
         htslib's `float` table bit-for-bit, then promotes to
         `float64` for the arithmetic exactly as C promotes
         `float`→`double`) closes the drift. A standalone fuzz over
         20k homopolymer-straddling alignments found exactly one
         integer-score divergence between the float64 and float32
         tables (141 vs 140); the upstream `probaln.c` self-test
         returns 140, matching the float32 path. That case is
         pinned as `TestProbalnGlocalHomopolymerFloatWidth` in
         `pkg/htsgo/baq/probaln_test.go`. With the fix, the live
         indel parity test now diffs the **full** `indel-AD.1`
         INDEL records (I16, QS, VDB included) byte-for-byte
         against the upstream binary
         (`TestMpileupIndelParity_Deletion_Live`), and the
         `annot-NMBZ.3.1` `chr16:75` column likewise agrees. No
         algorithmic behaviour changed — only the table's storage
         width was matched to htslib.

  One accepted divergence: `errmod_cal`'s downsampling of piles deeper
  than 255 reads uses Go's RNG rather than htslib's `drand48`, so
  byte-for-byte parity holds only at depth ≤255 (RNG byte-parity is not
  a project goal). All vendored `mpileup` fixtures are within that
  bound.

  **Parity watch-item (resolved) — QS-sum zero-break float comparison.**
  `bcfCallCombine` in `bam2bcf.go` breaks the allele-ordering loop on
  `qsum[ipos] == 0`, mirroring upstream `bcf_call_combine`
  (`bam2bcf.c:991-1001`, which tests a C `float`). The slice-4 goldens
  exercise this path heavily (87 multi-allelic SNP sites) and match
  byte-for-byte, so the C-`float`-vs-Go-`float64` width concern did not
  materialise on the upstream fixtures. The note is kept for awareness
  on future high-coverage inputs.
- **`phase` deterministic-DP phasing (DONE — byte parity).**
  Upstream `phase.c` has **no MCMC** (there is no `phase_core`): the
  CLI default path computes an `int8 path` via the `dynaprog`
  Viterbi recurrence and repairs chimeras in `fragphase`, fully
  deterministically. The Go port replicates that DP, `fragphase`,
  `genmask` and the khash/ksort emit order; the full `samtools
  phase` stream (PS/FL/M/EV) is byte-identical to upstream on the
  default LOD. `-b`/`-F`/`-A`/`-l`/`-e` all landed. Detail in the
  `phase` subsection below.
- **`targetcut` BAQ realignment with `-f` reference (DONE).** The
  HMM consensus mode is implemented (faithful port of
  `cut_target.c`, including the MAQ errmod port; see below). The
  per-record BAQ realignment that upstream's `read_aln` applies
  via `sam_prob_realn` when a `-f` reference is supplied is now
  wired: `targetcutHMM` runs `pkg/htsgo/baq.SamProbRealn` in
  apply+extend mode (`flag = 1<<1|1`) on every record that
  survives the read filter, on a per-chromosome cache of the
  reference fetched via `fasta.RandomAccess`. The stderr warning
  is gone; the BAQ-adjusted qualities feed `gencns` exactly as
  upstream feeds them into the pileup.
- **`tview`** — text (`-d T`) and HTML (`-d H`) modes implemented
  (byte-for-byte parity), plus the interactive `-d C` viewer (pure-Go
  Linux raw-mode termios, no ncurses; piped/non-Linux `-d C` exits with
  a clear message).

**`markdup` implemented features** (closed in the per-subcommand
sub-features PR; live-validated byte-for-byte against the upstream binary
in `subfeatures_live_oracle_test.go`):

- **`-d` optical-duplicate detection (DONE).** Each flagged duplicate is
  compared to the chosen original by the colon-delimited Illumina read-name
  tile coordinates (a port of `bam_markdup.c` `get_coordinates_colons` +
  `is_optical_duplicate`). Within `-d N` on both x and y axes the duplicate
  is counted as optical and tagged `dt:Z:SQ`; otherwise `dt:Z:LB`. The
  `--read-coords` regex variant and the `check_chain` whole-cluster
  re-expansion (which can promote a duplicate-of-a-duplicate to optical) are
  not reproduced — read names must use the colon Illumina layout.
- **`-s/--stats` and `-f FILE` (DONE).** The full duplicate-statistics block
  (READ / WRITTEN / EXCLUDED / EXAMINED / PAIRED / SINGLE / DUPLICATE
  PAIR|SINGLE|PAIR OPTICAL|SINGLE OPTICAL|NON PRIMARY|NON PRIMARY OPTICAL /
  DUPLICATE PRIMARY TOTAL / DUPLICATE TOTAL / ESTIMATED_LIBRARY_SIZE) is
  emitted, including the Picard-style library-size bisection solve. The CLI
  prepends the upstream `COMMAND:` line; byte-parity is asserted on the
  counter block (the COMMAND line echoes the verbatim argv and so differs by
  design between the two binaries).
- **`-S` supplementary/secondary duplicate marking (DONE).** By default
  non-primary records are emitted untouched (matching upstream). Under `-S`
  a non-primary record is flagged only when its primary duplicate carried an
  SA/XA aux tag or an unmapped mate, mirroring upstream's `add_duplicate`
  gate, and is counted under DUPLICATE NON PRIMARY. (Earlier the Go port
  over-marked all same-qname non-primary records by default — fixed here.)
- **`dt:Z:` "duplicate-type" aux tag (DONE under `-d`).** Written as SQ
  (optical) or LB (library) for every flagged duplicate when `-d` is active,
  in upstream's `do`-then-`dt` aux order.
- **EXCLUDED / EXAMINED accounting (corrected here).** The eligibility mask
  now matches upstream exactly (`SECONDARY|SUPPLEMENTARY|UNMAP|QCFAIL`): an
  unmapped primary counts as EXCLUDED, and a record that merely already
  carries the duplicate flag is EXAMINED (re-scored) rather than excluded.
  (The library-size "unable to calculate" stderr warnings upstream prints in
  degenerate cases are not reproduced.)

**`markdup` still-deferred features** (flag slots accepted on the CLI for
compat):

- Per-read-group statistics split (`--read-groups`) and barcode
  regex / barcode-tag keying. The single-group counters are byte-exact;
  multi-group output and barcode keying are not split. Fixture
  `reference_code/samtools/test/markdup/17_read_group.sam` remains a
  documented partial-parity skip.
- The `--json` stats variant (the text `-s` form is the validated target).

**`coverage` histogram mode (DONE).** The `-m` (UTF block histogram),
`-A` (ASCII `.`/`:` ramp), `-D` (plot summed depth per bin) and `-w N`
(bin count) modes are a faithful port of `coverage.c` `print_hist`: the
per-bin breadth/depth fill, the 10-row block-character ramp
(`round(8*(val-bin)/rowsize)-1`), the side-panel statistics
(reads/filtered/covered/percent/mean cov/baseQ/mapQ/bin width/max bin),
and the centered K/M/G/T x-axis labels. Without a TTY the bin count
defaults to 40 (upstream's terminal-width fallback). The tabular mode was
also corrected to upstream's exact `print_tabular_line` formatting (header
`meanbaseq`/`meanmapq`, `%g` coverage/meandepth, `%.3g` baseQ/mapQ), and
baseQ is now summed only at covered positions so `--min-depth` matches.
Live-validated in `subfeatures_live_oracle_test.go::TestLiveCoverageHistogram`
across `-m / -A / -D / -w / -Q / -q / --min-depth` and the tabular form.

**`targetcut` multi-region output (DONE — no separate "multi-window"
flag).** Upstream `cut_target.c` (getopt `f:Q:i:o:0:1:2:`) has no
window-size CLI option; the multi-region behaviour is the 2-state Viterbi
HMM emitting one consensus SAM record per identified covered/callable
segment, which the Go port already implements (see the `targetcut` BAQ
section). There is nothing additional to port beyond what landed in the
phase/targetcut PR.

**`markdup -l/--max-len` is a no-op-by-design.** Upstream uses `-l` solely
as the streaming buffer flush window in `bam_markdup.c:1949`; it does NOT
affect key construction or scoring. Our two-pass implementation buffers
per-bucket state in memory, so output is identical for any `-l` value.
The flag is accepted on the CLI and the option is preserved on
`MarkdupOptions.MaxLen` for forward compatibility if we ever move to a
single-pass streaming model.

**`stats` implemented sections.** The per-cycle quality matrices
**FFQ/LFQ**, the GC-content sections **GCF/GCL**, the ACGT-content
sections **GCC/GCT**, the indel sections **IC/ID**, the leading **CHK**
CRC32 checksum block, the **COV** coverage-distribution histogram and
the **GCD** GC-depth distribution are byte-faithful to upstream —
validated against the vendored `.sam` fixtures. CHK sums the per-record
CRC32 of read names, the BAM 4-bit-packed sequence and the quality
bytes; COV bins each reference position's M/=/X read depth via the
`-c MIN,MAX,STEP` option (default `1,1000,1`) and is emitted only for
coordinate-sorted input, matching upstream's `is_sorted` gating. COV
depth is accumulated in a bounded per-contig sliding window that is
flushed as records advance (mirroring upstream's `cov_rbuf` ring
buffer), so COV memory is O(longest read span), not O(genome).

GCD splits the reference into `--GC-depth`-wide segments (default
20000 bases) and records per-segment read depth plus GC content, then
at output sorts segments by GC and reports depth percentiles. Both
upstream code paths are ported: the default no-reference path
approximates GC content from the read sequences, while the
`-r/--ref-seq` reference path reads GC content from the indexed
reference FASTA (`fai_gc_content`). Like COV, GCD is emitted only for
coordinate-sorted input. The upstream `igcd`/`ngcd` indexing quirk —
where `gcd[0]` is an empty placeholder and the final segment is never
finalised before sorting — is replicated exactly for byte-parity.

For unsorted input the COV section is silently omitted, whereas
upstream aborts with `Expected coordinates in ascending order` — a
deliberate, friendlier divergence.

**`stats` reference-free option tails.** BWA-style quality trimming
(`-q/--trim-quality`) and the `-t/--target-regions` restriction are
implemented and validated against the upstream `stat/11` fixtures.
`-q` ports `bwa_trim_read` faithfully, including bwa's documented
off-by-one and the `BWA_MIN_RDLEN` (35) early return, and feeds the
`bases trimmed` SN counter. `-t` parses the upstream target-regions
format (`seq-name beg end`, 1-based inclusive — not BED), merges
overlapping intervals, restricts every counter to reads overlapping a
target interval, clips `bases mapped (cigar)` and the COV depth to the
target, and emits the `bases inside the target` and
`percentage of target genome with coverage > N` SN lines (the latter
threshold set by `-g/--cov-threshold`, default 0). The SN and COV
sections are byte-faithful to `11.stats.expected` /
`11.stats.g4.expected`.

**`stats` reference-statistics sections.** The **MPC**
mismatches-per-cycle section (emitted with `--ref-seq`) and the **RFS**
reference-statistics section (emitted with `--ref-stats`, plus its
`--ref-stats-chunk` companion) are implemented and byte-faithful to the
upstream `stat/` golden files. MPC ports `count_mismatches_per_cycle`
faithfully — including the cycle-index handling for soft/hard clips and
insertions, the reverse-strand mirroring, the N-base bucketing in
quality slot 0, and the documented `uint8` `qual+1` wrap that lands
mismatches of `*`-quality reads in the N column. RFS ports
`collect_refstats`: without `-t` it reports one row per `@SQ` header
entry, with `-t` it reports the merged target intervals as
`name:start-end` rows, and without `--ref-seq` the GC/N columns report
the `-1` lack-of-data sentinel. Validated against `test.pl` stats cases
1-8 (`-r test.fa`, MPC) and 16/17/19 (`--ref-stats`, RFS).

**`stats` per-fragment ACGT sections.** **FBC/LBC** (ACGT content per
cycle for first / last fragments) and **FTC/LTC** (the matching
A/C/G/T/N raw-counter totals) are implemented and byte-faithful to the
vendored `stat/` fixtures. They derive from the same per-fragment
cycle buffers GCC/GCT already accumulate. Note: despite the name these
are NOT barcode tables — earlier roadmap text mislabelled them.

**`stats` per-barcode sections (implemented).** The per-barcode
ACGT-content (`<tag>C`) and quality (`<tag>Q`) tables are a faithful
port of `collect_barcode_stats` (stats.c:773) and its output
(stats.c:1748). Collection is unconditional for the four fixed aux-tag
pairs upstream's `init_barcode_tags` installs — `BC`/`QT`, `CR`/`CY`,
`OX`/`BZ`, `RX`/`QX` — and a section is emitted only once a barcode for
its tag has been observed. The barcode separator, the two segment
columns either side of it, and the per-tag `max_qual` quality-column
count all match upstream. Note: this samtools version has **no**
`--barcodes` CLI option — barcodes are always collected — so no such
flag is wired; `--barcode-tag` / `--quality-tag` tag renaming likewise
does not exist here. Validated byte-for-byte against the exact
`test.pl` invocations `samtools stats 13_barcodes_ok.sam` and
`samtools stats 13_barcodes_ok_ox_bz.sam` (test.pl:3325-3326). The
malformed-barcode `expect_fail` fixtures (test.pl:3327-3329) are
exercised for graceful warn-and-skip behaviour but not byte-compared,
since their output deliberately diverges from any clean baseline.

**`stats` `--sparse` (corrected).** Upstream `-x/--sparse`
(stats.c:2170) only "suppresses outputting IS rows where there are no
insertions" — its sole effect is at stats.c:1796, thinning all-zero IS
rows. Earlier this code wrongly suppressed *every* histogram section.
It now emits every section unconditionally and honours `sparse` per-row
in the IS section only. The IS section itself was reworked: it now
emits the full `0..ibulk-1` row range (including zero rows in the
default mode) where `ibulk` mirrors upstream's last-non-zero /
99%-bulk-truncation logic, instead of the old observed-sizes-only map.
Insert sizes are now classified per-read from each record's own flags
and mate position and halved at output, matching upstream exactly
(this also fixed an inward/outward misclassification the prior
per-pair logic produced on multi-pair inputs).

**`stats` remaining tail** (also documented in `PARITY_VALIDATION.md`):

- Command-line positional region arguments **are now accepted**
  (`samtools stats in.bam chr1:100-200 chr2 ...`). They restrict every
  counter the same way `-t/--target-regions` does — routed through the
  same `isInRegions` overlap filter and driving the same "bases inside the
  target" / "percentage of target genome" SN lines (upstream's
  `sam_itr_regarray` + `replicate_regions` path). The CLI requires an
  indexed input (`.bai`/`.csi`) for positional regions, mirroring upstream's
  "Random alignment retrieval only works for indexed files" error; the
  restriction itself is a streaming overlap filter, so the full report is
  byte-identical to upstream across all sections. `-t` and positional
  regions are mutually exclusive (upstream's `if (!targets)` guard); `-t`
  wins when both are given. Validated live against the upstream binary in
  `stats_regions_upstream_test.go` (all sections compared, header `#`
  comment lines aside), with binary-free `TestUnitRegionsFromSpecs` /
  `TestUnitStatsIsInRegions` coverage of the parser + overlap predicate.
- `--remove-overlaps` is accepted as a no-op; single-record stats are
  unaffected by overlap removal for the counters emitted.

With `--sparse` corrected and the per-barcode sections implemented,
`samtools stats` reaches full 1:1 output parity for every section.

The output emits the byte-faithful **CHK** checksum block, **SN**
(Summary Numbers), the per-cycle and base-content sections
(**FFQ/LFQ/GCF/GCL/GCC/GCT/FBC/FTC/LBC/LTC/IC/ID**), the per-barcode
**`<tag>C`/`<tag>Q`** tables, the **MPC** mismatches-per-cycle matrix,
the **RL / MAPQ / IS** rollups, the **COV** coverage histogram, the
**GCD** GC-depth distribution and the **RFS** reference-statistics
section. `--sparse` thins only all-zero IS rows.

**Validation:** upstream fixtures from `reference_code/samtools/test/markdup/`
and `.../test/stat/` are vendored under
`tools/samtools/testdata/parity/{markdup,stat}/`. The byte-exact /
flag-exact / SN-byte cases are exercised in
`tools/samtools/pkg/samtools/markdup_test.go` and `stats_test.go`.

**`calmd` BAQ realignment** (`-r`, `-E`, `-A`, `-C`) — implemented:

- The probabilistic banded glocal forward-backward HMM
  (`probaln_glocal`) and the CIGAR-aware BAQ driver / MAPQ cap
  (`sam_prob_realn`, `sam_cap_mapq`) are ported faithfully into the
  shared `pkg/htsgo/baq` package — shared so `mpileup` can reuse it.
- **`-r`** computes the `BQ:Z` base-alignment-quality aux tag;
  **`-rE`** is extended-BAQ mode; **`-rA`** applies BAQ to the base
  qualities and writes a `ZQ:Z` tag. **`-C INT`** (threshold `> 10`)
  caps MAPQ via `sam_cap_mapq`.
- Validated byte-for-byte against htslib's own `realn0{1,2,3}*`
  golden fixtures: all 8 `test_realn` flag combinations plus the
  `FlagRedo` path pass exactly. `probaln_glocal` is additionally
  pinned against score / posterior-quality vectors captured from the
  upstream `probaln.c` self-test.

**`calmd -C` / `-e` / `-u` (DONE).** Despite an earlier note that these
were "deferred", the library and CLI implement them:

- **`-C INT`** caps MAPQ via `sam_cap_mapq` (threshold `> 10`), see the BAQ
  paragraph above.
- **`-e`** rewrites matching M-op SEQ bases to `=` (NM/MD still computed
  against the original base).
- **`-u`** emits uncompressed BAM (implies `-b`).

These are live-validated byte-for-byte against the upstream binary in
`subfeatures_live_oracle_test.go::TestLiveCalmdKnobs` (records identical
modulo the `@PG` line the Go port omits by design).

**`calmd` deferred features** (accepted as CLI flags, behaviour partial):

- **`-h` HASH_QNM** (hash-based query-name binarisation) — niche
  upstream-only optimisation; not implemented.
- **`-N` clear-MD/NM-bits**, **`--no-PG`** —
  CLI-accepted-and-ignored stubs.

**`calmd` implemented post-MD/NM transforms** (`bam_md.c` upstream
order — max-NM masking → write NM → write MD → DROP_TAG → BIN_QUAL):

- **`-d` DROP_TAG** — drops every aux tag except `RG`. Applied after
  the NM/MD fill, so the freshly-computed NM/MD are dropped too; only
  `RG` survives (records without `RG` keep no aux at all).
- **`-q` BIN_QUAL** — reduces base-quality resolution: each quality
  `>= 3` maps to `qual/10*10 + 7` (integer division); lower values
  unchanged.
- **`-n INT` max-NM** — for reads whose computed NM `>= INT`, masks
  every matching M/=/X base (SEQ `->` `=`, quality `-> 0`); the
  emitted NM/MD are unaffected.

**`import` deferred features**:

- **`--i1` / `--i2`** index-read inputs (the index-as-aux BC/QT shape).
  v1 wires `-0/-1/-2/-s` and the positional shapes; index files would
  attach a BC:Z and QT:Z tag computed from a separate index FASTQ. The
  parser scaffolding is in place; just needs a third walker.
- **`-i` CASAVA** parsing (extract barcode from CASAVA-style headers).
  We do parse the description tail for SAM aux fields directly, which
  covers the common case where the FASTQ was produced by `samtools
  fastq`; CASAVA-format input is a follow-up.
- **`--barcode-tag` / `--quality-tag`** renaming of the BC/QT tag pair.
  Not exposed in the v1 CLI (defaults to BC/QT).
- **`-O` / `--output-fmt`** and **`-@` / `--threads`** —
  CLI-accepted-and-ignored stubs. v1 picks output format from the
  output-path extension (`.sam` vs `.bam`) and is single-threaded.

**Note on `import --skipBamQ`.** There is no `--skipBamQ` flag in upstream
`samtools import` (getopt string `1:2:s:0:bhiT:r:R:o:O:u@:N`, long options
`--i1/--i2/--r1/--r2/--rg/--rg-line/--order/--barcode-tag/--quality-tag/
--name2/--no-PG`). The tool reads FASTQ, never BAM, so there is no "BAM
pair-position recovery" path to port — nothing to do here.

**Validation:** small hand-built fixtures live under
`tools/samtools/testdata/parity/{calmd,import}/` covering the four
calmd code paths (match, mismatch, deletion, insertion+softclip) and
the six import shapes (-0, -1/-2, -s, single positional, two
positionals, -T aux extraction, --order, -R/-r RG). The calmd BAQ
path additionally diffs the `-r` / `-rA` output against htslib's
vendored `realn01_exp*.sam` goldens, and `pkg/htsgo/baq` carries the
full `realn0{1,2,3}` golden corpus. The upstream `bam_md.c` /
`bam_import.c` regression cases now run as LIVE parity gates
(`TestParity_Calmd_UpstreamCorpus`, `TestParity_Import_UpstreamCorpus`):
both sides emit plain SAM (sidestepping the BGZF/libdeflate byte-identity
problem) and the streams are compared byte-for-byte modulo the `@PG`
line. Two port bugs were surfaced and fixed by these gates: (1) `import`
of a `/1`/`/2`-suffixed FASTQ record must set `FMUNMAP` (the htslib FASTQ
reader sets it from the suffix alone) — our `-s`/`-0` paths previously
omitted it; and (2) `calmd` must remove a *differing* MD/NM tag and
re-append it at the end of the aux list (leaving an *unchanged* tag in
place), matching `bam_md.c`'s `bam_aux_del`+`bam_aux_append` ordering.

**`phase` upstream-schema emit (DONE — byte parity, deterministic DP).**
The byte-faithful upstream `phase` text stream (CC banner + PS / FL /
M / EV / `//`) is ported across `phase_algo.go`, `phase_emit.go`,
`phase_frag.go`, `phase_khash.go`, `phase_ksort.go`, and
`phase_pileup.go`.

**There is no MCMC in upstream `phase.c`** — earlier roadmap notes
that referred to a "Markov-chain-Monte-Carlo `phase_core` loop" and a
"greedy same-vs-opposite vote / tied junctions emit label 0" were
inaccurate. Upstream's actual phasing is fully deterministic:

1. `phase()` (phase.c:401) builds the consensus `cns[]`, then calls
   `count_all` → `dynaprog` (phase.c:163) — a **Viterbi-style dynamic
   program** over `k`-bit local-haplotype states — to fill an
   `int8_t *path` giving each het's hap0 assignment.
2. `fragphase()` (phase.c:211) assigns each fragment to a haplotype
   vs. `path`, and (when `FLAG_FIX_CHIMERA`, the default) finds the
   best per-read flip point (`FLIP_PENALTY`/`FLIP_THRES`) to repair
   chimeric reads, emitting `YF:i:1` for flipped frags.
3. `genmask()` (phase.c:302) emits `FL` masked regions; the per-site
   `M0/M1/M2` lines and `EV` evidence reads follow.

The Go port replicates each step line-for-line, including the
**in-place Cuckoo-style khash kick-out rehash** (`phase_khash.go`,
`kh_resize` semantics) and the unstable `ks_introsort_rseq` so the
EV-line tie order over equal-`vpos` fragments matches upstream
byte-for-byte even after the fragment table grows past 16 buckets.
`UpstreamSchema` is the CLI default; `-b`/`-F`/`-A` route through
`dump_aln` with the in-tree `drand48` port; `-l FILE`/`-e` implement
`loadpos`/`FLAG_LIST_EXCL`.

**Empirical parity (default LOD):** byte-identical full-stream output
vs. the upstream binary across 300 randomized complex fixtures
(varied lengths, 4–12 hets, dense piles, chimeras, MAPQ=0 reads) and
the `-F`/`-A`/`-k`/`-Q`/`-D` flag matrix. Pinned by `TestLivePhase`
(simple) and `TestLivePhaseComplex` (chimera + FL + table-grow) in
`tools/samtools/pkg/samtools/`, plus binary-free
`TestUnitFragKhashKickoutLayout` / `TestUnitFragKhashLookupAfterGrow`.

> Note: a separate `phaseLegacyTSV` / `phaseHets` path (a greedy
> adjacent-het chainer emitting a simplified PS-label TSV) still exists
> behind `UpstreamSchema=false`, used only by the in-process v1 unit
> tests. It is NOT the CLI path and is not the upstream emit.

**`phase` het-CALLING at low LOD — the errmod/gl2cns LOD precision part
is FIXED.** Earlier this was filed as an "errmod-precision" matter: when
`-q` is lowered well below the default 37 the *set* of het sites called
diverged from upstream, and the suspicion was that the `errmod`/`gl2cns`
genotype-likelihood LOD differed at the margin. That suspicion was only
half right, and the errmod half is now resolved:

- **Root cause of the errmod-level divergence (FIXED).** Upstream
  `errmod_cal` declares the per-genotype likelihood sum `tmp1` as a C
  `float` and does `tmp1 += aux.bsum[k]` against a `double`, re-rounding
  the running sum to single precision at every step. The Go port had been
  accumulating in `float64` and rounding once at the end, which shifted a
  small fraction of float32 `q` entries by 1 ULP — enough to flip a few
  het calls at a low LOD threshold. The port now accumulates in `float32`
  with an explicit `float32(float64(tmp1)+bsum[k])` step, exactly
  mirroring upstream. (A second, latent 1-ULP source was also removed: the
  `-10/ln(10)` phred coefficient is now a runtime `double` division so the
  beta table is bit-identical to upstream rather than 1 ULP off from Go
  constant folding.) With these two fixes the errmod q-matrix is
  **bit-for-bit identical** to the vendored htslib `errmod.c` for every
  ≤255-deep pileup across the full quality range (verified over tens of
  thousands of random columns and pinned by the `TestUnit*` oracle
  goldens in `pkg/htsgo/errmod`). `samtools phase` now matches upstream
  byte-for-byte across a `-q` sweep from 1 to 37 on Phred-40 fixtures
  (`TestLivePhaseLowQ`).
- **Residual `phase`-specific gap (CLOSED 2026-06).** The low-`-q`
  marginal-base divergence is fixed. Root cause was a **flag-wiring
  bug**, not a pileup/errmod issue: upstream phase's getopt string is
  `"Q:eFq:k:b:l:D:A"`, where `-q` is the **min het Phred-LOD**
  (`g.min_varLOD`, default 37) and `-Q` is the min base quality
  (`g.min_baseQ`, default 13). Our CLI had bound `-q` to a spurious
  `minMAPQ` field (phase has **no** MAPQ CLI flag — MAPQ only enters via
  `min(baseQ, mapQ)` during het detection and the per-read
  `core.qual==0` skip in fragment build), and `runUpstreamPhase`
  **hardcoded** `minVarLOD: 37`. So `-q` never reached the variant-column
  test `(c&0xffff)>>2 >= g.min_varLOD` (phase.c:758); the het-LOD
  threshold was pinned at 37. With marginal Q15 bases a het column's
  errmod/gl2cns LOD is 21 (byte-identical to upstream, `c=00060057`),
  which is `< 37`, so every marginal het was dropped at *every* `-q` —
  emitting no `M` block at all (the doc's reported `M0` singleton was the
  pre-investigation symptom). Fix: added `PhaseOptions.MinVarLOD`
  (negative = upstream default 37; `>= 0` used verbatim so `-q 0` works),
  wired it through `runUpstreamPhase`, and bound the CLI `-q` to it. Now
  byte-parity across a full `-Q`×`-q` sweep straddling the LOD boundary
  (`TestLivePhaseMarginalQ`, plus binary-free `TestUnitAdmitPhaseBase` /
  `TestUnitIsPhaseVariantColumn`); the default-LOD oracles
  (`TestLivePhase`, `TestLivePhaseLowQ`, `TestLivePhaseComplex`) still
  pass unchanged.
- **libm boundary note (informational).** At the *double* level the beta
  table's interior entries use `exp`/`log1p`, and Go's `math.Exp` differs
  from glibc's `exp` by at most 1 ULP for a handful of arguments. That
  sub-ULP double noise is far below float32 resolution and is fully
  absorbed when bsum terms round into the float32 `q` output, so it never
  changes a result. It is a genuine libm transcendental last-ULP
  property, not a fixable algorithmic difference; see `docs/UPSTREAM_BUGS.md`.

**`targetcut` HMM consensus mode** (implemented). The Go port is now
a faithful translation of upstream `cut_target.c`: per-position
consensus via the MAQ revised error model (`errmod.c` is ported
in-tree as a shared package at `pkg/htsgo/errmod` — both
`samtools targetcut` and `bcftools mpileup` import the same
implementation, eliminating the duplicate ports that briefly
existed in PRs #199 and #216), followed by a 2-state Viterbi
over the per-chrom consensus track to segment "covered, callable"
regions from "no-info or uninformative" regions, then one SAM-
format consensus record is emitted per identified region in the
exact upstream printf shape
(`%s:%d-%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t<seq>\t<qual>\n`). The
emit-loop's "position 0 never participates in a region" quirk —
a consequence of upstream's backtrack loop running over the half-
open range `(0, l-1]` — is reproduced verbatim so we agree with
the C output by construction. CLI flags `-Q -i -0 -1 -2` carry
their upstream semantics; `-f` reference triggers per-record BAQ
realignment via `pkg/htsgo/baq.SamProbRealn` (apply+extend mode,
flag `1<<1|1`), matching upstream `cut_target.c::read_aln`. The
pre-port "aligned-slice FASTA" mode remains available behind
`--simple` (library: `TargetcutOptions.SimpleMode`); `-f` has no
effect in simple mode.

One implementation detail diverges from upstream by design: when
n > 255 bases pile at a single position the upstream `errmod_cal`
shuffles via `ks_shuffle` (drand48 state) and truncates to 255 to
fit its pre-computed coefficient table. We truncate deterministically
to the first 255 because the downstream gencns caps depth at 255
anyway and a drand48-dependent output would not be reproducible.
For coverages ≤255 (the overwhelming common case in practice) we
are byte-equivalent to upstream.

**Validation:** hand-built SAM fixtures in
`tools/samtools/pkg/samtools/phase_test.go` and
`tools/samtools/pkg/samtools/targetcut_test.go`. Phase tests cover
single-block chaining (consistent & label-flipping orderings),
ambiguous-label fall-back when reads don't bridge two hets, and
the MinMAPQ filter. Targetcut tests cover the simple-mode legacy
behaviour AND the HMM mode: a uniform-coverage region (one
emitted region with the expected SAM shape), an entirely-empty
chrom (no output), the upstream read filter (unmapped /
secondary / qcfail / dup all skipped), MinBaseQ pushing every
cell to "no info", a tuned `-i` entry penalty separating two
coverage blocks into two regions, and a majority-vote consensus
base check. A self-consistency test on the errmod port checks
that 10 identical 'A' observations yield the minimum (best)
homozygous score at q[A,A]. There is no upstream regression-test
fixture for either tool in `reference_code/samtools/test/`, so
expected values are hand-derived directly from the C source — the
project's accepted standard for tools with no upstream fixture.

**`consensus` bayesian mode** (implemented). Upstream `bam_consensus.c`
ships five modes — `simple` and four bayesian flavours (`bayesian_r`
aka "bayesian", `bayesian_m`, `bayesian_p`, `bayesian_116`). All five
are now implemented. The Gap5-derived posterior caller
(`calculate_consensus_gap5` / `_gap5m`), the localised-MAPQ NM-halo
adjustment (`nm_init` / `nm_local` / `poly_len`), and the bayesian
knobs (`-C/--cutoff`, `--P-het`, `--P-indel`, `--het-scale`,
`--adj-qual`, `--use-MQ`, `--adj-MQ`, `--NM-halo`, `--SC-cost`,
`--scale-MQ`, `--low-MQ`, `--high-MQ`, `--default-qual`,
`-p/--homopoly-fix`, `--homopoly-score`, `--homopoly-redux`) are
ported faithfully; the upstream `fast_exp` 0.1-resolution quantization
and the degree-3 `fast_log2` are reproduced so phred conversions are
byte-exact. The default invocation (`MODE_RECALL`, MQUAL + NM-adjust
on) is byte-for-byte parity with upstream on the
`reference_code/samtools/test/consensus/` corpus for the non-`-T`
fixtures (18q/19q/18p/19p/20p/21p on consen1, 30/31/32/40/41/42 on
consen1c).

**Indel calling** (insertions and deletions) is implemented and
byte-faithful for both `simple` and `bayesian` modes across
FASTA/FASTQ/pileup. Insertion-column rows (`nth>0`) are emitted by a
single per-mode dispatcher (`consensusInsertionColumns`) that the
FASTA/FASTQ and pileup emitters share. The insertion-column membership
rule is ported from upstream's pileup engine
(`consensus_pileup.c::get_next_base`): a read whose alignment terminates
exactly at the reference position is removed before the insertion column,
so it neither inserts nor pads there (`spansInsertionColumn`). Deletion
columns (`*`) honour `--show-del` for both the `nth==0` row and, crucially,
do NOT suppress the following `nth>0` insertion columns — upstream invokes
its emit callback independently per `nth`. The simple-mode gap bucket is
quality-weighted under `--use-qual` (`score[16] += 8*q`), the per-event
min-qual gate is applied to gap/pad events as well as bases, and the
pileup display columns (seq/qual/depth) are the RAW pileup column,
unfiltered by `--min-BQ` (which affects only the consensus call) —
all matching `bam_consensus.c`. A live parity sweep
(`TestConsensus_IndelUpstreamParity`) compares the Go port against the
freshly-built upstream binary byte-for-byte over an indel-rich fixture
across both modes, all three formats, and the `plain`/`mark-ins`/
`show-del`/`ambig`/`all-pos`/`no-show-ins`/`use-qual`/`min-bq` flag
variants.

Remaining indel-adjacent gap: **none** — closed.

- **`-a/--all-positions` in pileup format** (DONE). The placeholder
  `<chrom>\t<pos>\t0\t0\tN\t0\t*\t*` rows that upstream prints for
  deletion-only and zero-coverage reference positions are now emitted,
  reproducing upstream's `empty_pileup2` byte-for-byte — including its
  quirky duplicate rows at deletion sites (a suppressed `'*'` column does
  not advance `last_pos`, so each deletion position is re-filled by every
  following column: a D-bp deletion run yields the position emitted
  `D, D-1, … , 1` times). The fill is driven by a per-window
  `last_pos` cursor that mirrors `basic_pileup`/`empty_pileup2`
  (`bam_consensus.c:2202`/`2832-2842`): genuine pileup columns trigger a
  lazy gap fill back to the last emitted row, and a tail fill closes out
  the window. Leading/internal/trailing gaps, `-aa` empty contigs, and
  the `-l`/BED filter (placeholders are confined to selected positions)
  are all handled. Validated by the live
  `TestConsensus_AllPositionsUpstreamParity` sweep (simple + bayesian
  modes × `-a`/`-aa` × `--show-del` on/off, byte-for-byte against the
  freshly-built upstream binary) plus deterministic unit tests
  (`TestConsensus_AllPositionsPileup_*`).

Genuinely deferred sub-knobs (precisely scoped):

- **`-t/--qual-calibration` and `-X/--config`.** Accepted on the CLI
  but apply only the FLAT identity calibration table. The per-machine
  calibration tables (HiFi/HiSeq/ONT/Ultima) and the QUAL-file parser
  (`load_qcal`, `bam_consensus.c:672-736`) are a separable ~300-line
  table/parser block that does not affect the default invocation.
- **`-T/--reference` uncovered-base fill — DONE.** The reference FASTA is
  now loaded and used to fill every no-coverage / gap position that would
  otherwise be 'N', mirroring upstream `bam_consensus.c` (`update_ref` +
  `empty_pileup2` + the `basic_fasta` / `ref_or_Ns` gap fills). Scope of the
  substitution, byte-validated against the upstream binary
  (`consensus_reference_upstream_test.go`):
  - **pileup** (`--format pileup`): the no-coverage placeholder rows emit
    the reference base in the call column (`rseq ? rseq[i] : 'N'`); depth
    (0), quality (0) and the seq/qual (`* *`) columns are unchanged.
  - **FASTA**: internal gaps, and — under `-a` — leading/trailing gaps, are
    filled from the reference; the covered span keeps its computed call
    (a low-depth 'N' call stays 'N'). The default (no `-a`) still emits only
    the covered span.
  - **FASTQ**: as FASTA, plus the gap quality is `--ref-qual + '!'`
    (`opts->ref_qual`, default 0); covered positions keep their computed
    quality. `--ref-qual` is now honoured.
  - **`-aa` empty contigs**: the whole contig is filled from the reference
    (pileup rows and FASTA/FASTQ records alike).
  - A reference missing a contig, or a position past the contig end, falls
    back to 'N' exactly like upstream's `update_ref` <0 path. The substituted
    bytes preserve the FASTA's original case (soft-masked lowercase verbatim).
  Binary-free `TestUnitConsensusRefBase` /
  `TestUnitWriteEmptyPileupRowsRefSubstitution` cover the pure helpers.
- **Mate-overlap dedup.** `--ignore-overlaps` is accepted but is a
  no-op; v1 counts each mate independently in the pileup walker.
- **Threading.** `-@/--threads`, `-Z/--block-size`, and
  `--input-fmt-option` are accepted but ignored; v1 is single-pass
  and single-threaded.
- **Read-flag filtering.** `--rf/--incl-flags` and `--ff/--excl-flags`
  are accepted as text/int but ignored. v1's filter set is fixed
  (drop UNMAP|SECONDARY|QCFAIL|DUP, matching upstream's default
  `excl_flags`).
- **`--het-only`** restricts output to HETEROZYGOUS-called positions
  (homozygous and no-call positions become `N` in FASTA/FASTQ — with
  coordinates preserved — and are omitted entirely in pileup). This is a
  DELIBERATE divergence from upstream, which parses `--het-only` into its
  options struct but never reads it (a dead-option bug — the flag is inert
  upstream through samtools 1.22). We implement the intended
  heterozygous-only filtering the flag name promises. Het-ness is
  determined independently of `--ambig`. See
  [docs/UPSTREAM_BUGS.md#samtools-consensus-het-only](UPSTREAM_BUGS.md#samtools-consensus-het-only)
  for the bug write-up. Pinned by `TestConsensus_HetOnly_*` (unit) and the
  live `TestConsensus_HetOnlyUpstreamBug` (which proves upstream ignores
  the flag and our output correctly diverges).
- **`--verbosity`** is accepted and ignored.
- **Indel calling (REMAINING).** `samtools consensus` only emits
  substitution/base consensus; the gap5-derived posterior **indel** caller
  (the insertion/deletion-aware branch of `calculate_consensus_gap5`) is a
  larger separate algorithm and is **not** ported. This is the one
  consciously-left-out consensus item; the per-subcommand sub-features PR
  explicitly scoped it out (it is not a small/medium knob like the items
  above).

**`consensus` correctness model.** v1 mirrors upstream's
`calculate_consensus_simple` (`bam_consensus.c:1900-2006`) bit-for-bit
where it matters:

- One fraction gate only: `used_score < call_fract * tscore`
  (`bam_consensus.c:1988-1994`). There is **no** separate
  "min-fraction on the dominant base alone" gate; an earlier PR
  fabricated one and is corrected here.
- Heterozygous condition is `score2 >= het_fract * score1 && ambig`
  with no `score1 > 0` guard (`bam_consensus.c:1982`).
- `use_qual=0` by default (`bam_consensus.c:2984`) — bases score by
  frequency, not quality, until `-q/--use-qual` is set.
- `--show-del` is honoured in pileup mode too — rows whose call is
  `'*'` are suppressed when `--show-del no`
  (`bam_consensus.c:2244`).
- Insertion gating uses `MinCallFraction`, the same knob as the
  per-position gate.

**Validation:** table-driven hand-built SAM fixtures in
`tools/samtools/pkg/samtools/consensus_test.go` covering:
all-match FASTA/FASTQ/pileup, mixed at the 0.75 boundary,
`--show-del no` in pileup (no `*` rows),
`--show-del yes` in pileup (keeps `*` rows),
the canonical **50/30/20 + `-A`** fixture (must land `M`, not `N`),
frequency-only counting (low-Q vs high-Q gives identical output when
`UseQual=false`),
the `UseQual=true` flip (high-Q minority beats low-Q majority),
multi-contig, `-a` zero-fill, `--min-depth`, line-len wrapping,
insertion include/suppress + `--mark-ins`,
and the default-bayesian invocation emitting no fallback warning.
Bayesian-mode parity is additionally verified byte-for-byte against
the vendored upstream `reference_code/samtools/test/consensus/`
golden files (`TestConsensus_BayesianUpstreamParity`): twelve
fixtures spanning FASTQ, pileup, the `-C` cutoff, `-A` ambiguity,
`-a`/`-aa`, and the default MQUAL + NM-adjust path.
Coverage of the `pkg/samtools` package after this PR is ~80%.

### `bcftools`

**Status:** all 24 subcommands present (~96%). `view`, `index`, `stats`,
`query`, `concat`, `norm`, `call` (consensus + **full** multi-allelic),
`annotate`, `head`, `isec`, `merge`, `reheader`, `sort`, `convert`
(pass-through + GEN/HAP/TSV/gVCF modes), `mendelian`, `gtcheck`, `roh`,
`filter` (incl. `-M/--mask-file`), `consensus` (chain + iupac),
`mendelian2`, `polysomy`, `cnv` (incl. `--AF-file`), `csq` (slices 1–4),
and `mpileup` (SNP slices 1–4 + legacy `bam2bcf_indel` + `--indels-cns`).

All bcftools subcommands now have a real implementation in the Go port.
`gtcheck` is feature-complete: `-c/--cluster` clustering (this port's own
design — upstream leaves it an error stub), `--n-matches`,
`--distinctive-sites`, `-i/-e` filter expressions, and both `-O` output
containers (`t` text and `z` BGZF with an optional `z<N>` level) are all
implemented and live-oracle validated. `convert` covers the full upstream
mode set, including `--gvcf2vcf` (and its prefix-abbreviation `--gvcf`;
upstream `convert` has no separate VCF→gVCF blocking mode — that lives in
the `+gvcfz` plugin). `convert` PLINK exporters (`--plink`/`--tped`/`--plink-bed`, PLINK1 spec),
`csq -l/--local-csq`, `concat --ligate`, and `som` train/classify (upstream
`fwrite`-return write bug fixed) are implemented; samtools `tview` text/HTML
(`-d T`/`-d H`) and the interactive `-d C` viewer (pure-Go Linux raw-mode
termios, no ncurses) are all implemented.
`query %N_ALT` / `import --skipBamQ` are **not** upstream flags. The port now
rejects `%N_ALT` (and any undeclared bare tag) at header-validation time with
upstream's exact "no such tag defined in the VCF header: INFO/<TAG>" error
(`TestParityQuery_NAlleles`), rather than silently emitting ".".

Recently closed parity gaps, now asserted byte-for-byte against the upstream
binary (fixtures vendored under `tools/bcftools/testdata/parity/`):

- `concat -a`/`--allow-overlaps` matches upstream's contig-ordering
  heuristic exactly — the synced reader visits contigs in **first-seen-in-data
  order across the inputs** (not the merged-header `##contig` order), with
  ties broken reader-by-reader (`TestParityConcat_AllowOverlaps`).
- `concat -D`/`-d` now requires `-a`, erroring with upstream's
  "The -D option is supported only with -a" when used standalone, and the
  `-a -D`/`-a -d {exact|snps|indels|both|all}` cross-file de-duplication
  matches upstream's `BCF_SR_PAIR_*` collapse logic
  (`TestParityConcat_DedupRequiresA`, `TestParityConcat_DedupAllowOverlaps`).
- `norm -f` left-alignment and `norm -c {e|w|x|s}` (the upstream bitmask:
  warn / exclude / set-REF, with `e` exclusive) now match byte-for-byte,
  including the REF/ALT swap with genotype re-indexing for `-c s`
  (`TestParityNorm_LeftAlign`, `TestParityNorm_CheckRef*`). Note the port's
  `-c` letters were corrected: `s` is **set/fix** (not skip) and `x` is
  exclude, matching upstream `vcfnorm.c`.
- CSI indexing of an upstream-produced BCF is functionally parity-tested
  (`TestParityIndex_CSIForBCF`): our reader handles htslib's optional int64
  FORMAT descriptors, so region reads return identical records. Byte-identical
  `.csi`/`.tbi` output remains a deliberate non-target (BGZF framing differs,
  and our BCF CSI carries a small tabix-style aux block htslib omits for BCF).

Recently closed (2026-06-17 parity-pipeline bugfix wave), each asserted
byte-for-byte against the upstream binary (live-oracle tests, never `t.Skip`):

- `view -c/-C/-q/-Q` (and `-x/-X`) now recompute and append `INFO/AC` and
  `INFO/AN` even without a sample subset, matching `vcfview.c`'s `calc_ac`
  path (the non-subset path prefers pre-existing INFO/AC,AN, then GT;
  `-I/--no-update` still suppresses it). `TestView_ACANRecompute_UpstreamParity`.
- `norm -m+` join now matches `vcfnorm.c`: biallelics at the same position are
  bucketed by variant-type category (or all together for `-m+any`), merged via
  a port of htslib `merge_alleles` (common padded REF + allele-index map) and
  `merge_format_genotype` (so `0/1`+`0/1` joins to `2/1`, indels with
  differing-length REFs share a padded REF). `TestNorm_JoinMultiallelic_UpstreamParity`.
- `annotate -x` now strips the matching `##FILTER`/`##INFO`/`##FORMAT` header
  lines (bare `FILTER`/`INFO` drop all; `FORMAT`/`FMT` keep GT; `FILTER/NAME`
  rewrites an emptied record to PASS), fixing the header line ordering for
  `-x FILTER`. `TestAnnotate_Remove_UpstreamParity`.
- `merge --force-samples` added: duplicate sample names across inputs are
  de-duped by prefixing the clashing name from input *i* with `<i+1>:`
  (`A + A -> A, 2:A`), matching `vcfmerge.c merge_headers`.
  `TestMerge_ForceSamples_UpstreamParity`.
- `concat -a` overlap-merge record ordering now matches the synced reader
  (`bcf_sr_sort.c`): records at a shared position are grouped by REF>ALT and
  the groups emitted by descending pre-dedup count, ties by first-appearance.
  `TestConcat_OverlapOrder_UpstreamParity`. (The `-d snps/indels/both` collapse
  *model* — which records are dropped — remains a separate, pre-existing
  divergence from upstream's variant-set pairing.)
- `consensus` default het handling: when the VCF has samples and neither `-H`
  nor an allele pick is given, upstream applies IUPAC ambiguity codes
  (`consensus.c` `iupac_GTs`) across the `-s` sample or all samples; the port
  did this only with `-I`. Now matched. `TestConsensus_DefaultIUPAC_UpstreamParity`.
- `roh` ST/RG output is now ordered chromosome-major then sample-major (header
  order), matching upstream's per-chromosome synced-reader flush.
  `TestRoh_RecordOrder_UpstreamParity`.

**`csq` and the parity pipeline's synthetic GFF:** the pipeline's
`bcftools csq --force` matrix entry diverges, but the cause is the *fixture*,
not the consequence engine. The synthetic `annotations.gff3` has invalid GFF3
phase columns (`.`-phase exons, CDS phase ≠ len%3); upstream's GFF parser
detects this and skips/truncates the offending transcripts ("inconsistent
phase column ... skipping", indexing 689 CDSs not all 800), so it emits
different/fewer consequences (e.g. `intron|gene00000` where the transcript was
dropped). Our port does not yet validate the GFF phase column, so it keeps
those transcripts and calls `missense`. On clean, valid-phase GFFs — single
gene and overlapping opposite-strand genes alike — csq matches upstream
byte-for-byte, and the extensive per-tool csq parity suite passes. The
remaining narrow gap is GFF phase-column validation (the
`phase != len%3` / `.`-phase skip logic in upstream's GFF reader).

**Multi-threaded output compression** (`-@ / --threads N`) — DONE for the
output-writer subcommands. Like upstream (which calls htslib
`hts_set_threads` on the output file), the Go port now performs genuine
parallel BGZF compression of bgzipped output when `N > 1`, reusing the
shared block-parallel `pkg/htsgo/bgzf.MultiWriter` (the same writer used
by `bgzip` and samtools threading). The shared `openOutput` writer path
selects `MultiWriter` for `-O z` (VCF.gz) and `-O b` (BCF) whenever
`Threads > 1`, so the flag is wired end-to-end through **`view`,
`concat`, `norm`, `call`, `annotate`, `sort`, `merge`, `isec`,
`convert`, `reheader`, `mendelian`, `mendelian2`, `csq`,
`filter`/`vcffilter`, and `plugin`**, and (via its own writer)
**`mpileup`** output. As part of this change `-O z` now emits
BGZF (gzip-compatible) rather than plain `compress/gzip`, matching
upstream's `.vcf.gz` framing. Output decodes byte-identically regardless
of the thread count (every BGZF block is an independent gzip member);
compressed bytes may differ at block boundaries, so parity is asserted on
the **decoded** records. The pileup/call computation itself remains
single-threaded in v1 — only the output-compression stage is
parallelised, which is the dominant win for `-O z`/`-O b` workloads.

**Threading tail wired (follow-up complete):** the remaining
`openOutput` callers — `isec`, `convert`, `reheader`, `mendelian`,
`mendelian2`, `csq`, `filter`/`vcffilter`, and the `plugin` output path —
now route their bgzipped output through the same `bgzf.MultiWriter`
chokepoint when `Threads > 1`, so `-@/--threads` is honoured end-to-end
for every BGZF-framed-output bcftools subcommand. (mendelian2's
`-W`-only BGZF-VCF path was also folded onto `newBGZFOutput`, so it
threads and flushes the header into its own block too.) The remaining
subcommands — `gtcheck`, `roh`, `cnv`, `polysomy`, `consensus`, and
`query`/`stats`/`index` — either have no BGZF-framed record output or do
not accept `-@`.

Additionally, the shared `view`/`call`/etc. `openOutput` now flushes the
header into its own BGZF block for `-O z`/`-O b` (matching upstream
htslib's `bgzf_flush` after `vcf_hdr_write` / `bcf_hdr_write`), bringing
it in line with the `mpileup` writer and keeping tabix/.csi virtual
offsets clean.

**Validation:** live upstream-binary parity (`view_threads_test.go`, the
uniquely-named `upstreamBcftoolsThreads` sync.Once builder, never
`t.Skip`): our `-@ {2,4,8}` output decodes byte-identically to our `-@ 1`
output for both `-O z` and `-O b`, and the decoded records match the live
upstream `bcftools view -O z` / `-O b --threads 4` output on a
multi-block VCF fixture. The threading tail adds `thread_tail_test.go`
(uniquely-named `upstreamBcftoolsThreadTail` sync.Once builder, never
`t.Skip`): `convert`/`reheader`/`filter`/`isec`/`mendelian2` `-@ {2,4,8}`
output decodes byte-identically to `-@ 1` (and matches the live upstream
binary where an `-O z` surface exists), plus structural assertions that
`view -Oz`/`-Ob` now place the whole header in its own first BGZF block.
Race-clean (`go test -race`). The isolated
parallel-compression speedup (~3x at 4 threads) is demonstrated by
`BenchmarkMultiWriter` in `pkg/htsgo/bgzf`; end-to-end `view` throughput
moves little because VCF text (de)serialisation, not deflate, dominates
that path (`BenchmarkViewThreadsVCFGz`).

Closed in the #220–#225 wave: `view -x/--private` & `-X/--exclude-private`
(private-allele site filter), `stats -u/--user-tstv` (user-defined Ts/Tv
binning), and `csq -b/--brief-predictions` & `-C/--genetic-code` (standard
table 0).

The former "boulders" are now **closed**: mpileup indel calling (both the
legacy `bam2bcf_indel` path and `--indels-cns`), the `convert` GEN/HAP/TSV/
gVCF modes plus the PLINK exporters (`--plink`/`--tped`/`--plink-bed`,
implemented to the PLINK1 spec since upstream comments those options out),
and csq slice 4 (FORMAT/TBCSQ, `--unify-chr-names`, `--dump-gff`,
`-O b|u|z` non-text output **and `-O t` streaming-text output**) all
landed and are live-oracle validated (see the per-subcommand sections
below). The only bcftools item still open is `gtcheck`'s `-c/--cluster`
dendrogram + filter expressions. (`csq -l/--local-csq` is now
implemented; `csq -O t` is now implemented, matching upstream
`text_print_vcsq` byte-for-byte — the residual `csq -s -` /
sample-subsetting gap is unaffected.)

The plugin system (`bcftools plugin` / `bcftools +<name>`) is **done**,
but with a deliberate design divergence from upstream:

- **`+plugins`** — implemented as a **subprocess plugin system**.
  Upstream loads plugins as native shared objects (`.so`) via `dlopen`
  against a fixed C ABI. The Go port instead resolves `bcftools +<name>`
  to an ordinary *executable* found by name in the `BCFTOOLS_PLUGINS`
  colon-separated directory list, runs it as a child process, pipes the
  input VCF as uncompressed text to its stdin, and reads VCF back from
  its stdout. A plugin is therefore "a VCF-on-stdin to VCF-on-stdout
  filter" and can be written in any language — no C ABI, no version
  check, no rebuild against the host. `bcftools plugin -l`/`-lv` lists
  discoverable plugins; the host applies `-o`/`-O` formatting around the
  plugin's output; a non-zero plugin exit is surfaced as an error with
  its stderr. The contract is specified in `docs/PLUGIN_PROTOCOL.md`.
  The mechanism lives in `tools/bcftools/pkg/bcftools/plugin.go` with a
  reference example plugin under `tools/bcftools/plugins/example/`.
  **Now fully ported:** all 41 bundled upstream plugins (`+fill-tags`,
  `+split-vep`, `+setGT`, `+prune`, `+fixploidy`, `+vrfs`, ...) are
  reimplemented in pure Go and dispatched by `+<name>` ahead of the
  subprocess lookup, each driven to CLI-to-CLI byte-parity against the
  real upstream `bcftools` 1.23.1 binary (see the "FULLY CLOSED" summary
  at the top of this file). The subprocess protocol is retained as the
  **fallback** for user-supplied executables in any language. Upstream
  plugin sources (`plugins/*.c`) remain vendored under
  `reference_code/bcftools/` as the parity reference.

- **Native plugin region/target selection (`-r/-R/-t/-T`)** — **done**.
  A shared host-side filter (`tools/bcftools/pkg/bcftools/region_target.go`)
  consumes the region/target options out of each native plugin's argv and
  applies them before any record reaches the plugin's `Process`, so no plugin
  re-implements (or rejects) the flags. `-r`/`-R` is span-OVERLAP based, `-t`/`-T`
  is record-START based with `^` negation — replicating upstream 1.23.1's exact
  `-r` vs `-t` difference (an indel at POS=100 spanning 100..104 is included by
  `-r chr:102-102` but excluded by `-t chr:102-102`). `-R`/`-T` files use the
  synced-reader TSV format (`.bed` = 0-based; otherwise 1-based `chr,pos` or
  `chr,beg,end`). Supported by check-sparsity, remove-overlaps, prune,
  smpl-stats, indel-stats, contrast, guess-ploidy (only `-r/-R`; its `-t` is
  `--tag`), mendelian2, trio-stats, isecGT (applied to both input streams),
  split and scatter. Plugins opt in via a `RegionTargetCaps` capability so the
  letters that other plugins repurpose (tag2tag's `-r`=`--replace`/`-t`=`--tags`)
  are left untouched. check-sparsity keeps its own per-region report grouping and,
  per the fix-on-port policy, now accepts BED/TSV `-R` files (which upstream
  silently drops — `docs/UPSTREAM_BUGS.md#bcftools-check-sparsity-regions-file`)
  while keeping colon region-list lines byte-identical to upstream. Byte-validated
  vs the upstream binary in `native_plugin_region_target_oracle_test.go`.

- **Native `+prune` — all modes** — **done**. The native port
  (`native_plugin_prune.go` + `native_plugin_prune_ld.go`) now covers every
  upstream mode, byte-validated vs 1.23.1 in `native_plugin_prune_oracle_test.go`:
  - `-n/--nsites-per-win` in `1st`, `maxAF` (incl. the default maxAF without
    `--AF-tag`) and `rand` selection (`--random-seed` pins the in-tree drand48
    draw order — `native_drand48.go` reused);
  - `-m count=N` cluster removal and `-m R2=/LD=/RD=` (and bare-number r2)
    linkage-disequilibrium thresholding (hard drop), the LD math a byte-exact
    port of `_calc_r2_ld`/`vcfbuf_ld` (`+ - * /` and `sqrt` only);
  - `-f LABEL` soft-filtering (sets FILTER instead of dropping) and
    `-a count|r2|LD|RD` annotation (POS_* + value INFO tags, header lines);
  - `--keep-sites` (-k), `-i/-e` filtering, `--randomize-missing`, and the
    `-w` window in bp and site-count forms.
  Two surprising upstream behaviours are reproduced for parity (maxAF ranks by
  alt/ref, and the soft-filter header renders "within 0kb" via integer
  division); the pre-port code's claims that the `rand`/default-maxAF modes
  "cannot be matched byte-for-byte" were false and are corrected — see
  `docs/UPSTREAM_BUGS.md#bcftools-prune-maxaf-ranks-by-altref-not-allele-frequency`.

- **Native plugin output auto-indexing (`-W/--write-index[=FMT]`)** — **done**.
  A shared helper (`tools/bcftools/pkg/bcftools/native_plugin_writeindex.go`)
  parses `-W`/`--write-index` (bare or `=csi`/`=tbi`) and writes a CSI (default)
  or TBI index next to each indexable output, reusing the in-tree
  `pkg/htsgo/tabix` CSI/TBI writers. Plain-VCF/stdout outputs are non-indexable
  and reproduce upstream's exact error. Supported by contrast, isecGT,
  mendelian2, split and scatter (the multi-output plugins index every emitted
  file). The CSI writer now always emits the trailing `n_no_coor` field even when
  zero, matching htslib `hts_idx_save_core` byte-for-byte. Byte-validated vs the
  upstream binary (decoded index content == `bcftools index` over our data file)
  in `native_plugin_writeindex_oracle_test.go`.

- **Native stats-plugin curly-brace multi-threshold `-i/-e` expansion** —
  **done**. `smpl-stats`, `indel-stats` and `trio-stats` now support upstream's
  `-i 'EXPR{a,b,c}'` syntax: a shared helper
  (`tools/bcftools/pkg/bcftools/native_plugin_filter_expand.go`,
  `expandPluginFilterExpr`) replicates the C `parse_filters()` routine
  byte-for-byte — expanding each `{a,b,c}` list into one filter expression per
  element (braces replaced by the element), combining multiple `{...}` groups as
  a cartesian product in upstream's exact order, treating an empty `{}` as a
  collapse to the single default "all" filter, and erroring on an unmatched `{`.
  Each expanded threshold becomes its own `FLT*`/`SITE*` (and, for indel-stats,
  `SN*`/`DVAF*`/`DLEN*`/`DFRAC*`/`NFRAC*`) report section labelled by the
  expanded expression, with the per-filter `MERR`/`TRANSMITTED` debug lines (for
  trio-stats `-d`) streamed interleaved per record, and the stderr "Collecting
  data for N filtering expressions" note reporting the expanded count. The
  single-filter (no-brace) and no-filter paths stay byte-identical. Byte-validated
  vs the upstream binary in `TestNativePluginSmplStats`, `TestNativePluginIndelStats`
  and `TestNativePluginTrioStats`. (Note: a site-level `-e` brace expression hits
  the same pre-existing upstream NULL-`smpl_pass` segfault as a site-level `-e`
  without braces in smpl-stats/indel-stats; our port handles it robustly, so it is
  oracled only with FORMAT-level `-e` or site-level `-i`.)

- **Native stats-plugin remaining flags (`-o`, indel-stats `-p`, trio-stats
  `-a`, contrast `-f`)** — **done**. The four stats/contrast plugins now close
  every "not supported" mode, byte-validated vs the upstream binary 1.23.1:
  - **`-o/--output FILE`** for `smpl-stats`, `indel-stats` and `trio-stats`: the
    report is written to FILE instead of stdout via a shared
    `statsReportWriter` (`native_plugin_stats_common.go`), mirroring
    report_stats()'s `!output_fname||!strcmp("-") ? stdout : fopen(...)`. The
    bytes are byte-identical to the stdout form (the `CMD` line echoes the
    verbatim argv, including `-o`, in both). `TestNativePluginStatsOutputFile`
    drives both binaries with the SAME `-o` path and compares the FILE contents.
  - **indel-stats `-p/--ped FILE`** de-novo mode: the PED-resolved trios restrict
    the stats to de-novo indels in each child (the same Mendelian/`--alt2ref-DNM`
    DNM test as indel-stats.c), the SN* "number of samples" column reports the
    trio count, `npass`/`npass_gt` count DNM sites/genotypes, and the per-trio
    FORMAT/site filter is folded exactly as upstream. indel-stats.c's `parse_ped`
    is replicated including its *lack* of dedup (a trio listed twice is kept,
    unlike trio-stats.c). Byte-validated in `TestNativePluginIndelStatsPED` on a
    two-trio AD-bearing indel fixture (`trio_indels.vcf`/`.ped`). Fix-on-port:
    upstream aborts (`Incorrect GT allele`, exit 255) on a PED indel VCF lacking
    FORMAT/AD; our port skips the AD-derived DVAF/DFRAC/NFRAC contributions and
    still reports — see `docs/UPSTREAM_BUGS.md` and `TestUnitIndelStatsPEDNoADRobust`.
  - **trio-stats `-a/--alt-trios INT`**: the deferred singleton/doubleton
    (transmission-rate) accounting — a singleton/doubleton is counted only when
    its allele appears in at most `-a` alternate trios at the site, replicating
    alt_trios_reset/alt_trios_add and the final deferred loop (including the
    `-d transmitted` debug lines emitted from it). Byte-validated in
    `TestNativePluginTrioStatsAltTrios` on the two-trio `trio_multi.vcf`.
  - **contrast `-f/--max-allele-freq NUM`** rare-allele enrichment: the per-site
    VCF + PASSOC/FASSOC/NASSOC/NOVEL* output is unchanged, and the region-wide
    pooled minor-allele counts (folded over the records whose minor allele is at
    or below the `-f` threshold, ref/alt columns swapped when REF is the minor
    allele, exactly as contrast.c) feed the extra stderr
    `max_AC/PASSOC/FASSOC/NASSOC:` summary line (`%e` Fisher probability + `%f,%f`
    control/case non-REF fractions). An integer `-f` is a raw allele-count
    threshold; a float in `[0,1]` is scaled by the total sample count (floored,
    min 1). Byte-validated in `TestNativePluginContrastEnrichment` (both stdout
    and the full stderr summary). The `--regions-overlap`/`--targets-overlap`
    region-matching modes are now **supported** too (pos/record/variant overlap
    semantics implemented once in the shared `region_target.go` filter and
    honoured across contrast/mendelian2/scatter/split/trio-dnm3; byte-validated
    in the overlap-mode oracle suite). Contrast has no remaining unsupported
    options.
  - Binary-free `TestUnit*` tests cover the pure helpers: `parseIndelStatsPED`,
    `parseContrastMaxAC`, the contrast enrichment folding, the trio-stats
    deferred alt-trio accounting and the shared `statsReportWriter`.

- **Native `+split-vep` full selection/override surface** — **done**. The
  native split-vep port (`tools/bcftools/pkg/bcftools/native_plugin_splitvep*.go`)
  now closes the last five modes it previously rejected, byte-validated vs the
  upstream binary 1.23.1:
  - **`-s` transcript selectors** beyond all/worst: `primary` (CANONICAL=YES),
    `pick` (PICK=1), `mane` (MANE_SELECT!=""), and an arbitrary
    `<FIELD><OP><VALUE>` EXPRESSION with the `=`, `!=`, `~`, `!~` operators
    (`initSelectTrExpr`/`matchingTranscripts`, porting init_select_tr_expr /
    get_matching_transcript). The EXPRESSION-only/no-`-c` case defaults to
    `drop_sites=1` and reproduces the `-X` "no effect" error.
  - **PRN qualifier** (`:all`/`:worst`): `:worst` rewrites the printed
    Consequence to its single worst `&`-joined term (`csqRewriteWorst`, porting
    csq_rewrite_worst — including upstream's surprising exact-match ranking,
    documented at
    `docs/UPSTREAM_BUGS.md#bcftools-split-vep-prn-worst-exact-match`).
  - **`-g/--gene-list [+]FILE`** restrict vs prioritise modes and
    `--gene-list-fields` (`initGeneList`/`restrictCsqsToGenes`, porting
    init_gene_list / restrict_csqs_to_genes, including the two-pointer partition
    order).
  - **`-S/--severity -|FILE`** custom severity-scale override (the file-based
    scale re-orders worst-transcript selection and the `:term[+|-]` filter);
    `-S -`/`-S ?` print the default scale to stderr and exit non-zero.
  - **`--columns-types -|FILE`** regex type table replacing the built-in
    presets (drives both the `##INFO` Type and numeric re-parsing); `-`
    prints the default table to stderr and exits non-zero. A bad FILE errors
    only when an untyped column actually needs it, matching upstream's lazy
    `get_column_type`/`init_column2type`.

  Fix-on-port note: the `:csq` severity term lookup is case-sensitive against
  the lowercased scale keys, so `-s :MISSENSE` is rejected exactly as upstream
  (the earlier port lower-cased it and wrongly accepted it). Byte-validated in
  `TestNativePluginSplitVepSelectors`, `TestNativePluginSplitVepGeneList`,
  `TestNativePluginSplitVepSeverity` and `TestNativePluginSplitVepColumnsTypes`
  (fixture `tools/bcftools/testdata/parity/vep_select.vcf`).

- **Native `+setGT` binomial / random / read-depth modes** — **done**. The
  native setGT port (`tools/bcftools/pkg/bcftools/native_plugin_setgt.go`) now
  closes the last three target/new-gt modes it previously rejected, all
  byte-validated vs the upstream binary 1.23.1:
  - **`-t b:TAG CMP VAL`** — the two-tailed binomial selector over a FORMAT
    integer tag (typically `AD`) for each diploid heterozygous genotype. The
    p-value is `calc_binom_two_sided(AD[ia], AD[ib], 0.5)` indexed by the two GT
    alleles, computed through the in-tree `kfBetai` (htslib's `kf_betai`), so it
    is bit-exact — no libm `pow`/`lgamma` in the tail. `CMP` parsing mirrors
    `parse_binom_expr` exactly (`<`, `<=`, `>`, `>=`, `==`/`=`); the tag must be
    declared in the header.
  - **`-t r:FLOAT` + `-s/--seed INT`** — random selection of a `FLOAT` proportion
    of the targeted genotypes. **RNG finding:** setGT seeds htslib's
    `hts_srand48(rand_seed)` from `-s` (default 0) and *nothing else* — no
    `time()`/`getpid()`. `hts_drand48` is the deterministic POSIX/FreeBSD 48-bit
    LCG (`a=0x5DEECE66D, c=0xB`, seed low bits `0x330E`), so the stream is fully
    reproducible. The earlier code comment claiming the RNG "cannot be matched
    byte-for-byte" was **wrong**; the port reimplements the LCG
    (`native_drand48.go`, pinned to the canonical sequence by
    `TestDrand48KnownVectors`, cross-checked against glibc `drand48`) and is
    byte-identical to upstream for any fixed seed and across thread counts
    (random mode forces serial execution so the shared stream advances in input
    order). Used alone, `-t r` implicitly targets all genotypes (setGT.c:271).
  - **`-n X`** — set every allele of the genotype to the allele with the largest
    FORMAT/AD value for that sample (also usable inside a `c:` template, e.g.
    `-n c:0/X`); requires a FORMAT/AD header.

  Note: setGT has **no** `b:e`/`b:f` binomial subtype — `b:` is followed
  directly by `TAG CMP VAL` (a leading `e:`/`f:` is parsed as a tag name and
  errors as "tag not present"). Byte-validated in `TestNativePluginSetGTBinom`,
  `TestNativePluginSetGTRandom` and `TestNativePluginSetGTReadDepth` (fixture
  `tools/bcftools/testdata/parity/gt_plugins.vcf` plus an in-test skewed-AD
  fixture that exercises the genotype-changing binomial path).

- **Native `+fill-tags` — all modes** — **done**. The native fill-tags port
  (`tools/bcftools/pkg/bcftools/native_plugin_filltags*.go`) now closes the four
  modes it previously rejected, all byte-validated vs the upstream binary
  1.23.1:
  - **Every built-in tag.** `AN`, `AC`, `AC_Hom`, `AC_Het`, `AC_Hemi`, `AF`,
    `MAF`, `NS`, `HWE`, `ExcHet`, `END`, `TYPE`, `FORMAT/VAF`, `FORMAT/VAF1` and
    the `F_MISSING` expression — counted by the exact `process_fmt` BRANCH_INT
    het/hom/hemi/half classification (incl. `-d/--drop-missing`). `HWE`/`ExcHet`
    use the in-tree Wigginton 2005 exact test (`calcHWE`), and the sites-only
    `AF`-from-`AN,AC` path (`process_info_af`) is supported. `-t LIST` selection,
    `INFO/`/`FORMAT/` qualifiers and the `all` keyword match upstream; unknown
    tags error with the exact upstream message.
  - **`-S/--samples-file FILE` population grouping.** The file is
    `SAMPLE  GRP1[,GRP2,...]` per line (porting `parse_samples`): each distinct
    group becomes a population whose tags are suffixed `_GROUP`, plus the summary
    `ALL` population (empty suffix, appended last as `init_pops` does). Per-pop
    `##INFO` headers ("... in GROUP") and the per-pop tag values match upstream
    byte-for-byte; missing/duplicate samples warn as upstream does.
  - **Custom expression `TAG[:Number]=[int|integer|float](EXPR)`.** A
    self-contained evaluator (`native_plugin_filltags_expr.go`) ports the slice
    of filter.c that fill-tags exercises: INFO/FORMAT tag references, the
    aggregations `SUM`/`AVG`|`MEAN`/`MAX`/`MIN`/`MEDIAN`/`STDEV` and their
    per-sample `SMPL_*`/`sXXX` variants, arithmetic `+ - * /`, unary minus,
    `ABS`, `PHRED`, and the genotype reductions `F_MISSING`/`N_MISSING`/
    `F_PASS(COND)`/`N_PASS(COND)` (the condition compiled with the shared native
    filter engine; `-S` restricts the active-sample set). `int()`/`integer()`
    yields Integer with C `round()` (half-away) rounding, `float()`/bare yields
    Float; `:Number` sets a fixed count, else `Number=.`. The "Added by
    +fill-tags expression ..." header is reproduced verbatim (quotes escaped),
    one per population.
  - **`-l/--list-tags`.** Prints the exact upstream available-tag table to
    stderr and exits non-zero with no stdout (matching upstream's `error()`).
  Byte-validated in `TestNativePluginFillTagsPops` and
  `TestNativePluginFillTagsListTags` (fixtures
  `tools/bcftools/testdata/parity/filltags_pops.vcf`, `filltags_groups.txt`,
  `filltags_sites.vcf`), plus binary-free `TestUnitFillTags*` unit tests for the
  pure helpers. Scope note: the rarer filter.c statistical functions over GT
  index/subscript forms (`binom`/`fisher` on `FMT/AD`) are not part of this
  evaluator; an expression using an unsupported function returns a clear
  evaluation error rather than the former blanket "not supported".

- **Native plugin format / output-mode tails** — **done**. The last
  per-plugin format and output-container modes are closed, byte-validated vs the
  upstream binary 1.23.1 in `native_plugin_format_tails_oracle_test.go` with
  binary-free `TestUnit*` coverage of the pure helpers
  (`native_plugin_format_tails_unit_test.go`):
  - **`+ad-bias --clean-vcf`/`-c`** emits the VCF subset to only the ALT alleles
    whose Fisher p-value passes `-t` (and drops sites where no comparison
    passes). The allele subsetting is an in-tree text-model port of htslib's
    `bcf_remove_allele_set` (`remove_allele_set.go`): the ALT list, the
    `Number=A`/`R`/`G` INFO and FORMAT fields, and `FORMAT/GT` allele indices are
    all remapped/reindexed (removed alleles → missing, phasing preserved).
    **`+ad-bias -f`/`--format`** appends a `bcftools query`-style column
    (evaluated once per record via the shared `ParseFormatString`/`emitRecord`
    engine) to every `FT` report line, with the `[N-]User data:` header column;
    `-f`+`-c` errors as upstream. The `-c` short form takes no argument while the
    `--clean-vcf` long form consumes and ignores one (an upstream `getopt_long`
    quirk reproduced exactly).
  - **`+remove-overlaps -m 'min(QUAL)' --missing`** now supports both the scalar
    `0` default and the **`DP`** coverage heuristic (a missing-QUAL record is
    valued at `maxQUAL*INFO/DP/maxQUAL_DP` over its overlap window, in `float32`
    to match htslib), and **`-Ot`/`-Otz`** emit a plain / bgzip-framed
    `chr<TAB>pos` list instead of the VCF. The min(QUAL) resolution is a faithful
    port of the vcfbuf MARK_EXPR push/flush state machine (`minQualBuf`). One
    upstream ring-index off-by-one (a stale overlap mark leaking across windows
    at a deletion-window boundary) is corrected on port — see
    `docs/UPSTREAM_BUGS.md#bcftools-remove-overlaps-minqual-stale-mark`; the
    oracle fixture deliberately avoids that corner so the upstream comparison is
    self-consistent.
  - **`+tag2tag --LXX-to-XX`** (and the partial `--LPL-to-PL`, `--LAD-to-AD`)
    expand the localized FORMAT tags back to the standard `Number=G` `PL` and
    `Number=R` `AD`, porting `process_LXX`: per sample, `FORMAT/LAA` maps the
    localized indices to global allele indices, `dst[0]=src[0]` with the rest
    defaulting to `-d`/`--defaults` (`AD:.`/`PL:.` by default), and the LPL→PL
    expansion uses the `tmp_laa[j]*(tmp_laa[j]+1)/2 + tmp_laa[k]` genotype-index
    layout. `-r`/`--replace` drops the consumed source tags (LAA only when it is
    the last remaining one) and `-s`/`--skip-nalt` skips sites above the
    threshold (suppressing header removal). The reverse `--XX-to-LXX` direction
    is a `todo` upstream too and is rejected with the same restriction.
  - **`+guess-ploidy -g`/`--genome`** expands the `b37`/`b38`/`hg19`/`hg38`
    non-PAR-chrX presets (`X:2699521-154931043`, `X:2781480-155701381`, and the
    `chr`-prefixed variants) into the equivalent `-r CHR:BEG-END` before the
    shared host region filter runs (a new `argvRewriter` hook in
    `native_plugin.go`), matching upstream's `case 'g'` region shortcut.
  - **`+af-dist -p`/`-d`** bin lists may be read from a file (one boundary per
    line, gzip-transparent) in addition to the inline comma-separated form,
    matching `bin_init`/`hts_readlist`'s "a comma indicates a list, otherwise a
    file" decision.
  - **`+gvcfz` and `+frameshifts`** — **done** (previously registration stubs
    that errored cleanly from `Init`, deferred for needing the FORMAT/GT filter
    engine and the synced-reader region cursor respectively; both now exist
    in-tree). **`+gvcfz`** ports the gVCF block state machine: the `-g
    FILTER:EXPR[;…]` clauses each compile to a full bcftools filter expression
    (now with the FORMAT/GT predicates upstream's `filter_init`/`filter_test`
    cover; the filter lexer accepts single `&`/`|` alongside `&&`/`||`), the
    first matching group selects each record's block, consecutive same-group
    reference blocks (ALT `<NON_REF>`/`<*>`) merge into the first record with
    INFO/END extended to the block end, FORMAT/DP set to the min MIN_DP/DP,
    FORMAT/GQ\|RGQ to the min, and FORMAT/PL to the element-wise min; a real
    variant flushes the block; non-PASS group labels add a `##FILTER` line whose
    Description is the verbatim `-g` string (`"`→`'`). `-i/-e` pre-filter and
    `-a` are supported; `-o/-O/-W` are handled by the host pipeline.
    **`+frameshifts`** annotates INFO/OOF over an exon BED / region-list,
    porting htslib's `bcf_sr_regions_overlap` monotonic forward cursor
    (`exonCursor`: per-chromosome sorted+merged 0-based regions, re-seek on a
    backwards or new-chromosome query) plus the per-allele
    `bcf_set_variant_type` classification (the signed `var->n` and the
    `VCF_INDEL` bit). The default reproduces upstream's **dead-code** behaviour
    byte-for-byte — the per-allele guard `var[i].type != VCF_INDEL` is always
    true under modern htslib (`VCF_INDEL|VCF_INS`/`|VCF_DEL`), so every
    exon-overlapping indel allele gets `OOF=-1` and the mod-3 result is never
    produced (see `docs/UPSTREAM_BUGS.md#bcftools-frameshifts-oof-dead-code`).
    The corrected exon-trim + length-mod-3 computation is implemented as a pure
    helper and exposed via the opt-in `--fix-oof` flag. Both are byte-validated
    vs 1.23.1 across the multi-group/RGQ/-i/-e gvcfz cases and the
    plain-BED/tabixed-BED frameshifts cursor paths, with binary-free `TestUnit*`
    coverage of the grouping, the variant classification, the exon cursor, and
    the OOF computation.

- **Native `+fixref` `id`/`--use-id` mode** — **done** (previously rejected
  from Init). `native_plugin_fixref_id.go` ports MODE_USE_ID: the REF allele is
  determined from a separate dbSNP VCF/BCF keyed by the **ID (rsID) column**
  rather than from strand convention. Upstream consults the dbSNP file through a
  per-input-chromosome region-restricted synced reader
  (`bcf_sr_set_regions(sr, chr, 0)`) and builds an rsID→{pos,ref} hash map
  rebuilt on every chromosome change (skip non-SNPs / non-[ACGT] REF / missing
  `.` IDs; first-wins on duplicate IDs). We reproduce that map exactly by
  streaming the dbSNP VCF once per chromosome through `iohelper.OpenReader`
  (transparent BGZF/gzip), keeping only same-chromosome records. The orientation
  decision mirrors `dbsnp_check`: input REF already equals the dbSNP REF →
  `none`; input ALT equals the dbSNP REF → `swap` REF/ALT **and** flip every
  sample GT (0↔1, phase preserved); a missing/unknown ID or neither-allele match
  → unresolved (annotated `skip`, dropped under `-d/--discard`). The
  position-correction path (move `rec->pos`, re-fetch the forward REF, count a
  `fixed pos`, fatal on a dbSNP-vs-FASTA REF mismatch) and the one-shot
  "corrected position(s) results in unsorted VCF" warning are reproduced.
  Byte-validated vs 1.23.1 on BOTH the corrected VCF (stdout) and the stderr
  stats summary in `native_plugin_fixref_id_oracle_test.go` (match / swap+GT /
  unknown-ID / ID="." / `-d` discard / `-t` custom tag), with binary-free
  `TestUnitFixref*` coverage of the map builder and the orientation decision.
  **Fix-on-port:** upstream's synced reader refuses a plain (un-bgzipped /
  un-indexed) dbSNP VCF; our streaming reader accepts it too (a one-directional
  superset that never changes output on upstream-accepted inputs) — see
  `docs/UPSTREAM_BUGS.md#bcftools-fixref-id-plain-vcf`.

- **Native `+vrfs` (variant read frequency score)** — **done** (previously
  rejected from Init as "requires mpileup2 + BAM/CRAM + regidx"). The native
  port (`native_plugin_vrfs.go`) reads the alignment list, the sites file and
  the FASTA reference, piles up each BAM/CRAM (via `pkg/htsgo` `alnio`/`sam`/
  `fasta`), counts per-sample ref/alt supporting reads at each indexed site,
  bins the VAF (`nn2bin`), and emits the `SITE`/`MEAN`/`VAR2` profile
  **byte-for-byte** vs 1.23.1. The key parity enabler is that vrfs runs the
  mpileup2 engine in `LEGACY_MODE 1`, whose realignment step is a **stub** in
  upstream (`mpileup2/mpileup.c` `legacy_mplp_func` has `// realign` commented
  out): there is **no BAQ and no base-quality adjustment**, and the vrfs count
  loop never consults base quality, so the `MAX_BQ`/`DELTA_BQ`/`MIN_REALN_*`/
  `MAX_DP_PER_SAMPLE` knobs vrfs sets are no-ops for the count. The pileup
  therefore reduces to read-level flag filtering (drop unmapped/secondary/
  qcfail/dup; **no** MAPQ floor, **no** supplementary drop, **no** orphan/proper-
  pair filtering) plus a direct CIGAR walk, which the classic engine
  reproduces exactly — including htslib `bam_pileup1_t` `is_del`/`indel`
  semantics (a deleted column reads the post-deletion base; the last M base
  before an I/D op classifies as a generic indel). Supported: streaming and
  `-i/--use-index` modes (identical output), `-d/--min-depth`, `-n/--nbins`
  with hard-coded-profile rescaling, `-r/--recalc hc|data|file:PATH`,
  `-b/--batch I/N` and `k=N`, `-m/--merge-batches`/`-M/--merge-files`, and
  `-o`/`-O t|z[0-9]` output. The empty-profile `MEAN` line emits `-nan`
  (glibc `printf` of `0.0/0`), matched by a NaN special-case in the formatter.
  Byte-validated in `native_plugin_vrfs_oracle_test.go` (10 profiling cases +
  2 merge cases over upstream-samtools-built BAM fixtures); binary-free
  `TestUnitVrfs*` cover every pure helper.

Note on vendored reference source: `reference_code/bcftools` and
`reference_code/htslib` are now both vendored as submodules. Earlier
roadmap text in this section was written when bcftools internals were
unavailable and called several HMM-based features unportable for that
reason. That is no longer true — `vcfcnv.c` (CNV HMM), `vcfroh.c` (RoH
HMM), `polysomy.c` + `peakfit.c`, `HMM.c` (the shared Viterbi/Baum-Welch
core), `bam2bcf*.c` (mpileup MAQ model + indel caller) and
`reference_code/htslib/errmod.c` are all available as porting
references. The deferrals below are scope/effort calls, not
source-availability blockers; each is annotated accordingly.

Genuine algorithmic gaps (subcommands present but running a v1
heuristic in place of the upstream algorithm — full detail in the
per-subcommand option-tail sections below):

- **`cnv`** — full port: the upstream copy-number HMM (`vcfcnv.c` +
  `HMM.c`, both vendored). 4-state single-sample / 16-state paired
  Viterbi + forward-backward over BAF+LRR Gaussian emissions, with
  `--optimize` cell-fraction estimation and `--baum-welch` transition
  re-estimation. See the `cnv` option-tail section for the validation
  situation (no upstream golden exists).
- **`roh`** — full port: 2-state Viterbi + forward-backward HMM
  (`vcfroh.c` + `HMM.c`) with physical-distance- and genetic-map-scaled
  transitions, allele-frequency estimation and Baum-Welch
  Viterbi training, plus `-G/--GTs-only` hard-GT emission and `-O z`
  BGZF output. Validated byte-for-byte against the upstream
  `roh.1.*.out` goldens and the live binary. The only remaining gap is
  PL-based emission (`-G` hard-GT is the supported path).
- **`polysomy`** — DONE: faithful port of the upstream
  Gaussian-mixture peak fit (`polysomy.c` + `peakfit.c`). The GSL
  Levenberg-Marquardt solver is ported in-tree as pure Go
  (`peakfit_lm.go`); no third-party dependency was added. All
  algorithm knobs are live. See the per-subcommand section below for
  the validation situation (no upstream golden exists).
- **`gtcheck`** — full discordance engine ported (`vcfgtcheck.c`):
  the dosage-bitmask concordance model, the probability discordance
  score (default `-E 40`) and the integer-mismatch score (`-E 0`), the
  per-site `-log P(HWE)` column (allele frequency from INFO/AC,AN or
  counted from FORMAT/GT), GT **and** PL scoring with upstream's
  header-driven auto tag selection (query prefers PL, panel prefers
  GT), `-u GT,PL`/`-u PL,GT` mixing, cross-check (lower-triangle) and
  `-g` panel modes, explicit `-p/-P` pairs (sorted by sample index),
  `--n-matches` top-N trimming (incl. the cross-check `i>ism` half),
  the monoallelic-site skip and `--keep-refs`, and `--distinctive-sites`
  with the `hts_srand48(0)`/`lrand48` tie-break reproduced exactly. All
  scoring/output paths are verified byte-for-byte against the live
  upstream binary in `gtcheck_parity_test.go` (modulo the
  non-reproducible provenance header and timing line, which are
  stripped). Remaining deferrals: the `-c/--cluster` dendrogram (which
  upstream itself rejects as "to be implemented"), `-O z` compressed
  output, and filter expressions (`-i/-e gt:/qry:`).
- **`mpileup`** — the upstream MAQ genotype-likelihood model is fully
  ported for the SNP path (`errmod.c` → `errmod.go`; `bam2bcf.c`
  glfgen/combine/2bcf → `bam2bcf.go`): the multi-allelic PL grid, the
  `<*>` unseen allele, BCF output, BAQ realignment (`pkg/htsgo/baq`),
  the per-site bias annotations (VDB/SGB/RPBZ/MQBZ/BQBZ/MQSBZ/SCBZ),
  MQ0F and `MPLP_SMART_OVERLAPS` read-pair quality merging are all wired
  in and verified byte-for-byte against the upstream goldens (slices
  1-4 done). Only indel calling (`bam2bcf_indel.c`) remains deferred.
- **`call`** — consensus and biallelic multi-allelic calling are
  implemented; the full upstream multi-allelic `-m` grid over >2 ALTs
  pairs with the mpileup MAQ port.
- **`csq`** — haplotype-aware consequence engine (slices 1-4 done):
  the `hap_node_t` tree, `cds_translate`, compound consequences,
  `@pos` reference pointers and `-p/-n` are ported, and the
  INFO/BCSQ output matches upstream byte-for-byte on the targeted
  goldens. Slice 4 added the GFF/output tail: `FORMAT/TBCSQ`
  per-haplotype text expansion (`query -f'[%TBCSQ\n]'`,
  `expandTBCSQ`), `--unify-chr-names` contig reconciliation,
  `--dump-gff FILE` model dumping (byte-exact vs upstream `gff_dump`),
  and non-text `-O b|u|z` output via the in-tree BCF/BGZF writers —
  all validated byte-for-byte against the live upstream binary in
  `csq_slice4_test.go`. `-l/--local-csq` (per-record,
  non-haplotype-aware `test_cds_local`) is now ported in `csq_local.go`
  and validated byte-for-byte (INFO/BCSQ) against the live upstream
  binary; see the "csq full-parity slicing plan" below.

Option-tail status on `gtcheck`:

- `-u PL` / `-u GT,PL` / `-u PL,GT` — **implemented.** PL→dosage and
  PL→probability (`pl_to_dsg`/`pl_to_prob`) ported; tag selection is
  header-driven (query prefers PL, panel prefers GT) with the
  per-record `set_data` fallback.
- `-E`-weighted probability scoring — **implemented** and is the
  default (`-E 40`). `-E 0` selects the integer-mismatch path. The
  `[5]Average -log P(HWE)` column carries the real per-site HWE
  estimate (AF from INFO/AC,AN or counted from FORMAT/GT).
- `--n-matches INT` (incl. negative = sort-by-HWE) — **implemented**,
  including the cross-check half-triangle `i>ism` emission rule.
- `--distinctive-sites NUM[,MEM[,TMP]]` — **implemented.** Greedy block
  selection with the `hts_srand48(0)`/`lrand48` tie-break reproduced in
  pure Go. The `,MEM,TMP` external-sort suffix is parsed and ignored
  (the port sorts in memory). NOTE: upstream builds a `# DS` comment
  header but never writes it; we reproduce that quirk and emit only the
  `DS` data rows.
- `--keep-refs` and the monoallelic-site skip — **implemented.**
- `--cluster N,N` — deferred (upstream itself errors "to be
  implemented"); accepted-and-rejected with a PARITY_ROADMAP pointer.
- `-O z` — bgzip output; only tab-text (`-O t`) is emitted.
- Filter expressions (`-i/-e` with `gt:`/`qry:` prefixes) are parsed but
  not yet applied.
- Index-backed `-r/-R` seek (post-filter only).
- Multi-allelic input is rejected (matches upstream's
  `bcftools norm -m -` requirement).
- Validation: `gtcheck_parity_test.go` builds/locates the live upstream
  binary and asserts byte-for-byte equality across every mode above.

Option-tail gaps on `roh`:

- **The HMM is now the full upstream port.** `vcfroh.c` + `HMM.c` are
  ported in-tree (`tools/bcftools/pkg/bcftools/hmm.go`,
  `roh.go`, `roh_genmap.go`): a 2-state Viterbi decode plus a
  forward-backward posterior, transition probabilities scaled by the
  physical (and, with `-m`/`-M`, the genetic-map) inter-marker
  distance, allele-frequency estimation (`-e/--estimate-AF` from
  GT cohorts), Baum-Welch parameter re-estimation
  (`-V/--viterbi-training`) and overlapping-window buffering
  (`-b/--buffer-size`). RG/ST quality scores are forward-backward
  phred scores and match upstream byte-for-byte on the
  `roh.1.*.out` goldens.
- `-O z` — **DONE.** The text output is wrapped in a BGZF writer
  (`bgzf.NewWriter`), mirroring upstream's `bgzf_open(.., "wg")`; the
  framed bytes decode byte-identically to upstream. Bare `-O z` (no
  `s`/`r`) defaults the sections to `s`+`r` first, matching
  vcfroh.c:1246. Validated against the live binary in
  `TestSubfeatRohOzParity` (both `-Osrz` and bare `-Oz`).
- `-G/--GTs-only FLOAT` hard-GT emission — **DONE** (the supported
  emission path; the upstream test corpus uses `-G30`). PL-based
  emission scoring is the remaining gap: `-e PL,...` falls back to
  rejecting sites with no usable genotype rather than reading PLs.

Option-tail gaps on the wave-1 additions (PR #86):

- `annotate`: the deferred option-tail is now implemented and validated
  byte-for-byte against the live upstream binary (see
  `tools/bcftools/pkg/bcftools/annotate_advanced_test.go`):
  - `--set-id [+]<FORMAT>` macro expansion — the non-FORMAT subset of the
    `bcftools query` macro language (`%CHROM`, `%POS`, `%POS0`, `%END`,
    `%END0`, `%ID`, `%REF`, `%ALT`, `%FIRST_ALT`, `%QUAL`, `%FILTER`,
    `%TYPE`, `%INFO/TAG` and bare `%TAG`), `\t`/`\n`/`\<c>` escapes, and the
    leading-`+` "only if ID empty" prefix. `%TYPE` ports htslib's
    `bcf_set_variant_type` (case-sensitive single-base SNP/REF, MNP, INDEL,
    OTHER, BND, OVERLAP).
  - `--merge-logic <tag:logic>` for range (BEG/END, aka FROM/TO) tables —
    `first` (default), `append`, `append-missing`, `unique`, `sum`, `avg`,
    `min`, `max`. Integer-typed `avg` truncates like upstream. Repeated
    `--merge-logic` flags accumulate (comma-joined) as upstream does.
  - `--min-overlap <ann:vcf>` reciprocal-overlap thresholds.
  - `--pair-logic <exact|some|all|any|snps|indels|both|id>` for VCF sources,
    ported from htslib's `bcf_sr_sort` pairing-score table.
  - `--single-overlaps` (apply only the first overlapping row).
  - `--rename-annots <file>` (rename INFO/FORMAT/FILTER tags in the header
    and per-record).
  Still deferred: `-c CHROM,FROM,TO` BED-style annotation against arbitrary
  INFO array setters with `-i/-e` filter expressions on the annotation rows,
  and `--mark-sites`.
- `isec`: `--collapse some` (REF match + any-ALT-in-common) is approximated
  via strict tuple match; deeper semantics deferred.
- `merge`: pre-sort assumption is enforced; no automatic CHROM/POS sort.
- `reheader`: in-place rewrite (`-i`) currently emits to stdout — caller
  is responsible for the swap.
- `sort`: `-m/--max-mem` and `-T/--tmpdir` are accepted but always
  in-memory.

Option-tail progress on `stats`:

- `-u/--user-tstv TAG[:min:max:n]` — **implemented**. Collects 1st-ALT
  Ts/Tv counts stratified by a numeric INFO tag, binned into `n` buckets
  over `[min,max]` (defaults `0:1:100`), including the `TAG[idx]`
  multi-value index form. The `USR:TAG/idx` section is byte-for-byte
  identical to upstream `bcftools stats -u` for both Float (`%e` bin
  labels) and Integer (`%.0f` bin labels) tags — validated against the
  C binary built from `reference_code/bcftools` and committed as goldens
  under `tools/bcftools/testdata/parity/stats/`.
- Core sections **SN, AF, QUAL, IDD, ST, DP, PSC, PSI, HWE** are now
  byte-for-byte identical to upstream. AF is rebuilt from AC/AN (or
  `--af-tag`) with the singleton bin folded into bin 1 and the dynamic
  `(i-1)/(mAF-1)` labels; QUAL prints the upstream `0.1*(iqual-1)` value
  with one decimal; IDD emits `.` for the unset mean-VAF column; ST walks
  the `ref<<2|alt` codes; DP uses the idist `<min`/`>max` edge labels and
  per-genotype counts under `-s/-S`; PSC counts ref/het/hom on SNP-typed
  (or hom-ref) genotypes with `%.1f` mean depth, singletons, hap and
  missing columns; HWE emits the het-fraction quartile distribution per AF
  bin (only with `-s/-S`). Validated by the un-skipped
  `TestParityStats_*` cases against `reference_code/bcftools` goldens in
  `tools/bcftools/testdata/parity/`.

Option-tail gaps on the convert/mendelian PR:

- `convert`: v1 covers the pass-through round-trip
  (VCF↔BCF↔VCF.gz) with sample/region filtering and -i/-e expressions,
  **plus the full Oxford GEN/sample family**:
  - `-g/--gensample PREFIX|GEN,SAMPLE` — VCF/BCF → `.gen`(.gz)+`.samples`.
  - `-G/--gensample2vcf PREFIX|GEN,SAMPLE` — `.gen`+`.sample` → VCF/BCF.
  - `--tag {GT|PL|GP}` — FORMAT tag driving the `.gen` genotype
    probability triples (default `GT`); mirrors upstream's
    `process_gt_to_prob3` / `process_pl_to_prob3` / `process_gp_to_prob3`.
  - `--3N6` — the 3*N+6 column layout (leading bare-CHROM column).
  - `--sex FILE` — adds the `sex` column to the `.sample` file
    ("ID\\t[MF]"); every sample must have an entry (matches
    `init_sample2sex`).
  - `--vcf-ids` — VCF IDs in the `.gen` ID column. Import also handles
    the IMPUTE2 reference-panel case where the CHROM:POS_REF_ALT label
    sits in the second column and `--` in the first.
  - `--chrom` — deprecated upstream; v1 errors with the same
    "please use --3N6 instead" message and a non-zero exit.

  Validated byte-for-byte against the live upstream `bcftools` binary
  (built on demand in `convert_gen_test.go`) for both directions across
  GT/PL/GP, `--3N6`, `--vcf-ids`, and `--sex`.

  plus the **TSV→VCF** and **gVCF→VCF** import shapes:
  - `--tsv2vcf FILE` with `-c/--columns`, `-f/--fasta-ref`,
    `-s/-S` (the AA genotype path), and `--keep-duplicates` (a no-op
    here, matching upstream's tsv_to_vcf which never consults it).
    The column setters CHROM/POS/ID/REF/ALT/AA mirror
    `reference_code/bcftools/tsv2vcf.c` and `vcfconvert.c`
    (`tsv_to_vcf`); REF and the `##contig` header are filled from the
    faidx-indexed reference. Live-upstream byte-for-byte parity tests
    cover the REF/ALT path, the AA genotype path (incl. indel-skip and
    all-missing rows), and `--keep-duplicates`.
  - `--gvcf2vcf` (with `-f/--fasta-ref`) expands reference blocks
    (symbolic ALT `<*>`/`<X>`/`<NON_REF>` or REF-only plus `INFO/END`)
    into one per-site record, fetching each REF base from the
    reference and dropping `INFO/END`. The malformed-gVCF overlap
    clamp from `gvcf_next_line` is reproduced. `-i/-e` filters pass a
    failing record through verbatim, as upstream does. Live-upstream
    parity tests cover the basic expansion and the overlap clamp.

  **plus the full IMPUTE2 HAP/legend family**: `--hapsample`,
  `--hapsample2vcf`, `--haplegendsample`, `--haplegendsample2vcf` and
  the `--haploid2diploid` modifier (`tools/bcftools/pkg/bcftools/convert_hap.go`).
  The HAP exporters mirror `vcfconvert.c`'s `vcf_to_hapsample` /
  `vcf_to_haplegendsample` and `convert.c`'s `process_gt_to_hap[2]`
  byte-for-byte (BGZF `.hap.gz`/`.legend.gz` content, plain `.samples`,
  no-ALT / non-biallelic skip counters and the per-run summary line);
  the inverse `*2vcf` importers mirror `hapsample_to_vcf` /
  `haplegendsample_to_vcf` including the `rev_als` allele-orientation
  check and the synthetic END/GT/contig header. Validated against the
  live upstream binary in `convert_hap_test.go` (no goldens). The
  `--vcf-ids` modifier is still deferred, so the `--vcf-ids`-only hap
  ID format is not yet emitted.

  Still deferred (sibling follow-up): only `--gvcf` block-output
  pairing remains, hard-rejected in `checkConvertDeferred` with a
  roadmap pointer.
- `mendelian`: the v1 port detects Mendelian inconsistencies for
  PED-style trios (one or more `-t CHILD,FATHER,MOTHER` flags, or
  `-T trio-file`), emits `INFO/MERR` per record, and supports the
  five upstream modes `a|c|x|d|+`. The `+` mode is treated as a
  synonym for annotate (the upstream-only `##bcftools_PG=` provenance
  header is deferred). The `--rules FILE` ploidy specification is
  accepted but currently only the chrX heuristic is honoured; full
  per-contig ploidy override (PAR boundaries, mitochondrial
  haploidy) is a follow-up.

Option-tail gaps on `filter` (this PR, simple-mode):

- `--mask [^]REGION` and `-M/--mask-file [^]FILE` — **DONE.** The mask
  region set is loaded (`loadMaskFile` mirrors htslib `regidx_init`'s
  extension detection: `.bed`/`.bed.gz`/`.bed.bgz` → 0-based BED,
  otherwise 1-based `CHR,POS[,END]` tab), and a record's span is tested
  for overlap; a hit fails the filter so the record is tagged with the
  `-s/--soft-filter` name, mirroring vcffilter.c's `pass &= mpass`. The
  leading `^` negates the source. `--soft-filter` is required (matching
  vcffilter.c:656). When both `--mask` and `-M` are given the file mask
  wins (upstream stores both into one `mask_list`, last write wins).
  Validated byte-for-byte against the live upstream binary
  (`TestSubfeatFilterMaskParity`: BED, tab, negation, `--mask` string,
  `--mask-overlap 0`).
- `--mask-overlap 0|1|2` — **DONE** for modes 0 (POS only) and 1 (REF
  span, the default). Mode 2 (variant boundaries) is treated as 1; the
  spans coincide for non-symbolic records.
- `-W/--write-index[=FMT]` — accepted; v1 never auto-indexes outputs.
- `-v/--verbosity INT` — accepted; v1 ignores.
- `--regions-overlap` / `--targets-overlap` — accepted; v1 always
  uses POS-in-region semantics.
- `-g/--SnpGap INT:TYPE` — the `:TYPE` qualifier (indel|mnp|bnd|other|overlap)
  is parsed but always treated as "indel" in v1.
- BCF output (`-O b|u`) round-trips through the shared `pkg/htsgo/bcf`
  writer; CSI auto-indexing is the `-W` follow-up above.

Option-tail status on `query` format tokens:

- Bare `%INFO` (whole INFO column) and `%SAMPLE` (per-sample name inside
  a `[ ... ]` group) — **DONE.** `%INFO` reconstructs the original INFO
  column in record order via `vcf.Variant.InfoString()`, mirroring
  convert.c's `process_info` with a NULL key; `%SAMPLE` resolves the
  sample index to its header name. Char-typed FORMAT/INFO tags
  (`%FMT/<char-tag>` and bare char tags) already round-trip through the
  VCF string value. Validated byte-for-byte against the live binary in
  `TestSubfeatQueryTokensParity`.
- `%N_ALT` is **NOT** a `query -f` format token: upstream `convert.c`
  has no such tag and `bcftools query -f '%N_ALT'` errors with "no such
  tag defined in the VCF header: INFO/N_ALT". `N_ALT` exists only in the
  `-i/-e` filter-expression language (`filter.c:3384`), computed on the
  fly. There is therefore nothing to port here; adding a `%N_ALT` format
  token would diverge from upstream.

Status on `som` (self-organizing-map filtering) — **SHIPPED (upstream
write bug fixed)**:

- Upstream `bcftools som` (`vcfsom.c`) is a standalone train/classify
  tool. It is **broken in the vendored upstream**: `som_write_map`
  (`vcfsom.c:170`) checks `fwrite("SOMv1",5,1,fp)!=5`, but `fwrite`
  returns the element count (1), so the comparison is always true and
  `--train` calls `error()` and exits 255 after truncating the `.som`
  map to 5 bytes. Consequently `--classify` can never read a map and the
  subcommand is effectively dead upstream
  (`docs/UPSTREAM_BUGS.md#bcftools-som-write-map`).
- The Go port **registers `som` and fixes the write bug** so the
  `--train`→`--classify` pipeline works. `tools/bcftools/pkg/bcftools/som.go`
  ports the SOM math (BMU search, neighbourhood update, count
  normalisation, distance scoring) verbatim from `vcfsom.c` and writes a
  usable map. Three deliberate divergences from a byte-exact port: (a)
  the on-disk map format is our own clean, versioned binary format
  (magic `SOMGO1`) because upstream's `SOMv1` layout is unusable; (b) the
  SOM reads INFO annotations straight out of a VCF/BCF (the
  `-t/--training-annots` set, default `QUAL,MQ,MQ0F,BQB,MQB,RPB,SGB`,
  min/max-normalised) rather than a pre-extracted `annots.tab.gz`; and
  (c) weight initialisation uses Go's `math/rand` (deterministic per
  seed) instead of glibc's `random()`. No live oracle exists (upstream
  crashes), so the port is validated by train→classify and map-file
  round-trips plus hand-checkable BMU/update unit tests
  (`som_test.go`). Upstream's experimental `-f/--nfold` cross-validation
  and `-m/--merge` knobs are accepted-as-surface only (v1 trains a single
  map). See `docs/UPSTREAM_BUGS.md`.

Option-tail status on `mendelian2`:

- `--rules ASSEMBLY[?]` — predefined inheritance rules. **DONE.** The
  GRCh37 / GRCh38 tables (`rules_predefs[]` in `mendelian2.c`) are
  ported verbatim into `mendelian2_rules.go` and drive a per-site,
  per-sex (1X / 2X) ploidy + inheritance model that replaces the old
  chrX-only heuristic: PAR regions stay diploid, the male-specific X
  is haploid maternal, Y is haploid paternal for males / absent for
  females, MT is haploid maternal. `list` / `list?` print the
  catalogue and `GRCh38?` prints that assembly's detailed table, each
  exiting 255 like upstream. Live-binary parity: `-m c` count output
  is byte-for-byte identical to `bcftools +mendelian2 --rules ...`
  across GRCh37/GRCh38 × male/female (see
  `TestLiveParity_RulesCountMode`).
- `--rules-file FILE` — custom `SEX_ID CHROM:BEG-END INHERITED_FROM`
  rules file. **DONE.** Same parser as the built-in tables; count
  output matches the live upstream binary
  (`TestLiveParity_RulesFileCountMode`).
- `-W/--write-index[=csi|tbi]` — auto-index output. **DONE** for
  bgzipped VCF output (`-Oz -o FILE`): the in-tree CSI/TBI writer
  (`BuildIndex`) is invoked after the output is flushed, and the live
  upstream `bcftools` reads the resulting index for region queries
  (`TestLiveParity_WriteIndexEmitsValidIndex`). Output is emitted as
  BGZF (not plain gzip) when `-W` is requested so it is indexable.
  Indexing a non-bgzipped / stdout output is rejected. (Go's `flag`
  parser stops at the first positional, so flags must precede the
  input file — a pre-existing limitation shared by every bcftools
  subcommand here.)
- `-m a` VCF output — **DONE.** Annotates the full upstream per-site
  INFO quartet `MERR` (trios with a Mendelian error), `MGOOD`
  (evaluable + consistent trios), `MMISS` (trios with missing/unusable
  genotypes) and `MNORULE` (trios with no applicable inheritance rule),
  with the verbatim `##INFO` definitions and values — byte-validated
  vs upstream 1.23.1 (`TestNativePluginMendelian2Plugin`). `-m d` sets
  offending GTs to `./.` regardless of the site ploidy. Count-mode
  parity is exact.
- `--regions-overlap 0|1|2`, `--targets-overlap 0|1|2` — accepted;
  v1 always uses POS-in-region semantics.
- `-v/--verbosity INT`, `--no-version` — accepted; v1 ignores both.
- `-r/-R/-t/-T` — region / target post-filter is wired through the
  CLI but the BCF synced-reader region-jump path is not used in v1;
  filtering happens after the records are read.
- `sites_not_diploid` counter never goes up in v1 because our
  `vcf.Variant` decoder coerces non-diploid GTs to missing rather
  than tracking ploidy. Tracked alongside the broader BCF FORMAT
  reconstruction in `docs/UPSTREAM_BUGS.md`.
- The PED-row sort order is by child name (deterministic);
  upstream sorts by min(sample-index) for sequential VCF reads. The
  sort is a performance optimisation only — the set of reported
  trios and per-trio counters are identical.
- `-i/-e EXPR` are accepted at the library boundary but only
  applied as record-level filters in v1 (no per-sample mask). The
  `sites_fail` counter is therefore always 0 in v1.

Option-tail status on `polysomy` (Gaussian-mixture peak fit — DONE):

- **The algorithm is the upstream Gaussian-mixture peak fit.**
  `polysomy.c` is ported faithfully: each chromosome's BAF values are
  binned into an `--nbins`-bin histogram, smoothed, and the RR/RA/AA
  bands are isolated and per-segment normalised (`init_dist`). Three
  candidate fits — CN2 (one bounded Gaussian near 0.5), CN3 (two
  symmetric Gaussians near 1/3 and 2/3) and CN4 (a central peak plus
  two symmetric side peaks) — are run over the heterozygous band, and
  the lowest CN that passes `--fit-th` plus the symmetry / peak-size
  checks is chosen, with `--cn-penalty` as the tiebreaker
  (`fit_curves`).
- **The peak-fitting engine `peakfit.c` is ported in-tree as pure
  Go** (`peakfit.go`): the Gaussian / centre-bounded-Gaussian /
  exponential peak models, the residual objective `(model-y)/0.01`,
  the unscaled `Σ|model-y|` goodness measure, and the Monte-Carlo
  restart driver. The GSL `gsl_multifit_fdfsolver_lmsder` non-linear
  least-squares solver is replaced by an in-tree pure-Go
  Levenberg-Marquardt solver (`peakfit_lm.go`) — the normal-equations
  damping loop `(JᵀJ + λ·diag)·δ = -Jᵀr` with an analytic Jacobian, a
  λ up/down schedule, and convergence on the parameter / gradient
  deltas at tolerance 1e-8. `peakfit_lm.go` also ports glibc's
  `random()` so the Monte-Carlo restart stream after `srand(0)`
  matches upstream. **No third-party dependency was added** — CLAUDE.md
  scopes the one sanctioned numerical dep (gonum) narrowly and
  explicitly excludes stats-fitting tools like this.
- **All algorithm knobs are live:** `-b/--peak-size`,
  `-c/--cn-penalty`, `-f/--fit-th`, `-i/--include-aa`,
  `-m/--min-fraction`, `-p/--peak-symmetry`, `-n/--nbins`,
  `-S/--smooth`, `--ra-rr-scaling`, `--force-cn`.
- **Validation — no upstream golden exists.** `bcftools polysomy`
  produces only per-chromosome PNG plots and a `dist.dat` dump under
  `--output-dir`; `reference_code/bcftools/test/test.pl` has no
  `polysomy` invocation. Byte-for-byte parity therefore cannot be
  demonstrated and is **not claimed**. The port is validated instead
  with: (a) unit tests of the LM solver against analytic curves with
  known optima (a linear-in-parameters quadratic and a non-linear
  single Gaussian); (b) a check of the glibc `rand()` port against the
  published `srand(0)` reference sequence; (c) unit tests of each peak
  model (Gaussian / bounded-Gaussian / exp) recovering known
  parameters from clean synthetic samples; (d) hand-constructed BAF
  distributions for the canonical karyotypes — a clean diploid
  (single peak at 0.5 → CN2) and a clear trisomy (two peaks at 1/3 and
  2/3 → CN3), exercised both directly and end-to-end through the VCF
  reader. See `polysomy_test.go` and `peakfit_test.go`.
- `-o/--output-dir PATH` — accepted but ignored; the port writes the
  per-chromosome TSV to stdout (no PNG plots, no `dist.dat`).
- `--regions-overlap`, `--targets-overlap` — accepted; the port uses a
  chromosome-name post-filter (per-base interval filtering deferred).
- `-v/--verbosity` / `--verbose` — accepted; the port ignores it (no
  per-iteration fit trace).
- BAF source: upstream requires FORMAT/BAF; the port also accepts
  FORMAT/AD = REF,ALT as a fallback (synthesises BAF as
  `ALT / (REF + ALT)` at het sites).
- Per-record `-i/-e` are NOT in upstream `polysomy.c:main_polysomy`
  and we follow upstream's surface exactly (no invented flags).

Option-tail gaps on `cnv` (full HMM port):

- **The algorithm is the upstream HMM.** `vcfcnv.c` is ported
  faithfully: a 4-state HMM (CN0/CN1/CN2/CN3) for a single sample, or
  a 16-state HMM for the paired tumour/control mode (`-c`), swept per
  contig with Viterbi + forward-backward. Emission probabilities are
  the upstream joint BAF + LRR model — a truncated-Gaussian BAF peak
  mixture weighted by genotype frequencies fRR/fRA/fAA, combined with
  a per-state LRR Gaussian. The generic engine is the shared
  `hmm.go` port (also used by `roh`), reused unchanged. Every HMM
  tuning knob is now load-bearing: `-a/--aberrant` (CN3 BAF peak
  shift), `-b/--BAF-weight`, `-e/--err-prob`, `-l/--LRR-weight`,
  `-L/--LRR-smooth-win`, `-d/--BAF-dev`, `-k/--LRR-dev`,
  `-x/--xy-prob`, `-P/--same-prob`, `-O/--optimize` (iterated
  forward-backward cell-fraction estimation), and `-W/--baum-welch`
  (per-contig transition re-estimation).
- **Validation: no upstream golden exists.** bcftools `test/test.pl`
  contains no `cnv` invocation and `test/` ships no `cnv` fixtures
  (`cnv` output is plot-oriented). The Go port is therefore validated
  with hand-derived cases: a clean diploid region decodes to a single
  all-CN2 region; a long missing-BAF run in paired mode decodes to
  CN0; a het-band-split + positive-LRR run decodes to CN3; and
  unit tests pin the ported transition matrix (column-stochastic,
  the bad-xy-prob guard), the truncated-Gaussian `norm_cdf`, the
  emission model, the smoother and the initial-probability vector.
  Knob load-bearingness is asserted by tests that change `--err-prob`
  / `--xy-prob` and observe a different decode. Byte-for-byte parity
  against upstream is NOT claimed because upstream emits no
  comparable golden.
- `--AF-file` — **DONE.** Upstream's `vcfcnv.c` uses the AF file two
  ways and the port now mirrors both: (1) it acts as a targets filter
  (sites whose CHROM:POS is absent are dropped, `vcfcnv.c:27-31`,
  `:1429`), and (2) it recomputes the per-site genotype frequencies
  fRR/fRA/fAA under Hardy-Weinberg from each site's non-reference AF
  (`vcfcnv.c:735-739`) instead of the fixed defaults. The lookup
  mirrors `read_AF`: the AF is used only when the record's full allele
  vector (REF + all ALTs) matches a file entry, otherwise the listed
  site falls back to `nonref_af_dflt = 0.1` (`vcfcnv.c:1257`). When no
  AF file is given the port still uses the fixed defaults
  fRR/fRA/fAA = 0.76/0.14/0.098. Validated against the live binary in
  `TestSubfeatCNVAFFileParity` (RG region, quality and site/HET counts
  byte-identical) and a targets-filter / multiallelic-match unit test.
- `-o/--output-dir` — upstream writes per-sample / per-region plot
  data and several `.tab`/`.cn` files into this directory; this port
  streams a single summary TSV (the upstream `summary.tab` "RG"
  rows) to stdout regardless of the path (the flag is still required
  for CLI parity). The per-site `cn.<sample>.tab` and `dat.*.tab`
  files and the `CF` cell-fraction summary rows are not emitted.
- **Paired-mode control counts are dropped from `summary.tab`.** In
  paired (tumour/control) mode upstream writes per-sample summary
  files and a `summary.tab` whose "RG" rows carry four count columns
  — query `nSites`/`nHETs` and control `nSites`/`nHETs`
  (`vcfcnv.c:296-298,1095`). This port emits only the merged
  `summary.tab` "RG" rows with the query sample's `nSites`/`nHETs`;
  it computes the control counts internally but does not write them.
  This is a deliberate I/O-surface reduction, consistent with the
  single-stream `-o/--output-dir` simplification above.
- `-p/--plot-threshold` — accepted; this port emits no plots.
- `--regions-overlap` / `--targets-overlap` — accepted; always
  POS-in-region semantics.
- `-v/--verbosity` — accepted; ignored.
- Indel / non-SNP records — the BAF/LRR signals are typically
  per-marker SNP data; the port honours upstream's behaviour (treat
  each record as one marker regardless of REF/ALT).

Option-tail status on `csq` (slices 1-4 done, plus `-l/--local-csq`):

- **The engine IS haplotype-aware.** `bcftools csq` now phases
  variants per haplotype, walks the GFF transcripts, builds the
  `hap_node_t` tree and reports per-haplotype compound consequences
  in `INFO/BCSQ`. Indels, splice-site disruption, start/stop
  refinement and compound-het bookkeeping are all handled.
- `-p/--phase a|m|r|R|s` — load-bearing: the haplotype-construction
  modes are ported (`phaseRequire/Merge/AsIs/Skip/NonRef/DropGT`).
- `-i/--include` / `-e/--exclude` — accepted; not yet evaluated
  (slice 4). The expression evaluator already exists in
  `pkg/bcftools`; the wire-up is a small follow-up.
- `-s/--samples` / `-S/--samples-file` — accepted; the engine
  currently walks every header sample. Sample subsetting is slice 4.
- `-n/--ncsq INT` — parsed into the per-haplotype `FORMAT/BCSQ` cap
  (`ncsq2`), load-bearing for the bitmask emission and its
  `FORMAT/TBCSQ` text expansion (slice 4, `expandTBCSQ`).
- `-B/--trim-protein-seq INT` — **implemented**: amino-acid
  predictions in INFO/BCSQ are abbreviated to the first `INT`
  residues plus `..<index>` (ports `kprint_aa_prediction`).
- `-b/--brief-predictions` — **implemented**: upstream alias for
  `-B 1`. Validated byte-for-byte against the upstream binary
  (`tools/bcftools/testdata/parity/csq/`).
- `-C/--genetic-code INT|l` — **implemented** for the transcribed
  NCBI tables (`0, 1, 2, 3, 5`; `l` lists them). Codon translation
  uses the selected table; validated against the upstream binary.
  Additional tables can be added by appending to `gencodeTables`.
- `-l/--local-csq` — **implemented**: selects the per-record,
  non-haplotype-aware caller (`test_cds_local`, ported in
  `csq_local.go`). Each record's coding consequence is derived from its
  own ref/alt against the spliced reference, so compound consequences
  spanning several records are not joined (unlike the default
  haplotype-aware path). Validated byte-for-byte (INFO/BCSQ) against the
  live upstream binary (`TestCSQ_LocalUpstreamParity`).
- `--unify-chr-names 0|VCF,GFF,FAI` — **done (slice 4)**: the three
  comma-separated prefixes reconcile VCF/GFF/FASTA contig namespaces
  (`parseUnifyChrNames` / `unifyChrName`); `0` disables. Validated vs
  the live upstream binary.
- `--dump-gff FILE` — **done (slice 4)**: writes the parsed GFF model
  (genes/transcripts/CDS/UTR/exons) as a BGZF GFF3, byte-exact with
  upstream `gff_dump` on position-ordered inputs (`csq_dump.go`).
- `FORMAT/TBCSQ` — **done (slice 4)**: per-haplotype text expansion of
  the `FORMAT/BCSQ` bitmask (`query -f'[%TBCSQ\n]'`, `expandTBCSQ`),
  byte-for-byte vs upstream.
- `-O b|u|z|v|t` — **done**: VCF text (`v`), BGZF VCF (`z`) and BCF
  (`b`/`u`) via the in-tree writers (`openCSQOutput`); and the
  streaming tab-delimited text form (`t`, upstream `FT_TAB_TEXT`),
  which emits one `CSQ<TAB>sample<TAB>haplotype<TAB>chrom<TAB>pos<TAB>consequence`
  row per (sample, haplotype) consequence. It ports upstream's
  `text_stage`/`hap_stage_text`/`text_print_vcsq` path (the intron /
  non_coding tscript-level consequences are pushed to INFO/BCSQ but
  never text-staged, exactly as upstream, because their staged
  `vcf_ial` is 0), byte-validated vs 1.23.1 in `csq_text_test.go`
  (`TestCSQTextOracleParity`) plus binary-free `TestUnitCSQText*`
  coverage. The leading `#`-comment version/command provenance lines
  carry our build identity and are stripped before comparison.
  (Sample subsetting / `-s -` GT dropping is a separate open item.)
- `--threads`, `-v/--verbosity`, `-W/--write-index`, `--force`,
  `--no-version`, `-q/--quiet` — accepted; v1 ignores.
- The minimal GFF3 parser (`pkg/htsgo/gff`) understands `gene`,
  `mRNA` / `transcript`, `CDS`, and `exon` rows. Other feature
  types are silently skipped — fine for the v1 SNP classifier
  but the parser will need extension for splice-site / UTR work.

### csq full-parity slicing plan

Upstream `csq.c` is ~3994 lines. The Go port is sliced as follows.
Every upstream `csq` golden (`test/csq*.out`) is produced by the
haplotype engine and contains compound consequences (`103G>A+108T>A`),
reference pointers (`@107`), indels and frameshift on the *same* line,
so **no golden validates byte-for-byte until the haplotype engine
(slice 3) is complete**. Slices 1-2 were validated with hand-derived
unit tables; slice 3 unblocked the INFO/BCSQ goldens (now passing
byte-for-byte). The fixtures are vendored under
`tools/bcftools/testdata/csq/`.

- **Slice 1 — region classifier + SO-term completion (per-record, no
  haplotype tree). DONE.** The GFF3 model now carries `five_prime_UTR`
  / `three_prime_UTR` rows (explicit, or derived as exon-minus-CDS)
  and non-coding biotypes; per-record detection of `5_prime_utr`,
  `3_prime_utr`, `intron`, `non_coding` and the splice set
  (`splice_donor`, `splice_acceptor`, `splice_region`) is ported from
  `splice_init` + the SNP/MNP arm of `splice_csq` (the
  `N_SPLICE_DONOR=2` / `N_SPLICE_REGION_INTRON=8` /
  `N_SPLICE_REGION_EXON=3` boundary math, plus the 8bp exon-index
  padding from `gff.c`) and the `test_utr` / `test_splice` /
  `test_tscript` dispatch. The SO-term codon set is complete
  (`stop_gained`, `stop_lost`, `start_lost`, `stop_retained`,
  `coding_sequence`). Landed in `tools/bcftools/pkg/bcftools/
  csq_classify.go`.
- **Slice 2 — indel consequence classification (per-record). DONE.**
  `splice_csq_ins` / `splice_csq_del` / `splice_csq_mnp` /
  `splice_csq_complex` are ported: frameshift vs
  inframe-insertion/deletion, `feature_elongation` /
  `feature_truncation`, and indels at splice sites against a single
  transcript. Bundled with slice 1 in `csq_classify.go`.

  **PROVISIONAL: per-record indel frame bits.** `spliceCSQIns` /
  `spliceCSQDel` set the `csqFrameshift` vs `csqInframeIns` /
  `csqInframeDel` bit from the raw allele-length delta (`%3`).
  Upstream's `splice_csq_ins` / `splice_csq_del` do **not** set any
  frame bit at the splice layer — `hap_add_csq` recomputes frameshift
  vs inframe from the *translated* `dlen` once the haplotype is
  threaded through the spliced reference. The per-record bits are
  therefore an approximation that holds for a clean single-exon indel
  but is wrong whenever the indel spans an intron, partially overlaps
  the CDS boundary, or interacts with other variants on the same
  haplotype. The slice-3 engine MUST replace these per-record bits
  with the `hap_add_csq` `dlen`-based computation. To stop a
  CDS-internal indel being double-staged, `classifyTranscriptVariant`'s
  test_splice arm masks `spliceCSQNonSplice` (the provisional frame /
  elongation / truncation bits) off before deciding whether a splice
  hit occurred — see `csq_classify.go`.

Slices 1+2 shipped together as the per-record-classifier PR. The
`*`-upstream-stop prefix, `shifted_del_synonymous` start/stop refinement
and `inframe_altering` need the spliced reference / haplotype context
and are deferred to slice 3 with the engine. Validation is by
hand-derived unit tables (`csq_classify_test.go`): one case per SO
class — UTR5/UTR3, intron, splice donor/acceptor/region (fwd + rev
strand), stop_gained/stop_lost/start_lost, missense, synonymous,
inframe vs frameshift indel, splice-site indel — plus the kput_vcsq
SO-term precedence ordering.

- **Slice 3 — the haplotype-aware engine. DONE.** The haplotype tree
  (`hap_node_t`, `hap_init`, `hap_finalize`, `hap_add_csq`,
  `cds_translate`), the per-transcript padded + spliced reference
  build (`tscript_init_ref` / `tscript_splice_ref`), the full
  `splice_csq` family with `set_refalt` / `splice_build_hap` /
  `shifted_del_synonymous`, the `vbuf` / `pos2vbuf` position-clustered
  VCF buffer, `csq_push` / `csq_stage`, `kput_vcsq`, and the
  `-p/--phase {a|m|r|R|s}` haplotype-construction modes plus
  `-n/--ncsq` are ported in `csq_hap.go`, `csq_splice.go`,
  `csq_engine.go` and `csq_process.go`. This produces compound
  consequences (`103G>A+108T>A`), the `@pos` reference pointers, the
  `*`-upstream-stop prefix, and the true frameshift / inframe /
  elongation / truncation / start_retained / stop_retained calls from
  the translated `dlen`. The GFF3 CDS reading-frame phase is now
  trimmed off the 5' CDS exon at index-build time (mirroring `gff.c`),
  so the spliced CDS is frame-aligned for both the engine and the
  per-record classifier. *Passes byte-for-byte:* `csq.1.out`,
  `csq.oob-codon.out`, `csq.splice.issue-2543.1.out` — see
  `csq_golden_test.go::TestCSQGoldenINFO`.
- **Slice 4 — GFF/output tail. DONE.** The
  `FORMAT/TBCSQ` `bcftools query` expansion (`expandTBCSQ`, decoding
  the per-haplotype `FORMAT/BCSQ` bitmask into the `hap1\thap2`
  consequence list); `--unify-chr-names 0|VCF,GFF,FAI` (the 3-field
  rename spec, `parseUnifyChrNames` / `unifyChrName`); `--dump-gff
  FILE` (the parsed gene/transcript/CDS/UTR/exon model written as a
  BGZF GFF3, byte-exact with upstream `gff_dump`, `csq_dump.go`); and
  BCF / `-O b|u|z` output via the in-tree BCF/BGZF writers
  (`openCSQOutput`, sharing `bcftools view`'s writer). All four are
  validated **byte-for-byte against the live upstream binary** in
  `csq_slice4_test.go` (`TestCSQSlice4{DumpGFF,OutputFormats,TBCSQ,
  UnifyChrNames}`). The BCF writer gained unified-dictionary IDX
  dedup across INFO/FILTER/FORMAT (needed because `BCSQ` is both an
  INFO and a FORMAT tag) plus deterministic INFO emission ordering
  via `InfoOrder` (`pkg/htsgo/bcf`). **Still deferred:**
  the `-i/-e` filter wire-up. Also the `GF_NMD`/`NMD_transcript`
  branch of upstream `kput_vcsq`: `kputVcsq` currently omits
  NMD-transcript consequence emission.
  *Unblocks:* the `FORMAT/TBCSQ` and `--unify-chr-names` parity paths;
  the `csq.2.out` / `csq.3.out` / `csq.chr.out` goldens still require
  their (uncommitted) `.out` fixtures to run under
  `TestCSQGoldenINFO`.

Status: slices 1+2 (the per-record consequence logic), slice 3 (the
haplotype-aware engine) and slice 4 (the GFF/output tail) **DONE** —
see `csq_hap.go`, `csq_splice.go`, `csq_engine.go`, `csq_process.go`,
`csq_dump.go`, the `expandTBCSQ`/`unifyChrName` paths, and the tests
`csq_golden_test.go` + `csq_slice4_test.go`. The standalone v1
per-record classifier has been folded into the engine;
`csq_classify.go` now holds only the shared `csqStrings` table and the
`csq*` SO-term bit constants. The `bcftools csq` INFO/BCSQ output
matches upstream byte-for-byte on `csq.1.out`, `csq.oob-codon.out` and
`csq.splice.issue-2543.1.out`, and slice 4's `FORMAT/TBCSQ`,
`--dump-gff`, `--unify-chr-names` and `-O b|u|z` paths are validated
byte-for-byte against the live upstream binary. `-l/--local-csq`
(`test_cds_local`, the per-record non-haplotype-aware caller) is now
ported in `csq_local.go` and validated byte-for-byte (INFO/BCSQ)
against the live upstream binary (`TestCSQ_LocalUpstreamParity`).

Option-tail gaps on `mpileup` (SNP-only MAQ model; slices 1, 2 & 3 done):

- **The likelihood model IS the upstream MAQ model.** Slice 2 wired
  `bam2bcf.c::bcf_call_glfgen` / `bcf_call_combine` / `bcf_call2bcf`
  (ported in `bam2bcf.go`) onto the slice-1 errmod port. `mpileup`
  emits one BCF/VCF record per covered position with the `<*>`
  unseen allele, the full multi-allelic PL grid, and
  INFO/DP/I16/QS/MQ0F. The obsolete uniform-error binomial is gone.
- **BAQ recalibration IS wired (slice 3 done).** Reads are run
  through `pkg/htsgo/baq.SamProbRealn` (apply+extend mode) before their
  bases enter the pileup, gated by the ported `mplp_realn` /
  `MPLP_REALN_PARTIAL` column heuristic. `-B/--no-BAQ` disables it;
  `-E/--redo-BAQ` recomputes BAQ. Partial realignment is the default;
  `-D/--full-BAQ` clears `MPLP_REALN_PARTIAL` and forces full BAQ (every
  read realigned). The `<*>`-only `mpileup/*.out` golden records now
  byte-match (`TestMpileupBAQGoldens`).
- **No bias annotations (slice 4 TODO).** VDB / SGB / RPBZ / MQBZ /
  BQBZ / MQSBZ / SCBZ are not yet computed; records carry only
  INFO/DP/I16/QS/MQ0F. The `calc_vdb` / `calc_SegBias` /
  `calc_mwu_biasZ` machinery and the indel caller land in slice 4.
- **No indel calling.** The full upstream indel realigner
  (`bam2bcf_indel.c`) and the consensus indel mode
  (`bam2bcf_edlib.c`) are deferred. Every knob that drives the indel
  model — `-e/--ext-prob`, `-F/--gap-frac`, `-h/--tandem-qual`,
  `--indel-bias`, `--indel-size`, `-I/--skip-indels`,
  `-L/--max-idepth`, `-m/--min-ireads`, `-M/--max-read-len`,
  `--open-prob`, `--indels-cns`, `--indels-2.0`, `--no-indels-cns`,
  `--ar-prob`, `--ambig-reads / --ar`, `--del-bias`, `--poly-mqual`,
  `--no-poly-mqual`, `--score-vs-ref`, `--seqq-offset` — is accepted
  at the CLI for parity but inert in v1. The v1 emit path is
  equivalent to running upstream with `-I/--skip-indels` set.
- **The FORMAT/PL grid is multi-allelic.** Slice 2 emits the full
  upper-triangle `g[z++] = a[j]*5 + a[i]` grid of
  `n_alleles*(n_alleles+1)/2` values per sample, including the `<*>`
  unseen allele.
- **No BAI seek.** `-r/--regions` and `-R/--regions-file` are
  post-filters applied after a linear scan of every input BAM; the
  BAI-seek fast path lives in `pkg/htsgo/sam` but is not wired
  through `mpileup` in v1. Tracked as a follow-up — perf only, no
  output difference.
- **No per-read group filtering.** `-G/--read-groups` is parsed and
  stored; v1 includes every record whose @RG passes the standard
  filters. `-Z/--ignore-RG` is accepted but inert.
- **No gVCF blocking.** `-g/--gvcf` is accepted; one BCF/VCF record
  is emitted per covered reference position (gVCF range-blocking is a
  follow-up).
- **`-a/--annotate LIST` is parsed and partially honoured.**
  `parseFormatFlag` ports the upstream tag-list parser (mpileup.c:1141);
  tokens are accepted as bare names ("AD"), with FORMAT/INFO prefixes,
  and with the "-" prefix to clear bits. The upstream default
  (BQBZ/IDV/IMF/MQ0F/MQBZ/MQSBZ/RPBZ/SCBZ/SGB/VDB + FORMAT/AD) is the
  starting bitset; user tokens layer on top. Today FORMAT/AD, ADF, ADR
  and the bias INFO tags emit; the remaining `INFO/AD,ADF,ADR,SP,SCR`,
  `FORMAT/DP,DV,DPR,SP,SCR,QS,NMBZ,QM,DP4` tags land alongside the
  indel branch of 2bcf (slices 4e.3 and beyond).
- **`-O u|b` (BCF output) works.** Slice 2 wired BCF output through
  `pkg/htsgo/bcf` (`-O b` is BGZF-wrapped, `-O u` uncompressed);
  `-O v` (text VCF) is the default and `-O z` is gzipped VCF.
- **`--threads`, `-v/--verbosity`, `-W/--write-index`, `--no-version`,
  `-A/--count-orphans`, `-x/--ignore-overlaps`, `-d/--max-depth`,
  `-q/--min-MQ`, `-Q/--min-BQ`, `--max-bq`** — fully implemented in v1.
- `-X/--config STR` (presets like `1.12`, `2.1`, `ultima`,
  `pacbio-ccs-1.20`) — accepted; v1 ignores. Most presets toggle the
  indel-model knobs which are inert above.
- `-6/--illumina1.3+` — accepted; v1 ignores (input BAMs are
  Phred+33 across the board).
- `-C/--adjust-MQ INT` (MAPQ tail adjustment) — accepted; v1 ignores.
- The 5-flag mask flags (`--skip-any-unset`, `--skip-all-unset`,
  `--skip-any-set`, `--skip-all-set`, `--ls`) — accepted and stored;
  v1 honours only the standard `--ff` defaults (UNMAP, SECONDARY,
  QCFAIL, DUP, SUPPLEMENTARY) baked into `mpileupKeepRecord`.
- `--delta-BQ` — implemented (default 30): a base quality is capped at
  `neighbour_qual + delta` before the MAQ model sees it.
- `--seed` — accepted; ignored (no subsampling below the 255-read
  errmod cap).

Option-tail gaps on `consensus` (simple-mode):

- `--het-only` — **DONE (PR #220–#225 wave).** Implemented as a fix for
  an upstream dead-option bug (the upstream flag is parsed but never
  consulted); see `UPSTREAM_BUGS.md`. `--ignore-overlaps` landed earlier.
- `-c/--chain FILE` — **implemented.** Writes a UCSC-format liftover
  chain mapping reference to consensus coordinates alongside the
  consensus FASTA. The chain engine (`consensus_chain.go`) mirrors
  upstream's `init_chain` / `push_chain_gap` / `print_chain`, including
  the back-to-back gap merge and the 1-base block extension when REF
  and ALT share their leading base. Byte-for-byte parity against the
  live upstream binary is locked in by `consensus_chain_parity_test.go`
  (chain file AND consensus FASTA compared).
- `-H N` / `-H NpIu` / `-H I` — **implemented.** `-H N` selects the
  N-th haplotype slot of FORMAT/GT (resolved through the GT, matching
  upstream's `ialt = GT[haplotype-1]`, not a bare ALT index). `-H NpIu`
  applies the N-th haplotype for phased genotypes and an IUPAC
  ambiguity code for unphased ones; `-H I` applies IUPAC codes for all
  genotypes. The IUPAC encoder OR-s the per-position nucleotide
  bitmasks across the genotype's alleles (mirroring `iupac_set_allele`).
  Live-upstream parity is covered by `TestConsensusHaplotypeParity`.
- `--regions-overlap 0|1|2` — accepted; v1 ignores (no synced-reader
  region jump path).
- `-v/--verbosity INT` — accepted; v1 ignores.
- `-r/-R/-t/-T` — upstream `consensus` does NOT advertise these flags
  (only `--regions-overlap`); the v1 port follows upstream exactly.
- Multi-sample apply: `-s LIST` accepts a comma list but v1 honours only
  the first name (the other entries are silently ignored, mirroring
  upstream's single-sample focus when -H is unset).
- Complex MNP and SV (BND, breakend, `<DEL>` etc.) records are not yet
  applied to the reference. The current v1 covers SNPs and simple
  REF/ALT length-difference indels.
- Overlapping variants: first wins (left-to-right). Upstream emits a
  warning and folds them into a longer ALT; v1 doesn't.

Plus:

- **`bcftools view`** more flags: `--regions-overlap`, `--targets-overlap`,
  `--no-version`, `--write-index`, `--phased`.
  - **Implemented:** `-x/--private` and `-X/--exclude-private` (select /
    exclude sites whose non-reference alleles are exclusive to the sample
    subset). Mirrors upstream `vcfview.c`
    (`non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub`); applied after
    sample subsetting, gated on a `-s`/`-S` subset and a GT FORMAT field.
    Validated by table-driven unit tests plus a live upstream-parity test
    (`TestView_PrivateUpstreamParity`) that builds the upstream C binary
    from the vendored submodules and compares record selection in-process
    (`bcftools view -s S1,S2 -x/-X`); no committed golden snapshots. The
    input fixture lives under `tools/bcftools/testdata/parity/view/`.
    INFO/AC/AN recomputation after subsetting remains the separate
    documented gap (see the `view -s` note above), so the comparison
    blanks the INFO column.
- **CSI seek** for region queries: today we validate via the index then
  linear-scan. Real chunk-seek is the natural follow-up.

Subcommand-tail gaps on `bcftools call` (largely **CLOSED — salvaged
from PR #219**):

- **Full multi-allelic caller (`-m`) — DONE (salvaged from PR #219).**
  `callm.go` is a faithful port of `mcall.c`: EM allele-frequency
  estimation, per-site QUAL, the max-likelihood GT, and the INFO
  rewrite (AN/AC/DP4/MQ). The consensus caller (`-c`) is ported in
  `callc.go`. Live-oracle parity over a whole 4000+-site contig is
  asserted in `TestLiveCall` (`call -m`, `-v`, `-A`, BCF input,
  regions, and the items below).
- **`--ploidy GRCh37 / GRCh38` and `--ploidy-file` — DONE (salvaged).**
  `call_ploidy.go` builds the per-region, per-sex ploidy table from the
  predefined GRCh37/GRCh38 maps or a `--ploidy-file`; the default sex is
  F (matching `vcfcall.c` sample2sex init). Gated by
  `TestLiveCall/m_ploidy_grch37`, `/m_ploidy_grch38`, `/m_ploidy_file`.
- **`--gvcf` block-emit mode — DONE (salvaged).** `callm_gvcf.go` +
  `mpileup_gvcf.go` band consecutive REF-only records sharing a
  per-sample MIN_DP bin into a single `INFO/END`+`MIN_DP` record.
  Gated by `TestLiveCall/m_gvcf_0_5_10`, `/m_gvcf_5`.
- **`-C alleles` / `-T sites.tsv` constrain family — DONE (salvaged).**
  `call_constrain.go` loads the sites TSV and projects each record;
  `-C trio` mirrors upstream's own runtime "temporarily disabled"
  error. Gated by `TestLiveCall/m_C_alleles`,
  `/m_C_alleles_insert_missed`.
- **`-G` sample groups, `-V` skip-variants, `-*`/`-M` allele flags,
  `-F` prior-freqs — DONE (salvaged).** `call_groups.go` plus the
  CallOptions wiring in `call.go` / `cmd/bcftools/main.go`. Gated by
  `TestLiveCall/m_G_*`, `/m_V_*`, `/m_keep_*`, `/m_F_prior_freqs`.
- **Index-backed region queries** (`-r` reuses the post-filter path) —
  still uses the linear-scan post-filter; real CSI chunk-seek remains a
  follow-up.

**Validation:** live-oracle byte-parity against the vendored upstream
binary in `TestLiveCall` (built on demand from `reference_code/bcftools`,
no committed goldens).

### `mosdepth`

**Status:** 1 / 1 command, all flags. **`-d/--d4`, `-a/--fragment-mode`,
`-q/--quantize`, `-t/--threads`, `-m/--use-median`, CRAM input, and full
region (`--by`) output are all DONE** (byte-identical to the upstream v0.3.14
binary). No flags remain unimplemented.

Done:

- **Region (`--by`) mode — all output files byte-identical.** In `--by` mode
  we now emit `<prefix>.mosdepth.region.dist.txt` (cumulative region
  distribution; per-base for a BED `--by`, one rounded-mean entry per window
  for an integer `--by`), the `<chrom>_region` + `total_region`
  `summary.txt` rows (per-base aggregate over region-covered bases, in
  upstream's row order), and `per-base.bed.gz` (upstream keeps per-base in
  `--by` mode; suppressed only by `-n`). Zero-coverage references are omitted
  from `summary.txt` / `*.dist.txt` (upstream's BAM-index gate) but retained
  in `per-base.bed.gz` / `regions.bed.gz`; the thresholds region name for
  unnamed/window regions is the literal `unknown`. Validated byte-for-byte
  (BED + fixed window, every produced file) by `TestUpstream_By_AllFiles_Parity`,
  with binary-free unit coverage of the aggregation/cumulation. See
  `UPSTREAM_BUGS.md#mosdepth-region-mode`.

- **CRAM input** — `.cram` inputs are auto-detected (by the `CRAM` file
  magic) and decoded through `pkg/htsgo/alnio.NewReaderWithReference`, which
  dispatches to the in-tree `pkg/htsgo/cram` reader. `-f/--fasta` (and the
  samtools-style `--reference` alias) supplies the decode reference and
  `REF_CACHE` is honoured; an embedded-reference CRAM decodes with no `-f`.
  Depth depends only on alignment coordinates, so per-base / regions /
  summary / distribution / thresholds / quantized outputs are
  **byte-identical** to the equivalent BAM. Validated by
  `TestCRAM_InputMatchesBAM` (BAM-vs-CRAM byte parity across five option
  modes, with and without `-f`, on a samtools-transcoded `ovl.bam`) and
  `TestCRAM_InputMatchesUpstream` (byte-identical to the upstream `mosdepth`
  binary run on the same CRAM).
- **`.csi`** output — now emits a `.csi` (min_shift=14, depth=5, matching
  htslib's `tbx_index_build`) alongside each bgzipped BED output, replacing
  the earlier `.tbi`. Built via `pkg/htsgo/tabix.BuildCSIFromDataFile`.
- **`--mapq 0` fast-path** — when no MAPQ filter is in effect the record
  filter binds a MAPQ-free keep-predicate once, dropping the per-read MAPQ
  comparison from the hot loop. Verified byte-identical to the general path.

Missing:

- **Default-mode overlap-pair correction** — our default (non-fast) mode
  does not subtract double-counted depth where mate pairs overlap; output
  matches upstream's `--fast-mode` (see `UPSTREAM_BUGS.md`). Unchanged by
  this wave.

Implemented:

- **`-a/--fragment-mode`** — counts coverage across the whole template
  (fragment) between properly-paired mates rather than the aligned reads
  only. Only read1 of a proper, non-supplementary pair contributes; it
  covers `[min(read,mate) start, +|TLEN|)`. **Byte-identical** to the
  upstream `mosdepth` v0.3.14 binary on the `full-fragment-pairs` fixture
  (`TestUpstream_FragmentMode_Parity`).
- **`-q/--quantize SEGS`** — bins per-base depth into the user's
  `:`-separated segments and writes `<prefix>.quantized.bed.gz` (plus a
  `.csi`). Mirrors upstream's segment parsing (leading/trailing `:`
  prepend `0` / append the open-ended top bin), the `MOSDEPTH_Q*` label
  overrides, and the "skip depths outside the range" gap behaviour.
  **Byte-identical** to upstream across basic / leading-colon /
  trailing-colon / env-label specs (`TestUpstream_Quantize_Parity`).
- **`-t/--threads N`** — real parallelism for BGZF/BAM decode. A new
  in-tree parallel BGZF reader (`pkg/htsgo/bgzf.MultiReader`,
  symmetric to the existing `MultiWriter`) inflates blocks across N
  worker goroutines and reassembles the decoded byte stream in order, so
  the decoded bytes — and every output file — are **byte-identical** for
  any thread count. Threads < 2 falls back to the sequential reader.
  Verified by `TestThreads_OutputIdentical` (threads {1,2,4,8} identical;
  multi-block fixture spans 66 BGZF blocks) and
  `bgzf.TestMultiReader_MatchesSequential`.
- **`-m/--use-median`** — reports the per-region **median** depth instead
  of the mean in the `--by` regions output, mirroring upstream's
  `imean()` routing through `depthstat.CountStat`: a depth histogram
  (size 65536, top bucket folds in depths ≥ 65535) whose median is the
  first depth where the cumulative count reaches
  `stop_n = int(0.5 + n*0.5)` (round-half-up of n/2; even counts take the
  upper-middle value, not an average). Changes **only** the
  `regions.bed.gz` depth column — the summary, distribution, thresholds,
  quantized, and per-base outputs are untouched, exactly as upstream does.
  The histogram is built from the same `regionStats` sweep that computes
  the mean/threshold columns (single pass, identical depth profile).
  **Byte-identical** to upstream v0.3.14 on the `ovl.bam` MT region
  (`TestUpstream_UseMedian_Parity`), with a mean-vs-median divergence
  cross-check (`TestUpstream_UseMedian_DiffersFromMean`) and direct unit
  coverage of odd/even/empty/cap-fold cases (`TestRegionMedian_Unit`,
  `TestRegionMedian_CapFold`).
- **D4 output** (`-d/--d4`) — writes `<prefix>.per-base.d4` as a real D4
  framefile that is **byte-identical** to the upstream `mosdepth_d4`
  binary's output for the same BAM (same on-disk size). The track uses
  the upstream encoding: a 7-bit-packed primary table with the
  `SimpleRange{0,128}` dictionary, depths ≥ 128 clamped to the all-ones
  code 127 (matching upstream's per-base d4 writer, which leaves the
  secondary table empty), plus the `.metadata`, `.stab` and `.index`
  framefile members. No documented byte exceptions — the files match
  exactly, including the embedded JSON metadata.

**Validation:** `TestD4_UpstreamBinaryParity` downloads the official
`mosdepth_d4` release binary, runs it and our implementation on a fixture
BAM, and asserts the two `.per-base.d4` files are byte-for-byte equal
(verified: 3863 / 3863 bytes, identical). `.csi` validated structurally
and via in-tree round-trip query (`TestRunCsiReadable`,
`TestParity_IndexFiles_Csi`), plus an optional real-`tabix` read when the
binary is on `PATH` (`TestRunCsiReadableByRealTabix`). Fast-path
byte-identity proven by `TestMapqFastPathByteIdentical`.
`--fragment-mode`, `--quantize`, `-t/--threads`, and `--use-median` are
validated byte-for-byte against the upstream `mosdepth` v0.3.14 release
binary (`TestUpstream_FragmentMode_Parity`, `TestUpstream_Quantize_Parity`,
`TestThreads_OutputIdentical`, `TestUpstream_UseMedian_Parity`); the binary is fetched from the GitHub
release with retry/backoff and cached, overridable via `MOSDEPTH_BIN`.
Offline these fall back to internal-consistency assertions and log the
reduced tier (never a silent skip). The broader upstream functional-test
suite is still pending.

---

## Phase plan

### Phase 1 (this PR)

- Truth pass on `PORTING_STATUS.md`, `tools/README.md`, and `analysis/tool_ranking_2026.md`.
- Initial scaffold for `UPSTREAM_BUGS.md` and this file.

### Phase 2 — validated parity audits

For each tool that has an upstream test suite (samtools, bcftools, mosdepth,
vcftools, plus sickle and skewer which have small test corpora), set up
the same kind of validated-parity rig that bedtools got in PR #55:

- Submodule-init the upstream repo.
- Pull representative test cases (5-20 per subcommand) plus their
  expected outputs.
- Diff our output against expected; pass or `t.Skip()` with a documented
  reason.
- Fix small bugs we find inline; record upstream bugs in `UPSTREAM_BUGS.md`.

### Phase 3 — systematic gap closure

Tool-by-tool, dispatched as parallel agent waves of 3-4 PRs each:

1. **bedtools long tail** (~30 small subcommands, parallelisable).
2. **vcftools long tail** (~87 options, parallelisable).
3. **seqtk + fastp + sickle + skewer + prinseq closure** (small individual
   tools, all parallelisable).
4. **samtools mpileup + bcftools call + bcftools annotate** (the big
   variant-calling loop).
5. **`samtools markdup`/`idxstats`/`stats`/`merge`** (workhorse
   utilities).
6. **bcftools merge / isec / sort / annotate** (everyday set-ops).
7. **CRAM** (last; multi-week effort on its own).

Phase 1 lands first. Phases 2 and 3 are several sessions of parallel-agent
waves each.
