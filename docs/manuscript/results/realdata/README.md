# Real-data tests on GIAB NA12878 exome (merge + CRAM)

Inputs (real, multi-contig): two GRCh37 exome BAM lanes
(`project.NIST_…_1/2_NA12878.bwa.markDuplicates.bam`, 39.7 M + 40.5 M reads) and
the `hs37d5.fa.gz` reference, pulled from the **AWS S3 GIAB mirror**
(`s3.amazonaws.com/giab`, ~100 MB/s; the NCBI FTP/HTTPS mirrors were slow and
aria2c multi-segment over FTP corrupted the files). Both BAMs verified valid
(`quickcheck` + full decode). Run in the `linux/arm64` container with
upstream samtools 1.22.1 as the oracle. Log: `merge_cram_giab.log`.

These are **parity/functionality** findings (md5 of decoded records, provenance
excluded), not timings.

## Finding A — `samtools merge` records differ from upstream

`samtools merge` of the two lanes (80,232,730 reads):

| side | md5 of merged records |
|---|---|
| upstream | `be940de6f4cd57179277315fae8c2907` |
| ours | `8323aa2961ad5c482a7cce9e82906418` |

**MERGE PARITY: DIFF.** Root-caused on a small region (chr20:1–2 Mb, 71 797
merged reads). Two candidate causes — only the first is a real bug:

1. **RG-collision suffix is NOT a bug.** Both lanes carry `@RG ID:1`, so merge
   disambiguates the second. Upstream's suffix is `%08lX` of `lrand48()`
   (`bam_sort.c:408`) seeded by `(long)time(NULL)` (`bam_sort.c:1622`) — i.e.
   **time-random**; ours is too. Re-running both with `-s 1` makes the `@RG` IDs
   **identical** (`1-055424A4`), which also **proves our `lrand48` is byte-exact
   to htslib's**. So this difference is inherent to upstream (no `-s`), not a port
   defect.
2. **Aux-tag repositioning — the real bug.** With `-s 1` (suffix matched), every
   read still differs: upstream's `bam_translate` (`bam_sort.c:947`) does
   `bam_aux_del` + `bam_aux_append` for `RG` then `PG` on **every** read (even
   identity translations), moving both tags to the **end** of the aux list; ours
   updates `RG` in place and never repositions `PG`. So ours emits `…MD PG RG XG…`
   while upstream emits `…MD XG…XT RG PG`. Our plain `view -b` round-trip
   preserves tag order byte-for-byte (verified), so this is specific to merge.
   **FIXED.** `merge.go` now mirrors `bam_translate` — `bamTranslateAux`
   del+re-appends `RG` (via `rgTrans`) then `PG` (via a new `pgTrans`, with
   `PP`-chain translation) at each record's aux-list end. On the chr20 region
   (same seed) the merged **records are now byte-identical to upstream**
   (143 594 → 0 differing record-lines); regression tests added, merge + broader
   samtools tests green. A separate minor gap remains: the merged header groups
   `@RG`/`@PG` lines differently (content identical, line order differs) — does
   not affect record bytes.

Also observed (separate gaps): our `merge -R <region>` ignores the region (merges
the whole file), and our merge is markedly slower than upstream on large inputs.

## Finding B — our CRAM **decoder** fails on real data (rANS4x16 RLE)

Reference-based round-trip (`view -C -T hs37d5`), all writer/reader combos vs the
upstream/upstream baseline:

| write | read | cram MB | result |
|---|---|---|---|
| up | up | 2963 | PASS (baseline) |
| ours | up | 4276 | DIFF (non-empty, differs from baseline) |
| ours | ours | 4276 | DIFF (empty — decode failed) |
| up | ours | 2963 | DIFF (empty — decode failed) |

Every `read=ours` is **empty** because our CRAM decoder aborts with:

```text
samtools view: cram: container 1 slice 0: cram: slice external block
(content id 12): cram: decompressing rans4x16 block (content id 12):
rans4x16: RLE run expands past the declared output size 942308
```

- **Location:** `pkg/htsgo/cram/codec/rans4x16_transform.go:779` `htsRLEDecode`
  (the bounds check at line 805).
- **Significance:** this rANS4x16 RLE path passes the synthetic htscodecs
  conformance corpus but fails on real GIAB data — a real test-coverage gap.
  The header decodes fine; a data block (content id 12) overruns its declared
  output during RLE expansion.
- **Root-cause (diagnostic).** The RLE loop itself matches htscodecs exactly
  (both loop over every literal byte; bounds checks are arithmetically
  equivalent), so the bug is **upstream of the RLE step**. Instrumenting the
  failure: at `litIdx=626098/638416`, output was `942288/942308` — only 20 bytes
  of room left, yet **12 318 literal tokens still remained**, which cannot fit.
  So our rANS4x16 **order-1 literal decode produces a wrong stream** (wrong bytes
  and/or wrong length) for this block; the RLE then necessarily overruns.
  htscodecs requires `rle_len ≤ osz` and uses the same lengths, yet decodes the
  same CRAM fine — so the divergence is in our O1 rANS decode of this block.
  Fixing needs a byte-level comparison of our O1 rANS output vs htscodecs for the
  failing block (a focused codec session); a blind change here risks silent data
  corruption.
- **Encoder side:** `write=ours read=up` is non-empty but differs from the
  upstream baseline, and our CRAM is larger (4276 vs 2963 MB) — so our CRAM
  *encoder* also diverges from upstream (worse ratio + different records),
  though upstream can read it.

## Finding C — fastp adaptor trimming: works, comparable, not byte-identical

Our ported `fastp` vs upstream `fastp` on 1 M real exome read-pairs
(NIST7035 L001), default + `--detect_adapter_for_pe`:

- **Detection works, but PE trimming is broken on read 2.** Ours auto-detects
  the Nextera adapter `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`, but the
  **residual-adapter check** (reads still containing the adapter core after
  trimming) exposes a real bug:

  | read | raw | ours | upstream |
  |---|---|---|---|
  | R1 | 4.98 % (49 843) | **0.023 %** (224) | 0.003 % (25) |
  | R2 | 4.91 % (49 067) | **4.29 % (41 684)** | 0.002 % (22) |

  Our fastp clears R1 read-through adapter (4.98 → 0.02 %) but **barely trims R2**
  (4.91 → 4.29 %, only ~15 % removed); upstream clears both to ~0 %. So in PE
  overlap mode our adapter trimming is **one-sided** — it leaves adapters in
  read 2.

  **Fuzzy-detection cross-check (fastp 2×2).** Re-running fastp (which uses
  fuzzy k-mer overlap detection, not exact match) as the QC oracle over each
  trimmed set — residual reads it still detects+trims:

  | QC engine | on upstream-trimmed | on ours-trimmed |
  |---|---|---|
  | upstream fastp | 0 | **120 804** (1.68 M bases) |
  | our fastp | 0 | **0** |

  Two defects: (a) upstream's oracle finds **120 804 reads still carry adapter in
  our output** (~12 % of pairs) — trimming is incomplete; and (b) **our fastp
  reports 0 residual in its own output**, i.e. its adapter *detection* is blind
  to the same R2 read-through it fails to trim.

  **FIXED.** Our PE pipeline was missing upstream's overlap-based adapter trim
  (`AdapterTrimmer::trimByOverlapAnalysis`): it trimmed both mates using only the
  read-1 adapter sequence. Ported it (`overlap.go: trimByOverlapAnalysis`,
  wired into `processPairOnce`) so a read-through pair is trimmed to the insert
  length on **both** reads, removing the read-2 adapter regardless of sequence.
  After the fix, upstream's QC oracle finds **2 residual reads in our output**
  (== upstream's own 2); regression tests added; fastp tests green. A smaller
  read-count gap remains (ours retains ~2 % more pairs — a min-length /
  read-2-sequence-trim nuance), tracked separately.
- **CLI incompatibility (drop-in gap).** Our fastp's short flags differ from
  upstream's: ours uses `-I`=in1, `-O`=out1, `--in2`/`--out2` for read 2 and
  `--json`/`--html` (no `-i`/`-o`=out2 / `-j`/`-h`), whereas upstream uses
  `-i`/`-I`/`-o`/`-O` and `-j`/`-h`. A standard upstream fastp command scrambles
  I/O on our port. This breaks drop-in CLI compatibility and should be aligned.

## Finding D — QC/format tool sweep on real FASTQ

Ran each ported QC tool vs its upstream on the real NIST7035 exome reads:

| tool | result |
|---|---|
| **seqtk** | **7/7 PASS** byte-exact (`seq -A`, `seq`, `comp`, `seq -r`, `trimfq`, `fqchk`, `seq -q20 -n N`); CLI drop-in compatible. |
| **sickle** | **3/3 PASS** byte-exact (`pe` sanger: out1/out2/singles); CLI drop-in compatible. |
| **skewer** | adapter trimming **correct** (residual matches upstream exactly, R1+R2), but **CLI incompatible** (ours `skewer pe -i/-j/-o/-p`; upstream `skewer -m pe … -o prefix`), and ours keeps **2 pairs** upstream drops (a min-length-after-trim nuance). |
| **prinseq** | **CLI incompatible** — ours is subcommand-based (`prinseq filter -i/-o --fastq`); upstream is flat (`prinseq-lite.pl -fastq -out_good -min_len …`). |
| **fastp** | Finding C: PE R2 adapter trimming/detection broken; CLI incompatible. |

**Systematic gap:** `fastp`, `skewer`, `prinseq` are **not drop-in CLI replacements**
— they were redesigned with subcommands / different flag names, contrary to the
project's "drop-in POSIX CLI" goal. `seqtk` and `sickle` are compatible.

## Caveat — embed_ref / no_ref modes not validly tested

The `embed_ref` and `no_ref` round-trips produced **0-byte CRAMs for both
engines**, so their apparent "PASS" (empty == empty) is spurious — a bug in the
test harness's option invocation for those modes, **not** a real result. They
need the command corrected (and the decoder bug fixed) before they mean anything.
Only the reference-based mode above is valid.

## Finding E — `realparity` differential battery (H2a, whole-file)

A systematic differential-parity battery (`pipeline/cmd/realparity`) run on the
real GIAB inputs: `hs37d5.fa` reference + the NA12878 exome BAM (39.7 M reads) +
the HG002 GRCh37 benchmark VCF, each samtools/bcftools cell run on **our** port
and the **upstream** binary with provenance-stripped byte-exact comparison.

**Data-naming discovery (a real-world reference mismatch).** The BAM is
**hg19/UCSC-named** (`chrM`, `chr1`, …, `chr17_gl000203_random`) while `hs37d5`
and the VCF are **GRCh37-named** (`MT`, `1`, …). None of the BAM contigs match
the reference by name — which exercised exactly the edge cases a clean dataset
would not.

**Harness memory-safety (enabling the run).** The battery originally buffered
*both* sides' fully-decoded output (`CompareByteExact` on two `[]byte`), so a
multi-GB BAM decoded to ~15–20 GB of SAM ×2 and OOM-killed the 12 GB VM before
comparing. Refactored to a **streaming provenance-stripping md5 digest** (no
full-output buffers, no temp files; harness peak RSS ~12 MB), proven byte-for-
byte equivalent to `CompareByteExact` by a 200 k-input fuzz (commit `fd6e308`).
This first complete pass surfaced the bugs below.

**Bugs found and fixed** (each byte-exact vs the upstream oracle on real data,
independently reviewed, committed):

| # | bug | symptom on real data | commit |
|---|---|---|---|
| 1 | CRAM **reference-free decode** | `view_cram` decode aborted `fasta: contig "chrM" not in index` — our decoder ignored the slice `RR=0` (reference-not-required) flag and tried to fetch the absent contig; upstream encodes/decodes such slices reference-free | `3769fda` |
| 2 | CRAM **@SQ M5/UR** | encoding `-T ref` omitted the `M5`/`UR` tags upstream injects; UR must be added iff the @SQ has an M5 (computed or pre-existing), bare otherwise — verified against the upstream binary | `aa5c784`, `3769fda` |
| 3 | CRAM **inline MD/NM aux order** | decoded records differed from upstream purely in aux order — upstream pulls inline `MD`/`NM`/`RG` to the tail (`MD,NM,RG`), ours kept `MD`/`NM` inline | `d74a9b4` |
| 4 | `samtools stats` | GCD GC% used `float64` not upstream's `float` (`e3e98ea`); and the CIGAR/indel/NM-mismatch walk wasn't gated on `IS_UNMAPPED`, folding in BWA unmapped-with-CIGAR primary reads (`bases mapped (cigar)` +505, `mismatches` +20), plus integer-truncated average length | `e3e98ea`, `986d520` |
| 5 | `samtools depth` **OOM** | `depth -a -r chr20` (63 M positions) buffered its whole output → **11.4 GB, OOM-killed**; rewritten to stream like `bam2depth.c` → **72 MB**, byte-exact | `c0c573e` |
| 6 | `samtools faidx` **missing** | the subcommand did not exist (only `dict`); implemented `faidx`/`fqidx` byte-exact (index build + extract, plain + bgzipped), with a streaming bgzipped index build (hs37d5 `.fai`+`.gzi` at ~16 MB, not the ~3 GB genome) | `deb6e01` |
| 7 | `realparity` harness | `depth`/`bcftools query`/`bcftools stats` take `-r REGION`, not a positional region (which upstream reads as a filename) | `a7e7177` |

**Final verdict — `PASS=15 / DIVERGE=0 / ERROR=0`** (uncontended, reps=3,
`giabfinal/`). Every cell is byte-exact after provenance stripping on full-scale
real data: `view_sam`/`view_sam_header` (chr20), `flagstat` (39.7 M reads),
`idxstats`, `stats` (whole BAM), `depth -a` (chr20, 63 M lines), `quickcheck`,
`view_bam`, **`sort`** (whole 2.8 GB BAM, no OOM), `view_cram` (reference-free),
and `bcftools view`/`view_body`/`norm`/`stats`/`query`. (The first complete
re-run scored 14 PASS / 0 DIVERGE; its lone ERROR — the `depth` OOM — was then
fixed and this run confirms it.)

### Real-data performance (min over reps; `ratio = ours/upstream`)

Correctness is byte-exact everywhere; the timing/memory record surfaces
optimisation targets (not defects):

| cell | wall× | note |
|---|---|---|
| `samtools_quickcheck` | 0.44 | faster |
| `samtools_view_bam` | 0.76 | faster (whole 2.8 GB re-encode) |
| `bcftools_view` / `view_body` / `query` | 0.78–0.87 | faster (whole VCF) |
| `samtools_view_cram` | 1.04 | par wall (CPU 4.8× — our reference-free encode does more work) |
| `samtools_idxstats` | 1.29 | par |
| `samtools_flagstat` | 1.38 | slower |
| `samtools_sort` | 1.50 | slower (wall; CPU 2.9×, RSS 6.0 GiB vs 0.9 GiB) |
| `samtools_stats` | 1.56 | slower |
| `bcftools_norm` | 5.00 | **RSS 9.1 GiB vs 21 MiB (436×)** — buffers the whole VCF |
| `samtools_view` (region→SAM) | 12.7 | **slow** — region-query + SAM serialisation |
| `samtools_depth_a` | 27.7 | **slow** — correct + bounded (73 MB) but slow vs upstream pileup |

`bcftools_stats` (123×) is the empty `-r chr20` cell (the VCF is GRCh37-named
`20`, not `chr20`), so its ratio is noise.

> **Update (subsequently fixed, still byte-exact):** `bcftools norm`'s 9.1 GiB
> RSS is gone — it now streams with a bounded reorder window at **18 MB**
> (commit `4e53d2f`); and `samtools depth -a -r chr20` dropped from **~72 s to
> ~10 s** (≈ upstream) via an indexed BGZF seek + lean decode (`9766251`). The
> remaining `view` region→SAM speed is a lower-priority follow-up
> (`docs/PARITY_ROADMAP.md`).

## Status — real-data bugs found

The GIAB exome surfaced several genuine parity gaps. The bulk were then fixed
(each validated only on real data vs the upstream binary oracle, with an
independent review agent re-running the comparison).

**Fixed (validated vs upstream oracle on real data, committed):**

- **`bcftools norm`** multi-contig output ordering (see `../large_tier/`).
- **`samtools merge`** aux-tag reposition (Finding A) — records byte-identical to
  upstream. Minor `@RG`/`@PG` header line-grouping gap remains.
- **fastp PE read-2 adapter** (Finding C) — overlap-based PE adapter trim ported;
  residual 120 804 → 2 reads (== upstream).
- **CRAM rANS4x16 decoder** (Finding B) — the failure was the X_32 (32-way) RLE
  meta-coder selection, not the literal decode; our samtools now decodes the real
  GIAB CRAM byte-identical to upstream across all 39.7 M records.
- **fastp / prinseq / skewer CLI** (Finding D) — now accept upstream's drop-in CLI
  (`fastp -i/-I/-o/-O/-j/-h`; `prinseq-lite.pl` flat flags; `skewer -m pe -o prefix`).
  fastp also closed the ~2 % read-retention gap; prinseq + skewer-SE byte-exact.
- **`samtools merge -R`** region — implemented (indexed region query); a region
  merge is byte-identical to upstream.

**Also fixed (second multi-agent round, validated vs upstream oracle):**

- **skewer PE trimming** — ported skewer's PE overlap analysis + base
  error-correction; `-m pe` trimmed pairs byte-identical to upstream (the ~2
  dropped pairs are a faithful reproduction of upstream's own final-partial-block
  data loss).
- **`samtools merge`** `@RG`/`@PG` header line-grouping now matches upstream
  (0-diff header), records still byte-identical.
- **CRAM MD/NM regeneration** on decode was already implemented and correctly
  gated to external-reference CRAM; confirmed byte-exact vs upstream and pinned
  with a regression test (the earlier delta was an embed_ref test artifact).

**Remaining:** `samtools merge` is slower than upstream on large inputs
(performance, not correctness — deferred to the performance pass).

Downloads used the **AWS S3 GIAB mirror** throughout (fast + reliable byte-range).
Fixes were produced by a bounded multi-agent run (fix + independent review agent
per bug, ≤3 concurrent), validating only against the upstream binary on real data.
