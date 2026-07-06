# `samtools view` region / BED → SAM speed & RSS (C3 outlier close-out)

This note records the profiling and the fixes applied to the `samtools view`
**region-query / `-L <bed>` → SAM** path — the last C3 performance outlier
historically reported at ~12× slower than upstream with ~11× peak RSS.

Data: real GIAB HG002/GRCh38 **chr20** BAM
(`pipeline/.fixtures/realchr20/chr20.bam`, 2,054,615 reads, 240 MB) with its
`.bai`, and the co-located 10,192-interval `chr20.bed`. Oracle:
`bin/upstream/samtools` (htslib 1.23.1, `libdeflate=no`). Ours:
`bin/ours/samtools`. All timings warm-cache, best of 4 reps, `/usr/bin/time -l`.

## Headline: the 12× is already gone; the residual is inflate-bound

The BAM→SAM **direct-serialise fast path** (raw BGZF bytes →
`sam.BAMReader.WriteSAMBody`, no intermediate `sam.Record`) already covers the
region and BED-filtered cases via `viewIndexedChunksFast` /
`scanChunksFast` (indexed seek) and `viewStreamFast` / `fastSAMScan`
(linear scan). Measured today, **before** any change in this pass:

| query | ours (orig) | upstream | ratio | ours RSS | upstream RSS |
|---|---|---|---|---|---|
| `view chr20.bam chr20` (whole-chr, indexed) | 3.83 s | 1.55 s | 0.40× | 20 MB | 7 MB |
| `view -L chr20.bed chr20.bam` (linear + BED) | 3.86 s | 1.45 s | 0.38× | 21 MB | 7 MB |

So the outlier is **~2.5×**, not 12×, and RSS is **~3×**, not 11×.

## Profile (CPU + allocation)

`pipeline/bench/viewprof_main.go` (a `//go:build ignore` CPU-profile harness)
plus package benchmarks `BenchmarkViewRegionSAM` / `BenchmarkViewBEDSAM`
(`tools/samtools/pkg/samtools/view_region_bench_test.go`, `-benchmem`) located
the costs:

* **Wall / user CPU is dominated by DEFLATE inflate.** `/usr/bin/time` shows
  `sys` ≈ 0.17 s but `user` ≈ 3.7 s. The klauspost pure-Go inflater carries the
  bulk of that user time; upstream's zlib inflate (this htslib build has
  `libdeflate=no`, so it is *not* using the fast libdeflate path) is simply
  faster per byte. This cost is **irreducible without a faster in-tree
  inflater** — it is the documented floor for every BGZF-consuming Go tool here.
* **Allocation churn on the hot per-record loop**, which drives GC and RSS:
  - `sam.BAMReader.ReadSAMInto` allocated **~1 object per record** — a
    stack-local `[4]byte` block-size prefix that escaped to the heap when handed
    to `io.ReadFull` through the `io.Reader` interface (~2.9 M objects).
  - the BED fast filter called `IntervalTree.Query`, which allocates a
    `[]*Record` results slice *and* a query `*Record` **per alignment**
    (~8.2 M objects on the BED path) when all it needs is a yes/no answer.
  - the BGZF reader issued several tiny `read()` syscalls per 64 KB block
    (12-byte header, XLEN extra, payload, 8-byte footer) straight to the
    unbuffered `*os.File`.

## Fixes applied (all byte-exact, logic preserved)

1. **`pkg/htsgo/sam/bam_reader.go` + `bam_sam_fastpath.go`** — hoisted the
   4-byte block-size scratch to a reusable `BAMReader.sizeBuf` field so
   `ReadSAMInto` no longer heap-allocates once per record.
2. **`pkg/htsgo/bed/intervaltree.go`** — added an allocation-free,
   short-circuiting `IntervalTree.Overlaps(start, end) bool`; the view `-L`
   fast filter (`view_fastpath.go`) now calls it instead of
   `len(Query(...)) > 0`, dropping both the per-record query `*Record` and the
   results slice. Cross-checked against `Query` in a new unit test.
3. **`pkg/htsgo/bgzf/bgzf.go`** — interposed a 256 KiB `bufio.Reader` between
   the caller's `*os.File` and the BGZF block parser (skipped when the source
   already buffers, e.g. `*bytes.Reader` in tests), collapsing the per-block
   tiny reads into a handful of large sequential reads. Virtual-offset counting
   is unchanged (the `countingReader` still counts exactly the bytes consumed),
   so seek/chunk semantics and `RawRemaining` (samtools cat) are unaffected.

## After

Package benchmark (`-benchmem`, 5×), allocation & bytes per full-chr / BED pass:

| bench | before allocs | after allocs | before B/op | after B/op |
|---|---|---|---|---|
| `ViewRegionSAM` (whole-chr) | 2,178,607 | **92,898** (23× fewer) | 18.4 MB | 10.1 MB |
| `ViewBEDSAM` (`-L`) | 4,096,686 | **111,330** (37× fewer) | 31.3 MB | 7.7 MB |

CLI wall / RSS (warm, best of 4), before → after this pass:

| query | wall before → after | upstream | RSS before → after | upstream RSS |
|---|---|---|---|---|
| `view chr20.bam chr20` | 3.83 s → **3.78 s** | 1.55 s | 20 MB → **17 MB** | 7 MB |
| `view -L chr20.bed chr20.bam` | 3.86 s → **3.82 s** | 1.45 s | 21 MB → **18 MB** | 7 MB |

**Byte-exact vs upstream: PASS** for all three probes (whole-chr region, 1 Mb
indexed slice, `-L` BED) — verified with `cmp` on the full SAM output.

## Honest floor

Wall time moved only marginally (~1–2 %) because the path is **inflate-bound**,
not allocation- or serialise-bound: the SAM serialiser (`WriteSAMBody`), the BED
interval query and the region overlap test together are <10 % of samples; the
remaining ~2.5× wall gap is the klauspost-vs-zlib inflate throughput
difference. That is the **documented irreducible floor** for this path short of
adopting a faster inflater (a separate, larger decision — see the CRAM/BGZF
codec discussion in `CLAUDE.md`). The safe wins realised here are on the
**allocation / GC / RSS** axis: 23–37× fewer heap objects on the hot loop and a
~15 % peak-RSS reduction, with byte-identical output.
