# Large tier (16-contig) — a real bug found, plus the per-cell memory ceiling

Run on the local laptop container ([`../hardware.md`](../hardware.md)) with the
VM at 12 GB and `GOMEMLIMIT=8–9GiB GODEBUG=madvdontneed=1`. The large tier is the
first to use **16 contigs** (chr1…chr16), and that multi-contig coverage
immediately paid off.

## Finding: `bcftools norm` sorts output by contig *name*, not header order

Three large-tier cells DIVERGE — `bcftools_norm_check_ref`, `bcftools_norm_dedup`,
`bcftools_norm_split` — all with the same signature. Confirmed directly outside
the harness (`norm_contig_order_bug.txt`):

```text
ours:     chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr2 chr3 … chr9
upstream: chr1 chr2 chr3 … chr9 chr10 chr11 … chr16
```

The input VCF has its contigs in numeric order in **both** the `##contig` header
and the records. Upstream `norm` preserves that header (rid) order; **ours
re-sorts by contig name as a string**, so `chr10` sorts before `chr2`.

- **Root cause:** `tools/bcftools/pkg/bcftools/norm.go:1843` `sortVariants` does
  `variants[i].Chrom < variants[j].Chrom` (string compare) then `Pos`. It should
  order by the VCF header's contig index (rid), like upstream.
- **Why earlier tiers missed it:** smoke/small/medium use ≤8 contigs (chr1…chr8),
  where lexical and numeric order coincide. The bug only appears at **≥10
  contigs**, where `chr10` < `chr2` lexically. This is precisely the multi-contig
  behaviour the large tier exists to exercise.
- **Fix (applied):** `sortVariants` now takes the header contig order
  (`contigOrder(hdr)`, the existing helper from `concat.go`) and compares by
  contig index, not name — off-header contigs sort after declared ones,
  lexically, matching bcftools. After the fix, `bcftools norm -m -any` over the
  16-contig fixture is **byte-identical to upstream on every data row** (only the
  stripped `##bcftools_*` provenance lines differ). A regression test,
  `TestNormOrdersByHeaderContigNotLexical`, feeds chr10/chr2/chr1 out of order
  and asserts header (chr1,chr2,chr10) output order. All norm unit + parity tests
  pass; `gofmt`/`vet`/`build` clean.

This was a genuine fix-on-port parity defect, **not** an arm64/`libc++`/FP
artifact — and it is now corrected.

## The heavy cells were disk-bound, not memory-bound (resolved)

An earlier run reported `sam_mpileup` / `bcf_call` / `bcf_isec` as OOM at the
large tier. That diagnosis was wrong: the tools themselves are **bounded RAM**.
Direct peak-RSS measurement at large is **106 / 18 / 464 MB** (ours) vs
**42 / 10 / 11 MB** (upstream) — all exit 0. Two distinct disk problems, not a
memory one, were aborting the cells:

- **`mpileup` / `call` huge stdout.** Each writes a ~17–20 GB VCF; the bench was
  spilling it to a temp file on a container overlay with ~5 GB free, so it failed
  on `ENOSPC`. Fixed by streaming both sides' stdout to `/dev/null` (the
  `write()` cost is still counted symmetrically) — committed in
  *"fix bench harness disk-fill on huge stdout"*.
- **`isec` output prefix dir.** `bcftools isec -p DIR` writes real files; with
  `DIR` under a small-overlay `TMPDIR` it filled the disk. Fixed by pointing
  `TMPDIR` at the 205 GB host mount.

With both fixes, all three complete at large; their numbers are folded into the
figures (`figures/bench_oom_large.json`). `bcf_isec`'s memory has since been
fixed (streaming k-way merge + byte-bounded batch: RSS 39× → ~6.5×, multi-GB →
~80 MB on a many-contig human-scale corpus). Its earlier ~15× **wall** turned
out to be a benchmarking artifact: isec writes its `-p` output to disk, the bench
used the slow Docker bind mount, and our VCF / `sites.txt` writers were
under-buffered — fixed (256 KiB buffers), and on a fast disk isec is a steady
**~2.1–2.5×** (gap G6 is now just the minor multi-sample FORMAT over-decode).

> Separately, the **parity matrix** harness (`pipeline/runner.RunEntry`) still
> buffers each cell's entire ours+upstream stdout in RAM to byte-diff them, so it
> *would* OOM on the ~17 GB heavy cells. That is a harness limitation, not a tool
> one; it wants the **stream-comparing** approach `pipeline/cmd/realparity`
> already uses (write each side to a file, diff files). Tracked as gap G7.

## Large-tier performance (bench, reps=10, robust stats)

Run **serially per format-group** (uncontended) with `GOMEMLIMIT=8GiB`. Each cell
times our binary against upstream over 10 reps; `wall×` is the **median** ratio
(`ours/upstream`, `< 1.0` = we are faster) with its **95 % bootstrap CI** (H1a).
Per-group reports: `bench/{FASTQ,CRAM,BED}/bench.md`. All groups completed,
including the heavy `sam_mpileup` / `bcf_call` / `bcf_isec` cells once the bench
stopped spilling their huge output to a small-overlay temp dir (see above).

| cell | wall× (median) | 95% CI | note |
|---|---|---|---|
| `sickle_se` | 0.60 | [0.50, 0.69] | faster |
| `bed_intersect_pair` | 0.64 | [0.63, 0.66] | faster |
| `sam_view_bam2bam` | 0.70 | [0.69, 0.71] | faster |
| `bed_intersect_self` | 0.71 | [0.68, 0.72] | faster |
| `bed_coverage` | 0.75 | [0.72, 0.78] | faster |
| `sam_view_bam2cram` | 0.79 | [0.78, 0.81] | faster (RSS 11×) |
| `bed_genomecov` | 0.84 | [0.83, 0.86] | faster |
| `sam_sort` | 0.88 | [0.87, 0.92] | faster |
| `bed_sort` | 0.99 | [0.95, 1.01] | par |
| `bcf_stats` | 1.05 | [1.04, 1.07] | par |
| `sam_stats` | 1.12 | [1.10, 1.14] | par |
| `bcf_view` | 1.14 | [1.11, 1.15] | par |
| `sam_view_cram2bam` | 1.26 | [1.25, 1.27] | slower |
| `sam_flagstat` | 1.31 | [1.30, 1.33] | slower |
| `seqtk_seq` | 1.31 | [1.28, 1.53] | slower |
| `sam_depth` | 1.36 | [1.34, 1.38] | slower (RSS 11×) |
| `bcf_norm` | 1.42 | [1.37, 1.61] | slower (**RSS 48×** — buffers the VCF; see real-data #20) |
| `bcf_query` | 1.58 | [1.54, 1.59] | slower |
| `bed_merge` | 1.82 | [1.73, 1.88] | **slow** |
| `bcf_call` | 1.76 | [1.75, 1.77] | slower (RSS 1.6×) |
| `sam_mpileup` | 2.20 | [2.13, 2.27] | **slow** (RSS 2.6×) |
| `bcf_isec` | 2.48 | [2.44, 2.50] | slower (fast-disk; RSS ~6.5×). The earlier 15× was a slow-bind-mount + under-buffered-write artifact, now fixed. |

The large-tier picture matches medium: I/O-bound conversions, `bedtools`
intersect/coverage/genomecov, and `sickle` are **faster** than upstream; the
compute-heavy cells are slower. `bcf_isec`'s memory is bounded — the k-way
position-window merge with a byte-bounded batch dropped its peak RSS from 39× to
~6.5× upstream (multi-GB → ~80 MB on a many-contig human-scale corpus). Its
**wall** is a steady **~2.1–2.5×** across tiers (it had looked 15× at large only
because the bench wrote isec's `-p` output to the slow Docker bind mount with
under-buffered writers; both are now buffered and the figure uses fast-disk
numbers). `bcf_norm`'s 48× RSS and `sam_depth`/`sam_view_bam2cram`'s ~11× RSS are
the other memory-side optimisation targets (tracked in the real-data perf
follow-ups).

## Status

- Large fixtures generate fine (~22 GB, on the host-mounted disk).
- Multi-contig bug found and confirmed (above) — actionable.
- **Full large matrix completed**: the heavy `mpileup`/`call`/`isec` cells run
  once their huge output is streamed to `/dev/null` (bench) and `TMPDIR` points
  at the big host disk; the tools are bounded RAM, not OOM. The parity *matrix*
  harness still buffers outputs in RAM (gap G7) and wants the stream-comparing
  diff before it can byte-check those cells at large.
