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

## The per-cell memory ceiling (why the full large matrix doesn't finish here)

`GOMEMLIMIT` let the whole **medium** tier pass (it bounded cross-cell RSS
growth). At **large**, a few individual cells produce a single output bigger than
the 12 GB VM — and the harness buffers each cell's entire ours+upstream stdout in
RAM (`pipeline/runner.RunEntry`), so those cells OOM regardless of `GOMEMLIMIT`
(it cannot shrink a *live* buffer):

- `bcftools mpileup_heavy` / `call` over 192 Mbp,
- (likely) `bedtools genomecov -d` per-base over 192 Mbp,
- `samtools view` of the 2.5 M-read BAM to SAM.

These need either a bigger box (the runbook's 32–64 GB fat node) or a
**stream-comparing** harness (write each side to a temp file and diff files, the
way `pipeline/cmd/realparity` already does for its file-producing cells) instead
of buffering whole outputs in RAM. The lighter large cells (and the `norm`
DIVERGEs above) run fine here.

## Large-tier performance (bench, reps=10, robust stats)

Run **serially per format-group** (uncontended) with `GOMEMLIMIT=8GiB`. Each cell
times our binary against upstream over 10 reps; `wall×` is the **median** ratio
(`ours/upstream`, `< 1.0` = we are faster) with its **95 % bootstrap CI** (H1a).
Per-group reports: `bench/{FASTQ,CRAM,BED}/bench.md`. The FASTQ/CRAM/BED groups
completed; the BAM and VCF groups' light cells completed (numbers below) but
their **heavy compute cells OOM at large** for the same reason as the parity
matrix above — `bcf_call` runs over a ~20 GB intermediate mpileup VCF, etc. (see
the memory ceiling).

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
| `sam_mpileup`, `bcf_call`, `bcf_isec` | — | — | **OOM at large** (>12 GB output on this box) |

The large-tier picture matches medium: I/O-bound conversions, `bedtools`
intersect/coverage/genomecov, and `sickle` are **faster** than upstream; the
compute-heavy cells are slower; and the same handful of cells that exceed the
12 GB box's memory are reported as OOM rather than hidden. `bcf_norm`'s 48×
RSS and `sam_depth`/`sam_view_bam2cram`'s ~11× RSS are the memory-side
optimisation targets (tracked in the real-data perf follow-ups).

## Status

- Large fixtures generate fine (~22 GB, on the host-mounted disk).
- Multi-contig bug found and confirmed (above) — actionable.
- Full large matrix not completed on this box: the few >12 GB-output cells need
  the fat node or a stream-comparing harness.
