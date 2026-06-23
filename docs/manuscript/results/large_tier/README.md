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

## Status

- Large fixtures generate fine (~22 GB, on the host-mounted disk).
- Multi-contig bug found and confirmed (above) — actionable.
- Full large matrix not completed on this box: the few >12 GB-output cells need
  the fat node or a stream-comparing harness.
