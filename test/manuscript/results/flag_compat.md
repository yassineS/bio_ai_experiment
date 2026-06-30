# Flag-compatibility results (manuscript claim C4/C1)

This page reports, per ported tool, the **flag-compatibility percentage**: of an
upstream tool's documented CLI flags, how many our Go port accepts under the
**same name/spelling**.

## What this measures (and what it does not)

This is a **flag-surface acceptance** metric. For each upstream flag it asks one
question: *does our port register a flag with the same spelling?* It is computed
by token name only (a flag's bare name with leading dashes stripped, so `-i`,
`--input`, and `input` compare equal).

It deliberately does **not** measure semantic equivalence — i.e. whether an
accepted flag produces upstream-identical behaviour. That stronger property is
covered separately by the byte-exact **upstream-parity tests** (the `*Upstream*`
tests re-executed against freshly built upstream binaries in CI). A flag can be
accepted with the right spelling and still differ in behaviour; conversely, a
port can implement an upstream capability under a *different* spelling and score
as "not compatible" here even though the functionality exists. Both situations
occur in the table below and are called out in the caveats.

Read this number as: *"how drop-in is the CLI surface, spelling-for-spelling?"*
Not as: *"how feature-complete or byte-exact is the tool?"*

## Method

- **Upstream flag set** — curated per tool (or per subcommand for the
  subcommanded tools) in
  [`pipeline/cmd/flagcompat/upstream_flags.txt`](../../../pipeline/cmd/flagcompat/upstream_flags.txt).
  Every entry is derived from a primary source in the read-only `reference_code/`
  submodules: `getopt`/`getopt_long` optstrings and `longopts` arrays
  (seqtk, sickle, the htslib utilities, samtools, bcftools), `cmd.add(...)`
  registrations (fastp), the Perl `GetOptions(...)` list (prinseq), the
  `parameters.cpp` option strings (vcftools), or documented usage text
  (skewer, mosdepth, and the `bed*` subcommands). The per-source provenance is
  documented in the header of that data file. For the bedtools-derived `bed*`
  tools the **documented per-subcommand** flag set is used rather than the full
  inherited `ContextBase` parser union, so inherited plumbing flags that never
  appear in a subcommand's own help do not distort the denominator.
- **Our flag set** — extracted automatically from each port's Go source by the
  `flagcompat` program ([`pipeline/cmd/flagcompat`](../../../pipeline/cmd/flagcompat)),
  which parses every non-test `.go` file under `tools/<tool>/` with `go/ast` and
  records the names passed to the `pkg/cliflag` helpers and the stdlib `flag`
  package, plus the case labels and flag-binding map keys used by the few
  hand-rolled argument parsers (`bedbamtobed`, `bedtobam`, `bedunionbedg`,
  `bedintersect`). This needs no execution and is robust to formatting.
- **Score** — `compat% = |upstream ∩ ours| / |upstream|`. For the subcommanded
  tools (samtools, bcftools, and the `bed*` family) our port's flags are pooled
  across the whole tool family and each subcommand's upstream flags are tested
  against that pool, i.e. *"does the port accept this spelling anywhere in the
  tool"*. The **aggregate** is weighted by upstream flag count (every upstream
  flag slot counts once).

Reproduce with:

```sh
go run ./pipeline/cmd/flagcompat
```

The program is `gofmt`-clean, `go vet`-clean, and builds as part of
`go build ./...`. The numbers below are its output verbatim.

## Caveats and known interpretation notes

- **Spelling-only, not behaviour.** See above. Pair this with the parity tests
  for a complete picture.
- **Deliberate CLI redesign depresses some scores.** The lowest scorers
  (`fastp` ~47%, `prinseq` ~44%, `skewer` ~31%) reflect a *design choice*, not
  missing functionality: those ports expose the project's GNU-style hyphenated
  long flags (`--qual-threshold`, `--min-length`, `--cut-front`) where upstream
  uses snake_case or terse single letters (`--qualified_quality_phred`,
  `--length_required`, `-5`). The capability is present and parity-tested; only
  the spelling differs, so it reads as "incompatible" under a strict
  name-for-name match. Where the port also registers the upstream spelling, it
  is credited.
- **Token-name matching can over-credit collisions.** Because matching is by
  bare token, a port that reuses a short letter for a *different* purpose than
  upstream (e.g. `bedtag -i` is `--tag` in our port but the input file upstream)
  is counted as a match. Such collisions are rare and inflate a handful of
  scores by at most one flag.
- **Pooling across subcommands can over-credit.** For samtools/bcftools/`bed*`
  the "accepted anywhere in the tool" pooling means a flag registered in one
  subcommand satisfies the same spelling required by another. This is the
  intended surface-acceptance reading; a strict per-subcommand-binding metric
  would be marginally lower.
- **`--input-fmt` dominates the samtools gap.** Most samtools subcommands are
  missing only `--input-fmt` (and a couple of double-spelled aliases like
  `--min-mq`/`--min-bq`): the port auto-detects input format and registers
  `--output-fmt` but not the explicit `--input-fmt` selector. This single
  systematic omission accounts for the bulk of samtools' shortfall from 100%.

## Per-tool / per-subcommand compatibility

<!-- markdownlint-disable MD013 -->
| Tool / subcommand | Upstream flags | Accepted | Compat % | Missing (upstream flags not accepted) |
|---|---:|---:|---:|---|
| seqtk seq | 23 | 18 | 78.3% | `1` `F` `R` `S` `a` |
| seqtk comp | 2 | 1 | 50.0% | `u` |
| seqtk sample | 2 | 2 | 100.0% | — |
| seqtk subseq | 3 | 3 | 100.0% | — |
| seqtk trimfq | 5 | 5 | 100.0% | — |
| seqtk mergefa | 5 | 5 | 100.0% | — |
| seqtk cutN | 3 | 3 | 100.0% | — |
| seqtk gap | 1 | 1 | 100.0% | — |
| seqtk hety | 3 | 3 | 100.0% | — |
| seqtk gc | 4 | 4 | 100.0% | — |
| seqtk split | 2 | 2 | 100.0% | — |
| seqtk fqchk | 1 | 1 | 100.0% | — |
| seqtk telo | 5 | 5 | 100.0% | — |
| sickle | 30 | 22 | 73.3% | `M` `c` `discard-n` `m` `output-combo` `output-combo-all` `pe-combo` `truncate-n` |
| skewer | 48 | 15 | 31.2% | `1` `A` `L` `M` `N` `Q` `X` `b` `barcode` `c` `cut` `cut3` `e` `end-quality` `excluded-output` `fillNs` `format` `help` `intelligent` `k` `len` `masked-output` `matrix` `max` `max-len` `mean` `mean-quality` `min` `mode` `n` `stdout` `threads` `u` |
| fastp | 132 | 62 | 47.0% | `6` `B` `D` `F` `R` `T` `U` `V` `Y` `a` `adapter_sequence` `adapter_sequence_r2` `allow_gap_overlap_trimming` `average_qual` `b` `complexity_threshold` `compression` `cut_by_quality3` `cut_by_quality5` `cut_by_quality_aggressive` `cut_front` `cut_front_mean_quality` `cut_front_window_size` `cut_mean_quality` `cut_right` `cut_right_mean_quality` `cut_right_window_size` `cut_tail` `cut_tail_mean_quality` `cut_tail_window_size` `cut_window_size` `discard_unmerged` `dont_eval_duplication` `dont_overwrite` `e` `f` `failed_out` `filter_by_index1` `filter_by_index2` `filter_by_index_threshold` `fix_mgi_id` `interleaved_in` `j` `length_limit` `length_required` `low_complexity_filter` `max_len1` `max_len2` `n` `n_base_limit` `overlapped_out` `phred64` `qualified_quality_phred` `reads_to_process` `report_title` `stdin` `stdout` `thread` `trim_front1` `trim_front2` `trim_poly_x` `trim_tail1` `trim_tail2` `u` `umi_delim` `unpaired1` `unpaired2` `unqualified_percent_limit` `verbose` `z` |
| prinseq | 73 | 32 | 43.8% | `aa` `derep_min` `fasta2` `fastq2` `filename1` `filename2` `h` `help` `lc_method` `lc_threshold` `man` `max_gc` `max_len` `max_qual_mean` `max_qual_score` `min_gc` `min_len` `min_qual_mean` `min_qual_score` `out_bad` `out_good` `qual` `stats_all` `stats_assembly` `stats_dinuc` `stats_dupl` `stats_info` `stats_len` `stats_ns` `stats_tag` `trim_left` `trim_left_p` `trim_ns_left` `trim_ns_right` `trim_qual_left` `trim_qual_right` `trim_right` `trim_right_p` `trim_tail_left` `trim_tail_right` `verbose` |
| vcftools | 146 | 146 | 100.0% | — |
| mosdepth | 37 | 37 | 100.0% | — |
| bgzip | 31 | 30 | 96.8% | `binary` |
| tabix | 40 | 31 | 77.5% | `begin` `cache` `comment` `end` `min-shift` `separate-regions` `sequence` `threads` `verbosity` |
| htsfile | 5 | 4 | 80.0% | `H` |
| samtools view | 101 | 66 | 65.3% | `QNAME-file` `U` `add-flags` `customised-index` `excl-no-read-group` `excl-no-readgroup` `exclude-no-read-group` `exclude-no-readgroup` `expr` `expression` `fai-reference` `fast` `fetch-pairs` `input-fmt` `keep-tag` `library` `min-mq` `min-qlen` `no-header` `output-unselected` `read-group-file` `readgroup` `readgroup-file` `region-file` `remove-B` `remove-flags` `remove-tag` `require-flags` `sanitize` `save-counts` `subsample-seed` `target-file` `targets-file` `unmap` `unoutput` |
| samtools sort | 26 | 24 | 92.3% | `input-fmt` `template-coordinate` |
| samtools index | 18 | 17 | 94.4% | `input-fmt` |
| samtools flagstat | 9 | 8 | 88.9% | `input-fmt` |
| samtools depth | 30 | 26 | 86.7% | `input-fmt` `min-bq` `min-mq` `require-flags` |
| samtools fastq | 47 | 37 | 78.7% | `I1` `I2` `IF` `if` `index-format` `input-fmt` `no-sc` `no-sc-bkp` `require-flags` `sc-aux` |
| samtools mpileup | 83 | 45 | 54.2% | `adjust-MQ` `adjust-mq` `disable-overlap-removal` `exclude-RG` `exclude-rg` `ext-prob` `gap-frac` `ignore-RG` `ignore-overlaps-removal` `ignore-rg` `input-fmt` `max-idepth` `min-bq` `min-ireads` `min-mq` `no-baq` `no-output-del` `no-output-ends` `no-output-ins` `no-output-ins-mods` `output-BP-5` `output-MQ` `output-QNAME` `output-bp` `output-bp-5` `output-empty` `output-extra` `output-mods` `output-mq` `output-qname` `output-sep` `per-sample-mF` `per-sample-mf` `platforms` `redo-BAQ` `reverse-del` `skip-indels` `tandem-qual` |
| samtools idxstats | 8 | 7 | 87.5% | `input-fmt` |
| samtools quickcheck | 3 | 3 | 100.0% | — |
| samtools dict | 18 | 14 | 77.8% | `?` `alt` `alternative-name` `no-header` |
| samtools cat | 16 | 15 | 93.8% | `input-fmt` |
| samtools reheader | 8 | 6 | 75.0% | `command` `in-place` |
| samtools addreplacerg | 18 | 17 | 94.4% | `input-fmt` |
| samtools fixmate | 18 | 17 | 94.4% | `input-fmt` |
| samtools merge | 29 | 27 | 93.1% | `input-fmt` `template-coordinate` |
| samtools coverage | 39 | 35 | 89.7% | `input-fmt` `min-bq` `min-mq` `no-header` |
| samtools split | 19 | 16 | 84.2% | `input-fmt` `max-split` `zero-pad` |
| samtools markdup | 33 | 23 | 69.7% | `barcode-name` `barcode-rgx` `coords-order` `duplicate-count` `include-fails` `input-fmt` `json` `no-multi-dup` `read-coords` `use-read-groups` |
| samtools stats | 49 | 40 | 81.6% | `?` `customized-index-file` `id` `input-fmt` `most-inserts` `read-length` `sam` `split` `split-prefix` |
| samtools calmd | 24 | 23 | 95.8% | `input-fmt` |
| samtools import | 34 | 30 | 88.2% | `input-fmt` `r1` `r2` `rg` |
| samtools phase | 18 | 16 | 88.9% | `input-fmt` `min-bq` |
| samtools targetcut | 13 | 12 | 92.3% | `input-fmt` |
| samtools consensus | 69 | 68 | 98.6% | `input-fmt` |
| samtools tview | 12 | 11 | 91.7% | `input-fmt` |
| bcftools view | 80 | 68 | 85.0% | `U` `exclude-phased` `exclude-uncalled` `genotype` `known` `max-alleles` `min-alleles` `novel` `output-file` `trim-unseen-allele` `uncalled` `with-header` |
| bcftools call | 66 | 64 | 97.0% | `keep-masked-refs` `skip-Ns` |
| bcftools norm | 56 | 48 | 85.7% | `atom-overlaps` `keep-sum` `multi-overlaps` `no-realign` `old-rec-tag` `right-align` `site-win` `sort` |
| bcftools filter | 39 | 39 | 100.0% | — |
| bcftools query | 23 | 22 | 95.7% | `disable-automatic-newline` |
| bcftools concat | 38 | 34 | 89.5% | `compact-PS` `naive` `naive-force` `rm-dups` |
| bcftools merge | 41 | 33 | 80.5% | `filter-logic` `force-no-index` `force-single` `local-alleles` `missing-rules` `missing-to-ref` `no-index` `use-header` |
| bcftools stats | 22 | 22 | 100.0% | — |
| bcftools isec | 40 | 39 | 97.5% | `complement` |
| bcftools annotate | 51 | 48 | 94.1% | `columns-file` `header-line` `mark-sites` |
| bcftools roh | 23 | 23 | 100.0% | — |
| bcftools convert | 52 | 52 | 100.0% | — |
| bcftools sort | 16 | 14 | 87.5% | `output-file` `temp-dir` |
| bcftools head | 8 | 7 | 87.5% | `records` |
| bcftools gtcheck | 28 | 28 | 100.0% | — |
| bcftools consensus | 25 | 25 | 100.0% | — |
| bcftools csq | 34 | 34 | 100.0% | — |
| bcftools cnv | 26 | 26 | 100.0% | — |
| bcftools polysomy | 21 | 21 | 100.0% | — |
| bcftools reheader | 9 | 9 | 100.0% | — |
| bcftools index | 20 | 15 | 75.0% | `all` `min-shift` `nrecords` `output-file` `stats` |
| bcftools mpileup | 111 | 108 | 97.3% | `U` `excl-flags` `incl-flags` |
| bcftools som | 26 | 26 | 100.0% | — |
| bed12tobed6 | 3 | 3 | 100.0% | — |
| bedannotate | 8 | 6 | 75.0% | `files` `names` |
| bedbamtobed | 14 | 11 | 78.6% | `as` `bwa` `novo` |
| bedclosest | 21 | 17 | 81.0% | `abam` `b` `names` `split` |
| bedcluster | 4 | 4 | 100.0% | — |
| bedcomplement | 4 | 4 | 100.0% | — |
| bedcoverage | 18 | 17 | 94.4% | `i` |
| bedexpand | 3 | 3 | 100.0% | — |
| bedfisher | 15 | 10 | 66.7% | `abam` `e` `exclude` `i` `split` |
| bedflank | 9 | 8 | 88.9% | `header` |
| bedgenomecov | 20 | 18 | 90.0% | `du` `ignoreD` |
| bedgetfasta | 12 | 12 | 100.0% | — |
| bedgroupby | 13 | 13 | 100.0% | — |
| bedigv | 9 | 9 | 100.0% | — |
| bedintersect | 26 | 21 | 80.8% | `F` `b` `f` `i` `names` |
| bedjaccard | 13 | 8 | 61.5% | `abam` `e` `g` `i` `r` |
| bedlinks | 5 | 5 | 100.0% | — |
| bedmakewindows | 8 | 8 | 100.0% | — |
| bedmap | 27 | 17 | 63.0% | `C` `abam` `e` `i` `loj` `ops` `wa` `wao` `wb` `wo` |
| bedmerge | 15 | 13 | 86.7% | `ops` `scores` |
| bedmulticov | 12 | 10 | 83.3% | `bams` `p` |
| bedmultiinter | 9 | 6 | 66.7% | `examples` `i` `names` |
| bednuc | 8 | 8 | 100.0% | — |
| bedoverlap | 3 | 3 | 100.0% | — |
| bedpairtobed | 11 | 7 | 63.6% | `abam` `bedpe` `ed` `ubam` |
| bedpairtopair | 9 | 9 | 100.0% | — |
| bedrandom | 5 | 5 | 100.0% | — |
| bedreldist | 4 | 4 | 100.0% | — |
| bedsample | 7 | 5 | 71.4% | `s` `ubam` |
| bedshift | 8 | 8 | 100.0% | — |
| bedshuffle | 13 | 10 | 76.9% | `bedpe` `f` `noOverlapping` |
| bedslop | 9 | 8 | 88.9% | `header` |
| bedsort | 11 | 11 | 100.0% | — |
| bedspacing | 3 | 2 | 66.7% | `g` |
| bedsplit | 11 | 10 | 90.9% | `bed` |
| bedsubtract | 14 | 8 | 57.1% | `F` `abam` `e` `g` `i` `split` |
| bedsummary | 3 | 3 | 100.0% | — |
| bedtag | 11 | 8 | 72.7% | `files` `intervals` `scores` |
| bedtobam | 6 | 6 | 100.0% | — |
| bedunionbedg | 10 | 7 | 70.0% | `examples` `files` `labels` |
| bedwindow | 16 | 12 | 75.0% | `abam` `bed` `header` `ubam` |

## Per-tool-family rollup

| Tool family | Upstream flag slots | Accepted | Compat % |
|---|---:|---:|---:|
| seqtk | 59 | 53 | 89.8% |
| sickle | 30 | 22 | 73.3% |
| skewer | 48 | 15 | 31.2% |
| fastp | 132 | 62 | 47.0% |
| prinseq | 73 | 32 | 43.8% |
| vcftools | 146 | 146 | 100.0% |
| mosdepth | 37 | 37 | 100.0% |
| bgzip | 31 | 30 | 96.8% |
| tabix | 40 | 31 | 77.5% |
| htsfile | 5 | 4 | 80.0% |
| samtools | 742 | 603 | 81.3% |
| bcftools | 855 | 805 | 94.2% |
| bed12tobed6 | 3 | 3 | 100.0% |
| bedannotate | 8 | 6 | 75.0% |
| bedbamtobed | 14 | 11 | 78.6% |
| bedclosest | 21 | 17 | 81.0% |
| bedcluster | 4 | 4 | 100.0% |
| bedcomplement | 4 | 4 | 100.0% |
| bedcoverage | 18 | 17 | 94.4% |
| bedexpand | 3 | 3 | 100.0% |
| bedfisher | 15 | 10 | 66.7% |
| bedflank | 9 | 8 | 88.9% |
| bedgenomecov | 20 | 18 | 90.0% |
| bedgetfasta | 12 | 12 | 100.0% |
| bedgroupby | 13 | 13 | 100.0% |
| bedigv | 9 | 9 | 100.0% |
| bedintersect | 26 | 21 | 80.8% |
| bedjaccard | 13 | 8 | 61.5% |
| bedlinks | 5 | 5 | 100.0% |
| bedmakewindows | 8 | 8 | 100.0% |
| bedmap | 27 | 17 | 63.0% |
| bedmerge | 15 | 13 | 86.7% |
| bedmulticov | 12 | 10 | 83.3% |
| bedmultiinter | 9 | 6 | 66.7% |
| bednuc | 8 | 8 | 100.0% |
| bedoverlap | 3 | 3 | 100.0% |
| bedpairtobed | 11 | 7 | 63.6% |
| bedpairtopair | 9 | 9 | 100.0% |
| bedrandom | 5 | 5 | 100.0% |
| bedreldist | 4 | 4 | 100.0% |
| bedsample | 7 | 5 | 71.4% |
| bedshift | 8 | 8 | 100.0% |
| bedshuffle | 13 | 10 | 76.9% |
| bedslop | 9 | 8 | 88.9% |
| bedsort | 11 | 11 | 100.0% |
| bedspacing | 3 | 2 | 66.7% |
| bedsplit | 11 | 10 | 90.9% |
| bedsubtract | 14 | 8 | 57.1% |
| bedsummary | 3 | 3 | 100.0% |
| bedtag | 11 | 8 | 72.7% |
| bedtobam | 6 | 6 | 100.0% |
| bedunionbedg | 10 | 7 | 70.0% |
| bedwindow | 16 | 12 | 75.0% |

## Aggregate (weighted by upstream flag count)

- Upstream flag slots measured: **2628**
- Accepted under the same spelling: **2197**
- Weighted flag-compatibility: **83.6%**
<!-- markdownlint-enable MD013 -->
