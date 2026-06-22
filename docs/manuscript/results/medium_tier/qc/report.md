# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-22T13:22:21Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 102 | 90 | 0 | 2 | 10 | 0 |

## bgzip

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bgzip_decompress | vcf | ByteExact | PASS | 21 | 67 |  |  |
| bgzip_decompress_heavy | vcf | ByteExact | PASS | 16 | 16 | 1.03x |  |
| bgzip_compress | vcf_plain | BGZFDecoded | PASS | 29 | 100 |  |  |
| bgzip_reindex | vcf | ByteExact | PASS | 4 | 11 |  |  |

## fastp

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| fastp_cut_tail | fastq | ByteExact | PASS | 2998 | 797 |  |  |
| fastp_default_filter | fastq | ByteExact | PASS | 3537 | 557 |  |  |
| fastp_detect_adapter_pe_heavy | fastq_paired | ByteExact | PASS | 24273 | 3259 | 7.45x |  |

## htsfile

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| htsfile_identify | vcf | ByteExact | PASS | 11 | 4 |  |  |
| htsfile_copy | vcf | ByteExact | PASS | 81 | 43 |  |  |

## mosdepth

PASS 0 · SIMILAR 0 · DIVERGE 0 · SKIP 10 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| mosdepth_default | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_fast_mode | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_mapq20 | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_flag | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_by_bed_regions | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_by_window_regions | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_by_bed_thresholds | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_default_heavy | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_by_summary_region_rows | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |
| mosdepth_by_region_dist | bam | ByteExact | SKIP | 0 | 0 |  | upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosd... |

## prinseq

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| prinseq_min_len | fastq | ByteExact | PASS | 1510 | 2648 |  |  |
| prinseq_max_len | fastq | ByteExact | PASS | 1361 | 2874 |  |  |
| prinseq_trim_left | fastq | ByteExact | PASS | 1571 | 2664 |  |  |
| prinseq_trim_right | fastq | ByteExact | PASS | 1494 | 2792 |  |  |
| prinseq_trim_qual_right | fastq | ByteExact | PASS | 1640 | 7871 |  |  |
| prinseq_trim_qual_left | fastq | ByteExact | PASS | 1519 | 6819 |  |  |
| prinseq_min_qual_mean | fastq | ByteExact | PASS | 1426 | 6890 |  |  |
| prinseq_max_ns | fastq | ByteExact | PASS | 3192 | 2519 |  |  |

## seqtk

PASS 32 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| seqtk_seq_fq_base | fastq | ByteExact | PASS | 225 | 163 |  |  |
| seqtk_seq_fq_flagA | fastq | ByteExact | PASS | 127 | 92 |  |  |
| seqtk_seq_fq_flagr | fastq | ByteExact | PASS | 227 | 185 |  |  |
| seqtk_seq_fq_flagL_95 | fastq | ByteExact | PASS | 154 | 112 |  |  |
| seqtk_seq_fq_flagq_20 | fastq | ByteExact | PASS | 162 | 170 |  |  |
| seqtk_seq_fq_flagl_60 | fastq | ByteExact | PASS | 159 | 118 |  |  |
| seqtk_seq_fq_combo_rev_upper | fastq | ByteExact | PASS | 210 | 214 |  |  |
| seqtk_seq_fq_combo_qmask_n | fastq | ByteExact | PASS | 186 | 153 |  |  |
| seqtk_comp_fa | fasta | ByteExact | PASS | 133 | 81 |  |  |
| seqtk_comp_fq | fastq | ByteExact | PASS | 564 | 553 |  |  |
| seqtk_fqchk | fastq | ByteExact | PASS | 146 | 82 |  |  |
| seqtk_fqchk_q20 | fastq | ByteExact | PASS | 165 | 85 |  |  |
| seqtk_size_fq | fastq | ByteExact | PASS | 77 | 40 |  |  |
| seqtk_size_fa | fasta | ByteExact | PASS | 29 | 9 |  |  |
| seqtk_trimfq | fastq | ByteExact | PASS | 289 | 182 |  |  |
| seqtk_trimfq_q | fastq | ByteExact | PASS | 284 | 188 |  |  |
| seqtk_trimfq_be | fastq | ByteExact | PASS | 156 | 108 |  |  |
| seqtk_sample_count | fastq | ByteExact | PASS | 97 | 44 |  |  |
| seqtk_sample_frac | fastq | ByteExact | PASS | 110 | 56 |  |  |
| seqtk_hpc_fa | fasta | ByteExact | PASS | 78 | 54 |  |  |
| seqtk_hpc_fq | fastq | ByteExact | PASS | 259 | 203 |  |  |
| seqtk_gap_fa | fasta | ByteExact | PASS | 48 | 19 |  |  |
| seqtk_subseq_bed | fasta | ByteExact | PASS | 54 | 160 |  |  |
| seqtk_mergepe | fastq_paired | ByteExact | PASS | 382 | 326 |  |  |
| seqtk_dropse | fastq | ByteExact | PASS | 89 | 46 |  |  |
| seqtk_randbase | fasta | ByteExact | PASS | 134 | 60 |  |  |
| seqtk_telo | fasta | ByteExact | PASS | 45 | 10 |  |  |
| seqtk_listhet | fasta | ByteExact | PASS | 49 | 16 |  |  |
| seqtk_hety | fasta | ByteExact | PASS | 102 | 80 |  |  |
| seqtk_seq_fa | fasta | ByteExact | PASS | 59 | 28 |  |  |
| seqtk_seq_fq_to_fa_heavy | fastq | ByteExact | PASS | 123 | 93 | 1.31x |  |
| seqtk_cutN | fasta | ByteExact | PASS | 57 | 71 |  |  |

## sickle

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| sickle_se_base | fastq | ByteExact | PASS | 1348 | 1203 |  |  |
| sickle_se_q30 | fastq | ByteExact | PASS | 1295 | 1225 |  |  |
| sickle_se_l30 | fastq | ByteExact | PASS | 1358 | 1170 |  |  |
| sickle_se_no5prime | fastq | ByteExact | PASS | 1158 | 1256 |  |  |
| sickle_se_truncn | fastq | ByteExact | PASS | 1149 | 1151 |  |  |
| sickle_se_q30_l40 | fastq | ByteExact | PASS | 1085 | 1109 |  |  |
| sickle_pe_base | fastq_paired | ByteExact | PASS | 2209 | 2235 | 0.99x |  |
| sickle_se_cli_default_window | fastq | ByteExact | PASS | 1194 | 1158 |  |  |

## skewer

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| skewer_se_base | fastq | ByteExact | PASS | 2383 | 2078 |  |  |
| skewer_se_minlen30 | fastq | ByteExact | PASS | 2116 | 1937 |  |  |
| skewer_se_full_heavy | fastq | ByteExact | PASS | 2101 | 1934 | 1.09x |  |

## tabix

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| tabix_region_contig | vcf | ByteExact | PASS | 51 | 7 |  |  |
| tabix_region_range | vcf | ByteExact | PASS | 4 | 3 |  |  |
| tabix_region_chr2 | vcf | ByteExact | PASS | 36 | 6 |  |  |
| tabix_region_with_header | vcf | ByteExact | PASS | 36 | 6 |  |  |
| tabix_list_chroms | vcf | ByteExact | PASS | 5 | 4 |  |  |
| tabix_region_heavy | vcf | ByteExact | PASS | 38 | 6 | 6.36x |  |
| tabix_regions_bed | vcf | ByteExact | PASS | 42141 | 266 |  |  |

## vcftools

PASS 23 · SIMILAR 0 · DIVERGE 2 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| vcftools_freq | vcf_plain | ByteExact | PASS | 2339 | 164 |  |  |
| vcftools_counts | vcf_plain | ByteExact | PASS | 2088 | 78 |  |  |
| vcftools_freq2 | vcf_plain | ByteExact | PASS | 2265 | 89 |  |  |
| vcftools_depth | vcf_plain | ByteExact | PASS | 53 | 44 |  |  |
| vcftools_site_depth | vcf_plain | ByteExact | PASS | 1618 | 1517 |  |  |
| vcftools_site_mean_depth | vcf_plain | ByteExact | DIVERGE | 1995 | 1554 |  | output file ".ldepth.mean": first diff at line 2:   ours:     chr1	15	56	-nan   upstream: chr1	15	56	nan |
| vcftools_site_pi | vcf_plain | ByteExact | PASS | 2417 | 1558 |  |  |
| vcftools_window_pi | vcf_plain | ByteExact | PASS | 590 | 333 |  |  |
| vcftools_tstv_summary | vcf_plain | ByteExact | PASS | 54 | 35 |  |  |
| vcftools_missing_indv | vcf_plain | ByteExact | PASS | 54 | 43 |  |  |
| vcftools_missing_site | vcf_plain | ByteExact | PASS | 1522 | 1556 |  |  |
| vcftools_het | vcf_plain | ByteExact | PASS | 66 | 56 |  |  |
| vcftools_singletons | vcf_plain | ByteExact | PASS | 2333 | 2021 |  |  |
| vcftools_recode_heavy | vcf_plain | ByteExact | PASS | 105 | 172 |  |  |
| vcftools_het_multi | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_relatedness | vcf_multi_plain | ByteExact | PASS | 5 | 1 |  |  |
| vcftools_relatedness2 | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_freq_multi | vcf_multi_plain | ByteExact | PASS | 5 | 1 |  |  |
| vcftools_missing_indv_multi | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_window_pi_heavy | vcf_plain | ByteExact | PASS | 1272 | 783 | 1.62x |  |
| vcftools_geno_r2 | vcf_multi_plain | ByteExact | PASS | 4 | 1 | 2.84x |  |
| vcftools_hap_r2 | vcf_multi_plain | ByteExact | PASS | 1075 | 18809 | 0.06x |  |
| vcftools_matrix012 | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_lroh | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_hardy | vcf_plain | ByteExact | DIVERGE | 2392 | 1665 |  | output file ".hwe": first diff at line 2:   ours:     chr1	15	0/0/1	0.00/0.00/1.00	-nan	1.000000e+00	1.000000e+00	1.0000... |

