# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-23T12:05:29Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 400 | 383 | 4 | 2 | 11 | 0 |

## bcftools

PASS 66 · SIMILAR 1 · DIVERGE 0 · SKIP 1 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_base | vcf_plain | ByteExact | PASS | 53 | 43 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 46 | 42 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 29 | 29 |  |  |
| view_flagV_indels | vcf_plain | ByteExact | PASS | 53 | 38 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 52 | 36 |  |  |
| view_flagi_DP>40 | vcf_plain | ByteExact | PASS | 5 | 4 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 52 | 39 |  |  |
| view_flagO_v | vcf_plain | ByteExact | PASS | 53 | 44 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 47 | 40 |  |  |
| view_combo_indels_dp | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| view_combo_excl_snps_qual | vcf_plain | ByteExact | PASS | 29 | 29 |  |  |
| bcftools_view_r_region | vcf | ByteExact | PASS | 8 | 5 |  |  |
| bcftools_view_t_targets | vcf | ByteExact | PASS | 35 | 35 |  |  |
| bcftools_view_r_multi | vcf | ByteExact | PASS | 178 | 11 |  |  |
| bcftools_view_s_sample | vcf_multi_plain | ByteExact | PASS | 306 | 95 |  |  |
| bcftools_view_G_drop_gt | vcf_multi_plain | ByteExact | PASS | 249 | 79 |  |  |
| bcftools_view_c_minac | vcf_multi_plain | ByteExact | PASS | 230 | 120 |  |  |
| bcftools_view_C_maxac | vcf_multi_plain | ByteExact | PASS | 166 | 74 |  |  |
| bcftools_view_q_minaf | vcf_multi_plain | ByteExact | PASS | 248 | 112 |  |  |
| bcftools_view_Q_maxaf | vcf_multi_plain | ByteExact | PASS | 241 | 107 |  |  |
| bcftools_query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 33 | 21 |  |  |
| bcftools_query_info_dp | vcf_plain | ByteExact | PASS | 30 | 27 |  |  |
| bcftools_query_qual_filter | vcf_plain | ByteExact | PASS | 34 | 24 |  |  |
| bcftools_query_info_dp_filtered | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_info_dp_excluded | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_gt_multi | vcf_multi_plain | ByteExact | PASS | 147 | 113 |  |  |
| bcftools_query_list_samples | vcf_multi_plain | ByteExact | PASS | 6 | 4 |  |  |
| bcftools_query_heavy_full | vcf_plain | ByteExact | PASS | 37 | 31 | 1.17x |  |
| bcftools_norm_split | vcf_plain | ByteExact | PASS | 99 | 66 |  |  |
| bcftools_norm_dedup | vcf_plain | ByteExact | PASS | 77 | 57 |  |  |
| bcftools_norm_check_ref | vcf_plain | ByteExact | PASS | 95 | 66 |  |  |
| filter_base | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_flagi_QUAL>=30 | vcf_plain | ByteExact | PASS | 5 | 2 |  |  |
| filter_combo_soft_lowqual | vcf_plain | ByteExact | PASS | 4 | 3 |  |  |
| filter_combo_snpgap | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_sort | vcf_plain | ByteExact | PASS | 69 | 78 |  |  |
| bcftools_stats | vcf_plain | ByteExact | PASS | 28 | 26 |  |  |
| bcftools_stats_multi | vcf_multi_plain | ByteExact | PASS | 105 | 41 |  |  |
| bcftools_head | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_annotate_drop_af | vcf_plain | ByteExact | PASS | 75 | 42 |  |  |
| bcftools_annotate_drop_dp | vcf_plain | ByteExact | PASS | 71 | 44 |  |  |
| bcftools_gtcheck | vcf_multi_plain | ByteExact | PASS | 1294 | 125 |  |  |
| bcftools_mpileup | bam | ByteExact | PASS | 80824 | 28961 |  |  |
| bcftools_mpileup_heavy | bam | ByteExact | PASS | 72924 | 22683 | 3.21x |  |
| bcftools_plugin_fill_tags | vcf_multi_plain | ByteExact | PASS | 719 | 176 |  |  |
| bcftools_plugin_fill_tags_an_ac | vcf_multi_plain | ByteExact | PASS | 577 | 186 |  |  |
| bcftools_view_c_update_acan | vcf_multi_plain | ByteExact | PASS | 500 | 110 |  |  |
| bcftools_annotate_drop_filter | vcf_plain | ByteExact | PASS | 100 | 41 |  |  |
| bcftools_roh | vcf_multi_plain | ByteExact | PASS | 1004 | 376 |  |  |
| bcftools_consensus | vcf | ByteExact | PASS | 6857 | 91 |  |  |
| bcftools_concat | vcf | ByteExact | PASS | 190 | 118 |  |  |
| bcftools_norm_join | vcf_plain | ByteExact | PASS | 131 | 80 |  |  |
| bcftools_call | vcf_plain | Similarity | SIMILAR | 39344 | 24454 |  |  |
| bcftools_csq | vcf_plain | ByteExact | SKIP | 0 | 0 |  | bcftools csq --force over the multi-gene fixture is byte-identical to upstream EXCEPT for 3 records, and that residual i... |
| bcftools_isec | vcf | ByteExact | PASS | 4449 | 294 |  |  |
| bcftools_merge | vcf_plain | ByteExact | PASS | 335 | 149 |  |  |
| bcftools_convert | vcf_plain | ByteExact | PASS | 70 | 40 |  |  |
| bcftools_reheader | vcf_multi_plain | BGZFDecoded | PASS | 417 | 13 |  |  |
| view_base | vcf_plain | ByteExact | PASS | 57 | 45 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 45 | 37 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 30 | 29 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 48 | 36 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 46 | 36 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 45 | 35 |  |  |
| query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 32 | 21 |  |  |
| query_gt | vcf_plain | ByteExact | PASS | 30 | 34 |  |  |
| query_info_dp | vcf_plain | ByteExact | PASS | 31 | 27 |  |  |

## bed12tobed6

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bed12tobed6_score_dropped | bed12 | ByteExact | PASS | 40 | 79 |  |  |

## bedannotate

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedannotate_default_header_order | bed | ByteExact | PASS | 162 | 461 |  |  |

## bedbamtobed

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bam | ByteExact | PASS | 650 | 576 |  |  |
| flagsplit | bam | ByteExact | PASS | 664 | 708 |  |  |
| flaged | bam | ByteExact | PASS | 14 | 2 |  |  |
| flagcigar | bam | ByteExact | PASS | 675 | 649 |  |  |
| combo_split_cigar | bam | ByteExact | PASS | 655 | 671 |  |  |

## bedclosest

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 2938 | 141 |  |  |
| flagd | bed | ByteExact | PASS | 2790 | 180 |  |  |
| flagio | bed | ByteExact | PASS | 2862 | 60 |  |  |
| flagiu | bed | ByteExact | PASS | 5 | 2 |  |  |
| flagt_first | bed | ByteExact | PASS | 2672 | 50 |  |  |
| flagt_last | bed | ByteExact | PASS | 2712 | 47 |  |  |
| flagt_all | bed | ByteExact | PASS | 2749 | 125 |  |  |
| flags | bed | ByteExact | PASS | 3033 | 100 |  |  |
| flagN | bed | ByteExact | PASS | 3209 | 96 |  |  |
| combo_d_t_first | bed | ByteExact | PASS | 2714 | 58 |  |  |

## bedcluster

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 30 | 29 |  |  |
| flagd_50 | bed | ByteExact | PASS | 23 | 30 |  |  |
| flagd_0 | bed | ByteExact | PASS | 22 | 30 |  |  |
| bedcluster_s | bed | ByteExact | PASS | 31 | 44 |  |  |
| bedcluster_d50_s | bed | ByteExact | PASS | 34 | 42 |  |  |

## bedcomplement

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedcomplement_base | bed | ByteExact | PASS | 22 | 12 |  |  |

## bedcoverage

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 131 | 164 |  |  |
| flagcounts | bed | ByteExact | PASS | 101 | 140 |  |  |
| flagd | bed | ByteExact | PASS | 4427 | 19711 |  |  |
| flaghist | bed | ByteExact | PASS | 613 | 423 |  |  |
| flags | bed | ByteExact | PASS | 113 | 142 |  |  |
| flagS | bed | ByteExact | PASS | 89 | 127 |  |  |
| flagmean | bed | ByteExact | PASS | 120 | 174 |  |  |
| combo_counts_s | bed | ByteExact | PASS | 84 | 115 |  |  |

## bedexpand

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedexpand_c5 | bed | ByteExact | PASS | 27 | 36 |  |  |
| bedexpand_trailing_comma_col11 | bed12 | ByteExact | PASS | 43 | 64 |  |  |

## bedfisher

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedfisher_overlap_count | bed | ByteExact | PASS | 109 | 27 |  |  |

## bedflank

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 10 | 1 |  |  |
| flagb_50 | bed | ByteExact | PASS | 30 | 44 |  |  |
| flags | bed | ByteExact | PASS | 4 | 1 |  |  |
| combo_l_r | bed | ByteExact | PASS | 29 | 44 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 27 | 44 |  |  |
| combo_b_s | bed | ByteExact | PASS | 31 | 42 |  |  |

## bedgenomecov

PASS 5 · SIMILAR 3 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | Similarity | SIMILAR | 296 | 246 |  |  |
| flagbg | bed | ByteExact | PASS | 75 | 81 |  |  |
| flagbga | bed | ByteExact | PASS | 83 | 82 |  |  |
| flagd | bed | ByteExact | PASS | 861 | 7990 |  |  |
| flagdz | bed | ByteExact | PASS | 693 | 7953 |  |  |
| flagmax_5 | bed | Similarity | SIMILAR | 274 | 308 |  |  |
| flagstrand_+ | bed | Similarity | SIMILAR | 263 | 259 |  |  |
| combo_bg_strand | bed | ByteExact | PASS | 66 | 63 |  |  |

## bedgetfasta

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | fasta | ByteExact | PASS | 171 | 433 |  |  |
| flags | fasta | ByteExact | PASS | 322 | 573 |  |  |
| flagname | fasta | ByteExact | PASS | 140 | 415 |  |  |
| flagtab | fasta | ByteExact | PASS | 136 | 392 |  |  |
| combo_s_name | fasta | ByteExact | PASS | 318 | 587 |  |  |
| combo_s_tab | fasta | ByteExact | PASS | 414 | 585 |  |  |

## bedgroupby

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 17 | 1 |  |  |
| flago_mean | bed | ByteExact | PASS | 13 | 11 |  |  |
| flago_sum | bed | ByteExact | PASS | 11 | 11 |  |  |
| flago_min | bed | ByteExact | PASS | 12 | 11 |  |  |
| flago_max | bed | ByteExact | PASS | 11 | 11 |  |  |
| flago_count | bed | ByteExact | PASS | 10 | 9 |  |  |
| combo_g1_c5_mean | bed | ByteExact | PASS | 13 | 11 |  |  |
| combo_g1_c5_count | bed | ByteExact | PASS | 11 | 11 |  |  |

## bedigv

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedigv_base | bed | ByteExact | PASS | 40 | 110 |  |  |

## bedintersect

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| intersect_base | bed | ByteExact | PASS | 109 | 125 |  |  |
| intersect_flagc | bed | ByteExact | PASS | 54 | 73 |  |  |
| intersect_flagv | bed | ByteExact | PASS | 47 | 61 |  |  |
| intersect_flagu | bed | ByteExact | PASS | 51 | 69 |  |  |
| intersect_flagwa | bed | ByteExact | PASS | 95 | 127 |  |  |
| intersect_flagwb | bed | ByteExact | PASS | 165 | 216 |  |  |
| intersect_flags | bed | ByteExact | PASS | 73 | 92 |  |  |
| intersect_combo_wa_wb | bed | ByteExact | PASS | 122 | 244 | 0.50x |  |
| intersect_combo_u_s | bed | ByteExact | PASS | 47 | 68 |  |  |
| intersect_combo_c_s | bed | ByteExact | PASS | 47 | 69 |  |  |

## bedjaccard

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedjaccard_base | bed | ByteExact | PASS | 38 | 40 |  |  |
| bedjaccard_s | bed | ByteExact | PASS | 33 | 28 |  |  |
| bedjaccard_f50 | bed | ByteExact | PASS | 23 | 32 |  |  |

## bedlinks

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedlinks_base | bed | ByteExact | PASS | 70 | 262 |  |  |

## bedmakewindows

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmakewindows_g_w_winnum | bed | ByteExact | PASS | 14 | 11 |  |  |
| bedmakewindows_g_ws_winnum | bed | ByteExact | PASS | 23 | 22 |  |  |
| bedmakewindows_b_winnum | bed | ByteExact | PASS | 54 | 92 |  |  |
| bedmakewindows_b_srcwinnum | bed | ByteExact | PASS | 77 | 93 |  |  |
| bedmakewindows_b_n_winnum | bed | ByteExact | PASS | 50 | 100 |  |  |
| bedmakewindows_default_none | bed | ByteExact | PASS | 6 | 11 |  |  |
| bedmakewindows_i_src | bed | ByteExact | PASS | 8 | 12 |  |  |

## bedmap

PASS 11 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 73 | 65 |  |  |
| flago_mean | bed | ByteExact | PASS | 71 | 68 |  |  |
| flago_sum | bed | ByteExact | PASS | 57 | 66 |  |  |
| flago_min | bed | ByteExact | PASS | 62 | 66 |  |  |
| flago_max | bed | ByteExact | PASS | 77 | 68 |  |  |
| flago_count | bed | ByteExact | PASS | 51 | 44 |  |  |
| flago_median | bed | ByteExact | PASS | 59 | 73 |  |  |
| flago_stdev | bed | ByteExact | PASS | 60 | 78 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 76 | 69 |  |  |
| combo_c5_sum_s | bed | ByteExact | PASS | 57 | 81 |  |  |
| bedmap_collapse_tiebreak | bed | ByteExact | PASS | 67 | 48 |  |  |

## bedmerge

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 34 | 11 |  |  |
| flagd_50 | bed | ByteExact | PASS | 28 | 12 |  |  |
| flagd_0 | bed | ByteExact | PASS | 26 | 12 |  |  |
| flags | bed | ByteExact | PASS | 41 | 14 |  |  |
| flagS_+ | bed | ByteExact | PASS | 35 | 12 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 30 | 16 |  |  |
| combo_c5_count | bed | ByteExact | PASS | 31 | 12 |  |  |
| combo_d50_c5_sum | bed | ByteExact | PASS | 34 | 14 |  |  |
| combo_s_c5_mean | bed | ByteExact | PASS | 38 | 22 |  |  |
| bedmerge_collapse_tiebreak | bed | ByteExact | PASS | 37 | 13 |  |  |

## bedmulticov

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmulticov_one_bam | bam | ByteExact | PASS | 615 | 35310 |  |  |
| bedmulticov_mapq20 | bam | ByteExact | PASS | 730 | 34550 |  |  |

## bedmultiinter

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmultiinter_two | bed | ByteExact | PASS | 35 | 49 |  |  |
| bedmultiinter_two_names | bed | ByteExact | PASS | 32 | 49 |  |  |

## bednuc

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bednuc_base | fasta | ByteExact | PASS | 309 | 406 |  |  |

## bedoverlap

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedoverlap_base | bed | ByteExact | PASS | 14 | 6 |  |  |

## bedpairtobed

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedpairtobed_base | bed | ByteExact | PASS | 99 | 88 |  |  |

## bedpairtopair

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedpairtopair_base | bed | ByteExact | PASS | 35 | 46 |  |  |

## bedrandom

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedrandom_n50_l100_seed | bed | ByteExact | PASS | 10 | 1 |  |  |
| bedrandom_n100_l500_seed | bed | ByteExact | PASS | 4 | 2 |  |  |

## bedreldist

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedreldist_base | bed | ByteExact | PASS | 38 | 39 |  |  |
| bedreldist_detail | bed | ByteExact | PASS | 122 | 67 |  |  |

## bedsample

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsample_n50_seed | bed | ByteExact | PASS | 14 | 10 |  |  |
| bedsample_n200_seed | bed | ByteExact | PASS | 5 | 10 |  |  |

## bedshift

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 11 | 1 |  |  |
| flags_50 | bed | ByteExact | PASS | 16 | 27 |  |  |
| flags_50 | bed | ByteExact | PASS | 20 | 28 |  |  |
| combo_p_m | bed | ByteExact | PASS | 20 | 27 |  |  |
| combo_s_pct | bed | ByteExact | PASS | 19 | 27 |  |  |

## bedshuffle

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedshuffle_seed | bed | ByteExact | PASS | 29 | 29 |  |  |
| bedshuffle_seed_chrom | bed | ByteExact | PASS | 31 | 30 |  |  |

## bedslop

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 12 | 1 |  |  |
| flagb_50 | bed | ByteExact | PASS | 18 | 27 |  |  |
| flagb_100 | bed | ByteExact | PASS | 22 | 27 |  |  |
| flagpct | bed | ByteExact | PASS | 4 | 1 |  |  |
| flags | bed | ByteExact | PASS | 3 | 1 |  |  |
| combo_l_r | bed | ByteExact | PASS | 17 | 27 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 18 | 27 |  |  |
| combo_b_s | bed | ByteExact | PASS | 19 | 29 |  |  |

## bedsort

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsort_sizeA | bed | ByteExact | PASS | 57 | 55 |  |  |
| bedsort_chrThenSizeA | bed | ByteExact | PASS | 46 | 43 |  |  |
| bedsort_default_tiebreak | bed | ByteExact | PASS | 37 | 34 |  |  |
| bedsort_sizeD_tiebreak | bed | ByteExact | PASS | 48 | 57 |  |  |

## bedspacing

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedspacing_base | bed | ByteExact | PASS | 22 | 21 |  |  |

## bedsplit

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsplit_simple_n3 | bed | ByteExact | PASS | 36 | 35 |  |  |
| bedsplit_size_n3 | bed | ByteExact | PASS | 39 | 60 |  |  |

## bedsubtract

PASS 9 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 332 | 93 |  |  |
| flagA | bed | ByteExact | PASS | 318 | 57 |  |  |
| flagN | bed | ByteExact | PASS | 4 | 2 |  |  |
| flags | bed | ByteExact | PASS | 1332 | 81 |  |  |
| flagS | bed | ByteExact | PASS | 1389 | 87 |  |  |
| flagf_0_5 | bed | ByteExact | PASS | 322 | 84 |  |  |
| combo_A_s | bed | ByteExact | PASS | 1322 | 54 |  |  |
| combo_N_f | bed | ByteExact | PASS | 320 | 77 |  |  |
| bedsubtract_reciprocal_r_unsupported | bed | ByteExact | PASS | 320 | 82 |  |  |

## bedsummary

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsummary_format_and_missing_g | bed | ByteExact | PASS | 25 | 18 |  |  |

## bedtag

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtag_tag | bam | BAMDecoded | PASS | 1898 | 2820 |  |  |

## bedtobam

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtobam_decoded | bed | BAMDecoded | PASS | 49 | 48 |  |  |

## bedunionbedg

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedunionbedg_base | bed | ByteExact | PASS | 35 | 92 |  |  |

## bedwindow

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedwindow_v_w100 | bed | ByteExact | PASS | 199 | 83 |  |  |
| bedwindow_c_w100 | bed | ByteExact | PASS | 198 | 93 |  |  |
| bedwindow_v_lr | bed | ByteExact | PASS | 187 | 85 |  |  |
| bedwindow_join_w100 | bed | ByteExact | PASS | 310 | 198 |  |  |

## bgzip

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bgzip_decompress | vcf | ByteExact | PASS | 26 | 17 |  |  |
| bgzip_decompress_heavy | vcf | ByteExact | PASS | 18 | 16 | 1.12x |  |
| bgzip_compress | vcf_plain | BGZFDecoded | PASS | 29 | 99 |  |  |
| bgzip_reindex | vcf | ByteExact | PASS | 3 | 11 |  |  |

## fastp

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| fastp_cut_tail | fastq | ByteExact | PASS | 2995 | 692 |  |  |
| fastp_default_filter | fastq | ByteExact | PASS | 2849 | 444 |  |  |
| fastp_detect_adapter_pe_heavy | fastq_paired | ByteExact | PASS | 22963 | 2898 | 7.92x |  |

## htsfile

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| htsfile_identify | vcf | ByteExact | PASS | 13 | 4 |  |  |
| htsfile_copy | vcf | ByteExact | PASS | 88 | 48 |  |  |

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
| prinseq_min_len | fastq | ByteExact | PASS | 1485 | 2456 |  |  |
| prinseq_max_len | fastq | ByteExact | PASS | 1219 | 2207 |  |  |
| prinseq_trim_left | fastq | ByteExact | PASS | 1384 | 2473 |  |  |
| prinseq_trim_right | fastq | ByteExact | PASS | 1427 | 2446 |  |  |
| prinseq_trim_qual_right | fastq | ByteExact | PASS | 1576 | 7827 |  |  |
| prinseq_trim_qual_left | fastq | ByteExact | PASS | 1490 | 6787 |  |  |
| prinseq_min_qual_mean | fastq | ByteExact | PASS | 1484 | 6729 |  |  |
| prinseq_max_ns | fastq | ByteExact | PASS | 1401 | 2307 |  |  |

## samtools

PASS 72 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_bam_base | bam | ByteExact | PASS | 960 | 525 |  |  |
| view_bam_flagh | bam | ByteExact | PASS | 865 | 438 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 7 | 3 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 536 | 315 |  |  |
| view_bam_flagq_20 | bam | ByteExact | PASS | 881 | 432 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 787 | 403 |  |  |
| view_bam_flagf_x2 | bam | ByteExact | PASS | 529 | 309 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 805 | 367 |  |  |
| view_bam_flagF_x100 | bam | ByteExact | PASS | 860 | 427 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 697 | 379 |  |  |
| view_bam_flagL_bed | bam | ByteExact | PASS | 878 | 440 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 516 | 308 |  |  |
| view_bam_combo_count_L | bam | ByteExact | PASS | 564 | 319 |  |  |
| view_bam_combo_f2_F256_q20_count | bam | ByteExact | PASS | 529 | 310 |  |  |
| view_bam_combo_header_plus_body | bam | ByteExact | PASS | 830 | 408 |  |  |
| samtools_view_region_contig | bam | ByteExact | PASS | 107 | 53 |  |  |
| samtools_view_region_range | bam | ByteExact | PASS | 169 | 4 |  |  |
| samtools_view_region_count | bam | ByteExact | PASS | 140 | 43 |  |  |
| samtools_view_cram_body | cram | ByteExact | PASS | 2351 | 398 |  |  |
| samtools_view_cram_header | cram | ByteExact | PASS | 12 | 5 |  |  |
| samtools_view_cram_count | cram | ByteExact | PASS | 2031 | 43 |  |  |
| samtools_view_cram_q30 | cram | ByteExact | PASS | 2403 | 329 |  |  |
| samtools_view_cram_decode_sam_heavy | cram | ByteExact | PASS | 2359 | 309 | 7.63x |  |
| view_subsample_seed | bam | ByteExact | PASS | 700 | 367 |  |  |
| samtools_sort_name_sam | bam | ByteExact | PASS | 1327 | 606 |  |  |
| samtools_sort_byname_tag | bam | ByteExact | PASS | 1092 | 564 |  |  |
| samtools_sort_name_sam_heavy | bam | ByteExact | PASS | 1260 | 579 | 2.18x |  |
| samtools_sort_coord_sam | bam | ByteExact | PASS | 787 | 458 |  |  |
| samtools_flagstat | bam | ByteExact | PASS | 402 | 309 |  |  |
| samtools_idxstats | bam | ByteExact | PASS | 7 | 3 |  |  |
| samtools_stats | bam | ByteExact | PASS | 999 | 857 |  |  |
| samtools_quickcheck | bam | ByteExact | PASS | 6 | 4 |  |  |
| samtools_dict | fasta | ByteExact | PASS | 117 | 62 |  |  |
| samtools_stats_heavy | bam | ByteExact | PASS | 939 | 856 | 1.10x |  |
| depth_base | bam | ByteExact | PASS | 952 | 924 |  |  |
| depth_flaga | bam | ByteExact | PASS | 1092 | 955 |  |  |
| depth_flagr_chr1 | bam | ByteExact | PASS | 329 | 113 |  |  |
| depth_flagb_bed | bam | ByteExact | PASS | 1497 | 981 |  |  |
| depth_combo_all_region | bam | ByteExact | PASS | 243 | 6 |  |  |
| depth_combo_all_bed | bam | ByteExact | PASS | 1471 | 956 |  |  |
| samtools_coverage | bam | ByteExact | PASS | 12804 | 717 |  |  |
| samtools_coverage_region | bam | ByteExact | PASS | 7193 | 91 |  |  |
| samtools_depth_mapq_filter | bam | ByteExact | PASS | 1126 | 1042 |  |  |
| samtools_depth_baseq_filter | bam | ByteExact | PASS | 1034 | 865 |  |  |
| samtools_calmd | bam | ByteExact | PASS | 1134 | 522 |  |  |
| samtools_consensus | bam | ByteExact | PASS | 14565 | 2788 |  |  |
| samtools_consensus_region | bam | ByteExact | PASS | 929 | 8 |  |  |
| samtools_fastq | bam | ByteExact | PASS | 1149 | 608 |  |  |
| samtools_fastq_n | bam | ByteExact | PASS | 1069 | 614 |  |  |
| samtools_fastq_heavy | bam | ByteExact | PASS | 1001 | 575 | 1.74x |  |
| samtools_tview_text | bam | ByteExact | PASS | 661 | 80 |  |  |
| samtools_mpileup_pileup | bam | ByteExact | PASS | 21525 | 10331 |  |  |
| cat_concat | bam | BAMDecoded | PASS | 3688 | 221 |  |  |
| samtools_markdup | bam | BAMDecoded | PASS | 2596 | 2485 |  |  |
| samtools_fixmate | bam | BAMDecoded | PASS | 1956 | 2865 |  |  |
| samtools_addreplacerg | bam | BAMDecoded | PASS | 1845 | 495 |  |  |
| samtools_merge | bam | BAMDecoded | PASS | 2725 | 2830 |  |  |
| samtools_reheader | bam | BAMDecoded | PASS | 1813 | 49 |  |  |
| samtools_split | bam | BAMDecoded | PASS | 1975 | 2710 |  |  |
| samtools_import | fastq | BAMDecoded | PASS | 1532 | 233 |  |  |
| samtools_phase | bam | ByteExact | PASS | 10085 | 2115 |  |  |
| view_bam_base | bam | ByteExact | PASS | 922 | 585 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 9 | 4 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 530 | 310 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 831 | 414 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 830 | 385 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 750 | 390 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 543 | 313 |  |  |
| view_bam_combo_header_and_reads | bam | ByteExact | PASS | 913 | 453 |  |  |
| view_cram_body | cram | ByteExact | PASS | 2485 | 343 |  |  |
| view_cram_count | cram | ByteExact | PASS | 2006 | 45 |  |  |
| view_cram_decode_sam_heavy | cram | ByteExact | PASS | 2393 | 361 | 6.62x |  |

## seqtk

PASS 32 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| seqtk_seq_fq_base | fastq | ByteExact | PASS | 367 | 192 |  |  |
| seqtk_seq_fq_flagA | fastq | ByteExact | PASS | 153 | 129 |  |  |
| seqtk_seq_fq_flagr | fastq | ByteExact | PASS | 218 | 210 |  |  |
| seqtk_seq_fq_flagL_95 | fastq | ByteExact | PASS | 224 | 139 |  |  |
| seqtk_seq_fq_flagq_20 | fastq | ByteExact | PASS | 225 | 223 |  |  |
| seqtk_seq_fq_flagl_60 | fastq | ByteExact | PASS | 183 | 274 |  |  |
| seqtk_seq_fq_combo_rev_upper | fastq | ByteExact | PASS | 204 | 226 |  |  |
| seqtk_seq_fq_combo_qmask_n | fastq | ByteExact | PASS | 190 | 156 |  |  |
| seqtk_comp_fa | fasta | ByteExact | PASS | 126 | 80 |  |  |
| seqtk_comp_fq | fastq | ByteExact | PASS | 549 | 559 |  |  |
| seqtk_fqchk | fastq | ByteExact | PASS | 147 | 87 |  |  |
| seqtk_fqchk_q20 | fastq | ByteExact | PASS | 145 | 86 |  |  |
| seqtk_size_fq | fastq | ByteExact | PASS | 77 | 42 |  |  |
| seqtk_size_fa | fasta | ByteExact | PASS | 27 | 8 |  |  |
| seqtk_trimfq | fastq | ByteExact | PASS | 356 | 205 |  |  |
| seqtk_trimfq_q | fastq | ByteExact | PASS | 291 | 184 |  |  |
| seqtk_trimfq_be | fastq | ByteExact | PASS | 164 | 124 |  |  |
| seqtk_sample_count | fastq | ByteExact | PASS | 106 | 43 |  |  |
| seqtk_sample_frac | fastq | ByteExact | PASS | 112 | 55 |  |  |
| seqtk_hpc_fa | fasta | ByteExact | PASS | 71 | 53 |  |  |
| seqtk_hpc_fq | fastq | ByteExact | PASS | 322 | 204 |  |  |
| seqtk_gap_fa | fasta | ByteExact | PASS | 49 | 18 |  |  |
| seqtk_subseq_bed | fasta | ByteExact | PASS | 79 | 204 |  |  |
| seqtk_mergepe | fastq_paired | ByteExact | PASS | 428 | 422 |  |  |
| seqtk_dropse | fastq | ByteExact | PASS | 112 | 46 |  |  |
| seqtk_randbase | fasta | ByteExact | PASS | 135 | 68 |  |  |
| seqtk_telo | fasta | ByteExact | PASS | 31 | 10 |  |  |
| seqtk_listhet | fasta | ByteExact | PASS | 43 | 17 |  |  |
| seqtk_hety | fasta | ByteExact | PASS | 97 | 73 |  |  |
| seqtk_seq_fa | fasta | ByteExact | PASS | 59 | 27 |  |  |
| seqtk_seq_fq_to_fa_heavy | fastq | ByteExact | PASS | 161 | 108 | 1.49x |  |
| seqtk_cutN | fasta | ByteExact | PASS | 63 | 67 |  |  |

## sickle

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| sickle_se_base | fastq | ByteExact | PASS | 1188 | 1140 |  |  |
| sickle_se_q30 | fastq | ByteExact | PASS | 1088 | 1108 |  |  |
| sickle_se_l30 | fastq | ByteExact | PASS | 1117 | 1172 |  |  |
| sickle_se_no5prime | fastq | ByteExact | PASS | 1206 | 1107 |  |  |
| sickle_se_truncn | fastq | ByteExact | PASS | 1095 | 1101 |  |  |
| sickle_se_q30_l40 | fastq | ByteExact | PASS | 995 | 1059 |  |  |
| sickle_pe_base | fastq_paired | ByteExact | PASS | 2354 | 2261 | 1.04x |  |
| sickle_se_cli_default_window | fastq | ByteExact | PASS | 1101 | 1172 |  |  |

## skewer

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| skewer_se_base | fastq | ByteExact | PASS | 2008 | 1895 |  |  |
| skewer_se_minlen30 | fastq | ByteExact | PASS | 1954 | 1907 |  |  |
| skewer_se_full_heavy | fastq | ByteExact | PASS | 2065 | 1958 | 1.05x |  |

## tabix

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| tabix_region_contig | vcf | ByteExact | PASS | 46 | 6 |  |  |
| tabix_region_range | vcf | ByteExact | PASS | 4 | 3 |  |  |
| tabix_region_chr2 | vcf | ByteExact | PASS | 31 | 6 |  |  |
| tabix_region_with_header | vcf | ByteExact | PASS | 41 | 6 |  |  |
| tabix_list_chroms | vcf | ByteExact | PASS | 5 | 3 |  |  |
| tabix_region_heavy | vcf | ByteExact | PASS | 37 | 6 | 6.12x |  |
| tabix_regions_bed | vcf | ByteExact | PASS | 41504 | 287 |  |  |

## vcftools

PASS 23 · SIMILAR 0 · DIVERGE 2 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| vcftools_freq | vcf_plain | ByteExact | PASS | 1768 | 104 |  |  |
| vcftools_counts | vcf_plain | ByteExact | PASS | 2429 | 84 |  |  |
| vcftools_freq2 | vcf_plain | ByteExact | PASS | 2238 | 92 |  |  |
| vcftools_depth | vcf_plain | ByteExact | PASS | 61 | 44 |  |  |
| vcftools_site_depth | vcf_plain | ByteExact | PASS | 1971 | 1468 |  |  |
| vcftools_site_mean_depth | vcf_plain | ByteExact | DIVERGE | 2364 | 1551 |  | output file ".ldepth.mean": first diff at line 2:   ours:     chr1	15	56	-nan   upstream: chr1	15	56	nan |
| vcftools_site_pi | vcf_plain | ByteExact | PASS | 1788 | 1501 |  |  |
| vcftools_window_pi | vcf_plain | ByteExact | PASS | 588 | 329 |  |  |
| vcftools_tstv_summary | vcf_plain | ByteExact | PASS | 58 | 34 |  |  |
| vcftools_missing_indv | vcf_plain | ByteExact | PASS | 55 | 45 |  |  |
| vcftools_missing_site | vcf_plain | ByteExact | PASS | 2359 | 1528 |  |  |
| vcftools_het | vcf_plain | ByteExact | PASS | 74 | 61 |  |  |
| vcftools_singletons | vcf_plain | ByteExact | PASS | 2138 | 1966 |  |  |
| vcftools_recode_heavy | vcf_plain | ByteExact | PASS | 98 | 162 |  |  |
| vcftools_het_multi | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_relatedness | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_relatedness2 | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_freq_multi | vcf_multi_plain | ByteExact | PASS | 3 | 1 |  |  |
| vcftools_missing_indv_multi | vcf_multi_plain | ByteExact | PASS | 3 | 1 |  |  |
| vcftools_window_pi_heavy | vcf_plain | ByteExact | PASS | 1074 | 762 | 1.41x |  |
| vcftools_geno_r2 | vcf_multi_plain | ByteExact | PASS | 4 | 1 | 3.43x |  |
| vcftools_hap_r2 | vcf_multi_plain | ByteExact | PASS | 1115 | 18354 | 0.06x |  |
| vcftools_matrix012 | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_lroh | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_hardy | vcf_plain | ByteExact | DIVERGE | 2565 | 1659 |  | output file ".hwe": first diff at line 2:   ours:     chr1	15	0/0/1	0.00/0.00/1.00	-nan	1.000000e+00	1.000000e+00	1.0000... |

