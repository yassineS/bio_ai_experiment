# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-23T10:57:07Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 400 | 382 | 4 | 3 | 11 | 0 |

## bcftools

PASS 66 · SIMILAR 1 · DIVERGE 0 · SKIP 1 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_base | vcf_plain | ByteExact | PASS | 54 | 43 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 45 | 38 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 29 | 31 |  |  |
| view_flagV_indels | vcf_plain | ByteExact | PASS | 54 | 41 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 52 | 38 |  |  |
| view_flagi_DP>40 | vcf_plain | ByteExact | PASS | 6 | 4 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 47 | 45 |  |  |
| view_flagO_v | vcf_plain | ByteExact | PASS | 56 | 43 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 50 | 39 |  |  |
| view_combo_indels_dp | vcf_plain | ByteExact | PASS | 6 | 4 |  |  |
| view_combo_excl_snps_qual | vcf_plain | ByteExact | PASS | 31 | 29 |  |  |
| bcftools_view_r_region | vcf | ByteExact | PASS | 9 | 4 |  |  |
| bcftools_view_t_targets | vcf | ByteExact | PASS | 35 | 36 |  |  |
| bcftools_view_r_multi | vcf | ByteExact | PASS | 231 | 12 |  |  |
| bcftools_view_s_sample | vcf_multi_plain | ByteExact | PASS | 369 | 110 |  |  |
| bcftools_view_G_drop_gt | vcf_multi_plain | ByteExact | PASS | 289 | 89 |  |  |
| bcftools_view_c_minac | vcf_multi_plain | ByteExact | PASS | 232 | 129 |  |  |
| bcftools_view_C_maxac | vcf_multi_plain | ByteExact | PASS | 187 | 77 |  |  |
| bcftools_view_q_minaf | vcf_multi_plain | ByteExact | PASS | 230 | 109 |  |  |
| bcftools_view_Q_maxaf | vcf_multi_plain | ByteExact | PASS | 254 | 105 |  |  |
| bcftools_query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 32 | 21 |  |  |
| bcftools_query_info_dp | vcf_plain | ByteExact | PASS | 31 | 28 |  |  |
| bcftools_query_qual_filter | vcf_plain | ByteExact | PASS | 35 | 26 |  |  |
| bcftools_query_info_dp_filtered | vcf_plain | ByteExact | PASS | 7 | 3 |  |  |
| bcftools_query_info_dp_excluded | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_gt_multi | vcf_multi_plain | ByteExact | PASS | 145 | 112 |  |  |
| bcftools_query_list_samples | vcf_multi_plain | ByteExact | PASS | 6 | 4 |  |  |
| bcftools_query_heavy_full | vcf_plain | ByteExact | PASS | 37 | 30 | 1.22x |  |
| bcftools_norm_split | vcf_plain | ByteExact | PASS | 97 | 65 |  |  |
| bcftools_norm_dedup | vcf_plain | ByteExact | PASS | 83 | 55 |  |  |
| bcftools_norm_check_ref | vcf_plain | ByteExact | PASS | 96 | 68 |  |  |
| filter_base | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_flagi_QUAL>=30 | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_combo_soft_lowqual | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_combo_snpgap | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_sort | vcf_plain | ByteExact | PASS | 80 | 85 |  |  |
| bcftools_stats | vcf_plain | ByteExact | PASS | 28 | 25 |  |  |
| bcftools_stats_multi | vcf_multi_plain | ByteExact | PASS | 104 | 41 |  |  |
| bcftools_head | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_annotate_drop_af | vcf_plain | ByteExact | PASS | 77 | 41 |  |  |
| bcftools_annotate_drop_dp | vcf_plain | ByteExact | PASS | 69 | 40 |  |  |
| bcftools_gtcheck | vcf_multi_plain | ByteExact | PASS | 1467 | 125 |  |  |
| bcftools_mpileup | bam | ByteExact | PASS | 91607 | 28131 |  |  |
| bcftools_mpileup_heavy | bam | ByteExact | PASS | 88523 | 28220 | 3.14x |  |
| bcftools_plugin_fill_tags | vcf_multi_plain | ByteExact | PASS | 1495 | 274 |  |  |
| bcftools_plugin_fill_tags_an_ac | vcf_multi_plain | ByteExact | PASS | 627 | 167 |  |  |
| bcftools_view_c_update_acan | vcf_multi_plain | ByteExact | PASS | 376 | 162 |  |  |
| bcftools_annotate_drop_filter | vcf_plain | ByteExact | PASS | 105 | 55 |  |  |
| bcftools_roh | vcf_multi_plain | ByteExact | PASS | 1409 | 426 |  |  |
| bcftools_consensus | vcf | ByteExact | PASS | 7892 | 128 |  |  |
| bcftools_concat | vcf | ByteExact | PASS | 202 | 180 |  |  |
| bcftools_norm_join | vcf_plain | ByteExact | PASS | 110 | 88 |  |  |
| bcftools_call | vcf_plain | Similarity | SIMILAR | 46761 | 33279 |  |  |
| bcftools_csq | vcf_plain | ByteExact | SKIP | 0 | 0 |  | bcftools csq --force over the multi-gene fixture is byte-identical to upstream EXCEPT for 3 records, and that residual i... |
| bcftools_isec | vcf | ByteExact | PASS | 5479 | 325 |  |  |
| bcftools_merge | vcf_plain | ByteExact | PASS | 882 | 154 |  |  |
| bcftools_convert | vcf_plain | ByteExact | PASS | 84 | 46 |  |  |
| bcftools_reheader | vcf_multi_plain | BGZFDecoded | PASS | 508 | 19 |  |  |
| view_base | vcf_plain | ByteExact | PASS | 59 | 42 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 45 | 37 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 29 | 30 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 48 | 36 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 46 | 36 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 46 | 35 |  |  |
| query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 30 | 20 |  |  |
| query_gt | vcf_plain | ByteExact | PASS | 30 | 34 |  |  |
| query_info_dp | vcf_plain | ByteExact | PASS | 30 | 27 |  |  |

## bed12tobed6

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bed12tobed6_score_dropped | bed12 | ByteExact | PASS | 46 | 77 |  |  |

## bedannotate

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedannotate_default_header_order | bed | ByteExact | PASS | 141 | 440 |  |  |

## bedbamtobed

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bam | ByteExact | PASS | 667 | 600 |  |  |
| flagsplit | bam | ByteExact | PASS | 734 | 782 |  |  |
| flaged | bam | ByteExact | PASS | 10 | 3 |  |  |
| flagcigar | bam | ByteExact | PASS | 682 | 691 |  |  |
| combo_split_cigar | bam | ByteExact | PASS | 678 | 687 |  |  |

## bedclosest

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 3379 | 134 |  |  |
| flagd | bed | ByteExact | PASS | 2937 | 256 |  |  |
| flagio | bed | ByteExact | PASS | 2862 | 55 |  |  |
| flagiu | bed | ByteExact | PASS | 10 | 3 |  |  |
| flagt_first | bed | ByteExact | PASS | 2696 | 57 |  |  |
| flagt_last | bed | ByteExact | PASS | 2820 | 49 |  |  |
| flagt_all | bed | ByteExact | PASS | 3177 | 238 |  |  |
| flags | bed | ByteExact | PASS | 3059 | 82 |  |  |
| flagN | bed | ByteExact | PASS | 3233 | 93 |  |  |
| combo_d_t_first | bed | ByteExact | PASS | 2808 | 57 |  |  |

## bedcluster

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 32 | 29 |  |  |
| flagd_50 | bed | ByteExact | PASS | 30 | 30 |  |  |
| flagd_0 | bed | ByteExact | PASS | 31 | 32 |  |  |
| bedcluster_s | bed | ByteExact | PASS | 43 | 43 |  |  |
| bedcluster_d50_s | bed | ByteExact | PASS | 43 | 48 |  |  |

## bedcomplement

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedcomplement_base | bed | ByteExact | PASS | 21 | 12 |  |  |

## bedcoverage

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 144 | 174 |  |  |
| flagcounts | bed | ByteExact | PASS | 106 | 142 |  |  |
| flagd | bed | ByteExact | PASS | 7490 | 22777 |  |  |
| flaghist | bed | ByteExact | PASS | 705 | 519 |  |  |
| flags | bed | ByteExact | PASS | 130 | 146 |  |  |
| flagS | bed | ByteExact | PASS | 96 | 134 |  |  |
| flagmean | bed | ByteExact | PASS | 125 | 188 |  |  |
| combo_counts_s | bed | ByteExact | PASS | 83 | 150 |  |  |

## bedexpand

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedexpand_c5 | bed | ByteExact | PASS | 26 | 36 |  |  |
| bedexpand_trailing_comma_col11 | bed12 | ByteExact | PASS | 46 | 69 |  |  |

## bedfisher

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedfisher_overlap_count | bed | ByteExact | PASS | 147 | 39 |  |  |

## bedflank

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 30 | 7 |  |  |
| flagb_50 | bed | ByteExact | PASS | 56 | 43 |  |  |
| flags | bed | ByteExact | PASS | 11 | 3 |  |  |
| combo_l_r | bed | ByteExact | PASS | 34 | 50 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 34 | 52 |  |  |
| combo_b_s | bed | ByteExact | PASS | 33 | 43 |  |  |

## bedgenomecov

PASS 5 · SIMILAR 3 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | Similarity | SIMILAR | 397 | 264 |  |  |
| flagbg | bed | ByteExact | PASS | 120 | 84 |  |  |
| flagbga | bed | ByteExact | PASS | 110 | 83 |  |  |
| flagd | bed | ByteExact | PASS | 1330 | 8331 |  |  |
| flagdz | bed | ByteExact | PASS | 763 | 7706 |  |  |
| flagmax_5 | bed | Similarity | SIMILAR | 302 | 263 |  |  |
| flagstrand_+ | bed | Similarity | SIMILAR | 283 | 239 |  |  |
| combo_bg_strand | bed | ByteExact | PASS | 143 | 64 |  |  |

## bedgetfasta

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | fasta | ByteExact | PASS | 241 | 459 |  |  |
| flags | fasta | ByteExact | PASS | 386 | 602 |  |  |
| flagname | fasta | ByteExact | PASS | 178 | 431 |  |  |
| flagtab | fasta | ByteExact | PASS | 148 | 395 |  |  |
| combo_s_name | fasta | ByteExact | PASS | 332 | 605 |  |  |
| combo_s_tab | fasta | ByteExact | PASS | 312 | 576 |  |  |

## bedgroupby

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 10 | 1 |  |  |
| flago_mean | bed | ByteExact | PASS | 13 | 11 |  |  |
| flago_sum | bed | ByteExact | PASS | 14 | 11 |  |  |
| flago_min | bed | ByteExact | PASS | 13 | 11 |  |  |
| flago_max | bed | ByteExact | PASS | 14 | 11 |  |  |
| flago_count | bed | ByteExact | PASS | 11 | 11 |  |  |
| combo_g1_c5_mean | bed | ByteExact | PASS | 13 | 11 |  |  |
| combo_g1_c5_count | bed | ByteExact | PASS | 14 | 11 |  |  |

## bedigv

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedigv_base | bed | ByteExact | PASS | 42 | 108 |  |  |

## bedintersect

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| intersect_base | bed | ByteExact | PASS | 116 | 114 |  |  |
| intersect_flagc | bed | ByteExact | PASS | 48 | 69 |  |  |
| intersect_flagv | bed | ByteExact | PASS | 43 | 64 |  |  |
| intersect_flagu | bed | ByteExact | PASS | 49 | 69 |  |  |
| intersect_flagwa | bed | ByteExact | PASS | 74 | 115 |  |  |
| intersect_flagwb | bed | ByteExact | PASS | 136 | 215 |  |  |
| intersect_flags | bed | ByteExact | PASS | 60 | 89 |  |  |
| intersect_combo_wa_wb | bed | ByteExact | PASS | 110 | 242 | 0.46x |  |
| intersect_combo_u_s | bed | ByteExact | PASS | 46 | 74 |  |  |
| intersect_combo_c_s | bed | ByteExact | PASS | 50 | 69 |  |  |

## bedjaccard

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedjaccard_base | bed | ByteExact | PASS | 31 | 24 |  |  |
| bedjaccard_s | bed | ByteExact | PASS | 33 | 27 |  |  |
| bedjaccard_f50 | bed | ByteExact | PASS | 28 | 24 |  |  |

## bedlinks

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedlinks_base | bed | ByteExact | PASS | 71 | 241 |  |  |

## bedmakewindows

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmakewindows_g_w_winnum | bed | ByteExact | PASS | 18 | 15 |  |  |
| bedmakewindows_g_ws_winnum | bed | ByteExact | PASS | 16 | 24 |  |  |
| bedmakewindows_b_winnum | bed | ByteExact | PASS | 56 | 94 |  |  |
| bedmakewindows_b_srcwinnum | bed | ByteExact | PASS | 56 | 93 |  |  |
| bedmakewindows_b_n_winnum | bed | ByteExact | PASS | 52 | 96 |  |  |
| bedmakewindows_default_none | bed | ByteExact | PASS | 7 | 11 |  |  |
| bedmakewindows_i_src | bed | ByteExact | PASS | 7 | 11 |  |  |

## bedmap

PASS 11 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 71 | 70 |  |  |
| flago_mean | bed | ByteExact | PASS | 63 | 72 |  |  |
| flago_sum | bed | ByteExact | PASS | 60 | 70 |  |  |
| flago_min | bed | ByteExact | PASS | 68 | 74 |  |  |
| flago_max | bed | ByteExact | PASS | 89 | 99 |  |  |
| flago_count | bed | ByteExact | PASS | 59 | 47 |  |  |
| flago_median | bed | ByteExact | PASS | 89 | 93 |  |  |
| flago_stdev | bed | ByteExact | PASS | 352 | 85 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 150 | 83 |  |  |
| combo_c5_sum_s | bed | ByteExact | PASS | 234 | 70 |  |  |
| bedmap_collapse_tiebreak | bed | ByteExact | PASS | 63 | 59 |  |  |

## bedmerge

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 33 | 11 |  |  |
| flagd_50 | bed | ByteExact | PASS | 28 | 11 |  |  |
| flagd_0 | bed | ByteExact | PASS | 21 | 11 |  |  |
| flags | bed | ByteExact | PASS | 36 | 14 |  |  |
| flagS_+ | bed | ByteExact | PASS | 27 | 11 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 32 | 14 |  |  |
| combo_c5_count | bed | ByteExact | PASS | 24 | 12 |  |  |
| combo_d50_c5_sum | bed | ByteExact | PASS | 26 | 14 |  |  |
| combo_s_c5_mean | bed | ByteExact | PASS | 33 | 28 |  |  |
| bedmerge_collapse_tiebreak | bed | ByteExact | PASS | 32 | 13 |  |  |

## bedmulticov

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmulticov_one_bam | bam | ByteExact | PASS | 607 | 35280 |  |  |
| bedmulticov_mapq20 | bam | ByteExact | PASS | 627 | 35156 |  |  |

## bedmultiinter

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmultiinter_two | bed | ByteExact | PASS | 42 | 50 |  |  |
| bedmultiinter_two_names | bed | ByteExact | PASS | 39 | 50 |  |  |

## bednuc

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bednuc_base | fasta | ByteExact | PASS | 303 | 404 |  |  |

## bedoverlap

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedoverlap_base | bed | ByteExact | PASS | 20 | 7 |  |  |

## bedpairtobed

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedpairtobed_base | bed | ByteExact | PASS | 99 | 84 |  |  |

## bedpairtopair

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedpairtopair_base | bed | ByteExact | PASS | 39 | 47 |  |  |

## bedrandom

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedrandom_n50_l100_seed | bed | ByteExact | PASS | 11 | 2 |  |  |
| bedrandom_n100_l500_seed | bed | ByteExact | PASS | 4 | 2 |  |  |

## bedreldist

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedreldist_base | bed | ByteExact | PASS | 36 | 42 |  |  |
| bedreldist_detail | bed | ByteExact | PASS | 115 | 51 |  |  |

## bedsample

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsample_n50_seed | bed | ByteExact | PASS | 14 | 10 |  |  |
| bedsample_n200_seed | bed | ByteExact | PASS | 4 | 11 |  |  |

## bedshift

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 18 | 1 |  |  |
| flags_50 | bed | ByteExact | PASS | 17 | 27 |  |  |
| flags_50 | bed | ByteExact | PASS | 20 | 37 |  |  |
| combo_p_m | bed | ByteExact | PASS | 31 | 37 |  |  |
| combo_s_pct | bed | ByteExact | PASS | 26 | 29 |  |  |

## bedshuffle

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedshuffle_seed | bed | ByteExact | PASS | 28 | 29 |  |  |
| bedshuffle_seed_chrom | bed | ByteExact | PASS | 23 | 28 |  |  |

## bedslop

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 10 | 1 |  |  |
| flagb_50 | bed | ByteExact | PASS | 21 | 26 |  |  |
| flagb_100 | bed | ByteExact | PASS | 21 | 28 |  |  |
| flagpct | bed | ByteExact | PASS | 4 | 2 |  |  |
| flags | bed | ByteExact | PASS | 21 | 10 |  |  |
| combo_l_r | bed | ByteExact | PASS | 29 | 29 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 17 | 33 |  |  |
| combo_b_s | bed | ByteExact | PASS | 18 | 31 |  |  |

## bedsort

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsort_sizeA | bed | ByteExact | PASS | 66 | 59 |  |  |
| bedsort_chrThenSizeA | bed | ByteExact | PASS | 45 | 41 |  |  |
| bedsort_default_tiebreak | bed | ByteExact | PASS | 41 | 33 |  |  |
| bedsort_sizeD_tiebreak | bed | ByteExact | PASS | 49 | 66 |  |  |

## bedspacing

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedspacing_base | bed | ByteExact | PASS | 26 | 22 |  |  |

## bedsplit

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsplit_simple_n3 | bed | ByteExact | PASS | 37 | 34 |  |  |
| bedsplit_size_n3 | bed | ByteExact | PASS | 30 | 54 |  |  |

## bedsubtract

PASS 9 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 373 | 106 |  |  |
| flagA | bed | ByteExact | PASS | 341 | 76 |  |  |
| flagN | bed | ByteExact | PASS | 7 | 2 |  |  |
| flags | bed | ByteExact | PASS | 1670 | 151 |  |  |
| flagS | bed | ByteExact | PASS | 1679 | 174 |  |  |
| flagf_0_5 | bed | ByteExact | PASS | 331 | 100 |  |  |
| combo_A_s | bed | ByteExact | PASS | 1337 | 59 |  |  |
| combo_N_f | bed | ByteExact | PASS | 317 | 81 |  |  |
| bedsubtract_reciprocal_r_unsupported | bed | ByteExact | PASS | 315 | 90 |  |  |

## bedsummary

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsummary_format_and_missing_g | bed | ByteExact | PASS | 24 | 16 |  |  |

## bedtag

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtag_tag | bam | BAMDecoded | PASS | 2000 | 2890 |  |  |

## bedtobam

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtobam_decoded | bed | BAMDecoded | PASS | 55 | 52 |  |  |

## bedunionbedg

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedunionbedg_base | bed | ByteExact | PASS | 36 | 95 |  |  |

## bedwindow

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedwindow_v_w100 | bed | ByteExact | PASS | 207 | 85 |  |  |
| bedwindow_c_w100 | bed | ByteExact | PASS | 203 | 103 |  |  |
| bedwindow_v_lr | bed | ByteExact | PASS | 205 | 92 |  |  |
| bedwindow_join_w100 | bed | ByteExact | PASS | 290 | 222 |  |  |

## bgzip

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bgzip_decompress | vcf | ByteExact | PASS | 20 | 18 |  |  |
| bgzip_decompress_heavy | vcf | ByteExact | PASS | 16 | 15 | 1.05x |  |
| bgzip_compress | vcf_plain | BGZFDecoded | PASS | 27 | 101 |  |  |
| bgzip_reindex | vcf | ByteExact | PASS | 4 | 11 |  |  |

## fastp

PASS 2 · SIMILAR 0 · DIVERGE 1 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| fastp_cut_tail | fastq | ByteExact | PASS | 2867 | 1285 |  |  |
| fastp_default_filter | fastq | ByteExact | PASS | 2837 | 327 |  |  |
| fastp_detect_adapter_pe_heavy | fastq_paired | ByteExact | DIVERGE | 6 | 2779 | 0.00x | exit mismatch: ours_err=exit status 1 upstream_err=<nil> ours stderr: Error: must specify either:   Single-end: -i/--in1... |

## htsfile

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| htsfile_identify | vcf | ByteExact | PASS | 14 | 5 |  |  |
| htsfile_copy | vcf | ByteExact | PASS | 81 | 44 |  |  |

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
| prinseq_min_len | fastq | ByteExact | PASS | 1394 | 2415 |  |  |
| prinseq_max_len | fastq | ByteExact | PASS | 1281 | 2290 |  |  |
| prinseq_trim_left | fastq | ByteExact | PASS | 1466 | 2487 |  |  |
| prinseq_trim_right | fastq | ByteExact | PASS | 1395 | 2505 |  |  |
| prinseq_trim_qual_right | fastq | ByteExact | PASS | 1555 | 7584 |  |  |
| prinseq_trim_qual_left | fastq | ByteExact | PASS | 1536 | 6653 |  |  |
| prinseq_min_qual_mean | fastq | ByteExact | PASS | 1439 | 6744 |  |  |
| prinseq_max_ns | fastq | ByteExact | PASS | 1276 | 2271 |  |  |

## samtools

PASS 72 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_bam_base | bam | ByteExact | PASS | 1084 | 752 |  |  |
| view_bam_flagh | bam | ByteExact | PASS | 1025 | 456 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 7 | 4 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 524 | 308 |  |  |
| view_bam_flagq_20 | bam | ByteExact | PASS | 799 | 421 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 754 | 426 |  |  |
| view_bam_flagf_x2 | bam | ByteExact | PASS | 528 | 307 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 752 | 371 |  |  |
| view_bam_flagF_x100 | bam | ByteExact | PASS | 890 | 490 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 698 | 363 |  |  |
| view_bam_flagL_bed | bam | ByteExact | PASS | 872 | 592 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 534 | 381 |  |  |
| view_bam_combo_count_L | bam | ByteExact | PASS | 578 | 318 |  |  |
| view_bam_combo_f2_F256_q20_count | bam | ByteExact | PASS | 524 | 310 |  |  |
| view_bam_combo_header_plus_body | bam | ByteExact | PASS | 881 | 480 |  |  |
| samtools_view_region_contig | bam | ByteExact | PASS | 113 | 53 |  |  |
| samtools_view_region_range | bam | ByteExact | PASS | 161 | 4 |  |  |
| samtools_view_region_count | bam | ByteExact | PASS | 140 | 42 |  |  |
| samtools_view_cram_body | cram | ByteExact | PASS | 2428 | 386 |  |  |
| samtools_view_cram_header | cram | ByteExact | PASS | 7 | 4 |  |  |
| samtools_view_cram_count | cram | ByteExact | PASS | 2002 | 43 |  |  |
| samtools_view_cram_q30 | cram | ByteExact | PASS | 2289 | 337 |  |  |
| samtools_view_cram_decode_sam_heavy | cram | ByteExact | PASS | 2291 | 312 | 7.34x |  |
| view_subsample_seed | bam | ByteExact | PASS | 743 | 376 |  |  |
| samtools_sort_name_sam | bam | ByteExact | PASS | 1357 | 653 |  |  |
| samtools_sort_byname_tag | bam | ByteExact | PASS | 1162 | 584 |  |  |
| samtools_sort_name_sam_heavy | bam | ByteExact | PASS | 1326 | 587 | 2.26x |  |
| samtools_sort_coord_sam | bam | ByteExact | PASS | 757 | 447 |  |  |
| samtools_flagstat | bam | ByteExact | PASS | 402 | 307 |  |  |
| samtools_idxstats | bam | ByteExact | PASS | 7 | 3 |  |  |
| samtools_stats | bam | ByteExact | PASS | 954 | 857 |  |  |
| samtools_quickcheck | bam | ByteExact | PASS | 6 | 4 |  |  |
| samtools_dict | fasta | ByteExact | PASS | 117 | 61 |  |  |
| samtools_stats_heavy | bam | ByteExact | PASS | 952 | 854 | 1.11x |  |
| depth_base | bam | ByteExact | PASS | 883 | 953 |  |  |
| depth_flaga | bam | ByteExact | PASS | 875 | 858 |  |  |
| depth_flagr_chr1 | bam | ByteExact | PASS | 301 | 141 |  |  |
| depth_flagb_bed | bam | ByteExact | PASS | 1728 | 976 |  |  |
| depth_combo_all_region | bam | ByteExact | PASS | 240 | 6 |  |  |
| depth_combo_all_bed | bam | ByteExact | PASS | 1489 | 1097 |  |  |
| samtools_coverage | bam | ByteExact | PASS | 15616 | 945 |  |  |
| samtools_coverage_region | bam | ByteExact | PASS | 10497 | 97 |  |  |
| samtools_depth_mapq_filter | bam | ByteExact | PASS | 1201 | 1005 |  |  |
| samtools_depth_baseq_filter | bam | ByteExact | PASS | 1042 | 875 |  |  |
| samtools_calmd | bam | ByteExact | PASS | 1223 | 524 |  |  |
| samtools_consensus | bam | ByteExact | PASS | 14778 | 2811 |  |  |
| samtools_consensus_region | bam | ByteExact | PASS | 955 | 9 |  |  |
| samtools_fastq | bam | ByteExact | PASS | 1087 | 663 |  |  |
| samtools_fastq_n | bam | ByteExact | PASS | 1109 | 628 |  |  |
| samtools_fastq_heavy | bam | ByteExact | PASS | 1037 | 604 | 1.72x |  |
| samtools_tview_text | bam | ByteExact | PASS | 648 | 83 |  |  |
| samtools_mpileup_pileup | bam | ByteExact | PASS | 21670 | 10139 |  |  |
| cat_concat | bam | BAMDecoded | PASS | 3719 | 154 |  |  |
| samtools_markdup | bam | BAMDecoded | PASS | 2550 | 2498 |  |  |
| samtools_fixmate | bam | BAMDecoded | PASS | 1904 | 2873 |  |  |
| samtools_addreplacerg | bam | BAMDecoded | PASS | 1784 | 421 |  |  |
| samtools_merge | bam | BAMDecoded | PASS | 2777 | 2832 |  |  |
| samtools_reheader | bam | BAMDecoded | PASS | 1777 | 49 |  |  |
| samtools_split | bam | BAMDecoded | PASS | 1919 | 2695 |  |  |
| samtools_import | fastq | BAMDecoded | PASS | 1501 | 270 |  |  |
| samtools_phase | bam | ByteExact | PASS | 10143 | 2120 |  |  |
| view_bam_base | bam | ByteExact | PASS | 1022 | 653 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 9 | 4 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 526 | 309 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 885 | 411 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 738 | 416 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 793 | 408 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 531 | 308 |  |  |
| view_bam_combo_header_and_reads | bam | ByteExact | PASS | 898 | 592 |  |  |
| view_cram_body | cram | ByteExact | PASS | 2495 | 343 |  |  |
| view_cram_count | cram | ByteExact | PASS | 1936 | 43 |  |  |
| view_cram_decode_sam_heavy | cram | ByteExact | PASS | 2367 | 380 | 6.22x |  |

## seqtk

PASS 32 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| seqtk_seq_fq_base | fastq | ByteExact | PASS | 386 | 282 |  |  |
| seqtk_seq_fq_flagA | fastq | ByteExact | PASS | 210 | 140 |  |  |
| seqtk_seq_fq_flagr | fastq | ByteExact | PASS | 285 | 267 |  |  |
| seqtk_seq_fq_flagL_95 | fastq | ByteExact | PASS | 322 | 182 |  |  |
| seqtk_seq_fq_flagq_20 | fastq | ByteExact | PASS | 232 | 202 |  |  |
| seqtk_seq_fq_flagl_60 | fastq | ByteExact | PASS | 191 | 144 |  |  |
| seqtk_seq_fq_combo_rev_upper | fastq | ByteExact | PASS | 225 | 241 |  |  |
| seqtk_seq_fq_combo_qmask_n | fastq | ByteExact | PASS | 334 | 218 |  |  |
| seqtk_comp_fa | fasta | ByteExact | PASS | 127 | 80 |  |  |
| seqtk_comp_fq | fastq | ByteExact | PASS | 533 | 568 |  |  |
| seqtk_fqchk | fastq | ByteExact | PASS | 147 | 85 |  |  |
| seqtk_fqchk_q20 | fastq | ByteExact | PASS | 141 | 84 |  |  |
| seqtk_size_fq | fastq | ByteExact | PASS | 77 | 42 |  |  |
| seqtk_size_fa | fasta | ByteExact | PASS | 41 | 9 |  |  |
| seqtk_trimfq | fastq | ByteExact | PASS | 336 | 218 |  |  |
| seqtk_trimfq_q | fastq | ByteExact | PASS | 456 | 211 |  |  |
| seqtk_trimfq_be | fastq | ByteExact | PASS | 163 | 126 |  |  |
| seqtk_sample_count | fastq | ByteExact | PASS | 104 | 43 |  |  |
| seqtk_sample_frac | fastq | ByteExact | PASS | 112 | 55 |  |  |
| seqtk_hpc_fa | fasta | ByteExact | PASS | 78 | 52 |  |  |
| seqtk_hpc_fq | fastq | ByteExact | PASS | 295 | 209 |  |  |
| seqtk_gap_fa | fasta | ByteExact | PASS | 48 | 18 |  |  |
| seqtk_subseq_bed | fasta | ByteExact | PASS | 93 | 160 |  |  |
| seqtk_mergepe | fastq_paired | ByteExact | PASS | 713 | 756 |  |  |
| seqtk_dropse | fastq | ByteExact | PASS | 104 | 45 |  |  |
| seqtk_randbase | fasta | ByteExact | PASS | 143 | 60 |  |  |
| seqtk_telo | fasta | ByteExact | PASS | 31 | 9 |  |  |
| seqtk_listhet | fasta | ByteExact | PASS | 41 | 14 |  |  |
| seqtk_hety | fasta | ByteExact | PASS | 98 | 73 |  |  |
| seqtk_seq_fa | fasta | ByteExact | PASS | 51 | 27 |  |  |
| seqtk_seq_fq_to_fa_heavy | fastq | ByteExact | PASS | 123 | 82 | 1.49x |  |
| seqtk_cutN | fasta | ByteExact | PASS | 52 | 67 |  |  |

## sickle

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| sickle_se_base | fastq | ByteExact | PASS | 1116 | 1427 |  |  |
| sickle_se_q30 | fastq | ByteExact | PASS | 1131 | 1116 |  |  |
| sickle_se_l30 | fastq | ByteExact | PASS | 1189 | 1115 |  |  |
| sickle_se_no5prime | fastq | ByteExact | PASS | 1132 | 1096 |  |  |
| sickle_se_truncn | fastq | ByteExact | PASS | 1128 | 1111 |  |  |
| sickle_se_q30_l40 | fastq | ByteExact | PASS | 1126 | 1113 |  |  |
| sickle_pe_base | fastq_paired | ByteExact | PASS | 2238 | 2337 | 0.96x |  |
| sickle_se_cli_default_window | fastq | ByteExact | PASS | 1065 | 1139 |  |  |

## skewer

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| skewer_se_base | fastq | ByteExact | PASS | 2050 | 1954 |  |  |
| skewer_se_minlen30 | fastq | ByteExact | PASS | 1961 | 1952 |  |  |
| skewer_se_full_heavy | fastq | ByteExact | PASS | 1962 | 1907 | 1.03x |  |

## tabix

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| tabix_region_contig | vcf | ByteExact | PASS | 50 | 7 |  |  |
| tabix_region_range | vcf | ByteExact | PASS | 5 | 3 |  |  |
| tabix_region_chr2 | vcf | ByteExact | PASS | 33 | 7 |  |  |
| tabix_region_with_header | vcf | ByteExact | PASS | 44 | 6 |  |  |
| tabix_list_chroms | vcf | ByteExact | PASS | 4 | 3 |  |  |
| tabix_region_heavy | vcf | ByteExact | PASS | 37 | 6 | 5.64x |  |
| tabix_regions_bed | vcf | ByteExact | PASS | 41769 | 290 |  |  |

## vcftools

PASS 23 · SIMILAR 0 · DIVERGE 2 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| vcftools_freq | vcf_plain | ByteExact | PASS | 2137 | 99 |  |  |
| vcftools_counts | vcf_plain | ByteExact | PASS | 1915 | 80 |  |  |
| vcftools_freq2 | vcf_plain | ByteExact | PASS | 2273 | 90 |  |  |
| vcftools_depth | vcf_plain | ByteExact | PASS | 55 | 46 |  |  |
| vcftools_site_depth | vcf_plain | ByteExact | PASS | 2055 | 1504 |  |  |
| vcftools_site_mean_depth | vcf_plain | ByteExact | DIVERGE | 2424 | 1540 |  | output file ".ldepth.mean": first diff at line 2:   ours:     chr1	15	56	-nan   upstream: chr1	15	56	nan |
| vcftools_site_pi | vcf_plain | ByteExact | PASS | 2024 | 1530 |  |  |
| vcftools_window_pi | vcf_plain | ByteExact | PASS | 409 | 308 |  |  |
| vcftools_tstv_summary | vcf_plain | ByteExact | PASS | 58 | 44 |  |  |
| vcftools_missing_indv | vcf_plain | ByteExact | PASS | 63 | 44 |  |  |
| vcftools_missing_site | vcf_plain | ByteExact | PASS | 2114 | 1532 |  |  |
| vcftools_het | vcf_plain | ByteExact | PASS | 63 | 56 |  |  |
| vcftools_singletons | vcf_plain | ByteExact | PASS | 2527 | 1988 |  |  |
| vcftools_recode_heavy | vcf_plain | ByteExact | PASS | 101 | 169 |  |  |
| vcftools_het_multi | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_relatedness | vcf_multi_plain | ByteExact | PASS | 12 | 3 |  |  |
| vcftools_relatedness2 | vcf_multi_plain | ByteExact | PASS | 9 | 1 |  |  |
| vcftools_freq_multi | vcf_multi_plain | ByteExact | PASS | 5 | 1 |  |  |
| vcftools_missing_indv_multi | vcf_multi_plain | ByteExact | PASS | 7 | 1 |  |  |
| vcftools_window_pi_heavy | vcf_plain | ByteExact | PASS | 1428 | 774 | 1.84x |  |
| vcftools_geno_r2 | vcf_multi_plain | ByteExact | PASS | 4 | 1 | 3.55x |  |
| vcftools_hap_r2 | vcf_multi_plain | ByteExact | PASS | 1130 | 18459 | 0.06x |  |
| vcftools_matrix012 | vcf_multi_plain | ByteExact | PASS | 5 | 1 |  |  |
| vcftools_lroh | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_hardy | vcf_plain | ByteExact | DIVERGE | 2231 | 1654 |  | output file ".hwe": first diff at line 2:   ours:     chr1	15	0/0/1	0.00/0.00/1.00	-nan	1.000000e+00	1.000000e+00	1.0000... |

