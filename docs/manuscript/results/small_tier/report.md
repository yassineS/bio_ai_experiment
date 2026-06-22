# Parity pipeline report

- Scale: `small`
- Seed: `1`
- Generated: 2026-06-22T13:25:03Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 400 | 382 | 4 | 3 | 11 | 0 |

## bcftools

PASS 66 · SIMILAR 1 · DIVERGE 0 · SKIP 1 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_base | vcf_plain | ByteExact | PASS | 13 | 9 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 9 | 7 |  |  |
| view_flagV_indels | vcf_plain | ByteExact | PASS | 11 | 9 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_flagi_DP>40 | vcf_plain | ByteExact | PASS | 5 | 4 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_flagO_v | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_combo_indels_dp | vcf_plain | ByteExact | PASS | 5 | 4 |  |  |
| view_combo_excl_snps_qual | vcf_plain | ByteExact | PASS | 8 | 7 |  |  |
| bcftools_view_r_region | vcf | ByteExact | PASS | 9 | 4 |  |  |
| bcftools_view_t_targets | vcf | ByteExact | PASS | 9 | 8 |  |  |
| bcftools_view_r_multi | vcf | ByteExact | PASS | 61 | 6 |  |  |
| bcftools_view_s_sample | vcf_multi_plain | ByteExact | PASS | 35 | 16 |  |  |
| bcftools_view_G_drop_gt | vcf_multi_plain | ByteExact | PASS | 36 | 14 |  |  |
| bcftools_view_c_minac | vcf_multi_plain | ByteExact | PASS | 29 | 16 |  |  |
| bcftools_view_C_maxac | vcf_multi_plain | ByteExact | PASS | 23 | 13 |  |  |
| bcftools_view_q_minaf | vcf_multi_plain | ByteExact | PASS | 28 | 16 |  |  |
| bcftools_view_Q_maxaf | vcf_multi_plain | ByteExact | PASS | 29 | 15 |  |  |
| bcftools_query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 8 | 6 |  |  |
| bcftools_query_info_dp | vcf_plain | ByteExact | PASS | 8 | 7 |  |  |
| bcftools_query_qual_filter | vcf_plain | ByteExact | PASS | 9 | 7 |  |  |
| bcftools_query_info_dp_filtered | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_info_dp_excluded | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_gt_multi | vcf_multi_plain | ByteExact | PASS | 19 | 16 |  |  |
| bcftools_query_list_samples | vcf_multi_plain | ByteExact | PASS | 6 | 4 |  |  |
| bcftools_query_heavy_full | vcf_plain | ByteExact | PASS | 10 | 8 | 1.17x |  |
| bcftools_norm_split | vcf_plain | ByteExact | PASS | 21 | 13 |  |  |
| bcftools_norm_dedup | vcf_plain | ByteExact | PASS | 17 | 12 |  |  |
| bcftools_norm_check_ref | vcf_plain | ByteExact | PASS | 19 | 13 |  |  |
| filter_base | vcf_plain | ByteExact | PASS | 7 | 3 |  |  |
| filter_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_flagi_QUAL>=30 | vcf_plain | ByteExact | PASS | 4 | 3 |  |  |
| filter_combo_soft_lowqual | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_combo_snpgap | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_sort | vcf_plain | ByteExact | PASS | 16 | 16 |  |  |
| bcftools_stats | vcf_plain | ByteExact | PASS | 8 | 7 |  |  |
| bcftools_stats_multi | vcf_multi_plain | ByteExact | PASS | 19 | 10 |  |  |
| bcftools_head | vcf_plain | ByteExact | PASS | 6 | 3 |  |  |
| bcftools_annotate_drop_af | vcf_plain | ByteExact | PASS | 15 | 12 |  |  |
| bcftools_annotate_drop_dp | vcf_plain | ByteExact | PASS | 24 | 10 |  |  |
| bcftools_gtcheck | vcf_multi_plain | ByteExact | PASS | 120 | 19 |  |  |
| bcftools_mpileup | bam | ByteExact | PASS | 4515 | 1377 |  |  |
| bcftools_mpileup_heavy | bam | ByteExact | PASS | 4584 | 1433 | 3.20x |  |
| bcftools_plugin_fill_tags | vcf_multi_plain | ByteExact | PASS | 67 | 23 |  |  |
| bcftools_plugin_fill_tags_an_ac | vcf_multi_plain | ByteExact | PASS | 65 | 25 |  |  |
| bcftools_view_c_update_acan | vcf_multi_plain | ByteExact | PASS | 42 | 17 |  |  |
| bcftools_annotate_drop_filter | vcf_plain | ByteExact | PASS | 21 | 13 |  |  |
| bcftools_roh | vcf_multi_plain | ByteExact | PASS | 134 | 48 |  |  |
| bcftools_consensus | vcf | ByteExact | PASS | 189 | 14 |  |  |
| bcftools_concat | vcf | ByteExact | PASS | 33 | 24 |  |  |
| bcftools_norm_join | vcf_plain | ByteExact | PASS | 24 | 15 |  |  |
| bcftools_call | vcf_plain | Similarity | SIMILAR | 2399 | 1333 |  |  |
| bcftools_csq | vcf_plain | ByteExact | SKIP | 0 | 0 |  | bcftools csq --force over the multi-gene fixture is byte-identical to upstream EXCEPT for 3 records, and that residual i... |
| bcftools_isec | vcf | ByteExact | PASS | 608 | 40 |  |  |
| bcftools_merge | vcf_plain | ByteExact | PASS | 46 | 30 |  |  |
| bcftools_convert | vcf_plain | ByteExact | PASS | 24 | 15 |  |  |
| bcftools_reheader | vcf_multi_plain | BGZFDecoded | PASS | 49 | 9 |  |  |
| view_base | vcf_plain | ByteExact | PASS | 25 | 12 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 9 | 8 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 9 | 7 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 10 | 8 |  |  |
| query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 9 | 6 |  |  |
| query_gt | vcf_plain | ByteExact | PASS | 8 | 8 |  |  |
| query_info_dp | vcf_plain | ByteExact | PASS | 8 | 7 |  |  |

## bed12tobed6

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bed12tobed6_score_dropped | bed12 | ByteExact | PASS | 15 | 12 |  |  |

## bedannotate

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedannotate_default_header_order | bed | ByteExact | PASS | 52 | 105 |  |  |

## bedbamtobed

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bam | ByteExact | PASS | 100 | 78 |  |  |
| flagsplit | bam | ByteExact | PASS | 93 | 89 |  |  |
| flaged | bam | ByteExact | PASS | 5 | 2 |  |  |
| flagcigar | bam | ByteExact | PASS | 89 | 86 |  |  |
| combo_split_cigar | bam | ByteExact | PASS | 84 | 89 |  |  |

## bedclosest

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 179 | 32 |  |  |
| flagd | bed | ByteExact | PASS | 189 | 47 |  |  |
| flagio | bed | ByteExact | PASS | 118 | 9 |  |  |
| flagiu | bed | ByteExact | PASS | 4 | 2 |  |  |
| flagt_first | bed | ByteExact | PASS | 118 | 10 |  |  |
| flagt_last | bed | ByteExact | PASS | 119 | 9 |  |  |
| flagt_all | bed | ByteExact | PASS | 199 | 30 |  |  |
| flags | bed | ByteExact | PASS | 162 | 20 |  |  |
| flagN | bed | ByteExact | PASS | 165 | 28 |  |  |
| combo_d_t_first | bed | ByteExact | PASS | 121 | 10 |  |  |

## bedcluster

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 29 | 12 |  |  |
| flagd_50 | bed | ByteExact | PASS | 5 | 6 |  |  |
| flagd_0 | bed | ByteExact | PASS | 5 | 6 |  |  |
| bedcluster_s | bed | ByteExact | PASS | 8 | 7 |  |  |
| bedcluster_d50_s | bed | ByteExact | PASS | 6 | 7 |  |  |

## bedcomplement

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedcomplement_base | bed | ByteExact | PASS | 13 | 3 |  |  |

## bedcoverage

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 42 | 39 |  |  |
| flagcounts | bed | ByteExact | PASS | 27 | 34 |  |  |
| flagd | bed | ByteExact | PASS | 355 | 2535 |  |  |
| flaghist | bed | ByteExact | PASS | 110 | 94 |  |  |
| flags | bed | ByteExact | PASS | 26 | 32 |  |  |
| flagS | bed | ByteExact | PASS | 22 | 29 |  |  |
| flagmean | bed | ByteExact | PASS | 31 | 41 |  |  |
| combo_counts_s | bed | ByteExact | PASS | 24 | 28 |  |  |

## bedexpand

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedexpand_c5 | bed | ByteExact | PASS | 14 | 6 |  |  |
| bedexpand_trailing_comma_col11 | bed12 | ByteExact | PASS | 9 | 11 |  |  |

## bedfisher

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedfisher_overlap_count | bed | ByteExact | PASS | 19 | 6 |  |  |

## bedflank

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 9 | 1 |  |  |
| flagb_50 | bed | ByteExact | PASS | 5 | 6 |  |  |
| flags | bed | ByteExact | PASS | 3 | 1 |  |  |
| combo_l_r | bed | ByteExact | PASS | 6 | 7 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 6 | 7 |  |  |
| combo_b_s | bed | ByteExact | PASS | 6 | 7 |  |  |

## bedgenomecov

PASS 5 · SIMILAR 3 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | Similarity | SIMILAR | 37 | 21 |  |  |
| flagbg | bed | ByteExact | PASS | 11 | 14 |  |  |
| flagbga | bed | ByteExact | PASS | 12 | 13 |  |  |
| flagd | bed | ByteExact | PASS | 64 | 537 |  |  |
| flagdz | bed | ByteExact | PASS | 44 | 529 |  |  |
| flagmax_5 | bed | Similarity | SIMILAR | 23 | 19 |  |  |
| flagstrand_+ | bed | Similarity | SIMILAR | 22 | 18 |  |  |
| combo_bg_strand | bed | ByteExact | PASS | 10 | 13 |  |  |

## bedgetfasta

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | fasta | ByteExact | PASS | 53 | 64 |  |  |
| flags | fasta | ByteExact | PASS | 53 | 82 |  |  |
| flagname | fasta | ByteExact | PASS | 25 | 67 |  |  |
| flagtab | fasta | ByteExact | PASS | 40 | 60 |  |  |
| combo_s_name | fasta | ByteExact | PASS | 55 | 97 |  |  |
| combo_s_tab | fasta | ByteExact | PASS | 51 | 76 |  |  |

## bedgroupby

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 11 | 1 |  |  |
| flago_mean | bed | ByteExact | PASS | 5 | 3 |  |  |
| flago_sum | bed | ByteExact | PASS | 6 | 4 |  |  |
| flago_min | bed | ByteExact | PASS | 10 | 4 |  |  |
| flago_max | bed | ByteExact | PASS | 6 | 4 |  |  |
| flago_count | bed | ByteExact | PASS | 5 | 3 |  |  |
| combo_g1_c5_mean | bed | ByteExact | PASS | 5 | 3 |  |  |
| combo_g1_c5_count | bed | ByteExact | PASS | 5 | 4 |  |  |

## bedigv

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedigv_base | bed | ByteExact | PASS | 14 | 19 |  |  |

## bedintersect

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| intersect_base | bed | ByteExact | PASS | 38 | 32 |  |  |
| intersect_flagc | bed | ByteExact | PASS | 12 | 15 |  |  |
| intersect_flagv | bed | ByteExact | PASS | 11 | 14 |  |  |
| intersect_flagu | bed | ByteExact | PASS | 11 | 15 |  |  |
| intersect_flagwa | bed | ByteExact | PASS | 23 | 30 |  |  |
| intersect_flagwb | bed | ByteExact | PASS | 69 | 74 |  |  |
| intersect_flags | bed | ByteExact | PASS | 18 | 22 |  |  |
| intersect_combo_wa_wb | bed | ByteExact | PASS | 59 | 64 | 0.92x |  |
| intersect_combo_u_s | bed | ByteExact | PASS | 13 | 16 |  |  |
| intersect_combo_c_s | bed | ByteExact | PASS | 13 | 16 |  |  |

## bedjaccard

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedjaccard_base | bed | ByteExact | PASS | 18 | 5 |  |  |
| bedjaccard_s | bed | ByteExact | PASS | 6 | 6 |  |  |
| bedjaccard_f50 | bed | ByteExact | PASS | 7 | 6 |  |  |

## bedlinks

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedlinks_base | bed | ByteExact | PASS | 18 | 37 |  |  |

## bedmakewindows

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmakewindows_g_w_winnum | bed | ByteExact | PASS | 13 | 4 |  |  |
| bedmakewindows_g_ws_winnum | bed | ByteExact | PASS | 5 | 3 |  |  |
| bedmakewindows_b_winnum | bed | ByteExact | PASS | 10 | 15 |  |  |
| bedmakewindows_b_srcwinnum | bed | ByteExact | PASS | 11 | 15 |  |  |
| bedmakewindows_b_n_winnum | bed | ByteExact | PASS | 10 | 16 |  |  |
| bedmakewindows_default_none | bed | ByteExact | PASS | 4 | 2 |  |  |
| bedmakewindows_i_src | bed | ByteExact | PASS | 3 | 2 |  |  |

## bedmap

PASS 11 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 27 | 14 |  |  |
| flago_mean | bed | ByteExact | PASS | 19 | 15 |  |  |
| flago_sum | bed | ByteExact | PASS | 18 | 14 |  |  |
| flago_min | bed | ByteExact | PASS | 18 | 14 |  |  |
| flago_max | bed | ByteExact | PASS | 15 | 14 |  |  |
| flago_count | bed | ByteExact | PASS | 16 | 9 |  |  |
| flago_median | bed | ByteExact | PASS | 16 | 15 |  |  |
| flago_stdev | bed | ByteExact | PASS | 16 | 16 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 16 | 14 |  |  |
| combo_c5_sum_s | bed | ByteExact | PASS | 15 | 12 |  |  |
| bedmap_collapse_tiebreak | bed | ByteExact | PASS | 17 | 10 |  |  |

## bedmerge

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 14 | 3 |  |  |
| flagd_50 | bed | ByteExact | PASS | 5 | 3 |  |  |
| flagd_0 | bed | ByteExact | PASS | 6 | 3 |  |  |
| flags | bed | ByteExact | PASS | 7 | 3 |  |  |
| flagS_+ | bed | ByteExact | PASS | 6 | 3 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 5 | 3 |  |  |
| combo_c5_count | bed | ByteExact | PASS | 5 | 3 |  |  |
| combo_d50_c5_sum | bed | ByteExact | PASS | 6 | 3 |  |  |
| combo_s_c5_mean | bed | ByteExact | PASS | 7 | 3 |  |  |
| bedmerge_collapse_tiebreak | bed | ByteExact | PASS | 5 | 3 |  |  |

## bedmulticov

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmulticov_one_bam | bam | ByteExact | PASS | 97 | 3571 |  |  |
| bedmulticov_mapq20 | bam | ByteExact | PASS | 85 | 3719 |  |  |

## bedmultiinter

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedmultiinter_two | bed | ByteExact | PASS | 17 | 7 |  |  |
| bedmultiinter_two_names | bed | ByteExact | PASS | 7 | 7 |  |  |

## bednuc

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bednuc_base | fasta | ByteExact | PASS | 58 | 57 |  |  |

## bedoverlap

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedoverlap_base | bed | ByteExact | PASS | 18 | 2 |  |  |

## bedpairtobed

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedpairtobed_base | bed | ByteExact | PASS | 45 | 22 |  |  |

## bedpairtopair

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedpairtopair_base | bed | ByteExact | PASS | 16 | 8 |  |  |

## bedrandom

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedrandom_n50_l100_seed | bed | ByteExact | PASS | 10 | 1 |  |  |
| bedrandom_n100_l500_seed | bed | ByteExact | PASS | 3 | 1 |  |  |

## bedreldist

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedreldist_base | bed | ByteExact | PASS | 16 | 7 |  |  |
| bedreldist_detail | bed | ByteExact | PASS | 19 | 9 |  |  |

## bedsample

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsample_n50_seed | bed | ByteExact | PASS | 12 | 4 |  |  |
| bedsample_n200_seed | bed | ByteExact | PASS | 4 | 3 |  |  |

## bedshift

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 9 | 1 |  |  |
| flags_50 | bed | ByteExact | PASS | 4 | 5 |  |  |
| flags_50 | bed | ByteExact | PASS | 5 | 5 |  |  |
| combo_p_m | bed | ByteExact | PASS | 5 | 5 |  |  |
| combo_s_pct | bed | ByteExact | PASS | 5 | 5 |  |  |

## bedshuffle

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedshuffle_seed | bed | ByteExact | PASS | 13 | 5 |  |  |
| bedshuffle_seed_chrom | bed | ByteExact | PASS | 5 | 5 |  |  |

## bedslop

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 10 | 2 |  |  |
| flagb_50 | bed | ByteExact | PASS | 5 | 6 |  |  |
| flagb_100 | bed | ByteExact | PASS | 5 | 5 |  |  |
| flagpct | bed | ByteExact | PASS | 3 | 1 |  |  |
| flags | bed | ByteExact | PASS | 3 | 1 |  |  |
| combo_l_r | bed | ByteExact | PASS | 5 | 5 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 5 | 5 |  |  |
| combo_b_s | bed | ByteExact | PASS | 5 | 5 |  |  |

## bedsort

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsort_sizeA | bed | ByteExact | PASS | 25 | 9 |  |  |
| bedsort_chrThenSizeA | bed | ByteExact | PASS | 10 | 7 |  |  |
| bedsort_default_tiebreak | bed | ByteExact | PASS | 8 | 6 |  |  |
| bedsort_sizeD_tiebreak | bed | ByteExact | PASS | 9 | 9 |  |  |

## bedspacing

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedspacing_base | bed | ByteExact | PASS | 11 | 4 |  |  |

## bedsplit

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsplit_simple_n3 | bed | ByteExact | PASS | 21 | 10 |  |  |
| bedsplit_size_n3 | bed | ByteExact | PASS | 14 | 13 |  |  |

## bedsubtract

PASS 9 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 31 | 20 |  |  |
| flagA | bed | ByteExact | PASS | 21 | 15 |  |  |
| flagN | bed | ByteExact | PASS | 3 | 2 |  |  |
| flags | bed | ByteExact | PASS | 58 | 18 |  |  |
| flagS | bed | ByteExact | PASS | 62 | 18 |  |  |
| flagf_0_5 | bed | ByteExact | PASS | 21 | 18 |  |  |
| combo_A_s | bed | ByteExact | PASS | 57 | 14 |  |  |
| combo_N_f | bed | ByteExact | PASS | 22 | 18 |  |  |
| bedsubtract_reciprocal_r_unsupported | bed | ByteExact | PASS | 21 | 19 |  |  |

## bedsummary

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsummary_format_and_missing_g | bed | ByteExact | PASS | 19 | 4 |  |  |

## bedtag

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtag_tag | bam | BAMDecoded | PASS | 297 | 441 |  |  |

## bedtobam

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtobam_decoded | bed | BAMDecoded | PASS | 23 | 9 |  |  |

## bedunionbedg

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedunionbedg_base | bed | ByteExact | PASS | 12 | 7 |  |  |

## bedwindow

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedwindow_v_w100 | bed | ByteExact | PASS | 21 | 17 |  |  |
| bedwindow_c_w100 | bed | ByteExact | PASS | 13 | 19 |  |  |
| bedwindow_v_lr | bed | ByteExact | PASS | 13 | 17 |  |  |
| bedwindow_join_w100 | bed | ByteExact | PASS | 22 | 48 |  |  |

## bgzip

PASS 4 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bgzip_decompress | vcf | ByteExact | PASS | 8 | 5 |  |  |
| bgzip_decompress_heavy | vcf | ByteExact | PASS | 4 | 6 | 0.64x |  |
| bgzip_compress | vcf_plain | BGZFDecoded | PASS | 6 | 18 |  |  |
| bgzip_reindex | vcf | ByteExact | PASS | 2 | 4 |  |  |

## fastp

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| fastp_cut_tail | fastq | ByteExact | PASS | 404 | 368 |  |  |
| fastp_default_filter | fastq | ByteExact | PASS | 388 | 93 |  |  |
| fastp_detect_adapter_pe_heavy | fastq_paired | ByteExact | PASS | 3021 | 586 | 5.15x |  |

## htsfile

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| htsfile_identify | vcf | ByteExact | PASS | 17 | 4 |  |  |
| htsfile_copy | vcf | ByteExact | PASS | 13 | 9 |  |  |

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
| prinseq_min_len | fastq | ByteExact | PASS | 187 | 346 |  |  |
| prinseq_max_len | fastq | ByteExact | PASS | 180 | 340 |  |  |
| prinseq_trim_left | fastq | ByteExact | PASS | 191 | 377 |  |  |
| prinseq_trim_right | fastq | ByteExact | PASS | 188 | 357 |  |  |
| prinseq_trim_qual_right | fastq | ByteExact | PASS | 227 | 1090 |  |  |
| prinseq_trim_qual_left | fastq | ByteExact | PASS | 216 | 966 |  |  |
| prinseq_min_qual_mean | fastq | ByteExact | PASS | 197 | 920 |  |  |
| prinseq_max_ns | fastq | ByteExact | PASS | 174 | 327 |  |  |

## samtools

PASS 72 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_bam_base | bam | ByteExact | PASS | 136 | 63 |  |  |
| view_bam_flagh | bam | ByteExact | PASS | 126 | 62 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 5 | 4 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 76 | 43 |  |  |
| view_bam_flagq_20 | bam | ByteExact | PASS | 127 | 58 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 128 | 63 |  |  |
| view_bam_flagf_x2 | bam | ByteExact | PASS | 74 | 46 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 103 | 50 |  |  |
| view_bam_flagF_x100 | bam | ByteExact | PASS | 121 | 59 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 95 | 60 |  |  |
| view_bam_flagL_bed | bam | ByteExact | PASS | 159 | 59 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 74 | 42 |  |  |
| view_bam_combo_count_L | bam | ByteExact | PASS | 90 | 43 |  |  |
| view_bam_combo_f2_F256_q20_count | bam | ByteExact | PASS | 71 | 43 |  |  |
| view_bam_combo_header_plus_body | bam | ByteExact | PASS | 148 | 61 |  |  |
| samtools_view_region_contig | bam | ByteExact | PASS | 40 | 22 |  |  |
| samtools_view_region_range | bam | ByteExact | PASS | 47 | 4 |  |  |
| samtools_view_region_count | bam | ByteExact | PASS | 43 | 14 |  |  |
| samtools_view_cram_body | cram | ByteExact | PASS | 331 | 58 |  |  |
| samtools_view_cram_header | cram | ByteExact | PASS | 5 | 4 |  |  |
| samtools_view_cram_count | cram | ByteExact | PASS | 274 | 10 |  |  |
| samtools_view_cram_q30 | cram | ByteExact | PASS | 329 | 51 |  |  |
| samtools_view_cram_decode_sam_heavy | cram | ByteExact | PASS | 323 | 46 | 6.97x |  |
| view_subsample_seed | bam | ByteExact | PASS | 108 | 59 |  |  |
| samtools_sort_name_sam | bam | ByteExact | PASS | 151 | 77 |  |  |
| samtools_sort_byname_tag | bam | ByteExact | PASS | 139 | 74 |  |  |
| samtools_sort_name_sam_heavy | bam | ByteExact | PASS | 172 | 82 | 2.10x |  |
| samtools_sort_coord_sam | bam | ByteExact | PASS | 123 | 66 |  |  |
| samtools_flagstat | bam | ByteExact | PASS | 57 | 43 |  |  |
| samtools_idxstats | bam | ByteExact | PASS | 6 | 4 |  |  |
| samtools_stats | bam | ByteExact | PASS | 134 | 115 |  |  |
| samtools_quickcheck | bam | ByteExact | PASS | 5 | 4 |  |  |
| samtools_dict | fasta | ByteExact | PASS | 12 | 7 |  |  |
| samtools_stats_heavy | bam | ByteExact | PASS | 133 | 109 | 1.22x |  |
| depth_base | bam | ByteExact | PASS | 102 | 78 |  |  |
| depth_flaga | bam | ByteExact | PASS | 81 | 81 |  |  |
| depth_flagr_chr1 | bam | ByteExact | PASS | 54 | 22 |  |  |
| depth_flagb_bed | bam | ByteExact | PASS | 104 | 105 |  |  |
| depth_combo_all_region | bam | ByteExact | PASS | 42 | 4 |  |  |
| depth_combo_all_bed | bam | ByteExact | PASS | 73 | 96 |  |  |
| samtools_coverage | bam | ByteExact | PASS | 708 | 81 |  |  |
| samtools_coverage_region | bam | ByteExact | PASS | 547 | 24 |  |  |
| samtools_depth_mapq_filter | bam | ByteExact | PASS | 77 | 82 |  |  |
| samtools_depth_baseq_filter | bam | ByteExact | PASS | 89 | 81 |  |  |
| samtools_calmd | bam | ByteExact | PASS | 151 | 72 |  |  |
| samtools_consensus | bam | ByteExact | PASS | 1667 | 304 |  |  |
| samtools_consensus_region | bam | ByteExact | PASS | 168 | 5 |  |  |
| samtools_fastq | bam | ByteExact | PASS | 622 | 86 |  |  |
| samtools_fastq_n | bam | ByteExact | PASS | 158 | 82 |  |  |
| samtools_fastq_heavy | bam | ByteExact | PASS | 150 | 82 | 1.82x |  |
| samtools_tview_text | bam | ByteExact | PASS | 206 | 78 |  |  |
| samtools_mpileup_pileup | bam | ByteExact | PASS | 2754 | 1111 |  |  |
| cat_concat | bam | BAMDecoded | PASS | 445 | 14 |  |  |
| samtools_markdup | bam | BAMDecoded | PASS | 329 | 303 |  |  |
| samtools_fixmate | bam | BAMDecoded | PASS | 271 | 385 |  |  |
| samtools_addreplacerg | bam | BAMDecoded | PASS | 243 | 65 |  |  |
| samtools_merge | bam | BAMDecoded | PASS | 352 | 339 |  |  |
| samtools_reheader | bam | BAMDecoded | PASS | 231 | 11 |  |  |
| samtools_split | bam | BAMDecoded | PASS | 252 | 314 |  |  |
| samtools_import | fastq | BAMDecoded | PASS | 214 | 30 |  |  |
| samtools_phase | bam | ByteExact | PASS | 1896 | 296 |  |  |
| view_bam_base | bam | ByteExact | PASS | 113 | 68 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 6 | 3 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 71 | 41 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 121 | 55 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 101 | 58 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 114 | 49 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 77 | 43 |  |  |
| view_bam_combo_header_and_reads | bam | ByteExact | PASS | 115 | 58 |  |  |
| view_cram_body | cram | ByteExact | PASS | 320 | 45 |  |  |
| view_cram_count | cram | ByteExact | PASS | 262 | 10 |  |  |
| view_cram_decode_sam_heavy | cram | ByteExact | PASS | 311 | 44 | 7.06x |  |

## seqtk

PASS 32 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| seqtk_seq_fq_base | fastq | ByteExact | PASS | 52 | 21 |  |  |
| seqtk_seq_fq_flagA | fastq | ByteExact | PASS | 16 | 22 |  |  |
| seqtk_seq_fq_flagr | fastq | ByteExact | PASS | 39 | 24 |  |  |
| seqtk_seq_fq_flagL_95 | fastq | ByteExact | PASS | 23 | 18 |  |  |
| seqtk_seq_fq_flagq_20 | fastq | ByteExact | PASS | 32 | 28 |  |  |
| seqtk_seq_fq_flagl_60 | fastq | ByteExact | PASS | 24 | 26 |  |  |
| seqtk_seq_fq_combo_rev_upper | fastq | ByteExact | PASS | 36 | 43 |  |  |
| seqtk_seq_fq_combo_qmask_n | fastq | ByteExact | PASS | 38 | 28 |  |  |
| seqtk_comp_fa | fasta | ByteExact | PASS | 9 | 5 |  |  |
| seqtk_comp_fq | fastq | ByteExact | PASS | 75 | 74 |  |  |
| seqtk_fqchk | fastq | ByteExact | PASS | 21 | 12 |  |  |
| seqtk_fqchk_q20 | fastq | ByteExact | PASS | 22 | 12 |  |  |
| seqtk_size_fq | fastq | ByteExact | PASS | 12 | 6 |  |  |
| seqtk_size_fa | fasta | ByteExact | PASS | 3 | 1 |  |  |
| seqtk_trimfq | fastq | ByteExact | PASS | 67 | 26 |  |  |
| seqtk_trimfq_q | fastq | ByteExact | PASS | 49 | 33 |  |  |
| seqtk_trimfq_be | fastq | ByteExact | PASS | 24 | 22 |  |  |
| seqtk_sample_count | fastq | ByteExact | PASS | 16 | 6 |  |  |
| seqtk_sample_frac | fastq | ByteExact | PASS | 16 | 8 |  |  |
| seqtk_hpc_fa | fasta | ByteExact | PASS | 5 | 4 |  |  |
| seqtk_hpc_fq | fastq | ByteExact | PASS | 55 | 35 |  |  |
| seqtk_gap_fa | fasta | ByteExact | PASS | 4 | 2 |  |  |
| seqtk_subseq_bed | fasta | ByteExact | PASS | 22 | 31 |  |  |
| seqtk_mergepe | fastq_paired | ByteExact | PASS | 75 | 36 |  |  |
| seqtk_dropse | fastq | ByteExact | PASS | 17 | 6 |  |  |
| seqtk_randbase | fasta | ByteExact | PASS | 9 | 4 |  |  |
| seqtk_telo | fasta | ByteExact | PASS | 4 | 1 |  |  |
| seqtk_listhet | fasta | ByteExact | PASS | 5 | 1 |  |  |
| seqtk_hety | fasta | ByteExact | PASS | 8 | 5 |  |  |
| seqtk_seq_fa | fasta | ByteExact | PASS | 4 | 1 |  |  |
| seqtk_seq_fq_to_fa_heavy | fastq | ByteExact | PASS | 17 | 17 | 1.00x |  |
| seqtk_cutN | fasta | ByteExact | PASS | 5 | 4 |  |  |

## sickle

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| sickle_se_base | fastq | ByteExact | PASS | 165 | 157 |  |  |
| sickle_se_q30 | fastq | ByteExact | PASS | 156 | 150 |  |  |
| sickle_se_l30 | fastq | ByteExact | PASS | 154 | 177 |  |  |
| sickle_se_no5prime | fastq | ByteExact | PASS | 157 | 150 |  |  |
| sickle_se_truncn | fastq | ByteExact | PASS | 130 | 165 |  |  |
| sickle_se_q30_l40 | fastq | ByteExact | PASS | 140 | 150 |  |  |
| sickle_pe_base | fastq_paired | ByteExact | PASS | 332 | 316 | 1.05x |  |
| sickle_se_cli_default_window | fastq | ByteExact | PASS | 123 | 149 |  |  |

## skewer

PASS 3 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| skewer_se_base | fastq | ByteExact | PASS | 291 | 300 |  |  |
| skewer_se_minlen30 | fastq | ByteExact | PASS | 271 | 304 |  |  |
| skewer_se_full_heavy | fastq | ByteExact | PASS | 277 | 303 | 0.91x |  |

## tabix

PASS 7 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| tabix_region_contig | vcf | ByteExact | PASS | 22 | 5 |  |  |
| tabix_region_range | vcf | ByteExact | PASS | 4 | 3 |  |  |
| tabix_region_chr2 | vcf | ByteExact | PASS | 9 | 4 |  |  |
| tabix_region_with_header | vcf | ByteExact | PASS | 9 | 3 |  |  |
| tabix_list_chroms | vcf | ByteExact | PASS | 4 | 3 |  |  |
| tabix_region_heavy | vcf | ByteExact | PASS | 11 | 3 | 2.99x |  |
| tabix_regions_bed | vcf | ByteExact | PASS | 3504 | 68 |  |  |

## vcftools

PASS 22 · SIMILAR 0 · DIVERGE 3 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| vcftools_freq | vcf_plain | ByteExact | PASS | 283 | 15 |  |  |
| vcftools_counts | vcf_plain | ByteExact | PASS | 207 | 12 |  |  |
| vcftools_freq2 | vcf_plain | ByteExact | PASS | 229 | 14 |  |  |
| vcftools_depth | vcf_plain | ByteExact | PASS | 11 | 7 |  |  |
| vcftools_site_depth | vcf_plain | ByteExact | PASS | 209 | 268 |  |  |
| vcftools_site_mean_depth | vcf_plain | ByteExact | DIVERGE | 287 | 207 |  | output file ".ldepth.mean": first diff at line 2:   ours:     chr1	295	53	-nan   upstream: chr1	295	53	nan |
| vcftools_site_pi | vcf_plain | ByteExact | PASS | 212 | 204 |  |  |
| vcftools_window_pi | vcf_plain | ByteExact | PASS | 52 | 32 |  |  |
| vcftools_tstv_summary | vcf_plain | ByteExact | PASS | 11 | 6 |  |  |
| vcftools_missing_indv | vcf_plain | ByteExact | PASS | 11 | 7 |  |  |
| vcftools_missing_site | vcf_plain | ByteExact | PASS | 226 | 211 |  |  |
| vcftools_het | vcf_plain | ByteExact | PASS | 14 | 10 |  |  |
| vcftools_singletons | vcf_plain | ByteExact | PASS | 277 | 262 |  |  |
| vcftools_recode_heavy | vcf_plain | ByteExact | PASS | 17 | 23 |  |  |
| vcftools_het_multi | vcf_multi_plain | ByteExact | PASS | 3 | 1 |  |  |
| vcftools_relatedness | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_relatedness2 | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_freq_multi | vcf_multi_plain | ByteExact | PASS | 3 | 1 |  |  |
| vcftools_missing_indv_multi | vcf_multi_plain | ByteExact | PASS | 3 | 1 |  |  |
| vcftools_window_pi_heavy | vcf_plain | ByteExact | PASS | 95 | 80 | 1.19x |  |
| vcftools_geno_r2 | vcf_multi_plain | ByteExact | PASS | 3 | 1 | 3.47x |  |
| vcftools_hap_r2 | vcf_multi_plain | ByteExact | DIVERGE | 243 | 5004 | 0.05x | output file ".hap.ld": first diff at line 105:   ours:     chr1	705	859	24	3.69779e-33	1.38778e-17	1.11022e-16   upstrea... |
| vcftools_matrix012 | vcf_multi_plain | ByteExact | PASS | 4 | 1 |  |  |
| vcftools_lroh | vcf_multi_plain | ByteExact | PASS | 3 | 1 |  |  |
| vcftools_hardy | vcf_plain | ByteExact | DIVERGE | 219 | 218 |  | output file ".hwe": first diff at line 2:   ours:     chr1	295	0/0/1	0.00/0.00/1.00	-nan	1.000000e+00	1.000000e+00	1.000... |

