# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-23T02:44:59Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 68 | 66 | 1 | 0 | 1 | 0 |

## bcftools

PASS 66 · SIMILAR 1 · DIVERGE 0 · SKIP 1 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_base | vcf_plain | ByteExact | PASS | 64 | 42 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 47 | 37 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 28 | 28 |  |  |
| view_flagV_indels | vcf_plain | ByteExact | PASS | 49 | 40 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 51 | 36 |  |  |
| view_flagi_DP>40 | vcf_plain | ByteExact | PASS | 6 | 3 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 47 | 37 |  |  |
| view_flagO_v | vcf_plain | ByteExact | PASS | 47 | 38 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 51 | 38 |  |  |
| view_combo_indels_dp | vcf_plain | ByteExact | PASS | 6 | 3 |  |  |
| view_combo_excl_snps_qual | vcf_plain | ByteExact | PASS | 28 | 29 |  |  |
| bcftools_view_r_region | vcf | ByteExact | PASS | 8 | 5 |  |  |
| bcftools_view_t_targets | vcf | ByteExact | PASS | 36 | 41 |  |  |
| bcftools_view_r_multi | vcf | ByteExact | PASS | 244 | 12 |  |  |
| bcftools_view_s_sample | vcf_multi_plain | ByteExact | PASS | 317 | 101 |  |  |
| bcftools_view_G_drop_gt | vcf_multi_plain | ByteExact | PASS | 250 | 84 |  |  |
| bcftools_view_c_minac | vcf_multi_plain | ByteExact | PASS | 226 | 103 |  |  |
| bcftools_view_C_maxac | vcf_multi_plain | ByteExact | PASS | 170 | 72 |  |  |
| bcftools_view_q_minaf | vcf_multi_plain | ByteExact | PASS | 230 | 103 |  |  |
| bcftools_view_Q_maxaf | vcf_multi_plain | ByteExact | PASS | 219 | 109 |  |  |
| bcftools_query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 31 | 22 |  |  |
| bcftools_query_info_dp | vcf_plain | ByteExact | PASS | 30 | 27 |  |  |
| bcftools_query_qual_filter | vcf_plain | ByteExact | PASS | 35 | 24 |  |  |
| bcftools_query_info_dp_filtered | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_info_dp_excluded | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_gt_multi | vcf_multi_plain | ByteExact | PASS | 142 | 105 |  |  |
| bcftools_query_list_samples | vcf_multi_plain | ByteExact | PASS | 5 | 3 |  |  |
| bcftools_query_heavy_full | vcf_plain | ByteExact | PASS | 37 | 32 | 1.15x |  |
| bcftools_norm_split | vcf_plain | ByteExact | PASS | 105 | 66 |  |  |
| bcftools_norm_dedup | vcf_plain | ByteExact | PASS | 91 | 54 |  |  |
| bcftools_norm_check_ref | vcf_plain | ByteExact | PASS | 105 | 64 |  |  |
| filter_base | vcf_plain | ByteExact | PASS | 5 | 3 |  |  |
| filter_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 4 | 3 |  |  |
| filter_flagi_QUAL>=30 | vcf_plain | ByteExact | PASS | 4 | 2 |  |  |
| filter_combo_soft_lowqual | vcf_plain | ByteExact | PASS | 4 | 2 |  |  |
| filter_combo_snpgap | vcf_plain | ByteExact | PASS | 4 | 3 |  |  |
| bcftools_sort | vcf_plain | ByteExact | PASS | 79 | 79 |  |  |
| bcftools_stats | vcf_plain | ByteExact | PASS | 29 | 24 |  |  |
| bcftools_stats_multi | vcf_multi_plain | ByteExact | PASS | 101 | 40 |  |  |
| bcftools_head | vcf_plain | ByteExact | PASS | 6 | 3 |  |  |
| bcftools_annotate_drop_af | vcf_plain | ByteExact | PASS | 76 | 38 |  |  |
| bcftools_annotate_drop_dp | vcf_plain | ByteExact | PASS | 74 | 41 |  |  |
| bcftools_gtcheck | vcf_multi_plain | ByteExact | PASS | 1454 | 123 |  |  |
| bcftools_mpileup | bam | ByteExact | PASS | 83201 | 27361 |  |  |
| bcftools_mpileup_heavy | bam | ByteExact | PASS | 91736 | 27922 | 3.29x |  |
| bcftools_plugin_fill_tags | vcf_multi_plain | ByteExact | PASS | 811 | 172 |  |  |
| bcftools_plugin_fill_tags_an_ac | vcf_multi_plain | ByteExact | PASS | 622 | 290 |  |  |
| bcftools_view_c_update_acan | vcf_multi_plain | ByteExact | PASS | 380 | 275 |  |  |
| bcftools_annotate_drop_filter | vcf_plain | ByteExact | PASS | 160 | 68 |  |  |
| bcftools_roh | vcf_multi_plain | ByteExact | PASS | 2494 | 469 |  |  |
| bcftools_consensus | vcf | ByteExact | PASS | 7482 | 95 |  |  |
| bcftools_concat | vcf | ByteExact | PASS | 187 | 122 |  |  |
| bcftools_norm_join | vcf_plain | ByteExact | PASS | 98 | 69 |  |  |
| bcftools_call | vcf_plain | Similarity | SIMILAR | 43540 | 28295 |  |  |
| bcftools_csq | vcf_plain | ByteExact | SKIP | 0 | 0 |  | bcftools csq --force over the multi-gene fixture is byte-identical to upstream EXCEPT for 3 records, and that residual i... |
| bcftools_isec | vcf | ByteExact | PASS | 5098 | 294 |  |  |
| bcftools_merge | vcf_plain | ByteExact | PASS | 320 | 149 |  |  |
| bcftools_convert | vcf_plain | ByteExact | PASS | 83 | 57 |  |  |
| bcftools_reheader | vcf_multi_plain | BGZFDecoded | PASS | 439 | 14 |  |  |
| view_base | vcf_plain | ByteExact | PASS | 46 | 40 |  |  |
| view_flagv_snps | vcf_plain | ByteExact | PASS | 45 | 46 |  |  |
| view_flagv_indels | vcf_plain | ByteExact | PASS | 32 | 32 |  |  |
| view_flagi_QUAL>30 | vcf_plain | ByteExact | PASS | 49 | 39 |  |  |
| view_flage_QUAL<30 | vcf_plain | ByteExact | PASS | 47 | 40 |  |  |
| view_combo_snps_qual | vcf_plain | ByteExact | PASS | 48 | 42 |  |  |
| query_chrom_pos_ref_alt | vcf_plain | ByteExact | PASS | 35 | 23 |  |  |
| query_gt | vcf_plain | ByteExact | PASS | 33 | 37 |  |  |
| query_info_dp | vcf_plain | ByteExact | PASS | 32 | 27 |  |  |

