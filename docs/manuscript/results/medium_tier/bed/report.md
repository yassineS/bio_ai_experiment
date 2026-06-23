# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-22T14:40:39Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 112 | 109 | 3 | 0 | 0 | 0 |

## bedbamtobed

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bam | ByteExact | PASS | 794 | 664 |  |  |
| flagsplit | bam | ByteExact | PASS | 719 | 797 |  |  |
| flaged | bam | ByteExact | PASS | 7 | 2 |  |  |
| flagcigar | bam | ByteExact | PASS | 668 | 652 |  |  |
| combo_split_cigar | bam | ByteExact | PASS | 642 | 650 |  |  |

## bedclosest

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 3072 | 161 |  |  |
| flagd | bed | ByteExact | PASS | 3069 | 170 |  |  |
| flagio | bed | ByteExact | PASS | 2981 | 64 |  |  |
| flagiu | bed | ByteExact | PASS | 6 | 3 |  |  |
| flagt_first | bed | ByteExact | PASS | 2994 | 49 |  |  |
| flagt_last | bed | ByteExact | PASS | 2961 | 46 |  |  |
| flagt_all | bed | ByteExact | PASS | 3045 | 121 |  |  |
| flags | bed | ByteExact | PASS | 3005 | 84 |  |  |
| flagN | bed | ByteExact | PASS | 3336 | 95 |  |  |
| combo_d_t_first | bed | ByteExact | PASS | 2963 | 54 |  |  |

## bedcluster

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 25 | 34 |  |  |
| flagd_50 | bed | ByteExact | PASS | 26 | 32 |  |  |
| flagd_0 | bed | ByteExact | PASS | 23 | 30 |  |  |
| bedcluster_s | bed | ByteExact | PASS | 37 | 42 |  |  |
| bedcluster_d50_s | bed | ByteExact | PASS | 38 | 41 |  |  |

## bedcoverage

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 118 | 163 |  |  |
| flagcounts | bed | ByteExact | PASS | 107 | 135 |  |  |
| flagd | bed | ByteExact | PASS | 7248 | 24957 |  |  |
| flaghist | bed | ByteExact | PASS | 781 | 841 |  |  |
| flags | bed | ByteExact | PASS | 171 | 163 |  |  |
| flagS | bed | ByteExact | PASS | 136 | 166 |  |  |
| flagmean | bed | ByteExact | PASS | 154 | 183 |  |  |
| combo_counts_s | bed | ByteExact | PASS | 103 | 159 |  |  |

## bedflank

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 3 | 1 |  |  |
| flagb_50 | bed | ByteExact | PASS | 29 | 40 |  |  |
| flags | bed | ByteExact | PASS | 4 | 1 |  |  |
| combo_l_r | bed | ByteExact | PASS | 29 | 40 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 29 | 40 |  |  |
| combo_b_s | bed | ByteExact | PASS | 28 | 41 |  |  |

## bedgenomecov

PASS 5 · SIMILAR 3 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | Similarity | SIMILAR | 365 | 274 |  |  |
| flagbg | bed | ByteExact | PASS | 91 | 91 |  |  |
| flagbga | bed | ByteExact | PASS | 68 | 84 |  |  |
| flagd | bed | ByteExact | PASS | 1043 | 10506 |  |  |
| flagdz | bed | ByteExact | PASS | 2062 | 8089 |  |  |
| flagmax_5 | bed | Similarity | SIMILAR | 374 | 252 |  |  |
| flagstrand_+ | bed | Similarity | SIMILAR | 249 | 240 |  |  |
| combo_bg_strand | bed | ByteExact | PASS | 52 | 61 |  |  |

## bedgetfasta

PASS 6 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | fasta | ByteExact | PASS | 253 | 632 |  |  |
| flags | fasta | ByteExact | PASS | 295 | 554 |  |  |
| flagname | fasta | ByteExact | PASS | 136 | 468 |  |  |
| flagtab | fasta | ByteExact | PASS | 144 | 433 |  |  |
| combo_s_name | fasta | ByteExact | PASS | 315 | 580 |  |  |
| combo_s_tab | fasta | ByteExact | PASS | 319 | 575 |  |  |

## bedgroupby

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 2 | 3 |  |  |
| flago_mean | bed | ByteExact | PASS | 15 | 11 |  |  |
| flago_sum | bed | ByteExact | PASS | 14 | 11 |  |  |
| flago_min | bed | ByteExact | PASS | 17 | 11 |  |  |
| flago_max | bed | ByteExact | PASS | 12 | 11 |  |  |
| flago_count | bed | ByteExact | PASS | 11 | 10 |  |  |
| combo_g1_c5_mean | bed | ByteExact | PASS | 13 | 12 |  |  |
| combo_g1_c5_count | bed | ByteExact | PASS | 11 | 10 |  |  |

## bedintersect

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| intersect_base | bed | ByteExact | PASS | 75 | 107 |  |  |
| intersect_flagc | bed | ByteExact | PASS | 50 | 70 |  |  |
| intersect_flagv | bed | ByteExact | PASS | 43 | 63 |  |  |
| intersect_flagu | bed | ByteExact | PASS | 49 | 69 |  |  |
| intersect_flagwa | bed | ByteExact | PASS | 87 | 112 |  |  |
| intersect_flagwb | bed | ByteExact | PASS | 182 | 223 |  |  |
| intersect_flags | bed | ByteExact | PASS | 66 | 89 |  |  |
| intersect_combo_wa_wb | bed | ByteExact | PASS | 98 | 207 | 0.47x |  |
| intersect_combo_u_s | bed | ByteExact | PASS | 49 | 70 |  |  |
| intersect_combo_c_s | bed | ByteExact | PASS | 47 | 70 |  |  |

## bedmap

PASS 11 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 70 | 71 |  |  |
| flago_mean | bed | ByteExact | PASS | 110 | 97 |  |  |
| flago_sum | bed | ByteExact | PASS | 97 | 80 |  |  |
| flago_min | bed | ByteExact | PASS | 126 | 102 |  |  |
| flago_max | bed | ByteExact | PASS | 129 | 74 |  |  |
| flago_count | bed | ByteExact | PASS | 66 | 46 |  |  |
| flago_median | bed | ByteExact | PASS | 92 | 82 |  |  |
| flago_stdev | bed | ByteExact | PASS | 75 | 96 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 92 | 84 |  |  |
| combo_c5_sum_s | bed | ByteExact | PASS | 62 | 67 |  |  |
| bedmap_collapse_tiebreak | bed | ByteExact | PASS | 72 | 54 |  |  |

## bedmerge

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 27 | 16 |  |  |
| flagd_50 | bed | ByteExact | PASS | 37 | 12 |  |  |
| flagd_0 | bed | ByteExact | PASS | 24 | 11 |  |  |
| flags | bed | ByteExact | PASS | 34 | 14 |  |  |
| flagS_+ | bed | ByteExact | PASS | 27 | 13 |  |  |
| combo_c5_mean | bed | ByteExact | PASS | 28 | 15 |  |  |
| combo_c5_count | bed | ByteExact | PASS | 29 | 12 |  |  |
| combo_d50_c5_sum | bed | ByteExact | PASS | 31 | 16 |  |  |
| combo_s_c5_mean | bed | ByteExact | PASS | 40 | 21 |  |  |
| bedmerge_collapse_tiebreak | bed | ByteExact | PASS | 32 | 13 |  |  |

## bedshift

PASS 5 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 3 | 1 |  |  |
| flags_50 | bed | ByteExact | PASS | 21 | 28 |  |  |
| flags_50 | bed | ByteExact | PASS | 20 | 28 |  |  |
| combo_p_m | bed | ByteExact | PASS | 19 | 28 |  |  |
| combo_s_pct | bed | ByteExact | PASS | 19 | 28 |  |  |

## bedslop

PASS 8 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 3 | 1 |  |  |
| flagb_50 | bed | ByteExact | PASS | 19 | 28 |  |  |
| flagb_100 | bed | ByteExact | PASS | 19 | 27 |  |  |
| flagpct | bed | ByteExact | PASS | 3 | 1 |  |  |
| flags | bed | ByteExact | PASS | 3 | 1 |  |  |
| combo_l_r | bed | ByteExact | PASS | 18 | 28 |  |  |
| combo_b_pct | bed | ByteExact | PASS | 20 | 27 |  |  |
| combo_b_s | bed | ByteExact | PASS | 21 | 27 |  |  |

## bedsplit

PASS 2 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedsplit_simple_n3 | bed | ByteExact | PASS | 28 | 44 |  |  |
| bedsplit_size_n3 | bed | ByteExact | PASS | 35 | 66 |  |  |

## bedsubtract

PASS 9 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| base | bed | ByteExact | PASS | 350 | 96 |  |  |
| flagA | bed | ByteExact | PASS | 331 | 60 |  |  |
| flagN | bed | ByteExact | PASS | 4 | 2 |  |  |
| flags | bed | ByteExact | PASS | 1341 | 84 |  |  |
| flagS | bed | ByteExact | PASS | 1460 | 103 |  |  |
| flagf_0_5 | bed | ByteExact | PASS | 318 | 82 |  |  |
| combo_A_s | bed | ByteExact | PASS | 1375 | 53 |  |  |
| combo_N_f | bed | ByteExact | PASS | 321 | 78 |  |  |
| bedsubtract_reciprocal_r_unsupported | bed | ByteExact | PASS | 317 | 82 |  |  |

## bedtag

PASS 1 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| bedtag_tag | bam | BAMDecoded | PASS | 1953 | 2859 |  |  |

