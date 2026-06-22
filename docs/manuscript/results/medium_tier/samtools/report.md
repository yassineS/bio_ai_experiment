# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-22T13:05:58Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 72 | 72 | 0 | 0 | 0 | 0 |

## samtools

PASS 72 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| view_bam_base | bam | ByteExact | PASS | 1164 | 643 |  |  |
| view_bam_flagh | bam | ByteExact | PASS | 942 | 413 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 6 | 4 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 528 | 304 |  |  |
| view_bam_flagq_20 | bam | ByteExact | PASS | 980 | 425 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 854 | 403 |  |  |
| view_bam_flagf_x2 | bam | ByteExact | PASS | 521 | 301 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 691 | 349 |  |  |
| view_bam_flagF_x100 | bam | ByteExact | PASS | 825 | 406 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 688 | 354 |  |  |
| view_bam_flagL_bed | bam | ByteExact | PASS | 871 | 403 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 537 | 306 |  |  |
| view_bam_combo_count_L | bam | ByteExact | PASS | 566 | 316 |  |  |
| view_bam_combo_f2_F256_q20_count | bam | ByteExact | PASS | 562 | 314 |  |  |
| view_bam_combo_header_plus_body | bam | ByteExact | PASS | 815 | 464 |  |  |
| samtools_view_region_contig | bam | ByteExact | PASS | 112 | 58 |  |  |
| samtools_view_region_range | bam | ByteExact | PASS | 173 | 5 |  |  |
| samtools_view_region_count | bam | ByteExact | PASS | 147 | 42 |  |  |
| samtools_view_cram_body | cram | ByteExact | PASS | 2385 | 344 |  |  |
| samtools_view_cram_header | cram | ByteExact | PASS | 7 | 4 |  |  |
| samtools_view_cram_count | cram | ByteExact | PASS | 2109 | 45 |  |  |
| samtools_view_cram_q30 | cram | ByteExact | PASS | 2408 | 294 |  |  |
| samtools_view_cram_decode_sam_heavy | cram | ByteExact | PASS | 2358 | 312 | 7.55x |  |
| view_subsample_seed | bam | ByteExact | PASS | 795 | 367 |  |  |
| samtools_sort_name_sam | bam | ByteExact | PASS | 1359 | 591 |  |  |
| samtools_sort_byname_tag | bam | ByteExact | PASS | 1417 | 566 |  |  |
| samtools_sort_name_sam_heavy | bam | ByteExact | PASS | 1397 | 718 | 1.94x |  |
| samtools_sort_coord_sam | bam | ByteExact | PASS | 918 | 448 |  |  |
| samtools_flagstat | bam | ByteExact | PASS | 444 | 320 |  |  |
| samtools_idxstats | bam | ByteExact | PASS | 7 | 4 |  |  |
| samtools_stats | bam | ByteExact | PASS | 959 | 852 |  |  |
| samtools_quickcheck | bam | ByteExact | PASS | 6 | 4 |  |  |
| samtools_dict | fasta | ByteExact | PASS | 118 | 61 |  |  |
| samtools_stats_heavy | bam | ByteExact | PASS | 944 | 854 | 1.11x |  |
| depth_base | bam | ByteExact | PASS | 936 | 866 |  |  |
| depth_flaga | bam | ByteExact | PASS | 1332 | 900 |  |  |
| depth_flagr_chr1 | bam | ByteExact | PASS | 467 | 116 |  |  |
| depth_flagb_bed | bam | ByteExact | PASS | 1527 | 957 |  |  |
| depth_combo_all_region | bam | ByteExact | PASS | 293 | 6 |  |  |
| depth_combo_all_bed | bam | ByteExact | PASS | 1622 | 973 |  |  |
| samtools_coverage | bam | ByteExact | PASS | 13742 | 686 |  |  |
| samtools_coverage_region | bam | ByteExact | PASS | 7892 | 90 |  |  |
| samtools_depth_mapq_filter | bam | ByteExact | PASS | 1034 | 1003 |  |  |
| samtools_depth_baseq_filter | bam | ByteExact | PASS | 1080 | 805 |  |  |
| samtools_calmd | bam | ByteExact | PASS | 1103 | 498 |  |  |
| samtools_consensus | bam | ByteExact | PASS | 14245 | 2805 |  |  |
| samtools_consensus_region | bam | ByteExact | PASS | 884 | 7 |  |  |
| samtools_fastq | bam | ByteExact | PASS | 1023 | 582 |  |  |
| samtools_fastq_n | bam | ByteExact | PASS | 1097 | 560 |  |  |
| samtools_fastq_heavy | bam | ByteExact | PASS | 1099 | 571 | 1.92x |  |
| samtools_tview_text | bam | ByteExact | PASS | 653 | 81 |  |  |
| samtools_mpileup_pileup | bam | ByteExact | PASS | 21528 | 9924 |  |  |
| cat_concat | bam | BAMDecoded | PASS | 3514 | 97 |  |  |
| samtools_markdup | bam | BAMDecoded | PASS | 2604 | 2514 |  |  |
| samtools_fixmate | bam | BAMDecoded | PASS | 1956 | 2927 |  |  |
| samtools_addreplacerg | bam | BAMDecoded | PASS | 1770 | 429 |  |  |
| samtools_merge | bam | BAMDecoded | PASS | 2716 | 2805 |  |  |
| samtools_reheader | bam | BAMDecoded | PASS | 1780 | 45 |  |  |
| samtools_split | bam | BAMDecoded | PASS | 1945 | 2622 |  |  |
| samtools_import | fastq | BAMDecoded | PASS | 1491 | 192 |  |  |
| samtools_phase | bam | ByteExact | PASS | 10714 | 2127 |  |  |
| view_bam_base | bam | ByteExact | PASS | 852 | 435 |  |  |
| view_bam_flagH | bam | ByteExact | PASS | 7 | 4 |  |  |
| view_bam_flagc | bam | ByteExact | PASS | 556 | 352 |  |  |
| view_bam_flagq_30 | bam | ByteExact | PASS | 864 | 402 |  |  |
| view_bam_flagf_x10 | bam | ByteExact | PASS | 730 | 366 |  |  |
| view_bam_flagF_x10 | bam | ByteExact | PASS | 702 | 364 |  |  |
| view_bam_combo_count_q30 | bam | ByteExact | PASS | 540 | 304 |  |  |
| view_bam_combo_header_and_reads | bam | ByteExact | PASS | 817 | 440 |  |  |
| view_cram_body | cram | ByteExact | PASS | 2491 | 335 |  |  |
| view_cram_count | cram | ByteExact | PASS | 2006 | 45 |  |  |
| view_cram_decode_sam_heavy | cram | ByteExact | PASS | 2426 | 333 | 7.28x |  |

