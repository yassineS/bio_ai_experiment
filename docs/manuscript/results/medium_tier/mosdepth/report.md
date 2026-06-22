# Parity pipeline report

- Scale: `medium`
- Seed: `1`
- Generated: 2026-06-22T13:45:22Z

## Summary

| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |
|------:|-----:|--------:|--------:|-----:|------:|
| 10 | 10 | 0 | 0 | 0 | 0 |

## mosdepth

PASS 10 · SIMILAR 0 · DIVERGE 0 · SKIP 0 · ERROR 0

| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |
|-------|-------|---------|--------|---------:|-------------:|------:|--------|
| mosdepth_default | bam | ByteExact | PASS | 1345 | 554 |  |  |
| mosdepth_fast_mode | bam | ByteExact | PASS | 1317 | 533 |  |  |
| mosdepth_mapq20 | bam | ByteExact | PASS | 1285 | 518 |  |  |
| mosdepth_flag | bam | ByteExact | PASS | 1317 | 525 |  |  |
| mosdepth_by_bed_regions | bam | ByteExact | PASS | 12351 | 831 |  |  |
| mosdepth_by_window_regions | bam | ByteExact | PASS | 10009 | 603 |  |  |
| mosdepth_by_bed_thresholds | bam | ByteExact | PASS | 12481 | 975 | 12.79x |  |
| mosdepth_default_heavy | bam | ByteExact | PASS | 1326 | 552 | 2.40x |  |
| mosdepth_by_summary_region_rows | bam | ByteExact | PASS | 12484 | 817 |  |  |
| mosdepth_by_region_dist | bam | ByteExact | PASS | 12235 | 804 |  |  |

