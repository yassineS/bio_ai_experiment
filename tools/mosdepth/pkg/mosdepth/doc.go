// Package mosdepth is a pure-Go re-implementation of brentp's `mosdepth`
// per-base / per-region / per-window depth-of-coverage tool.
//
// The upstream tool is written in Nim and is widely used in clinical and
// research NGS pipelines, but its Nim toolchain is a real friction point in
// CI/conda environments. This package targets feature parity for the common
// modes:
//
//   - Per-base BED-gz output (default).
//   - Per-region summary over a user-supplied BED.
//   - Per-window summary over fixed-width windows (-b INT).
//   - Threshold proportions (-T 1,5,10) reporting the fraction of a region
//     at or above each integer depth.
//   - Global cumulative distribution file (.mosdepth.global.dist.txt).
//   - Per-chromosome summary file (.mosdepth.summary.txt).
//
// The runtime is single-threaded (the -t/--threads flag is accepted for CLI
// compatibility). Outputs are bgzipped using
// github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf and indexed
// with a .csi index built by
// github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix, matching upstream
// mosdepth, which writes a CSI (via htslib's tbx_index_build) alongside each
// bgzipped BED output. The CSI uses min_shift=14, depth=5 — the htslib
// default — so consumers that read CSI indexes (tabix, bcftools, htslib)
// work transparently.
//
// The algorithm is a single streaming pass over a coordinate-sorted BAM:
// for each record we accumulate +1/-1 events at the start and end of its
// reference-consuming runs, then sweep a rolling depth counter forward and
// emit per-position records as the cursor advances. This keeps memory
// bounded to one reference's-worth of events at a time.
package mosdepth
