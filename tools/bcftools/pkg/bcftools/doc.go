// Package bcftools is a pure-Go reimplementation of selected bcftools
// subcommands. The first slice covers `view`: VCF/BCF I/O with sample,
// region, and filter expression selection.
//
// Subsequent slices will add `query`, `stats`, `norm`, `concat`, and `merge`.
// BCF writing and .csi indexing are explicitly deferred — the view command
// reads BCF (via pkg/htsgo/bcf) but emits only VCF/VCF.gz so far.
package bcftools
