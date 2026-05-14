// Package samtools is a pure-Go reimplementation of selected samtools
// subcommands. The current slice covers `view` (filter/convert SAM↔BAM) and
// `flagstat` (alignment summary stats) on top of pkg/bioformats/sam.
//
// Future subcommands (sort, index, depth, fastq, mpileup) will land in
// follow-up PRs.
package samtools
