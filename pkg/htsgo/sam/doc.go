// Package sam implements readers and writers for the Sequence Alignment/Map
// formats — text SAM and its binary counterpart BAM — as defined by the
// hts-specs SAM/BAM specification.
//
// Both formats share the same logical record (header lines plus alignment
// records); the difference is purely on-disk encoding. SAM is tab-separated
// text, BAM is a binary stream wrapped in BGZF blocks. This package presents
// them through a common Header/Record model so that callers (notably
// tools/samtools) can mix and match.
//
// # Format auto-detection
//
// NewReader sniffs the first bytes of its input. If they look like a BGZF
// gzip member with the BAM magic ("BAM\1") just inside, it constructs a
// BAM reader on top of tools/bgzip/pkg/bgzip.Reader. Otherwise it falls back
// to the line-oriented text SAM reader.
//
// References:
//   - SAMv1 specification (https://samtools.github.io/hts-specs/SAMv1.pdf).
//   - htslib htslib/sam.h, htslib/sam.c.
package sam
