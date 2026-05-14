// Package bgzip implements the Blocked GNU Zip Format (BGZF) used by htslib
// and the samtools/bcftools/tabix ecosystem.
//
// A BGZF file is a concatenation of independent gzip member blocks. Each block
// has a small (≤64 KiB) uncompressed payload and carries the BGZF "BC" subfield
// in its gzip extra field, which records the total compressed length of the
// block minus one (BSIZE). The terminal block of a well-formed BGZF stream is
// a 28-byte gzip member representing an empty payload; readers use it to
// distinguish a complete stream from a truncated one.
//
// The package exposes a Writer that emits well-formed BGZF (including the EOF
// block on Close), a Reader that decodes a stream of BGZF blocks and verifies
// the BC subfield, and helpers (Scan, ReadGZI, WriteGZI, UncompressedOffsetAt,
// DecompressedSize) used by the bgzip CLI's -b/-s/-r subcommands and by tabix.
//
// References:
//   - SAM/BAM specification, section "The BGZF compression format".
//   - htslib htslib/bgzf.h, htslib/bgzf.c.
package bgzip
