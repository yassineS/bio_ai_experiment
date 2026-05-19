// Package bcf provides a pure-Go decoder for BCF, the binary VCF format
// produced by htslib's bcftools.
//
// # Scope of this first slice
//
// This package is intentionally read-only. It implements:
//
//   - Parsing the BCF magic header and the VCF-style text header that follows
//     it. Header dictionaries (CHROM, FILTER, INFO, FORMAT) are extracted from
//     the meta-information lines so that record bodies can be decoded into
//     plain VCF text via the existing pkg/htsgo/vcf types.
//   - Sequential record decoding including the "shared" portion (CHROM, POS,
//     RLEN, QUAL, FILTER, INFO, ID, ref/alt alleles) and the per-sample
//     "individual" portion (FORMAT fields for every sample).
//   - The htslib "typed" value encoding: int8 / int16 / int32 / float / char
//     scalars, vectors, missing-value sentinels, and the small-prefix vs.
//     length-prefixed variants used for fields longer than 14 entries.
//
// Out of scope for this slice (deferred to a follow-up PR):
//
//   - Writing BCF (encoder).
//   - .csi indexing of BCF for random access (today the BAM .bai code lives
//     under tools/samtools and tabix under tools/tabix; both move into htsgo
//     in PRs D and E. .csi for BCF will come with the next bcftools slice).
//   - GVCF-specific REF/ALT conventions beyond what the standard decoder
//     already handles.
//
// # BCF layout
//
// A BCF file is BGZF-wrapped (so callers should hand the decoder a BGZF-
// decompressed stream from pkg/htsgo/iohelper). The on-disk layout is:
//
//	magic[5] = "BCF\2\2"
//	l_text   uint32          // length of the text header including trailing NUL
//	text[l_text] byte         // the VCF-style header (##fileformat=... and #CHROM line)
//
//	repeated records:
//	    l_shared uint32
//	    l_indiv  uint32
//	    shared[l_shared] byte
//	    indiv[l_indiv]   byte
//
// The shared portion is fixed-prefix:
//
//	chrom    int32   // index into the CHROM dictionary, or -1
//	pos      int32   // 0-based on the wire (the decoder converts to 1-based)
//	rlen     int32   // length of REF
//	qual     float32 // missing = 0x7F800001
//	n_info   uint16
//	n_allele uint16
//	(n_sample : uint24, n_fmt : uint8) packed in one little-endian uint32
//	id       typed string
//	refalt   typed string array of n_allele entries (REF first, then ALTs)
//	filter   typed int vector (dictionary indices)
//	info     n_info repetitions of (key dict idx, typed value)
//
// The indiv portion is n_fmt repetitions of (key dict idx, typed array of
// n_sample entries).
//
// The "typed" encoding uses a single descriptor byte of the form
// `(descriptor << 4) | size_class`. descriptor 1=int8, 2=int16, 3=int32,
// 5=float, 7=char. size_class 0..14 is a literal count of entries that
// follow; size_class 15 means "the next typed integer holds the real count".
// Per-type missing values are 0x80, 0x8000, 0x80000000, 0x7F800001 for
// int8/int16/int32/float respectively; char missing is encoded as the
// zero-length string.
package bcf
