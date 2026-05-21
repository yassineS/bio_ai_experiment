// Package codec implements the custom compression codecs CRAM relies on
// that have no Go standard-library equivalent: the rANS (range
// Asymmetric Numeral System) entropy coders.
//
// CRAM blocks may be compressed with one of several methods. The ones
// covered (or planned) here:
//
//   - rANS 4x8  — CRAM v3.0's entropy coder. Order-0 and order-1.
//     Implemented in rans4x8.go. Pure Go, in-tree.
//   - rANS 4x16 — CRAM v3.1's entropy coder, with the v3.1 transform
//     bits (RLE, bit-packing, 4-way striping) and the 32-way coder.
//   - arith_dynamic — CRAM v3.1's adaptive range coder (compression
//     method 6). Order-0 and order-1 adaptive models plus the shared
//     PACK/RLE/STRIPE/CAT transform layer. Implemented in arith.go and
//     arith_transform.go. Decode is fully supported; encode is too,
//     except the rare X_EXT (bzip2) transform — Go has no standard-
//     library bzip2 encoder, so X_EXT encode returns an error.
//   - LZMA      — a rare optional per-block codec, implemented in
//     lzma.go; the one place a third-party dependency (ulikunitz/xz)
//     is sanctioned.
//
// gzip, bzip2 and RAW blocks are handled by the CRAM container layer
// directly via the Go standard library, not here.
//
// The rANS and arith_dynamic implementations are ports of the reference
// C in samtools/htscodecs (rANS_static.c / rANS_byte.h, c_range_coder.h /
// c_simple_model.h / arith_dynamic.c); the on-wire byte layout is
// identical, so output is interoperable with htslib and the htscodecs
// test corpus serves as the compliance oracle. See docs/CRAM_ROADMAP.md
// for the testing + compliance strategy.
package codec
