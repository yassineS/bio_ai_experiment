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
//     bits (RLE, bit-packing, 4-way striping). Planned (PR C2).
//   - LZMA      — a rare optional per-block codec. Planned (PR C-LZMA),
//     and the one place a third-party dependency (ulikunitz/xz) is
//     sanctioned. Until then an LZMA block decodes to an
//     "unsupported codec" error.
//
// gzip, bzip2 and RAW blocks are handled by the CRAM container layer
// directly via the Go standard library, not here.
//
// The rANS implementations are ports of the reference C in
// samtools/htscodecs (rANS_static.c / rANS_byte.h); the on-wire byte
// layout is identical, so output is interoperable with htslib and the
// htscodecs test corpus serves as the compliance oracle. See
// docs/CRAM_ROADMAP.md for the testing + compliance strategy.
package codec
