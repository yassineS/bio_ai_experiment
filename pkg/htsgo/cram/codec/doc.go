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
//     arith_transform.go. Both decode and encode are supported,
//     including the rare X_EXT (bzip2) transform: Go's standard library
//     has no bzip2 encoder, so an in-tree one lives in bzip2_encode.go.
//   - bzip2 — a pure-Go in-tree bzip2 *encoder* (bzip2_encode.go),
//     the counterpart to compress/bzip2's decoder. It backs the X_EXT
//     external codec above and CRAM bzip2 blocks (compression method 2);
//     its output is byte-compatible with compress/bzip2, the system
//     bzip2 and htslib/libbz2.
//   - fqzcomp   — CRAM v3.1's quality-score codec (compression
//     method 7). A parameterised context model driving the same
//     adaptive range coder as arith_dynamic. Implemented in
//     fqzcomp.go (decode + parameter machinery) and
//     fqzcomp_encode.go (the strategy-driven encoder). Both decode
//     and encode are byte-exact against the htscodecs corpus.
//   - name tokeniser — CRAM v3.1's read-name codec (compression
//     method 8). It splits each name into typed tokens, models each
//     name as a diff against an earlier one, and routes the per-token
//     streams through the rANS / arith sub-codecs above. Implemented
//     in nametok.go (decode) and nametok_encode.go (encode). Decode is
//     byte-exact against the full htscodecs tok3 vector set; the
//     encoder round-trips at every level (see NameTokEncode for the
//     byte-exactness caveat).
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
