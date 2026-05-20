// Package cram implements a structural parser for the CRAM file format
// used by the htslib/samtools ecosystem.
//
// CRAM is a tree-shaped binary container: a file definition, followed by
// a sequence of containers, each holding a number of blocks. The first
// container carries the SAM text header; every later container holds a
// compression-header block followed by one or more slices (a slice being
// a slice-header block plus its data blocks).
//
// This package covers the structural layer (the C3 milestone of
// docs/CRAM_ROADMAP.md), the data-series encoding layer (C4a) and the
// record-reconstruction layer (C4b). It parses the file/container/block
// tree, validates the per-structure CRC32 checksums that CRAM v3+
// embeds, and decompresses every block whose compression method is
// supported. On top of that it parses a container's compression header
// (the preservation, data-series and tag-encoding maps) and every slice
// header, and decodes a data series through its encoding: the CRAM
// encoding zoo — NULL, EXTERNAL, GOLOMB, HUFFMAN, BYTE_ARRAY_LEN,
// BYTE_ARRAY_STOP, BETA, SUBEXP, GOLOMB_RICE and GAMMA — reading from
// the slice's CORE bitstream and external blocks. Use ParseDataContainer
// to obtain a DataContainer and its Slices, then Slice.DrainSeries /
// Slice.DrainTag (or the fixed-count Decode*Series methods) for raw
// series access.
//
// Record reconstruction (C4b) sits on top: RecordReader walks the whole
// stream, decodes each slice's data series through the per-record
// traversal — bit flags, CRAM flags, reference id, read length and
// position, read group, read name, mate information, tags, and the
// read-feature list (mapped) or raw bases (unmapped) — and yields
// reconstructed *sam.Record values. It resolves the read-feature codes
// into SEQ, QUAL and CIGAR, links downstream mate pairs, and can emit
// the decoded stream as text SAM via WriteSAM. A reference-free CRAM
// file is fully recoverable; a reference-backed file fills the bases an
// external reference would supply with 'N' and reports NeedsReference.
//
// Integer encodings: CRAM uses two self-delimiting integer encodings,
// ITF-8 (a 1-5 byte 32-bit value) and LTF-8 (a 1-9 byte 64-bit value).
// Both are read by counting the leading 1-bits of the first byte. They
// are distinct from the varints used inside the rANS codec sub-package.
//
// Block compression: method 0 is raw, 1 is gzip, 2 is bzip2, 4 is
// rANS 4x8 and 5 is rANS 4x16 — all supported. Methods 3 (LZMA), 6
// (range/arithmetic), 7 (fqzcomp) and 8 (name tokeniser) are out of
// scope for C3 and decompressing such a block returns a clear
// "unsupported compression method" error.
//
// CRC32 validation: CRAM v3.0 and v3.1 append a 4-byte little-endian
// CRC32 (the standard IEEE polynomial, hash/crc32) to every container
// header and every block, covering all the bytes of that structure up
// to the CRC field. This parser validates every checksum as it walks
// the tree; a mismatch is reported as an error because it means a
// structure was mis-delineated.
//
// References:
//   - CRAM format specification v3.0 and v3.1 (hsformats.github.io).
//   - htslib cram/cram_io.c, cram/cram_structs.h.
package cram
