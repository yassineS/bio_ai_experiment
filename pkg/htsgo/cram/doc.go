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
// the decoded stream as text SAM via WriteSAM.
//
// Writing (C8): RecordWriter encodes a sam.Header plus a stream of
// sam.Records as a CRAM v3.0 file — the inverse of RecordReader.
// Construct one with NewRecordWriter (or CreateCRAM for a file path),
// feed records with Write, and finalise with Close; WriteCRAM is a
// one-shot convenience. The writer is deliberately simple: it produces
// reference-free CRAM (a mapped read's bases are stored literally as a
// base-stretch read feature rather than diffed against a reference), it
// encodes every data series with the EXTERNAL codec (integers as ITF-8,
// byte runs raw, byte arrays via BYTE_ARRAY_LEN) and never uses the CORE
// bitstream, and it gzip-compresses each data block when that shrinks
// it. Every record is encoded "detached" so it carries its own mate
// fields, which makes RNEXT/PNEXT/TLEN round-trip without the
// downstream-mate optimisation. One slice is written per container,
// capped at a fixed record count. The output is valid CRAM — it
// re-reads through RecordReader and through samtools — but is not
// optimised for compression ratio and is not byte-identical to
// samtools' own encoder. A record shape the simple writer cannot encode
// (an unknown reference, a CIGAR/SEQ length mismatch, an unsupported
// CIGAR operation) is rejected by Write with a clear error. The CIGAR
// match operators =, X and M all encode as a copy of the read bases, so
// a = or X run re-reads as M — expected CRAM lossiness, the same kind
// as an unmapped read losing its MAPQ and CIGAR.
//
// Reference resolution (C5) completes the reference-backed decode path.
// A reference-free CRAM file is fully recoverable on its own; a
// reference-backed file stores each mapped read's bases as a copy of a
// reference span plus the read features, so it needs the reference
// supplied. Call RecordReader.SetReference (or SetReferenceFASTA) with
// an indexed FASTA, or RecordReader.SetRefCache / UseRefCacheFromEnv to
// use the htslib REF_CACHE local cache, before the first Read. Each
// slice's reference span is fetched and its MD5 verified against the
// slice header — an MD5 mismatch is always a hard error, so the decoder
// never produces sequence against the wrong reference. Without a
// reference the bases an external reference would supply are filled with
// 'N' and NeedsReference reports it. The network REF_PATH URL-fetch
// mechanism is also supported (opt-in via the REF_PATH environment
// variable; see RefPathReference / UseRefPathFromEnv): a reference
// addressed only by MD5 is fetched from a URL endpoint — the EBI ENA
// registry by default — when every local source misses. An unresolvable
// reference is a clear error naming the missing MD5.
//
// CRAM index (.crai): ReadCRAI / OpenCRAI parse the gzip-compressed TSV
// index that accompanies a CRAM file. CRAIIndex.Query (and the
// region-typed QueryRegion) return the index entries whose alignment
// range overlaps a wanted region, giving a caller the container and
// slice byte offsets to seek to.
//
// CRAM v4.0 decode: the v4.0 draft is supported for decode (encode stays
// v3.0). v4 keeps the v3 container/block/slice tree but changes the
// variable-length integers and the data-series codec set. Every integer —
// container header, block header, slice header, the compression-header
// maps and the per-encoding parameters — is a big-endian uint7 LEB128
// varint (htscodecs varint.h, cram_io.c uint7_get_*) rather than ITF-8 /
// LTF-8, with signed fields zig-zag encoded; alignment coordinates widen
// to 64-bit. The integer data series move to the VARINT_UNSIGNED,
// VARINT_SIGNED and CONST_INT/CONST_BYTE codecs (EXTERNAL is restricted to
// byte data), the read group appears in the tag dictionary as an "RG*"
// placeholder, paired read names are deduplicated (the upstream mate's
// name is copied from the downstream mate), and the QO preservation flag
// drives reverse-strand quality reversal. The end-of-file marker is a
// distinct 31-byte sentinel. The transform codecs XPACK (bit-packing),
// XRLE (run-length) and XDELTA (delta) now decode: each wraps one or two
// sub-codecs (dispatched recursively through the same codec table) and
// reverses the transform — XPACK expands packed bytes through its reverse
// map (htscodecs hts_unpack, least-significant-symbol-first), XRLE expands
// runs from a literal and a length stream (hts_rle_decode), and XDELTA
// prefix-sums zig-zag deltas (value-at-a-time for the 32-bit integer path,
// 2-byte words for the block path). Version dispatch is threaded through
// FileDefinition.Major via a small
// intReader so the v2/v3 path is byte-for-byte unchanged.
//
// Integer encodings: CRAM v2/v3 use two self-delimiting integer
// encodings, ITF-8 (a 1-5 byte 32-bit value) and LTF-8 (a 1-9 byte 64-bit
// value), both read by counting the leading 1-bits of the first byte;
// CRAM v4 uses the uint7 varint described above. They are all distinct
// from the varints used inside the rANS codec sub-package.
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
//   - CRAM v4.0 draft (htslib cram_io.c / cram_decode.c major>=4 paths).
//   - htscodecs varint.h (the uint7 varint encoding v4 adopts).
//   - htslib cram/cram_io.c, cram/cram_structs.h.
package cram
