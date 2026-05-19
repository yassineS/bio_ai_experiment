// Package tabix implements htslib's generic position index for bgzipped
// tab-delimited files (VCF, BED, GFF, SAM, and custom layouts).
//
// A `.tbi` file records, for each reference sequence, two structures:
//
//   - A bin index (UCSC binning scheme) mapping bin numbers to a list of
//     "chunks" — pairs of 64-bit virtual offsets that bracket the byte range
//     in the bgzipped data file where records belonging to that bin live.
//   - A linear index — one virtual offset per 16,384 bp tile, giving the
//     earliest record that *might* overlap the tile. Used to trim chunks
//     during a region query.
//
// A virtual offset packs (compressed-file offset, uncompressed-byte offset)
// into a single uint64: `(coffset << 16) | uoffset`. With BGZF blocks capped
// at 64 KiB the low 16 bits suffice for the in-block byte offset.
//
// The on-disk format is described in the SAM specification ("The Tabix index
// format") and in Heng Li's "Tabix: fast retrieval of sequence features from
// generic TAB-delimited files" (Bioinformatics 27, 2011).
//
// The Index type supports three core operations:
//
//   - Build: stream a bgzipped data file once, record the bin/linear-index
//     contributions for every record, and return an in-memory Index.
//   - Write / Read: serialise / deserialise the bgzipped `.tbi` byte
//     representation.
//   - Query: given (chrom, beg, end), enumerate the byte ranges that may
//     contain matching records and yield each matching raw line.
//
// The implementation targets bytewise compatibility with htslib's `.tbi`
// files so that an index built by this package can be consumed by upstream
// tabix / bcftools / samtools and vice versa.
package tabix
