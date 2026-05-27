package libdeflate

import (
	"encoding/binary"
	"hash/crc32"
)

// Gzip magic numbers and field values used by libdeflate_gzip_compress.
const (
	gzipID1         = 0x1F
	gzipID2         = 0x8B
	gzipCMDeflate   = 8
	gzipMTimeNone   = 0
	gzipOSUnknown   = 0xFF
	gzipXFLFastest  = 4
	gzipXFLSlowest  = 2
	gzipHeaderLen   = 10
	gzipTrailerLen  = 8
	gzipMinOverhead = gzipHeaderLen + gzipTrailerLen
)

// writeGzipHeader appends the 10-byte gzip header used by
// libdeflate_gzip_compress to dst. FLG is always 0 (no FNAME, FCOMMENT,
// FEXTRA, FHCRC, FTEXT). MTIME is left at 0. XFL is derived from the
// compression level the same way libdeflate does it.
func writeGzipHeader(dst []byte, level int) []byte {
	var hdr [gzipHeaderLen]byte
	hdr[0] = gzipID1
	hdr[1] = gzipID2
	hdr[2] = gzipCMDeflate
	hdr[3] = 0 // FLG
	// MTIME (4 bytes, zero == unavailable).
	binary.LittleEndian.PutUint32(hdr[4:8], gzipMTimeNone)
	// XFL: libdeflate sets 4 for fastest (level<2), 2 for slowest (level>=8).
	switch {
	case level < 2:
		hdr[8] = gzipXFLFastest
	case level >= 8:
		hdr[8] = gzipXFLSlowest
	default:
		hdr[8] = 0
	}
	hdr[9] = gzipOSUnknown
	return append(dst, hdr[:]...)
}

// writeGzipTrailer appends the 8-byte gzip trailer: CRC-32 of the
// uncompressed input followed by the input length modulo 2^32, both
// little-endian.
func writeGzipTrailer(dst, src []byte) []byte {
	var trailer [gzipTrailerLen]byte
	binary.LittleEndian.PutUint32(trailer[0:4], crc32.ChecksumIEEE(src))
	binary.LittleEndian.PutUint32(trailer[4:8], uint32(len(src)))
	return append(dst, trailer[:]...)
}
