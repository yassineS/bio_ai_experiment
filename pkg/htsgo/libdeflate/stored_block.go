package libdeflate

import "encoding/binary"

// maxStoredBlockLen is the largest payload that fits in a single DEFLATE
// stored block (RFC 1951 §3.2.4: LEN is a 16-bit field).
const maxStoredBlockLen = 0xFFFF

// writeStoredBlocks emits one or more STORED blocks covering data. The
// final block has BFINAL=1 iff last is true; intermediate blocks
// always have BFINAL=0. Each block: 3-bit header (BFINAL | BTYPE=00),
// then padding to the next byte boundary, then LEN (2 bytes LE),
// ~LEN (2 bytes LE), then the raw bytes.
//
// libdeflate's deflate_compress_none takes this same path for inputs
// short enough to use the level-0/passthrough route, so we match its
// behavior of emitting an empty stored block for zero-length input.
func writeStoredBlocks(bw *bitWriter, data []byte, last bool) {
	if len(data) == 0 {
		writeStoredBlockChunk(bw, nil, last)
		return
	}
	for len(data) > 0 {
		n := len(data)
		if n > maxStoredBlockLen {
			n = maxStoredBlockLen
		}
		isLast := last && n == len(data)
		writeStoredBlockChunk(bw, data[:n], isLast)
		data = data[n:]
	}
}

// writeStoredBlockChunk emits a single STORED block of len(chunk) bytes.
func writeStoredBlockChunk(bw *bitWriter, chunk []byte, last bool) {
	var bfinal uint64
	if last {
		bfinal = 1
	}
	bw.writeBits(bfinal|(uint64(blockTypeStored)<<1), 3)
	bw.alignToByte()
	var hdr [4]byte
	binary.LittleEndian.PutUint16(hdr[0:2], uint16(len(chunk)))
	binary.LittleEndian.PutUint16(hdr[2:4], ^uint16(len(chunk)))
	bw.writeBytes(hdr[:])
	if len(chunk) > 0 {
		bw.writeBytes(chunk)
	}
}
