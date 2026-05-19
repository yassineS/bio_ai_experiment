package sam

import (
	"bufio"
	"io"
)

// NewReader returns a Reader for r, auto-detecting whether the input is SAM
// text or BGZF-wrapped BAM by sniffing the first bytes.
//
// The detection is conservative: BGZF (gzip + BC subfield), and a raw
// "BAM\1" magic (a BAM body that has already been decompressed, e.g. by
// pkg/htsgo/iohelper which transparently strips BGZF), both route to
// the BAM reader. Plain text falls through to the line-oriented SAM reader.
// If your stream is plain gzip-compressed SAM, wrap it in compress/gzip
// first.
func NewReader(r io.Reader) (Reader, error) {
	br := bufio.NewReader(r)
	head, _ := br.Peek(16)
	if looksLikeBGZF(head) {
		return NewBAMReader(br)
	}
	if len(head) >= 4 && head[0] == 'B' && head[1] == 'A' && head[2] == 'M' && head[3] == 0x01 {
		// A raw (already-decompressed) BAM body. Skip the BGZF layer in
		// BAMReader by exposing the magic-aware constructor.
		return newBAMReaderRaw(br)
	}
	return NewSAMReader(br)
}

// looksLikeBGZF reports whether b begins with a BGZF gzip header (the same
// pattern checked by pkg/htsgo/iohelper.bgzfSniff).
func looksLikeBGZF(b []byte) bool {
	if len(b) < 16 {
		return false
	}
	if b[0] != 0x1f || b[1] != 0x8b {
		return false
	}
	if b[2] != 0x08 {
		return false
	}
	if b[3]&0x04 == 0 {
		return false
	}
	xlen := uint16(b[10]) | uint16(b[11])<<8
	if xlen < 6 {
		return false
	}
	return b[12] == 'B' && b[13] == 'C' && b[14] == 0x02 && b[15] == 0x00
}
