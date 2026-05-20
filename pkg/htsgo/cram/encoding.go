package cram

import (
	"fmt"
)

// EncodingID identifies one of the CRAM data-series encodings (the
// "codec" of a data series, distinct from a block's compression
// method). Every data series in a CRAM compression header is described
// by an encoding ID followed by that encoding's parameter bytes.
type EncodingID int32

// The CRAM v3.0 data-series encoding identifiers.
const (
	EncodingNull          EncodingID = 0 // No data.
	EncodingExternal      EncodingID = 1 // Values verbatim from an external block.
	EncodingGolomb        EncodingID = 2 // Golomb code (CORE bitstream).
	EncodingHuffman       EncodingID = 3 // Canonical Huffman code (CORE bitstream).
	EncodingByteArrayLen  EncodingID = 4 // Length sub-encoding + values sub-encoding.
	EncodingByteArrayStop EncodingID = 5 // External bytes up to a stop byte.
	EncodingBeta          EncodingID = 6 // Fixed nbits binary with an offset.
	EncodingSubexp        EncodingID = 7 // Sub-exponential code with an offset.
	EncodingGolombRice    EncodingID = 8 // Golomb-Rice code (CORE bitstream).
	EncodingGamma         EncodingID = 9 // Elias gamma code with an offset.
)

// String returns the spec name of the encoding identifier.
func (e EncodingID) String() string {
	switch e {
	case EncodingNull:
		return "NULL"
	case EncodingExternal:
		return "EXTERNAL"
	case EncodingGolomb:
		return "GOLOMB"
	case EncodingHuffman:
		return "HUFFMAN"
	case EncodingByteArrayLen:
		return "BYTE_ARRAY_LEN"
	case EncodingByteArrayStop:
		return "BYTE_ARRAY_STOP"
	case EncodingBeta:
		return "BETA"
	case EncodingSubexp:
		return "SUBEXP"
	case EncodingGolombRice:
		return "GOLOMB_RICE"
	case EncodingGamma:
		return "GAMMA"
	default:
		return fmt.Sprintf("unknown(%d)", int32(e))
	}
}

// Encoding is a parsed CRAM data-series encoding: its identifier plus
// the parameters that govern it. Only the fields relevant to the
// encoding's ID are populated; the others are zero. A nil *Encoding,
// or one whose ID is EncodingNull, decodes to no values.
type Encoding struct {
	// ID is the encoding identifier.
	ID EncodingID

	// ExternalID is the content-id of the external block the values
	// come from. It applies to EXTERNAL and BYTE_ARRAY_STOP.
	ExternalID int32

	// StopByte is the delimiter byte that ends each value of a
	// BYTE_ARRAY_STOP series.
	StopByte byte

	// LenEnc and ValEnc are the length and values sub-encodings of a
	// BYTE_ARRAY_LEN encoding.
	LenEnc *Encoding
	ValEnc *Encoding

	// Offset is the constant added to every decoded value of a BETA,
	// SUBEXP, GAMMA, GOLOMB or GOLOMB_RICE encoding.
	Offset int32

	// NumBits is the fixed code width of a BETA encoding.
	NumBits int32

	// K is the order parameter of a SUBEXP or GOLOMB_RICE encoding.
	K int32

	// M is the divisor parameter of a GOLOMB encoding.
	M int32

	// Symbols is the alphabet of a HUFFMAN encoding (one entry per
	// code word).
	Symbols []int32

	// BitLengths is the canonical-code bit length of each Huffman
	// symbol, parallel to Symbols.
	BitLengths []int32

	// huffman is the lazily-built canonical Huffman decoder; nil until
	// the first decode through this encoding.
	huffman *huffmanTable
}

// parseEncoding reads one CRAM encoding from p starting at off: an
// ITF-8 codec id, an ITF-8 parameter-byte count, and that many bytes of
// codec-specific parameters. It returns the parsed encoding and the
// offset just past the parameter bytes.
func parseEncoding(p []byte, off int) (*Encoding, int, error) {
	id, n, err := itf8At(p, off)
	if err != nil {
		return nil, off, fmt.Errorf("cram: encoding id: %w", err)
	}
	off += n
	plen, n, err := itf8At(p, off)
	if err != nil {
		return nil, off, fmt.Errorf("cram: encoding parameter length: %w", err)
	}
	off += n
	if plen < 0 {
		return nil, off, fmt.Errorf("cram: encoding declares negative parameter length %d", plen)
	}
	end := off + int(plen)
	if end > len(p) || end < off {
		return nil, off, fmt.Errorf("cram: encoding parameter block (%d bytes) overruns the compression header", plen)
	}
	enc := &Encoding{ID: EncodingID(id)}
	// The parameters are parsed from the sub-slice p[off:end]; a
	// well-formed encoding consumes exactly that region. Sub-encodings
	// (BYTE_ARRAY_LEN) recurse into the same byte slice.
	body := p[:end]
	cur := off
	switch enc.ID {
	case EncodingNull:
		// No parameters.
	case EncodingExternal:
		enc.ExternalID, cur, err = readITF8Param(body, cur, "EXTERNAL content id")
	case EncodingByteArrayStop:
		if cur >= end {
			return nil, off, fmt.Errorf("cram: BYTE_ARRAY_STOP encoding missing stop byte")
		}
		enc.StopByte = body[cur]
		cur++
		enc.ExternalID, cur, err = readITF8Param(body, cur, "BYTE_ARRAY_STOP content id")
	case EncodingByteArrayLen:
		enc.LenEnc, cur, err = parseEncoding(body, cur)
		if err == nil {
			enc.ValEnc, cur, err = parseEncoding(body, cur)
		}
	case EncodingBeta:
		enc.Offset, cur, err = readITF8Param(body, cur, "BETA offset")
		if err == nil {
			enc.NumBits, cur, err = readITF8Param(body, cur, "BETA num bits")
		}
		if err == nil && (enc.NumBits < 0 || enc.NumBits > 32) {
			return nil, off, fmt.Errorf("cram: BETA encoding declares %d bits (must be 0..32)", enc.NumBits)
		}
	case EncodingSubexp:
		enc.Offset, cur, err = readITF8Param(body, cur, "SUBEXP offset")
		if err == nil {
			enc.K, cur, err = readITF8Param(body, cur, "SUBEXP k")
		}
		if err == nil && enc.K < 0 {
			return nil, off, fmt.Errorf("cram: SUBEXP encoding declares negative k %d", enc.K)
		}
	case EncodingGamma:
		enc.Offset, cur, err = readITF8Param(body, cur, "GAMMA offset")
	case EncodingGolomb:
		enc.Offset, cur, err = readITF8Param(body, cur, "GOLOMB offset")
		if err == nil {
			enc.M, cur, err = readITF8Param(body, cur, "GOLOMB M")
		}
		if err == nil && enc.M <= 0 {
			return nil, off, fmt.Errorf("cram: GOLOMB encoding declares non-positive M %d", enc.M)
		}
	case EncodingGolombRice:
		enc.Offset, cur, err = readITF8Param(body, cur, "GOLOMB_RICE offset")
		if err == nil {
			enc.K, cur, err = readITF8Param(body, cur, "GOLOMB_RICE log2m")
		}
		if err == nil && (enc.K < 0 || enc.K > 30) {
			return nil, off, fmt.Errorf("cram: GOLOMB_RICE encoding declares out-of-range log2m %d", enc.K)
		}
	case EncodingHuffman:
		var nsym int32
		nsym, cur, err = readITF8Param(body, cur, "HUFFMAN alphabet size")
		if err == nil {
			if nsym < 0 || cur+int(nsym) > end {
				return nil, off, fmt.Errorf("cram: HUFFMAN alphabet size %d overruns parameters", nsym)
			}
			enc.Symbols = make([]int32, nsym)
			for i := int32(0); i < nsym && err == nil; i++ {
				enc.Symbols[i], cur, err = readITF8Param(body, cur, "HUFFMAN symbol")
			}
		}
		if err == nil {
			var nlen int32
			nlen, cur, err = readITF8Param(body, cur, "HUFFMAN bit-length count")
			if err == nil {
				if nlen != nsym {
					return nil, off, fmt.Errorf("cram: HUFFMAN has %d symbols but %d bit lengths", nsym, nlen)
				}
				enc.BitLengths = make([]int32, nlen)
				for i := int32(0); i < nlen && err == nil; i++ {
					enc.BitLengths[i], cur, err = readITF8Param(body, cur, "HUFFMAN bit length")
				}
			}
		}
	default:
		return nil, off, fmt.Errorf("cram: unknown data-series encoding id %d", id)
	}
	if err != nil {
		return nil, off, err
	}
	// The parser advances cur within p[:end]; trailing parameter bytes
	// are tolerated (some writers pad) but cur must not exceed end.
	if cur > end {
		return nil, off, fmt.Errorf("cram: %s encoding parameters overran their %d-byte block", enc.ID, plen)
	}
	return enc, end, nil
}

// readITF8Param reads an ITF-8 integer at off within p, bounded by
// len(p), wrapping any error with what for context.
func readITF8Param(p []byte, off int, what string) (int32, int, error) {
	v, n, err := itf8At(p, off)
	if err != nil {
		return 0, off, fmt.Errorf("cram: %s: %w", what, err)
	}
	return v, off + n, nil
}
