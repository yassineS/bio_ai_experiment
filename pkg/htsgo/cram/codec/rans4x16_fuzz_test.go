package codec

import "testing"

// FuzzRANS4x16Decode exercises the decoder against arbitrary byte
// strings. The decoder parses an attacker-controllable format byte,
// raw-size varint and frequency table, so the contract under fuzz is
// simply: never panic, never hang, never allocate unboundedly — return
// an error or valid output. The maxRANSRawSize ceiling and the
// per-cursor bounds checks are what this target guards.
func FuzzRANS4x16Decode(f *testing.F) {
	for _, seed := range [][]byte{
		{},
		{0},
		{0x00, 0x00},
		{0x20, 0x00},
		{0x01, 0x00},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		mustEncode16([]byte("ACGTACGTACGTNNNN"), 0),
		mustEncode16([]byte("ACGTACGTACGTNNNN"), 1),
		mustEncode16(bytes512(), 0),
		mustEncode16(bytes512(), 1),
		// Transform-layer seeds: PACK, RLE, PACK+RLE and STRIPE, with
		// and without order-1, exercising the C2.2 decode paths.
		mustEncode16([]byte("ACGTACGTACGTNNNN"), 0x80),
		mustEncode16([]byte("ACGTACGTACGTNNNN"), 0x81),
		mustEncode16([]byte("AAAACCCCGGGGTTTT"), 0x40),
		mustEncode16([]byte("AAAACCCCGGGGTTTT"), 0x41),
		mustEncode16(bytes512(), 0xC0),
		mustEncode16(bytes512(), 0xC1),
		mustEncode16(bytes512(), 0x08),
		mustEncode16(bytes512(), 0x09),
		// X_32 (32-way) seeds: plain order-0/1 and combined with the
		// PACK/RLE/STRIPE transforms, exercising the C2.3 decode paths.
		mustEncode16([]byte("ACGTACGTACGTNNNN"), 0x04),
		mustEncode16([]byte("ACGTACGTACGTNNNN"), 0x05),
		mustEncode16(bytes512(), 0x04),
		mustEncode16(bytes512(), 0x05),
		mustEncode16(bytes512(), 0x84),
		mustEncode16(bytes512(), 0x85),
		mustEncode16(bytes512(), 0x44),
		mustEncode16(bytes512(), 0x45),
		mustEncode16(bytes512(), 0x0C),
		mustEncode16(bytes512(), 0x0D),
	} {
		f.Add(seed)
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		out, err := RANS4x16Decode(data)
		if err == nil && len(out) > maxRANSRawSize {
			t.Fatalf("decode returned %d bytes, above the safety ceiling", len(out))
		}
	})
}

func mustEncode16(in []byte, order int) []byte {
	out, err := RANS4x16Encode(in, order)
	if err != nil {
		panic(err)
	}
	return out
}
