package codec

import "testing"

// FuzzArithDecode exercises the arith_dynamic decoder against arbitrary
// byte strings. The decoder parses an attacker-controllable format byte,
// raw-size varint, PACK/RLE meta and adaptive range-coded payload, so the
// contract under fuzz is simply: never panic, never hang, never allocate
// unboundedly — return an error or valid output. The maxRANSRawSize
// ceiling, the maxStripeDepth cap and the per-cursor bounds checks are
// what this target guards.
func FuzzArithDecode(f *testing.F) {
	for _, seed := range [][]byte{
		{},
		{0},
		{0x00, 0x00},
		{0x20, 0x00},
		{0x01, 0x00},
		{0x10, 0x00},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		mustArithEncode([]byte("ACGTACGTACGTNNNN"), 0),
		mustArithEncode([]byte("ACGTACGTACGTNNNN"), 1),
		mustArithEncode(bytes512(), 0),
		mustArithEncode(bytes512(), 1),
		// Transform-layer seeds: PACK, RLE, PACK+RLE and STRIPE, with and
		// without order-1.
		mustArithEncode([]byte("ACGTACGTACGTNNNN"), 0x80),
		mustArithEncode([]byte("ACGTACGTACGTNNNN"), 0x81),
		mustArithEncode([]byte("AAAACCCCGGGGTTTT"), 0x40),
		mustArithEncode([]byte("AAAACCCCGGGGTTTT"), 0x41),
		mustArithEncode(bytes512(), 0xC0),
		mustArithEncode(bytes512(), 0xC1),
		mustArithEncode(bytes512(), 0x08),
		mustArithEncode(bytes512(), 0x09),
	} {
		f.Add(seed)
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		out, err := ArithDecode(data)
		if err == nil && len(out) > maxRANSRawSize {
			t.Fatalf("decode returned %d bytes, above the safety ceiling", len(out))
		}
	})
}

func mustArithEncode(in []byte, order int) []byte {
	out, err := ArithEncode(in, order)
	if err != nil {
		panic(err)
	}
	return out
}
