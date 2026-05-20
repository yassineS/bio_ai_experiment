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
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		mustEncode16([]byte("ACGTACGTACGTNNNN")),
		mustEncode16(bytes512()),
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

func mustEncode16(in []byte) []byte {
	out, err := RANS4x16Encode(in, 0)
	if err != nil {
		panic(err)
	}
	return out
}
