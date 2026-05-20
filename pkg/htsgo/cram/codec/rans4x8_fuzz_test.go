package codec

import "testing"

// FuzzRANS4x8Decode exercises the decoder against arbitrary byte
// strings. The decoder parses an attacker-controllable header and
// frequency table, so the contract under fuzz is simply: never panic,
// never hang, never allocate unboundedly — return an error or valid
// output. The maxRANSRawSize ceiling and the per-cursor bounds checks
// are what this target guards.
func FuzzRANS4x8Decode(f *testing.F) {
	// Seed with valid streams (encoder output) and a few malformed
	// shapes so the corpus starts from interesting points.
	for _, seed := range [][]byte{
		{},
		{0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0},
		mustEncode([]byte("ACGTACGTACGTNNNN"), 0),
		mustEncode([]byte("ACGTACGTACGTNNNN"), 1),
		mustEncode(bytes512(), 0),
		mustEncode(bytes512(), 1),
	} {
		f.Add(seed)
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		// The only requirement: no panic. A decode error is fine;
		// successful output is fine. We don't assert round-trip here
		// because arbitrary bytes won't be valid streams.
		out, err := RANS4x8Decode(data)
		if err == nil && len(out) > maxRANSRawSize {
			t.Fatalf("decode returned %d bytes, above the safety ceiling", len(out))
		}
	})
}

func mustEncode(in []byte, order int) []byte {
	out, err := RANS4x8Encode(in, order)
	if err != nil {
		panic(err)
	}
	return out
}

func bytes512() []byte {
	out := make([]byte, 512)
	for i := range out {
		out[i] = byte(i * 7 % 17)
	}
	return out
}
