package codec

import "testing"

// FuzzLZMADecode exercises the LZMA block decoder against arbitrary byte
// strings. The input is an attacker-controllable .xz container stream,
// so the contract under fuzz is simply: never panic, never hang, never
// allocate unboundedly — return an error or valid output. The
// maxLZMARawSize ceiling is what this target guards.
func FuzzLZMADecode(f *testing.F) {
	// Seed with valid streams (encoder output) and a few malformed
	// shapes so the corpus starts from interesting points.
	for _, seed := range [][]byte{
		{},
		{0},
		{0xFD, '7', 'z', 'X', 'Z', 0x00},
		mustLZMAEncode(nil),
		mustLZMAEncode([]byte("ACGTACGTACGTNNNN")),
		mustLZMAEncode(bytes512()),
	} {
		f.Add(seed)
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		// The only requirement: no panic. A decode error is fine;
		// successful output is fine. We don't assert round-trip here
		// because arbitrary bytes won't be valid streams.
		out, err := LZMADecode(data)
		if err == nil && len(out) > maxLZMARawSize {
			t.Fatalf("decode returned %d bytes, above the safety ceiling", len(out))
		}
	})
}

func mustLZMAEncode(in []byte) []byte {
	out, err := LZMAEncode(in)
	if err != nil {
		panic(err)
	}
	return out
}
