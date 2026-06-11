package fastp

import "testing"

// TestTrimPolyXUpstream_NoPanic guards the index hazards in the upstream
// poly-X walkback: an all-N tail (scan never breaks, pos==rlen) and a tail
// whose argmax base never appears literally must not panic and must return
// a no-trim result (newLen == len(seq)).
func TestTrimPolyXUpstream_NoPanic(t *testing.T) {
	cases := []struct {
		name string
		seq  string
	}{
		{"all_N", "NNNNNNNNNNNNNNNNNNNN"},
		{"N_tail", "GGGGGGGGGGNNNNNNNNNN"},
		{"empty", ""},
		{"short", "ACG"},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			defer func() {
				if r := recover(); r != nil {
					t.Fatalf("trimPolyXUpstream panicked on %q: %v", tc.seq, r)
				}
			}()
			newLen, trimmed := trimPolyXUpstream(tc.seq, 10)
			if newLen < 0 || newLen > len(tc.seq) {
				t.Fatalf("newLen %d out of range for seq len %d", newLen, len(tc.seq))
			}
			if trimmed < 0 {
				t.Fatalf("negative trimmed count %d", trimmed)
			}
		})
	}
}
