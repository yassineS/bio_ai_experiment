package bgzf

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// TestBGZFOracleUpstreamLibdeflate is the definitive byte-for-byte parity
// gate for the BGZF writer against genuine upstream htslib `bgzip`.
//
// The golden files under testdata/oracle/*.bgzf.golden were produced by a
// libdeflate-linked build of upstream htslib 1.23.1 (`bgzip -c -l 6`),
// the same configuration shipped by the htslib release binaries and the
// conda/debian packages. Each *.in is compressed by our pure-Go BGZF
// writer (which routes block payloads through the in-tree libdeflate port
// in pkg/htsgo/libdeflate) and the result must equal the upstream golden
// byte-for-byte.
//
// This locks in 1:1 parity for the whole BGZF layer — block framing, the
// BC subfield, the CRC32/ISIZE trailers, the multi-block split boundaries
// and the DEFLATE payload bytes themselves. The goldens are committed so
// the gate runs without the C toolchain; regenerate them with a
// libdeflate-linked `bgzip -c -l 6` if the fixtures ever change.
func TestBGZFOracleUpstreamLibdeflate(t *testing.T) {
	ins, err := filepath.Glob("testdata/oracle/*.in")
	if err != nil {
		t.Fatal(err)
	}
	if len(ins) == 0 {
		t.Skip("no oracle fixtures present")
	}
	for _, in := range ins {
		name := filepath.Base(in)
		t.Run(name, func(t *testing.T) {
			src, err := os.ReadFile(in)
			if err != nil {
				t.Fatal(err)
			}
			golden := in[:len(in)-len(".in")] + ".bgzf.golden"
			want, err := os.ReadFile(golden)
			if err != nil {
				t.Skipf("golden missing: %v", err)
			}
			var buf bytes.Buffer
			w := NewWriter(&buf) // default level 6, matching upstream bgzip -l 6
			if _, err := w.Write(src); err != nil {
				t.Fatal(err)
			}
			if err := w.Close(); err != nil {
				t.Fatal(err)
			}
			got := buf.Bytes()
			if !bytes.Equal(got, want) {
				n := len(got)
				if len(want) < n {
					n = len(want)
				}
				diff := -1
				for i := 0; i < n; i++ {
					if got[i] != want[i] {
						diff = i
						break
					}
				}
				t.Fatalf("BGZF output differs from upstream libdeflate bgzip: got %d bytes, want %d bytes, first diff at byte %d", len(got), len(want), diff)
			}
		})
	}
}
