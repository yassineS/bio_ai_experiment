package libdeflate

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// oracleDir is the on-disk location of the Slice 0 oracle fixtures.
const oracleDir = "testdata/oracle"

// readOracle returns the raw input and the libdeflate-reference gzip
// output for the named fixture.
func readOracle(t *testing.T, name string) (src, gz []byte) {
	t.Helper()
	src, err := os.ReadFile(filepath.Join(oracleDir, name+".bin"))
	if err != nil {
		t.Fatalf("read input: %v", err)
	}
	gz, err = os.ReadFile(filepath.Join(oracleDir, name+".gz"))
	if err != nil {
		t.Fatalf("read reference gz: %v", err)
	}
	return src, gz
}

func TestGzipCompress_OracleEmpty(t *testing.T) {
	src, want := readOracle(t, "empty")
	got, err := GzipCompress(src, 6)
	if err != nil {
		t.Fatalf("GzipCompress: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("byte mismatch:\ngot:  % x\nwant: % x", got, want)
	}
}

func TestGzipCompress_OracleSingleByte(t *testing.T) {
	src, want := readOracle(t, "single_byte")
	got, err := GzipCompress(src, 6)
	if err != nil {
		t.Fatalf("GzipCompress: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("byte mismatch:\ngot:  % x\nwant: % x", got, want)
	}
}

func TestGzipCompress_OracleRepeatedA(t *testing.T) {
	src, want := readOracle(t, "repeated_a")
	got, err := GzipCompress(src, 6)
	if err != nil {
		t.Fatalf("GzipCompress: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("byte mismatch:\ngot:  % x\nwant: % x", got, want)
	}
}

func TestGzipCompress_OracleRandom64K(t *testing.T) {
	src, want := readOracle(t, "random_64k")
	got, err := GzipCompress(src, 6)
	if err != nil {
		t.Fatalf("GzipCompress: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("byte mismatch: len got=%d want=%d, first diff at %d",
			len(got), len(want), firstDiff(got, want))
	}
}

func TestGzipCompress_OracleBGZFPayload(t *testing.T) {
	src, want := readOracle(t, "bgzf_payload")
	got, err := GzipCompress(src, 6)
	if err != nil {
		t.Fatalf("GzipCompress: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("byte mismatch: len got=%d want=%d, first diff at %d",
			len(got), len(want), firstDiff(got, want))
	}
}

// firstDiff returns the index of the first byte that differs between
// a and b, or -1 if they share the entire shorter prefix.
func firstDiff(a, b []byte) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return i
		}
	}
	if len(a) != len(b) {
		return n
	}
	return -1
}

func TestGzipDecompress_RoundTrip(t *testing.T) {
	for _, name := range []string{"empty", "single_byte", "repeated_a", "random_64k", "bgzf_payload"} {
		t.Run(name, func(t *testing.T) {
			src, _ := readOracle(t, name)
			gz, err := GzipCompress(src, 6)
			if err != nil {
				t.Fatalf("GzipCompress: %v", err)
			}
			gr, err := gzip.NewReader(bytes.NewReader(gz))
			if err != nil {
				t.Fatalf("gzip.NewReader: %v", err)
			}
			defer gr.Close()
			got, err := io.ReadAll(gr)
			if err != nil {
				t.Fatalf("read decompressed: %v", err)
			}
			if !bytes.Equal(got, src) {
				t.Fatalf("round trip mismatch:\ngot:  % x\nwant: % x", got, src)
			}
		})
	}
}

func TestGzipCompress_InvalidLevel(t *testing.T) {
	for _, lvl := range []int{-1, 0, 13, 99} {
		if _, err := GzipCompress(nil, lvl); err == nil {
			t.Errorf("level %d: want error, got nil", lvl)
		}
	}
}
