package libdeflate

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// splitDir holds the genuine-libdeflate block-split oracle fixtures: the
// 65280-byte payload that libdeflate splits into two DEFLATE blocks, and
// the gzip-wrapped reference outputs at levels 6 and 7.
const splitDir = "testdata/split"

// rawDeflateFromGzip extracts the raw DEFLATE stream from a gzip member
// produced by libdeflate-gzip: a fixed 10-byte RFC 1952 header (no extra
// fields are emitted by libdeflate-gzip on these fixtures) and an 8-byte
// CRC32 + ISIZE trailer.
func rawDeflateFromGzip(t *testing.T, gz []byte) []byte {
	t.Helper()
	if len(gz) < 18 {
		t.Fatalf("gzip member too short: %d bytes", len(gz))
	}
	return gz[10 : len(gz)-8]
}

// TestBlockSplitOracle verifies that DeflateCompress reproduces genuine
// libdeflate's two-block split decision byte-for-byte. Before the fix to
// minLensTable, our parser used min_len=3 (rather than 4) on this VCF
// payload, producing extra short matches that shifted the split point
// from offset 8798 to 11116 and changed every downstream byte.
func TestBlockSplitOracle(t *testing.T) {
	in, err := os.ReadFile(filepath.Join(splitDir, "block_splits.in"))
	if err != nil {
		t.Fatalf("read input: %v", err)
	}

	for _, level := range []int{6, 7} {
		gz, err := os.ReadFile(filepath.Join(splitDir,
			fileForLevel(level)))
		if err != nil {
			t.Fatalf("read level %d reference: %v", level, err)
		}
		want := rawDeflateFromGzip(t, gz)

		got, err := DeflateCompress(in, level)
		if err != nil {
			t.Fatalf("DeflateCompress level %d: %v", level, err)
		}
		if !bytes.Equal(got, want) {
			t.Fatalf("level %d byte mismatch: len got=%d want=%d, first diff at %d",
				level, len(got), len(want), firstDiff(got, want))
		}
	}
}

func fileForLevel(level int) string {
	switch level {
	case 6:
		return "block_splits.l6.gz"
	case 7:
		return "block_splits.l7.gz"
	default:
		return ""
	}
}
