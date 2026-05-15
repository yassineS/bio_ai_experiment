package fasta

import (
	"bytes"
	"encoding/binary"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

// makeBGZF returns the BGZF-compressed bytes of payload.
func makeBGZF(t *testing.T, payload []byte) []byte {
	t.Helper()
	var buf bytes.Buffer
	w := bgzip.NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("bgzip write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	return buf.Bytes()
}

// writeTempFile writes data to a fresh temp file under t.TempDir() and
// returns the path. The file is cleaned up automatically when the test
// exits.
func writeTempFile(t *testing.T, dir, name string, data []byte) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, data, 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

func TestIsBGZFFile(t *testing.T) {
	dir := t.TempDir()

	bgzfData := makeBGZF(t, []byte(">chr1\nACGT\n"))
	bgzfPath := writeTempFile(t, dir, "ref.fa.gz", bgzfData)
	got, err := isBGZFFile(bgzfPath)
	if err != nil {
		t.Fatalf("isBGZFFile(bgzf): %v", err)
	}
	if !got {
		t.Errorf("isBGZFFile on BGZF file: got false, want true")
	}

	plainPath := writeTempFile(t, dir, "plain.fa", []byte(">chr1\nACGT\n"))
	got, err = isBGZFFile(plainPath)
	if err != nil {
		t.Fatalf("isBGZFFile(plain): %v", err)
	}
	if got {
		t.Errorf("isBGZFFile on plain text: got true, want false")
	}

	// Tiny file (<4 bytes) → not BGZF.
	tinyPath := writeTempFile(t, dir, "tiny.bin", []byte{0x1f, 0x8b})
	got, err = isBGZFFile(tinyPath)
	if err != nil {
		t.Fatalf("isBGZFFile(tiny): %v", err)
	}
	if got {
		t.Errorf("isBGZFFile on truncated header: got true, want false")
	}
}

func TestLoadGZIBlockOffsets_RoundTrip(t *testing.T) {
	// Hand-build a .gzi with two entries beyond the implicit (0,0).
	var buf bytes.Buffer
	var hdr [8]byte
	binary.LittleEndian.PutUint64(hdr[:], 2)
	buf.Write(hdr[:])
	var entry [16]byte
	binary.LittleEndian.PutUint64(entry[0:8], 100)
	binary.LittleEndian.PutUint64(entry[8:16], 65280)
	buf.Write(entry[:])
	binary.LittleEndian.PutUint64(entry[0:8], 240)
	binary.LittleEndian.PutUint64(entry[8:16], 130560)
	buf.Write(entry[:])

	entries, err := readGZIBlockOffsets(&buf)
	if err != nil {
		t.Fatalf("readGZIBlockOffsets: %v", err)
	}
	if got, want := len(entries), 2; got != want {
		t.Fatalf("entries: got %d, want %d", got, want)
	}
	if entries[0].CompressedOffset != 100 || entries[0].UncompressedOffset != 65280 {
		t.Errorf("entry[0] = %+v, want {100, 65280}", entries[0])
	}
	if entries[1].CompressedOffset != 240 || entries[1].UncompressedOffset != 130560 {
		t.Errorf("entry[1] = %+v, want {240, 130560}", entries[1])
	}
}

func TestLoadGZIBlockOffsets_ShortFile(t *testing.T) {
	if _, err := readGZIBlockOffsets(bytes.NewReader([]byte{0x01, 0x02})); err == nil {
		t.Fatal("expected error on short .gzi header, got nil")
	}
	// Header says 2 entries, but only one is provided.
	var buf bytes.Buffer
	var hdr [8]byte
	binary.LittleEndian.PutUint64(hdr[:], 2)
	buf.Write(hdr[:])
	var entry [16]byte
	binary.LittleEndian.PutUint64(entry[0:8], 100)
	binary.LittleEndian.PutUint64(entry[8:16], 65280)
	buf.Write(entry[:])
	if _, err := readGZIBlockOffsets(&buf); err == nil {
		t.Fatal("expected error on missing entry, got nil")
	}
}

func TestDecompressBGZF_RoundTrip(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1\nACGT\nACGT\n>chr2\nNNNN\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))
	got, err := decompressBGZF(p)
	if err != nil {
		t.Fatalf("decompressBGZF: %v", err)
	}
	if !bytes.Equal(got, payload) {
		t.Errorf("payload mismatch: got %q want %q", got, payload)
	}
}

func TestOpenRandomAccessBGZF_Fetch(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nNNNNNNNNNN\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))

	ra, err := OpenRandomAccessBGZF(p)
	if err != nil {
		t.Fatalf("OpenRandomAccessBGZF: %v", err)
	}
	defer ra.Close()

	cases := []struct {
		name       string
		start, end int64
		want       string
	}{
		{"chr1", 0, 4, "ACGT"},
		{"chr1", 8, 12, "ACGT"}, // spans the line break in the file
		{"chr1", 0, 20, "ACGTACGTACGTACGTACGT"},
		{"chr2", 0, 10, "NNNNNNNNNN"},
	}
	for _, c := range cases {
		got, err := ra.Fetch(c.name, c.start, c.end)
		if err != nil {
			t.Errorf("Fetch(%s, %d, %d): %v", c.name, c.start, c.end, err)
			continue
		}
		if string(got) != c.want {
			t.Errorf("Fetch(%s, %d, %d) = %q, want %q", c.name, c.start, c.end, got, c.want)
		}
	}

	// Through-the-roof OpenRandomAccess should auto-route to the BGZF
	// implementation by sniffing the magic.
	ra2, err := OpenRandomAccess(p)
	if err != nil {
		t.Fatalf("OpenRandomAccess (auto BGZF): %v", err)
	}
	defer ra2.Close()
	seq, err := ra2.Fetch("chr1", 0, 4)
	if err != nil {
		t.Fatalf("OpenRandomAccess Fetch: %v", err)
	}
	if string(seq) != "ACGT" {
		t.Errorf("auto-BGZF Fetch = %q, want %q", seq, "ACGT")
	}
}

func TestOpenRandomAccessBGZF_RejectsPlainFASTA(t *testing.T) {
	dir := t.TempDir()
	p := writeTempFile(t, dir, "plain.fa", []byte(">chr1\nACGT\n"))
	if _, err := OpenRandomAccessBGZF(p); err == nil {
		t.Fatal("expected error opening plain FASTA via BGZF API, got nil")
	}
}

func TestOpenRandomAccessBGZFFullHeader(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1 assembled by consortium X\nACGT\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))
	ra, err := OpenRandomAccessBGZFFullHeader(p)
	if err != nil {
		t.Fatalf("OpenRandomAccessBGZFFullHeader: %v", err)
	}
	defer ra.Close()
	got, err := ra.Fetch("chr1 assembled by consortium X", 0, 4)
	if err != nil {
		t.Fatalf("Fetch full-header contig: %v", err)
	}
	if string(got) != "ACGT" {
		t.Errorf("Fetch = %q, want %q", got, "ACGT")
	}
}

func TestOpenRandomAccessBGZF_BadGZI(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1\nACGT\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))
	// Drop a broken .gzi next to it.
	if err := os.WriteFile(p+".gzi", []byte{0x01, 0x02}, 0o644); err != nil {
		t.Fatalf("write bad .gzi: %v", err)
	}
	if _, err := OpenRandomAccessBGZF(p); err == nil {
		t.Fatal("expected error on bad .gzi, got nil")
	}
}

// Sanity-check the BuildIndex{,FullHeader}Bytes exports.
func TestBuildIndexBytes(t *testing.T) {
	payload := []byte(">chr1\nACGT\n>chr2\nNNNN\n")
	idx, err := BuildIndexBytes(payload)
	if err != nil {
		t.Fatalf("BuildIndexBytes: %v", err)
	}
	if got := len(idx.Entries()); got != 2 {
		t.Errorf("entries: got %d, want 2", got)
	}
	idxF, err := BuildIndexFullHeaderBytes([]byte(">chr1 X\nACGT\n"))
	if err != nil {
		t.Fatalf("BuildIndexFullHeaderBytes: %v", err)
	}
	if _, ok := idxF.Get("chr1 X"); !ok {
		t.Errorf("full-header index missing %q", "chr1 X")
	}
}
