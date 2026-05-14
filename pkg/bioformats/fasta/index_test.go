package fasta

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestReadIndex(t *testing.T) {
	raw := "chr1\t100\t6\t60\t61\nchr2\t50\t117\t50\t51\n"
	idx, err := ReadIndex(strings.NewReader(raw))
	if err != nil {
		t.Fatalf("ReadIndex: %v", err)
	}
	if got := idx.Names(); len(got) != 2 || got[0] != "chr1" || got[1] != "chr2" {
		t.Fatalf("Names() = %v", got)
	}
	e, ok := idx.Get("chr2")
	if !ok {
		t.Fatalf("Get(chr2) missing")
	}
	if e.Length != 50 || e.Offset != 117 || e.LineBases != 50 || e.LineWidth != 51 {
		t.Fatalf("chr2 entry mismatch: %+v", e)
	}
	if _, ok := idx.Get("missing"); ok {
		t.Fatalf("Get(missing) should be false")
	}
}

func TestReadIndexBadInput(t *testing.T) {
	cases := []struct {
		name string
		raw  string
	}{
		{"too few cols", "chr1\t100\n"},
		{"bad length", "chr1\tNOPE\t6\t60\t61\n"},
		{"bad offset", "chr1\t100\tNOPE\t60\t61\n"},
		{"bad lineblen", "chr1\t100\t6\tNOPE\t61\n"},
		{"bad linelen", "chr1\t100\t6\t60\tNOPE\n"},
		{"zero lineblen", "chr1\t100\t6\t0\t1\n"},
		{"linelen < lineblen", "chr1\t100\t6\t60\t10\n"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := ReadIndex(strings.NewReader(tc.raw)); err == nil {
				t.Fatalf("expected error for %q", tc.raw)
			}
		})
	}
}

const tinyFasta = ">chr1\nACGTAC\nGTACGT\nAC\n>chr2\nTTTTTTTTTT\n"

func writeTempFasta(t *testing.T, body string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	return path
}

func TestBuildIndex(t *testing.T) {
	path := writeTempFasta(t, tinyFasta)
	idx, err := BuildIndex(path)
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	if names := idx.Names(); len(names) != 2 || names[0] != "chr1" || names[1] != "chr2" {
		t.Fatalf("Names() = %v", names)
	}
	chr1, _ := idx.Get("chr1")
	if chr1.Length != 14 {
		t.Fatalf("chr1 length = %d (want 14)", chr1.Length)
	}
	if chr1.LineBases != 6 || chr1.LineWidth != 7 {
		t.Fatalf("chr1 line geometry = (%d,%d), want (6,7)", chr1.LineBases, chr1.LineWidth)
	}
	chr2, _ := idx.Get("chr2")
	if chr2.Length != 10 {
		t.Fatalf("chr2 length = %d (want 10)", chr2.Length)
	}
}

func TestBuildIndexBadGeometry(t *testing.T) {
	// Two full-width lines surrounding a short one; the short one is fine
	// only if it's the final line. This case puts a short line in the middle.
	bad := ">c\nAAA\nA\nAAA\n"
	path := writeTempFasta(t, bad)
	if _, err := BuildIndex(path); err == nil {
		t.Fatalf("expected non-uniform line width error")
	}
}

func TestIndexWriteToAndSave(t *testing.T) {
	raw := "chr1\t100\t6\t60\t61\nchr2\t50\t117\t50\t51\n"
	idx, err := ReadIndex(strings.NewReader(raw))
	if err != nil {
		t.Fatalf("ReadIndex: %v", err)
	}
	var buf bytes.Buffer
	if _, err := idx.WriteTo(&buf); err != nil {
		t.Fatalf("WriteTo: %v", err)
	}
	if buf.String() != raw {
		t.Fatalf("roundtrip mismatch:\nwant %q\n got %q", raw, buf.String())
	}

	// Save round-trip.
	dir := t.TempDir()
	out := filepath.Join(dir, "out.fai")
	if err := idx.Save(out); err != nil {
		t.Fatalf("Save: %v", err)
	}
	got, err := os.ReadFile(out)
	if err != nil {
		t.Fatalf("read back: %v", err)
	}
	if string(got) != raw {
		t.Fatalf("save mismatch: %q vs %q", got, raw)
	}
}

func TestEntries(t *testing.T) {
	raw := "chr1\t100\t6\t60\t61\n"
	idx, err := ReadIndex(strings.NewReader(raw))
	if err != nil {
		t.Fatalf("ReadIndex: %v", err)
	}
	es := idx.Entries()
	if len(es) != 1 || es[0].Name != "chr1" {
		t.Fatalf("Entries() = %+v", es)
	}
}

func TestRandomAccessFetch(t *testing.T) {
	path := writeTempFasta(t, tinyFasta)
	ra, err := OpenRandomAccess(path)
	if err != nil {
		t.Fatalf("OpenRandomAccess: %v", err)
	}
	defer ra.Close()

	// Full chr1: ACGTACGTACGTAC (14 bases)
	full, err := ra.Fetch("chr1", 0, 14)
	if err != nil {
		t.Fatalf("Fetch full chr1: %v", err)
	}
	if string(full) != "ACGTACGTACGTAC" {
		t.Fatalf("full chr1 = %q", full)
	}

	// Sub-range straddling a line break.
	mid, err := ra.Fetch("chr1", 4, 9)
	if err != nil {
		t.Fatalf("Fetch mid: %v", err)
	}
	if string(mid) != "ACGTA" {
		t.Fatalf("chr1[4:9] = %q want ACGTA", mid)
	}

	// Empty range.
	empty, err := ra.Fetch("chr1", 5, 5)
	if err != nil {
		t.Fatalf("empty fetch: %v", err)
	}
	if len(empty) != 0 {
		t.Fatalf("empty fetch = %q", empty)
	}

	if l := ra.Length("chr1"); l != 14 {
		t.Fatalf("Length(chr1) = %d", l)
	}
	if l := ra.Length("missing"); l != -1 {
		t.Fatalf("Length(missing) = %d", l)
	}
}

func TestRandomAccessFetchErrors(t *testing.T) {
	path := writeTempFasta(t, tinyFasta)
	ra, err := OpenRandomAccess(path)
	if err != nil {
		t.Fatalf("OpenRandomAccess: %v", err)
	}
	defer ra.Close()

	if _, err := ra.Fetch("ghost", 0, 1); err == nil {
		t.Fatalf("expected missing-contig error")
	}
	if _, err := ra.Fetch("chr1", -1, 2); err == nil {
		t.Fatalf("expected out-of-bounds error (start<0)")
	}
	if _, err := ra.Fetch("chr1", 0, 999); err == nil {
		t.Fatalf("expected out-of-bounds error (end too big)")
	}
	if _, err := ra.Fetch("chr1", 3, 2); err == nil {
		t.Fatalf("expected end<start error")
	}
}

func TestRandomAccessUsesSidecarIndex(t *testing.T) {
	path := writeTempFasta(t, tinyFasta)
	// Drop a .fai sidecar that lies about lengths to prove it gets used.
	if err := os.WriteFile(path+".fai", []byte("chr1\t6\t6\t6\t7\n"), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}
	ra, err := OpenRandomAccess(path)
	if err != nil {
		t.Fatalf("OpenRandomAccess: %v", err)
	}
	defer ra.Close()
	if l := ra.Length("chr1"); l != 6 {
		t.Fatalf("sidecar length not honoured: %d", l)
	}
}

func TestRandomAccessWrappingReaderAt(t *testing.T) {
	// Exercise NewRandomAccess on an in-memory FASTA + index.
	body := []byte(">x\nACGTACGT\n")
	idx, err := ReadIndex(strings.NewReader("x\t8\t3\t8\t9\n"))
	if err != nil {
		t.Fatalf("read idx: %v", err)
	}
	ra := NewRandomAccess(bytes.NewReader(body), idx)
	defer ra.Close()
	got, err := ra.Fetch("x", 2, 6)
	if err != nil {
		t.Fatalf("Fetch: %v", err)
	}
	if string(got) != "GTAC" {
		t.Fatalf("got %q", got)
	}
}
