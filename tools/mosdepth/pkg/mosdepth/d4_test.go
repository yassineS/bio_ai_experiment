package mosdepth

import (
	"os"
	"path/filepath"
	"testing"
)

// TestD4WriterReaderRoundTrip exercises the low-level D4 writer/reader on a
// hand-built two-chromosome track, independent of the BAM pipeline.
func TestD4WriterReaderRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "track.d4")
	chroms := []d4Chrom{
		{Name: "chrA", Length: 5},
		{Name: "chrB", Length: 3},
	}
	w, err := newD4Writer(path, chroms)
	if err != nil {
		t.Fatalf("newD4Writer: %v", err)
	}
	// 100000 exceeds the SimpleRange{0,128} dictionary, so it is clamped to
	// the sentinel code 127 in the primary table — exactly as upstream
	// mosdepth's per-base D4 output caps over-dictionary depths.
	depthsA := []int32{0, 1, 2, 100000, 0}
	wantA := []int32{0, 1, 2, 127, 0}
	depthsB := []int32{7, 7, 0}
	if err := w.writeChrom("chrA", depthsA); err != nil {
		t.Fatalf("writeChrom chrA: %v", err)
	}
	if err := w.writeChrom("chrB", depthsB); err != nil {
		t.Fatalf("writeChrom chrB: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r, err := openD4Reader(path)
	if err != nil {
		t.Fatalf("openD4Reader: %v", err)
	}
	gotA, err := r.chromDepths("chrA")
	if err != nil {
		t.Fatalf("chromDepths chrA: %v", err)
	}
	gotB, err := r.chromDepths("chrB")
	if err != nil {
		t.Fatalf("chromDepths chrB: %v", err)
	}
	if !equalInt32(gotA, wantA) {
		t.Errorf("chrA: got %v, want %v", gotA, wantA)
	}
	if !equalInt32(gotB, depthsB) {
		t.Errorf("chrB: got %v, want %v", gotB, depthsB)
	}
	if _, err := r.chromDepths("missing"); err == nil {
		t.Errorf("expected error for missing chromosome")
	}
}

// TestD4WriterOrdering rejects out-of-order and wrong-length chromosome writes.
func TestD4WriterOrdering(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "track.d4")
	chroms := []d4Chrom{{Name: "chrA", Length: 2}, {Name: "chrB", Length: 2}}
	w, err := newD4Writer(path, chroms)
	if err != nil {
		t.Fatalf("newD4Writer: %v", err)
	}
	if err := w.writeChrom("chrB", []int32{1, 1}); err == nil {
		t.Errorf("expected out-of-order error")
	}
	if err := w.writeChrom("chrA", []int32{1}); err == nil {
		t.Errorf("expected wrong-length error")
	}
	// Closing before all chromosomes are written should error.
	if err := w.Close(); err == nil {
		t.Errorf("expected incomplete-close error")
	}
}

// TestD4ReaderBadMagic rejects a file that is not a D4 container.
func TestD4ReaderBadMagic(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "bad.d4")
	if err := os.WriteFile(path, []byte("not a d4 file at all really"), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if _, err := openD4Reader(path); err == nil {
		t.Errorf("expected bad-magic error")
	}
}

// TestD4DenseDepths confirms dense materialisation matches an emit() sweep.
func TestD4DenseDepths(t *testing.T) {
	a := newCovAccum(10)
	a.add(2, 5) // depth 1 over [2,5)
	a.add(3, 4) // +1 over [3,4) -> depth 2 there
	got := d4DenseDepths(a)
	want := []int32{0, 0, 1, 2, 1, 0, 0, 0, 0, 0}
	if !equalInt32(got, want) {
		t.Errorf("d4DenseDepths: got %v, want %v", got, want)
	}
}

func equalInt32(a, b []int32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
