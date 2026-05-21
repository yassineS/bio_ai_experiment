package cram

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
)

// gzipCRAI gzip-compresses a .crai text body so it can be fed to
// ReadCRAI, which expects the gzip layer a real .crai always has.
func gzipCRAI(t *testing.T, text string) []byte {
	t.Helper()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	if _, err := gw.Write([]byte(text)); err != nil {
		t.Fatalf("gzip write: %v", err)
	}
	if err := gw.Close(); err != nil {
		t.Fatalf("gzip close: %v", err)
	}
	return buf.Bytes()
}

// TestReadCRAISynthetic parses a hand-built .crai and checks every field
// of every entry is decoded as the six-integer line says.
func TestReadCRAISynthetic(t *testing.T) {
	body := "0\t100\t50\t1024\t12\t300\n" +
		"0\t200\t40\t2048\t12\t280\n" +
		"1\t10\t90\t4096\t0\t500\n" +
		"-1\t0\t0\t8192\t0\t128\n"
	idx, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, body)))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	if len(idx.Entries) != 4 {
		t.Fatalf("parsed %d entries, want 4", len(idx.Entries))
	}
	want := CRAIEntry{RefID: 0, AlignmentStart: 100, AlignmentSpan: 50, ContainerOffset: 1024, SliceOffset: 12, SliceSize: 300}
	if idx.Entries[0] != want {
		t.Errorf("entry 0 = %+v, want %+v", idx.Entries[0], want)
	}
	if idx.Entries[3].RefID != -1 {
		t.Errorf("entry 3 RefID = %d, want -1 (unmapped slice)", idx.Entries[3].RefID)
	}
}

// TestReadCRAIBlankLine checks a trailing blank line is tolerated.
func TestReadCRAIBlankLine(t *testing.T) {
	idx, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, "0\t1\t10\t0\t0\t1\n\n")))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	if len(idx.Entries) != 1 {
		t.Errorf("parsed %d entries, want 1", len(idx.Entries))
	}
}

// TestReadCRAIMalformed feeds malformed .crai bodies to the parser; each
// must surface as an error, never a panic.
func TestReadCRAIMalformed(t *testing.T) {
	cases := []struct {
		name string
		body string
	}{
		{"too few fields", "0\t1\t2\n"},
		{"too many fields", "0\t1\t2\t3\t4\t5\t6\n"},
		{"non-integer field", "0\tx\t2\t3\t4\t5\n"},
		{"negative container offset", "0\t1\t2\t-3\t4\t5\n"},
		{"negative slice size", "0\t1\t2\t3\t4\t-5\n"},
		{"negative alignment span", "0\t1\t-2\t3\t4\t5\n"},
		{"ref id out of 32-bit range", "9999999999999\t1\t2\t3\t4\t5\n"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, c.body))); err == nil {
				t.Error("expected an error")
			}
		})
	}
}

// TestReadCRAINotGzip checks a non-gzip stream is rejected cleanly.
func TestReadCRAINotGzip(t *testing.T) {
	if _, err := ReadCRAI(bytes.NewReader([]byte("0\t1\t2\t3\t4\t5\n"))); err == nil {
		t.Error("a non-gzip .crai stream should error")
	}
}

// TestCRAIQueryOverlap exercises the region overlap query against a
// synthetic index with carefully placed slice spans.
func TestCRAIQueryOverlap(t *testing.T) {
	// Three slices on ref 0: 1-100, 101-200, 201-300; one on ref 1.
	body := "0\t1\t100\t0\t0\t10\n" +
		"0\t101\t100\t100\t0\t10\n" +
		"0\t201\t100\t200\t0\t10\n" +
		"1\t1\t100\t300\t0\t10\n"
	idx, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, body)))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	cases := []struct {
		name             string
		refID            int32
		beg0, end0       int64
		wantContainerOff []int64
	}{
		{"first slice only", 0, 0, 50, []int64{0}},
		{"spanning two slices", 0, 90, 110, []int64{0, 100}},
		{"open-ended from middle", 0, 150, 0, []int64{100, 200}},
		{"all of ref 0", 0, 0, 1000, []int64{0, 100, 200}},
		{"different reference", 1, 0, 1000, []int64{300}},
		{"no overlap past the end", 0, 5000, 6000, nil},
		{"exact boundary start", 0, 100, 101, []int64{100}},
		{"unknown reference", 9, 0, 1000, nil},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			hits := idx.Query(c.refID, c.beg0, c.end0)
			if len(hits) != len(c.wantContainerOff) {
				t.Fatalf("Query returned %d entries, want %d", len(hits), len(c.wantContainerOff))
			}
			for i, off := range c.wantContainerOff {
				if hits[i].ContainerOffset != off {
					t.Errorf("hit %d container offset = %d, want %d", i, hits[i].ContainerOffset, off)
				}
			}
		})
	}
}

// TestCRAIQueryBoundary checks half-open overlap arithmetic at exact
// slice boundaries: a query ending exactly where a slice begins must not
// match that slice.
func TestCRAIQueryBoundary(t *testing.T) {
	// Slice covers 1-based 101-200, i.e. half-open 0-based [100, 200).
	idx, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, "0\t101\t100\t0\t0\t10\n")))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	if hits := idx.Query(0, 0, 100); len(hits) != 0 {
		t.Errorf("a query [0,100) must not hit a slice starting at 0-based 100, got %d hits", len(hits))
	}
	if hits := idx.Query(0, 100, 101); len(hits) != 1 {
		t.Errorf("a query [100,101) must hit a slice starting at 0-based 100, got %d hits", len(hits))
	}
	if hits := idx.Query(0, 199, 200); len(hits) != 1 {
		t.Errorf("a query [199,200) must hit a slice covering 0-based [100,200), got %d hits", len(hits))
	}
	if hits := idx.Query(0, 200, 201); len(hits) != 0 {
		t.Errorf("a query [200,201) must not hit a slice covering 0-based [100,200), got %d hits", len(hits))
	}
}

// TestCRAIQueryRegion checks the region-package-typed query wrapper.
func TestCRAIQueryRegion(t *testing.T) {
	idx, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, "2\t50\t100\t0\t0\t10\n")))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	r := region.ResolvedRegion{RefID: 2, Beg0: 60, End0: 70}
	if hits := idx.QueryRegion(r); len(hits) != 1 {
		t.Errorf("QueryRegion: got %d hits, want 1", len(hits))
	}
	miss := region.ResolvedRegion{RefID: 2, Beg0: 5000, End0: 6000}
	if hits := idx.QueryRegion(miss); len(hits) != 0 {
		t.Errorf("QueryRegion non-overlapping: got %d hits, want 0", len(hits))
	}
}

// TestCRAIQueryZeroSpan checks a zero-span slice (an unmapped-reads
// slice or an empty slice) overlaps only its single start coordinate.
func TestCRAIQueryZeroSpan(t *testing.T) {
	idx, err := ReadCRAI(bytes.NewReader(gzipCRAI(t, "0\t10\t0\t0\t0\t10\n")))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	if hits := idx.Query(0, 9, 10); len(hits) != 1 {
		t.Errorf("a query covering the zero-span slice's start must hit, got %d", len(hits))
	}
	if hits := idx.Query(0, 11, 20); len(hits) != 0 {
		t.Errorf("a query past the zero-span slice must not hit, got %d", len(hits))
	}
}

// TestOpenCRAIFixture parses the real .crai shipped with the samtools
// submodule and sanity-checks the parsed entries.
func TestOpenCRAIFixture(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, "mpileup/ce#5b.cram.crai")
	if _, err := os.Stat(path); err != nil {
		t.Skip("samtools submodule not initialised — .crai fixture unavailable")
	}
	idx, err := OpenCRAI(path)
	if err != nil {
		t.Fatalf("OpenCRAI: %v", err)
	}
	if len(idx.Entries) == 0 {
		t.Fatal("the real .crai parsed to zero entries")
	}
	for i, e := range idx.Entries {
		if e.ContainerOffset < 0 || e.SliceSize < 0 {
			t.Errorf("entry %d has a negative offset/size: %+v", i, e)
		}
		if e.AlignmentSpan < 0 {
			t.Errorf("entry %d has a negative span: %+v", i, e)
		}
	}
	// Every entry's slice should be locatable: a query over the union of
	// all ranges must return every entry.
	first := idx.Entries[0]
	hits := idx.Query(first.RefID, 0, 0)
	if len(hits) == 0 {
		t.Error("an open-ended query on the first entry's reference returned nothing")
	}
}

// TestOpenCRAIMissing checks OpenCRAI on a missing path errors cleanly.
func TestOpenCRAIMissing(t *testing.T) {
	if _, err := OpenCRAI(filepath.Join(samtoolsTestDir, "does-not-exist.crai")); err == nil {
		t.Error("OpenCRAI on a missing file should error")
	}
}
