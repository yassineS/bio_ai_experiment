package giab

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestParseBED_Membership(t *testing.T) {
	rs, err := ParseBED(strings.NewReader("chr1\t100\t200\nchr2\t0\t50\n"))
	if err != nil {
		t.Fatal(err)
	}
	if rs.Count() != 2 {
		t.Fatalf("count: %d", rs.Count())
	}
	// BED [100,200) on chr1 -> 0-based 100..199 -> VCF POS 101..200.
	cases := []struct {
		chrom string
		pos   int
		want  bool
	}{
		{"chr1", 100, false}, // 0-based 99, < 100, outside
		{"chr1", 101, true},  // 0-based 100, inside
		{"chr1", 200, true},  // 0-based 199, inside
		{"chr1", 201, false}, // 0-based 200, == end, outside
		{"chr2", 1, true},    // 0-based 0, inside [0,50)
		{"chr2", 50, true},   // 0-based 49, inside
		{"chr2", 51, false},  // 0-based 50, outside
		{"chr3", 10, false},  // absent chrom
	}
	for _, c := range cases {
		if got := rs.Contains(c.chrom, c.pos); got != c.want {
			t.Errorf("Contains(%s,%d)=%v want %v", c.chrom, c.pos, got, c.want)
		}
	}
}

func TestParseBED_MergesOverlaps(t *testing.T) {
	rs, err := ParseBED(strings.NewReader("chr1\t100\t200\nchr1\t150\t300\nchr1\t1000\t1100\n"))
	if err != nil {
		t.Fatal(err)
	}
	if rs.Count() != 2 {
		t.Fatalf("expected 2 merged intervals, got %d", rs.Count())
	}
	if !rs.Contains("chr1", 250) { // inside the merged [100,300)
		t.Fatal("merged interval should contain pos 250")
	}
}

func TestParseBED_SkipsHeaders(t *testing.T) {
	in := "track name=foo\nbrowser position chr1\n# comment\nchr1\t10\t20\n"
	rs, err := ParseBED(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	if rs.Count() != 1 {
		t.Fatalf("header lines should be skipped, count=%d", rs.Count())
	}
}

func TestParseBEDFile_Gzip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "strat.bed.gz")
	var buf bytes.Buffer
	gz := gzip.NewWriter(&buf)
	gz.Write([]byte("chr1\t0\t100\n"))
	gz.Close()
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	rs, err := ParseBEDFile(path)
	if err != nil {
		t.Fatalf("ParseBEDFile gz: %v", err)
	}
	if !rs.Contains("chr1", 50) {
		t.Fatal("gzipped BED not parsed")
	}
}

func TestRegionSet_EmptyNil(t *testing.T) {
	var rs *RegionSet
	if !rs.Empty() {
		t.Fatal("nil RegionSet should be empty")
	}
	if rs.Contains("chr1", 1) {
		t.Fatal("nil RegionSet contains nothing")
	}
	if rs.Count() != 0 {
		t.Fatal("nil count is 0")
	}
}
