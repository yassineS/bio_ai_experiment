package tabix

import (
	"os"
	"path/filepath"
	"testing"
)

// writeFile writes text to a fresh file named name inside t.TempDir() and
// returns its absolute path.
func writeFile(t *testing.T, name, text string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(text), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return path
}

func TestLoadTargetsTabCoordinates(t *testing.T) {
	// Tab format: 1-based inclusive. "chr1 100 150" => 0-based [99, 150).
	path := writeFile(t, "t.txt", "chr1\t100\t150\nchr2\t10\t10\n")
	tg, err := LoadTargets(path)
	if err != nil {
		t.Fatalf("LoadTargets: %v", err)
	}
	cases := []struct {
		chrom    string
		beg, end int
		want     bool
	}{
		{"chr1", 99, 100, true},   // first base of the interval
		{"chr1", 149, 150, true},  // last base (1-based 150)
		{"chr1", 150, 151, false}, // just past the inclusive end
		{"chr1", 98, 99, false},   // just before the start
		{"chr2", 9, 10, true},     // single-position target
		{"chr2", 10, 11, false},
		{"chrX", 0, 1000, false}, // unknown chromosome
	}
	for _, c := range cases {
		if got := tg.Overlaps(c.chrom, c.beg, c.end); got != c.want {
			t.Errorf("Overlaps(%s,%d,%d)=%v want %v", c.chrom, c.beg, c.end, got, c.want)
		}
	}
}

func TestLoadTargetsBEDCoordinates(t *testing.T) {
	// BED format: 0-based half-open. "chr1 149 200" => [149, 200).
	path := writeFile(t, "t.bed", "chr1\t149\t200\n")
	tg, err := LoadTargets(path)
	if err != nil {
		t.Fatalf("LoadTargets: %v", err)
	}
	cases := []struct {
		beg, end int
		want     bool
	}{
		{149, 150, true},  // first base
		{199, 200, true},  // last base (exclusive end 200)
		{200, 201, false}, // exclusive end excluded
		{148, 149, false},
	}
	for _, c := range cases {
		if got := tg.Overlaps("chr1", c.beg, c.end); got != c.want {
			t.Errorf("Overlaps(chr1,%d,%d)=%v want %v", c.beg, c.end, got, c.want)
		}
	}
}

func TestTargetsBoundaryAndAdjacent(t *testing.T) {
	// Adjacent intervals [99,150) and [199,250) (tab, 1-based 100-150,200-250).
	path := writeFile(t, "t.txt", "chr1\t100\t150\nchr1\t200\t250\n")
	tg, err := LoadTargets(path)
	if err != nil {
		t.Fatalf("LoadTargets: %v", err)
	}
	// A record sitting exactly in the gap between the two intervals is out.
	if tg.Overlaps("chr1", 150, 151) {
		t.Errorf("record at 150 should fall in the gap between adjacent intervals")
	}
	// Record straddling the boundary into the second interval overlaps.
	if !tg.Overlaps("chr1", 198, 200) {
		t.Errorf("record straddling into the second interval should overlap")
	}
	// A record that spans across the gap touches both intervals -> overlap.
	if !tg.Overlaps("chr1", 100, 260) {
		t.Errorf("spanning record should overlap")
	}
}

func TestTargetsChromOnlyLine(t *testing.T) {
	path := writeFile(t, "t.txt", "chr2\n")
	tg, err := LoadTargets(path)
	if err != nil {
		t.Fatalf("LoadTargets: %v", err)
	}
	if !tg.Overlaps("chr2", 0, 1) {
		t.Errorf("chrom-only line should select the whole chromosome")
	}
	if !tg.Overlaps("chr2", 1<<29, 1<<29+1) {
		t.Errorf("chrom-only line should cover far positions")
	}
	if tg.Overlaps("chr1", 0, 1) {
		t.Errorf("chrom-only chr2 line must not select chr1")
	}
}

func TestTargetsSkipsCommentsAndBlankLines(t *testing.T) {
	path := writeFile(t, "t.txt", "# comment\n\nchr1\t10\t20\n")
	tg, err := LoadTargets(path)
	if err != nil {
		t.Fatalf("LoadTargets: %v", err)
	}
	if !tg.Overlaps("chr1", 9, 10) {
		t.Errorf("interval after comment/blank lines should be loaded")
	}
}

func TestParseTargetInterval(t *testing.T) {
	cases := []struct {
		name     string
		fields   []string
		bed      bool
		wantOK   bool
		wantBeg  int
		wantEnd  int
		chromOnl bool
	}{
		{"tab range", []string{"c", "100", "150"}, false, true, 99, 150, false},
		{"tab single col2 only", []string{"c", "100"}, false, true, 99, 100, false},
		{"tab zero start invalid", []string{"c", "0", "5"}, false, false, 0, 0, false},
		{"bed range", []string{"c", "149", "200"}, true, true, 149, 200, false},
		{"bed start only", []string{"c", "10"}, true, true, 10, 11, false},
		{"bad begin", []string{"c", "x", "5"}, false, false, 0, 0, false},
		{"chrom only", []string{"c"}, false, true, 0, maxTargetCoor, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			iv, ok := parseTargetInterval(c.fields, c.bed)
			if ok != c.wantOK {
				t.Fatalf("ok=%v want %v", ok, c.wantOK)
			}
			if !ok {
				return
			}
			if iv.beg != c.wantBeg || iv.end != c.wantEnd {
				t.Errorf("got [%d,%d) want [%d,%d)", iv.beg, iv.end, c.wantBeg, c.wantEnd)
			}
		})
	}
}

func TestLoadTargetsMissingFile(t *testing.T) {
	if _, err := LoadTargets(filepath.Join(t.TempDir(), "nope.txt")); err == nil {
		t.Fatalf("expected error for missing targets file")
	}
}
