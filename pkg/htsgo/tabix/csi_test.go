package tabix

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

func TestCSINewDefaults(t *testing.T) {
	c := NewCSI(0, 0)
	if c.MinShift != 14 || c.Depth != 5 {
		t.Fatalf("default (MinShift, Depth) = (%d, %d); want (14, 5)", c.MinShift, c.Depth)
	}
}

func TestCSIBinMath(t *testing.T) {
	c := NewCSI(14, 5)

	// BAI's tile size is 1<<14 = 16 KiB; the smallest bin level has 4681 base.
	// A zero-length range collapses to bin 0 in this implementation.
	if got := c.Reg2bin(0, 1); got == 0 {
		t.Fatal("Reg2bin(0,1) should not be 0; want a smallest-level bin")
	}

	// A range spanning the full chromosome must collapse to bin 0.
	if got := c.Reg2bin(0, c.MaxPos()); got != 0 {
		t.Fatalf("Reg2bin(0, MaxPos) = %d; want 0", got)
	}

	bins := c.Reg2bins(0, 1<<15)
	if len(bins) == 0 || bins[0] != 0 {
		t.Fatalf("Reg2bins(0, 32768) should start with bin 0; got %v", bins)
	}
}

func TestCSIRoundTrip(t *testing.T) {
	c := NewCSI(14, 5)
	c.Names = []string{"chr1", "chr2"}
	c.AddRecord(0, 100, 200, MakeVOffset(0, 10), MakeVOffset(0, 50))
	c.AddRecord(0, 5000, 5100, MakeVOffset(0, 60), MakeVOffset(0, 90))
	c.AddRecord(1, 1, 2, MakeVOffset(100, 0), MakeVOffset(100, 25))
	c.NoCoor = 7

	var buf bytes.Buffer
	if err := c.Write(&buf); err != nil {
		t.Fatalf("Write: %v", err)
	}
	got, err := ReadCSI(&buf)
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if got.MinShift != c.MinShift || got.Depth != c.Depth {
		t.Errorf("MinShift/Depth lost: got (%d,%d) want (%d,%d)", got.MinShift, got.Depth, c.MinShift, c.Depth)
	}
	if len(got.Refs) != len(c.Refs) {
		t.Fatalf("refs: got %d want %d", len(got.Refs), len(c.Refs))
	}
	if got.NoCoor != 7 {
		t.Errorf("NoCoor: got %d want 7", got.NoCoor)
	}
}

func TestCSIWriteReadFile(t *testing.T) {
	c := NewCSI(14, 5)
	cfg := Config{Format: 2, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}
	c.SetAuxFromTabix(cfg, []string{"chr1"})
	c.AddRecord(0, 0, 1024, MakeVOffset(0, 0), MakeVOffset(0, 64))

	dir := t.TempDir()
	path := filepath.Join(dir, "test.csi")
	if err := c.WriteFile(path); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("expected file at %s, got %v", path, err)
	}
	got, err := ReadCSIFile(path)
	if err != nil {
		t.Fatalf("ReadCSIFile: %v", err)
	}
	if got.MinShift != c.MinShift {
		t.Errorf("MinShift round-trip: got %d want %d", got.MinShift, c.MinShift)
	}
	if len(got.Names) == 0 || got.Names[0] != "chr1" {
		t.Errorf("Names[0]: got %v want [\"chr1\", ...]", got.Names)
	}
}

func TestCSIRegionChunksAndChromID(t *testing.T) {
	c := NewCSI(14, 5)
	c.Names = []string{"chr1", "chr2"}
	beg, end := int64(100), int64(200)
	c.AddRecord(0, beg, end, MakeVOffset(0, 10), MakeVOffset(0, 60))

	chunks := c.RegionChunks(0, beg, end)
	if len(chunks) == 0 {
		t.Fatal("RegionChunks returned nothing for in-range query")
	}
	if c2 := c.RegionChunks(0, 1_000_000, 1_000_100); c2 != nil {
		t.Errorf("RegionChunks for empty region returned %v; want nil", c2)
	}
	if c2 := c.RegionChunks(-1, 0, 10); c2 != nil {
		t.Errorf("RegionChunks(-1) returned %v; want nil", c2)
	}

	if got := ChromIDInCSI(c, "chr2"); got != 1 {
		t.Errorf("ChromIDInCSI(\"chr2\") = %d; want 1", got)
	}
	if got := ChromIDInCSI(c, "missing"); got != -1 {
		t.Errorf("ChromIDInCSI(missing) = %d; want -1", got)
	}
}

func TestCSITabixAuxRoundTrip(t *testing.T) {
	c := NewCSI(14, 5)
	cfg := Config{Format: 2, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}
	c.SetAuxFromTabix(cfg, []string{"chr1", "chr2", "chr3"})

	got := parseCSIAuxNames(c.Aux)
	want := []string{"chr1", "chr2", "chr3"}
	if len(got) != len(want) {
		t.Fatalf("parseCSIAuxNames got %v want %v", got, want)
	}
	for i := range got {
		if got[i] != want[i] {
			t.Errorf("name %d: got %q want %q", i, got[i], want[i])
		}
	}

	if names := parseCSIAuxNames(nil); names != nil {
		t.Errorf("parseCSIAuxNames(nil) = %v; want nil", names)
	}
}

func TestCSIBuildFromDataFile(t *testing.T) {
	// Build a tiny bgzipped VCF-like file and index it via BuildCSIFromDataFile.
	dir := t.TempDir()
	dataPath := filepath.Join(dir, "data.vcf.gz")

	f, err := os.Create(dataPath)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	bw := bgzip.NewWriter(f)
	rows := []string{
		"##fileformat=VCFv4.3\n",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		"chr1\t100\t.\tA\tT\t.\tPASS\t.\n",
		"chr1\t250\t.\tG\tC\t.\tPASS\t.\n",
		"chr2\t10\t.\tA\tT\t.\tPASS\t.\n",
	}
	for _, r := range rows {
		if _, err := bw.Write([]byte(r)); err != nil {
			t.Fatalf("bgzip write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	f.Close()

	cfg := Config{Format: 2, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}
	c, err := BuildCSIFromDataFile(dataPath, cfg, 14)
	if err != nil {
		t.Fatalf("BuildCSIFromDataFile: %v", err)
	}
	if len(c.Refs) < 2 {
		t.Errorf("CSI refs: got %d want >=2", len(c.Refs))
	}
	if len(c.Names) < 2 || c.Names[0] != "chr1" || c.Names[1] != "chr2" {
		t.Errorf("CSI names: got %v want [chr1 chr2 ...]", c.Names)
	}
}
