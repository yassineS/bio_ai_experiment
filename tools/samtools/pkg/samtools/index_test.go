package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// sortedSAM is a tiny coordinate-sorted dataset that exercises bin/linear
// behaviour across two chromosomes plus an unmapped record.
const sortedSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:1000
r1	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	200	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr2	50	60	5M	*	0	0	ACGTA	IIIII
u1	4	*	0	0	*	*	0	0	*	*
`

// samToBAM converts text SAM bytes to a BGZF-wrapped BAM byte slice for
// use as test input to Index.
func samToBAM(t *testing.T, text string) []byte {
	t.Helper()
	r, err := sam.NewSAMReader(strings.NewReader(text))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	var bam bytes.Buffer
	bw := sam.NewBAMWriter(&bam)
	if err := bw.WriteHeader(r.Header()); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("SAM Read: %v", err)
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("BAM Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("BAM Close: %v", err)
	}
	return bam.Bytes()
}

func TestIndexBuildsBAI(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	var idxBuf bytes.Buffer
	if err := Index(bytes.NewReader(bamBytes), &idxBuf, IndexOptions{}); err != nil {
		t.Fatalf("Index: %v", err)
	}
	idx, err := bam.ReadBAI(bytes.NewReader(idxBuf.Bytes()))
	if err != nil {
		t.Fatalf("bam.ReadBAI: %v", err)
	}
	if len(idx.Refs) != 2 {
		t.Fatalf("Refs len: got %d, want 2", len(idx.Refs))
	}
	if idx.NoCoor != 1 {
		t.Errorf("NoCoor: got %d, want 1 (u1)", idx.NoCoor)
	}
	// Both refs should have at least one bin and a meta pseudo-bin.
	for i, ref := range idx.Refs {
		var sawMeta, sawData bool
		for _, b := range ref.Bins {
			if b.BinID == bam.MetaBin {
				sawMeta = true
			} else {
				sawData = true
			}
		}
		if !sawMeta {
			t.Errorf("ref %d: missing meta pseudo-bin", i)
		}
		if !sawData {
			t.Errorf("ref %d: missing data bin", i)
		}
	}
	// Meta counts on ref 0: 2 mapped, 0 unmapped.
	mapped, unmapped, ok := idx.MetaCounts(0)
	if !ok || mapped != 2 || unmapped != 0 {
		t.Errorf("ref 0 meta counts: ok=%v mapped=%d unmapped=%d, want ok=true 2 0", ok, mapped, unmapped)
	}
}

func TestIndexCSIRejected(t *testing.T) {
	if err := Index(nil, &bytes.Buffer{}, IndexOptions{SelectCSI: true}); err != ErrCSIUnsupported {
		t.Errorf("expected ErrCSIUnsupported, got %v", err)
	}
}

func TestIndexFileRoundTrip(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if err := IndexFile(bamPath, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile: %v", err)
	}
	baiPath := bamPath + ".bai"
	raw, err := os.ReadFile(baiPath)
	if err != nil {
		t.Fatalf("read bai: %v", err)
	}
	if !bytes.Equal(raw[:4], bam.BAIMagic[:]) {
		t.Errorf("bai magic: got % x, want % x", raw[:4], bam.BAIMagic)
	}
}

func TestIndexFileCustomOutput(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	baiPath := filepath.Join(dir, "custom.bai")
	if err := IndexFile(bamPath, baiPath, IndexOptions{}); err != nil {
		t.Fatalf("IndexFile: %v", err)
	}
	if _, err := os.Stat(baiPath); err != nil {
		t.Errorf("custom bai not created: %v", err)
	}
}

func TestIndexFileMissingInput(t *testing.T) {
	if err := IndexFile("/nonexistent/path.bam", "", IndexOptions{}); err == nil {
		t.Error("expected error on missing input")
	}
}

// TestIndexThenRegionQuery is the end-to-end test: build a BAI from a BAM,
// then use it through the region-query path to fetch a chromosome slice.
func TestIndexThenRegionQuery(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if err := IndexFile(bamPath, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile: %v", err)
	}
	// Query chr1:50-150 should return r1 only.
	var out bytes.Buffer
	n, err := ViewFile(bamPath, &out, ViewOptions{
		Regions: []string{"chr1:50-150"},
		Count:   true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile: %v", err)
	}
	if n != 1 {
		t.Errorf("chr1:50-150 count: got %d, want 1", n)
	}
	// Query an unknown chromosome yields zero matches without error.
	var out2 bytes.Buffer
	n2, err := ViewFile(bamPath, &out2, ViewOptions{
		Regions: []string{"chrUnknown:1-1000"},
		Count:   true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile chrUnknown: %v", err)
	}
	if n2 != 0 {
		t.Errorf("chrUnknown count: got %d, want 0", n2)
	}
}

// TestIndexFileCSIRejected confirms that IndexFile surfaces the CSI
// deferral error before opening the input.
func TestIndexFileCSIRejected(t *testing.T) {
	if err := IndexFile("anything", "", IndexOptions{SelectCSI: true}); err != ErrCSIUnsupported {
		t.Errorf("expected ErrCSIUnsupported, got %v", err)
	}
}

// TestViewFileNoRegionsStreams covers the no-regions path of ViewFile,
// which delegates straight to View.
func TestViewFileNoRegionsStreams(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := ViewFile(bamPath, &out, ViewOptions{Count: true}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile: %v", err)
	}
	if n != 4 {
		t.Errorf("count: got %d, want 4", n)
	}
}

// TestViewFileMissingFile surfaces a clean error when the input path is
// missing.
func TestViewFileMissingFile(t *testing.T) {
	if _, err := ViewFile("/nope/nonexistent.bam", &bytes.Buffer{}, ViewOptions{}, io.Discard); err == nil {
		t.Error("expected error for missing file")
	}
}

// TestViewFileFallsBackWithoutIndex confirms that ViewFile with regions on
// a BAM that has no .bai sibling falls back to a linear scan + warning.
func TestViewFileFallsBackWithoutIndex(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	var warn bytes.Buffer
	var out bytes.Buffer
	n, err := ViewFile(bamPath, &out, ViewOptions{
		Regions: []string{"chr1:50-150"},
		Count:   true,
	}, &warn)
	if err != nil {
		t.Fatalf("ViewFile: %v", err)
	}
	if n != 1 {
		t.Errorf("fallback count: got %d, want 1", n)
	}
	if !strings.Contains(warn.String(), "no index") {
		t.Errorf("expected 'no index' warning, got %q", warn.String())
	}
}
