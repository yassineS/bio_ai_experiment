package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"reflect"
	"sort"
	"strconv"
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

// TestIndexBuildsCSI confirms that Index with SelectCSI emits a
// BGZF-compressed .csi index that round-trips through bam.ReadCSI and
// carries the same pseudo-bin metadata the BAI path does.
func TestIndexBuildsCSI(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	var idxBuf bytes.Buffer
	if err := Index(bytes.NewReader(bamBytes), &idxBuf, IndexOptions{SelectCSI: true}); err != nil {
		t.Fatalf("Index CSI: %v", err)
	}
	idx, err := bam.ReadCSI(bytes.NewReader(idxBuf.Bytes()))
	if err != nil {
		t.Fatalf("bam.ReadCSI: %v", err)
	}
	if len(idx.CSI.Refs) != 2 {
		t.Fatalf("Refs len: got %d, want 2", len(idx.CSI.Refs))
	}
	if idx.CSI.NoCoor != 1 {
		t.Errorf("NoCoor: got %d, want 1 (u1)", idx.CSI.NoCoor)
	}
	mapped, unmapped, ok := idx.MetaCounts(0)
	if !ok || mapped != 2 || unmapped != 0 {
		t.Errorf("ref 0 meta counts: ok=%v mapped=%d unmapped=%d, want ok=true 2 0", ok, mapped, unmapped)
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

// TestIndexFileCSI confirms that IndexFile with SelectCSI writes a
// sibling <bam>.csi file and that region queries against it agree with
// the .bai path for chromosomes within BAI's range.
func TestIndexFileCSI(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if err := IndexFile(bamPath, "", IndexOptions{SelectCSI: true}); err != nil {
		t.Fatalf("IndexFile CSI: %v", err)
	}
	csiPath := bamPath + ".csi"
	if _, err := os.Stat(csiPath); err != nil {
		t.Fatalf("csi not created: %v", err)
	}
	idx, err := bam.ReadCSIFile(csiPath)
	if err != nil {
		t.Fatalf("ReadCSIFile: %v", err)
	}
	if len(idx.CSI.Refs) != 2 {
		t.Errorf("Refs len: got %d, want 2", len(idx.CSI.Refs))
	}
	// Region query through the .csi index returns the same record as the
	// .bai path does for chr1:50-150 (see TestIndexThenRegionQuery).
	var out bytes.Buffer
	n, err := ViewFile(bamPath, &out, ViewOptions{
		Regions: []string{"chr1:50-150"},
		Count:   true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile via csi: %v", err)
	}
	if n != 1 {
		t.Errorf("chr1:50-150 count via csi: got %d, want 1", n)
	}
}

// regionQNAMEs runs an indexed region query against the BAM at bamPath and
// returns the sorted set of QNAMEs of the emitted records. ViewFile without
// Count writes headerless SAM, so the QNAME is the first tab-separated
// field of every output line.
func regionQNAMEs(t *testing.T, bamPath, reg string) []string {
	t.Helper()
	var out bytes.Buffer
	if _, err := ViewFile(bamPath, &out, ViewOptions{Regions: []string{reg}}, io.Discard); err != nil {
		t.Fatalf("ViewFile %s (%s): %v", bamPath, reg, err)
	}
	var names []string
	for _, line := range strings.Split(strings.TrimRight(out.String(), "\n"), "\n") {
		if line == "" || strings.HasPrefix(line, "@") {
			continue
		}
		names = append(names, line[:strings.IndexByte(line, '\t')])
	}
	sort.Strings(names)
	return names
}

// TestCSIvsBAIRegionParity builds both a .bai and a .csi for the same BAM
// and confirms region queries return the identical set of record QNAMEs
// for every region within BAI's coordinate range. Comparing QNAMEs (not
// just counts) catches a shared bug that yields a same-sized but wrong
// record set.
func TestCSIvsBAIRegionParity(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	regions := []string{"chr1:50-150", "chr1:150-250", "chr1:1-1000", "chr2:1-1000", "chr2:40-60"}

	dirBAI := t.TempDir()
	baiBAM := filepath.Join(dirBAI, "in.bam")
	if err := os.WriteFile(baiBAM, bamBytes, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := IndexFile(baiBAM, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile bai: %v", err)
	}

	dirCSI := t.TempDir()
	csiBAM := filepath.Join(dirCSI, "in.bam")
	if err := os.WriteFile(csiBAM, bamBytes, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := IndexFile(csiBAM, "", IndexOptions{SelectCSI: true}); err != nil {
		t.Fatalf("IndexFile csi: %v", err)
	}

	for _, reg := range regions {
		baiNames := regionQNAMEs(t, baiBAM, reg)
		csiNames := regionQNAMEs(t, csiBAM, reg)
		if !reflect.DeepEqual(baiNames, csiNames) {
			t.Errorf("region %s: bai QNAMEs %v != csi QNAMEs %v", reg, baiNames, csiNames)
		}
	}
}

// TestCSILargeCoordinate exercises a position beyond the BAI 2^29 bp
// ceiling: a record at ~600 Mbp cannot be addressed by BAI but is
// indexed and queried correctly by a CSI built with a larger min_shift.
//
// Note the depth-5 default CSI has the same 2^29 reach as BAI by design
// (it is the BAI-equivalent scheme); the larger coordinate range comes
// from raising min_shift (samtools index -m), which is what this test
// exercises.
func TestCSILargeCoordinate(t *testing.T) {
	// Position 600,000,001 (1-based) is past the BAI 2^29 (~536 Mbp) limit.
	const bigPos = 600000001
	largeSAM := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chrBig\tLN:700000000\n" +
		"rBig\t0\tchrBig\t" + strconv.Itoa(bigPos) + "\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n"
	bamBytes := samToBAM(t, largeSAM)
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "big.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatal(err)
	}
	// min_shift 15 with depth 5 addresses 1<<30 bp (~1.07 Gbp), comfortably
	// past the BAI ceiling and past the 600 Mbp record.
	if err := IndexFile(bamPath, "", IndexOptions{SelectCSI: true, CSIMinShift: 15}); err != nil {
		t.Fatalf("IndexFile csi: %v", err)
	}
	idx, err := bam.ReadCSIFile(bamPath + ".csi")
	if err != nil {
		t.Fatalf("ReadCSIFile: %v", err)
	}
	if idx.MaxPos() <= 1<<29 {
		t.Errorf("CSI MaxPos %d should exceed the BAI 2^29 ceiling", idx.MaxPos())
	}
	// A region query straddling the large coordinate finds the record.
	var out bytes.Buffer
	reg := "chrBig:599999990-600000010"
	n, err := ViewFile(bamPath, &out, ViewOptions{Regions: []string{reg}, Count: true}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile %s: %v", reg, err)
	}
	if n != 1 {
		t.Errorf("%s count: got %d, want 1", reg, n)
	}
	// A region well clear of the record returns nothing.
	var out2 bytes.Buffer
	n2, err := ViewFile(bamPath, &out2, ViewOptions{Regions: []string{"chrBig:1-1000"}, Count: true}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile clear region: %v", err)
	}
	if n2 != 0 {
		t.Errorf("chrBig:1-1000 count: got %d, want 0", n2)
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
