package tabix

import (
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// writeBGZF writes the given uncompressed text to a fresh bgzipped file
// inside t.TempDir() and returns the absolute path.
func writeBGZF(t *testing.T, name, text string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	defer f.Close()
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(text)); err != nil {
		t.Fatalf("bgzip write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	return path
}

func TestPresetConfigs(t *testing.T) {
	gff, err := PresetConfig(PresetGFF)
	if err != nil {
		t.Fatalf("gff preset: %v", err)
	}
	if gff.Format != FormatGeneric || gff.ColSeq != 1 || gff.ColBeg != 4 || gff.ColEnd != 5 || gff.Meta != '#' {
		t.Errorf("gff preset wrong: %+v", gff)
	}

	bed, _ := PresetConfig(PresetBED)
	if !bed.ZeroBased() {
		t.Errorf("bed preset should be zero-based")
	}
	if bed.ColSeq != 1 || bed.ColBeg != 2 || bed.ColEnd != 3 {
		t.Errorf("bed columns wrong: %+v", bed)
	}

	sam, _ := PresetConfig(PresetSAM)
	if sam.FormatBase() != FormatSAM || sam.ColSeq != 3 || sam.ColBeg != 4 || sam.Meta != '@' {
		t.Errorf("sam preset wrong: %+v", sam)
	}

	vcf, _ := PresetConfig(PresetVCF)
	if vcf.FormatBase() != FormatVCF || vcf.ColSeq != 1 || vcf.ColBeg != 2 || vcf.ColEnd != 0 {
		t.Errorf("vcf preset wrong: %+v", vcf)
	}

	if _, err := PresetConfig("bogus"); !errors.Is(err, ErrBadPreset) {
		t.Errorf("expected ErrBadPreset, got %v", err)
	}
}

func TestBuildAndQueryVCF(t *testing.T) {
	vcf := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=chr1,length=1000000>",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tA\tT\t.\tPASS\t.",
		"chr1\t200\t.\tG\tC\t.\tPASS\t.",
		"chr1\t300\t.\tT\tA\t.\tPASS\t.",
		"chr1\t1000\t.\tC\tG\t.\tPASS\t.",
		"chr2\t100\t.\tA\tT\t.\tPASS\t.",
		"chr2\t250\t.\tA\tG\t.\tPASS\t.",
	}, "\n") + "\n"
	path := writeBGZF(t, "x.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	if len(idx.Names) != 2 {
		t.Fatalf("expected 2 chroms, got %d", len(idx.Names))
	}
	if idx.Names[0] != "chr1" || idx.Names[1] != "chr2" {
		t.Errorf("chrom order wrong: %v", idx.Names)
	}

	// Query chr1:150-250 (1-based -> our internal 0-based half-open is
	// [149, 250)). Records at POS=200 should match; 100 and 300 should not.
	got, err := idx.QueryBytes(path, "chr1", 149, 250)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("expected 1 record, got %d: %v", len(got), got)
	}
	if !bytes.Contains(got[0], []byte("\t200\t")) {
		t.Errorf("wrong record returned: %s", got[0])
	}

	// Query an empty region.
	empty, err := idx.QueryBytes(path, "chr1", 400, 500)
	if err != nil {
		t.Fatalf("empty query: %v", err)
	}
	if len(empty) != 0 {
		t.Errorf("expected empty result, got %d", len(empty))
	}

	// Query an unknown chrom.
	unk, err := idx.QueryBytes(path, "chrZ", 0, 1000)
	if err != nil {
		t.Fatalf("unknown chrom query: %v", err)
	}
	if len(unk) != 0 {
		t.Errorf("expected empty result for unknown chrom")
	}

	// Query the whole of chr2.
	chr2, err := idx.QueryBytes(path, "chr2", 0, 1000)
	if err != nil {
		t.Fatalf("chr2 query: %v", err)
	}
	if len(chr2) != 2 {
		t.Errorf("chr2 expected 2 records, got %d", len(chr2))
	}
}

func TestBuildBED(t *testing.T) {
	bed := strings.Join([]string{
		"chr1\t100\t200\tregion1",
		"chr1\t300\t400\tregion2",
		"chr2\t0\t50\tregion3",
	}, "\n") + "\n"
	path := writeBGZF(t, "x.bed.gz", bed)
	cfg, _ := PresetConfig(PresetBED)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	if !idx.Config.ZeroBased() {
		t.Errorf("BED config should be zero-based")
	}
	// Region overlapping the first record.
	got, err := idx.QueryBytes(path, "chr1", 150, 175)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("expected 1 record, got %d", len(got))
	}
	if !bytes.Contains(got[0], []byte("region1")) {
		t.Errorf("wrong record: %s", got[0])
	}
}

func TestBuildGFF(t *testing.T) {
	gff := strings.Join([]string{
		"##gff-version 3",
		"chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1",
		"chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=g2",
	}, "\n") + "\n"
	path := writeBGZF(t, "x.gff.gz", gff)
	cfg, _ := PresetConfig(PresetGFF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	got, err := idx.QueryBytes(path, "chr1", 50, 250)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("expected 1, got %d", len(got))
	}
	if !bytes.Contains(got[0], []byte("ID=g1")) {
		t.Errorf("wrong record: %s", got[0])
	}
}

func TestBuildSAM(t *testing.T) {
	// SAM-style: chrom at col 3, POS at col 4 (1-based).
	sam := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r2\t0\tchr1\t500\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
	}, "\n") + "\n"
	path := writeBGZF(t, "x.sam.gz", sam)
	cfg, _ := PresetConfig(PresetSAM)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	if len(idx.Names) != 1 || idx.Names[0] != "chr1" {
		t.Errorf("chrom list wrong: %v", idx.Names)
	}
	got, err := idx.QueryBytes(path, "chr1", 90, 110)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 || !bytes.HasPrefix(got[0], []byte("r1\t")) {
		t.Errorf("expected r1, got %q", got)
	}
}

func TestRoundTripWriteRead(t *testing.T) {
	vcf := strings.Join([]string{
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tA\tT\t.\tPASS\t.",
		"chr1\t200\t.\tG\tC\t.\tPASS\t.",
		"chr2\t100\t.\tA\tT\t.\tPASS\t.",
	}, "\n") + "\n"
	path := writeBGZF(t, "rt.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	tbiPath := path + ".tbi"
	if err := idx.WriteFile(tbiPath); err != nil {
		t.Fatalf("write tbi: %v", err)
	}
	reloaded, err := ReadFile(tbiPath)
	if err != nil {
		t.Fatalf("read tbi: %v", err)
	}
	if len(reloaded.Names) != len(idx.Names) {
		t.Fatalf("names lost in round-trip")
	}
	for i, n := range idx.Names {
		if reloaded.Names[i] != n {
			t.Errorf("name %d: got %q want %q", i, reloaded.Names[i], n)
		}
	}
	if reloaded.Config != idx.Config {
		t.Errorf("config drift: got %+v want %+v", reloaded.Config, idx.Config)
	}
	// Re-query through the reloaded index.
	got, err := reloaded.QueryBytes(path, "chr1", 150, 250)
	if err != nil {
		t.Fatalf("query after reload: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("reload query: expected 1, got %d", len(got))
	}
}

func TestQueryAcrossLinearTiles(t *testing.T) {
	// Span more than 16384 bp to populate multiple linear-index entries.
	var b strings.Builder
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for pos := 1; pos < 200000; pos += 5000 {
		b.WriteString("chr1\t")
		b.WriteString(itoa(pos))
		b.WriteString("\t.\tA\tT\t.\tPASS\t.\n")
	}
	path := writeBGZF(t, "big.vcf.gz", b.String())
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	if len(idx.Refs) != 1 {
		t.Fatalf("expected 1 ref")
	}
	if len(idx.Refs[0].Linear) < 12 {
		t.Errorf("expected at least 12 linear tiles, got %d", len(idx.Refs[0].Linear))
	}
	// Query a region in the middle.
	got, err := idx.QueryBytes(path, "chr1", 90000, 110000)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) == 0 {
		t.Errorf("expected records, got 0")
	}
}

func TestEmptyFile(t *testing.T) {
	// A bgzipped file with only header lines.
	vcf := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	path := writeBGZF(t, "empty.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build empty: %v", err)
	}
	if len(idx.Names) != 0 {
		t.Errorf("expected 0 chroms, got %d", len(idx.Names))
	}
}

func TestSingleRecord(t *testing.T) {
	vcf := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t42\t.\tA\tT\t.\tPASS\t.\n"
	path := writeBGZF(t, "one.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	got, err := idx.QueryBytes(path, "chr1", 0, 100)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("expected 1, got %d", len(got))
	}
}

func TestRecordAtChromBoundary(t *testing.T) {
	// A record at position 1 (the lowest legal 1-based POS).
	vcf := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tT\t.\tPASS\t.\n"
	path := writeBGZF(t, "boundary.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	got, err := idx.QueryBytes(path, "chr1", 0, 5)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Errorf("expected 1, got %d", len(got))
	}
}

func TestReadBadMagic(t *testing.T) {
	bad := bytes.NewReader([]byte("NOPE----"))
	if _, err := Read(bad); !errors.Is(err, ErrBadMagic) {
		t.Errorf("expected ErrBadMagic, got %v", err)
	}
}

func TestReadTruncated(t *testing.T) {
	// Magic but nothing after.
	short := bytes.NewReader([]byte("TBI\x01"))
	_, err := Read(short)
	if !errors.Is(err, ErrTruncated) {
		t.Errorf("expected ErrTruncated, got %v", err)
	}
}

func TestQueryMissingIndex(t *testing.T) {
	// Calling QueryBytes with a stale data path is harmless; calling
	// ReadFile on a non-existent path is the actual error case.
	if _, err := ReadFile("/nonexistent.tbi"); err == nil {
		t.Error("expected error for missing .tbi")
	}
}

// TestHexFixture builds a known-good tiny .tbi by hand and feeds it through
// Read to make sure the parser handles the canonical bytewise layout.
func TestHexFixture(t *testing.T) {
	// Construct the bytes: 1 ref ("c"), VCF preset, one bin (4681) with
	// one chunk, one linear entry.
	var buf bytes.Buffer
	buf.Write(Magic[:])
	wInt32 := func(v int32) { binary.Write(&buf, binary.LittleEndian, v) }
	wU32 := func(v uint32) { binary.Write(&buf, binary.LittleEndian, v) }
	wU64 := func(v uint64) { binary.Write(&buf, binary.LittleEndian, v) }
	wInt32(1)         // n_ref
	wInt32(FormatVCF) // format
	wInt32(1)         // col_seq
	wInt32(2)         // col_beg
	wInt32(0)         // col_end
	wInt32('#')       // meta
	wInt32(0)         // skip
	wInt32(2)         // l_nm
	buf.WriteString("c")
	buf.WriteByte(0)
	wInt32(1)  // n_bin
	wU32(4681) // bin id
	wInt32(1)  // n_chunk
	wU64(uint64(MakeVOffset(0, 0)))
	wU64(uint64(MakeVOffset(0, 50)))
	wInt32(1) // n_intv
	wU64(uint64(MakeVOffset(0, 0)))

	idx, err := Read(&buf)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	if len(idx.Names) != 1 || idx.Names[0] != "c" {
		t.Errorf("name wrong: %v", idx.Names)
	}
	if idx.Config.FormatBase() != FormatVCF {
		t.Errorf("format wrong: %d", idx.Config.Format)
	}
	if len(idx.Refs) != 1 || len(idx.Refs[0].Bins) != 1 || idx.Refs[0].Bins[0].ID != 4681 {
		t.Errorf("bin wrong: %+v", idx.Refs)
	}
	if len(idx.Refs[0].Linear) != 1 {
		t.Errorf("linear wrong: %d", len(idx.Refs[0].Linear))
	}
}

func TestValidateConfigRejectsBad(t *testing.T) {
	bad := Config{ColSeq: 0, ColBeg: 1}
	if _, err := Build("nonexistent", bad); !errors.Is(err, ErrBadConfig) {
		t.Errorf("expected ErrBadConfig, got %v", err)
	}
}

func TestRegionChunksMerges(t *testing.T) {
	// Construct two bins with overlapping chunks; RegionChunks should
	// merge them into one.
	idx := NewIndex(Config{ColSeq: 1, ColBeg: 2, Meta: '#'})
	idx.Names = []string{"chr1"}
	idx.Refs = []RefIndex{{
		Bins: []Bin{
			{ID: uint32(Reg2bin(0, 1)), Chunks: []Chunk{
				{Beg: MakeVOffset(0, 0), End: MakeVOffset(0, 100)},
			}},
			// Add the same chunk under bin 0 to force overlap.
			{ID: 0, Chunks: []Chunk{
				{Beg: MakeVOffset(0, 50), End: MakeVOffset(0, 200)},
			}},
		},
		Linear: []VOffset{0},
	}}
	chunks, refID, err := idx.RegionChunks("chr1", 0, 1)
	if err != nil {
		t.Fatalf("region chunks: %v", err)
	}
	if refID != 0 {
		t.Errorf("wrong refID: %d", refID)
	}
	if len(chunks) != 1 {
		t.Fatalf("expected 1 merged chunk, got %d: %v", len(chunks), chunks)
	}
	if chunks[0].End != MakeVOffset(0, 200) {
		t.Errorf("merge end wrong: %x", chunks[0].End)
	}
}

func TestReadFileBadFile(t *testing.T) {
	// Write a file that's BGZF-valid but contains garbage payload.
	path := filepath.Join(t.TempDir(), "bogus.tbi")
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	bw := bgzip.NewWriter(f)
	bw.Write([]byte("not-a-tabix"))
	bw.Close()
	f.Close()
	if _, err := ReadFile(path); !errors.Is(err, ErrBadMagic) {
		t.Errorf("expected ErrBadMagic, got %v", err)
	}
}

func TestChromsCopy(t *testing.T) {
	idx := NewIndex(Config{ColSeq: 1, ColBeg: 2})
	idx.Names = []string{"a", "b"}
	got := idx.Chroms()
	got[0] = "mutated"
	if idx.Names[0] != "a" {
		t.Errorf("Chroms() should return a copy; got mutation through it")
	}
}

func TestSplitNullStrings(t *testing.T) {
	out := splitNullStrings([]byte("a\x00bb\x00ccc\x00"))
	if len(out) != 3 || out[0] != "a" || out[1] != "bb" || out[2] != "ccc" {
		t.Errorf("split wrong: %v", out)
	}
}

// itoa is a tiny strconv-free integer formatter for test fixture
// construction. Keeps the tests independent of the strconv package's exact
// behaviour to make TestQueryAcrossLinearTiles less fragile.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	negative := n < 0
	if negative {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if negative {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

// Ensure imports used only in some paths still compile.
var _ = io.EOF
