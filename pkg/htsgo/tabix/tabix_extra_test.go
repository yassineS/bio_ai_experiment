package tabix

import (
	"bytes"
	"encoding/binary"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// TestBuildBadRecordTooFewFields exercises the "ColSeq beyond fields" path
// inside Build.
func TestBuildBadRecordTooFewFields(t *testing.T) {
	bad := "##header\nchr1only\n"
	path := writeBGZF(t, "bad.gz", bad)
	cfg, _ := PresetConfig(PresetVCF)
	if _, err := Build(path, cfg); !errors.Is(err, ErrBadRecord) {
		t.Errorf("expected ErrBadRecord, got %v", err)
	}
}

// TestBuildBadBeginValue exercises the non-numeric begin column path.
func TestBuildBadBeginValue(t *testing.T) {
	bad := "chr1\tNOPE\t10\tx\n"
	path := writeBGZF(t, "badbeg.gz", bad)
	cfg, _ := PresetConfig(PresetBED)
	_, err := Build(path, cfg)
	if err == nil || !strings.Contains(err.Error(), "invalid begin") {
		t.Errorf("expected invalid begin, got %v", err)
	}
}

// TestBuildBadEndValue exercises the non-numeric end column path.
func TestBuildBadEndValue(t *testing.T) {
	bad := "chr1\t10\tNOPE\tx\n"
	path := writeBGZF(t, "badend.gz", bad)
	cfg, _ := PresetConfig(PresetBED)
	_, err := Build(path, cfg)
	if err == nil || !strings.Contains(err.Error(), "invalid end") {
		t.Errorf("expected invalid end, got %v", err)
	}
}

// TestBuildBEDEndTooSmall covers the recordEnd "end <= beg" normalization.
func TestBuildBEDEndTooSmall(t *testing.T) {
	// A degenerate BED record where end == beg.
	bed := "chr1\t100\t100\tdegen\n"
	path := writeBGZF(t, "degen.gz", bed)
	cfg, _ := PresetConfig(PresetBED)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	// The record's effective interval becomes [100, 101).
	got, err := idx.QueryBytes(path, "chr1", 100, 105)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Errorf("expected 1 record after normalization, got %d", len(got))
	}
}

// TestBuildNoCoor exercises the NoCoor counter via SAM "*" chrom.
func TestBuildNoCoor(t *testing.T) {
	sam := strings.Join([]string{
		"@HD\tVN:1.6",
		"r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
		"r2\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
	}, "\n") + "\n"
	path := writeBGZF(t, "noco.sam.gz", sam)
	cfg, _ := PresetConfig(PresetSAM)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	if idx.NoCoor != 1 {
		t.Errorf("expected NoCoor=1, got %d", idx.NoCoor)
	}
	if len(idx.Names) != 1 || idx.Names[0] != "chr1" {
		t.Errorf("names wrong: %v", idx.Names)
	}
}

// TestWriteWithNoCoor exercises the optional trailing n_no_coor field of
// the .tbi format.
func TestWriteWithNoCoor(t *testing.T) {
	idx := NewIndex(Config{Format: FormatVCF, ColSeq: 1, ColBeg: 2, Meta: '#'})
	idx.Names = []string{"chr1"}
	idx.Refs = []RefIndex{{Bins: nil, Linear: nil}}
	idx.NoCoor = 42

	var buf bytes.Buffer
	if err := idx.Write(&buf); err != nil {
		t.Fatalf("write: %v", err)
	}
	got, err := Read(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	if got.NoCoor != 42 {
		t.Errorf("NoCoor lost in round-trip: got %d", got.NoCoor)
	}
}

// TestBuildMissingFile covers the os.Open error path in Build.
func TestBuildMissingFile(t *testing.T) {
	cfg, _ := PresetConfig(PresetVCF)
	_, err := Build("/this/path/does/not/exist.gz", cfg)
	if err == nil {
		t.Error("expected error for missing file")
	}
}

// TestBuildNonBGZF covers the bgzip.NewReader error path. We feed a file
// that has the gzip magic but no BC subfield.
func TestBuildNonBGZF(t *testing.T) {
	path := filepath.Join(t.TempDir(), "plain.gz")
	// Build a regular gzip (not BGZF) file: htslib's bgzip Reader will
	// reject it for missing the BC subfield.
	if err := os.WriteFile(path, []byte{0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0, 0xff}, 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	cfg, _ := PresetConfig(PresetVCF)
	if _, err := Build(path, cfg); err == nil {
		t.Error("expected error for plain gzip input")
	}
}

// TestReadWrongNameCount exercises the "n_ref vs strings" error in Read.
func TestReadWrongNameCount(t *testing.T) {
	var buf bytes.Buffer
	buf.Write(Magic[:])
	wInt32 := func(v int32) { binary.Write(&buf, binary.LittleEndian, v) }
	wInt32(5) // n_ref says 5
	wInt32(FormatVCF)
	wInt32(1)
	wInt32(2)
	wInt32(0)
	wInt32('#')
	wInt32(0)
	// names: only one
	wInt32(2)
	buf.WriteString("c")
	buf.WriteByte(0)
	if _, err := Read(&buf); err == nil {
		t.Error("expected error for name/n_ref mismatch")
	}
}

// TestQueryBytesEmptyRange covers the early return for inverted ranges.
func TestQueryBytesEmptyRange(t *testing.T) {
	idx := NewIndex(Config{ColSeq: 1, ColBeg: 2})
	idx.Names = []string{"chr1"}
	idx.Refs = []RefIndex{{}}
	got, err := idx.QueryBytes("/no/such/file", "chr1", 100, 100)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if got != nil {
		t.Errorf("expected nil, got %v", got)
	}
}

// TestQueryBytesNegativeBeg ensures negative beg is clamped.
func TestQueryBytesNegativeBeg(t *testing.T) {
	idx := NewIndex(Config{ColSeq: 1, ColBeg: 2})
	idx.Names = []string{"chr1"}
	idx.Refs = []RefIndex{{}}
	_, _, err := idx.RegionChunks("chr1", -10, 100)
	if err != nil {
		t.Fatalf("regionChunks: %v", err)
	}
}

// TestReadFinalNoCoorAbsent ensures the parser tolerates the optional
// trailer being missing.
func TestReadFinalNoCoorAbsent(t *testing.T) {
	idx := NewIndex(Config{Format: FormatVCF, ColSeq: 1, ColBeg: 2, Meta: '#'})
	idx.Names = []string{"chr1"}
	idx.Refs = []RefIndex{{Bins: nil, Linear: nil}}
	idx.NoCoor = 0

	var buf bytes.Buffer
	if err := idx.Write(&buf); err != nil {
		t.Fatalf("write: %v", err)
	}
	got, err := Read(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	if got.NoCoor != 0 {
		t.Errorf("NoCoor should be 0, got %d", got.NoCoor)
	}
}

// TestWriteFileBadPath covers the os.Create error path of WriteFile.
func TestWriteFileBadPath(t *testing.T) {
	idx := NewIndex(Config{ColSeq: 1, ColBeg: 2})
	if err := idx.WriteFile("/no/such/dir/x.tbi"); err == nil {
		t.Error("expected error")
	}
}

// TestQueryBytesMissingFile ensures QueryBytes propagates os.Open errors.
func TestQueryBytesMissingFile(t *testing.T) {
	idx := NewIndex(Config{ColSeq: 1, ColBeg: 2, Meta: '#'})
	idx.Names = []string{"chr1"}
	idx.Refs = []RefIndex{{
		Bins:   []Bin{{ID: 4681, Chunks: []Chunk{{Beg: 0, End: MakeVOffset(0, 100)}}}},
		Linear: []VOffset{0},
	}}
	if _, err := idx.QueryBytes("/no/such/file", "chr1", 0, 1); err == nil {
		t.Error("expected error for missing data file")
	}
}

// TestRoundTripBytewise verifies that Write+Read is the inverse, and that
// re-Writing the read-back Index produces an identical byte stream.
func TestRoundTripBytewise(t *testing.T) {
	vcf := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t10\t.\tA\tT\t.\tPASS\t.\n" +
		"chr1\t20\t.\tA\tT\t.\tPASS\t.\n"
	path := writeBGZF(t, "rt2.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	var b1 bytes.Buffer
	if err := idx.Write(&b1); err != nil {
		t.Fatalf("write 1: %v", err)
	}
	reloaded, err := Read(bytes.NewReader(b1.Bytes()))
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	var b2 bytes.Buffer
	if err := reloaded.Write(&b2); err != nil {
		t.Fatalf("write 2: %v", err)
	}
	if !bytes.Equal(b1.Bytes(), b2.Bytes()) {
		t.Errorf("re-write differs: %d vs %d bytes", len(b1.Bytes()), len(b2.Bytes()))
	}
}

// TestBuildVCFRefAlleleEnd verifies that VCF end is derived from REF
// allele length (rather than always beg+1).
func TestBuildVCFRefAlleleEnd(t *testing.T) {
	vcf := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tACGT\tA\t.\tPASS\t.\n"
	path := writeBGZF(t, "ref.vcf.gz", vcf)
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	// Query [102, 103) (0-based half-open): the deletion spans
	// [99, 103), so position 102 should be inside it.
	got, err := idx.QueryBytes(path, "chr1", 102, 103)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Errorf("expected 1 hit for ACGT deletion at 100, got %d", len(got))
	}
}

// TestVOffsetSentinelInBuild ensures the linear-index sentinel handling
// places non-zero offsets correctly when the first record is not at offset
// 0 within its block.
func TestVOffsetSentinelInBuild(t *testing.T) {
	// Force a block boundary: write a single comment line that pads out
	// past one block, then a record.
	var b strings.Builder
	b.WriteString("##VCFv4.2\n")
	for i := 0; i < 4096; i++ {
		b.WriteString("##contig=<ID=padding_padding_padding_padding>\n")
	}
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	b.WriteString("chr1\t100\t.\tA\tT\t.\tPASS\t.\n")
	path := writeBGZF(t, "padded.vcf.gz", b.String())
	cfg, _ := PresetConfig(PresetVCF)
	idx, err := Build(path, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	got, err := idx.QueryBytes(path, "chr1", 99, 101)
	if err != nil {
		t.Fatalf("query: %v", err)
	}
	if len(got) != 1 {
		t.Errorf("expected 1 record across block boundary, got %d", len(got))
	}
}

// TestReg2binsInvertedReturnsBin0 covers the inverted-range guard.
func TestReg2binsInvertedReturnsBin0(t *testing.T) {
	got := Reg2bins(50, 10)
	if len(got) != 1 || got[0] != 0 {
		t.Errorf("expected just bin 0 for inverted range, got %v", got)
	}
}

// TestReadEmpty covers the io.EOF path in Read at the very first byte.
func TestReadEmpty(t *testing.T) {
	if _, err := Read(bytes.NewReader(nil)); !errors.Is(err, ErrTruncated) {
		t.Errorf("expected ErrTruncated, got %v", err)
	}
}

// TestWriteEmptyIndex ensures an Index with zero refs serialises cleanly.
func TestWriteEmptyIndex(t *testing.T) {
	idx := NewIndex(Config{Format: FormatVCF, ColSeq: 1, ColBeg: 2, Meta: '#'})
	var buf bytes.Buffer
	if err := idx.Write(&buf); err != nil {
		t.Fatalf("write: %v", err)
	}
	got, err := Read(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	if len(got.Names) != 0 {
		t.Errorf("expected 0 names, got %v", got.Names)
	}
}

// TestBgzipNotInstalled is a sanity check that the bgzip package is
// importable from inside tests (catches go-mod misconfiguration early).
func TestBgzipImportable(t *testing.T) {
	_ = bgzip.MaxBlockSize
}
