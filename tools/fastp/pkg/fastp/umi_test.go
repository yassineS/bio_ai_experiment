package fastp

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestValidUMILocation(t *testing.T) {
	good := []string{UMILocRead1, UMILocRead2, UMILocPerRead, UMILocIndex1, UMILocIndex2, UMILocPerIndex}
	for _, s := range good {
		if !ValidUMILocation(s) {
			t.Errorf("%q should be valid", s)
		}
	}
	for _, s := range []string{"", "READ1", "i7", "foo"} {
		if ValidUMILocation(s) {
			t.Errorf("%q should be invalid", s)
		}
	}
}

func TestApplyUMIRead1(t *testing.T) {
	rec := &fastq.Record{
		ID:          "r1",
		Description: "r1",
		Sequence:    []byte("ATCGATCGAAAACCCCGGGG"),
		Quality:     []byte("IIIIIIIIIIIIIIIIIIII"),
	}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocRead1
	opts.UMILen = 6

	out, _ := applyUMI(rec, nil, opts, stats)
	if got := out.ID; got != "r1:UMI_ATCGAT" {
		t.Fatalf("id: want %q, got %q", "r1:UMI_ATCGAT", got)
	}
	if got := string(out.Sequence); got != "CGAAAACCCCGGGG" {
		t.Fatalf("seq after UMI removal: want CGAAAACCCCGGGG, got %s", got)
	}
	if stats.UMIProcessed != 1 {
		t.Fatalf("UMIProcessed: want 1, got %d", stats.UMIProcessed)
	}
}

func TestApplyUMIRead1WithSkipAndPrefix(t *testing.T) {
	rec := &fastq.Record{
		ID:       "r1",
		Sequence: []byte("AAAAAATTGGGGCCCC"),
		Quality:  []byte("IIIIIIIIIIIIIIII"),
	}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocRead1
	opts.UMILen = 6 // AAAAAA
	opts.UMISkip = 2
	opts.UMIPrefix = "X"

	out, _ := applyUMI(rec, nil, opts, stats)
	if got := out.ID; got != "r1:UMI_XAAAAAA" {
		t.Fatalf("id: want r1:UMI_XAAAAAA, got %s", got)
	}
	// 6 UMI bases + 2 skipped bases removed from the 16-byte seq.
	if got := string(out.Sequence); got != "GGGGCCCC" {
		t.Fatalf("seq: want GGGGCCCC, got %s", got)
	}
	if got := string(out.Quality); got != "IIIIIIII" {
		t.Fatalf("qual: want IIIIIIII, got %s", got)
	}
}

func TestApplyUMITooShort(t *testing.T) {
	rec := &fastq.Record{
		ID:       "r1",
		Sequence: []byte("ACG"),
		Quality:  []byte("III"),
	}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocRead1
	opts.UMILen = 6

	out, _ := applyUMI(rec, nil, opts, stats)
	if out.ID != "r1" {
		t.Fatalf("short read should be untouched, got %s", out.ID)
	}
	if stats.UMIProcessed != 0 {
		t.Fatalf("short read should not increment counter, got %d", stats.UMIProcessed)
	}
}

func TestApplyUMIRead2(t *testing.T) {
	r1 := &fastq.Record{ID: "r1", Sequence: []byte("AAAAAAAAAA"), Quality: []byte("IIIIIIIIII")}
	r2 := &fastq.Record{ID: "r2", Sequence: []byte("CCCCCCCCCC"), Quality: []byte("IIIIIIIIII")}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocRead2
	opts.UMILen = 4

	r1Out, r2Out := applyUMI(r1, r2, opts, stats)
	if r1Out.ID != "r1" {
		t.Fatalf("r1 should be untouched: %s", r1Out.ID)
	}
	if r2Out.ID != "r2:UMI_CCCC" {
		t.Fatalf("r2 id: want r2:UMI_CCCC, got %s", r2Out.ID)
	}
	if string(r2Out.Sequence) != "CCCCCC" {
		t.Fatalf("r2 seq after UMI removal: want CCCCCC, got %s", string(r2Out.Sequence))
	}
	if stats.UMIProcessed != 1 {
		t.Fatalf("UMIProcessed: want 1, got %d", stats.UMIProcessed)
	}
}

func TestApplyUMIPerRead(t *testing.T) {
	r1 := &fastq.Record{ID: "r1", Sequence: []byte("AAAAAAAAAA"), Quality: []byte("IIIIIIIIII")}
	r2 := &fastq.Record{ID: "r2", Sequence: []byte("CCCCCCCCCC"), Quality: []byte("IIIIIIIIII")}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocPerRead
	opts.UMILen = 4

	r1Out, r2Out := applyUMI(r1, r2, opts, stats)
	if r1Out.ID != "r1:UMI_AAAA_CCCC" {
		t.Fatalf("r1 id: want r1:UMI_AAAA_CCCC, got %s", r1Out.ID)
	}
	if r2Out.ID != "r2:UMI_AAAA_CCCC" {
		t.Fatalf("r2 id: want r2:UMI_AAAA_CCCC, got %s", r2Out.ID)
	}
	if stats.UMIProcessed != 2 {
		t.Fatalf("UMIProcessed: want 2, got %d", stats.UMIProcessed)
	}
}

func TestApplyUMIIndex1(t *testing.T) {
	// Illumina-style header: ID + space + "1:N:0:ATCACG+CGATGT"
	r1 := &fastq.Record{
		ID:          "M00001:1:000",
		Description: "M00001:1:000 1:N:0:ATCACG+CGATGT",
		Sequence:    []byte("ACGTACGTAC"),
		Quality:     []byte("IIIIIIIIII"),
	}
	r2 := &fastq.Record{
		ID:          "M00001:1:000",
		Description: "M00001:1:000 2:N:0:ATCACG+CGATGT",
		Sequence:    []byte("ACGTACGTAC"),
		Quality:     []byte("IIIIIIIIII"),
	}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocIndex1

	r1Out, r2Out := applyUMI(r1, r2, opts, stats)
	if !strings.HasSuffix(r1Out.ID, ":UMI_ATCACG") {
		t.Fatalf("r1 id: want suffix :UMI_ATCACG, got %s", r1Out.ID)
	}
	if !strings.HasSuffix(r2Out.ID, ":UMI_ATCACG") {
		t.Fatalf("r2 id: want suffix :UMI_ATCACG, got %s", r2Out.ID)
	}
	// Sequence is untouched in index modes.
	if string(r1Out.Sequence) != "ACGTACGTAC" {
		t.Fatalf("seq should be untouched, got %s", string(r1Out.Sequence))
	}
	if stats.UMIProcessed != 2 {
		t.Fatalf("UMIProcessed: want 2 (PE), got %d", stats.UMIProcessed)
	}
}

func TestApplyUMIIndex2AndPerIndex(t *testing.T) {
	r1 := &fastq.Record{
		ID:          "M00001",
		Description: "M00001 1:N:0:ATCACG+CGATGT",
		Sequence:    []byte("ACGT"),
		Quality:     []byte("IIII"),
	}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = UMILocIndex2
	r1Out, _ := applyUMI(r1, nil, opts, stats)
	if !strings.HasSuffix(r1Out.ID, ":UMI_CGATGT") {
		t.Fatalf("index2 id: want suffix :UMI_CGATGT, got %s", r1Out.ID)
	}

	opts.UMILoc = UMILocPerIndex
	stats = &ProcessStats{}
	r1Out, _ = applyUMI(r1, nil, opts, stats)
	if !strings.HasSuffix(r1Out.ID, ":UMI_ATCACG_CGATGT") {
		t.Fatalf("per_index id: want suffix :UMI_ATCACG_CGATGT, got %s", r1Out.ID)
	}
}

func TestApplyUMIDefaultLoc(t *testing.T) {
	// No UMILoc set: SE defaults to read1, PE defaults to per_read.
	r := &fastq.Record{ID: "x", Sequence: []byte("AAAACCCCGG"), Quality: []byte("IIIIIIIIII")}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILen = 4
	out, _ := applyUMI(r, nil, opts, stats)
	if !strings.HasSuffix(out.ID, ":UMI_AAAA") {
		t.Fatalf("SE default should fall back to read1; got %s", out.ID)
	}
}

func TestProcessSingleEndUMIEndToEnd(t *testing.T) {
	input := strings.Join([]string{
		"@r1", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII",
		"@r2", "GGGGCCCCAAAATTTTAAAA", "+", "IIIIIIIIIIIIIIIIIIII",
		"",
	}, "\n")
	var out bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MinLength = 1
	opts.LengthRequired = 1
	opts.QualThreshold = 0
	opts.UMI = true
	opts.UMILoc = UMILocRead1
	opts.UMILen = 6

	stats, err := ProcessSingleEnd(strings.NewReader(input), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	if stats.UMIProcessed != 2 {
		t.Fatalf("UMIProcessed: want 2, got %d", stats.UMIProcessed)
	}
	written := out.String()
	if !strings.Contains(written, ":UMI_ACGTAC") {
		t.Fatalf("output should contain UMI tag, got:\n%s", written)
	}
	if !strings.Contains(written, ":UMI_GGGGCC") {
		t.Fatalf("output should contain UMI tag for second read, got:\n%s", written)
	}
}

func TestLegacyUMIFlagsStillWork(t *testing.T) {
	// Legacy callers that only set UMILength/UMILocation/UMISkip should
	// still produce the same name suffix.
	r := &fastq.Record{ID: "rid", Sequence: []byte("ACGTACGTACGT"), Quality: []byte("IIIIIIIIIIII")}
	stats := &ProcessStats{}
	opts := DefaultProcessOptions()
	opts.UMILength = 4
	opts.UMILocation = UMILocRead1

	out, _ := extractUMI(r, nil, opts, stats)
	if !strings.HasSuffix(out.ID, ":UMI_ACGT") {
		t.Fatalf("legacy umi-length flag: want suffix :UMI_ACGT, got %s", out.ID)
	}
}
