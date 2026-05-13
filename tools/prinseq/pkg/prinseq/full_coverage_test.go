package prinseq

import (
	"bufio"
	"bytes"
	"errors"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// bufioScannerFor wraps a string in a bufio.Scanner identical to the one used
// internally by the package.
func bufioScannerFor(s string) *bufio.Scanner {
	sc := bufio.NewScanner(strings.NewReader(s))
	sc.Buffer(make([]byte, 0, 1024*1024), 1024*1024)
	return sc
}

// --- Helpers --------------------------------------------------------------

// errReader is an io.Reader that returns a synthetic error on the first read.
type errReader struct{}

func (errReader) Read(p []byte) (int, error) { return 0, errors.New("synthetic read error") }

// errWriter is an io.Writer that errors on every Write.
type errWriter struct{ failAfter int }

func (e *errWriter) Write(p []byte) (int, error) {
	if e.failAfter <= 0 {
		return 0, errors.New("synthetic write error")
	}
	if len(p) > e.failAfter {
		n := e.failAfter
		e.failAfter = 0
		return n, errors.New("synthetic write error mid-stream")
	}
	e.failAfter -= len(p)
	return len(p), nil
}

// --- prinseq.go: empty-line skipping in FASTA parsers --------------------

func TestFastaStatsSkipsBlankLines(t *testing.T) {
	// blank lines between records exercise the `if len(line) == 0 { continue }`
	// branch in both calculateFastaStats and calculateEnhancedFastaStats.
	input := ">a\n\nACGT\n\n>b\n\nNNAA\n\n"
	stats, err := CalculateStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateStats: %v", err)
	}
	if stats.NumReads != 2 {
		t.Errorf("NumReads = %d, want 2", stats.NumReads)
	}
	if stats.TotalBases != 8 {
		t.Errorf("TotalBases = %d, want 8", stats.TotalBases)
	}

	es, err := CalculateEnhancedStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats: %v", err)
	}
	if es.NumReads != 2 {
		t.Errorf("enhanced NumReads = %d, want 2", es.NumReads)
	}
}

// --- prinseq.go: processEnhancedSequence N-counting branch ---------------

func TestEnhancedStatsCountsNs(t *testing.T) {
	// processEnhancedSequence's `case 'N','n'` is only hit when the input has Ns.
	input := ">a\nACGTNNNN\n"
	es, err := CalculateEnhancedStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if es.NumNs != 4 {
		t.Errorf("enhanced NumNs = %d, want 4", es.NumNs)
	}
}

// --- prinseq.go: scanner.Err propagation ---------------------------------

func TestCalculateStatsScannerError(t *testing.T) {
	// A reader that yields some valid bytes then errors mid-stream so that
	// bufio.Scanner returns the error from scanner.Err(). The default
	// scanner buffer is 1MB so an unterminated huge line triggers
	// bufio.ErrTooLong.
	huge := make([]byte, 2*1024*1024)
	huge[0] = '>'
	for i := 1; i < len(huge); i++ {
		huge[i] = 'A'
	}
	// no newline at the end => one giant token > buffer => scanner.Err = ErrTooLong
	_, err := CalculateStats(bytes.NewReader(huge), false)
	if err == nil {
		t.Errorf("expected scanner error from huge FASTA line")
	}

	// Enhanced variant.
	_, err = CalculateEnhancedStats(bytes.NewReader(huge), false)
	if err == nil {
		t.Errorf("expected enhanced scanner error from huge FASTA line")
	}

	// FASTQ scanner error path.
	hugeQ := make([]byte, 2*1024*1024)
	hugeQ[0] = '@'
	for i := 1; i < len(hugeQ); i++ {
		hugeQ[i] = 'A'
	}
	if _, err := CalculateStats(bytes.NewReader(hugeQ), true); err == nil {
		t.Errorf("expected scanner error from huge FASTQ line")
	}
	if _, err := CalculateEnhancedStats(bytes.NewReader(hugeQ), true); err == nil {
		t.Errorf("expected enhanced scanner error from huge FASTQ line")
	}
}

// --- prinseq.go: trimSequence with quality strings -----------------------

func TestTrimSequenceWithQuality(t *testing.T) {
	// All four trim modes (TrimLeft, TrimRight, TrimLeftP, TrimRightP) take
	// a separate qual branch (`if len(qualBytes) > 0`). Provide explicit
	// quality strings to cover those branches.
	tests := []struct {
		name     string
		seq      string
		qual     string
		opts     FilterOptions
		wantSeq  string
		wantQual string
	}{
		{
			name:     "fixed left with qual",
			seq:      "ACGTACGT",
			qual:     "ABCDEFGH",
			opts:     FilterOptions{TrimLeft: 3},
			wantSeq:  "TACGT",
			wantQual: "DEFGH",
		},
		{
			name:     "fixed right with qual",
			seq:      "ACGTACGT",
			qual:     "ABCDEFGH",
			opts:     FilterOptions{TrimRight: 3},
			wantSeq:  "ACGTA",
			wantQual: "ABCDE",
		},
		{
			name:     "percent left 25 with qual",
			seq:      "ACGTACGT",
			qual:     "ABCDEFGH",
			opts:     FilterOptions{TrimLeftP: 25}, // 8*25/100=2
			wantSeq:  "GTACGT",
			wantQual: "CDEFGH",
		},
		{
			name:     "percent right 25 with qual",
			seq:      "ACGTACGT",
			qual:     "ABCDEFGH",
			opts:     FilterOptions{TrimRightP: 25}, // 8*25/100=2
			wantSeq:  "ACGTAC",
			wantQual: "ABCDEF",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotSeq, gotQual := trimSequence(tt.seq, tt.qual, tt.opts)
			if gotSeq != tt.wantSeq {
				t.Errorf("seq = %q, want %q", gotSeq, tt.wantSeq)
			}
			if gotQual != tt.wantQual {
				t.Errorf("qual = %q, want %q", gotQual, tt.wantQual)
			}
		})
	}
}

// --- prinseq.go: trimPoly* with quality strings (the `if len(qual) > 0`
// branches inside trimPolyNRight/trimPolyATLeft/trimPolyATRight) ---------

func TestTrimPolyWithQuality(t *testing.T) {
	tests := []struct {
		name     string
		seq      string
		qual     string
		opts     FilterOptions
		wantSeq  string
		wantQual string
	}{
		{
			name:     "polyN right with qual",
			seq:      "ACGTNN",
			qual:     "IIIIII",
			opts:     FilterOptions{TrimNsRight: 2},
			wantSeq:  "ACGT",
			wantQual: "IIII",
		},
		{
			name:     "polyAT left with qual",
			seq:      "AATCG",
			qual:     "IIIII",
			opts:     FilterOptions{TrimTailLeft: 2},
			wantSeq:  "CG",
			wantQual: "II",
		},
		{
			name:     "polyAT right with qual",
			seq:      "GCAAA",
			qual:     "IIIII",
			opts:     FilterOptions{TrimTailRight: 2},
			wantSeq:  "GC",
			wantQual: "II",
		},
		{
			name:     "polyN left with qual",
			seq:      "NNACGT",
			qual:     "IIIIII",
			opts:     FilterOptions{TrimNsLeft: 2},
			wantSeq:  "ACGT",
			wantQual: "IIII",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotSeq, gotQual := trimSequence(tt.seq, tt.qual, tt.opts)
			if gotSeq != tt.wantSeq {
				t.Errorf("seq=%q want %q", gotSeq, tt.wantSeq)
			}
			if gotQual != tt.wantQual {
				t.Errorf("qual=%q want %q", gotQual, tt.wantQual)
			}
		})
	}
}

// --- prinseq.go: calculateFastqStats wrapper -----------------------------

func TestCalculateFastqStatsWrapperDirect(t *testing.T) {
	// calculateFastqStats is a thin wrapper that hard-codes offset=33.
	// It's not reachable via the public API (which always goes through
	// calculateFastqStatsWithOffset), so we exercise it directly to make
	// the dead-code cost explicit.
	scanner := bufioScannerFor("@a\nACGT\n+\nIIII\n")
	stats := &Stats{MinLength: -1}
	out, err := calculateFastqStats(scanner, stats)
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if out.NumReads != 1 {
		t.Errorf("NumReads = %d, want 1", out.NumReads)
	}
	if d := out.AvgQuality - 40.0; d > 0.001 || d < -0.001 {
		t.Errorf("AvgQuality = %v, want ~40 (proves Phred+33 offset)", out.AvgQuality)
	}
}

// --- prinseq.go: trimSequence corner cases for percentage trimming -------

func TestTrimSequencePercentageEdgeCases(t *testing.T) {
	// Percentage that rounds down to 0 should NOT trim (inner if-statement skips).
	gotSeq, _ := trimSequence("ACGTACGT", "", FilterOptions{TrimLeftP: 5}) // 8*5/100=0
	if gotSeq != "ACGTACGT" {
		t.Errorf("0-percent trim should be no-op, got %q", gotSeq)
	}

	// Percentage that rounds up to the entire length should also NOT trim.
	gotSeq, _ = trimSequence("ACGT", "", FilterOptions{TrimLeftP: 100})
	if gotSeq != "ACGT" {
		t.Errorf("100%% left trim should be no-op (guarded), got %q", gotSeq)
	}

	gotSeq, _ = trimSequence("ACGTACGT", "", FilterOptions{TrimRightP: 5})
	if gotSeq != "ACGTACGT" {
		t.Errorf("0-percent right trim should be no-op, got %q", gotSeq)
	}

	gotSeq, _ = trimSequence("ACGT", "", FilterOptions{TrimRightP: 100})
	if gotSeq != "ACGT" {
		t.Errorf("100%% right trim should be no-op (guarded), got %q", gotSeq)
	}

	// Fixed left/right trim that exceeds sequence length is a no-op.
	gotSeq, _ = trimSequence("ACGT", "", FilterOptions{TrimLeft: 10})
	if gotSeq != "ACGT" {
		t.Errorf("over-length TrimLeft should be no-op, got %q", gotSeq)
	}
	gotSeq, _ = trimSequence("ACGT", "", FilterOptions{TrimRight: 10})
	if gotSeq != "ACGT" {
		t.Errorf("over-length TrimRight should be no-op, got %q", gotSeq)
	}
}

// --- prinseq.go: trimPolyN/AT below threshold no-op branches -------------

func TestTrimPolyBelowThreshold(t *testing.T) {
	// trimPolyNRight: 1 trailing N but minLen=3 -> no-op.
	gotSeq, _ := trimSequence("ACGTN", "", FilterOptions{TrimNsRight: 3})
	if gotSeq != "ACGTN" {
		t.Errorf("polyN right below threshold should be no-op, got %q", gotSeq)
	}
	// trimPolyATRight: 1 trailing A but minLen=3 -> no-op.
	gotSeq, _ = trimSequence("CCCCA", "", FilterOptions{TrimTailRight: 3})
	if gotSeq != "CCCCA" {
		t.Errorf("polyAT right below threshold should be no-op, got %q", gotSeq)
	}
	// trimPolyATLeft: 1 leading T but minLen=3 -> no-op.
	gotSeq, _ = trimSequence("TCCCC", "", FilterOptions{TrimTailLeft: 3})
	if gotSeq != "TCCCC" {
		t.Errorf("polyAT left below threshold should be no-op, got %q", gotSeq)
	}
}

// --- prinseq.go: calculateAvgQualityScoreWithOffset empty quality --------

func TestCalculateAvgQualityScoreEmpty(t *testing.T) {
	if v := calculateAvgQualityScore(""); v != 0.0 {
		t.Errorf("empty qual avg = %v, want 0", v)
	}
}

// --- prinseq.go: getQualityEncoding default branch -----------------------

func TestGetQualityEncodingDefault(t *testing.T) {
	// sanger -> Phred33, anything else (incl. empty) -> Phred33,
	// illumina -> Phred64. The default-branch return covers both
	// "" and "sanger" inputs.
	if enc := getQualityEncoding("sanger"); enc != getQualityEncoding("") {
		t.Errorf("sanger and empty should produce the same encoding")
	}
	if enc := getQualityEncoding("illumina"); enc == getQualityEncoding("") {
		t.Errorf("illumina should differ from default")
	}
}

// --- prinseq.go: shouldFilterDuplicate Derep==0 branch -------------------

func TestShouldFilterDuplicateDerepZero(t *testing.T) {
	seen := make(map[string]int)
	if shouldFilterDuplicate("ACGT", seen, FilterOptions{Derep: 0}) {
		t.Errorf("with Derep=0 nothing is a duplicate")
	}
}

// --- prinseq.go: reverseComplement of base outside complement table ------

func TestReverseComplementUnknownBase(t *testing.T) {
	// 'X' has no entry in the complement map -> falls through to the
	// else branch and is copied verbatim (in reverse position).
	got := reverseComplement("ACGTX")
	want := "XACGT"
	if got != want {
		t.Errorf("reverseComplement(ACGTX) = %q, want %q", got, want)
	}
	// Also test ambiguity code 'R' for the same branch.
	got = reverseComplement("RR")
	if got != "RR" {
		t.Errorf("reverseComplement(RR) = %q, want %q", got, "RR")
	}
}

// --- prinseq.go: derep mode 4 actually triggers a drop -------------------

func TestDerepRevcompTriggersDrop(t *testing.T) {
	// Sequence ACGTA stores revcomp TACGT. Sequence TACGT stores revcomp ACGTA.
	// Neither collides. Need a sequence that is its own revcomp ("ACGT" -> "ACGT")
	// so the second occurrence's revcomp counter reaches the threshold.
	input := ">a\nACGT\n>b\nACGT\n"
	var out bytes.Buffer
	if err := Filter(strings.NewReader(input), &out, false, FilterOptions{Derep: 4, DerepMin: 2}); err != nil {
		t.Fatalf("err: %v", err)
	}
	// Both a and b store the same revcomp ACGT -> 'b' should be dropped.
	if n := strings.Count(out.String(), ">"); n != 1 {
		t.Errorf("derep revcomp kept %d records, want 1", n)
	}
}

// --- prinseq.go: filter*/filterPaired* writer errors ----------------------

// hugeFastaRecord produces a FASTA record large enough to overflow the
// internal bufio buffer used by fasta.Writer (default 4 KiB). Returning the
// record forces a real Write to the underlying writer during writer.Write(),
// which is where the inline "if err := writer.Write(...); err != nil" paths
// live (rather than only at Flush time).
func hugeFastaRecord() string {
	return ">big\n" + strings.Repeat("ACGT", 4096) + "\n"
}

func hugeFastqRecord() string {
	seq := strings.Repeat("ACGT", 4096)
	qual := strings.Repeat("I", len(seq))
	return "@big\n" + seq + "\n+\n" + qual + "\n"
}

func TestFilterFastaWriterError(t *testing.T) {
	err := Filter(strings.NewReader(hugeFastaRecord()), &errWriter{}, false, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from filter when writer fails")
	}
	// Also cover the Flush() error path with a small record.
	err = Filter(strings.NewReader(">a\nACGT\n"), &errWriter{}, false, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from filter Flush when writer fails")
	}
}

func TestFilterFastqWriterError(t *testing.T) {
	err := Filter(strings.NewReader(hugeFastqRecord()), &errWriter{}, true, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from fastq filter when writer fails")
	}
	err = Filter(strings.NewReader("@a\nACGT\n+\nIIII\n"), &errWriter{}, true, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from fastq Flush when writer fails")
	}
}

func TestFilterFastaOutBadWriterError(t *testing.T) {
	// Big "drop" record so badWriter.Write actually flushes to errWriter.
	input := ">drop\n" + strings.Repeat("A", 16384) + "\n"
	err := Filter(strings.NewReader(input), &bytes.Buffer{}, false, FilterOptions{
		MinLen: 20000, // drop
		OutBad: &errWriter{},
	})
	if err == nil {
		t.Errorf("expected error when bad-output writer Write fails")
	}
	// Small record exercises the Flush path on badWriter.
	err = Filter(strings.NewReader(">drop\nA\n"), &bytes.Buffer{}, false, FilterOptions{
		MinLen: 5,
		OutBad: &errWriter{},
	})
	if err == nil {
		t.Errorf("expected error when bad-output Flush fails")
	}
}

func TestFilterFastqOutBadWriterError(t *testing.T) {
	seq := strings.Repeat("A", 16384)
	qual := strings.Repeat("I", len(seq))
	input := "@drop\n" + seq + "\n+\n" + qual + "\n"
	err := Filter(strings.NewReader(input), &bytes.Buffer{}, true, FilterOptions{
		MinLen: 20000,
		OutBad: &errWriter{},
	})
	if err == nil {
		t.Errorf("expected error when fastq bad-output writer Write fails")
	}
}

func TestFilterFastaDerepBadWriterError(t *testing.T) {
	big := strings.Repeat("ACGT", 4096)
	input := ">a\n" + big + "\n>b\n" + big + "\n"
	err := Filter(strings.NewReader(input), &bytes.Buffer{}, false, FilterOptions{
		Derep:    1,
		DerepMin: 2,
		OutBad:   &errWriter{},
	})
	if err == nil {
		t.Errorf("expected error when derep bad-output writer Write fails (fasta)")
	}
}

func TestFilterFastqDerepBadWriterError(t *testing.T) {
	seq := strings.Repeat("ACGT", 4096)
	qual := strings.Repeat("I", len(seq))
	input := "@a\n" + seq + "\n+\n" + qual + "\n@b\n" + seq + "\n+\n" + qual + "\n"
	err := Filter(strings.NewReader(input), &bytes.Buffer{}, true, FilterOptions{
		Derep:    1,
		DerepMin: 2,
		OutBad:   &errWriter{},
	})
	if err == nil {
		t.Errorf("expected error when derep bad-output writer Write fails (fastq)")
	}
}

func TestFilterPairedFastaWriterError(t *testing.T) {
	in1 := hugeFastaRecord()
	in2 := hugeFastaRecord()
	// writer1 errors -> Write to writer1 returns error
	if err := FilterPaired(strings.NewReader(in1), strings.NewReader(in2),
		&errWriter{}, &bytes.Buffer{}, false, FilterOptions{}); err == nil {
		t.Errorf("expected error from paired fasta writer1 Write")
	}

	// writer2 errors -> Write to writer2 returns error (after writer1 succeeds)
	if err := FilterPaired(strings.NewReader(in1), strings.NewReader(in2),
		&bytes.Buffer{}, &errWriter{}, false, FilterOptions{}); err == nil {
		t.Errorf("expected error from paired fasta writer2 Write")
	}

	// Small records exercise the Flush paths.
	if err := FilterPaired(strings.NewReader(">a\nACGT\n"), strings.NewReader(">a\nACGT\n"),
		&errWriter{}, &bytes.Buffer{}, false, FilterOptions{}); err == nil {
		t.Errorf("expected error from paired fasta writer1 Flush")
	}
}

func TestFilterPairedFastqWriterError(t *testing.T) {
	in1 := hugeFastqRecord()
	in2 := hugeFastqRecord()
	if err := FilterPaired(strings.NewReader(in1), strings.NewReader(in2),
		&errWriter{}, &bytes.Buffer{}, true, FilterOptions{}); err == nil {
		t.Errorf("expected error from paired fastq writer1 Write")
	}

	if err := FilterPaired(strings.NewReader(in1), strings.NewReader(in2),
		&bytes.Buffer{}, &errWriter{}, true, FilterOptions{}); err == nil {
		t.Errorf("expected error from paired fastq writer2 Write")
	}

	// Small records buffer internally inside the fastq writer; the per-record
	// Write returns nil and the Flush() error path is taken instead.
	if err := FilterPaired(strings.NewReader("@a\nACGT\n+\nIIII\n"),
		strings.NewReader("@a\nACGT\n+\nIIII\n"),
		&errWriter{}, &bytes.Buffer{}, true, FilterOptions{}); err == nil {
		t.Errorf("expected error from paired fastq writer1 Flush")
	}
}

// --- prinseq.go: filterPaired{Fastq,Fasta} reader error paths ------------

func TestFilterPairedFastaReaderError(t *testing.T) {
	// reader1 errors immediately -> err1 returned.
	err := FilterPaired(errReader{}, strings.NewReader(">a\nACGT\n"),
		&bytes.Buffer{}, &bytes.Buffer{}, false, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from paired fasta reader1")
	}
	// reader2 errors immediately -> err2 returned.
	err = FilterPaired(strings.NewReader(">a\nACGT\n"), errReader{},
		&bytes.Buffer{}, &bytes.Buffer{}, false, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from paired fasta reader2")
	}
}

func TestFilterPairedFastqReaderError(t *testing.T) {
	err := FilterPaired(errReader{}, strings.NewReader("@a\nACGT\n+\nIIII\n"),
		&bytes.Buffer{}, &bytes.Buffer{}, true, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from paired fastq reader1")
	}
	err = FilterPaired(strings.NewReader("@a\nACGT\n+\nIIII\n"), errReader{},
		&bytes.Buffer{}, &bytes.Buffer{}, true, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from paired fastq reader2")
	}
}

func TestFilterPairedFastaShouldFilter(t *testing.T) {
	// One pair fails MinLen so the `continue` branch in filterPairedFasta is taken.
	in1 := ">a\nACGTACGTACGT\n>b\nAC\n"
	in2 := ">a\nTTTTTTTTTTTT\n>b\nGG\n"
	var o1, o2 bytes.Buffer
	err := FilterPaired(strings.NewReader(in1), strings.NewReader(in2),
		&o1, &o2, false, FilterOptions{MinLen: 5})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if !strings.Contains(o1.String(), ">a") || strings.Contains(o1.String(), ">b") {
		t.Errorf("paired fasta out1 wrong: %q", o1.String())
	}
	if !strings.Contains(o2.String(), ">a") || strings.Contains(o2.String(), ">b") {
		t.Errorf("paired fasta out2 wrong: %q", o2.String())
	}
}

func TestFilterPairedFastqMismatchedCounts(t *testing.T) {
	in1 := "@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n"
	in2 := "@a\nACGT\n+\nIIII\n"
	err := FilterPaired(strings.NewReader(in1), strings.NewReader(in2),
		&bytes.Buffer{}, &bytes.Buffer{}, true, FilterOptions{})
	if err == nil {
		t.Errorf("expected error for mismatched fastq pair counts")
	}
}

// --- prinseq.go: filter*  reader error paths -----------------------------

func TestFilterFastaReaderError(t *testing.T) {
	err := Filter(errReader{}, &bytes.Buffer{}, false, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from filter when reader errors")
	}
}

func TestFilterFastqReaderError(t *testing.T) {
	err := Filter(errReader{}, &bytes.Buffer{}, true, FilterOptions{})
	if err == nil {
		t.Errorf("expected error from fastq filter when reader errors")
	}
}

// --- api.go: error paths -------------------------------------------------

func TestAPIHandleStatsBodyReadError(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost, "/api/stats", errReader{})
	rr := httptest.NewRecorder()
	s.handleStats(rr, req)
	if rr.Code != http.StatusBadRequest {
		t.Errorf("status = %d, want 400 for body read error", rr.Code)
	}
}

func TestAPIHandleFilterBodyReadError(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost, "/api/filter?format=fastq", errReader{})
	rr := httptest.NewRecorder()
	s.handleFilter(rr, req)
	if rr.Code != http.StatusBadRequest {
		t.Errorf("status = %d, want 400 for filter body read error", rr.Code)
	}
}

func TestAPIHandleFilterBadInput(t *testing.T) {
	s := NewAPIServer(":0")
	// Invalid FASTQ -> Filter returns an error -> 500.
	req := httptest.NewRequest(http.MethodPost,
		"/api/filter?format=fastq", strings.NewReader("not fastq at all"))
	rr := httptest.NewRecorder()
	s.handleFilter(rr, req)
	if rr.Code != http.StatusInternalServerError {
		t.Errorf("status = %d, want 500 for bad fastq filter", rr.Code)
	}
}

func TestAPIHandleBenchmarkBodyReadError(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost, "/api/benchmark", errReader{})
	rr := httptest.NewRecorder()
	s.handleBenchmark(rr, req)
	if rr.Code != http.StatusBadRequest {
		t.Errorf("status = %d, want 400 for benchmark body read error", rr.Code)
	}
}

func TestAPIHandleBenchmarkBadInput(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost,
		"/api/benchmark?format=fastq", strings.NewReader("not fastq"))
	rr := httptest.NewRecorder()
	s.handleBenchmark(rr, req)
	if rr.Code != http.StatusInternalServerError {
		t.Errorf("status = %d, want 500 for bad fastq benchmark", rr.Code)
	}
}

func TestAPIHandleReportBodyReadError(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost, "/api/report", errReader{})
	rr := httptest.NewRecorder()
	s.handleReport(rr, req)
	if rr.Code != http.StatusBadRequest {
		t.Errorf("status = %d, want 400 for report body read error", rr.Code)
	}
}

func TestAPIHandleReportBadStats(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost,
		"/api/report?format=fastq", strings.NewReader("not fastq"))
	rr := httptest.NewRecorder()
	s.handleReport(rr, req)
	if rr.Code != http.StatusInternalServerError {
		t.Errorf("status = %d, want 500 for bad fastq report", rr.Code)
	}
}

func TestAPIHandleReportGenerationError(t *testing.T) {
	// An empty-sequence FASTQ produces NaN GCContent which causes
	// json.MarshalIndent inside GenerateHTMLReport to error. Exercises
	// the "Error generating report" branch in handleReport.
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost,
		"/api/report?format=fastq", strings.NewReader("@a\n\n+\n\n"))
	rr := httptest.NewRecorder()
	s.handleReport(rr, req)
	if rr.Code != http.StatusInternalServerError {
		t.Errorf("status = %d, want 500 for report-generation failure", rr.Code)
	}
	if !strings.Contains(rr.Body.String(), "Error generating report") {
		t.Errorf("expected 'Error generating report' in body, got %q", rr.Body.String())
	}
}

func TestAPIHandleGraphBodyReadError(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost, "/api/graph", errReader{})
	rr := httptest.NewRecorder()
	s.handleGraph(rr, req)
	if rr.Code != http.StatusBadRequest {
		t.Errorf("status = %d, want 400 for graph body read error", rr.Code)
	}
}

func TestAPIHandleGraphBadStats(t *testing.T) {
	s := NewAPIServer(":0")
	req := httptest.NewRequest(http.MethodPost,
		"/api/graph?format=fastq", strings.NewReader("not fastq"))
	rr := httptest.NewRecorder()
	s.handleGraph(rr, req)
	if rr.Code != http.StatusInternalServerError {
		t.Errorf("status = %d, want 500 for bad fastq graph stats", rr.Code)
	}
}

func TestAPIHandleGraphGenerationError(t *testing.T) {
	s := NewAPIServer(":0")
	// FASTA input + type=quality -> no QualityDistribution -> GenerateGraph errors.
	req := httptest.NewRequest(http.MethodPost,
		"/api/graph?format=fasta&type=quality", strings.NewReader(">a\nACGT\n"))
	rr := httptest.NewRecorder()
	s.handleGraph(rr, req)
	if rr.Code != http.StatusInternalServerError {
		t.Errorf("status = %d, want 500 for graph generation failure", rr.Code)
	}
}

// --- api.go: Start ----- ------------------------------------------------

func TestAPIStartBindError(t *testing.T) {
	// Asking Start to bind to an invalid address forces ListenAndServe
	// to return an error synchronously. This covers the Start function
	// without leaving a goroutine running.
	s := NewAPIServer("invalid-address-without-port")
	if err := s.Start(); err == nil {
		t.Errorf("expected Start to fail on invalid address")
	}
}

// --- benchmark.go: read-error paths --------------------------------------

func TestBenchmarkReadErrors(t *testing.T) {
	if _, _, err := BenchmarkStats(errReader{}, false); err == nil {
		t.Errorf("expected BenchmarkStats to surface read error")
	}
	if _, err := BenchmarkFilter(errReader{}, false, FilterOptions{}); err == nil {
		t.Errorf("expected BenchmarkFilter to surface read error")
	}
	if _, err := RunBenchmarkSuite(errReader{}, false); err == nil {
		t.Errorf("expected RunBenchmarkSuite to surface read error")
	}
}

// --- benchmark.go: enhanced-stats and filter inner errors ---------------

func TestBenchmarkStatsInvalidFastq(t *testing.T) {
	// Provide invalid FASTQ so the underlying CalculateEnhancedStats errors.
	_, _, err := BenchmarkStats(strings.NewReader("not fastq"), true)
	if err == nil {
		t.Errorf("expected BenchmarkStats to return error for invalid fastq")
	}
}

func TestBenchmarkFilterInvalidFastq(t *testing.T) {
	_, err := BenchmarkFilter(strings.NewReader("not fastq"), true, FilterOptions{})
	if err == nil {
		t.Errorf("expected BenchmarkFilter to return error for invalid fastq")
	}
}

// --- graph.go: writer errors at each SVG sub-section ---------------------

func TestGenerateSVGWriterErrors(t *testing.T) {
	stats, err := CalculateEnhancedStats(strings.NewReader(
		"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n!!!!\n"), true)
	if err != nil {
		t.Fatalf("stats: %v", err)
	}

	// Fail on the very first Fprintf (the SVG header).
	if err := GenerateSVG(stats, &errWriter{failAfter: 0}); err == nil {
		t.Errorf("expected error when SVG header write fails")
	}
	// Empirical write boundaries (cumulative bytes written when each Write
	// returns to the underlying writer):
	//   369 (xml+svg header), 469 (length title), 621 (length axes),
	//   695 (length bar), 702 (length </g>),
	//   810 (quality title), 962 (quality axes), 1036/1112 (quality bars),
	//   1119 (quality </g>), 1226 (positional title), 1379 (positional axes),
	//   1388 (positional M), 1400/1412/1424 (positional L), 1441 (class=line),
	//   1448 (positional </g>), 1455 (svg close).
	// To make the Fprintf landing at offset O fail, choose a budget that
	// allows the previous Write to succeed but is smaller than O.
	for _, budget := range []int{
		// header, length-svg sub-writes
		100, 400, 500, 650, 700,
		// quality-svg sub-writes
		800, 900, 1000, 1100, 1115,
		// positional-svg sub-writes (including 1425 to fail the
		// `class="line"` Fprintf and 1442 to fail the closing </g>)
		1200, 1300, 1383, 1395, 1407, 1419, 1425, 1442,
	} {
		if err := GenerateSVG(stats, &errWriter{failAfter: budget}); err == nil {
			t.Errorf("expected error at SVG write budget %d", budget)
		}
	}
}

// --- graph.go: positional-quality bar drawn (line 240 loop body) ---------

func TestPositionalQualityNonzeroBars(t *testing.T) {
	// Need positional quality values with a real spread so the loop
	// `for j := 0; j < barWidth; j++` body executes (barWidth > 0).
	// Read positions cycle 'I' (Q40) and '!' (Q0) so positions 1,3,5,...
	// average higher than 2,4,...
	input := "@a\nACGTACGTAC\n+\nI!I!I!I!I!\n" +
		"@b\nACGTACGTAC\n+\nI!I!I!I!I!\n"
	stats, err := CalculateEnhancedStats(strings.NewReader(input), true)
	if err != nil {
		t.Fatalf("stats: %v", err)
	}
	var buf bytes.Buffer
	if err := GenerateGraph(stats, GraphTypePositional, &buf); err != nil {
		t.Fatalf("graph err: %v", err)
	}
	// The non-zero bar character is the UTF-8 full block (U+2588).
	if !strings.Contains(buf.String(), "█") {
		t.Errorf("expected at least one positional-quality bar character, got %q",
			buf.String())
	}
}

// --- report.go: error paths ----------------------------------------------

func TestGenerateHTMLReportJSONMarshalError(t *testing.T) {
	// json.MarshalIndent fails on NaN/Inf floats. Inject one through
	// PositionalQuality to drive the "error marshaling stats" return.
	stats := &Stats{
		MinLength:           1,
		LengthDistribution:  map[int]int{},
		QualityDistribution: map[int]int{},
		BaseComposition:     map[string]int{},
		Dinucleotides:       map[string]int{},
		PositionalQuality:   []float64{nanFloat()},
	}
	err := GenerateHTMLReport(stats, &bytes.Buffer{})
	if err == nil {
		t.Errorf("expected json marshaling error for NaN PositionalQuality")
	}
	if !strings.Contains(err.Error(), "marshaling stats") {
		t.Errorf("expected 'marshaling stats' wrapping, got %v", err)
	}
}

// nanFloat returns a NaN as float64 without pulling in math.NaN at module
// init time (so the value is unambiguous in the test source).
func nanFloat() float64 {
	zero := 0.0
	return zero / zero
}

func TestGenerateHTMLReportSVGError(t *testing.T) {
	stats, err := CalculateEnhancedStats(strings.NewReader(">a\nACGT\n"), false)
	if err != nil {
		t.Fatalf("stats: %v", err)
	}
	// Writer fails immediately, so GenerateSVG errors and the report wrapper
	// returns the wrapped "error generating SVG" message.
	err = GenerateHTMLReport(stats, &errWriter{failAfter: 0})
	if err == nil {
		t.Errorf("expected error from GenerateHTMLReport when SVG write fails")
	}
}

func TestGenerateHTMLReportBodyWriteError(t *testing.T) {
	// Stats with AvgQuality > 0 forces the optional middle Fprintf section
	// (so we get to cover both error returns there). Use enough capacity in
	// failAfter to get past the SVG section, then fail in the main body.
	stats, err := CalculateEnhancedStats(strings.NewReader(
		"@a\nACGT\n+\nIIII\n"), true)
	if err != nil {
		t.Fatalf("stats: %v", err)
	}
	if stats.AvgQuality <= 0 {
		t.Fatalf("setup: AvgQuality should be > 0")
	}
	// Fail after enough bytes that the SVG completes but the main HTML
	// body Fprintf fails.
	err = GenerateHTMLReport(stats, &errWriter{failAfter: 4000})
	if err == nil {
		t.Errorf("expected error from GenerateHTMLReport when HTML write fails")
	}

	// Likewise, fail during the optional "Average Quality Score" Fprintf.
	// Try a slightly larger budget so we get past the first body Fprintf
	// but stop inside the second.
	for _, budget := range []int{4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5200, 5400, 5800, 6200} {
		if err := GenerateHTMLReport(stats, &errWriter{failAfter: budget}); err == nil {
			t.Errorf("expected error from GenerateHTMLReport at budget %d", budget)
		}
	}
}

// --- parallel.go: error paths --------------------------------------------

func TestFilterFilesParallelOutputDirError(t *testing.T) {
	// Create a regular file then point outputDir at a path *under* it so
	// MkdirAll cannot succeed (since the parent isn't a directory).
	dir := t.TempDir()
	regular := filepath.Join(dir, "notadir")
	if err := os.WriteFile(regular, []byte("x"), 0644); err != nil {
		t.Fatal(err)
	}
	badOut := filepath.Join(regular, "subdir")
	err := FilterFilesParallel([]string{filepath.Join(dir, "any.fasta")}, badOut, false, FilterOptions{}, 1)
	if err == nil {
		t.Errorf("expected error creating output dir under a regular file")
	}
}

func TestFilterSingleFileOutputCreateError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fasta")
	if err := os.WriteFile(in, []byte(">a\nACGT\n"), 0644); err != nil {
		t.Fatal(err)
	}
	// Put a *directory* where the output file should go so os.Create fails.
	outDir := filepath.Join(dir, "outroot")
	if err := os.MkdirAll(outDir, 0755); err != nil {
		t.Fatal(err)
	}
	conflict := filepath.Join(outDir, filepath.Base(in))
	if err := os.MkdirAll(conflict, 0755); err != nil {
		t.Fatal(err)
	}
	// FilterFilesParallel will join outDir + base(in) = conflict (a dir) -> Create fails.
	if err := FilterFilesParallel([]string{in}, outDir, false, FilterOptions{}, 1); err == nil {
		t.Errorf("expected error when output path is a directory")
	}
}

func TestBatchProcessOutputDirError(t *testing.T) {
	dir := t.TempDir()
	regular := filepath.Join(dir, "notadir")
	if err := os.WriteFile(regular, []byte("x"), 0644); err != nil {
		t.Fatal(err)
	}
	badOut := filepath.Join(regular, "subdir")
	_, err := BatchProcess(BatchProcessConfig{
		InputFiles: []string{filepath.Join(dir, "any.fasta")},
		OutputDir:  badOut,
	})
	if err == nil {
		t.Errorf("expected error creating BatchProcess output dir under a file")
	}
}

func TestBatchProcessStatsError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	// Invalid FASTQ -> CalculateEnhancedStats errors.
	if err := os.WriteFile(in, []byte("not fastq"), 0644); err != nil {
		t.Fatal(err)
	}
	results, err := BatchProcess(BatchProcessConfig{
		InputFiles: []string{in},
		OutputDir:  filepath.Join(dir, "out"),
		IsFastq:    true,
	})
	if err != nil {
		t.Fatalf("BatchProcess: %v", err)
	}
	if len(results) != 1 || results[0].Error == nil {
		t.Errorf("expected stats error result")
	}
	if results[0].Error != nil && !strings.Contains(results[0].Error.Error(), "stats") {
		t.Errorf("expected stats-error wrapping, got %v", results[0].Error)
	}
}

func TestBatchProcessReportCreateError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "input.fasta")
	if err := os.WriteFile(in, []byte(">a\nACGT\n"), 0644); err != nil {
		t.Fatal(err)
	}
	outDir := filepath.Join(dir, "out")
	if err := os.MkdirAll(outDir, 0755); err != nil {
		t.Fatal(err)
	}
	// Force os.Create for the .html path to fail by pre-creating a directory
	// at exactly the expected report path.
	reportPath := filepath.Join(outDir, filepath.Base(in)+".html")
	if err := os.MkdirAll(reportPath, 0755); err != nil {
		t.Fatal(err)
	}
	results, err := BatchProcess(BatchProcessConfig{
		InputFiles:     []string{in},
		OutputDir:      outDir,
		GenerateReport: true,
	})
	if err != nil {
		t.Fatalf("BatchProcess: %v", err)
	}
	if len(results) != 1 || results[0].Error == nil {
		t.Fatalf("expected report-create error, got %+v", results)
	}
	if !strings.Contains(results[0].Error.Error(), "creating report") {
		t.Errorf("expected 'creating report' wrapping, got %v", results[0].Error)
	}
}

func TestBatchProcessFilterCreateError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "input.fasta")
	if err := os.WriteFile(in, []byte(">a\nACGT\n"), 0644); err != nil {
		t.Fatal(err)
	}
	outDir := filepath.Join(dir, "out")
	if err := os.MkdirAll(outDir, 0755); err != nil {
		t.Fatal(err)
	}
	// Force the filtered output Create to fail.
	filteredPath := filepath.Join(outDir, "filtered_"+filepath.Base(in))
	if err := os.MkdirAll(filteredPath, 0755); err != nil {
		t.Fatal(err)
	}
	results, err := BatchProcess(BatchProcessConfig{
		InputFiles: []string{in},
		OutputDir:  outDir,
	})
	if err != nil {
		t.Fatalf("BatchProcess: %v", err)
	}
	if len(results) != 1 || results[0].Error == nil {
		t.Fatalf("expected filter-create error, got %+v", results)
	}
	if !strings.Contains(results[0].Error.Error(), "creating output") {
		t.Errorf("expected 'creating output' wrapping, got %v", results[0].Error)
	}
}

// TestBatchProcessFilterRunError makes the inline Filter call fail by
// supplying a FilterOptions.OutBad that always errors and an input that
// will be rejected (sent to OutBad). Covers the "error filtering" branch
// in processBatchFile.
func TestBatchProcessFilterRunError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fasta")
	// One short sequence -> rejected by MinLen -> sent to OutBad which errors.
	if err := os.WriteFile(in, []byte(">drop\nAC\n"), 0644); err != nil {
		t.Fatal(err)
	}
	outDir := filepath.Join(dir, "out")
	results, err := BatchProcess(BatchProcessConfig{
		InputFiles: []string{in},
		OutputDir:  outDir,
		FilterOpts: FilterOptions{
			MinLen: 5,
			OutBad: &errWriter{}, // explodes when Filter routes a record here
		},
	})
	if err != nil {
		t.Fatalf("BatchProcess: %v", err)
	}
	if len(results) != 1 || results[0].Error == nil {
		t.Fatalf("expected filter error result, got %+v", results)
	}
	if !strings.Contains(results[0].Error.Error(), "filtering") {
		t.Errorf("expected 'filtering' wrapping, got %v", results[0].Error)
	}
}

// TestBatchProcessReportGenerateError drives processBatchFile through the
// "error generating report" branch by giving it an empty-sequence FASTQ.
// CalculateEnhancedStats then yields stats with GCContent == NaN, and
// json.MarshalIndent inside GenerateHTMLReport returns an error.
func TestBatchProcessReportGenerateError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "empty.fastq")
	// One FASTQ record with an empty sequence -> NaN GCContent in stats.
	if err := os.WriteFile(in, []byte("@a\n\n+\n\n"), 0644); err != nil {
		t.Fatal(err)
	}
	outDir := filepath.Join(dir, "out")
	results, err := BatchProcess(BatchProcessConfig{
		InputFiles:     []string{in},
		OutputDir:      outDir,
		IsFastq:        true,
		GenerateReport: true,
	})
	if err != nil {
		t.Fatalf("BatchProcess: %v", err)
	}
	if len(results) != 1 || results[0].Error == nil {
		t.Fatalf("expected an error result, got %+v", results)
	}
	if !strings.Contains(results[0].Error.Error(), "generating report") {
		t.Errorf("expected 'generating report' wrapping, got %v",
			results[0].Error)
	}
}

func TestBatchProcessFilterError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "input.fastq")
	// Valid FASTQ for stats, but we'll set IsFastq=false so Filter runs as
	// FASTA on a FASTQ source. That's lenient; we instead force Filter
	// to fail with a malformed input that nonetheless calculates stats.
	// Strategy: provide a valid FASTQ for stats, then truncate during filter
	// by switching IsFastq -- but Filter is invoked with config.IsFastq.
	//
	// Simpler: use a FASTA input but make the filtered_output path a
	// non-writable file via permissions. We already cover Create error
	// above; instead cover "error filtering" via an unreadable-mid-stream
	// trick: write a header but no body so the fastq reader errors.
	if err := os.WriteFile(in, []byte("@a\nACGT\n+\nII\n"), 0644); err != nil {
		// quality length (2) != sequence length (4) -> filterFastq errors mid-stream
		t.Fatal(err)
	}
	outDir := filepath.Join(dir, "out")
	results, err := BatchProcess(BatchProcessConfig{
		InputFiles: []string{in},
		OutputDir:  outDir,
		IsFastq:    true,
	})
	if err != nil {
		t.Fatalf("BatchProcess: %v", err)
	}
	if len(results) != 1 {
		t.Fatalf("got %d results", len(results))
	}
	// Mismatched FASTQ should fail somewhere in the pipeline. Stats might
	// catch it first; either way an error should be reported. The "error
	// filtering" path is covered when stats happens to succeed first, so
	// we just assert that *some* error is produced.
	if results[0].Error == nil {
		t.Errorf("expected an error result from malformed FASTQ")
	}
}

// --- prinseq.go: incomplete enhanced fastq record -------------------------

func TestCalculateEnhancedFastqIncompleteRecord(t *testing.T) {
	// Only 3 of 4 lines provided -> "incomplete FASTQ record".
	input := "@a\nACGT\n+\n"
	_, err := CalculateEnhancedStats(strings.NewReader(input), true)
	if err == nil {
		t.Errorf("expected error for incomplete enhanced FASTQ record")
	}
}

func TestCalculateEnhancedFastqHeaderInvalid(t *testing.T) {
	// Header line doesn't start with '@'
	input := "X\nACGT\n+\nIIII\n"
	if _, err := CalculateEnhancedStats(strings.NewReader(input), true); err == nil {
		t.Errorf("expected error for bad enhanced FASTQ header")
	}
}

func TestCalculateEnhancedFastqPlusInvalid(t *testing.T) {
	input := "@a\nACGT\nX\nIIII\n"
	if _, err := CalculateEnhancedStats(strings.NewReader(input), true); err == nil {
		t.Errorf("expected error for bad enhanced FASTQ plus line")
	}
}

func TestCalculateEnhancedFastqQualLenMismatch(t *testing.T) {
	input := "@a\nACGTACGT\n+\nIII\n"
	if _, err := CalculateEnhancedStats(strings.NewReader(input), true); err == nil {
		t.Errorf("expected error for enhanced FASTQ quality length mismatch")
	}
}
