package skewer

import (
	"bytes"
	"compress/gzip"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// errWriter is an io.Writer that returns errAfter once nWritten >= afterBytes.
// Used to exercise error paths in fastq.Writer.Write / fastq.Writer.Flush.
type errWriter struct {
	afterBytes int
	written    int
	errAfter   error
}

func (e *errWriter) Write(p []byte) (int, error) {
	if e.written >= e.afterBytes {
		return 0, e.errAfter
	}
	canWrite := e.afterBytes - e.written
	if canWrite > len(p) {
		canWrite = len(p)
	}
	e.written += canWrite
	if canWrite < len(p) {
		return canWrite, e.errAfter
	}
	return canWrite, nil
}

// bigSeq returns a sequence of length n made of repeated "ACGT".
func bigSeq(n int) string {
	pat := "ACGT"
	var b strings.Builder
	for b.Len() < n {
		b.WriteString(pat)
	}
	return b.String()[:n]
}

// makeFASTQ builds a FASTQ string from records of (id, seq, qual).
func makeFASTQ(records ...[3]string) string {
	var b strings.Builder
	for _, r := range records {
		b.WriteString("@" + r[0] + "\n")
		b.WriteString(r[1] + "\n")
		b.WriteString("+\n")
		b.WriteString(r[2] + "\n")
	}
	return b.String()
}

func qstring(n int, c byte) string {
	return strings.Repeat(string(c), n)
}

func TestTrimSingleEndEndToEnd(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	core := "ACGTACGTACGTACGTACGT" // 20 bp
	in := makeFASTQ(
		[3]string{"keep", core + adapter, qstring(len(core)+len(adapter), 'I')},
		[3]string{"noadapter", core, qstring(len(core), 'I')},
		[3]string{"drop", "ACGT" + adapter, qstring(4+len(adapter), 'I')}, // 4 bp after trim, below min length
	)

	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.TotalReads != 3 {
		t.Errorf("TotalReads = %d, want 3", stats.TotalReads)
	}
	if stats.AdapterFound3 != 2 {
		t.Errorf("AdapterFound3 = %d, want 2", stats.AdapterFound3)
	}
	if stats.DiscardedReads != 1 {
		t.Errorf("DiscardedReads = %d, want 1", stats.DiscardedReads)
	}
	if stats.TrimmedReads != 1 {
		t.Errorf("TrimmedReads = %d, want 1 (only 'keep' written and shorter)", stats.TrimmedReads)
	}
	res := out.String()
	if strings.Contains(res, adapter) {
		t.Errorf("output still contains adapter: %q", res)
	}
	// keep read sequence should be exactly core
	if !strings.Contains(res, "@keep\n"+core+"\n") {
		t.Errorf("expected trimmed 'keep' read with sequence %q, got %q", core, res)
	}
	// noadapter read should be passed through unchanged
	if !strings.Contains(res, "@noadapter\n"+core+"\n") {
		t.Errorf("expected untrimmed 'noadapter' read, got %q", res)
	}
}

func TestTrimSingleEndEmptyInput(t *testing.T) {
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	stats, err := TrimSingleEnd(strings.NewReader(""), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd empty: %v", err)
	}
	if stats.TotalReads != 0 {
		t.Errorf("TotalReads = %d, want 0", stats.TotalReads)
	}
	if out.Len() != 0 {
		t.Errorf("expected empty output, got %q", out.String())
	}
}

func TestTrimSingleEndReadEntirelyAdapter(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	in := makeFASTQ([3]string{"alladapter", adapter, qstring(len(adapter), 'I')})
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 1

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.DiscardedReads != 1 {
		t.Errorf("DiscardedReads = %d, want 1", stats.DiscardedReads)
	}
	if out.Len() != 0 {
		t.Errorf("expected empty output, got %q", out.String())
	}
}

func TestTrimSingleEndAdapterLongerThanRead(t *testing.T) {
	// adapter longer than the read; only a prefix of the adapter could match.
	adapter := "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	read := "ACGTACGTACGT" + "AGATCGGAAGAGC" // 12 + 13 = 25 bp, adapter prefix matches at pos 12
	in := makeFASTQ([3]string{"r1", read, qstring(len(read), 'I')})
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 5

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.AdapterFound3 != 1 {
		t.Errorf("AdapterFound3 = %d, want 1", stats.AdapterFound3)
	}
	if !strings.Contains(out.String(), "@r1\nACGTACGTACGT\n") {
		t.Errorf("expected read trimmed to %q, got %q", "ACGTACGTACGT", out.String())
	}
}

func TestTrimSingleEndQualityTrimming(t *testing.T) {
	// 20 bp read, last 8 bases low quality ('#' = Q2 in Phred33), threshold 20.
	seq := "ACGTACGTACGTACGTACGT"
	qual := qstring(12, 'I') + qstring(8, '#')
	in := makeFASTQ([3]string{"q1", seq, qual})
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 20
	opts.MinLength = 5

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.TrimmedReads != 1 {
		t.Errorf("TrimmedReads = %d, want 1", stats.TrimmedReads)
	}
	if !strings.Contains(out.String(), "@q1\nACGTACGTACGT\n") {
		t.Errorf("expected quality-trimmed read %q, got %q", "ACGTACGTACGT", out.String())
	}
}

func TestTrimSingleEndQualityTrimmingBothEnds(t *testing.T) {
	// low quality at both ends.
	seq := "ACGTACGTACGTACGTACGT"
	qual := qstring(4, '#') + qstring(10, 'I') + qstring(6, '#')
	in := makeFASTQ([3]string{"q1", seq, qual})
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 20
	opts.MinLength = 5

	if _, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts); err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	want := seq[4:14]
	if !strings.Contains(out.String(), "@q1\n"+want+"\n") {
		t.Errorf("expected %q, got %q", want, out.String())
	}
}

func TestTrimSingleEndQualityTrimsBelowMinLength(t *testing.T) {
	// Entire read is low quality -> trimmed to nothing -> discarded.
	seq := "ACGTACGTACGTACGTACGT"
	qual := qstring(20, '#')
	in := makeFASTQ([3]string{"q1", seq, qual})
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 20
	opts.MinLength = 5

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.DiscardedReads != 1 {
		t.Errorf("DiscardedReads = %d, want 1", stats.DiscardedReads)
	}
	if out.Len() != 0 {
		t.Errorf("expected empty output, got %q", out.String())
	}
}

func TestTrimByQuality(t *testing.T) {
	tests := []struct {
		name               string
		qual               string
		threshold          int
		start, end         int
		wantStart, wantEnd int
	}{
		{"no trim", "IIIIII", 20, 0, 6, 0, 6},
		{"trim 3prime", "III###", 20, 0, 6, 0, 3},
		{"trim 5prime", "###III", 20, 0, 6, 3, 6},
		{"trim both", "##II##", 20, 0, 6, 2, 4},
		// When every base is below threshold, 3' trimming collapses end down
		// to start; the result is an empty range (start == end), reported as
		// the original start value.
		{"all low", "######", 20, 0, 6, 0, 0},
		{"sub-slice offset", "II###", 20, 2, 7, 2, 4},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gs, ge := trimByQuality([]byte(tt.qual), tt.threshold, tt.start, tt.end)
			if gs != tt.wantStart || ge != tt.wantEnd {
				t.Errorf("trimByQuality() = (%d,%d), want (%d,%d)", gs, ge, tt.wantStart, tt.wantEnd)
			}
		})
	}
}

func TestImprovedFindAdapterMismatchRates(t *testing.T) {
	adapter := "AGATCGGAAGAGC" // 13 bp
	// one mismatch in the adapter region at position 12+5
	seqOneErr := "ACGTACGTACGT" + "AGATCAGAAGAGC"
	tests := []struct {
		name       string
		seq        string
		minOverlap int
		errorRate  float64
		want       int
	}{
		{"exact", "ACGTACGTACGT" + adapter, 3, 0.1, 12},
		{"one error within rate", seqOneErr, 3, 0.1, 12}, // maxErrors = floor(13*0.1)=1
		{"one error rejected at zero rate", seqOneErr, 3, 0.0, -1},
		{"min overlap too long", "AGA", 5, 0.1, -1},
		{"empty adapter", "ACGT", 3, 0.1, -1},
		{"partial overlap at very end", "ACGTACGTACGTAGAT", 4, 0.1, 12},
		{"partial overlap below overlap", "ACGTACGTACGTAGA", 4, 0.1, -1},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := adapter
			if tt.name == "empty adapter" {
				a = ""
			}
			got := improvedFindAdapter(tt.seq, a, tt.minOverlap, tt.errorRate)
			if got != tt.want {
				t.Errorf("improvedFindAdapter(%q, %q, %d, %g) = %d, want %d", tt.seq, a, tt.minOverlap, tt.errorRate, got, tt.want)
			}
		})
	}
}

func TestDetectAdapters(t *testing.T) {
	// "AGATCGGAAGAGC" is first in commonAdapters, so detection deterministically
	// picks it when reads carry the TruSeq adapter at the 3' end.
	adapter := "AGATCGGAAGAGC"
	// 10 reads, all carrying the adapter at the 3' end.
	var recs [][3]string
	for i := 0; i < 10; i++ {
		recs = append(recs, [3]string{"r", "ACGTACGTACGTACGT" + adapter, qstring(16+len(adapter), 'I')})
	}
	in := makeFASTQ(recs...)
	reader := fastq.NewReader(strings.NewReader(in), fastq.Phred33)
	opts, err := detectAdapters(reader, 100)
	if err != nil {
		t.Fatalf("detectAdapters: %v", err)
	}
	if opts.Adapter3 != adapter {
		t.Errorf("detected Adapter3 = %q, want %q", opts.Adapter3, adapter)
	}
}

func TestDetectAdaptersEmpty(t *testing.T) {
	reader := fastq.NewReader(strings.NewReader(""), fastq.Phred33)
	opts, err := detectAdapters(reader, 100)
	if err != nil {
		t.Fatalf("detectAdapters: %v", err)
	}
	if opts.Adapter3 != "" || opts.Adapter5 != "" {
		t.Errorf("expected no adapters detected, got %q / %q", opts.Adapter3, opts.Adapter5)
	}
}

func TestTrimSingleEndAutoDetect(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	var recs [][3]string
	for i := 0; i < 20; i++ {
		recs = append(recs, [3]string{"r", "ACGTACGTACGTACGTACGT" + adapter, qstring(20+len(adapter), 'I')})
	}
	in := makeFASTQ(recs...)
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.AutoDetect = true
	opts.AutoDetectReads = 50
	opts.MinLength = 10

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd auto-detect: %v", err)
	}
	if stats.DetectedAdapter3 != adapter {
		t.Errorf("DetectedAdapter3 = %q, want %q", stats.DetectedAdapter3, adapter)
	}
	if stats.AdapterFound3 != 20 {
		t.Errorf("AdapterFound3 = %d, want 20", stats.AdapterFound3)
	}
	if strings.Contains(out.String(), adapter) {
		t.Errorf("output still contains adapter")
	}
}

func TestTrimPairedEndWithSingleOutput(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	// pair 1: both reads survive; pair 2: read2 trimmed below min length -> goes to singles
	in1 := makeFASTQ(
		[3]string{"p1", "ACGTACGTACGTACGTACGT" + adapter, qstring(20+len(adapter), 'I')},
		[3]string{"p2", "ACGTACGTACGTACGTACGT", qstring(20, 'I')},
	)
	in2 := makeFASTQ(
		[3]string{"p1", "TGCATGCATGCATGCATGCA" + adapter, qstring(20+len(adapter), 'I')},
		[3]string{"p2", "AAA" + adapter, qstring(3+len(adapter), 'I')},
	)
	var o1, o2, os bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10

	stats, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, &os, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if stats.TotalReads != 4 {
		t.Errorf("TotalReads = %d, want 4", stats.TotalReads)
	}
	if stats.DiscardedReads != 1 {
		t.Errorf("DiscardedReads = %d, want 1", stats.DiscardedReads)
	}
	// Only pair p1 has both mates passing, so the paired outputs contain p1
	// only; p2's surviving read1 is written to the singles file.
	if !strings.Contains(o1.String(), "@p1\n") || strings.Contains(o1.String(), "@p2\n") {
		t.Errorf("o1 contents wrong: %q", o1.String())
	}
	if !strings.Contains(o2.String(), "@p1\n") || strings.Contains(o2.String(), "@p2\n") {
		t.Errorf("o2 contents wrong: %q", o2.String())
	}
	// singles should contain the surviving p2 read1.
	if !strings.Contains(os.String(), "@p2\n") {
		t.Errorf("singles missing p2 read1: %q", os.String())
	}
}

func TestTrimPairedEndNoSingleOutputBothDiscarded(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	in1 := makeFASTQ([3]string{"p1", "AA" + adapter, qstring(2+len(adapter), 'I')})
	in2 := makeFASTQ([3]string{"p1", "CC" + adapter, qstring(2+len(adapter), 'I')})
	var o1, o2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10

	stats, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if stats.DiscardedReads != 2 {
		t.Errorf("DiscardedReads = %d, want 2", stats.DiscardedReads)
	}
	if o1.Len() != 0 || o2.Len() != 0 {
		t.Errorf("expected empty paired outputs")
	}
}

func TestTrimPairedEndMismatchedLengths(t *testing.T) {
	// Second input shorter than first -> loop stops at EOF for reader2.
	in1 := makeFASTQ(
		[3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')},
		[3]string{"p2", "ACGTACGTACGTACGTACGT", qstring(20, 'I')},
	)
	in2 := makeFASTQ([3]string{"p1", "TGCATGCATGCATGCATGCA", qstring(20, 'I')})
	var o1, o2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.MinLength = 10
	stats, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if stats.TotalReads != 2 {
		t.Errorf("TotalReads = %d, want 2 (only first pair processed)", stats.TotalReads)
	}
}

func TestTrimPairedEndWithUMI(t *testing.T) {
	in1 := makeFASTQ([3]string{"p1", "AAAACCCCACGTACGTACGTACGT", qstring(24, 'I')})
	in2 := makeFASTQ([3]string{"p1", "TGCATGCATGCATGCATGCA", qstring(20, 'I')})
	var o1, o2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.UMILength = 8
	opts.UMIPosition = "5prime"
	opts.MinLength = 10

	stats, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if stats.UMIStats == nil || stats.UMIStats.TotalUMIs != 1 {
		t.Fatalf("expected 1 UMI, got %+v", stats.UMIStats)
	}
	if !strings.Contains(o1.String(), "UMI:AAAACCCC") {
		t.Errorf("o1 missing UMI tag: %q", o1.String())
	}
	if !strings.Contains(o2.String(), "UMI:AAAACCCC") {
		t.Errorf("o2 missing UMI tag: %q", o2.String())
	}
}

func TestTrimSingleEndPhred64(t *testing.T) {
	// 'h' = Phred64 quality 40.
	seq := "ACGTACGTACGTACGTACGT"
	in := makeFASTQ([3]string{"r1", seq, qstring(20, 'h')})
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.MinLength = 10
	if _, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred64, opts); err != nil {
		t.Fatalf("TrimSingleEnd Phred64: %v", err)
	}
	if !strings.Contains(out.String(), "@r1\n"+seq+"\n") {
		t.Errorf("expected unchanged read, got %q", out.String())
	}
}

func TestTrimSingleEndBadFASTQ(t *testing.T) {
	// Malformed FASTQ should produce a read error.
	in := "@bad\nACGT\n" // truncated record
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	_, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err == nil {
		t.Errorf("expected error for malformed FASTQ")
	}
}

func TestExtractUMIDisabledAndTooLong(t *testing.T) {
	rec := &fastq.Record{ID: "r", Description: "r", Sequence: []byte("ACGTACGT"), Quality: []byte("IIIIIIII")}
	// disabled
	if _, umi := extractUMI(rec, TrimOptions{UMILength: 0}); umi != "" {
		t.Errorf("expected no UMI when disabled")
	}
	// too long
	if _, umi := extractUMI(rec, TrimOptions{UMILength: 100, UMIPosition: "5prime"}); umi != "" {
		t.Errorf("expected no UMI when length >= read length")
	}
}

func TestProcessBatchSingleEnd(t *testing.T) {
	dir := t.TempDir()
	adapter := "AGATCGGAAGAGC"
	content := makeFASTQ(
		[3]string{"r1", "ACGTACGTACGTACGTACGT" + adapter, qstring(20+len(adapter), 'I')},
		[3]string{"r2", "ACGTACGTACGTACGTACGT", qstring(20, 'I')},
	)
	var jobs []BatchJob
	for i := 0; i < 3; i++ {
		in := filepath.Join(dir, "in"+string(rune('0'+i))+".fastq")
		if err := os.WriteFile(in, []byte(content), 0644); err != nil {
			t.Fatal(err)
		}
		jobs = append(jobs, BatchJob{InputFile: in, OutputFile: filepath.Join(dir, "out"+string(rune('0'+i))+".fastq"), Index: i})
	}
	// Add one job pointing at a missing file.
	jobs = append(jobs, BatchJob{InputFile: filepath.Join(dir, "missing.fastq"), OutputFile: filepath.Join(dir, "out_missing.fastq"), Index: 99})

	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10

	results, err := ProcessBatch(jobs, fastq.Phred33, opts, 2)
	if err != nil {
		t.Fatalf("ProcessBatch: %v", err)
	}
	if len(results) != 4 {
		t.Fatalf("got %d results, want 4", len(results))
	}
	var ok, fail int
	for _, r := range results {
		if r.Error != nil {
			fail++
		} else {
			ok++
			if r.Stats.TotalReads != 2 {
				t.Errorf("expected 2 reads, got %d", r.Stats.TotalReads)
			}
		}
	}
	if ok != 3 || fail != 1 {
		t.Errorf("ok=%d fail=%d, want 3/1", ok, fail)
	}
	// Verify one output file content.
	data, err := os.ReadFile(filepath.Join(dir, "out0.fastq"))
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(data), adapter) {
		t.Errorf("output still has adapter: %q", string(data))
	}
}

func TestProcessBatchZeroWorkers(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte(makeFASTQ([3]string{"r1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	jobs := []BatchJob{{InputFile: in, OutputFile: filepath.Join(dir, "out.fastq")}}
	results, err := ProcessBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 0)
	if err != nil {
		t.Fatalf("ProcessBatch: %v", err)
	}
	if len(results) != 1 || results[0].Error != nil {
		t.Fatalf("unexpected result: %+v", results)
	}
}

func TestProcessBatchGzip(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq.gz")
	f, err := os.Create(in)
	if err != nil {
		t.Fatal(err)
	}
	gw := gzip.NewWriter(f)
	if _, err := gw.Write([]byte(makeFASTQ([3]string{"r1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')}))); err != nil {
		t.Fatal(err)
	}
	gw.Close()
	f.Close()

	out := filepath.Join(dir, "out.fastq.gz")
	jobs := []BatchJob{{InputFile: in, OutputFile: out}}
	opts := DefaultTrimOptions()
	opts.MinLength = 10
	results, err := ProcessBatch(jobs, fastq.Phred33, opts, 1)
	if err != nil {
		t.Fatalf("ProcessBatch: %v", err)
	}
	if results[0].Error != nil {
		t.Fatalf("job failed: %v", results[0].Error)
	}
	if results[0].Stats.TotalReads != 1 {
		t.Errorf("expected 1 read, got %d", results[0].Stats.TotalReads)
	}
	// Read gzip output back.
	of, err := os.Open(out)
	if err != nil {
		t.Fatal(err)
	}
	defer of.Close()
	gr, err := gzip.NewReader(of)
	if err != nil {
		t.Fatal(err)
	}
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(gr); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(buf.String(), "@r1\n") {
		t.Errorf("gzip output missing read: %q", buf.String())
	}
}

func TestProcessPairedBatch(t *testing.T) {
	dir := t.TempDir()
	adapter := "AGATCGGAAGAGC"
	c1 := makeFASTQ(
		[3]string{"p1", "ACGTACGTACGTACGTACGT" + adapter, qstring(20+len(adapter), 'I')},
		[3]string{"p2", "ACGTACGTACGTACGTACGT", qstring(20, 'I')},
	)
	c2 := makeFASTQ(
		[3]string{"p1", "TGCATGCATGCATGCATGCA" + adapter, qstring(20+len(adapter), 'I')},
		[3]string{"p2", "TT" + adapter, qstring(2+len(adapter), 'I')},
	)
	in1 := filepath.Join(dir, "r1.fastq")
	in2 := filepath.Join(dir, "r2.fastq")
	os.WriteFile(in1, []byte(c1), 0644)
	os.WriteFile(in2, []byte(c2), 0644)

	jobs := []BatchPairedJob{
		{
			InputFile1:   in1,
			InputFile2:   in2,
			OutputFile1:  filepath.Join(dir, "o1.fastq"),
			OutputFile2:  filepath.Join(dir, "o2.fastq"),
			OutputSingle: filepath.Join(dir, "single.fastq"),
			Index:        0,
		},
		// failing job: missing input1
		{
			InputFile1:  filepath.Join(dir, "nope1.fastq"),
			InputFile2:  in2,
			OutputFile1: filepath.Join(dir, "x1.fastq"),
			OutputFile2: filepath.Join(dir, "x2.fastq"),
			Index:       1,
		},
		// failing job: missing input2
		{
			InputFile1:  in1,
			InputFile2:  filepath.Join(dir, "nope2.fastq"),
			OutputFile1: filepath.Join(dir, "y1.fastq"),
			OutputFile2: filepath.Join(dir, "y2.fastq"),
			Index:       2,
		},
	}
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10

	results, err := ProcessPairedBatch(jobs, fastq.Phred33, opts, 2)
	if err != nil {
		t.Fatalf("ProcessPairedBatch: %v", err)
	}
	if len(results) != 3 {
		t.Fatalf("got %d results, want 3", len(results))
	}
	var ok, fail int
	for _, r := range results {
		if r.Error != nil {
			fail++
		} else {
			ok++
			if r.Stats.TotalReads != 4 {
				t.Errorf("expected 4 reads, got %d", r.Stats.TotalReads)
			}
		}
	}
	if ok != 1 || fail != 2 {
		t.Errorf("ok=%d fail=%d, want 1/2", ok, fail)
	}
	// single output should have p2 read1 survivor.
	data, err := os.ReadFile(filepath.Join(dir, "single.fastq"))
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(data), "@p2\n") {
		t.Errorf("single output missing survivor: %q", string(data))
	}
}

func TestProcessPairedBatchZeroWorkers(t *testing.T) {
	dir := t.TempDir()
	in1 := filepath.Join(dir, "r1.fastq")
	in2 := filepath.Join(dir, "r2.fastq")
	os.WriteFile(in1, []byte(makeFASTQ([3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})), 0644)
	os.WriteFile(in2, []byte(makeFASTQ([3]string{"p1", "TGCATGCATGCATGCATGCA", qstring(20, 'I')})), 0644)
	jobs := []BatchPairedJob{{
		InputFile1:  in1,
		InputFile2:  in2,
		OutputFile1: filepath.Join(dir, "o1.fastq"),
		OutputFile2: filepath.Join(dir, "o2.fastq"),
	}}
	results, err := ProcessPairedBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 0)
	if err != nil {
		t.Fatalf("ProcessPairedBatch: %v", err)
	}
	if len(results) != 1 || results[0].Error != nil {
		t.Fatalf("unexpected result: %+v", results)
	}
}

func TestToJSONWithUMIStats(t *testing.T) {
	s := &TrimStats{
		TotalReads: 10,
		UMIStats:   &UMIStats{TotalUMIs: 10, UniqueUMIs: 3, UMIDistribution: map[string]int{"AAAA": 5}},
	}
	out := s.ToJSON()
	if !strings.Contains(out, "umi_stats") || !strings.Contains(out, "unique_umis") {
		t.Errorf("JSON missing UMI stats: %s", out)
	}
}

func TestToHTMLWithDetectedAndUMI(t *testing.T) {
	s := &TrimStats{
		TotalReads:       100,
		TrimmedReads:     50,
		AdapterFound3:    50,
		AdapterFound5:    5,
		DiscardedReads:   2,
		TotalBases:       2000,
		TrimmedBases:     400,
		DetectedAdapter3: "AGATCGGAAGAGC",
		DetectedAdapter5: "GTTCAGAGTTCTACAGTCCGACGATC",
		UMIStats:         &UMIStats{TotalUMIs: 98, UniqueUMIs: 40},
	}
	html := s.ToHTML()
	for _, want := range []string{"Detected 3' Adapter", "Detected 5' Adapter", "Total UMIs", "Unique UMIs", "AGATCGGAAGAGC"} {
		if !strings.Contains(html, want) {
			t.Errorf("HTML missing %q", want)
		}
	}
}

func TestTrimSingleEndProgressReporting(t *testing.T) {
	var recs [][3]string
	for i := 0; i < 5; i++ {
		recs = append(recs, [3]string{"r", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
	}
	in := makeFASTQ(recs...)
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.ProgressReport = true
	opts.ProgressInterval = 2
	opts.MinLength = 10
	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.TotalReads != 5 {
		t.Errorf("TotalReads = %d, want 5", stats.TotalReads)
	}
}

func TestTrimRecordAdapterShorterThanMinOverlap(t *testing.T) {
	// minOverlap larger than sequence => findAdapter returns -1, no trimming.
	rec := &fastq.Record{ID: "r", Description: "r", Sequence: []byte("AGA"), Quality: []byte("III")}
	opts := DefaultTrimOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.MinOverlap = 10
	opts.MinLength = 1
	trimmed := trimRecord(rec, opts, &TrimStats{})
	if string(trimmed.Sequence) != "AGA" {
		t.Errorf("expected unchanged, got %q", string(trimmed.Sequence))
	}
}

// TestImprovedFindAdapterShorterThanMinOverlap forces compareLen < minOverlap
// on the first iteration by passing an adapter that is shorter than minOverlap
// but a sequence that is long enough for the outer loop to run. The function
// must therefore continue past every iteration and return -1.
func TestImprovedFindAdapterShorterThanMinOverlap(t *testing.T) {
	// adapter shorter than minOverlap; seq long enough for outer loop to iterate.
	got := improvedFindAdapter("AAAAAAAAAA", "AGA", 5, 0.1)
	if got != -1 {
		t.Errorf("improvedFindAdapter(adapter shorter than minOverlap) = %d, want -1", got)
	}
}

// TestDetectAdaptersFromReadsPicks5Prime verifies that when common adapters
// appear at the start (position < 10) of most reads, detectAdaptersFromReads
// reports them as the detected 5' adapter and not just the 3' adapter.
func TestDetectAdaptersFromReadsPicks5Prime(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	// Each read places the adapter at position 0 (5'-end) followed by random-ish
	// content. The same string also matches as a "3' adapter" candidate at
	// position 0 because detectAdaptersFromReads checks both ends with the same
	// scanner — that is the documented behaviour. We assert both Adapter3 and
	// Adapter5 come back populated.
	var recs []*fastq.Record
	for i := 0; i < 20; i++ {
		recs = append(recs, &fastq.Record{
			ID:       "r",
			Sequence: []byte(adapter + "TTTTTTTTTTTTTTTTTTTT"),
			Quality:  []byte(strings.Repeat("I", len(adapter)+20)),
		})
	}
	opts := detectAdaptersFromReads(recs)
	if opts.Adapter5 != adapter {
		t.Errorf("Adapter5 = %q, want %q", opts.Adapter5, adapter)
	}
	if opts.Adapter3 != adapter {
		t.Errorf("Adapter3 = %q, want %q", opts.Adapter3, adapter)
	}
}

// TestTrimSingleEndAutoDetect5Prime covers the auto-detect branch where the
// detected 5' adapter (not just 3') gets propagated into stats and opts.
func TestTrimSingleEndAutoDetect5Prime(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	var recs [][3]string
	for i := 0; i < 20; i++ {
		seq := adapter + "ACGTACGTACGTACGTACGT"
		recs = append(recs, [3]string{"r", seq, qstring(len(seq), 'I')})
	}
	in := makeFASTQ(recs...)
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.AutoDetect = true
	opts.AutoDetectReads = 30
	opts.MinLength = 5

	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd auto-detect 5': %v", err)
	}
	if stats.DetectedAdapter5 != adapter {
		t.Errorf("DetectedAdapter5 = %q, want %q", stats.DetectedAdapter5, adapter)
	}
	if stats.AdapterFound5 != 20 {
		t.Errorf("AdapterFound5 = %d, want 20", stats.AdapterFound5)
	}
}

// TestTrimSingleEndAutoDetectZeroReads exercises the maxReads<=0 default
// branch: AutoDetectReads=0 should fall back to the 1000 default and still
// auto-detect from however many reads are actually in the input.
func TestTrimSingleEndAutoDetectZeroReads(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	var recs [][3]string
	for i := 0; i < 12; i++ {
		recs = append(recs, [3]string{"r", "ACGTACGTACGTACGTACGT" + adapter, qstring(20+len(adapter), 'I')})
	}
	in := makeFASTQ(recs...)
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.AutoDetect = true
	opts.AutoDetectReads = 0 // forces default-1000 branch
	opts.MinLength = 10
	stats, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd: %v", err)
	}
	if stats.DetectedAdapter3 != adapter {
		t.Errorf("DetectedAdapter3 = %q, want %q", stats.DetectedAdapter3, adapter)
	}
	if stats.AdapterFound3 != 12 {
		t.Errorf("AdapterFound3 = %d, want 12", stats.AdapterFound3)
	}
}

// TestTrimSingleEndAutoDetectReadError ensures that a malformed record
// encountered during the auto-detect pre-read produces an error to the caller.
func TestTrimSingleEndAutoDetectReadError(t *testing.T) {
	// One full record then a truncated one — the second read errors before EOF.
	in := makeFASTQ([3]string{"r1", "ACGTACGT", qstring(8, 'I')}) + "@bad\nACGT\n"
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.AutoDetect = true
	opts.AutoDetectReads = 50
	_, err := TrimSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err == nil {
		t.Errorf("expected error for malformed FASTQ during auto-detect")
	} else if !strings.Contains(err.Error(), "error reading FASTQ") {
		t.Errorf("unexpected error message: %v", err)
	}
}

// TestDetectAdaptersReadError exercises detectAdapters' error path: a malformed
// FASTQ stream should surface the read error rather than silently returning
// empty options.
func TestDetectAdaptersReadError(t *testing.T) {
	in := "@bad\nACGT\n" // truncated record
	reader := fastq.NewReader(strings.NewReader(in), fastq.Phred33)
	_, err := detectAdapters(reader, 100)
	if err == nil {
		t.Errorf("expected error for malformed FASTQ")
	}
}

// TestTrimSingleEndWriteError covers the writer.Write error path. We make the
// underlying writer fail after the first byte, then write a single record
// whose total serialized size exceeds the default 4 KiB bufio.Writer buffer
// so that the buffered Write attempts a flush and returns the error to skewer.
func TestTrimSingleEndWriteError(t *testing.T) {
	seqLen := 8 * 1024 // 8 KiB > default bufio.Writer buffer (4 KiB)
	seq := bigSeq(seqLen)
	in := makeFASTQ([3]string{"r1", seq, qstring(seqLen, 'I')})
	w := &errWriter{afterBytes: 16, errAfter: errors.New("disk full")}
	opts := DefaultTrimOptions()
	opts.MinLength = 5
	_, err := TrimSingleEnd(strings.NewReader(in), w, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
	if !strings.Contains(err.Error(), "error writing FASTQ") {
		t.Errorf("expected wrapped write error, got %v", err)
	}
}

// TestTrimSingleEndFlushError covers the writer.Flush error path. A single
// short record fits entirely in the bufio buffer, so the error only surfaces
// when skewer flushes the writer at the end of the run.
func TestTrimSingleEndFlushError(t *testing.T) {
	in := makeFASTQ([3]string{"r1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
	w := &errWriter{afterBytes: 0, errAfter: errors.New("flush boom")}
	opts := DefaultTrimOptions()
	opts.MinLength = 5
	_, err := TrimSingleEnd(strings.NewReader(in), w, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected flush error, got nil")
	}
	if !strings.Contains(err.Error(), "error flushing output") {
		t.Errorf("expected wrapped flush error, got %v", err)
	}
}

// TestTrimPairedEndReader1Error ensures that an error from reader1 (other than
// EOF) is propagated. We give reader1 a malformed record but make sure reader2
// still has at least one good record so the EOF short-circuit at line ~327
// does not fire before err1 is checked.
func TestTrimPairedEndReader1Error(t *testing.T) {
	in1 := "@bad\nACGT\n"                                                       // malformed -> non-EOF error
	in2 := makeFASTQ([3]string{"p1", "TGCATGCATGCATGCATGCA", qstring(20, 'I')}) // valid record
	var o1, o2 bytes.Buffer
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, nil, fastq.Phred33, DefaultTrimOptions())
	if err == nil {
		t.Fatalf("expected reader1 error")
	}
	if !strings.Contains(err.Error(), "error reading first input") {
		t.Errorf("expected wrapped reader1 error, got %v", err)
	}
}

// TestTrimPairedEndReader2Error mirrors TestTrimPairedEndReader1Error for
// reader2: reader1 must yield a valid first record, then reader2's malformed
// record is what surfaces.
func TestTrimPairedEndReader2Error(t *testing.T) {
	in1 := makeFASTQ([3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
	in2 := "@bad\nTGCA\n"
	var o1, o2 bytes.Buffer
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, nil, fastq.Phred33, DefaultTrimOptions())
	if err == nil {
		t.Fatalf("expected reader2 error")
	}
	if !strings.Contains(err.Error(), "error reading second input") {
		t.Errorf("expected wrapped reader2 error, got %v", err)
	}
}

// TestTrimPairedEndWriter1Error exercises the writer1.Write error path. The
// pair-aware loop only writes to writer1 when both reads pass MinLength, so
// the seqs must be long enough to survive trimming and large enough to push
// through the bufio buffer.
func TestTrimPairedEndWriter1Error(t *testing.T) {
	seq := bigSeq(8 * 1024)
	in1 := makeFASTQ([3]string{"p1", seq, qstring(len(seq), 'I')})
	in2 := makeFASTQ([3]string{"p1", seq, qstring(len(seq), 'I')})
	w1 := &errWriter{afterBytes: 16, errAfter: errors.New("w1 fail")}
	var o2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.MinLength = 5
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), w1, &o2, nil, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected writer1 error")
	}
	if !strings.Contains(err.Error(), "error writing first output") {
		t.Errorf("expected wrapped writer1 error, got %v", err)
	}
}

// TestTrimPairedEndWriter2Error mirrors TestTrimPairedEndWriter1Error for
// writer2. writer1 must succeed so we reach the writer2.Write call.
func TestTrimPairedEndWriter2Error(t *testing.T) {
	seq := bigSeq(8 * 1024)
	in1 := makeFASTQ([3]string{"p1", seq, qstring(len(seq), 'I')})
	in2 := makeFASTQ([3]string{"p1", seq, qstring(len(seq), 'I')})
	var o1 bytes.Buffer
	w2 := &errWriter{afterBytes: 16, errAfter: errors.New("w2 fail")}
	opts := DefaultTrimOptions()
	opts.MinLength = 5
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, w2, nil, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected writer2 error")
	}
	if !strings.Contains(err.Error(), "error writing second output") {
		t.Errorf("expected wrapped writer2 error, got %v", err)
	}
}

// TestTrimPairedEndOnlyRead1Survives covers the (!pass1 -> DiscardedReads) +
// (pass2 -> writerSingle.Write) branches in the paired loop's singles handler.
// pair p1 has read1 too short after trimming and read2 long enough to survive.
func TestTrimPairedEndOnlyRead2Survives(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	in1 := makeFASTQ([3]string{"p1", "AA" + adapter, qstring(2+len(adapter), 'I')}) // read1 trimmed to 2bp, < MinLength
	in2 := makeFASTQ([3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})     // read2 untouched, 20bp
	var o1, o2, os bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10
	stats, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, &os, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if stats.DiscardedReads != 1 {
		t.Errorf("DiscardedReads = %d, want 1 (read1 discarded)", stats.DiscardedReads)
	}
	if o1.Len() != 0 || o2.Len() != 0 {
		t.Errorf("paired outputs should be empty, got o1=%q o2=%q", o1.String(), o2.String())
	}
	if !strings.Contains(os.String(), "@p1\n") {
		t.Errorf("singles should contain surviving read2: %q", os.String())
	}
	if !strings.Contains(os.String(), "ACGTACGTACGTACGTACGT") {
		t.Errorf("singles content wrong: %q", os.String())
	}
}

// TestTrimPairedEndSingleWriteErrorPass1 covers the writerSingle.Write error
// path for the surviving read1 (pass1, !pass2).
func TestTrimPairedEndSingleWriteErrorPass1(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	seq := bigSeq(8 * 1024)
	in1 := makeFASTQ([3]string{"p1", seq, qstring(len(seq), 'I')})                  // pass1
	in2 := makeFASTQ([3]string{"p1", "AA" + adapter, qstring(2+len(adapter), 'I')}) // !pass2
	var o1, o2 bytes.Buffer
	ws := &errWriter{afterBytes: 16, errAfter: errors.New("single fail")}
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, ws, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected single-output write error")
	}
	if !strings.Contains(err.Error(), "error writing single output") {
		t.Errorf("expected wrapped single write error, got %v", err)
	}
}

// TestTrimPairedEndSingleWriteErrorPass2 covers the writerSingle.Write error
// path for the surviving read2 (!pass1, pass2).
func TestTrimPairedEndSingleWriteErrorPass2(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	seq := bigSeq(8 * 1024)
	in1 := makeFASTQ([3]string{"p1", "AA" + adapter, qstring(2+len(adapter), 'I')}) // !pass1
	in2 := makeFASTQ([3]string{"p1", seq, qstring(len(seq), 'I')})                  // pass2
	var o1, o2 bytes.Buffer
	ws := &errWriter{afterBytes: 16, errAfter: errors.New("single fail2")}
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, ws, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected single-output write error for pass2")
	}
	if !strings.Contains(err.Error(), "error writing single output") {
		t.Errorf("expected wrapped single write error, got %v", err)
	}
}

// TestTrimPairedEndFlushError1 covers writer1.Flush error path. errWriter
// fails immediately, so the buffered records only surface the error at flush.
func TestTrimPairedEndFlushError1(t *testing.T) {
	in1 := makeFASTQ([3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
	in2 := makeFASTQ([3]string{"p1", "TGCATGCATGCATGCATGCA", qstring(20, 'I')})
	w1 := &errWriter{afterBytes: 0, errAfter: errors.New("flush1")}
	var o2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.MinLength = 5
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), w1, &o2, nil, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected flush1 error")
	}
	if !strings.Contains(err.Error(), "error flushing first output") {
		t.Errorf("expected wrapped flush1 error, got %v", err)
	}
}

// TestTrimPairedEndFlushError2 covers writer2.Flush error path. writer1's
// flush must succeed first to reach writer2.Flush.
func TestTrimPairedEndFlushError2(t *testing.T) {
	in1 := makeFASTQ([3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
	in2 := makeFASTQ([3]string{"p1", "TGCATGCATGCATGCATGCA", qstring(20, 'I')})
	var o1 bytes.Buffer
	w2 := &errWriter{afterBytes: 0, errAfter: errors.New("flush2")}
	opts := DefaultTrimOptions()
	opts.MinLength = 5
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, w2, nil, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected flush2 error")
	}
	if !strings.Contains(err.Error(), "error flushing second output") {
		t.Errorf("expected wrapped flush2 error, got %v", err)
	}
}

// TestTrimPairedEndFlushErrorSingle covers writerSingle.Flush error path.
// One pair's read2 fails so the singles writer receives data but errors when
// flushed at the end.
func TestTrimPairedEndFlushErrorSingle(t *testing.T) {
	adapter := "AGATCGGAAGAGC"
	in1 := makeFASTQ([3]string{"p1", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
	in2 := makeFASTQ([3]string{"p1", "AA" + adapter, qstring(2+len(adapter), 'I')})
	var o1, o2 bytes.Buffer
	ws := &errWriter{afterBytes: 1 << 20, errAfter: errors.New("flushsingle")}
	// afterBytes huge => buffered Write succeeds, only Flush surfaces error...
	// but we actually need flush to fail. Use a writer that fails after enough
	// bytes to buffer the record but fail on the trailing flush. The fastq
	// bufio buffer is 4 KiB and the record is ~36 bytes, so we set afterBytes
	// to 30 so the first chunk goes through then the flush at end fails.
	ws = &errWriter{afterBytes: 30, errAfter: errors.New("flushsingle")}
	opts := DefaultTrimOptions()
	opts.Adapter3 = adapter
	opts.MinLength = 10
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, ws, fastq.Phred33, opts)
	if err == nil {
		t.Fatalf("expected flush-single error")
	}
	if !strings.Contains(err.Error(), "error flushing single output") {
		t.Errorf("expected wrapped flush-single error, got %v", err)
	}
}

// TestTrimPairedEndProgressReporting exercises the paired-end progress
// reporting branch (ProgressReport=true, readCount % ProgressInterval == 0).
func TestTrimPairedEndProgressReporting(t *testing.T) {
	var recs1, recs2 [][3]string
	for i := 0; i < 4; i++ {
		recs1 = append(recs1, [3]string{"p", "ACGTACGTACGTACGTACGT", qstring(20, 'I')})
		recs2 = append(recs2, [3]string{"p", "TGCATGCATGCATGCATGCA", qstring(20, 'I')})
	}
	in1 := makeFASTQ(recs1...)
	in2 := makeFASTQ(recs2...)
	var o1, o2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.ProgressReport = true
	opts.ProgressInterval = 2 // readCount increments by 2 per pair
	opts.MinLength = 5
	stats, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &o1, &o2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if stats.TotalReads != 8 {
		t.Errorf("TotalReads = %d, want 8", stats.TotalReads)
	}
}

// TestProcessSingleJobOpenWriterError forces iohelper.OpenWriter to fail by
// pointing the output at a path whose parent directory does not exist.
func TestProcessSingleJobOpenWriterError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	// /nonexistent_xyz_skewer/... will fail to create; relies on no such path.
	bad := filepath.Join(dir, "no", "such", "dir", "out.fastq")
	jobs := []BatchJob{{InputFile: in, OutputFile: bad}}
	results, err := ProcessBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 1)
	if err != nil {
		t.Fatalf("ProcessBatch: %v", err)
	}
	if len(results) != 1 || results[0].Error == nil {
		t.Fatalf("expected job error, got %+v", results)
	}
	if !strings.Contains(results[0].Error.Error(), "error creating output file") {
		t.Errorf("expected output-create error, got %v", results[0].Error)
	}
}

// TestProcessSingleJobTrimError feeds processSingleJob a malformed FASTQ so
// TrimSingleEnd returns an error and that error is reported through the
// BatchResult instead of as the OpenWriter error.
func TestProcessSingleJobTrimError(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte("@bad\nACGT\n"), 0644); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(dir, "out.fastq")
	jobs := []BatchJob{{InputFile: in, OutputFile: out}}
	results, err := ProcessBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 1)
	if err != nil {
		t.Fatalf("ProcessBatch: %v", err)
	}
	if results[0].Error == nil {
		t.Fatalf("expected processing error")
	}
	if !strings.Contains(results[0].Error.Error(), "error processing file") {
		t.Errorf("expected processing-file error, got %v", results[0].Error)
	}
}

// TestProcessPairedJobOpenWriter1Error covers the OpenWriter(output1) error in
// processPairedJob.
func TestProcessPairedJobOpenWriter1Error(t *testing.T) {
	dir := t.TempDir()
	in1 := filepath.Join(dir, "in1.fastq")
	in2 := filepath.Join(dir, "in2.fastq")
	if err := os.WriteFile(in1, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(in2, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	jobs := []BatchPairedJob{{
		InputFile1:  in1,
		InputFile2:  in2,
		OutputFile1: filepath.Join(dir, "no", "such", "dir", "o1.fastq"),
		OutputFile2: filepath.Join(dir, "o2.fastq"),
	}}
	results, err := ProcessPairedBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 1)
	if err != nil {
		t.Fatalf("ProcessPairedBatch: %v", err)
	}
	if results[0].Error == nil {
		t.Fatalf("expected job error, got %+v", results[0])
	}
	if !strings.Contains(results[0].Error.Error(), "error creating output file") ||
		!strings.Contains(results[0].Error.Error(), "o1.fastq") {
		t.Errorf("expected wrapped o1 create error, got %v", results[0].Error)
	}
}

// TestProcessPairedJobOpenWriter2Error: output1 is fine, output2 fails to
// create.
func TestProcessPairedJobOpenWriter2Error(t *testing.T) {
	dir := t.TempDir()
	in1 := filepath.Join(dir, "in1.fastq")
	in2 := filepath.Join(dir, "in2.fastq")
	if err := os.WriteFile(in1, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(in2, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	jobs := []BatchPairedJob{{
		InputFile1:  in1,
		InputFile2:  in2,
		OutputFile1: filepath.Join(dir, "o1.fastq"),
		OutputFile2: filepath.Join(dir, "no", "such", "dir", "o2.fastq"),
	}}
	results, err := ProcessPairedBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 1)
	if err != nil {
		t.Fatalf("ProcessPairedBatch: %v", err)
	}
	if results[0].Error == nil {
		t.Fatalf("expected job error")
	}
	if !strings.Contains(results[0].Error.Error(), "o2.fastq") {
		t.Errorf("expected o2 create error, got %v", results[0].Error)
	}
}

// TestProcessPairedJobOpenSingleError: output1+2 fine, OutputSingle fails.
func TestProcessPairedJobOpenSingleError(t *testing.T) {
	dir := t.TempDir()
	in1 := filepath.Join(dir, "in1.fastq")
	in2 := filepath.Join(dir, "in2.fastq")
	if err := os.WriteFile(in1, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(in2, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	jobs := []BatchPairedJob{{
		InputFile1:   in1,
		InputFile2:   in2,
		OutputFile1:  filepath.Join(dir, "o1.fastq"),
		OutputFile2:  filepath.Join(dir, "o2.fastq"),
		OutputSingle: filepath.Join(dir, "no", "such", "dir", "single.fastq"),
	}}
	results, err := ProcessPairedBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 1)
	if err != nil {
		t.Fatalf("ProcessPairedBatch: %v", err)
	}
	if results[0].Error == nil {
		t.Fatalf("expected job error")
	}
	if !strings.Contains(results[0].Error.Error(), "single output file") {
		t.Errorf("expected single-output create error, got %v", results[0].Error)
	}
}

// TestProcessPairedJobTrimError forces TrimPairedEnd to fail (malformed input)
// after all files have been opened, covering the error wrap on line ~211.
func TestProcessPairedJobTrimError(t *testing.T) {
	dir := t.TempDir()
	in1 := filepath.Join(dir, "in1.fastq")
	in2 := filepath.Join(dir, "in2.fastq")
	// Both readers need at least one bad record AND the partner must yield a
	// valid record for the error to surface (otherwise the EOF short-circuit
	// fires first). Easiest: identical malformed input on both sides so the
	// first reader hits a non-EOF parse error.
	if err := os.WriteFile(in1, []byte("@bad\nACGT\n+\nXY\n"), 0644); err != nil { // length mismatch
		t.Fatal(err)
	}
	if err := os.WriteFile(in2, []byte(makeFASTQ([3]string{"r", "ACGT", qstring(4, 'I')})), 0644); err != nil {
		t.Fatal(err)
	}
	jobs := []BatchPairedJob{{
		InputFile1:  in1,
		InputFile2:  in2,
		OutputFile1: filepath.Join(dir, "o1.fastq"),
		OutputFile2: filepath.Join(dir, "o2.fastq"),
	}}
	results, err := ProcessPairedBatch(jobs, fastq.Phred33, DefaultTrimOptions(), 1)
	if err != nil {
		t.Fatalf("ProcessPairedBatch: %v", err)
	}
	if results[0].Error == nil {
		t.Fatalf("expected processing error")
	}
	if !strings.Contains(results[0].Error.Error(), "error processing files") {
		t.Errorf("expected processing-files error, got %v", results[0].Error)
	}
}

// Sanity: confirm that errWriter actually returns errors when expected so
// failures in the I/O-error tests above are not silently absorbed by the test
// framework.
func TestErrWriterBehaviour(t *testing.T) {
	w := &errWriter{afterBytes: 4, errAfter: errors.New("after4")}
	n, err := w.Write([]byte("abcdef"))
	if n != 4 || err == nil {
		t.Fatalf("Write = (%d,%v), want (4, error)", n, err)
	}
	_, err = w.Write([]byte("x"))
	if err == nil {
		t.Fatalf("second Write should also return error")
	}
	// io.Writer interface conformance check.
	var _ io.Writer = w
}
