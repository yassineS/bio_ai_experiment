package skewer

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

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
	out, err := s.ToJSON()
	if err != nil {
		t.Fatalf("ToJSON: %v", err)
	}
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
