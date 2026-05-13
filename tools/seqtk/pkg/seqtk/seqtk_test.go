package seqtk

import (
	"bytes"
	"compress/gzip"
	"os"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestCalculateFastaStats(t *testing.T) {
	fasta := `>seq1
ACGTACGT
>seq2
GCGCGCGCGC
>seq3
ATATATAT
`
	r := strings.NewReader(fasta)
	stats, err := CalculateFastaStats(r)
	if err != nil {
		t.Fatalf("CalculateFastaStats failed: %v", err)
	}

	if stats.NumSequences != 3 {
		t.Errorf("Expected 3 sequences, got %d", stats.NumSequences)
	}

	if stats.TotalBases != 26 {
		t.Errorf("Expected 26 total bases, got %d", stats.TotalBases)
	}

	if stats.MinLength != 8 {
		t.Errorf("Expected min length 8, got %d", stats.MinLength)
	}

	if stats.MaxLength != 10 {
		t.Errorf("Expected max length 10, got %d", stats.MaxLength)
	}

	expectedGC := (10.0 + 4.0) / 26.0 * 100.0 // 10 GC from seq2, 4 from seq1
	if diff := stats.GCContent - expectedGC; diff > 0.01 || diff < -0.01 {
		t.Errorf("Expected GC content %.2f%%, got %.2f%%", expectedGC, stats.GCContent)
	}
}

func TestCalculateFastqStats(t *testing.T) {
	fastqData := `@read1
ACGTACGT
+
IIIIIIII
@read2
GCGCGCGC
+
IIIIIIII
`
	r := strings.NewReader(fastqData)
	stats, err := CalculateFastqStats(r, fastq.Phred33)
	if err != nil {
		t.Fatalf("CalculateFastqStats failed: %v", err)
	}

	if stats.NumSequences != 2 {
		t.Errorf("Expected 2 sequences, got %d", stats.NumSequences)
	}

	if stats.TotalBases != 16 {
		t.Errorf("Expected 16 total bases, got %d", stats.TotalBases)
	}
}

func TestConvertFastqToFasta(t *testing.T) {
	fastqData := `@read1 description
ACGTACGT
+
IIIIIIII
@read2
GCGCGCGC
+
IIIIIIII
`
	r := strings.NewReader(fastqData)
	var buf bytes.Buffer

	err := ConvertFastqToFasta(r, &buf, fastq.Phred33)
	if err != nil {
		t.Fatalf("ConvertFastqToFasta failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, ">read1 description") {
		t.Error("Output doesn't contain expected header")
	}
	if !strings.Contains(output, "ACGTACGT") {
		t.Error("Output doesn't contain expected sequence")
	}
	if strings.Contains(output, "IIIIIIII") {
		t.Error("Output shouldn't contain quality scores")
	}
}

func TestReverseComplementFasta(t *testing.T) {
	fasta := `>seq1
ACGT
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	err := ReverseComplement(r, &buf, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("ReverseComplement failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "ACGT") {
		t.Error("Expected reverse complement ACGT")
	}
}

func TestSampleFasta(t *testing.T) {
	fasta := `>seq1
ACGT
>seq2
GCTA
>seq3
TGCA
>seq4
CATG
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	err := Sample(r, &buf, 0.5, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("Sample failed: %v", err)
	}

	output := buf.String()
	// Should contain at least one sequence
	if !strings.Contains(output, ">seq") {
		t.Error("Output should contain at least one sequence")
	}
}

func TestTrimQuality(t *testing.T) {
	// Quality scores: ! = 0, I = 40 (Phred+33)
	fastqData := `@read1
ACGTACGT
+
!!!!IIII
`
	r := strings.NewReader(fastqData)
	var buf bytes.Buffer

	err := TrimQuality(r, &buf, 30, fastq.Phred33)
	if err != nil {
		t.Fatalf("TrimQuality failed: %v", err)
	}

	output := buf.String()
	// Should only keep high-quality portion
	if !strings.Contains(output, "@read1") {
		t.Error("Output should contain trimmed read")
	}
}

func TestGetFileType(t *testing.T) {
	// Create temp FASTA file
	fastaContent := []byte(">seq1\nACGT\n")
	fastaFile := "/tmp/test.fasta"
	if err := os.WriteFile(fastaFile, fastaContent, 0644); err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	defer os.Remove(fastaFile)

	isFastq, err := GetFileType(fastaFile)
	if err != nil {
		t.Fatalf("GetFileType failed: %v", err)
	}
	if isFastq {
		t.Error("FASTA file detected as FASTQ")
	}

	// Create temp FASTQ file
	fastqContent := []byte("@seq1\nACGT\n+\nIIII\n")
	fastqFile := "/tmp/test.fastq"
	if err := os.WriteFile(fastqFile, fastqContent, 0644); err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	defer os.Remove(fastqFile)

	isFastq, err = GetFileType(fastqFile)
	if err != nil {
		t.Fatalf("GetFileType failed: %v", err)
	}
	if !isFastq {
		t.Error("FASTQ file not detected as FASTQ")
	}
}

func TestDecompressReader(t *testing.T) {
	testData := "test data content"

	// Test gzip
	var gzBuf bytes.Buffer
	gzWriter := gzip.NewWriter(&gzBuf)
	gzWriter.Write([]byte(testData))
	gzWriter.Close()

	reader, err := DecompressReader(&gzBuf, "test.gz")
	if err != nil {
		t.Fatalf("DecompressReader failed for gzip: %v", err)
	}

	var output bytes.Buffer
	output.ReadFrom(reader)
	if output.String() != testData {
		t.Errorf("Decompressed data doesn't match: got %q, want %q", output.String(), testData)
	}

	// Test plain file
	plainReader := strings.NewReader(testData)
	reader, err = DecompressReader(plainReader, "test.txt")
	if err != nil {
		t.Fatalf("DecompressReader failed for plain file: %v", err)
	}

	output.Reset()
	output.ReadFrom(reader)
	if output.String() != testData {
		t.Errorf("Plain data doesn't match: got %q, want %q", output.String(), testData)
	}
}

func TestCompressWriter(t *testing.T) {
	testData := "test data content"

	// Test gzip
	var gzBuf bytes.Buffer
	writer, err := CompressWriter(&gzBuf, "test.gz")
	if err != nil {
		t.Fatalf("CompressWriter failed for gzip: %v", err)
	}

	writer.Write([]byte(testData))
	writer.Close()

	// Verify compressed data can be read back
	reader, err := gzip.NewReader(&gzBuf)
	if err != nil {
		t.Fatalf("Failed to create gzip reader: %v", err)
	}

	var output bytes.Buffer
	output.ReadFrom(reader)
	if output.String() != testData {
		t.Errorf("Compressed/decompressed data doesn't match: got %q, want %q", output.String(), testData)
	}

	// Test plain file - CompressWriter returns nil for non-compressed files
	var plainBuf bytes.Buffer
	writer, err = CompressWriter(&plainBuf, "test.txt")
	if err != nil {
		t.Fatalf("CompressWriter failed for plain file: %v", err)
	}

	// For plain files, writer should be nil
	if writer != nil {
		t.Error("CompressWriter should return nil for plain files")
	}

	// Write directly to buffer for plain files
	plainBuf.Write([]byte(testData))

	if plainBuf.String() != testData {
		t.Errorf("Plain data doesn't match: got %q, want %q", plainBuf.String(), testData)
	}
}

func TestOpenInputOutput(t *testing.T) {
	// Test with compressed file
	testData := "@read1\nACGT\n+\nIIII\n"
	tmpFile := "/tmp/test_input.fastq.gz"

	// Create compressed test file
	file, err := os.Create(tmpFile)
	if err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}

	gzWriter := gzip.NewWriter(file)
	gzWriter.Write([]byte(testData))
	gzWriter.Close()
	file.Close()
	defer os.Remove(tmpFile)

	// Test OpenInput
	input, err := OpenInput(tmpFile)
	if err != nil {
		t.Fatalf("OpenInput failed: %v", err)
	}
	defer input.Close()

	var buf bytes.Buffer
	buf.ReadFrom(input)
	if buf.String() != testData {
		t.Errorf("Input data doesn't match: got %q, want %q", buf.String(), testData)
	}

	// Test OpenOutput
	outputFile := "/tmp/test_output.fastq.gz"
	defer os.Remove(outputFile)

	output, err := OpenOutput(outputFile)
	if err != nil {
		t.Fatalf("OpenOutput failed: %v", err)
	}

	output.Write([]byte(testData))
	output.Close()

	// Verify output file
	input2, err := OpenInput(outputFile)
	if err != nil {
		t.Fatalf("Failed to reopen output file: %v", err)
	}
	defer input2.Close()

	buf.Reset()
	buf.ReadFrom(input2)
	if buf.String() != testData {
		t.Errorf("Output data doesn't match: got %q, want %q", buf.String(), testData)
	}
}

func TestGetFileTypeCompressed(t *testing.T) {
	// Create compressed FASTQ file
	fastqContent := []byte("@seq1\nACGT\n+\nIIII\n")
	fastqFile := "/tmp/test_compressed.fastq.gz"

	file, err := os.Create(fastqFile)
	if err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}

	gzWriter := gzip.NewWriter(file)
	gzWriter.Write(fastqContent)
	gzWriter.Close()
	file.Close()
	defer os.Remove(fastqFile)

	isFastq, err := GetFileType(fastqFile)
	if err != nil {
		t.Fatalf("GetFileType failed: %v", err)
	}
	if !isFastq {
		t.Error("Compressed FASTQ file not detected as FASTQ")
	}

	// Create compressed FASTA file
	fastaContent := []byte(">seq1\nACGT\n")
	fastaFile := "/tmp/test_compressed.fasta.gz"

	file, err = os.Create(fastaFile)
	if err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}

	gzWriter = gzip.NewWriter(file)
	gzWriter.Write(fastaContent)
	gzWriter.Close()
	file.Close()
	defer os.Remove(fastaFile)

	isFastq, err = GetFileType(fastaFile)
	if err != nil {
		t.Fatalf("GetFileType failed: %v", err)
	}
	if isFastq {
		t.Error("Compressed FASTA file detected as FASTQ")
	}
}

func TestFilter(t *testing.T) {
	// Test length filtering
	fasta := `>seq1
ACGT
>seq2
ACGTACGTACGT
>seq3
ACGTACGTACGTACGT
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	opts := FilterOptions{
		MinLength: 8,
		MaxLength: 12,
	}

	err := Filter(r, &buf, opts, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	output := buf.String()
	if strings.Contains(output, "seq1") {
		t.Error("seq1 (length 4) should be filtered out")
	}
	if !strings.Contains(output, "seq2") {
		t.Error("seq2 (length 12) should be included")
	}
	if strings.Contains(output, "seq3") {
		t.Error("seq3 (length 16) should be filtered out")
	}
}

func TestFilterPattern(t *testing.T) {
	fasta := `>chr1_seq1
ACGT
>chr2_seq2
GCTA
>chr1_seq3
TGCA
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	opts := FilterOptions{
		Pattern: "chr1",
	}

	err := Filter(r, &buf, opts, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "chr1_seq1") {
		t.Error("chr1_seq1 should be included")
	}
	if strings.Contains(output, "chr2_seq2") {
		t.Error("chr2_seq2 should be filtered out")
	}
	if !strings.Contains(output, "chr1_seq3") {
		t.Error("chr1_seq3 should be included")
	}
}

func TestFilterCombined(t *testing.T) {
	fasta := `>chr1_short
ACG
>chr1_long
ACGTACGTACGT
>chr2_medium
ACGTACGT
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	opts := FilterOptions{
		MinLength: 8,
		Pattern:   "chr1",
	}

	err := Filter(r, &buf, opts, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	output := buf.String()
	if strings.Contains(output, "chr1_short") {
		t.Error("chr1_short should be filtered out (too short)")
	}
	if !strings.Contains(output, "chr1_long") {
		t.Error("chr1_long should be included")
	}
	if strings.Contains(output, "chr2_medium") {
		t.Error("chr2_medium should be filtered out (wrong pattern)")
	}
}

const subseqFASTA = `>chr1 first contig
ACGTACGTAA
>chr2
TTTTGGGGCC
>chr3
NNNNN
`

func TestSubseqNameList(t *testing.T) {
	tests := []struct {
		name    string
		spec    string
		want    string
		wantNot []string
		lineLen int
	}{
		{
			name: "single name emits full record with original header",
			spec: "chr2\n",
			want: ">chr2\nTTTTGGGGCC\n",
		},
		{
			name: "in-input order, not spec order",
			// chr2 listed before chr1, but chr1 appears first in the input.
			spec:    "chr2\nchr1\n",
			want:    ">chr1 first contig\nACGTACGTAA\n>chr2\nTTTTGGGGCC\n",
			wantNot: []string{">chr2\nTTTTGGGGCC\n>chr1"},
		},
		{
			name: "trailing whitespace after name is ignored",
			spec: "chr2\textra columns here\n",
			want: ">chr2\nTTTTGGGGCC\n",
		},
		{
			name:    "unknown name warns and is skipped",
			spec:    "chrX\nchr3\n",
			want:    ">chr3\nNNNNN\n",
			wantNot: []string{"chrX"},
		},
		{
			name:    "line wrapping",
			spec:    "chr1\n",
			lineLen: 4,
			want:    ">chr1 first contig\nACGT\nACGT\nAA\n",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var buf bytes.Buffer
			err := Subseq(strings.NewReader(subseqFASTA), strings.NewReader(tc.spec), &buf, tc.lineLen)
			if err != nil {
				t.Fatalf("Subseq failed: %v", err)
			}
			got := buf.String()
			if got != tc.want {
				t.Errorf("output = %q, want %q", got, tc.want)
			}
			for _, n := range tc.wantNot {
				if strings.Contains(got, n) {
					t.Errorf("output %q should not contain %q", got, n)
				}
			}
		})
	}
}

func TestSubseqBED(t *testing.T) {
	tests := []struct {
		name string
		spec string
		want string
	}{
		{
			name: "single region: chrom:start+1-end header and correct slice",
			spec: "chr1\t1\t4\n",
			want: ">chr1:2-4\nCGT\n",
		},
		{
			name: "multiple regions on same chrom, in BED order",
			spec: "chr1\t0\t2\nchr1\t4\t8\n",
			want: ">chr1:1-2\nAC\n>chr1:5-8\nACGT\n",
		},
		{
			name: "regions across multiple chroms",
			spec: "chr1\t0\t3\nchr2\t6\t10\n",
			want: ">chr1:1-3\nACG\n>chr2:7-10\nGGCC\n",
		},
		{
			name: "end is clamped to sequence length",
			spec: "chr3\t1\t100\n",
			want: ">chr3:2-5\nNNNN\n",
		},
		{
			name: "comment, track and browser lines are ignored",
			spec: "#comment\ntrack name=foo\nbrowser position chr1\nchr1\t2\t5\n",
			want: ">chr1:3-5\nGTA\n",
		},
		{
			name: "extra BED columns are ignored",
			spec: "chr2\t0\t4\tfeatureName\t100\t+\n",
			want: ">chr2:1-4\nTTTT\n",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var buf bytes.Buffer
			err := Subseq(strings.NewReader(subseqFASTA), strings.NewReader(tc.spec), &buf, 0)
			if err != nil {
				t.Fatalf("Subseq failed: %v", err)
			}
			if got := buf.String(); got != tc.want {
				t.Errorf("output = %q, want %q", got, tc.want)
			}
		})
	}
}

func TestSubseqBEDOutOfRange(t *testing.T) {
	// start >= len(seq): the region is skipped entirely (warning to stderr).
	var buf bytes.Buffer
	err := Subseq(strings.NewReader(subseqFASTA), strings.NewReader("chr3\t10\t20\n"), &buf, 0)
	if err != nil {
		t.Fatalf("Subseq failed: %v", err)
	}
	if got := buf.String(); got != "" {
		t.Errorf("out-of-range region should produce no output, got %q", got)
	}
}

func TestSubseqFormatAutodetect(t *testing.T) {
	tests := []struct {
		name      string
		spec      string
		wantIsBED bool
	}{
		{name: "three int-ish fields => BED", spec: "chr1\t0\t10\n", wantIsBED: true},
		{name: "space-separated three fields => BED", spec: "chr1 0 10\n", wantIsBED: true},
		{name: "extra columns still BED", spec: "chr1\t0\t10\tname\n", wantIsBED: true},
		{name: "two fields => name list", spec: "chr1\t0\n", wantIsBED: false},
		{name: "non-integer second field => name list", spec: "chr1\tfoo\t10\n", wantIsBED: false},
		{name: "non-integer third field => name list", spec: "chr1\t0\tbar\n", wantIsBED: false},
		{name: "single column => name list", spec: "chr1\n", wantIsBED: false},
		{name: "comments skipped before detection", spec: "# header\nchr1\t0\t10\n", wantIsBED: true},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			spec, err := parseSubseqSpec(strings.NewReader(tc.spec))
			if err != nil {
				t.Fatalf("parseSubseqSpec failed: %v", err)
			}
			if spec.isBED != tc.wantIsBED {
				t.Errorf("isBED = %v, want %v", spec.isBED, tc.wantIsBED)
			}
		})
	}
}

func TestMergePE(t *testing.T) {
	r1 := "@a/1\nAAAA\n+\nIIII\n@b/1\nCCCC\n+\nIIII\n@c/1\nGGGG\n+\nIIII\n"
	r2 := "@a/2\nTTTT\n+\nIIII\n@b/2\nGGGG\n+\nIIII\n@c/2\nAAAA\n+\nIIII\n"
	want := "@a/1\nAAAA\n+\nIIII\n" +
		"@a/2\nTTTT\n+\nIIII\n" +
		"@b/1\nCCCC\n+\nIIII\n" +
		"@b/2\nGGGG\n+\nIIII\n" +
		"@c/1\nGGGG\n+\nIIII\n" +
		"@c/2\nAAAA\n+\nIIII\n"

	var buf bytes.Buffer
	if err := MergePE(strings.NewReader(r1), strings.NewReader(r2), &buf); err != nil {
		t.Fatalf("MergePE failed: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("MergePE output = %q, want %q", got, want)
	}

	// Mismatched counts: in1 shorter.
	short1 := "@a/1\nAAAA\n+\nIIII\n"
	long2 := "@a/2\nTTTT\n+\nIIII\n@b/2\nGGGG\n+\nIIII\n"
	buf.Reset()
	err := MergePE(strings.NewReader(short1), strings.NewReader(long2), &buf)
	if err == nil {
		t.Fatal("expected error for mismatched counts (in1 shorter), got nil")
	}
	if !strings.Contains(err.Error(), "in1 is shorter") {
		t.Errorf("error %q should mention 'in1 is shorter'", err.Error())
	}

	// Mismatched counts: in2 shorter.
	long1 := "@a/1\nAAAA\n+\nIIII\n@b/1\nCCCC\n+\nIIII\n"
	short2 := "@a/2\nTTTT\n+\nIIII\n"
	buf.Reset()
	err = MergePE(strings.NewReader(long1), strings.NewReader(short2), &buf)
	if err == nil {
		t.Fatal("expected error for mismatched counts (in2 shorter), got nil")
	}
	if !strings.Contains(err.Error(), "in2 is shorter") {
		t.Errorf("error %q should mention 'in2 is shorter'", err.Error())
	}
}

func TestMergePEFasta(t *testing.T) {
	r1 := ">a/1\nAAAA\n>b/1\nCCCC\n"
	r2 := ">a/2\nTTTT\n>b/2\nGGGG\n"

	var buf bytes.Buffer
	if err := MergePE(strings.NewReader(r1), strings.NewReader(r2), &buf); err != nil {
		t.Fatalf("MergePE (fasta) failed: %v", err)
	}
	got := buf.String()
	// Order: a/1, a/2, b/1, b/2.
	wantOrder := []string{">a/1", ">a/2", ">b/1", ">b/2"}
	prev := -1
	for _, h := range wantOrder {
		idx := strings.Index(got, h)
		if idx < 0 {
			t.Fatalf("output missing header %q; got %q", h, got)
		}
		if idx <= prev {
			t.Fatalf("headers out of order in output %q (expected %v)", got, wantOrder)
		}
		prev = idx
	}
	// Sanity-check the sequences are present too.
	for _, s := range []string{"AAAA", "TTTT", "CCCC", "GGGG"} {
		if !strings.Contains(got, s) {
			t.Errorf("output missing sequence %q; got %q", s, got)
		}
	}
}

func TestCutN_BasicSplit(t *testing.T) {
	in := ">chr1\nACGNNNTGCANNNNG\n"
	want := ">chr1:1-3\nACG\n>chr1:7-10\nTGCA\n>chr1:15-15\nG\n"

	var buf bytes.Buffer
	if err := CutN(strings.NewReader(in), &buf, CutNOptions{MinN: 3}); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("CutN output = %q, want %q", got, want)
	}
}

func TestCutN_MultipleRuns(t *testing.T) {
	// Leading NN (kept, < 3), then a real cut at NNN, then an internal short
	// run (NN, kept, < 3), then another real cut at NNNNNN.
	// Sequence (1-based):
	//  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
	//  N  N  A  C  G  T  N  N  N  C  N  N  T  T  N  N  N  N  N  N
	// Internal NN at 11-12 is < 3, so it is NOT cut. Trailing NNNNNN (15-20)
	// is cut, leaving nothing after it.
	in := ">s\nNNACGTNNNCNNTTNNNNNN\n"
	want := ">s:1-6\nNNACGT\n>s:10-14\nCNNTT\n"

	var buf bytes.Buffer
	if err := CutN(strings.NewReader(in), &buf, CutNOptions{MinN: 3}); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("CutN output = %q, want %q", got, want)
	}
}

func TestCutN_AllNs(t *testing.T) {
	in := ">empty\nNNNNNNNN\n"

	var buf bytes.Buffer
	if err := CutN(strings.NewReader(in), &buf, CutNOptions{MinN: 3}); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	if got := buf.String(); got != "" {
		t.Errorf("all-N record should produce no output, got %q", got)
	}
}

func TestCutN_NoCutsRecordUnchanged(t *testing.T) {
	// No N-run, and no -g flag => emit unchanged with original header (no :start-end suffix).
	in := ">solo\nACGTACGT\n"
	want := ">solo\nACGTACGT\n"

	var buf bytes.Buffer
	if err := CutN(strings.NewReader(in), &buf, CutNOptions{MinN: 3}); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("CutN output = %q, want %q", got, want)
	}

	// Two short Ns (< 3) also count as "no qualifying run".
	in2 := ">solo\nACGNNACGT\n"
	want2 := ">solo\nACGNNACGT\n"
	buf.Reset()
	if err := CutN(strings.NewReader(in2), &buf, CutNOptions{MinN: 3}); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	if got := buf.String(); got != want2 {
		t.Errorf("CutN output = %q, want %q", got, want2)
	}
}

func TestCutN_GapsToStderr(t *testing.T) {
	in := ">chr1\nACGNNNTGCANNNNG\n"
	wantOut := ">chr1:1-3\nACG\n>chr1:7-10\nTGCA\n>chr1:15-15\nG\n"
	// 0-based half-open BED: NNN at [3,6), NNNN at [10,14).
	wantGaps := "chr1\t3\t6\tN\nchr1\t10\t14\tN\n"

	var out, gaps bytes.Buffer
	opts := CutNOptions{MinN: 3, EmitGaps: true, GapsW: &gaps}
	if err := CutN(strings.NewReader(in), &out, opts); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	if got := out.String(); got != wantOut {
		t.Errorf("CutN stdout = %q, want %q", got, wantOut)
	}
	if got := gaps.String(); got != wantGaps {
		t.Errorf("CutN gaps = %q, want %q", got, wantGaps)
	}
}

func TestCutN_FastqInput(t *testing.T) {
	// FASTQ input is auto-detected by the leading '@'; output is always FASTA.
	in := "@r1\nACGNNNTGCA\n+\nIIIIIIIIII\n"
	want := ">r1:1-3\nACG\n>r1:7-10\nTGCA\n"

	var buf bytes.Buffer
	if err := CutN(strings.NewReader(in), &buf, CutNOptions{MinN: 3}); err != nil {
		t.Fatalf("CutN failed: %v", err)
	}
	got := buf.String()
	if got != want {
		t.Errorf("CutN(FASTQ) output = %q, want %q", got, want)
	}
	if strings.Contains(got, "@") || strings.Contains(got, "IIII") {
		t.Errorf("FASTQ header/quality leaked into FASTA output: %q", got)
	}
}

func TestSubseqFastqInputFastaOutput(t *testing.T) {
	fastqData := "@read1 some desc\nACGTACGTAA\n+\nIIIIIIIIII\n@read2\nTTTTGGGGCC\n+\n##########\n"

	// Name list.
	var buf bytes.Buffer
	if err := Subseq(strings.NewReader(fastqData), strings.NewReader("read2\n"), &buf, 0); err != nil {
		t.Fatalf("Subseq failed: %v", err)
	}
	if got, want := buf.String(), ">read2\nTTTTGGGGCC\n"; got != want {
		t.Errorf("name-list FASTQ->FASTA: got %q, want %q", got, want)
	}

	// BED region.
	buf.Reset()
	if err := Subseq(strings.NewReader(fastqData), strings.NewReader("read1\t0\t4\n"), &buf, 0); err != nil {
		t.Fatalf("Subseq failed: %v", err)
	}
	if got, want := buf.String(), ">read1:1-4\nACGT\n"; got != want {
		t.Errorf("BED FASTQ->FASTA: got %q, want %q", got, want)
	}
	if strings.Contains(buf.String(), "I") || strings.Contains(buf.String(), "@") {
		t.Errorf("FASTQ quality/header leaked into FASTA output: %q", buf.String())
	}
}
