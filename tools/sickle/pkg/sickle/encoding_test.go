package sickle

import (
	"bufio"
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// makeFastq builds a minimal FASTQ string with one record where seq and qual
// have the same length (qual is repeated to fill seq).
func makeFastq(id, seq, qualChar string) string {
	q := strings.Repeat(qualChar, len(seq))
	return "@" + id + "\n" + seq + "\n+\n" + q + "\n"
}

func TestDetectEncoding(t *testing.T) {
	tests := []struct {
		name      string
		input     string
		wantName  string
		wantEnc   fastq.QualityEncoding
		wantErr   bool
		ambiguous bool
	}{
		{
			name: "clear sanger Phred+33",
			// '!'..'I' = ASCII 33..73 → Phred+33 quality 0..40
			input:    makeFastq("r1", "ACGTACGT", "!\"#$%&'(") + makeFastq("r2", "ACGTACGT", "IIIIIIII"),
			wantName: "sanger",
			wantEnc:  fastq.Phred33,
		},
		{
			name: "clear illumina Phred+64",
			// '@'..'h' = ASCII 64..104 → Phred+64 quality 0..40
			input:    makeFastq("r1", "ACGTACGT", "ABCDEFGH") + makeFastq("r2", "ACGTACGT", "hhhhhhhh"),
			wantName: "illumina",
			wantEnc:  fastq.Phred64,
		},
		{
			name: "all ! is sanger min boundary",
			// min=33, max=33 → matches min<64 && max<=73 branch.
			input:    makeFastq("r1", "ACGT", "!!!!"),
			wantName: "sanger",
			wantEnc:  fastq.Phred33,
		},
		{
			name: "all h is illumina max boundary",
			// min=104, max=104 → matches min>=64 && max<=104 branch.
			input:    makeFastq("r1", "ACGT", "hhhh"),
			wantName: "illumina",
			wantEnc:  fastq.Phred64,
		},
		{
			name: "solexa min in [59,64)",
			// ';' = 59, 'h' = 104 — min in solexa range, max within illumina max.
			input:    makeFastq("r1", "ACGTACGT", ";<=>?@AB") + makeFastq("r2", "ACGTACGT", "hhhhhhhh"),
			wantName: "solexa",
			wantEnc:  fastq.Phred64,
		},
		{
			name: "Illumina 1.8+ where Phred+33 reaches above ASCII 73",
			// 'J' = 74, '#' = 35 → min<64 && max>73, still Phred+33.
			input:    makeFastq("r1", "ACGTACGT", "##JJJJJJ"),
			wantName: "sanger",
			wantEnc:  fastq.Phred33,
		},
		{
			name:    "invalid: quality byte below 33",
			input:   "@r1\nACGT\n+\n\x10\x10\x10\x10\n",
			wantErr: true,
		},
		{
			name:    "empty input",
			input:   "",
			wantErr: true,
		},
		{
			name:    "malformed header",
			input:   "not-a-fastq\nACGT\n+\nIIII\n",
			wantErr: true,
		},
		{
			name: "ambiguous: high min, very high max → fallback sanger",
			// Force the default branch: min>=64 but max>104 (e.g. 'z'=122).
			input:     makeFastq("r1", "ACGT", "zzzz"),
			wantName:  "sanger",
			wantEnc:   fastq.Phred33,
			ambiguous: true,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			br := bufio.NewReaderSize(strings.NewReader(tc.input), 64*1024)
			res, err := DetectEncoding(br)
			if tc.wantErr {
				if err == nil {
					t.Fatalf("expected error, got result %+v", res)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if res.Name != tc.wantName {
				t.Errorf("Name: got %q, want %q (min=%d max=%d sampled=%d)", res.Name, tc.wantName, res.MinByte, res.MaxByte, res.Sampled)
			}
			if res.Encoding != tc.wantEnc {
				t.Errorf("Encoding: got %v, want %v", res.Encoding, tc.wantEnc)
			}
			if res.Ambiguous != tc.ambiguous {
				t.Errorf("Ambiguous: got %v, want %v", res.Ambiguous, tc.ambiguous)
			}
			// The peek must NOT consume the underlying reader: we should still
			// be able to read the first byte ('@') from br.
			b, perr := br.Peek(1)
			if perr != nil || len(b) == 0 || b[0] != '@' {
				t.Errorf("expected first byte still readable, got %q err=%v", b, perr)
			}
		})
	}
}

func TestEncodingFromName(t *testing.T) {
	tests := []struct {
		name    string
		want    fastq.QualityEncoding
		wantErr bool
	}{
		{"sanger", fastq.Phred33, false},
		{"phred33", fastq.Phred33, false},
		{"illumina", fastq.Phred64, false},
		{"phred64", fastq.Phred64, false},
		{"solexa", fastq.Phred64, false},
		{"unknown", fastq.Phred33, true},
		{"", fastq.Phred33, true},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, err := EncodingFromName(tc.name)
			if (err != nil) != tc.wantErr {
				t.Fatalf("err=%v wantErr=%v", err, tc.wantErr)
			}
			if got != tc.want {
				t.Errorf("got %v, want %v", got, tc.want)
			}
		})
	}
}

func TestTrimSingleEndAutoDetect(t *testing.T) {
	// Build a sanger FASTQ where each read is mostly high-qual ('I'=40, threshold 30)
	// then a low-qual tail ('#'=2). Auto-detected as sanger should trim the tail.
	sangerFastq := makeFastq("read1", "ACGTACGTACGTACGTACGTACGT", "I") +
		makeFastq("read2", "ACGTACGTACGTACGTACGTACGT", "I")
	// Make the input have a clearly-Phred+33 fingerprint by including a low byte.
	sangerFastq = "@hint\nACGT\n+\n!!!!\n" + sangerFastq

	br := bufio.NewReaderSize(strings.NewReader(sangerFastq), 64*1024)
	res, err := DetectEncoding(br)
	if err != nil {
		t.Fatalf("DetectEncoding(sanger) failed: %v", err)
	}
	if res.Name != "sanger" || res.Encoding != fastq.Phred33 {
		t.Fatalf("expected sanger/Phred33, got %s/%v", res.Name, res.Encoding)
	}

	// Continue trimming from the same buffered reader — peeked bytes should
	// still be available so all 3 records are processed.
	var out bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 30
	opts.LengthThreshold = 4
	stats, err := TrimSingleEnd(br, &out, res.Encoding, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}
	if stats.TotalReads != 3 {
		t.Errorf("expected 3 total reads, got %d", stats.TotalReads)
	}
	if !strings.Contains(out.String(), "read1") || !strings.Contains(out.String(), "read2") {
		t.Errorf("expected both high-qual reads in output, got: %s", out.String())
	}

	// Now an Illumina/Phred+64 input. 'h'=104 → quality 40; '!' isn't valid
	// here, so use '@'=64 for low quality.
	illuminaFastq := makeFastq("read1", "ACGTACGTACGTACGTACGTACGT", "h") +
		makeFastq("read2", "ACGTACGTACGTACGTACGTACGT", "h")
	br2 := bufio.NewReaderSize(strings.NewReader(illuminaFastq), 64*1024)
	res2, err := DetectEncoding(br2)
	if err != nil {
		t.Fatalf("DetectEncoding(illumina) failed: %v", err)
	}
	if res2.Name != "illumina" || res2.Encoding != fastq.Phred64 {
		t.Fatalf("expected illumina/Phred64, got %s/%v", res2.Name, res2.Encoding)
	}

	var out2 bytes.Buffer
	stats2, err := TrimSingleEnd(br2, &out2, res2.Encoding, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd(illumina) failed: %v", err)
	}
	if stats2.TotalReads != 2 {
		t.Errorf("expected 2 total reads, got %d", stats2.TotalReads)
	}
	if !strings.Contains(out2.String(), "read1") {
		t.Errorf("expected read1 in illumina output, got: %s", out2.String())
	}
}

func TestPairedAutoDetectFromR1(t *testing.T) {
	// Sanger R1 + sanger R2, same record count. Detection runs ONCE on R1 and
	// the result is applied to both files — exactly what the CLI does. We
	// give R1 a low-byte "hint" record (qual '!') so the encoding range is
	// unambiguously Phred+33; R2 has the same encoding but isn't sampled.
	r1 := "@det/1\nACGT\n+\n!!!!\n" +
		makeFastq("p1/1", "ACGTACGTACGTACGTACGT", "I") +
		makeFastq("p2/1", "ACGTACGTACGTACGTACGT", "I")
	r2 := "@det/2\nTGCA\n+\nIIII\n" +
		makeFastq("p1/2", "TGCATGCATGCATGCATGCA", "I") +
		makeFastq("p2/2", "TGCATGCATGCATGCATGCA", "I")

	br1 := bufio.NewReaderSize(strings.NewReader(r1), 64*1024)
	br2 := bufio.NewReaderSize(strings.NewReader(r2), 64*1024)

	res, err := DetectEncoding(br1)
	if err != nil {
		t.Fatalf("DetectEncoding(R1) failed: %v", err)
	}
	if res.Name != "sanger" || res.Encoding != fastq.Phred33 {
		t.Fatalf("expected sanger/Phred33 from R1, got %s/%v", res.Name, res.Encoding)
	}

	// Apply detected encoding to BOTH br1 and br2 (without re-detecting on R2).
	// br1's peeked bytes must still be available for the trimmer to read.
	var out1, out2, outS bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 30
	opts.LengthThreshold = 4
	stats, err := TrimPairedEnd(br1, br2, &out1, &out2, &outS, res.Encoding, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd failed: %v", err)
	}
	if stats.TotalReads != 6 {
		t.Errorf("expected 6 total reads (3 pairs), got %d", stats.TotalReads)
	}
	if !strings.Contains(out1.String(), "p1/1") || !strings.Contains(out2.String(), "p1/2") {
		t.Errorf("expected paired output to contain both mates of p1; got R1=%q R2=%q", out1.String(), out2.String())
	}
	if !strings.Contains(out1.String(), "det/1") || !strings.Contains(out2.String(), "det/2") {
		t.Errorf("expected the leading detection-hint pair to also be processed; got R1=%q R2=%q", out1.String(), out2.String())
	}
}
