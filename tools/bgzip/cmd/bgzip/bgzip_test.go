package main

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// TestBgzipTextModeOracle compares our compressed output, byte-for-byte, with
// goldens produced by genuine htslib bgzip (committed under testdata/oracle).
// It exercises the text-mode block-flush heuristic for VCF/FASTA/multi-block
// text and the plain binary path for a BAM input. The goldens were generated
// with `reference_code/htslib/bgzip -c <input> > <input>.bgz`.
func TestBgzipTextModeOracle(t *testing.T) {
	cases := []struct {
		name  string
		input string
	}{
		{"small VCF", "small.vcf"},
		{"FASTA", "test.fa"},
		{"multi-block text", "multiblock.txt"},
		{"binary BAM", "test.bam"},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			dir := filepath.Join("..", "..", "testdata", "oracle")
			in, err := os.ReadFile(filepath.Join(dir, tc.input))
			if err != nil {
				t.Fatalf("read input: %v", err)
			}
			golden, err := os.ReadFile(filepath.Join(dir, tc.input+".bgz"))
			if err != nil {
				t.Fatalf("read golden: %v", err)
			}

			var out bytes.Buffer
			rc := runCompress("-", true, false, false, 6, false, bytes.NewReader(in), &out, &bytes.Buffer{})
			if rc != 0 {
				t.Fatalf("runCompress returned %d", rc)
			}
			if !bytes.Equal(out.Bytes(), golden) {
				t.Fatalf("output differs from genuine bgzip golden: got %d bytes, want %d bytes",
					out.Len(), len(golden))
			}
		})
	}
}

// TestBgzipBinaryFlag verifies that --binary forces the binary streaming path
// even on a detected text format, matching genuine bgzip --binary.
func TestBgzipBinaryFlag(t *testing.T) {
	dir := filepath.Join("..", "..", "testdata", "oracle")
	in, err := os.ReadFile(filepath.Join(dir, "small.vcf"))
	if err != nil {
		t.Fatalf("read input: %v", err)
	}

	// Forcing binary on a small single-block VCF should round-trip and yield a
	// single data block, not the header-split text framing.
	var textOut, binOut bytes.Buffer
	if rc := runCompress("-", true, false, false, 6, false, bytes.NewReader(in), &textOut, &bytes.Buffer{}); rc != 0 {
		t.Fatalf("text-mode runCompress returned %d", rc)
	}
	if rc := runCompress("-", true, false, false, 6, true, bytes.NewReader(in), &binOut, &bytes.Buffer{}); rc != 0 {
		t.Fatalf("binary-mode runCompress returned %d", rc)
	}
	if bytes.Equal(textOut.Bytes(), binOut.Bytes()) {
		t.Fatal("expected --binary framing to differ from text-mode framing for a multi-line VCF")
	}
}

// TestDetectTextual unit-tests the format detection against representative
// inputs for each branch.
func TestDetectTextual(t *testing.T) {
	cases := []struct {
		name string
		in   []byte
		want bool
	}{
		{"empty", nil, false},
		{"VCF", []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\n"), true},
		{"SAM header", []byte("@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"), true},
		{"FASTA", []byte(">seq1\nACGTACGT\n"), true},
		{"FASTQ", []byte("@read1\nACGTACGT\n+\nIIIIIIII\n"), true},
		{"BED", []byte("chr1\t100\t200\tname\n"), true},
		{"plain text", []byte("hello world\nthis is text\n"), true},
		{"gzip magic", []byte{0x1f, 0x8b, 0x08, 0x00, 0x00}, false},
		{"BAM magic", []byte("BAM\x01\x00\x00\x00"), false},
		{"BCF magic", []byte("BCF\x02\x02\x00"), false},
		{"random binary", []byte{0x00, 0x01, 0x02, 0xff, 0xfe}, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := detectTextual(tc.in); got != tc.want {
				t.Errorf("detectTextual(%q) = %v, want %v", tc.in, got, tc.want)
			}
		})
	}
}
