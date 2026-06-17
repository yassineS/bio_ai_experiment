package htsfile

import (
	"bytes"
	"compress/gzip"
	"os"
	"strings"
	"testing"
)

// TestClassifyPayload covers the per-format payload heuristics
// directly (no compression wrapping). Each case is a minimal
// synthetic prefix sufficient to trigger the matcher.
func TestClassifyPayload(t *testing.T) {
	tests := []struct {
		name    string
		prefix  []byte
		want    Payload
		version string
	}{
		{"bam-magic", []byte("BAM\x01... rest"), PayloadBAM, "1"},
		{"bcf-2.2", []byte("BCF\x02\x02 ... rest"), PayloadBCF, "2.2"},
		{"bcf-2.1", []byte("BCF\x02\x01 ... rest"), PayloadBCF, "2.1"},
		{"cram-3.0", []byte("CRAM\x03\x00 ..."), PayloadCRAM, "3.0"},
		{"vcf-header", []byte("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"), PayloadVCF, "4.2"},
		{"sam-hd", []byte("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\n"), PayloadSAM, "1.6"},
		{"sam-sq-only", []byte("@SQ\tSN:chr1\tLN:1000\n"), PayloadSAM, ""},
		{"fasta", []byte(">chr1\nACGTACGT\n"), PayloadFASTA, ""},
		{"fastq", []byte("@read1\nACGT\n+\n!!!!\n"), PayloadFASTQ, ""},
		{"gff", []byte("##gff-version 3\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n"), PayloadGFF, "3"},
		{"bed", []byte("chr1\t100\t200\tfeat1\n"), PayloadBED, ""},
		{"bed-with-track-line", []byte("track name=foo\nchr1\t100\t200\n"), PayloadBED, ""},
		{"plain-text", []byte("hello world\nthis is just text.\n"), PayloadText, ""},
		{"empty", []byte{}, PayloadUnknown, ""},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			f := classifyPayload(tt.prefix)
			if f.Payload != tt.want {
				t.Errorf("payload: got %s want %s", f.Payload, tt.want)
			}
			if f.Version != tt.version {
				t.Errorf("version: got %q want %q", f.Version, tt.version)
			}
		})
	}
}

// TestDetectCompression verifies the magic-byte sniff for plain
// vs gzip vs BGZF.
func TestDetectCompression(t *testing.T) {
	plain := []byte("hello world this is plain text\n")
	if got := detectCompression(plain); got != CompressionPlain {
		t.Errorf("plain: got %s", got)
	}

	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	_, _ = zw.Write([]byte("payload"))
	_ = zw.Close()
	if got := detectCompression(buf.Bytes()); got != CompressionGzip {
		t.Errorf("gzip: got %s", got)
	}

	// BGZF: gzip header with FEXTRA set and a BC subfield. We
	// hand-craft a 18-byte prefix matching what bgzf.NewWriter would
	// emit.
	bgzfHeader := []byte{
		0x1f, 0x8b, // magic
		0x08,       // CM = deflate
		0x04,       // FLG = FEXTRA
		0, 0, 0, 0, // MTIME
		0,          // XFL
		0xff,       // OS
		0x06, 0x00, // XLEN = 6 (extra subfields total length, LE)
		'B', 'C', // SI1 SI2
		0x02, 0x00, // SLEN = 2 (LE)
		0x10, 0x00, // BSIZE-1 = 16 (LE)
	}
	if got := detectCompression(bgzfHeader); got != CompressionBGZF {
		t.Errorf("bgzf: got %s", got)
	}
}

// TestIdentifyReader_GzippedVCF wraps a minimal VCF in gzip and
// verifies the sniffer correctly reports VCF + gzip-compressed.
func TestIdentifyReader_GzippedVCF(t *testing.T) {
	const body = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write([]byte(body)); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	f, err := IdentifyReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatal(err)
	}
	if f.Payload != PayloadVCF || f.Version != "4.3" {
		t.Errorf("payload/version: got %s %s, want VCF 4.3", f.Payload, f.Version)
	}
	if f.Compression != CompressionGzip {
		t.Errorf("compression: got %s want gzip", f.Compression)
	}
}

// TestDescribe verifies the one-line htsfile-style summary string for
// a couple of common cases.
func TestDescribe(t *testing.T) {
	cases := []struct {
		f    Format
		want string
	}{
		{
			Format{Compression: CompressionBGZF, Payload: PayloadBAM, Version: "1"},
			"BAM version 1 compressed sequence data",
		},
		{
			Format{Compression: CompressionGzip, Payload: PayloadVCF, Version: "4.2"},
			"VCF version 4.2 gzip-compressed variant calling data",
		},
		{
			Format{Compression: CompressionPlain, Payload: PayloadFASTA},
			"FASTA sequence text",
		},
		{
			Format{Compression: CompressionPlain, Payload: PayloadBED},
			"BED genomic region text",
		},
	}
	for _, tt := range cases {
		t.Run(tt.want, func(t *testing.T) {
			if got := tt.f.Describe(); got != tt.want {
				t.Errorf("got %q want %q", got, tt.want)
			}
		})
	}
}

// TestIdentifyPlain_FASTQDoesNotMisidentifyAsSAM ensures a FASTQ
// record beginning with '@' is correctly routed to FASTQ rather than
// SAM (the SAM detector matches @HD/@SQ/@RG/@PG/@CO specifically).
func TestIdentifyPlain_FASTQDoesNotMisidentifyAsSAM(t *testing.T) {
	const body = "@read1\nACGT\n+\n!!!!\n"
	f, err := IdentifyReader(strings.NewReader(body))
	if err != nil {
		t.Fatal(err)
	}
	if f.Payload != PayloadFASTQ {
		t.Errorf("payload: got %s want FASTQ", f.Payload)
	}
}

// TestIdentify_PathHappyPath exercises Identify(path) end-to-end with
// a temporary file. The IdentifyReader path is covered separately by
// the other tests; this one specifically pins the file-open +
// identify path used by the CLI.
func TestIdentify_PathHappyPath(t *testing.T) {
	dir := t.TempDir()
	path := dir + "/sample.vcf"
	body := []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	if err := os.WriteFile(path, body, 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := Identify(path)
	if err != nil {
		t.Fatal(err)
	}
	if f.Payload != PayloadVCF || f.Version != "4.2" {
		t.Errorf("got %s %s, want VCF 4.2", f.Payload, f.Version)
	}
}

// TestIdentify_MissingFile pins the error path when the file at path
// doesn't exist. The CLI surfaces the error to stderr and exits 1; the
// library API just returns a wrapped error.
func TestIdentify_MissingFile(t *testing.T) {
	if _, err := Identify("/nope/does-not-exist.bam"); err == nil {
		t.Fatal("expected error for missing file")
	}
}

// TestLooksLikeBED_RejectsNonIntegerCoordinates verifies that a
// 3-column file whose 2nd/3rd columns aren't non-negative integers
// does NOT trigger the BED detector (regression-test the heuristic's
// strictness).
func TestLooksLikeBED_RejectsNonIntegerCoordinates(t *testing.T) {
	prefix := []byte("chr1\tabc\t200\tfeat1\n")
	f := classifyPayload(prefix)
	if f.Payload == PayloadBED {
		t.Errorf("non-integer coords should not classify as BED; got %s", f.Payload)
	}
}
