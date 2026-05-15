package bcftools

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// convertVCFFixture builds a deterministic multi-sample VCF body for the
// convert tests. Pre-sorted by (CHROM, POS) so the streaming pass-through
// produces stable output.
func convertVCFFixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	10	.	A	T	.	PASS	DP=10	GT	0/1	0/0	1/1
chr1	20	.	G	C	.	PASS	DP=20	GT	0/0	0/1	0/0
chr2	5	.	C	A	.	PASS	DP=5	GT	1/1	0/0	0/1
`
}

func TestConvert_RoundTripVCFPassthrough(t *testing.T) {
	var out bytes.Buffer
	n, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Convert: %v", err)
	}
	if n != 3 {
		t.Errorf("want 3 records emitted, got %d", n)
	}
	for _, want := range []string{
		"##fileformat=VCFv4.2",
		"S1\tS2\tS3",
		"chr1\t10",
		"chr2\t5",
		"DP=20",
	} {
		if !strings.Contains(out.String(), want) {
			t.Errorf("output missing %q:\n%s", want, out.String())
		}
	}
}

func TestConvert_VCFToGzipAndBack(t *testing.T) {
	dir := t.TempDir()
	gzPath := filepath.Join(dir, "fixture.vcf.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := Convert(strings.NewReader(convertVCFFixture()), f, ConvertOptions{
		OutputFormat: OutputVCFGz,
	}); err != nil {
		t.Fatalf("Convert gz: %v", err)
	}
	f.Close()

	// Read the gzipped output back through Convert and check we recovered the same body.
	r, err := os.Open(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	gr, err := gzip.NewReader(r)
	if err != nil {
		t.Fatalf("gzip read: %v", err)
	}
	defer gr.Close()
	var back bytes.Buffer
	if _, err := Convert(gr, &back, ConvertOptions{OutputFormat: OutputVCF}); err != nil {
		t.Fatalf("Convert from gz: %v", err)
	}
	for _, want := range []string{"chr1\t10", "chr1\t20", "chr2\t5"} {
		if !strings.Contains(back.String(), want) {
			t.Errorf("round-trip missing %q:\n%s", want, back.String())
		}
	}
}

func TestConvert_VCFToBCFRoundTrip(t *testing.T) {
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "fixture.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := Convert(strings.NewReader(convertVCFFixture()), f, ConvertOptions{
		OutputFormat: OutputBCF,
	}); err != nil {
		t.Fatalf("Convert bcf: %v", err)
	}
	f.Close()

	var back bytes.Buffer
	n, err := ConvertFile(bcfPath, &back, ConvertOptions{OutputFormat: OutputVCF})
	if err != nil {
		t.Fatalf("ConvertFile from bcf: %v", err)
	}
	if n != 3 {
		t.Errorf("want 3, got %d", n)
	}
	for _, want := range []string{"chr1\t10", "chr2\t5", "S1\tS2\tS3"} {
		if !strings.Contains(back.String(), want) {
			t.Errorf("BCF round-trip missing %q:\n%s", want, back.String())
		}
	}
}

func TestConvert_SamplesFilter(t *testing.T) {
	var out bytes.Buffer
	_, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		Samples:      []string{"S2"},
	})
	if err != nil {
		t.Fatalf("Convert: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "\tS2\n") {
		t.Errorf("want only S2 in column header line; got:\n%s", body)
	}
	if strings.Contains(body, "\tS1\t") || strings.Contains(body, "\tS3\t") {
		t.Errorf("want S1/S3 stripped; got:\n%s", body)
	}
	// First data record: original GT for S2 was 0/0.
	if !strings.Contains(body, "chr1\t10\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/0\n") {
		t.Errorf("expected S2-only column for chr1:10; got:\n%s", body)
	}
}

func TestConvert_RegionsPostFilter(t *testing.T) {
	var out bytes.Buffer
	_, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		Regions:      []string{"chr1:1-15"},
	})
	if err != nil {
		t.Fatalf("Convert: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "chr1\t10") {
		t.Errorf("want chr1:10 record kept:\n%s", body)
	}
	if strings.Contains(body, "chr1\t20") {
		t.Errorf("did not expect chr1:20 in region chr1:1-15:\n%s", body)
	}
	if strings.Contains(body, "chr2\t5") {
		t.Errorf("did not expect chr2:5 in region chr1:1-15:\n%s", body)
	}
}

func TestConvert_TargetsPostFilter(t *testing.T) {
	var out bytes.Buffer
	_, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		Targets:      []string{"chr2"},
	})
	if err != nil {
		t.Fatalf("Convert: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "chr2\t5") {
		t.Errorf("want chr2 record:\n%s", body)
	}
	if strings.Contains(body, "chr1") && strings.Contains(body, "PASS") {
		// chr1 metadata lines are fine; the assertion is no data row.
		for _, line := range strings.Split(body, "\n") {
			if strings.HasPrefix(line, "chr1\t") {
				t.Errorf("did not expect chr1 data row %q under -t chr2", line)
			}
		}
	}
}

func TestConvert_IncludeExclude(t *testing.T) {
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		IncludeExpr:  "INFO/DP>5",
	}); err != nil {
		t.Fatalf("Convert -i: %v", err)
	}
	body := out.String()
	if strings.Contains(body, "chr2\t5") {
		t.Errorf("DP=5 record should be excluded by INFO/DP>5:\n%s", body)
	}
	if !strings.Contains(body, "chr1\t20") {
		t.Errorf("DP=20 record should be retained by INFO/DP>5:\n%s", body)
	}

	out.Reset()
	if _, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		ExcludeExpr:  "INFO/DP<10",
	}); err != nil {
		t.Fatalf("Convert -e: %v", err)
	}
	body = out.String()
	if strings.Contains(body, "chr2\t5") {
		t.Errorf("DP=5 record should be dropped by -e INFO/DP<10:\n%s", body)
	}
}

func TestConvert_ForceSamplesGate(t *testing.T) {
	// Without --force-samples, a missing requested sample is an error.
	var out bytes.Buffer
	_, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		Samples:      []string{"S1", "NOPE"},
	})
	if err == nil {
		t.Fatal("expected error for missing sample without --force-samples")
	}

	// With --force-samples, the missing name is silently dropped.
	out.Reset()
	if _, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		Samples:      []string{"S1", "NOPE"},
		ForceSamples: true,
	}); err != nil {
		t.Fatalf("Convert --force-samples: %v", err)
	}
	if !strings.Contains(out.String(), "\tS1\n") {
		t.Errorf("want lone S1 sample after force-samples:\n%s", out.String())
	}
}

func TestConvert_BadRegion(t *testing.T) {
	var out bytes.Buffer
	_, err := Convert(strings.NewReader(convertVCFFixture()), &out, ConvertOptions{
		OutputFormat: OutputVCF,
		Regions:      []string{"chr1:not-a-num"},
	})
	if err == nil {
		t.Error("expected error parsing bad region")
	}
}

func TestConvertFile_MissingInput(t *testing.T) {
	var out bytes.Buffer
	if _, err := ConvertFile("/no/such/file.vcf", &out, ConvertOptions{}); err == nil {
		t.Error("expected error for missing input file")
	}
}

func TestConvertFile_SamplesFromFile(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(convertVCFFixture()), 0644); err != nil {
		t.Fatal(err)
	}
	sf := filepath.Join(dir, "samples.txt")
	if err := os.WriteFile(sf, []byte("S2\nS3\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ConvertFile(in, &out, ConvertOptions{
		OutputFormat: OutputVCF,
		SamplesFile:  sf,
	}); err != nil {
		t.Fatalf("ConvertFile: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "\tS2\tS3\n") {
		t.Errorf("want S2\\tS3 sample columns:\n%s", body)
	}
	if strings.Contains(body, "\tS1\t") {
		t.Errorf("S1 should be filtered out:\n%s", body)
	}
}

// Ensure the read path tolerates pipes (no seeking).
func TestConvert_StreamingReader(t *testing.T) {
	r, w := io.Pipe()
	go func() {
		_, _ = io.WriteString(w, convertVCFFixture())
		_ = w.Close()
	}()
	var out bytes.Buffer
	if _, err := Convert(r, &out, ConvertOptions{OutputFormat: OutputVCF}); err != nil {
		t.Fatalf("Convert from pipe: %v", err)
	}
	if !strings.Contains(out.String(), "chr1\t10") {
		t.Errorf("want first record:\n%s", out.String())
	}
}
