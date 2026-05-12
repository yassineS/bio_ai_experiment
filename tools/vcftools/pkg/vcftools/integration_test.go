package vcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Integration tests for vcftools

const testVCF = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	G	30	PASS	DP=60	GT:DP	0/0:20	0/1:20	1/1:20
chr1	200	.	C	T	25	PASS	DP=45	GT:DP	0/0:15	0/0:15	0/1:15
chr1	300	.	G	A	40	PASS	DP=90	GT:DP	0/1:30	1/1:30	1/1:30
chr2	100	.	T	C	35	PASS	DP=75	GT:DP	0/1:25	0/1:25	1/1:25
chr2	200	.	A	AT	28	PASS	DP=60	GT:DP	0/0:20	0/0:20	0/1:20
`

func TestIntegration_BasicFiltering(t *testing.T) {
	reader := strings.NewReader(testVCF)

	params := &Params{
		OutPrefix: filepath.Join(t.TempDir(), "test"),
		Chr:       "chr1",
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}
}

func TestIntegration_FrequencyCalculation(t *testing.T) {
	reader := strings.NewReader(testVCF)
	tmpDir := t.TempDir()

	params := &Params{
		OutPrefix: filepath.Join(tmpDir, "test"),
		Freq:      true,
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// Check output file exists
	freqFile := filepath.Join(tmpDir, "test.frq")
	if _, err := os.Stat(freqFile); os.IsNotExist(err) {
		t.Errorf("Expected frequency file %s to exist", freqFile)
	}

	// Read and verify content
	content, err := os.ReadFile(freqFile)
	if err != nil {
		t.Fatalf("Failed to read frequency file: %v", err)
	}

	if !strings.Contains(string(content), "CHROM") {
		t.Error("Frequency file should contain header")
	}
}

func TestIntegration_RecodeVCF(t *testing.T) {
	reader := strings.NewReader(testVCF)
	tmpDir := t.TempDir()

	params := &Params{
		OutPrefix:     filepath.Join(tmpDir, "test"),
		Recode:        true,
		RecodeInfoAll: true,
		RemoveIndels:  true,
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// Check output file exists
	recodeFile := filepath.Join(tmpDir, "test.recode.vcf")
	if _, err := os.Stat(recodeFile); os.IsNotExist(err) {
		t.Errorf("Expected recoded VCF file %s to exist", recodeFile)
	}

	// Read and verify content
	content, err := os.ReadFile(recodeFile)
	if err != nil {
		t.Fatalf("Failed to read recoded file: %v", err)
	}

	vcfStr := string(content)

	// Should have header
	if !strings.Contains(vcfStr, "##fileformat") {
		t.Error("Recoded VCF should contain fileformat header")
	}

	// Should not have indels (chr2:200 is A->AT)
	if strings.Contains(vcfStr, "chr2\t200") {
		t.Error("Recoded VCF should not contain indels when using --remove-indels")
	}

	// Should have SNPs
	if !strings.Contains(vcfStr, "chr1\t100") {
		t.Error("Recoded VCF should contain SNPs")
	}
}

func TestIntegration_MultipleStatistics(t *testing.T) {
	reader := strings.NewReader(testVCF)
	tmpDir := t.TempDir()

	params := &Params{
		OutPrefix:     filepath.Join(tmpDir, "test"),
		Freq:          true,
		SiteMeanDepth: true,
		MissingSite:   true,
		Hardy:         true,
		TsTvSummary:   true,
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// Check all output files exist
	expectedFiles := []string{
		"test.frq",
		"test.ldepth.mean",
		"test.lmiss",
		"test.hwe",
		"test.TsTv.summary",
	}

	for _, filename := range expectedFiles {
		path := filepath.Join(tmpDir, filename)
		if _, err := os.Stat(path); os.IsNotExist(err) {
			t.Errorf("Expected file %s to exist", filename)
		}
	}
}

func TestIntegration_SampleFiltering(t *testing.T) {
	reader := strings.NewReader(testVCF)
	tmpDir := t.TempDir()

	params := &Params{
		OutPrefix:     filepath.Join(tmpDir, "test"),
		IndvList:      []string{"sample1", "sample3"},
		Recode:        true,
		RecodeInfoAll: true,
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// Check output file exists and has correct samples
	recodeFile := filepath.Join(tmpDir, "test.recode.vcf")
	content, err := os.ReadFile(recodeFile)
	if err != nil {
		t.Fatalf("Failed to read recoded file: %v", err)
	}

	vcfStr := string(content)

	// Should have sample1 and sample3
	if !strings.Contains(vcfStr, "sample1") || !strings.Contains(vcfStr, "sample3") {
		t.Error("Recoded VCF should contain sample1 and sample3")
	}

	// Should not have sample2
	lines := strings.Split(vcfStr, "\n")
	for _, line := range lines {
		if strings.HasPrefix(line, "#CHROM") {
			if strings.Contains(line, "sample2") {
				t.Error("Recoded VCF should not contain sample2")
			}
		}
	}
}

func TestIntegration_AlleleFrequencyFilter(t *testing.T) {
	reader := strings.NewReader(testVCF)

	var buf bytes.Buffer
	params := &Params{
		OutPrefix: "test",
		UseStdout: true,
		Maf:       0.4, // Should filter out low MAF sites
		Recode:    true,
	}

	// Since we're using stdout, we need to capture it differently
	// For now, just verify it doesn't error
	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	_ = buf // Placeholder for potential stdout capture
}

func TestIntegration_QualityFilter(t *testing.T) {
	reader := strings.NewReader(testVCF)
	tmpDir := t.TempDir()

	params := &Params{
		OutPrefix: filepath.Join(tmpDir, "test"),
		MinQ:      30, // Should filter sites with QUAL < 30
		Freq:      true,
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// Check frequency file
	freqFile := filepath.Join(tmpDir, "test.frq")
	content, err := os.ReadFile(freqFile)
	if err != nil {
		t.Fatalf("Failed to read frequency file: %v", err)
	}

	// Should only have sites with QUAL >= 30
	// chr1:100 (QUAL=30), chr1:300 (QUAL=40), chr2:100 (QUAL=35)
	// Should NOT have chr1:200 (QUAL=25), chr2:200 (QUAL=28)
	lines := strings.Split(string(content), "\n")
	dataLines := 0
	for _, line := range lines {
		if strings.HasPrefix(line, "chr") {
			dataLines++
		}
	}

	// Should have 3 sites (chr1:100, chr1:300, chr2:100)
	if dataLines != 3 {
		t.Errorf("Expected 3 sites with QUAL >= 30, got %d", dataLines)
	}
}

func TestIntegration_PositionRange(t *testing.T) {
	reader := strings.NewReader(testVCF)
	tmpDir := t.TempDir()

	params := &Params{
		OutPrefix: filepath.Join(tmpDir, "test"),
		Chr:       "chr1",
		FromBp:    150,
		ToBp:      250,
		Freq:      true,
	}

	err := Run(reader, params)
	if err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// Check frequency file - should only have chr1:200
	freqFile := filepath.Join(tmpDir, "test.frq")
	content, err := os.ReadFile(freqFile)
	if err != nil {
		t.Fatalf("Failed to read frequency file: %v", err)
	}

	lines := strings.Split(string(content), "\n")
	dataLines := 0
	for _, line := range lines {
		if strings.HasPrefix(line, "chr") {
			dataLines++
			if !strings.Contains(line, "chr1\t200") {
				t.Errorf("Unexpected line in filtered output: %s", line)
			}
		}
	}

	if dataLines != 1 {
		t.Errorf("Expected 1 site in range [150, 250], got %d", dataLines)
	}
}

func TestIntegration_NewStatistics(t *testing.T) {
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "test")

	params := &Params{
		OutPrefix:    prefix,
		SitePi:       true,
		WindowPi:     1000,
		WindowPiStep: 500,
		TsTvByCount:  true,
		Depth:        true,
	}
	if err := Run(strings.NewReader(testVCF), params); err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	for _, suffix := range []string{".sites.pi", ".windowed.pi", ".TsTv.count", ".idepth"} {
		b, err := os.ReadFile(prefix + suffix)
		if err != nil {
			t.Fatalf("expected output file %s: %v", prefix+suffix, err)
		}
		if !strings.Contains(string(b), "\t") {
			t.Errorf("%s looks malformed:\n%s", prefix+suffix, string(b))
		}
	}

	// Each of the 3 samples appears at 5 sites with DP 20,15,30,25,20 -> mean 22.
	idepth, err := os.ReadFile(prefix + ".idepth")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(idepth), "22.00000") {
		t.Errorf(".idepth should report mean depth 22, got:\n%s", string(idepth))
	}

	// chr1 site 100 is balanced biallelic over 6 chromosomes: pi = 0.6.
	sitesPi, err := os.ReadFile(prefix + ".sites.pi")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(sitesPi), "chr1\t100\t0.600000") {
		t.Errorf(".sites.pi should report 0.6 for chr1:100, got:\n%s", string(sitesPi))
	}
}
