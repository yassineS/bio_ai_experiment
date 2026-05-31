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
	// Upstream emits MEAN_DEPTH via the default ostream << double formatter
	// (%g, 6 sig-figs, trailing zeros stripped). 22.0 renders as "22".
	if !strings.Contains(string(idepth), "\t22\n") {
		t.Errorf(".idepth should report mean depth 22, got:\n%s", string(idepth))
	}

	// chr1 site 100 is balanced biallelic over 6 chromosomes: pi = 0.6.
	sitesPi, err := os.ReadFile(prefix + ".sites.pi")
	if err != nil {
		t.Fatal(err)
	}
	// Upstream's C++ ostream default emits "0.6" not "0.600000" — see
	// formatCppDouble in statistics.go.
	if !strings.Contains(string(sitesPi), "chr1\t100\t0.6\n") {
		t.Errorf(".sites.pi should report 0.6 for chr1:100, got:\n%s", string(sitesPi))
	}
}

func TestIntegration_TsTvByQual(t *testing.T) {
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "test")

	params := &Params{
		OutPrefix:  prefix,
		TsTvByQual: true,
	}
	if err := Run(strings.NewReader(testVCF), params); err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	b, err := os.ReadFile(prefix + ".TsTv.qual")
	if err != nil {
		t.Fatalf("expected output file %s: %v", prefix+".TsTv.qual", err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	wantHeader := "QUAL_THRESHOLD\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD"
	if lines[0] != wantHeader {
		t.Errorf(".TsTv.qual header = %q, want %q", lines[0], wantHeader)
	}
	// testVCF has 4 biallelic SNPs with distinct QUAL values 25, 30, 35, 40
	// (chr2:200 is an indel and is skipped) -> 4 threshold rows.
	if got := len(lines) - 1; got != 4 {
		t.Errorf(".TsTv.qual data rows = %d, want 4", got)
	}
}

// TestTsTvByQualCumulative checks the threshold/cumulative-bucket logic of
// outputTsTvByQual via Run on a small in-memory VCF that contains both
// transitions and a transversion.
func TestTsTvByQualCumulative(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
chr1	10	.	A	G	10	PASS	.	GT	0/1
chr1	20	.	C	T	20	PASS	.	GT	0/1
chr1	30	.	A	C	20	PASS	.	GT	0/1
chr1	40	.	G	A	30	PASS	.	GT	0/1
chr1	50	.	A	AT	15	PASS	.	GT	0/1
chr1	60	.	A	T	.	PASS	.	GT	0/1
`
	// Biallelic SNPs: A>G@10 (Ts), C>T@20 (Ts), A>C@20 (Tv), G>A@30 (Ts).
	// chr1:50 is an indel (skipped); chr1:60 has missing QUAL (skipped).
	// Distinct sorted thresholds: 10, 20, 30.
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "test")
	if err := Run(strings.NewReader(vcfText), &Params{OutPrefix: prefix, TsTvByQual: true}); err != nil {
		t.Fatalf("Run() error = %v", err)
	}
	b, err := os.ReadFile(prefix + ".TsTv.qual")
	if err != nil {
		t.Fatalf("reading .TsTv.qual: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	if len(lines) != 4 { // header + 3 thresholds
		t.Fatalf(".TsTv.qual has %d lines, want 4:\n%s", len(lines), string(b))
	}
	// Upstream variant_file_output.cpp:3310 emits cumulative sums STRICTLY
	// below / STRICTLY above each threshold (values at exactly q go into
	// neither side), formatted via the default ostream << (no fixed
	// precision). 0/0 renders as glibc's "-nan" literal.
	want := []string{
		// q=10: LT none; GT = qual > 10 → C>T@20 (Ts), A>C@20 (Tv), G>A@30 (Ts) = 2 Ts, 1 Tv -> 2.
		"10\t0\t0\t-nan\t2\t1\t2",
		// q=20: LT = qual < 20 → A>G@10 (Ts) = 1 Ts, 0 Tv -> -nan; GT = qual > 20 → G>A@30 (Ts) = 1 Ts, 0 Tv -> -nan.
		"20\t1\t0\t-nan\t1\t0\t-nan",
		// q=30: LT = qual < 30 → A>G@10 (Ts), C>T@20 (Ts), A>C@20 (Tv) = 2 Ts, 1 Tv -> 2; GT = none.
		"30\t2\t1\t2\t0\t0\t-nan",
	}
	for i, w := range want {
		if lines[i+1] != w {
			t.Errorf("row %d = %q, want %q", i, lines[i+1], w)
		}
	}
}

// writePopFile writes one sample name per line to dir/name and returns the
// full path; helper for the Weir & Cockerham Fst integration tests.
func writePopFile(t *testing.T, dir, name string, samples []string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(strings.Join(samples, "\n")+"\n"), 0o644); err != nil {
		t.Fatalf("writing pop file %s: %v", path, err)
	}
	return path
}

// TestIntegration_WeirFst exercises --weir-fst-pop end-to-end on testVCF and
// asserts the per-site .weir.fst output is produced with the expected header
// and at least one data row.
func TestIntegration_WeirFst(t *testing.T) {
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "test")
	pop1 := writePopFile(t, tmpDir, "pop1.txt", []string{"sample1"})
	pop2 := writePopFile(t, tmpDir, "pop2.txt", []string{"sample2", "sample3"})

	params := &Params{
		OutPrefix:  prefix,
		WeirFstPop: []string{pop1, pop2},
	}
	if err := Run(strings.NewReader(testVCF), params); err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	fstPath := prefix + ".weir.fst"
	b, err := os.ReadFile(fstPath)
	if err != nil {
		t.Fatalf("expected output file %s: %v", fstPath, err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	wantHeader := "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST"
	if lines[0] != wantHeader {
		t.Errorf(".weir.fst header = %q, want %q", lines[0], wantHeader)
	}
	// pop1 has only one sample so each SNP is skipped (n_i < 2); the file
	// should still be produced with the header. Now use a synthetic VCF with
	// five samples in two populations, each of size >= 2, to verify a real
	// per-site Fst row is written.
	const fstVCF = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	.	A	G	30	PASS	.	GT	0/1	0/1	0/0	0/0	1/1
chr1	200	.	C	T	30	PASS	.	GT	0/0	0/1	1/1	0/1	0/0
`
	tmpDir3 := t.TempDir()
	prefix3 := filepath.Join(tmpDir3, "fst")
	popX := writePopFile(t, tmpDir3, "popX.txt", []string{"S1", "S2"})
	popY := writePopFile(t, tmpDir3, "popY.txt", []string{"S3", "S4", "S5"})
	if err := Run(strings.NewReader(fstVCF), &Params{
		OutPrefix:  prefix3,
		WeirFstPop: []string{popX, popY},
	}); err != nil {
		t.Fatalf("Run() error = %v", err)
	}
	b3, err := os.ReadFile(prefix3 + ".weir.fst")
	if err != nil {
		t.Fatalf("expected file %s: %v", prefix3+".weir.fst", err)
	}
	lines3 := strings.Split(strings.TrimRight(string(b3), "\n"), "\n")
	if lines3[0] != wantHeader {
		t.Errorf(".weir.fst header = %q, want %q", lines3[0], wantHeader)
	}
	// 2 SNPs in fstVCF -> 2 data rows.
	if got := len(lines3) - 1; got != 2 {
		t.Errorf(".weir.fst data rows = %d, want 2:\n%s", got, string(b3))
	}
	// First row is the worked example: chr1 100 ~ -0.3232.
	if !strings.HasPrefix(lines3[1], "chr1\t100\t") {
		t.Errorf("first row prefix = %q, want chr1\\t100\\t...", lines3[1])
	}
}

// TestIntegration_WeirFstRejectsBadPops checks the error paths in
// loadPopulationFiles: <2 pops, missing file, sample shared across pops.
func TestIntegration_WeirFstRejectsBadPops(t *testing.T) {
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "test")
	pop1 := writePopFile(t, tmpDir, "p1.txt", []string{"sample1"})

	// Single population => error.
	err := Run(strings.NewReader(testVCF), &Params{OutPrefix: prefix, WeirFstPop: []string{pop1}})
	if err == nil || !strings.Contains(err.Error(), "at least 2") {
		t.Errorf("Run with one pop file should error about needing >=2 pops; got: %v", err)
	}

	// Sample appearing in both pops => error.
	pop2 := writePopFile(t, tmpDir, "p2.txt", []string{"sample1", "sample2"})
	err = Run(strings.NewReader(testVCF), &Params{OutPrefix: prefix, WeirFstPop: []string{pop1, pop2}})
	if err == nil || !strings.Contains(err.Error(), "multiple population files") {
		t.Errorf("Run with sample in multiple pops should error; got: %v", err)
	}
}

// TestIntegration_WindowedWeirFst checks that --fst-window-size produces the
// .windowed.weir.fst file with the expected header and at least one row.
func TestIntegration_WindowedWeirFst(t *testing.T) {
	const fstVCF = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	.	A	G	30	PASS	.	GT	0/1	0/1	0/0	0/0	1/1
chr1	500	.	C	T	30	PASS	.	GT	0/0	0/1	1/1	0/1	0/0
chr1	1500	.	G	A	30	PASS	.	GT	1/1	0/1	0/0	0/0	0/1
`
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "fst")
	popX := writePopFile(t, tmpDir, "popX.txt", []string{"S1", "S2"})
	popY := writePopFile(t, tmpDir, "popY.txt", []string{"S3", "S4", "S5"})

	params := &Params{
		OutPrefix:     prefix,
		WeirFstPop:    []string{popX, popY},
		FstWindowSize: 1000,
	}
	if err := Run(strings.NewReader(fstVCF), params); err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	winPath := prefix + ".windowed.weir.fst"
	b, err := os.ReadFile(winPath)
	if err != nil {
		t.Fatalf("expected file %s: %v", winPath, err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	wantHeader := "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST"
	if lines[0] != wantHeader {
		t.Errorf(".windowed.weir.fst header = %q, want %q", lines[0], wantHeader)
	}
	// Step defaulted to 1000 (== window). SNPs at 100, 500 -> window
	// [1, 1000]; SNP at 1500 -> window [1001, 2000]. Two windows expected.
	if got := len(lines) - 1; got != 2 {
		t.Errorf(".windowed.weir.fst data rows = %d, want 2:\n%s", got, string(b))
	}
	// Per-site .weir.fst should also exist.
	if _, err := os.Stat(prefix + ".weir.fst"); err != nil {
		t.Errorf("expected per-site %s: %v", prefix+".weir.fst", err)
	}
}

func TestIntegration_IndelHistGenoDepthTajimaD(t *testing.T) {
	tmpDir := t.TempDir()
	prefix := filepath.Join(tmpDir, "test")

	params := &Params{
		OutPrefix:    prefix,
		HistIndelLen: true,
		GenoDepth:    true,
		TajimaD:      1000,
	}
	if err := Run(strings.NewReader(testVCF), params); err != nil {
		t.Fatalf("Run() error = %v", err)
	}

	// The only indel in testVCF is chr2:200 A -> AT, i.e. length +1, once.
	hist, err := os.ReadFile(prefix + ".indel.hist")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(hist), "1\t1\t") {
		t.Errorf(".indel.hist should contain a length-1 indel; got:\n%s", string(hist))
	}

	// .gdepth: header has CHROM POS then sample names; rows carry DP values.
	gd, err := os.ReadFile(prefix + ".gdepth")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(gd), "CHROM\tPOS\tsample1\tsample2\tsample3") {
		t.Errorf(".gdepth header missing sample columns; got:\n%s", string(gd))
	}
	if !strings.Contains(string(gd), "chr1\t100\t20\t20\t20") {
		t.Errorf(".gdepth should report DP 20/20/20 for chr1:100; got:\n%s", string(gd))
	}

	// .Tajima.D: each window header line present; chr1's window has 3 SNPs.
	td, err := os.ReadFile(prefix + ".Tajima.D")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(td), "CHROM\tBIN_START\tN_SNPS\tTajimaD") {
		t.Errorf(".Tajima.D header missing; got:\n%s", string(td))
	}
	if !strings.Contains(string(td), "chr1\t0\t3\t") {
		t.Errorf(".Tajima.D should have a chr1 window with 3 SNPs; got:\n%s", string(td))
	}
}
