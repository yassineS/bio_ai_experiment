package vcftools

// Tests for `--bcf` (BCF input) and `--contigs` (supplemental contig
// declarations for BCF header construction). PR-G chunk 2 follow-up
// to the wave-21 `--recode-bcf` PR.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestRun_BCFInput_Roundtrip exercises the full --recode-bcf →
// --bcf pipeline: write a BCF via --recode-bcf, read it back through
// the port's BCF input path, recode to VCF, and assert the resulting
// VCF carries every variant + per-sample field from the source.
func TestRun_BCFInput_Roundtrip(t *testing.T) {
	tmp := t.TempDir()

	// Step 1: VCF → BCF via --recode-bcf.
	bcfPrefix := filepath.Join(tmp, "step1")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:     bcfPrefix,
		RecodeBCF:     true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("Run --recode-bcf: %v", err)
	}

	// Step 2: BCF → VCF via --bcf + --recode.
	vcfPrefix := filepath.Join(tmp, "step2")
	if err := Run(io.NopCloser(strings.NewReader("")), &Params{
		OutPrefix:     vcfPrefix,
		BCFInputFile:  bcfPrefix + ".recode.bcf",
		Recode:        true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("Run --bcf --recode: %v", err)
	}

	body, err := os.ReadFile(vcfPrefix + ".recode.vcf")
	if err != nil {
		t.Fatal(err)
	}
	out := string(body)

	// Every data line from the original fixture (modulo INFO key
	// ordering, which our writer can't pin) must appear in the round-
	// tripped output. Use position+REF+ALT as a stable handle and check
	// the per-sample tail (genotype + scalars) verbatim.
	wantRows := []string{
		// chr1:100 — full GT:DP:GQ:PL trail for both samples.
		"chr1\t100\t.\tA\tG\t30\tPASS\t",
		"0/0:8:60:0,15,200",
		"0/1:12:80:50,0,180",
		// chr1:200
		"chr1\t200\t.\tC\tT\t25\tPASS\t",
		"0/1:7:55:30,0,150",
		"1/1:9:70:200,15,0",
	}
	for _, w := range wantRows {
		if !strings.Contains(out, w) {
			t.Errorf("output missing %q:\n%s", w, out)
		}
	}
}

// TestRun_BCFInput_ComposesWithFilters verifies that --bcf input
// flows through the normal vcftools filter pipeline. We use --chr +
// --from-bp on a 2-site BCF fixture and assert the expected single
// surviving variant.
func TestRun_BCFInput_ComposesWithFilters(t *testing.T) {
	tmp := t.TempDir()
	bcfPrefix := filepath.Join(tmp, "in")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:     bcfPrefix,
		RecodeBCF:     true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("Run --recode-bcf: %v", err)
	}
	outPrefix := filepath.Join(tmp, "filtered")
	if err := Run(io.NopCloser(strings.NewReader("")), &Params{
		OutPrefix:     outPrefix,
		BCFInputFile:  bcfPrefix + ".recode.bcf",
		Chr:           "chr1",
		FromBp:        150,
		Recode:        true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("Run --bcf --chr --from-bp: %v", err)
	}
	body, err := os.ReadFile(outPrefix + ".recode.vcf")
	if err != nil {
		t.Fatal(err)
	}
	out := string(body)
	if !strings.Contains(out, "chr1\t200\t") {
		t.Errorf("expected chr1:200 in output:\n%s", out)
	}
	if strings.Contains(out, "chr1\t100\t") {
		t.Errorf("--from-bp 150 should drop chr1:100; got:\n%s", out)
	}
}

// TestRun_ContigsFile_AddsContigLines verifies that --contigs prepends
// ##contig=<ID=X> meta lines to a header that lacks contig declarations.
// Mirrors upstream's gating on `meta_data.has_contigs == false`.
func TestRun_ContigsFile_AddsContigLines(t *testing.T) {
	tmp := t.TempDir()
	contigsPath := filepath.Join(tmp, "contigs.txt")
	if err := os.WriteFile(contigsPath, []byte("chr1\nchr2\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	const noContigsFixture = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
chr1	100	.	A	G	30	PASS	DP=10	GT	0/0
chr2	250	.	C	T	30	PASS	DP=15	GT	0/1
`
	outPrefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(noContigsFixture), &Params{
		OutPrefix:     outPrefix,
		ContigsFile:   contigsPath,
		Recode:        true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("Run --contigs --recode: %v", err)
	}
	body, err := os.ReadFile(outPrefix + ".recode.vcf")
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{"##contig=<ID=chr1>", "##contig=<ID=chr2>"} {
		if !bytes.Contains(body, []byte(want)) {
			t.Errorf("output missing %q:\n%s", want, body)
		}
	}
}

// TestRun_ContigsFile_NoOpWhenHeaderAlreadyHasContigs verifies the
// upstream gating: if the source header already declares any
// `##contig=` line, --contigs is silently ignored (matching upstream
// variant_file.cpp:154's `if (has_contigs == false)`).
func TestRun_ContigsFile_NoOpWhenHeaderAlreadyHasContigs(t *testing.T) {
	tmp := t.TempDir()
	contigsPath := filepath.Join(tmp, "contigs.txt")
	if err := os.WriteFile(contigsPath, []byte("chr99\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	// recodeBCFFixture (from recode_bcf_test.go) already has
	// ##contig=<ID=chr1>/chr2 declarations.
	outPrefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(recodeBCFFixture), &Params{
		OutPrefix:     outPrefix,
		ContigsFile:   contigsPath,
		Recode:        true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatal(err)
	}
	body, err := os.ReadFile(outPrefix + ".recode.vcf")
	if err != nil {
		t.Fatal(err)
	}
	if bytes.Contains(body, []byte("##contig=<ID=chr99>")) {
		t.Errorf("--contigs should be ignored when header has contigs; got chr99 line:\n%s", body)
	}
}

// TestAugmentHeaderContigs_AcceptsMetaInfoForm checks that --contigs
// lines starting with `##contig=<` are kept verbatim (not re-wrapped).
func TestAugmentHeaderContigs_AcceptsMetaInfoForm(t *testing.T) {
	tmp := t.TempDir()
	contigsPath := filepath.Join(tmp, "contigs.txt")
	body := "##contig=<ID=chrM,length=16569>\n##contig=<ID=chrX,length=155270560>\n"
	if err := os.WriteFile(contigsPath, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	const noContigsFixture = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
chrM	100	.	A	G	30	PASS	DP=10	GT	0/0
`
	outPrefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(noContigsFixture), &Params{
		OutPrefix:     outPrefix,
		ContigsFile:   contigsPath,
		Recode:        true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatal(err)
	}
	out, err := os.ReadFile(outPrefix + ".recode.vcf")
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{
		"##contig=<ID=chrM,length=16569>",
		"##contig=<ID=chrX,length=155270560>",
	} {
		if !bytes.Contains(out, []byte(want)) {
			t.Errorf("output missing %q:\n%s", want, out)
		}
	}
}
