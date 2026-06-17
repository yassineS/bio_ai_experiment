package bcftools

// Live upstream-binary parity for two divergences the parity pipeline flagged
// on its synthetic fixtures:
//
//  1. consensus default het handling — upstream applies IUPAC ambiguity codes
//     (consensus.c iupac_GTs) when the VCF carries samples and neither -H nor
//     an allele pick is given, across the selected sample (-s) or all samples.
//     Our port previously emitted a resolved base.
//
//  2. roh ST/RG record ordering — upstream emits the tables chromosome-major
//     then sample-major (header order), not sample-major.
//
// Both run the genuine upstream binary and our port on the same bgzipped+
// indexed fixtures and assert byte-equality (data rows, provenance stripped).
// t.Fatalf (never t.Skip) when a binary is unavailable.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

const consHetFasta = ">chr1\nACGTACGTACGTACGTACGT\n"

// consHetVCF carries het and hom genotypes across two samples so the IUPAC
// default folds C+A->M and A+G->R.
const consHetVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	2	.	C	A	.	.	.	GT	0/0	1/1
chr1	5	.	A	G	.	.	.	GT	0/0	0/1
chr1	9	.	A	C	.	.	.	GT	0/1	0/1
`

// writeFastaFaidx writes a FASTA and indexes it with the vendored samtools
// faidx (upstream consensus needs the .fai). Returns the FASTA path.
func writeFastaFaidx(t *testing.T, htslibDir, dir, name, content string) string {
	t.Helper()
	fa := filepath.Join(dir, name)
	if err := os.WriteFile(fa, []byte(content), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	// bcftools/htslib ship no faidx; use the sibling samtools binary under
	// reference_code/samtools (htslibDir is reference_code/htslib).
	samtools := filepath.Join(filepath.Dir(htslibDir), "samtools", "samtools")
	if out, err := exec.Command(samtools, "faidx", fa).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx %s: %v: %s", fa, err, out)
	}
	return fa
}

// rohOrderVCF: 3 samples across 2 chromosomes, enough sites for the HMM to
// emit ST rows, so the chromosome-major/sample-major ordering is exercised.
var rohOrderVCF = func() string {
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##contig=<ID=chr1,length=100000>\n")
	b.WriteString("##contig=<ID=chr2,length=100000>\n")
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n")
	gts := []string{"0/0", "0/1", "1/1"}
	for _, chrom := range []string{"chr1", "chr2"} {
		for i := 0; i < 30; i++ {
			pos := 100 + i*137
			g := func(j int) string { return gts[(i+j)%3] }
			b.WriteString(chrom + "\t")
			b.WriteString(itoaTest(pos) + "\t.\tA\tG\t.\t.\t.\tGT\t")
			b.WriteString(g(0) + "\t" + g(1) + "\t" + g(2) + "\n")
		}
	}
	return b.String()
}()

func itoaTest(n int) string {
	if n == 0 {
		return "0"
	}
	var d []byte
	for n > 0 {
		d = append([]byte{byte('0' + n%10)}, d...)
		n /= 10
	}
	return string(d)
}

// TestConsensus_DefaultIUPAC_UpstreamParity checks the iupac_GTs default
// against the upstream binary: all-samples (no -s) and single-sample (-s).
func TestConsensus_DefaultIUPAC_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)
	tools := upstreamBcftoolsConsensus(t)
	htslibDir := filepath.Dir(tools.bgzip)

	dir := t.TempDir()
	fa := writeFastaFaidx(t, htslibDir, dir, "ref.fa", consHetFasta)
	gz := writeMergeInput(t, htslibDir, dir, "het.vcf", consHetVCF)

	cases := [][]string{
		{"consensus", "-f", fa, gz},             // all samples -> IUPAC
		{"consensus", "-s", "S2", "-f", fa, gz}, // one sample -> IUPAC
	}
	for i, args := range cases {
		t.Run([]string{"all-samples", "one-sample"}[i], func(t *testing.T) {
			want := runBin(t, live, args...)
			got := runBin(t, ours, args...)
			if !bytes.Equal(want, got) {
				t.Fatalf("consensus IUPAC mismatch\n--- upstream ---\n%s\n--- ours ---\n%s", want, got)
			}
		})
	}
}

// TestRoh_RecordOrder_UpstreamParity checks the chromosome-major/sample-major
// ST/RG ordering on a multi-sample, multi-chromosome fixture.
func TestRoh_RecordOrder_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)
	tools := upstreamBcftoolsConsensus(t)
	htslibDir := filepath.Dir(tools.bgzip)

	dir := t.TempDir()
	gz := writeMergeInput(t, htslibDir, dir, "roh.vcf", rohOrderVCF)

	args := []string{"roh", "-G30", "--AF-dflt", "0.4", gz}
	want := stripProvenanceBytes(runBin(t, live, args...))
	got := stripProvenanceBytes(runBin(t, ours, args...))
	if !bytes.Equal(want, got) {
		// Show only the first divergence to keep the failure readable.
		t.Fatalf("roh record ordering mismatch vs upstream\n--- upstream (head) ---\n%.800s\n--- ours (head) ---\n%.800s", want, got)
	}
}
