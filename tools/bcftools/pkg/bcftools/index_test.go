package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// writeBCFForIndex builds a tiny BGZF-wrapped BCF on disk by piping a small
// VCF through View. The resulting file is what the index code path consumes.
func writeBCFForIndex(t *testing.T, vcfText string) string {
	t.Helper()
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		f.Close()
		t.Fatalf("View(-O b): %v", err)
	}
	f.Close()
	return bcfPath
}

// TestBuildIndexBCFCSI builds a CSI index for a freshly-emitted BCF file.
func TestBuildIndexBCFCSI(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr1	200	b	C	G	30	PASS	DP=20
chr2	50	c	G	A	30	PASS	DP=30
`
	bcfPath := writeBCFForIndex(t, vcfText)
	out, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	if _, err := os.Stat(out); err != nil {
		t.Fatalf("output index missing: %v", err)
	}
	if !strings.HasSuffix(out, ".csi") {
		t.Errorf("expected .csi suffix, got %s", out)
	}

	csi, err := tabix.ReadCSIFile(out)
	if err != nil {
		t.Fatalf("tabix.ReadCSIFile: %v", err)
	}
	// htslib computes n_lvls from the longest contig: with both contigs at
	// length 1000, hts_adjust_csi_settings yields depth 0 (1000+256 fits in
	// hts_bin_maxpos(14,0)=16384). Confirmed against genuine
	// `bcftools index -c`, which writes min_shift=14, n_lvls=0, l_aux=0.
	if csi.MinShift != 14 || csi.Depth != 0 {
		t.Errorf("unexpected params: %+v", csi)
	}
	if len(csi.Aux) != 0 {
		t.Errorf("BCF CSI must have no aux block, got %d bytes", len(csi.Aux))
	}
	// Should have at least 2 refs (chr1, chr2).
	if len(csi.Refs) < 2 {
		t.Errorf("got %d refs, want >= 2", len(csi.Refs))
	}
}

// TestBuildIndexMinShift exercises --csi-min-shift on a sparse dataset.
func TestBuildIndexMinShift(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
`
	bcfPath := writeBCFForIndex(t, vcfText)
	out, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, MinShift: 16})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	csi, err := tabix.ReadCSIFile(out)
	if err != nil {
		t.Fatal(err)
	}
	if csi.MinShift != 16 {
		t.Errorf("min_shift: got %d, want 16", csi.MinShift)
	}
}

// TestBuildIndexNoForceFails surfaces the "already exists" guard.
func TestBuildIndexNoForceFails(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI}); err != nil {
		t.Fatal(err)
	}
	// Second attempt without -f should error.
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI}); err == nil {
		t.Fatal("expected error on existing index without -f")
	}
	// With -f it should succeed.
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatalf("BuildIndex(force): %v", err)
	}
}

// TestBuildIndexTBIRejectsBCF asserts the .tbi route refuses BCF input.
func TestBuildIndexTBIRejectsBCF(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexTBI}); err == nil {
		t.Fatal("expected error: --tbi on BCF input")
	}
}

// TestBuildIndexVCFGzCSI exercises CSI generation on a bgzipped VCF.
func TestBuildIndexVCFGzCSI(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
chr1	150	b	C	G	30	PASS	.
`
	dir := t.TempDir()
	bgzPath := filepath.Join(dir, "y.vcf.gz")
	f, err := os.Create(bgzPath)
	if err != nil {
		t.Fatal(err)
	}
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(vcfText)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()

	out, err := BuildIndex(bgzPath, IndexOptions{Format: IndexCSI})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	csi, err := tabix.ReadCSIFile(out)
	if err != nil {
		t.Fatal(err)
	}
	if len(csi.Refs) == 0 {
		t.Fatal("expected at least one ref")
	}
	if csi.Names[0] != "chr1" {
		t.Errorf("names[0]: %q", csi.Names[0])
	}
}

// TestBuildIndexVCFGzTBI exercises the .tbi path on bgzipped VCF.
func TestBuildIndexVCFGzTBI(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
`
	dir := t.TempDir()
	bgzPath := filepath.Join(dir, "z.vcf.gz")
	f, err := os.Create(bgzPath)
	if err != nil {
		t.Fatal(err)
	}
	bw := bgzip.NewWriter(f)
	bw.Write([]byte(vcfText))
	bw.Close()
	f.Close()

	out, err := BuildIndex(bgzPath, IndexOptions{Format: IndexTBI})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	if !strings.HasSuffix(out, ".tbi") {
		t.Errorf("expected .tbi suffix, got %s", out)
	}
}

// TestViewRegionsWithCSI confirms the BCF .csi fast path: build a BCF + CSI,
// then query a region via ViewFile and verify the correct record returns.
func TestViewRegionsWithCSI(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr1	500	b	C	G	30	PASS	DP=20
chr1	900	c	G	A	30	PASS	DP=30
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ViewFile(bcfPath, &out, ViewOptions{Regions: []string{"chr1:50-200"}}, nil); err != nil {
		t.Fatalf("ViewFile: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "\ta\t") {
		t.Errorf("expected record 'a' in:\n%s", got)
	}
	if strings.Contains(got, "\tb\t") || strings.Contains(got, "\tc\t") {
		t.Errorf("records b/c should not be in chr1:50-200 result:\n%s", got)
	}
}

// TestLooksLikeBCF covers BCF detection for both BGZF-wrapped and raw input.
func TestLooksLikeBCF(t *testing.T) {
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(`##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	T	.	PASS	.
`), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	got, err := looksLikeBCF(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if !got {
		t.Error("BGZF-wrapped BCF not detected")
	}

	// Plain VCF should not match.
	vcfPath := filepath.Join(dir, "x.vcf")
	if err := os.WriteFile(vcfPath, []byte("##fileformat=VCFv4.2\n"), 0644); err != nil {
		t.Fatal(err)
	}
	if got, _ := looksLikeBCF(vcfPath); got {
		t.Error("plain VCF should not detect as BCF")
	}
}

// oracleDir is the path (relative to this package) to the committed genuine
// upstream CSI oracles produced by `bcftools index -c` / `tabix --csi -p vcf`.
const oracleDir = "../../../../pkg/htsgo/tabix/testdata/bcf_csi"

// TestCSIDepthAdjustBCF asserts our BCF .csi is byte-identical to one produced
// by genuine `bcftools index -c`. It exercises the depth computation
// (hts_adjust_csi_settings: depth 0 for short contigs) and the l_aux=0 (no aux
// block) behaviour htslib uses for BCF indexes.
func TestCSIDepthAdjustBCF(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(oracleDir, "in.bcf")
	bcfPath := filepath.Join(dir, "in.bcf")
	copyFileForIndex(t, src, bcfPath)

	out, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	got, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	want, err := os.ReadFile(filepath.Join(oracleDir, "expected.csi"))
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("BCF .csi not byte-equal to genuine `bcftools index -c`: got %d bytes, want %d", len(got), len(want))
	}
}

// TestCSIDepthAdjustVCFGz asserts our VCF.gz .csi is byte-identical to one
// produced by genuine `tabix --csi -p vcf`. The depth here is computed from the
// ##contig length (10000), and the tabix conf aux block is retained.
func TestCSIDepthAdjustVCFGz(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(oracleDir, "vcf_in.vcf.gz")
	vcfPath := filepath.Join(dir, "in.vcf.gz")
	copyFileForIndex(t, src, vcfPath)

	out, err := BuildIndex(vcfPath, IndexOptions{Format: IndexCSI, Force: true})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	got, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	want, err := os.ReadFile(filepath.Join(oracleDir, "vcf_expected.csi"))
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("VCF.gz .csi not byte-equal to genuine `tabix --csi -p vcf`: got %d bytes, want %d", len(got), len(want))
	}
}

// copyFileForIndex copies src to dst for the oracle round-trip tests.
func copyFileForIndex(t *testing.T, src, dst string) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read %s: %v", src, err)
	}
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatalf("write %s: %v", dst, err)
	}
}
