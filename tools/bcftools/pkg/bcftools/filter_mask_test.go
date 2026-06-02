package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const maskTestVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	T	.	PASS	DP=10	GT	0/1
chr1	250	.	C	G	.	PASS	DP=20	GT	0/1
chr1	300	.	GAT	G	.	PASS	DP=15;INDEL	GT	0/1
chr2	100	.	A	T	.	PASS	DP=10	GT	0/1
`

// TestVCFFilterMaskRegion: --mask chr1:1-200 tags only chr1:100 with the
// soft-filter ID; everything else stays PASS.
func TestVCFFilterMaskRegion(t *testing.T) {
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(maskTestVCF), &out, VCFFilterOptions{
		SoftFilter: "MASKED",
		MaskRegion: "chr1:1-200",
		NoVersion:  true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	s := out.String()
	if !strings.Contains(s, "chr1\t100\t.\tA\tT\t.\tMASKED") {
		t.Errorf("chr1:100 not masked:\n%s", s)
	}
	for _, want := range []string{"chr1\t250\t.\tC\tG\t.\tPASS", "chr1\t300\t.\tGAT\tG\t.\tPASS", "chr2\t100\t.\tA\tT\t.\tPASS"} {
		if !strings.Contains(s, want) {
			t.Errorf("expected PASS line missing %q in:\n%s", want, s)
		}
	}
	if !strings.Contains(s, `##FILTER=<ID=MASKED,Description="Record masked by region">`) {
		t.Errorf("missing mask FILTER header:\n%s", s)
	}
}

// TestVCFFilterMaskNegate: `^chr1:1-200` masks everything *outside* the
// region — chr1:250, chr1:300, chr2:100 all get the MASKED tag.
func TestVCFFilterMaskNegate(t *testing.T) {
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(maskTestVCF), &out, VCFFilterOptions{
		SoftFilter: "MASKED",
		MaskRegion: "^chr1:1-200",
		NoVersion:  true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	s := out.String()
	for _, want := range []string{"chr1\t250\t.\tC\tG\t.\tMASKED", "chr1\t300\t.\tGAT\tG\t.\tMASKED", "chr2\t100\t.\tA\tT\t.\tMASKED"} {
		if !strings.Contains(s, want) {
			t.Errorf("expected MASKED line missing %q in:\n%s", want, s)
		}
	}
	if !strings.Contains(s, "chr1\t100\t.\tA\tT\t.\tPASS") {
		t.Errorf("chr1:100 should have stayed PASS:\n%s", s)
	}
}

// TestVCFFilterMaskFile: BED end is half-open in BED but inclusive in
// the resolved mask; "0\t100" should cover POS=1..100 inclusive.
func TestVCFFilterMaskFile(t *testing.T) {
	dir := t.TempDir()
	bed := filepath.Join(dir, "m.bed")
	if err := os.WriteFile(bed, []byte("chr1\t0\t100\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(maskTestVCF), &out, VCFFilterOptions{
		SoftFilter: "MASKED",
		MaskFile:   bed,
		NoVersion:  true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	s := out.String()
	if !strings.Contains(s, "chr1\t100\t.\tA\tT\t.\tMASKED") {
		t.Errorf("BED end inclusive: chr1:100 should be masked:\n%s", s)
	}
}

// TestVCFFilterMaskRequiresSoftFilter: upstream's exact error text.
func TestVCFFilterMaskRequiresSoftFilter(t *testing.T) {
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(maskTestVCF), &out, VCFFilterOptions{
		MaskRegion: "chr1:1-200",
		NoVersion:  true,
	})
	if err == nil {
		t.Fatal("expected error for --mask without -s")
	}
	want := "The option --soft-filter is required with --mask and --mask-file options"
	if !strings.Contains(err.Error(), want) {
		t.Errorf("error text mismatch:\n got: %v\nwant substring: %q", err, want)
	}
}

// TestVCFFilterMaskOverlap exercises the three --mask-overlap modes
// against a BED region that touches the chr1:300 INDEL (REF=GAT, span
// [300,302] under overlap=1).
func TestVCFFilterMaskOverlap(t *testing.T) {
	dir := t.TempDir()
	bed := filepath.Join(dir, "m.bed")
	// BED 300-302 -> 1-based inclusive [301,302].
	if err := os.WriteFile(bed, []byte("chr1\t300\t302\n"), 0644); err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		overlap int
		want    bool // chr1:300 GAT masked?
		note    string
	}{
		{0, false, "POS=300 outside [301,302]"},
		{1, true, "record span [300,302] overlaps [301,302]"},
		{2, true, "variant span [300,302] overlaps [301,302]"},
	}
	for _, tc := range cases {
		t.Run(tc.note, func(t *testing.T) {
			var out bytes.Buffer
			if _, err := VCFFilter(strings.NewReader(maskTestVCF), &out, VCFFilterOptions{
				SoftFilter:  "MASKED",
				MaskFile:    bed,
				MaskOverlap: tc.overlap,
				NoVersion:   true,
			}); err != nil {
				t.Fatalf("VCFFilter: %v", err)
			}
			masked := strings.Contains(out.String(), "chr1\t300\t.\tGAT\tG\t.\tMASKED")
			if masked != tc.want {
				t.Errorf("overlap=%d: chr1:300 masked=%v, want %v\n%s", tc.overlap, masked, tc.want, out.String())
			}
		})
	}
}
