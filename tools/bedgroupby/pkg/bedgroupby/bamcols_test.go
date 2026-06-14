package bedgroupby

import (
	"bytes"
	"strings"
	"testing"
)

// samFixture is a tiny SAM-text alignment used to exercise the BAM/SAM input
// path without needing a binary BAM. It reproduces the column layout upstream
// `bedtools groupby` groups over: each alignment becomes
// QNAME, FLAG, RNAME, 0-based start, MAPQ, CIGAR(op-then-length), ...
const samFixture = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
read	0	chr1	11	60	3M2I5M	*	0	0	ACGTTACGTA	IIIIIIIIII
read	16	chr1	51	30	4M1D6M	*	0	0	ACGTACGTAC	IIIIIIIIII
read	0	chr1	101	40	5M	*	0	0	ACGTA	IIIII
`

// TestSAMColumns checks the per-record rendering matches the exact columns
// upstream bedtools' BamRecord::getField exposes: col4 is the 0-based start
// and col6 is the CIGAR with the op character before its length ("3M2I5M" ->
// "M3I2M5"). These were confirmed against reference_code/bedtools.
func TestSAMColumns(t *testing.T) {
	src, err := newLineSource(strings.NewReader(samFixture))
	if err != nil {
		t.Fatalf("newLineSource: %v", err)
	}
	want := []string{
		"read\t0\tchr1\t10\t60\tM3I2M5\t.\t-1\t0\tACGTTACGTA\tIIIIIIIIII",
		"read\t16\tchr1\t50\t30\tM4D1M6\t.\t-1\t0\tACGTACGTAC\tIIIIIIIIII",
		"read\t0\tchr1\t100\t40\tM5\t.\t-1\t0\tACGTA\tIIIII",
	}
	var got []string
	for {
		line, ok, err := src.next()
		if err != nil {
			t.Fatalf("next: %v", err)
		}
		if !ok {
			break
		}
		got = append(got, line)
	}
	if len(got) != len(want) {
		t.Fatalf("got %d lines, want %d: %q", len(got), len(want), got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("line %d:\n got %q\nwant %q", i, got[i], want[i])
		}
	}
}

// TestGroup_SAMInput drives the full grouping engine over SAM text, mirroring
// the upstream BAM test (`groupby -g 1,3 -c 4 -o mean`): group by read name
// and chrom, mean of the 0-based starts. All three reads share name+chrom, so
// the mean is (10+50+100)/3.
func TestGroup_SAMInput(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Group(strings.NewReader(samFixture), &buf, Options{
		GroupCols: []int{1, 3},
		AggCols:   []int{4},
		Ops:       []string{"mean"},
	}); err != nil {
		t.Fatalf("Group: %v", err)
	}
	got := strings.TrimSpace(buf.String())
	// (10+50+100)/3 = 53.3333333
	if !strings.HasPrefix(got, "read\tchr1\t53.33") {
		t.Errorf("got %q", got)
	}
}

// TestGroup_SAMInput_DistinctCigar groups on the chrom column and collapses
// the CIGAR column (6), confirming the op-then-length rendering survives the
// grouping engine end to end.
func TestGroup_SAMInput_DistinctCigar(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Group(strings.NewReader(samFixture), &buf, Options{
		GroupCols: []int{3},
		AggCols:   []int{6},
		Ops:       []string{"collapse"},
	}); err != nil {
		t.Fatalf("Group: %v", err)
	}
	got := strings.TrimSpace(buf.String())
	want := "chr1\tM3I2M5,M4D1M6,M5"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}
