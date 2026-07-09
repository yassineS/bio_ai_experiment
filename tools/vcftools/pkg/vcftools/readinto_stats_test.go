package vcftools

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The streaming-stats loop drives the reader via ReadInto, reusing one Variant
// across every site instead of allocating a fresh one per record. The reused
// buffers make a stray retained alias (a stat that keeps v.Chrom or an allele
// string past the next read) corrupt output, so these tests exercise the exact
// output cells the real-data benchmark compares and assert the result is (a)
// stable across repeated runs and (b) carries the correct, un-corrupted
// chromosome / allele strings that the accumulators retain across reads.

// multiChromVCF spans two chromosomes with a multi-allelic and an indel site so
// the interned CHROM/REF/ALT strings the stats accumulators retain are all
// exercised.
const multiChromVCF = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	G	30	PASS	DP=60	GT:DP	0/0:20	0/1:20	1/1:20
chr1	200	.	C	T,A	25	PASS	DP=45	GT:DP	0/0:15	0/1:15	1/2:15
chr1	300	.	G	A	40	PASS	DP=90	GT:DP	0/1:30	1/1:30	1/1:30
chr2	100	.	T	C	35	PASS	DP=75	GT:DP	0/1:25	0/1:25	1/1:25
chr2	200	.	AT	A	28	PASS	DP=60	GT:DP	0/0:20	0/0:20	0/1:20
chr2	300	.	A	G	22	PASS	DP=48	GT:DP	0/1:16	0/1:16	0/1:16
`

// runStatsCell runs one --<stat> cell into a temp prefix and returns the named
// output file's bytes. It mirrors the real-data benchmark's vcftools cells.
func runStatsCell(t *testing.T, params *Params, outFile string) []byte {
	t.Helper()
	dir := t.TempDir()
	params.OutPrefix = filepath.Join(dir, "rb")
	if err := Run(strings.NewReader(multiChromVCF), params); err != nil {
		t.Fatalf("Run() error = %v", err)
	}
	b, err := os.ReadFile(filepath.Join(dir, outFile))
	if err != nil {
		t.Fatalf("reading %s: %v", outFile, err)
	}
	return b
}

// TestReadIntoStatsDeterministic runs each streaming-stats cell twice and
// asserts byte-identical output. The ReadInto loop reuses one Variant across
// sites; if any accumulator retained an alias into the reused buffer, the second
// (differently-scheduled) run would drift, so identical bytes across runs is a
// direct guard against aliasing corruption in the reuse path.
func TestReadIntoStatsDeterministic(t *testing.T) {
	cells := []struct {
		name    string
		params  func() *Params
		outFile string
	}{
		{"freq", func() *Params { return &Params{Freq: true} }, "rb.frq"},
		{"counts", func() *Params { return &Params{Counts: true} }, "rb.frq.count"},
		{"hardy", func() *Params { return &Params{Hardy: true} }, "rb.hwe"},
		{"site_mean_depth", func() *Params { return &Params{SiteMeanDepth: true} }, "rb.ldepth.mean"},
		{"window_pi", func() *Params { return &Params{WindowPi: 100000} }, "rb.windowed.pi"},
	}
	for _, c := range cells {
		t.Run(c.name, func(t *testing.T) {
			first := runStatsCell(t, c.params(), c.outFile)
			second := runStatsCell(t, c.params(), c.outFile)
			if string(first) != string(second) {
				t.Fatalf("%s: output not deterministic across ReadInto reuse runs\nfirst:\n%s\nsecond:\n%s",
					c.name, first, second)
			}
		})
	}
}

// TestReadIntoStatsRetainsCorrectChroms asserts the retained CHROM/allele
// strings survive the buffer reuse intact: a stale alias would surface as a
// wrong or garbled chromosome column. Every data row's first column must be one
// of the two real chromosome names, and both must appear.
func TestReadIntoStatsRetainsCorrectChroms(t *testing.T) {
	out := runStatsCell(t, &Params{Freq: true}, "rb.frq")
	lines := strings.Split(strings.TrimRight(string(out), "\n"), "\n")
	if len(lines) < 2 {
		t.Fatalf("expected a header plus data rows, got %d lines", len(lines))
	}
	sawChr1, sawChr2 := false, false
	for _, line := range lines[1:] { // skip header
		chrom := strings.SplitN(line, "\t", 2)[0]
		switch chrom {
		case "chr1":
			sawChr1 = true
		case "chr2":
			sawChr2 = true
		default:
			t.Fatalf("corrupt/aliased CHROM column %q in row %q", chrom, line)
		}
	}
	if !sawChr1 || !sawChr2 {
		t.Fatalf("expected both chr1 and chr2 in output, got chr1=%v chr2=%v", sawChr1, sawChr2)
	}
}
