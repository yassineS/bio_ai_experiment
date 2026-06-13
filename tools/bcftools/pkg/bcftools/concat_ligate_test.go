package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// The header shared by all ligate fixtures. Upstream's `view -Oz` step adds a
// `##FILTER=<ID=PASS,...>` line and non-reproducible ##bcftools_* provenance
// lines; stripLigateProvenance removes both so the comparison is on the
// reproducible payload (the records and the PS/PQ-augmented FORMAT/sample data).
const ligateHdr = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=20>\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"

const ligateColumns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"

// writeLigateFixture writes a VCF fixture to dir/name and returns its path.
func writeLigateFixture(t *testing.T, dir, name, body string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(ligateHdr+ligateColumns+body), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return path
}

// bgzipIndexLigate compresses src to src.gz and writes a .tbi index using the
// upstream bcftools binary (concat --ligate requires indexed input). It returns
// the path to the .gz file.
func bgzipIndexLigate(t *testing.T, bin, src string) string {
	t.Helper()
	gz := src + ".gz"
	if out, err := exec.Command(bin, "view", "-Oz", "-o", gz, src).CombinedOutput(); err != nil {
		t.Fatalf("bgzip %s: %v\n%s", src, err, out)
	}
	if out, err := exec.Command(bin, "index", "-f", "-t", gz).CombinedOutput(); err != nil {
		t.Fatalf("index %s: %v\n%s", gz, err, out)
	}
	return gz
}

// stripLigateProvenance removes lines that upstream injects but that are not
// reproducible (the ##bcftools_* provenance lines and the ##FILTER=<ID=PASS>
// line added by the bgzip view step). The remaining text is identical between
// the Go port and upstream.
func stripLigateProvenance(s string) string {
	var keep []string
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "##bcftools_") {
			continue
		}
		if strings.HasPrefix(line, "##FILTER=<ID=PASS") {
			continue
		}
		keep = append(keep, line)
	}
	return strings.Join(keep, "\n")
}

// runGoLigate runs the Go ConcatFiles in ligate mode over plain-text inputs.
func runGoLigate(t *testing.T, opts ConcatOptions, paths ...string) string {
	t.Helper()
	var buf bytes.Buffer
	opts.Ligate = true
	if _, err := ConcatFiles(paths, &buf, opts); err != nil {
		t.Fatalf("Go ConcatFiles ligate: %v", err)
	}
	return stripLigateProvenance(buf.String())
}

// runUpstreamLigate runs upstream `bcftools concat -l` over the bgzipped+indexed
// versions of the given plain-text inputs.
func runUpstreamLigate(t *testing.T, bin string, extra []string, plainPaths ...string) string {
	t.Helper()
	args := []string{"concat", "-l"}
	args = append(args, extra...)
	for _, p := range plainPaths {
		args = append(args, bgzipIndexLigate(t, bin, p))
	}
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream concat -l: %v\n%s", err, stderr.String())
	}
	return stripLigateProvenance(stdout.String())
}

// TestConcatLigateParityNoSwap covers a clean overlap where the trailing chunk
// is already in phase with the leading one: no sample should be swapped.
func TestConcatLigateParityNoSwap(t *testing.T) {
	bin := upstreamBcftools(t)
	dir := t.TempDir()
	c1 := writeLigateFixture(t, dir, "c1.vcf",
		"20\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t200\t.\tC\tT\t.\t.\t.\tGT\t1|0\t0|1\n"+
			"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t1|0\t0|1\n")
	c2 := writeLigateFixture(t, dir, "c2.vcf",
		"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t1|0\t0|1\n"+
			"20\t500\t.\tA\tT\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t600\t.\tG\tC\t.\t.\t.\tGT\t1|0\t0|1\n")

	want := runUpstreamLigate(t, bin, nil, c1, c2)
	got := runGoLigate(t, ConcatOptions{MinPQ: 30}, c1, c2)
	if want != got {
		t.Fatalf("ligate no-swap mismatch:\n--- upstream ---\n%s\n--- go ---\n%s", want, got)
	}
}

// TestConcatLigateParitySwap covers an overlap where the trailing chunk's
// haplotypes are flipped relative to the leading chunk; upstream must swap them
// back and the swap must propagate to the trailing chunk's later records.
func TestConcatLigateParitySwap(t *testing.T) {
	bin := upstreamBcftools(t)
	dir := t.TempDir()
	c1 := writeLigateFixture(t, dir, "c1.vcf",
		"20\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t200\t.\tC\tT\t.\t.\t.\tGT\t1|0\t0|1\n"+
			"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t1|0\t0|1\n")
	// c2 overlap (300,400) has S1 flipped; later sites 500,600 are in the
	// flipped frame and must be un-flipped on output.
	c2 := writeLigateFixture(t, dir, "c2swap.vcf",
		"20\t300\t.\tG\tA\t.\t.\t.\tGT\t1|0\t1|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t0|1\t0|1\n"+
			"20\t500\t.\tA\tT\t.\t.\t.\tGT\t1|0\t1|0\n"+
			"20\t600\t.\tG\tC\t.\t.\t.\tGT\t0|1\t0|1\n")

	want := runUpstreamLigate(t, bin, nil, c1, c2)
	got := runGoLigate(t, ConcatOptions{MinPQ: 30}, c1, c2)
	if want != got {
		t.Fatalf("ligate swap mismatch:\n--- upstream ---\n%s\n--- go ---\n%s", want, got)
	}
}

// TestConcatLigateParityLowPQ covers a noisy overlap where one sample's votes
// disagree, producing a phase quality below the default min-PQ of 30 and thus
// starting a new phase set for that sample at the boundary.
func TestConcatLigateParityLowPQ(t *testing.T) {
	bin := upstreamBcftools(t)
	dir := t.TempDir()
	p1 := writeLigateFixture(t, dir, "p1.vcf",
		"20\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|1\n"+
			"20\t200\t.\tC\tT\t.\t.\t.\tGT\t0|1\t0|1\n"+
			"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t0|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t0|1\t0|1\n"+
			"20\t500\t.\tA\tT\t.\t.\t.\tGT\t0|1\t0|1\n")
	// Overlap (300,400,500): S1 consistent (match,match), S2 mixed (match at
	// 300, mismatch at 400 and 500) -> nmatch=1,nmism=2 -> swap + low PQ (8).
	p2 := writeLigateFixture(t, dir, "p2.vcf",
		"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t0|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t500\t.\tA\tT\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t600\t.\tG\tC\t.\t.\t.\tGT\t0|1\t0|1\n")

	want := runUpstreamLigate(t, bin, nil, p1, p2)
	got := runGoLigate(t, ConcatOptions{MinPQ: 30}, p1, p2)
	if want != got {
		t.Fatalf("ligate low-PQ mismatch:\n--- upstream ---\n%s\n--- go ---\n%s", want, got)
	}

	// With -q 5 the low PQ (8) no longer breaks the phase set.
	wantQ5 := runUpstreamLigate(t, bin, []string{"-q", "5"}, p1, p2)
	gotQ5 := runGoLigate(t, ConcatOptions{MinPQ: 5}, p1, p2)
	if wantQ5 != gotQ5 {
		t.Fatalf("ligate low-PQ -q5 mismatch:\n--- upstream ---\n%s\n--- go ---\n%s", wantQ5, gotQ5)
	}
}

// TestConcatLigateParityThreeChunks covers a chain of three overlapping chunks,
// exercising the carry-forward of the resolved overlap into the next pairwise
// ligation.
func TestConcatLigateParityThreeChunks(t *testing.T) {
	bin := upstreamBcftools(t)
	dir := t.TempDir()
	c1 := writeLigateFixture(t, dir, "t1.vcf",
		"20\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t200\t.\tC\tT\t.\t.\t.\tGT\t1|0\t0|1\n"+
			"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t1|0\t0|1\n")
	c2 := writeLigateFixture(t, dir, "t2.vcf",
		"20\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|1\n"+
			"20\t400\t.\tT\tC\t.\t.\t.\tGT\t1|0\t0|1\n"+
			"20\t500\t.\tA\tT\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t600\t.\tG\tC\t.\t.\t.\tGT\t1|0\t0|1\n")
	c3 := writeLigateFixture(t, dir, "t3.vcf",
		"20\t500\t.\tA\tT\t.\t.\t.\tGT\t0|1\t1|0\n"+
			"20\t600\t.\tG\tC\t.\t.\t.\tGT\t1|0\t0|1\n"+
			"20\t700\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1|0\n")

	want := runUpstreamLigate(t, bin, nil, c1, c2, c3)
	got := runGoLigate(t, ConcatOptions{MinPQ: 30}, c1, c2, c3)
	if want != got {
		t.Fatalf("ligate three-chunk mismatch:\n--- upstream ---\n%s\n--- go ---\n%s", want, got)
	}
}

// --- focused unit tests for the vote/swap/PQ helpers ---

// mkLigateVariant builds a minimal phased variant for the unit tests.
func mkLigateVariant(chrom string, pos int, ref, alt string, gts ...string) *vcf.Variant {
	v := &vcf.Variant{Chrom: chrom, Pos: pos, Ref: ref, Alt: []string{alt}, Format: []string{"GT"}}
	for i, g := range gts {
		v.Samples = append(v.Samples, vcf.Sample{Name: "S" + ligateItoa(i+1), Data: map[string]string{"GT": g}})
	}
	return v
}

func ligateItoa(i int) string { return string(rune('0' + i)) }

func TestLigateVoteHetMatch(t *testing.T) {
	st := &ligateState{nsmpl: 1, swapPhase: []int{0}, nmatch: []int{0}, nmism: []int{0}}
	a := mkLigateVariant("20", 100, "A", "G", "0|1")
	b := mkLigateVariant("20", 100, "A", "G", "0|1")
	st.vote(a, b)
	if st.nmatch[0] != 1 || st.nmism[0] != 0 {
		t.Fatalf("expected match: nmatch=%d nmism=%d", st.nmatch[0], st.nmism[0])
	}
}

func TestLigateVoteHetMismatch(t *testing.T) {
	st := &ligateState{nsmpl: 1, swapPhase: []int{0}, nmatch: []int{0}, nmism: []int{0}}
	a := mkLigateVariant("20", 100, "A", "G", "0|1")
	b := mkLigateVariant("20", 100, "A", "G", "1|0")
	st.vote(a, b)
	if st.nmatch[0] != 0 || st.nmism[0] != 1 {
		t.Fatalf("expected mismatch: nmatch=%d nmism=%d", st.nmatch[0], st.nmism[0])
	}
}

func TestLigateVoteSkipsHomUnphasedMissing(t *testing.T) {
	st := &ligateState{nsmpl: 1, swapPhase: []int{0}, nmatch: []int{0}, nmism: []int{0}}
	cases := [][2]string{
		{"1|1", "0|1"}, // homozygous a
		{"0|1", "0/1"}, // unphased b
		{"0|1", ".|1"}, // missing b
		{"0/1", "0|1"}, // unphased a
	}
	for _, c := range cases {
		st.nmatch[0], st.nmism[0] = 0, 0
		st.vote(mkLigateVariant("20", 100, "A", "G", c[0]), mkLigateVariant("20", 100, "A", "G", c[1]))
		if st.nmatch[0] != 0 || st.nmism[0] != 0 {
			t.Fatalf("expected skip for %v: nmatch=%d nmism=%d", c, st.nmatch[0], st.nmism[0])
		}
	}
}

func TestLigateVoteUnderSwap(t *testing.T) {
	// With swap_phase[0]=1, an as-is match counts as a mismatch and vice versa.
	st := &ligateState{nsmpl: 1, swapPhase: []int{1}, nmatch: []int{0}, nmism: []int{0}}
	a := mkLigateVariant("20", 100, "A", "G", "0|1")
	b := mkLigateVariant("20", 100, "A", "G", "0|1") // as-is match
	st.vote(a, b)
	if st.nmism[0] != 1 || st.nmatch[0] != 0 {
		t.Fatalf("under swap, as-is match should be a mismatch: nmatch=%d nmism=%d", st.nmatch[0], st.nmism[0])
	}
}

func TestLigateDecideSwapAndPQ(t *testing.T) {
	st := &ligateState{nsmpl: 2, swapPhase: []int{0, 0}, nmatch: []int{2, 1}, nmism: []int{0, 2}, phaseQual: make([]int32, 2)}
	st.decideSwap()
	if st.swapPhase[0] != 0 {
		t.Fatalf("sample0 nmatch>=nmism so no swap, got %d", st.swapPhase[0])
	}
	if st.swapPhase[1] != 1 {
		t.Fatalf("sample1 nmism>nmatch so swap, got %d", st.swapPhase[1])
	}
	if st.phaseQual[0] != 99 {
		t.Fatalf("sample0 one-sided votes -> PQ 99, got %d", st.phaseQual[0])
	}
	// sample1: f=1/3 -> 99*(0.7 + f*ln f + (1-f)*ln(1-f))/0.7
	if st.phaseQual[1] != 8 {
		t.Fatalf("sample1 PQ from entropy formula expected 8, got %d", st.phaseQual[1])
	}
}

func TestLigatePhaseUpdateFlips(t *testing.T) {
	st := &ligateState{nsmpl: 2, swapPhase: []int{1, 0}}
	rec := mkLigateVariant("20", 100, "A", "G", "0|1", "1|0")
	st.phaseUpdate(rec)
	if got := rec.Samples[0].Data["GT"]; got != "1|0" {
		t.Fatalf("swapped sample should flip 0|1 -> 1|0, got %s", got)
	}
	if got := rec.Samples[1].Data["GT"]; got != "1|0" {
		t.Fatalf("unswapped sample should be untouched, got %s", got)
	}
}

func TestLigateBuildPairsAlignment(t *testing.T) {
	a := []*vcf.Variant{mkLigateVariant("20", 300, "G", "A", "0|1"), mkLigateVariant("20", 400, "T", "C", "1|0")}
	b := []*vcf.Variant{mkLigateVariant("20", 300, "G", "A", "0|1"), mkLigateVariant("20", 400, "T", "C", "1|0")}
	pairs := buildPairs(a, b)
	if len(pairs) != 2 {
		t.Fatalf("expected 2 paired sites, got %d", len(pairs))
	}
	for i, p := range pairs {
		if p.a == nil || p.b == nil {
			t.Fatalf("pair %d should have both sides", i)
		}
	}
}

func TestLigatePQFormula(t *testing.T) {
	// Spot-check a few (nmatch,nmism) -> PQ values against the entropy formula.
	cases := []struct {
		match, mism int
		want        int32
	}{
		{2, 0, 99},
		{0, 2, 99},
		{1, 1, 0},
		{1, 2, 8},
		{3, 1, 19},
	}
	for _, c := range cases {
		st := &ligateState{nsmpl: 1, swapPhase: []int{0}, nmatch: []int{c.match}, nmism: []int{c.mism}, phaseQual: make([]int32, 1)}
		st.decideSwap()
		if st.phaseQual[0] != c.want {
			t.Fatalf("PQ(%d,%d): want %d got %d", c.match, c.mism, c.want, st.phaseQual[0])
		}
	}
}
