package bedgenomecov

import (
	"bytes"
	"strings"
	"testing"
)

// genomeOf parses a chrom-sizes-style string for tests.
func genomeOf(t *testing.T, s string) *GenomeSize {
	t.Helper()
	g, err := ReadGenome(strings.NewReader(s))
	if err != nil {
		t.Fatalf("ReadGenome: %v", err)
	}
	return g
}

func runOf(t *testing.T, bed, genome string, opts Options) string {
	t.Helper()
	g := genomeOf(t, genome)
	var buf bytes.Buffer
	if err := Run(strings.NewReader(bed), g, &buf, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return buf.String()
}

func TestReadGenome(t *testing.T) {
	cases := []struct {
		name    string
		in      string
		wantErr bool
		order   []string
	}{
		{"ok", "chr1\t10\nchr2\t5\n", false, []string{"chr1", "chr2"}},
		{"spaces", "chr1 10\nchr2 5\n", false, []string{"chr1", "chr2"}},
		{"comments", "# hello\nchr1\t10\n", false, []string{"chr1"}},
		{"empty", "", true, nil},
		{"missing size", "chr1\n", true, nil},
		{"bad size", "chr1\tabc\n", true, nil},
		{"zero size", "chr1\t0\n", true, nil},
		{"duplicate", "chr1\t10\nchr1\t20\n", true, nil},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			g, err := ReadGenome(strings.NewReader(c.in))
			if (err != nil) != c.wantErr {
				t.Fatalf("err=%v wantErr=%v", err, c.wantErr)
			}
			if !c.wantErr {
				if len(g.Order) != len(c.order) {
					t.Fatalf("order=%v want %v", g.Order, c.order)
				}
				for i, chrom := range c.order {
					if g.Order[i] != chrom {
						t.Errorf("order[%d]=%s want %s", i, g.Order[i], chrom)
					}
				}
			}
		})
	}
}

func TestHistogramBasic(t *testing.T) {
	// chr1 length 10. Two intervals: [0,3) and [2,5).
	// Coverage: positions 0,1 depth 1, 2 depth 2, 3,4 depth 1, 5..9 depth 0.
	// => depth 0 -> 5 bases, depth 1 -> 4 bases, depth 2 -> 1 base.
	out := runOf(t, "chr1\t0\t3\nchr1\t2\t5\n", "chr1\t10\n", Options{})
	want := "chr1\t0\t5\t10\t0.5\n" +
		"chr1\t1\t4\t10\t0.4\n" +
		"chr1\t2\t1\t10\t0.1\n" +
		"genome\t0\t5\t10\t0.5\n" +
		"genome\t1\t4\t10\t0.4\n" +
		"genome\t2\t1\t10\t0.1\n"
	if out != want {
		t.Errorf("histogram mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestRunBAM_SAMInput(t *testing.T) {
	// SAM text is accepted by RunBAM (sam.NewReader auto-detects SAM/BAM).
	// Genome comes from the @SQ header. r2 is a spliced 5M10N5M read.
	sam := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:100\n" +
		"r1\t0\tchr1\t11\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t0\tchr1\t16\t60\t5M10N5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"

	// Without -split: r1 covers [10,20), r2 covers [15,30) (whole span).
	var nosplit bytes.Buffer
	if err := RunBAM(strings.NewReader(sam), &nosplit, Options{Mode: ModeBedGraphAll}); err != nil {
		t.Fatalf("RunBAM: %v", err)
	}
	ns := nosplit.String()
	if !strings.Contains(ns, "chr1\t15\t20\t2\n") {
		t.Errorf("no-split: expected depth-2 overlap [15,20):\n%s", ns)
	}

	// With -split: r2 contributes only its blocks [15,20) and [30,35), so the
	// intron+gap [20,30) drops to depth 0.
	var split bytes.Buffer
	if err := RunBAM(strings.NewReader(sam), &split, Options{Mode: ModeBedGraphAll, Split: true}); err != nil {
		t.Fatalf("RunBAM split: %v", err)
	}
	sp := split.String()
	if !strings.Contains(sp, "chr1\t20\t30\t0\n") {
		t.Errorf("split: expected the intron [20,30) at depth 0:\n%s", sp)
	}
	if !strings.Contains(sp, "chr1\t30\t35\t1\n") {
		t.Errorf("split: expected r2 second block [30,35) at depth 1:\n%s", sp)
	}
}

// pairedSAM is a single proper pair: the +strand mate at pos 11 with TLEN 30,
// and its -strand mate at pos 31 with TLEN -30.
const pairedSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:100\n" +
	"p1\t99\tchr1\t11\t60\t10M\t=\t31\t30\tACGTACGTAC\tIIIIIIIIII\n" +
	"p1\t147\tchr1\t31\t60\t10M\t=\t11\t-30\tACGTACGTAC\tIIIIIIIIII\n"

func TestRunBAM_PairedCoverage(t *testing.T) {
	// -pc: the fragment spans [10,40) (pos0 10 + TLEN 30), counted once.
	var buf bytes.Buffer
	if err := RunBAM(strings.NewReader(pairedSAM), &buf, Options{Mode: ModeBedGraphAll, PairedCoverage: true}); err != nil {
		t.Fatalf("RunBAM -pc: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, "chr1\t10\t40\t1\n") {
		t.Errorf("-pc: expected single fragment [10,40) at depth 1:\n%s", out)
	}
}

func TestRunBAM_FragmentSize(t *testing.T) {
	// -fs 20: the +read covers [10,30), the -read covers [20,40) (anchored at
	// its 3' end), so [20,30) is depth 2.
	var buf bytes.Buffer
	if err := RunBAM(strings.NewReader(pairedSAM), &buf, Options{Mode: ModeBedGraphAll, FragmentSize: 20}); err != nil {
		t.Fatalf("RunBAM -fs: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, "chr1\t20\t30\t2\n") {
		t.Errorf("-fs 20: expected the fragment overlap [20,30) at depth 2:\n%s", out)
	}
}

func TestRunBAM_PCandFSExclusive(t *testing.T) {
	if err := RunBAM(strings.NewReader(pairedSAM), &bytes.Buffer{},
		Options{PairedCoverage: true, FragmentSize: 20}); err == nil {
		t.Fatal("expected error combining -pc and -fs")
	}
}

func TestRunBAM_NoSQHeader(t *testing.T) {
	if err := RunBAM(strings.NewReader("@HD\tVN:1.6\n"), &bytes.Buffer{}, Options{}); err == nil {
		t.Fatal("expected error when the header has no @SQ entries")
	}
}

func TestSplitBED12Blocks(t *testing.T) {
	// One BED12 record, 3 blocks of 10 at starts 0/20/40. Without -split the
	// whole [0,50) span is covered; with -split only the blocks are.
	rec := "chr1\t0\t50\tx\t0\t+\t0\t0\t0\t3\t10,10,10,\t0,20,40,\n"
	genome := "chr1\t60\n"

	whole := runOf(t, rec, genome, Options{Mode: ModeBedGraphAll})
	if !strings.Contains(whole, "chr1\t0\t50\t1\n") {
		t.Errorf("no-split should cover the whole span:\n%s", whole)
	}

	split := runOf(t, rec, genome, Options{Mode: ModeBedGraphAll, Split: true})
	for _, blk := range []string{"chr1\t0\t10\t1\n", "chr1\t20\t30\t1\n", "chr1\t40\t50\t1\n"} {
		if !strings.Contains(split, blk) {
			t.Errorf("split missing block %q in:\n%s", blk, split)
		}
	}
	// The inter-block gaps must be depth 0 under -split.
	if !strings.Contains(split, "chr1\t10\t20\t0\n") {
		t.Errorf("split should leave the inter-block gap at depth 0:\n%s", split)
	}
}

func TestHistogramMaxDepth(t *testing.T) {
	// depth 5 covered at pos 0 (5 intervals); cap to 2.
	bed := "chr1\t0\t1\nchr1\t0\t1\nchr1\t0\t1\nchr1\t0\t1\nchr1\t0\t1\n"
	out := runOf(t, bed, "chr1\t3\n", Options{MaxDepth: 2})
	want := "chr1\t0\t2\t3\t0.666667\n" +
		"chr1\t2\t1\t3\t0.333333\n" +
		"genome\t0\t2\t3\t0.666667\n" +
		"genome\t2\t1\t3\t0.333333\n"
	if out != want {
		t.Errorf("max-depth mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestBedGraphNoZero(t *testing.T) {
	out := runOf(t, "chr1\t0\t3\nchr1\t2\t5\n", "chr1\t10\n", Options{Mode: ModeBedGraph})
	want := "chr1\t0\t2\t1\n" +
		"chr1\t2\t3\t2\n" +
		"chr1\t3\t5\t1\n"
	if out != want {
		t.Errorf("bg mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestBedGraphWithZero(t *testing.T) {
	out := runOf(t, "chr1\t0\t3\n", "chr1\t6\n", Options{Mode: ModeBedGraphAll})
	want := "chr1\t0\t3\t1\n" +
		"chr1\t3\t6\t0\n"
	if out != want {
		t.Errorf("bga mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestPerBase(t *testing.T) {
	out := runOf(t, "chr1\t0\t2\n", "chr1\t3\n", Options{Mode: ModePerBase})
	want := "chr1\t1\t1\n" +
		"chr1\t2\t1\n" +
		"chr1\t3\t0\n"
	if out != want {
		t.Errorf("per-base mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestPerBaseNonZero(t *testing.T) {
	// Upstream `-dz` reports ZERO-based positions (offset 0), unlike `-d`
	// (1-based, offset 1). bedtools genomecov -i ... -dz on [0,2) emits
	// `chr1 0 1` / `chr1 1 1`.
	out := runOf(t, "chr1\t0\t2\n", "chr1\t3\n", Options{Mode: ModePerBaseNonZero})
	want := "chr1\t0\t1\n" +
		"chr1\t1\t1\n"
	if out != want {
		t.Errorf("per-base nz mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestStrandFilter(t *testing.T) {
	// Only + intervals should be counted.
	bed := "chr1\t0\t3\ta\t0\t+\nchr1\t2\t5\tb\t0\t-\n"
	out := runOf(t, bed, "chr1\t10\n", Options{Mode: ModeBedGraph, Strand: "+"})
	want := "chr1\t0\t3\t1\n"
	if out != want {
		t.Errorf("strand+ mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}

	out = runOf(t, bed, "chr1\t10\n", Options{Mode: ModeBedGraph, Strand: "-"})
	want = "chr1\t2\t5\t1\n"
	if out != want {
		t.Errorf("strand- mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestStrandValidation(t *testing.T) {
	g := genomeOf(t, "chr1\t10\n")
	var buf bytes.Buffer
	err := Run(strings.NewReader(""), g, &buf, Options{Strand: "x"})
	if err == nil {
		t.Fatal("expected error for invalid strand")
	}
}

func TestFivePrimeAndThreePrimeIncompatible(t *testing.T) {
	g := genomeOf(t, "chr1\t10\n")
	var buf bytes.Buffer
	err := Run(strings.NewReader(""), g, &buf, Options{FivePrime: true, ThreePrime: true})
	if err == nil {
		t.Fatal("expected error for combining -5 and -3")
	}
}

func TestFivePrime(t *testing.T) {
	// + interval contributes position 0; - interval contributes position 4 (end-1).
	bed := "chr1\t0\t3\ta\t0\t+\nchr1\t2\t5\tb\t0\t-\n"
	out := runOf(t, bed, "chr1\t10\n", Options{Mode: ModePerBase, FivePrime: true})
	// Pos 1 has depth 1 (5'+), pos 5 has depth 1 (5'-).
	want := "chr1\t1\t1\nchr1\t2\t0\nchr1\t3\t0\nchr1\t4\t0\nchr1\t5\t1\nchr1\t6\t0\nchr1\t7\t0\nchr1\t8\t0\nchr1\t9\t0\nchr1\t10\t0\n"
	if out != want {
		t.Errorf("five-prime mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestThreePrime(t *testing.T) {
	// + interval contributes its last base (index 2 for [0,3)); - interval
	// contributes its first base (index 2 for [2,5)); both land on the same
	// base, depth 2. Reported with -dz, which is ZERO-based (offset 0), so the
	// position is 2 — matching `bedtools genomecov ... -dz -3`.
	bed := "chr1\t0\t3\ta\t0\t+\nchr1\t2\t5\tb\t0\t-\n"
	out := runOf(t, bed, "chr1\t5\n", Options{Mode: ModePerBaseNonZero, ThreePrime: true})
	want := "chr1\t2\t2\n"
	if out != want {
		t.Errorf("three-prime mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestScaleBedGraph(t *testing.T) {
	out := runOf(t, "chr1\t0\t3\n", "chr1\t3\n", Options{Mode: ModeBedGraph, Scale: 0.5})
	want := "chr1\t0\t3\t0.5\n"
	if out != want {
		t.Errorf("scale mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestScalePerBase(t *testing.T) {
	out := runOf(t, "chr1\t0\t1\n", "chr1\t2\n", Options{Mode: ModePerBase, Scale: 2.0})
	want := "chr1\t1\t2\nchr1\t2\t0\n"
	if out != want {
		t.Errorf("scale mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestTrackline(t *testing.T) {
	out := runOf(t, "chr1\t0\t3\n", "chr1\t3\n",
		Options{Mode: ModeBedGraph, TrackLine: true, TrackOpts: "name=test"})
	want := "track type=bedGraph name=test\nchr1\t0\t3\t1\n"
	if out != want {
		t.Errorf("trackline mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}

	// Trackline only added in bedGraph modes (silently ignored elsewhere).
	out = runOf(t, "chr1\t0\t3\n", "chr1\t3\n", Options{Mode: ModeHistogram, TrackLine: true})
	if strings.Contains(out, "track") {
		t.Errorf("histogram should not get a trackline:\n%s", out)
	}
}

func TestIntervalClampedToGenome(t *testing.T) {
	// Interval extends past end-of-chrom; should be clamped.
	out := runOf(t, "chr1\t8\t20\n", "chr1\t10\n", Options{Mode: ModeBedGraph})
	want := "chr1\t8\t10\t1\n"
	if out != want {
		t.Errorf("clamp mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestUnknownChromosomeIgnored(t *testing.T) {
	// chrX not in genome file -> silently skipped.
	out := runOf(t, "chrX\t0\t3\nchr1\t0\t1\n", "chr1\t2\n", Options{Mode: ModeBedGraph})
	want := "chr1\t0\t1\t1\n"
	if out != want {
		t.Errorf("unknown chrom should be ignored.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestEmptyInput(t *testing.T) {
	out := runOf(t, "", "chr1\t3\n", Options{Mode: ModeBedGraphAll})
	want := "chr1\t0\t3\t0\n"
	if out != want {
		t.Errorf("empty input mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}

func TestNilGenomeError(t *testing.T) {
	var buf bytes.Buffer
	if err := Run(strings.NewReader(""), nil, &buf, Options{}); err == nil {
		t.Fatal("expected error for nil genome")
	}
	g := &GenomeSize{}
	if err := Run(strings.NewReader(""), g, &buf, Options{}); err == nil {
		t.Fatal("expected error for empty genome")
	}
}

func TestBadBedRecord(t *testing.T) {
	g := genomeOf(t, "chr1\t10\n")
	var buf bytes.Buffer
	err := Run(strings.NewReader("chr1\tnan\t3\n"), g, &buf, Options{})
	if err == nil {
		t.Fatal("expected error for malformed input")
	}
}

func TestSortedKeysAndFormatFraction(t *testing.T) {
	keys := sortedKeys(map[int]int{3: 1, 1: 1, 2: 1})
	if keys[0] != 1 || keys[1] != 2 || keys[2] != 3 {
		t.Errorf("sortedKeys=%v", keys)
	}
	if s := formatFraction(0.5); s != "0.5" {
		t.Errorf("formatFraction(0.5)=%s", s)
	}
	if s := formatFraction(1.0 / 3.0); s != "0.333333" {
		t.Errorf("formatFraction(1/3)=%s", s)
	}
}

func TestMultipleChromosomes(t *testing.T) {
	bed := "chr1\t0\t2\nchr2\t1\t3\n"
	out := runOf(t, bed, "chr1\t5\nchr2\t5\n", Options{Mode: ModeBedGraph})
	want := "chr1\t0\t2\t1\n" +
		"chr2\t1\t3\t1\n"
	if out != want {
		t.Errorf("multi-chrom mismatch.\nwant:\n%s\ngot:\n%s", want, out)
	}
}
