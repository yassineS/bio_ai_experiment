package bedcoverage

import (
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Binary-free unit tests for the two confirmed bedcoverage divergences:
//
//  1. the overlap-fraction (-f/-F/-r/-e) pass/skip decision, including the rule
//     that -split suppresses fraction filtering entirely; and
//  2. the verbatim BED12 blockSizes/blockStarts echo.
//
// They use synthetic *bed.Record values and call the pure helpers directly, so
// they pass with an UNPOPULATED reference_code/bedtools submodule.

// TestUnitFractionPassNonSplit checks -f/-F/-r/-e semantics on a single B
// against a single A in the non-split path (fractionPass is the decision).
func TestUnitFractionPassNonSplit(t *testing.T) {
	// A = [0,100); B = [90,110): overlap 10 -> fracA=0.10, fracB=0.50.
	a := plainB("chr1", 0, 100)
	b := plainB("chr1", 90, 110)
	cases := []struct {
		name string
		opts Options
		want bool
	}{
		{"no thresholds", Options{}, true},
		{"f pass", Options{FractionA: 0.05}, true},
		{"f fail", Options{FractionA: 0.5}, false},
		{"F pass", Options{FractionB: 0.5}, true},
		{"F fail", Options{FractionB: 0.6}, false},
		{"fF AND both pass", Options{FractionA: 0.05, FractionB: 0.5}, true},
		{"fF AND one fails", Options{FractionA: 0.5, FractionB: 0.5}, false},
		{"e OR one passes", Options{FractionA: 0.5, FractionB: 0.5, Either: true}, true},
		{"e OR both fail", Options{FractionA: 0.5, FractionB: 0.6, Either: true}, false},
		// -r forces the B-side threshold to equal FractionA: fracA=0.1 < 0.5 fails.
		{"r reject", Options{FractionA: 0.5, FractionB: 0.5, Reciprocal: true}, false},
		// -r where the A fraction holds on both sides: A=[0,20) over B=[0,20) full.
		{"r accept", Options{FractionA: 1.0, Reciprocal: true}, true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			la, lb := a, b
			if tc.name == "r accept" {
				la, lb = plainB("chr1", 0, 20), plainB("chr1", 0, 20)
			}
			if got := fractionPass(la, lb, tc.opts); got != tc.want {
				t.Fatalf("fractionPass = %v, want %v", got, tc.want)
			}
		})
	}
}

// TestUnitSplitSuppressesFractionFilter is the core of divergence #1: under
// -split, selectOverlapping must NOT drop a B for failing -f/-F, because
// upstream coverage keeps the always-populated overlap set. A 1bp B over a
// 50bp block fails -f 0.5 in the non-split path but is kept under -split.
func TestUnitSplitSuppressesFractionFilter(t *testing.T) {
	// Query: chr1 100-400, blocks 100-150 and 350-400 (two 50bp blocks).
	q := blockedQuery("chr1", 100, 400, [2]int{0, 50}, [2]int{250, 50})
	// 1bp B inside block1 -> fracA(over block) = 0.02, way below -f 0.5.
	tree := bed.NewIntervalTree([]*bed.Record{plainB("chr1", 100, 101)})

	nonSplit := selectOverlapping(q, tree, Options{FractionA: 0.5})
	if len(nonSplit) != 0 {
		t.Fatalf("non-split -f 0.5: kept %d, want 0 (fraction filter active)", len(nonSplit))
	}

	for _, opts := range []Options{
		{Split: true, FractionA: 0.5},
		{Split: true, FractionA: 1.0},
		{Split: true, FractionB: 1.0},
		{Split: true, FractionA: 1.0, FractionB: 1.0},
		{Split: true, FractionA: 1.0, Reciprocal: true},
		{Split: true, FractionA: 1.0, FractionB: 1.0, Either: true},
	} {
		got := selectOverlapping(q, tree, opts)
		if len(got) != 1 {
			t.Fatalf("split %+v: kept %d, want 1 (fractions ignored under -split)", opts, len(got))
		}
	}
}

// TestUnitVerbatimBlockEcho is the core of divergence #2: recordColumns must
// echo a BED12 record's blockSizes/blockStarts exactly as read (raw text wins),
// preserving or omitting the trailing comma; BAM-style records with no raw text
// fall back to the trailing-comma form.
func TestUnitVerbatimBlockEcho(t *testing.T) {
	base := bed.Record{
		Chrom: "chr1", ChromStart: 100, ChromEnd: 400,
		Name: "q", Strand: "+", ThickStart: 100, ThickEnd: 400,
		BlockCount: 2, BlockSizes: []int{50, 50}, BlockStarts: []int{0, 250},
	}

	noComma := base
	noComma.RawBlockSizes = "50,50"
	noComma.RawBlockStarts = "0,250"
	cols := recordColumns(&noComma)
	if got := strings.Join(cols, "\t"); !strings.HasSuffix(got, "50,50\t0,250") {
		t.Fatalf("no-comma echo = %q, want block cols 50,50 / 0,250 verbatim", got)
	}

	withComma := base
	withComma.RawBlockSizes = "50,50,"
	withComma.RawBlockStarts = "0,250,"
	cols = recordColumns(&withComma)
	if got := strings.Join(cols, "\t"); !strings.HasSuffix(got, "50,50,\t0,250,") {
		t.Fatalf("with-comma echo = %q, want trailing comma preserved", got)
	}

	// No raw text (BAM-derived): fall back to trailing-comma rendering.
	bam := base
	bam.ItemRGB = "0,0,0"
	cols = recordColumns(&bam)
	if got := strings.Join(cols, "\t"); !strings.HasSuffix(got, "50,50,\t0,250,") {
		t.Fatalf("BAM-fallback echo = %q, want trailing-comma form", got)
	}
}

// TestUnitBlockField checks the verbatim/fallback selection in isolation.
func TestUnitBlockField(t *testing.T) {
	if got := blockField("50,50", []int{50, 50}); got != "50,50" {
		t.Fatalf("blockField(raw no comma) = %q, want verbatim 50,50", got)
	}
	if got := blockField("50,50,", []int{50, 50}); got != "50,50," {
		t.Fatalf("blockField(raw comma) = %q, want verbatim 50,50,", got)
	}
	if got := blockField("", []int{50, 50}); got != "50,50," {
		t.Fatalf("blockField(no raw) = %q, want trailing-comma fallback 50,50,", got)
	}
}
