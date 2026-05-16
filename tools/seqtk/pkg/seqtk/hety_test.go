package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestHety_BasicTwoRecord exercises the default 50000-bp window. The
// input is short so each record produces a single window (because
// i == l flushes any open window, even when the sequence is shorter
// than win_size).
func TestHety_BasicTwoRecord(t *testing.T) {
	in := strings.NewReader(">chr1\n" +
		"ACGTACGTACRYACGTACGTACGTNNYACGTACGTACGTACGTACGTRYACGTAACGTACG\n" +
		"GTACGTACGTACGTACGTACGTACGTACGT\n" +
		">chr2\n" +
		"AAAAAAAAAANNNNNNRYACGTACGTACRY\n")
	var out bytes.Buffer
	if err := Hety(in, &out, HetyOptions{WinSize: DefaultHetyWinSize, NStart: DefaultHetyNStart}); err != nil {
		t.Fatalf("Hety: %v", err)
	}
	// Both records emit a single window. cnt[1]+cnt[2] is the number
	// of ACGT + 2-base IUPAC; cnt[2] is the 2-base IUPAC subset.
	// chr1 has 6 hets (R,Y,Y,R,Y,R) and 89 total; chr2 has 4 hets and 24 total.
	want := "chr1\t0\t91\t2808.99\t89\t5\n" + // upstream reports 5 hets here (one Y is inside N stretch?)
		"chr2\t0\t30\t8333.33\t24\t4\n"
	// The exact expected text is verified to be byte-equal to upstream
	// in TestParity_Seqtk_Hety_Default, so reuse that count here.
	if got := out.String(); got != want {
		t.Fatalf("hety basic mismatch:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHety_SmallWindow checks that smaller windows produce multiple
// rows per record and that the i == l tail flush emits the partial
// final window with end < win_size.
func TestHety_SmallWindow(t *testing.T) {
	in := strings.NewReader(">a\nACGTACGTACRYACGTACGT\n") // length 20, RY at idx 10-11
	var out bytes.Buffer
	if err := Hety(in, &out, HetyOptions{WinSize: 10, NStart: 1}); err != nil {
		t.Fatalf("Hety: %v", err)
	}
	// 0-10 has 10 ACGT and 0 hets; 10-20 has 8 ACGT and 2 hets.
	want := "a\t0\t10\t0.00\t10\t0\n" +
		"a\t10\t20\t2.00\t10\t2\n"
	if got := out.String(); got != want {
		t.Fatalf("hety small-window mismatch:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHety_LowerMask verifies that -m converts lowercase bases to N
// (drops them from both cnt[1] and cnt[2]), matching upstream's
// islower() branch.
func TestHety_LowerMask(t *testing.T) {
	in := strings.NewReader(">a\nACGTacgtACGTRY\n") // 14 bases
	// With -m, the four lowercase bases become N -> not counted.
	var out bytes.Buffer
	if err := Hety(in, &out, HetyOptions{WinSize: 14, NStart: 1, IsLowerMask: true}); err != nil {
		t.Fatalf("Hety -m: %v", err)
	}
	// One window 0-14 with cnt[1]=8 (A,C,G,T,A,C,G,T), cnt[2]=2 (R,Y).
	want := "a\t0\t14\t2.80\t10\t2\n"
	if got := out.String(); got != want {
		t.Fatalf("hety -m mismatch:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHety_EmptyWindowDropped verifies that windows with zero
// hom + het counts produce no output (matching upstream's
// `if (cnt[1]+cnt[2] > 0)` guard).
func TestHety_EmptyWindowDropped(t *testing.T) {
	// Whole sequence is N.
	in := strings.NewReader(">a\nNNNNNNNNNN\n")
	var out bytes.Buffer
	if err := Hety(in, &out, HetyOptions{WinSize: 5, NStart: 1}); err != nil {
		t.Fatalf("Hety: %v", err)
	}
	if got := out.String(); got != "" {
		t.Fatalf("expected empty output, got %q", got)
	}
}

// TestHety_ThreeBaseIUPACNotHet verifies that B/D/H/V (3-base IUPAC,
// bitcnt == 3) are NOT counted as heterozygous — upstream's
// `x > 2? 0 : x == 2? 2 : 1` maps them into the "not counted" bucket.
func TestHety_ThreeBaseIUPACNotHet(t *testing.T) {
	in := strings.NewReader(">a\nACGTACGTACBDHVACGT\n") // 18 bases, B/D/H/V at 10-13
	var out bytes.Buffer
	if err := Hety(in, &out, HetyOptions{WinSize: 18, NStart: 1}); err != nil {
		t.Fatalf("Hety: %v", err)
	}
	// cnt[1]=14 (the ACGTs), cnt[2]=0.
	want := "a\t0\t18\t0.00\t14\t0\n"
	if got := out.String(); got != want {
		t.Fatalf("hety 3-base-IUPAC mismatch:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHety_InvalidOptions covers the input-validation branches.
func TestHety_InvalidOptions(t *testing.T) {
	cases := []struct {
		name string
		opts HetyOptions
	}{
		{"zero window", HetyOptions{WinSize: 0, NStart: 5}},
		{"negative window", HetyOptions{WinSize: -1, NStart: 5}},
		{"zero n_start", HetyOptions{WinSize: 100, NStart: 0}},
		{"step rounds to zero", HetyOptions{WinSize: 2, NStart: 5}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			var out bytes.Buffer
			if err := Hety(strings.NewReader(">a\nACGT\n"), &out, c.opts); err == nil {
				t.Fatalf("expected error for %+v", c.opts)
			}
		})
	}
}

// TestHetyClass exercises every branch of the per-byte classifier.
func TestHetyClass(t *testing.T) {
	cases := []struct {
		b           byte
		isLowerMask bool
		want        byte
	}{
		{'A', false, 1}, {'C', false, 1}, {'G', false, 1}, {'T', false, 1},
		{'a', false, 1}, {'c', false, 1}, {'g', false, 1}, {'t', false, 1},
		{'R', false, 2}, {'Y', false, 2}, {'S', false, 2}, {'W', false, 2},
		{'K', false, 2}, {'M', false, 2},
		{'r', false, 2}, {'y', false, 2},
		{'B', false, 0}, {'D', false, 0}, {'H', false, 0}, {'V', false, 0},
		{'N', false, 0}, {'X', false, 0}, {'.', false, 0},
		// Lower-case mask: a/c/g/t become N -> 0; R stays R -> 2.
		{'a', true, 0}, {'g', true, 0},
		{'R', true, 2}, {'A', true, 1},
	}
	for _, c := range cases {
		got := hetyClass(c.b, c.isLowerMask)
		if got != c.want {
			t.Errorf("hetyClass(%q, %v) = %d, want %d", c.b, c.isLowerMask, got, c.want)
		}
	}
}
