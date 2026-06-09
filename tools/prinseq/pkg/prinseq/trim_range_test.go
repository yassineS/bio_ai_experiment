package prinseq

// Unit tests for the sliding-window quality trim, trim_to_len, and the
// range_len / range_gc filters. These do not depend on the upstream Perl
// oracle (see trim_range_parity_test.go for the live byte-for-byte checks)
// and lock in the exact edge semantics ported from prinseq-lite.pl.

import (
	"bytes"
	"testing"
)

func TestCheckRange(t *testing.T) {
	cases := []struct {
		ranges string
		val    int
		want   bool
	}{
		{"10-20", 10, true},
		{"10-20", 20, true},
		{"10-20", 9, false},
		{"10-20", 21, false},
		{"1-100", 50, true},
		{"50-50", 50, true},
		{"50-50", 49, false},
		// Multiple ranges use AND semantics (upstream's checkRange returns
		// false on the first range the value falls outside of). With
		// disjoint ranges nothing passes.
		{"10-20,30-40", 15, false},
		{"10-20,30-40", 35, false},
		// Overlapping ranges: the value must satisfy both.
		{"10-30,20-40", 25, true},
		{"10-30,20-40", 15, false},
		// A bare number (no '-') gives an implicit upper bound of 0.
		{"5", 0, false},  // 0 < 5
		{"0-5", 0, true}, // explicit lower bound 0
		{"0-0", 0, true},
	}
	for _, c := range cases {
		if got := checkRange(c.ranges, c.val); got != c.want {
			t.Errorf("checkRange(%q, %d) = %v, want %v", c.ranges, c.val, got, c.want)
		}
	}
}

func TestTrimQualityWindow_Mean(t *testing.T) {
	// Quality scores (Phred+33): a 4-base low block then high. With a
	// window of 4 and a mean rule (lt 20), the first window mean is low so
	// we advance; once the window mean climbs above 20 we stop.
	seq := []byte("ACGTACGTAC")
	// scores: 5 5 5 5 40 40 40 40 40 40  ('&' = 5, 'I' = 40)
	qual := []byte("&&&&IIIIII")
	opts := FilterOptions{TrimQualL: 20, TrimQualWindow: 4, TrimQualStep: 1, TrimQualType: "mean", QualType: "sanger"}
	gotSeq, gotQual := trimQualityWindow(seq, qual, opts)
	// Window means: [0..3]=5 (<20 trim), [1..4]=(5+5+5+40)/4=13.75 (<20
	// trim), [2..5]=(5+5+40+40)/4=22.5 (>=20 stop). So trim 2 bases.
	if string(gotSeq) != "GTACGTAC" {
		t.Errorf("seq = %q, want %q", gotSeq, "GTACGTAC")
	}
	if string(gotQual) != "&&IIIIII" {
		t.Errorf("qual = %q, want %q", gotQual, "&&IIIIII")
	}
}

func TestTrimQualityWindow_RuleGT(t *testing.T) {
	// rule "gt" trims while score > threshold. With threshold 10 and a
	// per-base window, the leading high-quality bases are trimmed until a
	// base at/below 10 is found.
	seq := []byte("ACGTAC")
	qual := []byte("IIII&&") // 40 40 40 40 5 5
	opts := FilterOptions{TrimQualL: 10, TrimQualWindow: 1, TrimQualStep: 1, TrimQualRule: "gt", QualType: "sanger"}
	gotSeq, _ := trimQualityWindow(seq, qual, opts)
	if string(gotSeq) != "AC" {
		t.Errorf("seq = %q, want %q", gotSeq, "AC")
	}
}

func TestTrimToLen_ViaFilter(t *testing.T) {
	in := []byte(">s1\nACGTACGTACGT\n>s2\nACGT\n")
	var out bytes.Buffer
	opts := FilterOptions{TrimToLen: 8}
	if err := Filter(bytes.NewReader(in), &out, false, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	// s1 (12 bp) hard-trimmed to 8; s2 (4 bp) left untouched (length <= 8).
	want := ">s1\nACGTACGT\n>s2\nACGT\n"
	if out.String() != want {
		t.Errorf("trim_to_len output:\n%q\nwant:\n%q", out.String(), want)
	}
}

func TestRangeLen_ViaFilter(t *testing.T) {
	in := []byte(">s1\nACGTACGTACGT\n>s2\nACGT\n>s3\nACGTACGT\n")
	var out bytes.Buffer
	// Keep only reads 8..12 long: drops s2 (4 bp), keeps s1 (12) and s3 (8).
	opts := FilterOptions{RangeLen: "8-12"}
	if err := Filter(bytes.NewReader(in), &out, false, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	want := ">s1\nACGTACGTACGT\n>s3\nACGTACGT\n"
	if out.String() != want {
		t.Errorf("range_len output:\n%q\nwant:\n%q", out.String(), want)
	}
}

func TestRangeGC_ViaFilter(t *testing.T) {
	in := []byte(">s1\nACGT\n>s2\nGGGG\n>s3\nAAAA\n")
	var out bytes.Buffer
	// GC%: s1=50, s2=100, s3=0. Keep 40..60 -> only s1.
	opts := FilterOptions{RangeGC: "40-60"}
	if err := Filter(bytes.NewReader(in), &out, false, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	want := ">s1\nACGT\n"
	if out.String() != want {
		t.Errorf("range_gc output:\n%q\nwant:\n%q", out.String(), want)
	}
}

func TestRangeGC_IntegerTruncation(t *testing.T) {
	// A 3-base read with 1 GC base has GC% = 33.33 -> truncated to 33.
	// A range of 33-33 must therefore keep it while 34-34 rejects it.
	in := []byte(">s1\nACA\n")
	for _, tc := range []struct {
		rng  string
		keep bool
	}{
		{"33-33", true},
		{"34-34", false},
		{"0-32", false},
	} {
		var out bytes.Buffer
		if err := Filter(bytes.NewReader(in), &out, false, FilterOptions{RangeGC: tc.rng}); err != nil {
			t.Fatalf("Filter: %v", err)
		}
		got := out.Len() > 0
		if got != tc.keep {
			t.Errorf("range_gc %s: kept=%v, want %v (out=%q)", tc.rng, got, tc.keep, out.String())
		}
	}
}

func TestZeroLengthDropped(t *testing.T) {
	// A read whose sequence is empty after trimming must be dropped, not
	// emitted as an empty record (upstream's zero_length filter).
	if !shouldFilterSequence("", "", FilterOptions{}) {
		t.Error("empty sequence should be filtered out")
	}
	// And a non-empty read must survive the same call.
	if shouldFilterSequence("ACGT", "", FilterOptions{}) {
		t.Error("non-empty sequence should not be filtered by the zero-length rule")
	}
}
