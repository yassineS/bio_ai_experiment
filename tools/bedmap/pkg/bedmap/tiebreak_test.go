package bedmap

import (
	"bytes"
	"testing"
)

// TestUnitMapCollapseInputOrder pins that bedmap collapses an A record's
// overlapping B values in B-file (input) order on equal (chrom, start) keys —
// NOT by chromEnd. This guards the tie-break fix without invoking the binary.
func TestUnitMapCollapseInputOrder(t *testing.T) {
	a := "chr1\t0\t200\tA1\n"
	// Equal-start B records, ends out of order (100,50,75,60 -> input a,b,c,d).
	b := "chr1\t10\t100\ta\nchr1\t10\t50\tb\nchr1\t10\t75\tc\nchr1\t10\t60\td\n"

	for _, op := range []string{"collapse", "distinct"} {
		var out bytes.Buffer
		n, err := Map(bytes.NewReader([]byte(a)), bytes.NewReader([]byte(b)), &out,
			Options{Columns: []int{4}, Ops: []string{op}})
		if err != nil {
			t.Fatalf("Map(%s): %v", op, err)
		}
		if n != 1 {
			t.Fatalf("Map(%s): n=%d, want 1", op, n)
		}
		want := "chr1\t0\t200\tA1\ta,b,c,d\n"
		if out.String() != want {
			t.Fatalf("%s order.\nwant: %q\ngot:  %q", op, want, out.String())
		}
	}
}

// TestUnitMapCollapseInputOrderMixedStrand checks the order is preserved when a
// strand filter selects a subset of the equal-start B records.
func TestUnitMapCollapseInputOrderMixedStrand(t *testing.T) {
	a := "chr1\t0\t200\tA1\t0\t+\n"
	b := "chr1\t10\t100\ta\t1\t+\n" +
		"chr1\t10\t50\tb\t2\t-\n" +
		"chr1\t10\t75\tc\t3\t+\n" +
		"chr1\t10\t60\td\t4\t-\n"
	var out bytes.Buffer
	if _, err := Map(bytes.NewReader([]byte(a)), bytes.NewReader([]byte(b)), &out,
		Options{Columns: []int{4}, Ops: []string{"collapse"}, SameStrand: true}); err != nil {
		t.Fatalf("Map: %v", err)
	}
	// Only the '+' B records (a, c) survive, in input order.
	want := "chr1\t0\t200\tA1\t0\t+\ta,c\n"
	if out.String() != want {
		t.Fatalf("strand collapse order.\nwant: %q\ngot:  %q", want, out.String())
	}
}
