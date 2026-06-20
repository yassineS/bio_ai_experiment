package difffuzz

import (
	"bytes"
	"testing"
)

// seedVCF is a small valid VCF used to seed the mutation tests.
var seedVCF = []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n" +
	"chr1\t200\t.\tG\tC\t60\tPASS\tDP=20\n")

func TestMutatorDeterministic(t *testing.T) {
	// Same seed -> identical input sequence.
	a := NewMutator(FormatVCF, seedVCF, 42)
	b := NewMutator(FormatVCF, seedVCF, 42)
	for i := 0; i < 50; i++ {
		ia := a.Next(i)
		ib := b.Next(i)
		if ia.Origin != ib.Origin {
			t.Fatalf("iter %d origin mismatch: %q vs %q", i, ia.Origin, ib.Origin)
		}
		if !bytes.Equal(ia.Data, ib.Data) {
			t.Fatalf("iter %d data mismatch for same seed", i)
		}
	}
}

func TestMutatorDifferentSeeds(t *testing.T) {
	a := NewMutator(FormatVCF, seedVCF, 1)
	b := NewMutator(FormatVCF, seedVCF, 2)
	same := 0
	for i := 0; i < 50; i++ {
		if bytes.Equal(a.Next(i).Data, b.Next(i).Data) {
			same++
		}
	}
	// Structured/raw outputs depend on RNG; different seeds should diverge often.
	if same > 25 {
		t.Fatalf("different seeds produced too many identical inputs: %d/50", same)
	}
}

func TestMutatorAllOriginsExercised(t *testing.T) {
	m := NewMutator(FormatVCF, seedVCF, 7)
	seen := map[Origin]bool{}
	for i := 0; i < 8; i++ {
		seen[m.Next(i).Origin] = true
	}
	for _, o := range []Origin{OriginMutation, OriginStructured, OriginRaw} {
		if !seen[o] {
			t.Errorf("origin %q never produced in 8 iterations", o)
		}
	}
}

func TestRawBytesTargetHasNoStructured(t *testing.T) {
	m := NewMutator(FormatRawBytes, nil, 3)
	for i := 0; i < 12; i++ {
		in := m.Next(i)
		if in.Origin == OriginStructured {
			t.Fatalf("RawBytes target should never produce structured inputs, got one at %d", i)
		}
	}
}

func TestBoundaryValueReplacesNumeric(t *testing.T) {
	m := NewMutator(FormatBED, nil, 1)
	in := []byte("chr1\t100\t200\tname\t0\t+\n")
	out := m.boundaryValue(in)
	// A boundary value should land in one of the numeric columns.
	if bytes.Equal(out, in) {
		// Acceptable only if no numeric field existed, but this input has them.
		t.Fatalf("boundaryValue did not change a numeric field: %q", out)
	}
	if !bytes.HasPrefix(out, []byte("chr1\t")) {
		t.Fatalf("boundaryValue corrupted structure: %q", out)
	}
}

func TestDuplicateAndReorderRecords(t *testing.T) {
	m := NewMutator(FormatBED, nil, 9)
	in := []byte("a\nb\nc\n")
	dup := m.duplicateRecords(in)
	if bytes.Count(dup, []byte("\n")) <= bytes.Count(in, []byte("\n")) {
		t.Fatalf("duplicateRecords did not add a line: %q", dup)
	}
	reordered := m.reorderRecords([]byte("a\nb\nc\nd\n"))
	if len(reordered) != len("a\nb\nc\nd\n") {
		t.Fatalf("reorderRecords changed length: %q", reordered)
	}
}

func TestTruncateShrinks(t *testing.T) {
	m := NewMutator(FormatVCF, nil, 5)
	in := bytes.Repeat([]byte("x"), 100)
	out := m.truncate(in)
	if len(out) >= len(in) {
		t.Fatalf("truncate did not shrink: %d >= %d", len(out), len(in))
	}
}

func TestStructuredGeneratorsParseShape(t *testing.T) {
	cases := []struct {
		format Format
		prefix string
	}{
		{FormatVCF, "##fileformat"},
		{FormatSAM, "@HD"},
		{FormatFASTA, ">"},
		{FormatFASTQ, "@"},
	}
	for _, c := range cases {
		m := NewMutator(c.format, nil, 11)
		got := m.structured()
		if !bytes.HasPrefix(got, []byte(c.prefix)) {
			t.Errorf("structured %s: want prefix %q, got %q", c.format, c.prefix, head(got))
		}
	}
	// BED has no header; just ensure it has tab-separated columns.
	bedM := NewMutator(FormatBED, nil, 11)
	bed := bedM.structured()
	if !bytes.Contains(bed, []byte("\t")) {
		t.Errorf("structured BED has no tabs: %q", head(bed))
	}
}

func head(b []byte) string {
	if len(b) > 40 {
		return string(b[:40])
	}
	return string(b)
}

func TestSplitJoinRoundTrip(t *testing.T) {
	for _, in := range [][]byte{
		[]byte("a\nb\nc\n"),
		[]byte("a\nb\nc"),
		[]byte(""),
		[]byte("\n"),
		[]byte("single"),
	} {
		got := joinLines(splitKeepEmpty(in))
		if !bytes.Equal(got, in) {
			t.Errorf("round-trip failed: %q -> %q", in, got)
		}
	}
}
