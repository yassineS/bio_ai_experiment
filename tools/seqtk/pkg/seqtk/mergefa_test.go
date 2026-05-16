package seqtk

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// TestMergefa_DefaultMode covers the OR-merge default: A+G -> R,
// C+G -> S, T+T -> T (homozygous same), etc. The case rule (uppercase
// only when both inputs are uppercase) is exercised by the first base
// of the second record (lowercase 'a' on both sides => 'a').
//
// Note: only the leading 'a' is lowercase on both inputs; positions
// 1..3 are uppercase on both inputs, so upstream and this port both
// keep them uppercase. We verified the expected output by piping the
// same fixtures through `reference_code/seqtk/seqtk mergefa`.
func TestMergefa_DefaultMode(t *testing.T) {
	a := ">x\nACGT\n>y\naCGT\n"
	b := ">x\nAGGT\n>y\naCGT\n"
	want := ">x\nASGT\n>y\naCGT\n"
	var out, warn bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, &warn, MergefaOptions{}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("want %q, got %q", want, got)
	}
	// Counter line: chr1 has 3 same (A,G,T) + 1 diff (C vs G);
	// chr2 has 3 same (C,G,T) — the leading 'a' isn't uppercase
	// and is skipped by the bucket logic.
	if !strings.Contains(warn.String(), "(same,diff,hom-het,het-hom,het-het)=(6,1,0,0,0)") {
		t.Errorf("expected counter line in warn: %q", warn.String())
	}
}

// TestMergefa_IntersectMode covers `-i`: A+G -> empty -> 'x' (lowercase
// X), C+G -> empty -> 'x', T+T -> T (intersect == T).
func TestMergefa_IntersectMode(t *testing.T) {
	a := ">x\nACGT\n"
	b := ">x\nAGGT\n"
	want := ">x\nAxGT\n"
	var out bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, io.Discard, MergefaOptions{Intersect: true}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("want %q, got %q", want, got)
	}
}

// TestMergefa_MaskMode covers `-m`: like intersect, but additionally
// lowercases positions where either side is N.
func TestMergefa_MaskMode(t *testing.T) {
	a := ">x\nACNT\n"
	b := ">x\nACGT\n"
	want := ">x\nACgT\n"
	var out bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, io.Discard, MergefaOptions{Mask: true}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("want %q, got %q", want, got)
	}
}

// TestMergefa_HaploidMode covers `-h`: heterozygous merges (b0>1 or
// b1>1, i.e. either input is itself an IUPAC ambiguity code) become
// lowercase. With pure ACGT inputs, b0=b1=1 always so haploid is a
// no-op — the test feeds an R (=A/G) in one input to trigger it.
func TestMergefa_HaploidMode(t *testing.T) {
	a := ">x\nARGT\n"
	b := ">x\nAAGT\n"
	want := ">x\nArGT\n"
	var out bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, io.Discard, MergefaOptions{Haploid: true}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("want %q, got %q", want, got)
	}
}

// TestMergefa_IMConflictRejected verifies that combining -i and -m
// returns an error, matching upstream's early-exit check.
func TestMergefa_IMConflictRejected(t *testing.T) {
	a := ">x\nACGT\n"
	b := ">x\nACGT\n"
	var out bytes.Buffer
	err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, io.Discard, MergefaOptions{Intersect: true, Mask: true})
	if err == nil {
		t.Fatal("expected error for -i and -m combined")
	}
}

// TestMergefa_QualityLowering verifies that low-quality FASTQ bases are
// lowercased before merging — both sides independently per upstream's
// `seq->qual.l` guard.
func TestMergefa_QualityLowering(t *testing.T) {
	// '!' is Phred 0, 'I' is Phred 40. With -q 20, both bases of
	// chr2 become lowercase before merging; their case is preserved
	// across the OR-merge.
	a := "@chr2\nAAAA\n+\n!!!!\n"
	b := "@chr2\nGGGG\n+\nIIII\n"
	want := ">chr2\nrrrr\n"
	var out bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, io.Discard, MergefaOptions{Quality: 20}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("want %q, got %q", want, got)
	}
}

// TestMergefa_NameMismatchWarns confirms the per-record name-mismatch
// stderr warning fires.
func TestMergefa_NameMismatchWarns(t *testing.T) {
	a := ">x\nACGT\n"
	b := ">y\nACGT\n"
	var out, warn bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, &warn, MergefaOptions{}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if !strings.Contains(warn.String(), "Different sequence names") {
		t.Errorf("expected name-mismatch warning: %q", warn.String())
	}
}

// TestMergefa_LengthMismatchWarns confirms the per-record length
// mismatch stderr warning fires and the shorter length is used.
func TestMergefa_LengthMismatchWarns(t *testing.T) {
	a := ">x\nACGTACGT\n"
	b := ">x\nACGT\n"
	wantOut := ">x\nACGT\n"
	var out, warn bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, &warn, MergefaOptions{}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	if got := out.String(); got != wantOut {
		t.Errorf("want %q, got %q", wantOut, got)
	}
	if !strings.Contains(warn.String(), "Unequal sequence length") {
		t.Errorf("expected length-mismatch warning: %q", warn.String())
	}
}

// TestMergefa_RandHetIsDeterministic verifies that the -r path produces
// reproducible output for a given seed (project policy:
// reproducibility within our tool; byte-parity with upstream lrand48
// is explicitly NOT a goal — see docs/PARITY_ROADMAP.md#rng-policy).
func TestMergefa_RandHetIsDeterministic(t *testing.T) {
	a := ">x\nARGT\n>y\nMYRK\n"
	b := ">x\nAAGT\n>y\nAAAA\n"
	var out1, out2 bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out1, io.Discard, MergefaOptions{RandHet: true, Seed: 42}); err != nil {
		t.Fatalf("mergefaImpl 1: %v", err)
	}
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out2, io.Discard, MergefaOptions{RandHet: true, Seed: 42}); err != nil {
		t.Fatalf("mergefaImpl 2: %v", err)
	}
	if !bytes.Equal(out1.Bytes(), out2.Bytes()) {
		t.Errorf("seeded -r mode is not deterministic:\nout1: %q\nout2: %q", out1.String(), out2.String())
	}
}

// TestMergefa_LineWrapAt60 confirms output is wrapped at 60 bases per
// line, matching upstream.
func TestMergefa_LineWrapAt60(t *testing.T) {
	a := ">x\n" + strings.Repeat("A", 125) + "\n"
	b := ">x\n" + strings.Repeat("G", 125) + "\n"
	var out bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, io.Discard, MergefaOptions{}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	want := ">x\n" + strings.Repeat("R", 60) + "\n" + strings.Repeat("R", 60) + "\n" + strings.Repeat("R", 5) + "\n"
	if got := out.String(); got != want {
		t.Errorf("wrap mismatch:\nwant %q\ngot  %q", want, got)
	}
}

// TestMergefa_CounterAccounting verifies the (same,diff,hom-het,
// het-hom,het-het) counts in the stderr summary line, exercising
// each bucket once.
func TestMergefa_CounterAccounting(t *testing.T) {
	// pos 0: A,A -> same (bucket 0)
	// pos 1: A,G -> diff (bucket 1)
	// pos 2: A,R -> hom-het (bucket 2; b0=1, b1=2 since R=A/G)
	// pos 3: R,A -> het-hom (bucket 3; b0=2, b1=1)
	// pos 4: R,Y -> het-het (bucket 4; b0=2, b1=2)
	a := ">x\nAAARR\n"
	b := ">x\nAGRAY\n"
	var out, warn bytes.Buffer
	if err := mergefaImpl(strings.NewReader(a), strings.NewReader(b), &out, &warn, MergefaOptions{}); err != nil {
		t.Fatalf("mergefaImpl: %v", err)
	}
	wantSummary := "(same,diff,hom-het,het-hom,het-het)=(1,1,1,1,1)"
	if !strings.Contains(warn.String(), wantSummary) {
		t.Errorf("counter mismatch: want %q in %q", wantSummary, warn.String())
	}
}
