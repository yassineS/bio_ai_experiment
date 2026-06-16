package vcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// TestParity_Derived_Counts verifies --derived --counts byte-for-byte
// against upstream. The fixture covers:
//
//   - 1:100 AA=A, REF=A   → no swap (A:3 C:3).
//   - 1:200 AA=T, REF=G   → swap   (T:3 G:3).
//   - 1:300 AA=N          → AA doesn't match any allele → dropped.
//   - 1:400 AA=a, REF=A   → upstream uppercases AA before match → no swap.
//   - 1:500 AA=.          → missing AA → dropped.
//   - 2:100 AA=?          → upstream sentinel for missing → dropped.
//
// Ported from variant_file_output.cpp:67-159 (output_frequency, the
// `derived` branch).
func TestParity_Derived_Counts(t *testing.T) {
	prefix := runVcftoolsParity(t, "derived_fixture.vcf", &Params{
		Counts:  true,
		Derived: true,
	})
	got := readFileBytes(t, prefix+".frq.count")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "derived.expected.frq.count"))
	if !bytes.Equal(got, want) {
		t.Errorf(".frq.count mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Derived_Freq verifies --derived --freq. The data lines now
// match upstream byte-for-byte, including float formatting: the port emits
// the same C++ `defaultfloat` precision-6 form upstream uses (e.g. "0.5",
// not "0.500000"). Allele identity and ordering are upstream-exact.
func TestParity_Derived_Freq(t *testing.T) {
	prefix := runVcftoolsParity(t, "derived_fixture.vcf", &Params{
		Freq:    true,
		Derived: true,
	})
	got := readFileBytes(t, prefix+".frq")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "derived.expected.frq"))
	if !bytes.Equal(got, want) {
		t.Errorf(".frq mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestDerived_NoFreqIsNoOp pins that --derived alone (without --freq or
// --counts) produces no .frq / .frq.count file: it's a modifier flag.
// Mirrors upstream behaviour (parameters.cpp:201 only sets a boolean;
// the reorder logic lives inside output_frequency, which only runs when
// --freq / --counts is set).
func TestDerived_NoFreqIsNoOp(t *testing.T) {
	tmp := t.TempDir()
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "derived_fixture.vcf"))
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer in.Close()
	prefix := filepath.Join(tmp, "out")
	if err := Run(in, &Params{OutPrefix: prefix, Derived: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if _, err := os.Stat(prefix + ".frq"); !os.IsNotExist(err) {
		t.Errorf("--derived alone produced .frq (should be no-op without --freq/--counts)")
	}
	if _, err := os.Stat(prefix + ".frq.count"); !os.IsNotExist(err) {
		t.Errorf("--derived alone produced .frq.count (should be no-op without --freq/--counts)")
	}
}

// TestDerived_DropsSitesWithoutMatchingAA pins that when --derived is
// supplied, sites where INFO/AA is absent, ".", "?", or does not match
// any allele are dropped, and the surviving rows (including multi-allelic
// ones) are emitted with the ancestral allele first — byte-for-byte
// against LIVE upstream. On sample.vcf the survivors are 20:1110696
// (AA=T, multi-allelic), 20:1230237 (AA=T, monomorphic), and 20:1234567
// (AA=G, multi-allelic).
func TestDerived_DropsSitesWithoutMatchingAA(t *testing.T) {
	vcf := fixtureVCF(t, "sample.vcf")
	up := readFileBytes(t, runUpstream(t, vcf, "--counts", "--derived")+".frq.count")
	got := readFileBytes(t, runGo(t, vcf, &Params{Counts: true, Derived: true})+".frq.count")
	if !bytes.Equal(got, up) {
		t.Errorf(".frq.count mismatch\nupstream:\n%s\ngot:\n%s", up, got)
	}
}
