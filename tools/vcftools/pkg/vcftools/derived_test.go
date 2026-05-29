package vcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
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

// TestParity_Derived_Freq verifies --derived --freq byte-for-byte against
// genuine VCFtools 0.1.18 (minimal %g float formatting, e.g. A:0.5 not
// A:0.500000). Allele identity and ordering are upstream-byte-for-byte.
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
// any REF/ALT are dropped from --counts output. sample.vcf has only
// two biallelic-with-AA sites:
//
//   - 20:1230237 REF=T ALT=. AA=T → AA matches REF, kept, no swap.
//   - 1230237 is technically degenerate ("biallelic" with ALT=".") but
//     the port emits it under plain --counts; under --derived it must
//     also be emitted (AA matches REF) with the leading column REF=T.
//
// Every other biallelic SNP in sample.vcf has no AA and must be
// dropped under --derived.
func TestDerived_DropsSitesWithoutMatchingAA(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Counts:  true,
		Derived: true,
	})
	lines := readFileLines(t, prefix+".frq.count")
	if len(lines) < 1 || !strings.HasPrefix(lines[0], "CHROM\tPOS\t") {
		t.Fatalf("missing header row: %q", lines)
	}
	// Expect exactly one data row: 20:1230237.
	var dataRows []string
	for _, ln := range lines[1:] {
		if strings.TrimSpace(ln) == "" {
			continue
		}
		dataRows = append(dataRows, ln)
	}
	if len(dataRows) != 1 {
		t.Fatalf("expected 1 derived-kept biallelic site on sample.vcf, got %d: %v", len(dataRows), dataRows)
	}
	if !strings.HasPrefix(dataRows[0], "20\t1230237\t") {
		t.Errorf("expected derived row at 20:1230237, got %q", dataRows[0])
	}
}
