package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// buildSeqtk compiles the seqtk binary into a temp dir once per test and
// returns its path. It is used by the CLI-level flag tests below.
func buildSeqtk(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "seqtk")
	cmd := exec.Command("go", "build", "-o", bin, ".")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build seqtk: %v\n%s", err, out)
	}
	return bin
}

// TestRandbaseSeedShortAndLong asserts that the randbase command honours both
// the POSIX short form (-s) and the GNU long form (--seed) of the seed flag,
// and that a fixed seed produces deterministic, identical output across both
// forms. This guards the cliflag.Int64Var wiring of the seed flag.
func TestRandbaseSeedShortAndLong(t *testing.T) {
	bin := buildSeqtk(t)

	const fasta = ">s\nRYSWKM\n"
	in := filepath.Join(t.TempDir(), "ambig.fa")
	if err := os.WriteFile(in, []byte(fasta), 0o644); err != nil {
		t.Fatal(err)
	}

	run := func(seedFlag string) []byte {
		cmd := exec.Command(bin, "randbase", seedFlag, "42", in)
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("seqtk randbase %s 42: %v\n%s", seedFlag, err, stderr.String())
		}
		return stdout.Bytes()
	}

	short := run("-s")
	long := run("--seed")

	if len(short) == 0 {
		t.Fatal("empty output from -s form")
	}
	if !bytes.Equal(short, long) {
		t.Errorf("-s and --seed produced different output:\n-s:     %q\n--seed: %q", short, long)
	}
}

// TestSeqBundledParsesAndEquivalent asserts that a bundled short-flag cluster
// mixing boolean switches and a value-concat ("-r6l30": -r -6 -l 30) parses
// through cliflag.Parse and yields output byte-identical to the fully expanded
// canonical form. Upstream seqtk's `seq` redefines -l/-L semantics relative to
// our port, so this is an in-binary bundled==canonical equivalence assertion
// (the parsing-rollout invariant) rather than an upstream byte-parity check.
func TestSeqBundledParsesAndEquivalent(t *testing.T) {
	bin := buildSeqtk(t)

	// Three FASTQ records of differing length so the -l 30 length filter is
	// observable; quality lines are present so -6 has a defined meaning.
	const fq = "@a\n" +
		"ACGTACGTACGTACGTACGTACGTACGTACGT\n" + // 32 bp, kept
		"+\n" +
		"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
		"@b\n" +
		"ACGTACGT\n" + // 8 bp, dropped by -l 30
		"+\n" +
		"IIIIIIII\n"
	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(fq), 0o644); err != nil {
		t.Fatal(err)
	}

	run := func(args ...string) []byte {
		full := append([]string{"seq"}, append(args, in)...)
		cmd := exec.Command(bin, full...)
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("seqtk seq %v: %v\n%s", args, err, stderr.String())
		}
		return stdout.Bytes()
	}

	bundled := run("-r6l30")
	canonical := run("-r", "-6", "-l", "30")

	if len(bundled) == 0 {
		t.Fatal("empty output from bundled form")
	}
	if !bytes.Equal(bundled, canonical) {
		t.Errorf("bundled -r6l30 vs canonical differ:\nbundled:   %q\ncanonical: %q", bundled, canonical)
	}
}
