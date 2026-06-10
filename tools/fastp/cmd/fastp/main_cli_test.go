package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// buildFastpCLI compiles our fastp binary into a temp dir once per test and
// returns its path. It backs the POSIX-bundling CLI test below.
func buildFastpCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "fastp")
	cmd := exec.Command("go", "build", "-o", bin, ".")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build fastp: %v\n%s", err, out)
	}
	return bin
}

// fastpFixture is a single 32 bp read with a low-quality 3' tail so that the
// -5/-3 quality-cut switches produce an observable trim.
const fastpFixture = "@r1\n" +
	"ACGTACGTACGTACGTACGTACGTACGTACGT\n" +
	"+\n" +
	"IIIIIIIIIIIIIIIIIIIIIIIIIIII####\n"

// TestFastpBundledParsesAndEquivalent asserts that the bundled boolean cluster
// "-53" (cut-front + cut-tail) parses through cliflag.Parse and yields output
// byte-identical to the expanded "-5 -3" form. fastp's single-letter flags were
// deliberately reassigned relative to upstream (our -x is adapter3, upstream's
// is trim_poly_x), so this is an in-binary bundled==canonical equivalence
// assertion of the parsing rollout rather than an upstream byte-parity check.
func TestFastpBundledParsesAndEquivalent(t *testing.T) {
	bin := buildFastpCLI(t)
	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(fastpFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	run := func(cutArg string) []byte {
		out := filepath.Join(t.TempDir(), "out.fq")
		cmd := exec.Command(bin, "-i", in, "-o", out, cutArg)
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("fastp %s: %v\n%s", cutArg, err, stderr.String())
		}
		got, err := os.ReadFile(out)
		if err != nil {
			t.Fatal(err)
		}
		return got
	}

	bundled := run("-53")

	// Canonical: two separate tokens.
	outC := filepath.Join(t.TempDir(), "outc.fq")
	cmd := exec.Command(bin, "-i", in, "-o", outC, "-5", "-3")
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("fastp -5 -3: %v\n%s", err, stderr.String())
	}
	canonical, err := os.ReadFile(outC)
	if err != nil {
		t.Fatal(err)
	}

	if len(bundled) == 0 {
		t.Fatal("empty output from bundled -53 form")
	}
	if !bytes.Equal(bundled, canonical) {
		t.Errorf("bundled -53 vs canonical -5 -3 differ:\nbundled:   %q\ncanonical: %q", bundled, canonical)
	}
}
