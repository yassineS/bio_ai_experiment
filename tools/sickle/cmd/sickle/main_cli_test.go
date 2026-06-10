package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// buildSickleCLI compiles our sickle binary into a temp dir once per test and
// returns its path. It backs the POSIX-bundling CLI tests below.
func buildSickleCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "sickle")
	cmd := exec.Command("go", "build", "-o", bin, ".")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build sickle: %v\n%s", err, out)
	}
	return bin
}

// upstreamSickle returns the path to the upstream sickle binary built from the
// reference_code submodule, or "" when it has not been built. The submodule is
// rooted four levels up from tools/sickle/cmd/sickle.
func upstreamSickle() string {
	p := filepath.Join("..", "..", "..", "..", "reference_code", "sickle", "sickle")
	if _, err := os.Stat(p); err == nil {
		return p
	}
	return ""
}

// sickleSEFixture is a small single-end FASTQ whose second read has a low-quality
// 3' tail (so trimming actually changes the output) and a passing first read.
// Reads are 30 bp so they survive the default length threshold after trimming.
const sickleSEFixture = "@r1\n" +
	"ACGTACGTACGTACGTACGTACGTACGTAC\n" +
	"+\n" +
	"IIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
	"@r2\n" +
	"ACGTACGTACGTACGTACGTACGTACGTAC\n" +
	"+\n" +
	"IIIIIIIIIIIIIIIIIIIIIIII######\n"

// TestSickleSEBundledParsesAndEquivalent asserts that the bundled short-flag
// form "-xn" (no-fiveprime + trunc-n) plus value-concat "-q20" parses through
// cliflag.Parse and produces output identical to the fully-expanded canonical
// form "-x -n -q 20". This exercises POSIX getopt bundling end-to-end.
func TestSickleSEBundledParsesAndEquivalent(t *testing.T) {
	bin := buildSickleCLI(t)
	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(sickleSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	run := func(args ...string) []byte {
		full := append([]string{"se", "-f", in, "-t", "sanger", "--quiet"}, args...)
		cmd := exec.Command(bin, full...)
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("sickle se %v: %v\n%s", args, err, stderr.String())
		}
		return stdout.Bytes()
	}

	bundled := run("-xnq20") // -x -n -q 20, fully bundled + value-concat
	canonical := run("-x", "-n", "-q", "20")

	if len(bundled) == 0 {
		t.Fatal("empty output from bundled form")
	}
	if !bytes.Equal(bundled, canonical) {
		t.Errorf("bundled vs canonical differ:\nbundled:   %q\ncanonical: %q", bundled, canonical)
	}
}

// TestSickleSEBundledMatchesUpstream compares our bundled-flag output against
// the upstream sickle binary for the same logical command line. Upstream sickle
// uses getopt_long, so it accepts the same "-xnq 20" cluster; we feed it the
// expanded equivalent it understands and require byte-identical trimming.
func TestSickleSEBundledMatchesUpstream(t *testing.T) {
	up := upstreamSickle()
	if up == "" {
		t.Fatalf("upstream sickle binary not found; build it with `make` in reference_code/sickle")
	}
	bin := buildSickleCLI(t)

	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(sickleSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}
	outOurs := filepath.Join(t.TempDir(), "ours.fq")
	outUp := filepath.Join(t.TempDir(), "up.fq")

	// Ours: bundled form through cliflag.Parse.
	ours := exec.Command(bin, "se", "-f", in, "-t", "sanger", "--quiet", "-o", outOurs, "-xnq20")
	if out, err := ours.CombinedOutput(); err != nil {
		t.Fatalf("ours sickle se: %v\n%s", err, out)
	}
	// Upstream: equivalent bundled cluster (getopt accepts -xnq 20).
	upCmd := exec.Command(up, "se", "-f", in, "-t", "sanger", "--quiet", "-o", outUp, "-xnq", "20")
	if out, err := upCmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream sickle se: %v\n%s", err, out)
	}

	gotOurs, err := os.ReadFile(outOurs)
	if err != nil {
		t.Fatal(err)
	}
	gotUp, err := os.ReadFile(outUp)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(gotOurs, gotUp) {
		t.Errorf("our bundled output differs from upstream:\nours:     %q\nupstream: %q", gotOurs, gotUp)
	}
}

// TestSickleSEQuietShortAlias asserts the -z short alias for --quiet (added for
// upstream retro-compat) parses and silences the stats banner like upstream.
func TestSickleSEQuietShortAlias(t *testing.T) {
	bin := buildSickleCLI(t)
	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(sickleSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command(bin, "se", "-f", in, "-t", "sanger", "-z")
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("sickle se -z: %v\n%s", err, stderr.String())
	}
	if bytes.Contains(stderr.Bytes(), []byte("records")) {
		t.Errorf("-z did not suppress statistics; stderr=%q", stderr.String())
	}
}
