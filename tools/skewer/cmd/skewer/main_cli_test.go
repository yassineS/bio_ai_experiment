package main

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// buildSkewerCLI compiles our skewer binary into a temp dir once per test and
// returns its path. It backs the POSIX-bundling CLI tests below.
func buildSkewerCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "skewer")
	cmd := exec.Command("go", "build", "-o", bin, ".")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build skewer: %v\n%s", err, out)
	}
	return bin
}

// skewerFixture is a single 32 bp read used to exercise the SE pipeline.
const skewerFixture = "@r1\n" +
	"ACGTACGTACGTACGTACGTACGTACGTACGT\n" +
	"+\n" +
	"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"

// TestSkewerBundledParsesAndEquivalent asserts that the bundled cluster
// "-al30" (auto-detect boolean + value-concat min-length 30) parses through
// cliflag.Parse and yields output byte-identical to the expanded "-a -l 30"
// form. Upstream skewer cannot be built offline with modern g++ (its old C++
// trips the libstdc++ comparator static_assert), so we validate against the
// canonical expansion within our own binary; the source-derived expectation is
// that POSIX getopt clustering ("-al30" == "-a -l 30") leaves trimming output
// unchanged. This is the in-binary equivalence tier.
func TestSkewerBundledParsesAndEquivalent(t *testing.T) {
	bin := buildSkewerCLI(t)
	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(skewerFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	run := func(args ...string) []byte {
		full := append([]string{"se", "-i", in, "--quiet"}, args...)
		cmd := exec.Command(bin, full...)
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("skewer se %v: %v\n%s", args, err, stderr.String())
		}
		return stdout.Bytes()
	}

	bundled := run("-al30")
	canonical := run("-a", "-l", "30")

	if len(bundled) == 0 {
		t.Fatal("empty output from bundled -al30 form")
	}
	if !bytes.Equal(bundled, canonical) {
		t.Errorf("bundled -al30 vs canonical -a -l 30 differ:\nbundled:   %q\ncanonical: %q", bundled, canonical)
	}
}

// TestSkewerCompressShortAlias asserts the -z/--compress retro-compat flag
// forces gzip output (matching upstream skewer's -z) and that the gzip stream
// decompresses to the same bytes as the uncompressed run.
func TestSkewerCompressShortAlias(t *testing.T) {
	bin := buildSkewerCLI(t)
	in := filepath.Join(t.TempDir(), "in.fq")
	if err := os.WriteFile(in, []byte(skewerFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	plainPath := filepath.Join(t.TempDir(), "plain.fq")
	gzPath := filepath.Join(t.TempDir(), "comp.fq")

	if out, err := exec.Command(bin, "se", "-i", in, "--quiet", "-o", plainPath).CombinedOutput(); err != nil {
		t.Fatalf("skewer se plain: %v\n%s", err, out)
	}
	if out, err := exec.Command(bin, "se", "-i", in, "--quiet", "-z", "-o", gzPath).CombinedOutput(); err != nil {
		t.Fatalf("skewer se -z: %v\n%s", err, out)
	}

	plain, err := os.ReadFile(plainPath)
	if err != nil {
		t.Fatal(err)
	}
	gzBytes, err := os.ReadFile(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	gr, err := gzip.NewReader(bytes.NewReader(gzBytes))
	if err != nil {
		t.Fatalf("-z output is not a valid gzip stream: %v", err)
	}
	decoded, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("decompress -z output: %v", err)
	}
	if !bytes.Equal(decoded, plain) {
		t.Errorf("-z decompressed output differs from plain run:\ngzip:  %q\nplain: %q", decoded, plain)
	}
}
