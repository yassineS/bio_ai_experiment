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

// upstreamSEFixture is a couple of reads with a clear 3' TruSeq adapter so the
// upstream-style CLI exercises real trimming (not just a pass-through).
const upstreamSEFixture = "@s1\n" +
	"ACGTACGTACGTACGTACGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n" +
	"+\n" +
	"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
	"@s2\n" +
	"GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG\n" +
	"+\n" +
	"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"

// TestSkewerUpstreamSE_NamingAndContent asserts the upstream positional CLI
// (`skewer [options] <reads.fastq> -o <base>`) accepts the drop-in form, exits
// 0, writes <base>-trimmed.fastq, and produces output byte-identical to the
// equivalent `skewer se` subcommand. This is the in-binary equivalence tier:
// upstream skewer's C++ does not always build offline, so we anchor on the
// subcommand whose trimming logic is separately validated byte-exact against
// upstream in pkg/skewer/parity_test.go.
func TestSkewerUpstreamSE_NamingAndContent(t *testing.T) {
	bin := buildSkewerCLI(t)
	dir := t.TempDir()
	in := filepath.Join(dir, "reads.fq")
	if err := os.WriteFile(in, []byte(upstreamSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	base := filepath.Join(dir, "out")
	x := "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	if out, err := exec.Command(bin, "-m", "tail", "-x", x, "-l", "8", "-o", base, "--quiet", in).CombinedOutput(); err != nil {
		t.Fatalf("skewer upstream SE: %v\n%s", err, out)
	}

	got, err := os.ReadFile(base + "-trimmed.fastq")
	if err != nil {
		t.Fatalf("expected upstream SE output %s-trimmed.fastq: %v", base, err)
	}

	subOut := filepath.Join(dir, "sub.fq")
	if out, err := exec.Command(bin, "se", "-i", in, "-x", x, "-l", "8", "-o", subOut, "--quiet").CombinedOutput(); err != nil {
		t.Fatalf("skewer se subcommand: %v\n%s", err, out)
	}
	want, err := os.ReadFile(subOut)
	if err != nil {
		t.Fatal(err)
	}

	if !bytes.Equal(got, want) {
		t.Errorf("upstream-CLI SE output differs from se subcommand:\nupstream:\n%s\nsubcommand:\n%s", got, want)
	}
}

// TestSkewerUpstreamPE_NamingAndContent asserts the upstream positional CLI
// for paired-end (`skewer -m pe <r1> <r2> -o <base>`) accepts the drop-in form,
// exits 0, and writes <base>-trimmed-pair1.fastq / <base>-trimmed-pair2.fastq
// with content byte-identical to the equivalent `skewer pe` subcommand (with
// the upstream default mate adapters supplied explicitly).
func TestSkewerUpstreamPE_NamingAndContent(t *testing.T) {
	bin := buildSkewerCLI(t)
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	if err := os.WriteFile(r1, []byte(upstreamSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(r2, []byte(upstreamSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	x := "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	y := "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"

	base := filepath.Join(dir, "out")
	if out, err := exec.Command(bin, "-m", "pe", "-l", "8", "-o", base, "--quiet", r1, r2).CombinedOutput(); err != nil {
		t.Fatalf("skewer upstream PE: %v\n%s", err, out)
	}

	gotP1, err := os.ReadFile(base + "-trimmed-pair1.fastq")
	if err != nil {
		t.Fatalf("expected %s-trimmed-pair1.fastq: %v", base, err)
	}
	gotP2, err := os.ReadFile(base + "-trimmed-pair2.fastq")
	if err != nil {
		t.Fatalf("expected %s-trimmed-pair2.fastq: %v", base, err)
	}

	subP1 := filepath.Join(dir, "sub1.fq")
	subP2 := filepath.Join(dir, "sub2.fq")
	if out, err := exec.Command(bin, "pe", "-i", r1, "-j", r2, "-o", subP1, "-p", subP2, "-x", x, "-y", y, "-l", "8", "--quiet").CombinedOutput(); err != nil {
		t.Fatalf("skewer pe subcommand: %v\n%s", err, out)
	}
	wantP1, err := os.ReadFile(subP1)
	if err != nil {
		t.Fatal(err)
	}
	wantP2, err := os.ReadFile(subP2)
	if err != nil {
		t.Fatal(err)
	}

	if !bytes.Equal(gotP1, wantP1) {
		t.Errorf("upstream-CLI PE pair1 differs from pe subcommand:\nupstream:\n%s\nsubcommand:\n%s", gotP1, wantP1)
	}
	if !bytes.Equal(gotP2, wantP2) {
		t.Errorf("upstream-CLI PE pair2 differs from pe subcommand:\nupstream:\n%s\nsubcommand:\n%s", gotP2, wantP2)
	}
}

// TestSkewerUpstreamDefaultPrefix asserts that omitting -o derives the output
// base from the input filename with the last extension stripped, matching
// upstream's `<reads>-trimmed.fastq` default-naming behaviour.
func TestSkewerUpstreamDefaultPrefix(t *testing.T) {
	bin := buildSkewerCLI(t)
	dir := t.TempDir()
	in := filepath.Join(dir, "sample.fq")
	if err := os.WriteFile(in, []byte(upstreamSEFixture), 0o644); err != nil {
		t.Fatal(err)
	}

	if out, err := exec.Command(bin, "-m", "tail", "-l", "8", "--quiet", in).CombinedOutput(); err != nil {
		t.Fatalf("skewer upstream default-prefix: %v\n%s", err, out)
	}
	if _, err := os.Stat(filepath.Join(dir, "sample-trimmed.fastq")); err != nil {
		t.Fatalf("expected default-prefix output sample-trimmed.fastq: %v", err)
	}
}
