package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// Live byte-for-byte parity tests for the DEFAULT (dynamic) sliding window.
//
// Upstream sickle (reference_code/sickle/src/sliding.c) has no -w flag: it
// always computes window_size = int(0.1 * read_length) per read, falling back
// to the full read length when that rounds to 0. Our Go port reproduces this
// when -w is unset (or <= 0). These tests run BOTH binaries with NO -w flag on
// reads of *varying* lengths — so 0.1*length (and therefore the trimmed result)
// differs per read — and require byte-identical output.
//
// They build the upstream binary's presence via upstreamSickle() (built by
// `make` in reference_code/sickle); the harness t.Fatalf's if it is absent so
// the parity assertion is never silently skipped.

// runSickle runs the given sickle binary's `se`/`pe` subcommand, writing the
// trimmed records to files, and returns their contents. stderr stats are
// silenced with --quiet so only the trimmed FASTQ is compared.
func runSickleBinarySE(t *testing.T, bin, inPath string, extraArgs ...string) []byte {
	t.Helper()
	outPath := filepath.Join(t.TempDir(), "out.fq")
	args := append([]string{"se", "-f", inPath, "-t", "sanger", "--quiet", "-o", outPath}, extraArgs...)
	cmd := exec.Command(bin, args...)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("%s se %v: %v\n%s", bin, args, err, out)
	}
	data, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatal(err)
	}
	return data
}

// makeQual builds a quality string of length n. The first hi bases are at the
// high score (Phred 'I' = Q40); the remainder alternates between just-above and
// just-below the q20 threshold so the trimmed length depends on the window
// width (and therefore on 0.1*length).
func makeQual(n, hi int) string {
	var b strings.Builder
	for i := 0; i < n; i++ {
		switch {
		case i < hi:
			b.WriteByte('I') // Q40
		case i%2 == 0:
			b.WriteByte('5') // Q20 (== threshold, passes)
		default:
			b.WriteByte('0') // Q15 (< threshold)
		}
	}
	return b.String()
}

func makeSeq(n int) string {
	const bases = "ACGT"
	var b strings.Builder
	for i := 0; i < n; i++ {
		b.WriteByte(bases[i%4])
	}
	return b.String()
}

// varyingLengthFASTQ returns a FASTQ with reads of several different lengths so
// that int(0.1*length) — the dynamic default window — differs across reads:
//
//	len 25  -> window 2
//	len 40  -> window 4
//	len 60  -> window 6
//	len 105 -> window 10
//	len 8   -> window int(0.8)=0 -> fallback to full length (8)
//	len 33  -> window 3
//
// Each read mixes a high-quality 5' region with an oscillating near-threshold
// 3' tail so the sliding window's averaging actually moves the 3' cut.
func varyingLengthFASTQ() string {
	type rec struct {
		name   string
		length int
		hi     int
	}
	recs := []rec{
		{"r_len25", 25, 12},
		{"r_len40", 40, 18},
		{"r_len60", 60, 25},
		{"r_len105", 105, 40},
		{"r_len8", 8, 6},
		{"r_len33", 33, 15},
	}
	var b strings.Builder
	for _, r := range recs {
		b.WriteString("@" + r.name + "\n")
		b.WriteString(makeSeq(r.length) + "\n")
		b.WriteString("+\n")
		b.WriteString(makeQual(r.length, r.hi) + "\n")
	}
	return b.String()
}

// TestParityDynamicWindowSE_VaryingLengths drives upstream sickle and the Go
// port with NO -w flag over reads of varying lengths and requires byte-for-byte
// identical SE output. This exercises the per-read int(0.1*length) default
// window directly.
func TestParityDynamicWindowSE_VaryingLengths(t *testing.T) {
	up := upstreamSickle()
	if up == "" {
		t.Fatalf("upstream sickle binary not found; build it with `make` in reference_code/sickle")
	}
	ours := buildSickleCLI(t)

	in := filepath.Join(t.TempDir(), "varying.fq")
	if err := os.WriteFile(in, []byte(varyingLengthFASTQ()), 0o644); err != nil {
		t.Fatal(err)
	}

	// Default window (no -w), q20 l10.
	gotUp := runSickleBinarySE(t, up, in, "-q", "20", "-l", "10")
	gotOurs := runSickleBinarySE(t, ours, in, "-q", "20", "-l", "10")
	if !bytes.Equal(gotUp, gotOurs) {
		t.Errorf("default-window SE output differs from upstream.\nupstream:\n%s\nours:\n%s", gotUp, gotOurs)
	}

	// A second threshold profile (q25 l5) to vary the cut points further while
	// still using the dynamic per-read window.
	gotUp2 := runSickleBinarySE(t, up, in, "-q", "25", "-l", "5")
	gotOurs2 := runSickleBinarySE(t, ours, in, "-q", "25", "-l", "5")
	if !bytes.Equal(gotUp2, gotOurs2) {
		t.Errorf("default-window SE (q25 l5) differs from upstream.\nupstream:\n%s\nours:\n%s", gotUp2, gotOurs2)
	}
}

// TestParityDynamicWindowPE_VaryingLengths is the paired-end counterpart: both
// mates use the per-read dynamic window (mates can have different lengths, so
// different windows), and the R1/R2/singleton outputs must all byte-match
// upstream with no -w flag.
func TestParityDynamicWindowPE_VaryingLengths(t *testing.T) {
	up := upstreamSickle()
	if up == "" {
		t.Fatalf("upstream sickle binary not found; build it with `make` in reference_code/sickle")
	}
	ours := buildSickleCLI(t)

	// R1 and R2 deliberately use different per-read lengths so the dynamic
	// window differs between mates of the same pair.
	r1 := "@p1\n" + makeSeq(60) + "\n+\n" + makeQual(60, 25) + "\n" +
		"@p2\n" + makeSeq(40) + "\n+\n" + makeQual(40, 5) + "\n" +
		"@p3\n" + makeSeq(105) + "\n+\n" + makeQual(105, 50) + "\n"
	r2 := "@p1\n" + makeSeq(33) + "\n+\n" + makeQual(33, 15) + "\n" +
		"@p2\n" + makeSeq(25) + "\n+\n" + makeQual(25, 12) + "\n" +
		"@p3\n" + makeSeq(80) + "\n+\n" + makeQual(80, 30) + "\n"

	dir := t.TempDir()
	r1Path := filepath.Join(dir, "r1.fq")
	r2Path := filepath.Join(dir, "r2.fq")
	if err := os.WriteFile(r1Path, []byte(r1), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(r2Path, []byte(r2), 0o644); err != nil {
		t.Fatal(err)
	}

	runPE := func(bin string) (o1, o2, os3 []byte) {
		d := t.TempDir()
		out1 := filepath.Join(d, "o1.fq")
		out2 := filepath.Join(d, "o2.fq")
		outS := filepath.Join(d, "s.fq")
		cmd := exec.Command(bin, "pe",
			"-f", r1Path, "-r", r2Path, "-t", "sanger", "--quiet",
			"-o", out1, "-p", out2, "-s", outS, "-q", "20", "-l", "10")
		if out, err := cmd.CombinedOutput(); err != nil {
			t.Fatalf("%s pe: %v\n%s", bin, err, out)
		}
		read := func(p string) []byte {
			b, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			return b
		}
		return read(out1), read(out2), read(outS)
	}

	u1, u2, us := runPE(up)
	g1, g2, gs := runPE(ours)
	if !bytes.Equal(u1, g1) {
		t.Errorf("default-window PE R1 differs.\nupstream:\n%s\nours:\n%s", u1, g1)
	}
	if !bytes.Equal(u2, g2) {
		t.Errorf("default-window PE R2 differs.\nupstream:\n%s\nours:\n%s", u2, g2)
	}
	if !bytes.Equal(us, gs) {
		t.Errorf("default-window PE singletons differ.\nupstream:\n%s\nours:\n%s", us, gs)
	}
}
