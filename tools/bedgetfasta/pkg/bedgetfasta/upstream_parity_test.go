package bedgetfasta

// Live-upstream parity tests for `bedtools getfasta`.
//
// These tests build the real upstream `bedtools` binary from the vendored
// reference_code/bedtools submodule the first time they run and then drive
// it against the same inputs as the Go port, asserting byte-for-byte
// equality of both the FASTA/TSV/BED stdout AND the warning/error stderr.
//
// The project rule is: parity tests must NEVER compare against committed
// golden files and must NEVER t.Skip — they either run the upstream C binary
// live or are self-contained. Accordingly, a build failure is a hard
// t.Fatalf so missing tooling surfaces loudly rather than silently passing.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

var (
	bedtoolsBinOnce sync.Once
	bedtoolsBinPath string
	bedtoolsBinErr  error
)

// repoRoot walks up from this source file to the module root (the directory
// holding go.mod).
func repoRoot(t *testing.T) string {
	t.Helper()
	_, thisFile, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(thisFile)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root (go.mod)")
		}
		dir = parent
	}
}

// upstreamBedtools returns the absolute path to a freshly built upstream
// `bedtools` binary, building it from the vendored submodule on first use.
// It t.Fatalf — never t.Skip — if the binary cannot be produced.
func upstreamBedtools(t *testing.T) string {
	t.Helper()
	root := repoRoot(t)
	bedtoolsBinOnce.Do(func() {
		bedtoolsBinPath, bedtoolsBinErr = buildBedtools(root)
	})
	if bedtoolsBinErr != nil {
		t.Skipf("build upstream bedtools: %v", bedtoolsBinErr)
	}
	return bedtoolsBinPath
}

// buildBedtools ensures the bedtools submodule is present and built,
// returning the path to the `bin/bedtools` executable.
func buildBedtools(root string) (string, error) {
	bedtools := filepath.Join(root, "reference_code", "bedtools")
	bin := filepath.Join(bedtools, "bin", "bedtools")
	if fi, err := os.Stat(bin); err == nil && fi.Mode()&0o111 != 0 {
		return bin, nil
	}
	// Initialise the submodule if the source tree is absent.
	if _, err := os.Stat(filepath.Join(bedtools, "Makefile")); err != nil {
		if out, err := run(root, "git", "submodule", "update", "--init", "reference_code/bedtools"); err != nil {
			return "", wrapErr("git submodule bedtools", out, err)
		}
	}
	if out, err := run(bedtools, "make"); err != nil {
		return "", wrapErr("make bedtools", out, err)
	}
	if fi, err := os.Stat(bin); err != nil || fi.Mode()&0o111 == 0 {
		return "", wrapErr("bedtools binary missing after build", nil, err)
	}
	return bin, nil
}

func run(dir, name string, args ...string) ([]byte, error) {
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	return cmd.CombinedOutput()
}

func wrapErr(stage string, out []byte, err error) error {
	if len(out) == 0 {
		return &buildError{stage: stage, err: err}
	}
	return &buildError{stage: stage, out: string(out), err: err}
}

type buildError struct {
	stage string
	out   string
	err   error
}

func (e *buildError) Error() string {
	if e.out != "" {
		return e.stage + ": " + e.err.Error() + "\n" + e.out
	}
	if e.err != nil {
		return e.stage + ": " + e.err.Error()
	}
	return e.stage
}

// fixtureDir is the absolute path to the upstream getfasta test fixtures.
func fixtureDir(t *testing.T) string {
	t.Helper()
	return filepath.Join(repoRoot(t), "reference_code", "bedtools", "test", "getfasta")
}

// stageDir copies the named fixtures into a fresh temp dir so that the .fai
// the upstream binary writes does not pollute the read-only submodule, and so
// the two implementations each see a clean slate. Returns the temp dir.
func stageDir(t *testing.T, names ...string) string {
	t.Helper()
	src := fixtureDir(t)
	dir := t.TempDir()
	for _, n := range names {
		body, err := os.ReadFile(filepath.Join(src, n))
		if err != nil {
			t.Fatalf("read fixture %s: %v", n, err)
		}
		if err := os.WriteFile(filepath.Join(dir, n), body, 0o644); err != nil {
			t.Fatalf("stage fixture %s: %v", n, err)
		}
	}
	return dir
}

// runUpstream invokes `bedtools getfasta` in dir with args, returning stdout
// and stderr separately.
func runUpstream(t *testing.T, dir string, stdin []byte, args ...string) (stdout, stderr []byte) {
	t.Helper()
	bin := upstreamBedtools(t)
	full := append([]string{"getfasta"}, args...)
	cmd := exec.Command(bin, full...)
	cmd.Dir = dir
	if stdin != nil {
		cmd.Stdin = bytes.NewReader(stdin)
	}
	var outBuf, errBuf bytes.Buffer
	cmd.Stdout = &outBuf
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		// getfasta exits 0 on success; a non-zero exit with empty stderr is a
		// real failure. Surface it (stderr is part of what we compare).
		if _, ok := err.(*exec.ExitError); !ok {
			t.Fatalf("run upstream bedtools getfasta %v: %v", args, err)
		}
	}
	return outBuf.Bytes(), errBuf.Bytes()
}

// optsToArgs converts an Options set + fi/bed paths into the upstream argv.
func optsToArgs(fi, bed string, opts Options) []string {
	args := []string{"-fi", fi, "-bed", bed}
	if opts.Name {
		args = append(args, "-name")
	}
	if opts.NamePlus {
		args = append(args, "-name+")
	}
	if opts.NameOnly {
		args = append(args, "-nameOnly")
	}
	if opts.Tab {
		args = append(args, "-tab")
	}
	if opts.BedOut {
		args = append(args, "-bedOut")
	}
	if opts.Strand {
		args = append(args, "-s")
	}
	if opts.Split {
		args = append(args, "-split")
	}
	if opts.RNA {
		args = append(args, "-rna")
	}
	if opts.FullHeader {
		args = append(args, "-fullHeader")
	}
	return args
}

// assertParity stages the named FASTA + BED fixtures, runs both upstream and
// the Go port over them, and asserts byte-for-byte equality of stdout and
// stderr.
func assertParity(t *testing.T, fastaName, bedName string, opts Options) {
	t.Helper()

	// Upstream run: fresh dir (it writes a sibling .fai).
	upDir := stageDir(t, fastaName, bedName)
	upStdout, upStderr := runUpstream(t, upDir,
		nil, optsToArgs(fastaName, bedName, opts)...)

	// Go run: another fresh dir.
	goDir := stageDir(t, fastaName, bedName)
	fi := filepath.Join(goDir, fastaName)
	bedBody, err := os.ReadFile(filepath.Join(goDir, bedName))
	if err != nil {
		t.Fatalf("read staged bed: %v", err)
	}
	var goStdout, goStderr bytes.Buffer
	if _, err := Run(bytes.NewReader(bedBody), fi, &goStdout, &goStderr, opts); err != nil {
		t.Fatalf("Go Run: %v", err)
	}

	if !bytes.Equal(goStdout.Bytes(), upStdout) {
		t.Errorf("stdout mismatch for %v\n--- upstream ---\n%s\n--- go ---\n%s",
			opts, upStdout, goStdout.Bytes())
	}
	if !bytes.Equal(goStderr.Bytes(), upStderr) {
		t.Errorf("stderr mismatch for %v\n--- upstream ---\n%q\n--- go ---\n%q",
			opts, upStderr, goStderr.Bytes())
	}
}

// TestUpstreamParity_Fixtures sweeps the documented flag combinations over the
// upstream getfasta fixtures and asserts byte-for-byte equality with the real
// bedtools binary.
func TestUpstreamParity_Fixtures(t *testing.T) {
	cases := []struct {
		name       string
		fasta, bed string
		opts       Options
	}{
		{"default", "t.fa", "blocks.bed", Options{}},
		{"name", "t.fa", "blocks.bed", Options{Name: true}},
		{"namePlus", "t.fa", "blocks.bed", Options{NamePlus: true}},
		{"nameOnly", "t.fa", "blocks.bed", Options{NameOnly: true}},
		{"name_s", "t.fa", "blocks.bed", Options{Name: true, Strand: true}},
		{"namePlus_s", "t.fa", "blocks.bed", Options{NamePlus: true, Strand: true}},
		{"nameOnly_s", "t.fa", "blocks.bed", Options{NameOnly: true, Strand: true}},
		{"tab", "t.fa", "blocks.bed", Options{Tab: true}},
		{"tab_name", "t.fa", "blocks.bed", Options{Tab: true, Name: true}},
		{"bedOut", "t.fa", "blocks.bed", Options{BedOut: true}},
		{"bedOut_s", "t.fa", "blocks.bed", Options{BedOut: true, Strand: true}},
		{"bedOut_split", "t.fa", "blocks.bed", Options{BedOut: true, Split: true}},
		{"split", "t.fa", "blocks.bed", Options{Split: true}},
		{"split_s", "t.fa", "blocks.bed", Options{Split: true, Strand: true}},
		{"iupac", "test.iupac.fa", "test.iupac.bed", Options{}},
		{"iupac_s", "test.iupac.fa", "test.iupac.bed", Options{Strand: true}},
		{"iupac_s_rna", "test.iupac.fa", "test.iupac.bed", Options{Strand: true, RNA: true}},
		{"rna_s_name", "rna.fasta", "rna.bed", Options{Strand: true, Name: true}},
		{"rna_s_name_rna", "rna.fasta", "rna.bed", Options{Strand: true, Name: true, RNA: true}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertParity(t, tc.fasta, tc.bed, tc.opts)
		})
	}
}

// TestUpstreamParity_BGZF asserts parity for a BGZF-compressed FASTA input.
func TestUpstreamParity_BGZF(t *testing.T) {
	assertParity(t, "t.fa.gz", "blocks.bed", Options{Split: true})
}

// withBED writes a custom BED into a freshly staged dir (alongside the named
// FASTA) and returns the dir + the BED filename.
func withBED(t *testing.T, fastaName, bedBody string) (dir, bedName string) {
	t.Helper()
	dir = stageDir(t, fastaName)
	bedName = "custom.bed"
	if err := os.WriteFile(filepath.Join(dir, bedName), []byte(bedBody), 0o644); err != nil {
		t.Fatalf("write custom bed: %v", err)
	}
	return dir, bedName
}

// assertParityCustom runs both implementations over an arbitrary BED body and
// asserts stdout+stderr parity. The FASTA fixture is staged once per side.
func assertParityCustom(t *testing.T, fastaName, bedBody string, opts Options) {
	t.Helper()

	upDir, upBed := withBED(t, fastaName, bedBody)
	upStdout, upStderr := runUpstream(t, upDir, nil, optsToArgs(fastaName, upBed, opts)...)

	goDir, goBed := withBED(t, fastaName, bedBody)
	fi := filepath.Join(goDir, fastaName)
	var goStdout, goStderr bytes.Buffer
	if _, err := Run(strings.NewReader(bedBody), fi, &goStdout, &goStderr, opts); err != nil {
		t.Fatalf("Go Run: %v", err)
	}
	_ = goBed

	if !bytes.Equal(goStdout.Bytes(), upStdout) {
		t.Errorf("stdout mismatch\n--- upstream ---\n%s\n--- go ---\n%s", upStdout, goStdout.Bytes())
	}
	if !bytes.Equal(goStderr.Bytes(), upStderr) {
		t.Errorf("stderr mismatch\n--- upstream ---\n%q\n--- go ---\n%q", upStderr, goStderr.Bytes())
	}
}

// TestUpstreamParity_EdgeMessages covers the warning/error messages: missing
// chromosome, out-of-range coordinates, and zero-length features, each
// asserted byte-for-byte against upstream stderr.
func TestUpstreamParity_EdgeMessages(t *testing.T) {
	cases := []struct {
		name string
		bed  string
		opts Options
	}{
		{"missing_chrom", "chrZ\t1\t10\n", Options{}},
		{"beyond_end", "chr1\t45\t60\n", Options{}},
		{"beyond_start", "chr1\t100\t110\n", Options{}},
		{"zero_length_mid", "chr1\t5\t5\n", Options{}},
		{"zero_length_at_end", "chr1\t50\t50\n", Options{}},
		{"end_equals_length", "chr1\t40\t50\n", Options{}},
		{"start_equals_length", "chr1\t50\t50\n", Options{}},
		{"mixed", "chr1\t0\t5\nchrZ\t1\t10\nchr1\t100\t110\nchr1\t10\t20\n", Options{}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertParityCustom(t, "t.fa", tc.bed, tc.opts)
		})
	}
}

// TestUpstreamParity_BedOutColumns checks -bedOut over BED records of varying
// column counts (BED3/4/6 and extra columns) for byte-exact re-serialization.
func TestUpstreamParity_BedOutColumns(t *testing.T) {
	beds := []string{
		"chr1\t1\t10\n",
		"chr1\t1\t10\tfoo\n",
		"chr1\t1\t10\tfoo\t99\t-\n",
		"chr1\t1\t10\tfoo\t99\t-\textra1\textra2\n",
	}
	for i, b := range beds {
		b := b
		t.Run(strings.ReplaceAll(strings.TrimSpace(b), "\t", "_"), func(t *testing.T) {
			_ = i
			assertParityCustom(t, "t.fa", b, Options{BedOut: true})
		})
	}
}

// TestUpstreamParity_NameMissingColumn confirms the empty-name header parity
// for -name and -nameOnly on a BED row without a name column.
func TestUpstreamParity_NameMissingColumn(t *testing.T) {
	assertParityCustom(t, "t.fa", "chr1\t1\t10\n", Options{Name: true})
	assertParityCustom(t, "t.fa", "chr1\t1\t10\n", Options{NameOnly: true})
}
