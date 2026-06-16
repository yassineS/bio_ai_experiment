package bedclosest

// Live-upstream parity tests for the `closest` sortAndNaming sub-suite
// (upstream test-sort-and-naming.sh, ids closest.t01-t23). They build the real
// upstream `bedtools` binary once (via the shared sync.Once builder in
// gaps_parity_test.go) and compare its `closest` stdout, stderr, AND exit code
// byte-for-byte against this port's cross-file chromosome sort-order /
// naming-convention validation (validate.go). They t.Fatalf (never t.Skip) so a
// missing/unbuildable submodule is a hard failure, matching the parity-rig
// policy.
//
// The upstream messages echo each file's command-line name, so both upstream
// and the port are driven with the bare fixture basenames (upstream is run with
// its working directory set to the sortAndNaming fixture dir).

import (
	"bytes"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// sortNamingDir returns the absolute path to the vendored sortAndNaming
// fixtures.
func sortNamingDir(t *testing.T) string {
	t.Helper()
	return filepath.Join(gapsRepoRoot(t), "reference_code", "bedtools", "test", "closest", "sortAndNaming")
}

// runUpstreamClosestRaw runs upstream `bedtools closest` with the given bare
// argument filenames, with its working directory set to dir, and returns its
// stdout, stderr, and exit code.
func runUpstreamClosestRaw(t *testing.T, bt, dir, aFile string, bFiles []string, flags ...string) (stdout, stderr []byte, exit int) {
	t.Helper()
	args := []string{"closest", "-a", aFile, "-b"}
	args = append(args, bFiles...)
	args = append(args, flags...)
	cmd := exec.Command(bt, args...)
	cmd.Dir = dir
	var outBuf, errBuf bytes.Buffer
	cmd.Stdout = &outBuf
	cmd.Stderr = &errBuf
	err := cmd.Run()
	exit = 0
	if err != nil {
		var ee *exec.ExitError
		if errors.As(err, &ee) {
			exit = ee.ExitCode()
		} else {
			t.Fatalf("upstream closest %v failed to run: %v", args, err)
		}
	}
	return outBuf.Bytes(), errBuf.Bytes(), exit
}

// runOurClosestRaw drives this port's ClosestMulti over the named fixtures with
// validation enabled, returning its stdout, stderr (warn writer), and the exit
// code it would map to (1 on a validation ERROR, else 0).
func runOurClosestRaw(t *testing.T, dir, aFile string, bFiles []string, opts Options) (stdout, stderr []byte, exit int) {
	t.Helper()
	read := func(name string) io.Reader {
		data, err := os.ReadFile(filepath.Join(dir, name))
		if err != nil {
			t.Fatalf("read fixture %s: %v", name, err)
		}
		return bytes.NewReader(data)
	}
	readers := make([]io.Reader, len(bFiles))
	for i, b := range bFiles {
		readers[i] = read(b)
	}
	var out, warn bytes.Buffer
	opts.WarnWriter = &warn
	opts.QueryName = aFile
	opts.DBNames = bFiles
	_, err := ClosestMulti(read(aFile), readers, &out, opts)
	exit = 0
	if err != nil {
		exit = 1
	}
	return out.Bytes(), warn.Bytes(), exit
}

// TestParity_Closest_SortAndNaming reproduces every case of the upstream
// test-sort-and-naming.sh script, asserting byte-for-byte parity of stdout,
// stderr, and exit code with the live binary.
func TestParity_Closest_SortAndNaming(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)
	dir := sortNamingDir(t)

	type snCase struct {
		name string
		a    string
		b    []string
	}
	// The whole closest.t01-t23 matrix. The mixed lexicographic/numeric and
	// chr/non-chr naming permutations exercise every branch of the validator:
	// out-of-order ERROR, lexico-disproven-while-assumed ERROR, db-chrom-not-in
	// -query ERROR, chr-prefix WARNING, and leading-zero WARNING, plus the clean
	// (exit 0) permutations.
	cases := []snCase{
		{"t01", "sq1.bed", []string{"sdb1.bed"}},
		{"t02", "q1a_num.bed", []string{"db1_num.bed", "db2_numBackwards.bed"}},
		{"t03", "q1a_num.bed", []string{"db1_num.bed", "db2_num.bed", "db3_numBackwards.bed"}},
		{"t04", "q1_num.bed", []string{"db1_num.bed", "db2_num.bed"}},
		{"t05", "q1_num.bed", []string{"db1_noChr.bed"}},
		{"t06", "q1_num.bed", []string{"db1_leadingZero.txt"}},
		{"t07", "q1_gls.bed", []string{"q1_gls.bed"}},
		{"t08", "alpha_all.bed", []string{"alpha_all.bed"}},
		{"t09", "alpha_all.bed", []string{"alpha_missing.bed"}},
		{"t10", "alpha_all.bed", []string{"num_all.bed"}},
		{"t11", "alpha_all.bed", []string{"num_missing.bed"}},
		{"t12", "alpha_missing.bed", []string{"alpha_all.bed"}},
		{"t13", "alpha_missing.bed", []string{"alpha_missing.bed"}},
		{"t14", "alpha_missing.bed", []string{"num_all.bed"}},
		{"t15", "alpha_missing.bed", []string{"num_missing.bed"}},
		{"t16", "num_all.bed", []string{"alpha_all.bed"}},
		{"t17", "num_all.bed", []string{"alpha_missing.bed"}},
		{"t18", "num_all.bed", []string{"num_all.bed"}},
		{"t19", "num_all.bed", []string{"num_missing.bed"}},
		{"t20", "num_missing.bed", []string{"alpha_all.bed"}},
		{"t21", "num_missing.bed", []string{"alpha_missing.bed"}},
		{"t22", "num_missing.bed", []string{"num_all.bed"}},
		{"t23", "num_missing.bed", []string{"num_missing.bed"}},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			wantOut, wantErr, wantExit := runUpstreamClosestRaw(t, bt, dir, c.a, c.b)
			gotOut, gotErr, gotExit := runOurClosestRaw(t, dir, c.a, c.b, Options{})
			if !bytes.Equal(gotOut, wantOut) {
				t.Fatalf("%s stdout mismatch\nupstream:\n%q\nours:\n%q", c.name, wantOut, gotOut)
			}
			if !bytes.Equal(gotErr, wantErr) {
				t.Fatalf("%s stderr mismatch\nupstream:\n%q\nours:\n%q", c.name, wantErr, gotErr)
			}
			if gotExit != wantExit {
				t.Fatalf("%s exit mismatch: upstream=%d ours=%d", c.name, wantExit, gotExit)
			}
		})
	}
}
