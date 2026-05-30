package seqtk

// Live-binary oracle sweep across every seqtk subcommand the Go port
// implements. Each row in the table runs an identical CLI invocation
// through (a) the genuine upstream binary at reference_code/seqtk/seqtk
// (built from the v1.5-r133 submodule), and (b) our Go binary built
// fresh in TestMain. Outputs are compared byte-for-byte.
//
// When the upstream binary is absent (e.g. CI without submodules), the
// whole test t.Skip's so existing parity coverage still passes. When a
// specific fixture for a parity-relevant invocation isn't yet checked
// into testdata/parity/, the affected subtest t.Skip's with a TODO.
// When the fixture IS present but our port diverges, t.Errorf surfaces
// the gap (we never silently swallow real divergence).
//
// Goal: extend the live-oracle proof-of-concept demonstrated by
// seq_live_test.go (`TestSeqtkSeqLiveOracle`) to all 25 subcommands.

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// ourBin is the path to our freshly-built Go seqtk binary, populated by
// TestMain before any live-oracle subtest runs. Empty string means the
// build failed; tests will skip rather than crash.
var ourBin string

// TestMain builds tools/seqtk/cmd/seqtk into a temp dir once for the
// entire package test run so the live-oracle table can exec it.
func TestMain(m *testing.M) {
	tmp, err := os.MkdirTemp("", "seqtk-live-bin-")
	if err != nil {
		fmt.Fprintf(os.Stderr, "live-oracle: mkdtemp: %v\n", err)
		os.Exit(m.Run())
	}
	defer os.RemoveAll(tmp)

	bin := filepath.Join(tmp, "seqtk")
	cmd := exec.Command("go", "build", "-o", bin, "github.com/yassineS/bio_ai_experiment/tools/seqtk/cmd/seqtk")
	if out, err := cmd.CombinedOutput(); err != nil {
		fmt.Fprintf(os.Stderr, "live-oracle: go build failed: %v\n%s\n", err, out)
	} else {
		ourBin = bin
	}

	os.Exit(m.Run())
}

// parityDir returns the absolute path to the seqtk parity-fixture dir.
// It is computed relative to the test binary's working directory
// (tools/seqtk/pkg/seqtk).
func parityDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity"))
	if err != nil {
		t.Fatalf("parityDir abs: %v", err)
	}
	return abs
}

// runBin execs the binary at `bin` with the given args, returning
// stdout, stderr, and any execution error. The args slice is taken
// verbatim. CWD is set to `cwd` so relative paths in args resolve there.
func runBin(bin string, cwd string, args ...string) (stdout, stderr []byte, err error) {
	cmd := exec.Command(bin, args...)
	cmd.Dir = cwd
	var so, se bytes.Buffer
	cmd.Stdout = &so
	cmd.Stderr = &se
	err = cmd.Run()
	return so.Bytes(), se.Bytes(), err
}

// liveCase describes one subcommand invocation. `args` is the seqtk
// argv tail (everything after "seqtk") fed to our Go binary.
// `liveArgs`, when non-nil, overrides the argv for the upstream binary
// — needed where upstream's flag syntax diverges from our port's (e.g.
// `seqtk hrun in.fa 2` vs our `seqtk hrun -l 2 in.fa`, or fq2fa which
// upstream lacks entirely — there it's `seqtk seq -A`). `compareStderr`
// asks the test to additionally diff stderr (used by e.g. `telo` which
// writes its summary there, and `mergefa`).
type liveCase struct {
	name          string
	args          []string
	liveArgs      []string
	compareStderr bool
	// fixtureRequired lists fixture filenames (relative to the parity
	// dir) that must exist for the subtest to run. If any are missing
	// the subtest t.Skip's with a TODO.
	fixtureRequired []string
	// skipReason, if non-empty, is logged via t.Skip without running
	// the subcommand. Used for subcommands whose port intentionally
	// diverges from upstream in ways already captured elsewhere.
	skipReason string
}

// liveCases enumerates the full sweep across the 25 subcommands.
// Fixtures live under tools/seqtk/testdata/parity/ — see
// liveCase.fixtureRequired for the inputs each row consumes.
func liveCases() []liveCase {
	return []liveCase{
		// --- comp ---
		{name: "comp_small_fa", args: []string{"comp", "small.fa"}, fixtureRequired: []string{"small.fa"}},
		{name: "comp_small_fq", args: []string{"comp", "small.fq"}, fixtureRequired: []string{"small.fq"}},
		{name: "comp_nruns_fa", args: []string{"comp", "nruns.fa"}, fixtureRequired: []string{"nruns.fa"}},
		{name: "comp_u_fa", args: []string{"comp", "comp_u.fa"}, fixtureRequired: []string{"comp_u.fa"}},

		// --- fq2fa --- upstream lacks this subcommand; the documented
		// upstream equivalent is `seqtk seq -A`. Compare our `fq2fa`
		// against that.
		{name: "fq2fa_small", args: []string{"fq2fa", "small.fq"}, liveArgs: []string{"seq", "-A", "small.fq"}, fixtureRequired: []string{"small.fq"}},
		{name: "fq2fa_p64", args: []string{"fq2fa", "p64.fq"}, liveArgs: []string{"seq", "-A", "p64.fq"}, fixtureRequired: []string{"p64.fq"}},

		// --- size ---
		{name: "size_small_fa", args: []string{"size", "small.fa"}, fixtureRequired: []string{"small.fa"}},
		{name: "size_small_fq", args: []string{"size", "small.fq"}, fixtureRequired: []string{"small.fq"}},
		{name: "size_nruns_fa", args: []string{"size", "nruns.fa"}, fixtureRequired: []string{"nruns.fa"}},

		// --- gap --- (BED rows)
		{name: "gap_small_l3", args: []string{"gap", "-l", "3", "gap_small.fa"}, fixtureRequired: []string{"gap_small.fa"}},
		{name: "gap_small_l5", args: []string{"gap", "-l", "5", "gap_small.fa"}, fixtureRequired: []string{"gap_small.fa"}},
		{name: "gap_nruns_l1", args: []string{"gap", "-l", "1", "nruns.fa"}, fixtureRequired: []string{"nruns.fa"}},

		// --- gc --- (BED rows)
		{name: "gc_small_default", args: []string{"gc", "gc_small.fa"}, fixtureRequired: []string{"gc_small.fa"}},
		{name: "gc_small_f07_l10", args: []string{"gc", "-f", "0.7", "-l", "10", "gc_small.fa"}, fixtureRequired: []string{"gc_small.fa"}},
		{name: "gc_small_w_f07_l10", args: []string{"gc", "-w", "-f", "0.7", "-l", "10", "gc_small.fa"}, fixtureRequired: []string{"gc_small.fa"}},

		// --- hpc ---
		{name: "hpc_homo", args: []string{"hpc", "homo.fa"}, fixtureRequired: []string{"homo.fa"}},
		{name: "hpc_small", args: []string{"hpc", "small.fa"}, fixtureRequired: []string{"small.fa"}},

		// --- hrun --- upstream uses a positional `[minLen]` argument
		// while our port uses `-l INT`. liveArgs encodes that
		// translation.
		{name: "hrun_nruns_default", args: []string{"hrun", "nruns.fa"}, fixtureRequired: []string{"nruns.fa"}},
		{name: "hrun_nruns_l2", args: []string{"hrun", "-l", "2", "nruns.fa"}, liveArgs: []string{"hrun", "nruns.fa", "2"}, fixtureRequired: []string{"nruns.fa"}},
		{name: "hrun_nruns_l3", args: []string{"hrun", "-l", "3", "nruns.fa"}, liveArgs: []string{"hrun", "nruns.fa", "3"}, fixtureRequired: []string{"nruns.fa"}},

		// --- hety ---
		{name: "hety_basic_default", args: []string{"hety", "hety_basic.fa"}, fixtureRequired: []string{"hety_basic.fa"}},
		{name: "hety_basic_w30", args: []string{"hety", "-w", "30", "hety_basic.fa"}, fixtureRequired: []string{"hety_basic.fa"}},
		{name: "hety_basic_w30_t3", args: []string{"hety", "-w", "30", "-t", "3", "hety_basic.fa"}, fixtureRequired: []string{"hety_basic.fa"}},
		{name: "hety_basic_w30_m", args: []string{"hety", "-w", "30", "-m", "hety_basic.fa"}, fixtureRequired: []string{"hety_basic.fa"}},
		{name: "hety_lowercase_w6_t1_m", args: []string{"hety", "-w", "6", "-t", "1", "-m", "hety_lowercase.fa"}, fixtureRequired: []string{"hety_lowercase.fa"}},

		// --- listhet ---
		{name: "listhet_small", args: []string{"listhet", "small.fa"}, fixtureRequired: []string{"small.fa"}},
		{name: "listhet_ambig", args: []string{"listhet", "ambig.fa"}, fixtureRequired: []string{"ambig.fa"}},
		{name: "listhet_hety_basic", args: []string{"listhet", "hety_basic.fa"}, fixtureRequired: []string{"hety_basic.fa"}},
		{name: "listhet_hety_lowercase", args: []string{"listhet", "hety_lowercase.fa"}, fixtureRequired: []string{"hety_lowercase.fa"}},

		// --- kfreq ---
		{name: "kfreq_small_AC", args: []string{"kfreq", "AC", "kfreq_small.fa"}, fixtureRequired: []string{"kfreq_small.fa"}},
		{name: "kfreq_small_AAAA", args: []string{"kfreq", "AAAA", "kfreq_small.fa"}, fixtureRequired: []string{"kfreq_small.fa"}},
		{name: "kfreq_small_ACGT", args: []string{"kfreq", "ACGT", "kfreq_small.fa"}, fixtureRequired: []string{"kfreq_small.fa"}},
		{name: "kfreq_mixed_AA", args: []string{"kfreq", "AA", "kfreq_mixed.fa"}, fixtureRequired: []string{"kfreq_mixed.fa"}},
		{name: "kfreq_mixed_ACGT", args: []string{"kfreq", "ACGT", "kfreq_mixed.fa"}, fixtureRequired: []string{"kfreq_mixed.fa"}},
		{name: "kfreq_mixed_CCCTAA", args: []string{"kfreq", "CCCTAA", "kfreq_mixed.fa"}, fixtureRequired: []string{"kfreq_mixed.fa"}},
		{name: "kfreq_mixed_CCGG", args: []string{"kfreq", "CCGG", "kfreq_mixed.fa"}, fixtureRequired: []string{"kfreq_mixed.fa"}},
		{name: "kfreq_edge_AC", args: []string{"kfreq", "AC", "kfreq_edge.fa"}, fixtureRequired: []string{"kfreq_edge.fa"}},

		// --- telo --- (writes BED to stdout, summary to stderr)
		{name: "telo_basic_default", args: []string{"telo", "telo_basic.fa"}, compareStderr: true, fixtureRequired: []string{"telo_basic.fa"}},
		{name: "telo_basic_mTTAGGG", args: []string{"telo", "-m", "TTAGGG", "telo_basic.fa"}, compareStderr: true, fixtureRequired: []string{"telo_basic.fa"}},
		{name: "telo_complex_default", args: []string{"telo", "telo_complex.fa"}, compareStderr: true, fixtureRequired: []string{"telo_complex.fa"}},
		{name: "telo_complex_s100", args: []string{"telo", "-s", "100", "telo_complex.fa"}, compareStderr: true, fixtureRequired: []string{"telo_complex.fa"}},
		{name: "telo_complex_p2d500", args: []string{"telo", "-p", "2", "-d", "500", "telo_complex.fa"}, compareStderr: true, fixtureRequired: []string{"telo_complex.fa"}},
		{name: "telo_edge_default", args: []string{"telo", "telo_edge.fa"}, compareStderr: true, fixtureRequired: []string{"telo_edge.fa"}},

		// --- cutN ---
		{name: "cutN_n4_nruns", args: []string{"cutN", "-n", "4", "nruns.fa"}, fixtureRequired: []string{"nruns.fa"}},
		{name: "cutN_n100_nruns", args: []string{"cutN", "-n", "100", "nruns.fa"}, fixtureRequired: []string{"nruns.fa"}},
		{name: "cutN_n5_gap_small", args: []string{"cutN", "-n", "5", "gap_small.fa"}, fixtureRequired: []string{"gap_small.fa"}},

		// --- dropse ---
		{name: "dropse_input_fq", args: []string{"dropse", "dropse_input.fq"}, fixtureRequired: []string{"dropse_input.fq"}},
		{name: "dropse_input_fa", args: []string{"dropse", "dropse_input.fa"}, fixtureRequired: []string{"dropse_input.fa"}},

		// --- famask ---
		// Argument order: famask <src> <mask>
		{name: "famask_main", args: []string{"famask", "famask_src.fa", "famask_mask.fa"}, fixtureRequired: []string{"famask_src.fa", "famask_mask.fa"}},
		{name: "famask_simple", args: []string{"famask", "famask_simple_src.fa", "famask_simple_mask.fa"}, fixtureRequired: []string{"famask_simple_src.fa", "famask_simple_mask.fa"}},

		// --- mergefa --- (writes summary line to stderr)
		{name: "mergefa_default", args: []string{"mergefa", "mergefa_a.fa", "mergefa_b.fa"}, compareStderr: true, fixtureRequired: []string{"mergefa_a.fa", "mergefa_b.fa"}},
		{name: "mergefa_i", args: []string{"mergefa", "-i", "mergefa_a.fa", "mergefa_b.fa"}, compareStderr: true, fixtureRequired: []string{"mergefa_a.fa", "mergefa_b.fa"}},
		{name: "mergefa_m", args: []string{"mergefa", "-m", "mergefa_a.fa", "mergefa_b.fa"}, compareStderr: true, fixtureRequired: []string{"mergefa_a.fa", "mergefa_b.fa"}},
		{name: "mergefa_h", args: []string{"mergefa", "-h", "mergefa_a.fa", "mergefa_b.fa"}, compareStderr: true, fixtureRequired: []string{"mergefa_a.fa", "mergefa_b.fa"}},
		{name: "mergefa_q20", args: []string{"mergefa", "-q", "20", "mergefa_a.fq", "mergefa_b.fq"}, compareStderr: true, fixtureRequired: []string{"mergefa_a.fq", "mergefa_b.fq"}},
		{name: "mergefa_long", args: []string{"mergefa", "mergefa_long_a.fa", "mergefa_long_b.fa"}, compareStderr: true, fixtureRequired: []string{"mergefa_long_a.fa", "mergefa_long_b.fa"}},

		// --- mergepe ---
		{name: "mergepe_pe", args: []string{"mergepe", "pe1.fq", "pe2.fq"}, fixtureRequired: []string{"pe1.fq", "pe2.fq"}},

		// --- mutfa ---
		{name: "mutfa_ref", args: []string{"mutfa", "ref.fa", "muts.tsv"}, fixtureRequired: []string{"ref.fa", "muts.tsv"}},

		// --- randbase --- upstream takes no flags and uses drand48
		// with implicit seed 0. Our CLI time-seeds by default; force
		// our `-s 0` to match upstream's implicit-0 path while
		// invoking upstream bare.
		{name: "randbase_ambig_s0", args: []string{"randbase", "-s", "0", "ambig.fa"}, liveArgs: []string{"randbase", "ambig.fa"}, fixtureRequired: []string{"ambig.fa"}},

		// --- rename ---
		{name: "rename_small_fa_noprefix", args: []string{"rename", "small.fa"}, fixtureRequired: []string{"small.fa"}},
		{name: "rename_small_fa_prefix", args: []string{"rename", "small.fa", "PX"}, fixtureRequired: []string{"small.fa"}},
		{name: "rename_small_fq_noprefix", args: []string{"rename", "small.fq"}, fixtureRequired: []string{"small.fq"}},
		{name: "rename_pairs_fa", args: []string{"rename", "rename_pairs.fa", "SAMPLE_"}, fixtureRequired: []string{"rename_pairs.fa"}},
		{name: "rename_pairs_fq", args: []string{"rename", "rename_pairs.fq", "SAMPLE_"}, fixtureRequired: []string{"rename_pairs.fq"}},

		// --- sample --- default seed (11), streaming kr_drand < frac
		{name: "sample_sample20_f0.3", args: []string{"sample", "sample20.fq", "0.3"}, fixtureRequired: []string{"sample20.fq"}},
		{name: "sample_sample20_f0.5", args: []string{"sample", "sample20.fq", "0.5"}, fixtureRequired: []string{"sample20.fq"}},

		// --- split --- handled specially (multi-file output); see runSplitOracle.
		// Listed here so the report covers it.
		{name: "split_handled_separately", args: nil, skipReason: "split is exercised via TestSeqtkSplitLiveOracle (multi-file diff)"},

		// --- subseq ---
		{name: "subseq_names", args: []string{"subseq", "names.fa", "namelist.txt"}, fixtureRequired: []string{"names.fa", "namelist.txt"}},
		{name: "subseq_bed", args: []string{"subseq", "bed.fa", "regions.bed"}, fixtureRequired: []string{"bed.fa", "regions.bed"}},

		// --- trimfq --- default Modified Mott trim (-q 0.05 -l 30)
		{name: "trimfq_default_p64", args: []string{"trimfq", "p64.fq"}, fixtureRequired: []string{"p64.fq"}},
		{name: "trimfq_borders", args: []string{"trimfq", "trimfq_borders.fq"}, fixtureRequired: []string{"trimfq_borders.fq"}},

		// --- fqchk ---
		{name: "fqchk_mixed_default", args: []string{"fqchk", "fqchk_mixed.fq"}, fixtureRequired: []string{"fqchk_mixed.fq"}},
		{name: "fqchk_mixed_q0", args: []string{"fqchk", "-q", "0", "fqchk_mixed.fq"}, fixtureRequired: []string{"fqchk_mixed.fq"}},
		{name: "fqchk_mixed_q30", args: []string{"fqchk", "-q", "30", "fqchk_mixed.fq"}, fixtureRequired: []string{"fqchk_mixed.fq"}},
		{name: "fqchk_small_default", args: []string{"fqchk", "small.fq"}, fixtureRequired: []string{"small.fq"}},
		{name: "fqchk_small_q0", args: []string{"fqchk", "-q", "0", "small.fq"}, fixtureRequired: []string{"small.fq"}},

		// --- seq --- exercised exhaustively by seq_live_test.go; spot-check
		// one invocation here so this sweep covers all 25 named subcommands.
		{name: "seq_small_fa_passthrough", args: []string{"seq", "small.fa"}, fixtureRequired: []string{"small.fa"}},
	}
}

// TestSeqtkAllSubcommandsLive is the table-driven live-oracle sweep.
// Each row drives both the genuine upstream binary and our freshly
// built Go binary with identical args from the parity-fixture
// directory and compares stdout (and optionally stderr) byte-for-byte.
func TestSeqtkAllSubcommandsLive(t *testing.T) {
	live := findGenuineSeqtk()
	if live == "" {
		t.Skip("genuine seqtk binary not found; skipping live-oracle sweep")
	}
	if ourBin == "" {
		t.Skip("our seqtk binary failed to build in TestMain; skipping live-oracle sweep")
	}
	// findGenuineSeqtk returns a path that's only valid relative to
	// the package's CWD (tools/seqtk/pkg/seqtk). Resolve to absolute
	// so we can run the binary from anywhere.
	liveAbs, err := filepath.Abs(live)
	if err != nil {
		t.Fatalf("abs(live): %v", err)
	}

	cwd := parityDir(t)
	for _, tc := range liveCases() {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			if tc.skipReason != "" {
				t.Skip(tc.skipReason)
			}
			for _, f := range tc.fixtureRequired {
				if _, err := os.Stat(filepath.Join(cwd, f)); err != nil {
					t.Skipf("TODO: parity fixture %q missing: %v", f, err)
				}
			}
			runOne(t, liveAbs, ourBin, cwd, tc)
		})
	}
}

// runOne executes a single liveCase against both binaries and reports.
func runOne(t *testing.T, live, ours, cwd string, tc liveCase) {
	t.Helper()
	liveInvoke := tc.liveArgs
	if liveInvoke == nil {
		liveInvoke = tc.args
	}
	liveOut, liveErr, err1 := runBin(live, cwd, liveInvoke...)
	ourOut, ourErr, err2 := runBin(ours, cwd, tc.args...)

	// Exit-code parity: both should agree on success/failure. We don't
	// require equal exit codes (Go's flag package uses 2 on misuse,
	// upstream uses 1), only that both succeed or both fail.
	liveOK := err1 == nil
	ourOK := err2 == nil
	if liveOK != ourOK {
		t.Errorf("exit-status divergence: live success=%v (err=%v, stderr=%q), ours success=%v (err=%v, stderr=%q)",
			liveOK, err1, string(liveErr), ourOK, err2, string(ourErr))
		return
	}

	if !bytes.Equal(ourOut, liveOut) {
		t.Errorf("stdout differs (%d vs %d bytes)\n--- live ---\n%s\n--- ours ---\n%s",
			len(liveOut), len(ourOut), truncate(liveOut, 2000), truncate(ourOut, 2000))
	}
	if tc.compareStderr {
		if !bytes.Equal(ourErr, liveErr) {
			t.Errorf("stderr differs\n--- live ---\n%s\n--- ours ---\n%s",
				truncate(liveErr, 2000), truncate(ourErr, 2000))
		}
	}
}

// truncate returns a copy of b cut at n bytes for less-noisy diffs.
func truncate(b []byte, n int) []byte {
	if len(b) <= n {
		return b
	}
	cut := make([]byte, 0, n+20)
	cut = append(cut, b[:n]...)
	cut = append(cut, []byte("\n...[truncated]")...)
	return cut
}

// TestSeqtkSplitLiveOracle drives `seqtk split -n N prefix in.fa` through
// both binaries into separate temp dirs and diffs every output file. The
// upstream binary insists on writing into its CWD, so each invocation
// runs in its own tempdir with the fixture symlinked in.
func TestSeqtkSplitLiveOracle(t *testing.T) {
	live := findGenuineSeqtk()
	if live == "" {
		t.Skip("genuine seqtk binary not found")
	}
	if ourBin == "" {
		t.Skip("our seqtk binary failed to build in TestMain")
	}
	liveAbs, err := filepath.Abs(live)
	if err != nil {
		t.Fatalf("abs(live): %v", err)
	}

	cases := []struct {
		name    string
		input   string
		n       string
		lineLen string // empty means omit
	}{
		{name: "fa_n2", input: "small.fa", n: "2"},
		{name: "fa_n3", input: "small.fa", n: "3"},
		{name: "fa_n2_l5", input: "small.fa", n: "2", lineLen: "5"},
		{name: "fq_n2", input: "small.fq", n: "2"},
		{name: "fq_n2_l4", input: "small.fq", n: "2", lineLen: "4"},
	}

	parity := parityDir(t)
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			src := filepath.Join(parity, tc.input)
			if _, err := os.Stat(src); err != nil {
				t.Skipf("fixture %s missing: %v", tc.input, err)
			}
			liveDir := t.TempDir()
			ourDir := t.TempDir()
			if err := copyFile(src, filepath.Join(liveDir, tc.input)); err != nil {
				t.Fatalf("stage live fixture: %v", err)
			}
			if err := copyFile(src, filepath.Join(ourDir, tc.input)); err != nil {
				t.Fatalf("stage ours fixture: %v", err)
			}

			args := []string{"split", "-n", tc.n}
			if tc.lineLen != "" {
				args = append(args, "-l", tc.lineLen)
			}
			args = append(args, "out", tc.input)

			if _, se, err := runBin(liveAbs, liveDir, args...); err != nil {
				t.Fatalf("live split failed: %v stderr=%q", err, string(se))
			}
			if _, se, err := runBin(ourBin, ourDir, args...); err != nil {
				t.Fatalf("our split failed: %v stderr=%q", err, string(se))
			}

			liveFiles := listOutputs(t, liveDir, "out.")
			ourFiles := listOutputs(t, ourDir, "out.")
			if !stringSlicesEqual(liveFiles, ourFiles) {
				t.Errorf("output file set differs: live=%v ours=%v", liveFiles, ourFiles)
				return
			}
			for _, name := range liveFiles {
				lb, err := os.ReadFile(filepath.Join(liveDir, name))
				if err != nil {
					t.Fatalf("read live %s: %v", name, err)
				}
				ob, err := os.ReadFile(filepath.Join(ourDir, name))
				if err != nil {
					t.Fatalf("read ours %s: %v", name, err)
				}
				if !bytes.Equal(lb, ob) {
					t.Errorf("split %s file %s differs (%d vs %d bytes)",
						tc.name, name, len(lb), len(ob))
				}
			}
		})
	}
}

// copyFile copies src to dst, creating parent dirs as needed.
func copyFile(src, dst string) error {
	in, err := os.Open(src)
	if err != nil {
		return err
	}
	defer in.Close()
	if err := os.MkdirAll(filepath.Dir(dst), 0o755); err != nil {
		return err
	}
	out, err := os.Create(dst)
	if err != nil {
		return err
	}
	defer out.Close()
	if _, err := io.Copy(out, in); err != nil {
		return err
	}
	return nil
}

// listOutputs returns sorted basenames in dir that begin with prefix.
func listOutputs(t *testing.T, dir, prefix string) []string {
	t.Helper()
	ents, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("read dir %s: %v", dir, err)
	}
	var out []string
	for _, e := range ents {
		if e.IsDir() {
			continue
		}
		if strings.HasPrefix(e.Name(), prefix) {
			out = append(out, e.Name())
		}
	}
	sort.Strings(out)
	return out
}

// stringSlicesEqual reports whether a and b contain the same strings
// in the same order. Used by the split oracle.
func stringSlicesEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
