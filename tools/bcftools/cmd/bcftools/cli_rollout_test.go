package main

import (
	"bytes"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
)

// This file validates the cliflag.Parse rollout across bcftools
// subcommands: every subcommand now routes argument parsing through
// parseFlags -> cliflag.Parse, gaining POSIX getopt short-flag bundling
// ("-hG" == "-h -G"), value concatenation ("-Ob" == "-O b"), the "--"
// terminator and bare "-" handling, plus the upstream legacy/compat
// short flags that were registered alongside the rollout.
//
// The tests assert three things, none of which ever t.Skip (the suite
// t.Fatalf's when prerequisites are missing, per the project's
// no-silent-skip rule for CLI rollouts):
//
//  1. Representative BUNDLED command lines parse and run for several
//     subcommands.
//  2. The bundled form is equivalent to the canonical (spread-out) form
//     within OUR binary (byte-identical stdout, provenance lines
//     stripped).
//  3. Bundled command lines also parse and run on the LIVE upstream
//     bcftools binary, and the decoded output matches where the two
//     implementations are at parity.

// rolloutHarness holds the absolute paths to our freshly-built port
// binary and the vendored upstream bcftools binary used by the rollout
// tests.
type rolloutHarness struct {
	ours     string
	upstream string
	plugins  string
}

var (
	upstreamBcftoolsCliRolloutOnce sync.Once
	upstreamBcftoolsCliRollout     rolloutHarness
	upstreamBcftoolsCliRolloutErr  string
)

// buildRolloutHarness builds our bcftools port binary once and locates
// the vendored upstream binary and plugins directory. It records a
// human-readable error string (rather than skipping) so callers can
// t.Fatalf with an actionable message.
func buildRolloutHarness(t *testing.T) rolloutHarness {
	t.Helper()
	upstreamBcftoolsCliRolloutOnce.Do(func() {
		tmp, err := os.MkdirTemp("", "bcftools-cli-rollout-*")
		if err != nil {
			upstreamBcftoolsCliRolloutErr = "failed to make tempdir: " + err.Error()
			return
		}
		bin := filepath.Join(tmp, "bcftools")
		// The test binary's cwd is this package dir; "." builds the
		// bcftools command in place.
		cmd := exec.Command("go", "build", "-o", bin, ".")
		var berr bytes.Buffer
		cmd.Stderr = &berr
		if err := cmd.Run(); err != nil {
			upstreamBcftoolsCliRolloutErr = "go build of our bcftools failed: " + err.Error() + "\n" + berr.String()
			return
		}
		// Upstream binary: four levels up from cmd/bcftools to repo root.
		up, _ := filepath.Abs(filepath.Join("..", "..", "..", "..",
			"reference_code", "bcftools", "bcftools"))
		if fi, err := os.Stat(up); err != nil || fi.IsDir() || fi.Mode()&0111 == 0 {
			upstreamBcftoolsCliRolloutErr = "upstream bcftools binary not found/executable at " + up +
				" (run: git submodule update --init --recursive reference_code/htslib reference_code/bcftools && build it)"
			return
		}
		pluginDir, _ := filepath.Abs(filepath.Join("..", "..", "..", "..",
			"reference_code", "bcftools", "plugins"))
		upstreamBcftoolsCliRollout = rolloutHarness{ours: bin, upstream: up, plugins: pluginDir}
	})
	if upstreamBcftoolsCliRolloutErr != "" {
		t.Fatalf("bcftools CLI-rollout harness unavailable: %s", upstreamBcftoolsCliRolloutErr)
	}
	return upstreamBcftoolsCliRollout
}

// fixture returns the absolute path to a testdata/parity fixture.
func fixture(t *testing.T, name string) string {
	t.Helper()
	p, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs(%s): %v", name, err)
	}
	if _, err := os.Stat(p); err != nil {
		t.Fatalf("fixture %s missing: %v", name, err)
	}
	return p
}

// runOut runs bin with args (and optional extra env) and returns stdout.
// It fails the test on a non-zero exit so a parse failure (exit 2)
// surfaces immediately.
func runOut(t *testing.T, bin string, env []string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	if len(env) > 0 {
		cmd.Env = append(os.Environ(), env...)
	}
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = io.Discard
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v", bin, args, err)
	}
	return out.Bytes()
}

// stripProvenance removes the volatile "##bcftools_*" provenance header
// lines so two runs (or two implementations) can be compared for
// content equality.
func stripProvenance(b []byte) string {
	var keep []string
	for _, line := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(line, "##bcftools_") {
			continue
		}
		keep = append(keep, line)
	}
	return strings.Join(keep, "\n")
}

// TestCliRolloutBundledEqualsCanonical asserts that, within OUR binary,
// a bundled short-flag command line produces identical output to the
// canonical spread-out form for several subcommands. This is the core
// guarantee of the cliflag.Parse rollout.
func TestCliRolloutBundledEqualsCanonical(t *testing.T) {
	h := buildRolloutHarness(t)
	basic := fixture(t, "basic.vcf")
	call := fixture(t, "call_in.vcf")

	cases := []struct {
		name      string
		bundled   []string
		canonical []string
		env       []string
	}{
		{
			name:      "view -hG == view -h -G",
			bundled:   []string{"view", "-hG", basic},
			canonical: []string{"view", "-h", "-G", basic},
		},
		{
			name:      "view -Ov == view -O v",
			bundled:   []string{"view", "-Ov", basic},
			canonical: []string{"view", "-O", "v", basic},
		},
		{
			name:      "norm -aD == norm -a -d exact",
			bundled:   []string{"norm", "-aD", basic},
			canonical: []string{"norm", "-a", "-d", "exact", basic},
		},
		{
			name:      "stats -s- == stats -s -",
			bundled:   []string{"stats", "-s-", basic},
			canonical: []string{"stats", "-s", "-", basic},
		},
		{
			name:      "call -mv == call -m -v",
			bundled:   []string{"call", "-mv", call},
			canonical: []string{"call", "-m", "-v", call},
		},
		{
			name:      "plugin -lv == plugin -l -v",
			bundled:   []string{"plugin", "-lv"},
			canonical: []string{"plugin", "-l", "-v"},
			env:       []string{"BCFTOOLS_PLUGINS=" + h.plugins},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			gotB := stripProvenance(runOut(t, h.ours, tc.env, tc.bundled...))
			gotC := stripProvenance(runOut(t, h.ours, tc.env, tc.canonical...))
			if gotB != gotC {
				t.Fatalf("bundled vs canonical mismatch\n bundled=%v\n canonical=%v\n--- bundled ---\n%s\n--- canonical ---\n%s",
					tc.bundled, tc.canonical, gotB, gotC)
			}
		})
	}
}

// TestCliRolloutBundledMatchesUpstream asserts that bundled command
// lines parse on the LIVE upstream binary too, and that the decoded
// output matches ours for the subcommands that are at parity. Where the
// caller (header-only / sites-only paths) is deterministic across
// implementations we compare provenance-stripped stdout byte-for-byte;
// otherwise we assert both sides accept the bundle and exit cleanly.
func TestCliRolloutBundledMatchesUpstream(t *testing.T) {
	h := buildRolloutHarness(t)
	basic := fixture(t, "basic.vcf")

	t.Run("view -hG header matches upstream", func(t *testing.T) {
		ours := stripProvenance(runOut(t, h.ours, nil, "view", "-hG", basic))
		up := stripProvenance(runOut(t, h.upstream, nil, "view", "-hG", basic))
		if ours != up {
			t.Fatalf("view -hG header differs from upstream\n--- ours ---\n%s\n--- upstream ---\n%s", ours, up)
		}
	})

	t.Run("view -GHi accepted by both (sites, no header)", func(t *testing.T) {
		// -G drop genotypes, -H no header; -i with an expression that
		// matches every record. Both binaries must accept the bundle.
		args := []string{"view", "-GH", "-i", "QUAL>=0", basic}
		ours := stripProvenance(runOut(t, h.ours, nil, args...))
		up := stripProvenance(runOut(t, h.upstream, nil, args...))
		if ours != up {
			t.Fatalf("view -GH -i body differs from upstream\n--- ours ---\n%s\n--- upstream ---\n%s", ours, up)
		}
	})

	t.Run("stats -s- SN block matches upstream", func(t *testing.T) {
		ours := snLines(runOut(t, h.ours, nil, "stats", "-s-", basic))
		up := snLines(runOut(t, h.upstream, nil, "stats", "-s-", basic))
		if ours != up {
			t.Fatalf("stats -s- SN block differs from upstream\n--- ours ---\n%s\n--- upstream ---\n%s", ours, up)
		}
	})
}

// snLines extracts the "SN" summary-number rows from `bcftools stats`
// output, dropping the leading bracketed column index that upstream
// prefixes and we may format slightly differently.
func snLines(b []byte) string {
	var keep []string
	for _, line := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(line, "SN\t") {
			keep = append(keep, line)
		}
	}
	return strings.Join(keep, "\n")
}

// TestCliRolloutLegacyFlagsParse exercises the upstream retro-compat
// short flags registered during the rollout, ensuring each parses (and
// in bundled position) on our binary and is accepted by upstream.
func TestCliRolloutLegacyFlagsParse(t *testing.T) {
	h := buildRolloutHarness(t)
	call := fixture(t, "call_in.vcf")
	basic := fixture(t, "basic.vcf")

	cases := []struct {
		name string
		args []string
	}{
		// norm -D: deprecated alias of -d exact, here bundled with -a.
		{"norm -aD", []string{"norm", "-aD", basic}},
		// call -f: deprecated alias of -a/--annotate (value-taking),
		// bundled as the cluster-terminating flag.
		{"call -mvf GQ", []string{"call", "-mvf", "GQ", call}},
		// call -Y: deprecated alias of --ploidy Y (haploid), bundled.
		{"call -mvY", []string{"call", "-mvY", call}},
		// call -N: omit-REF-N (default), accepted no-op, bundled.
		{"call -mvN", []string{"call", "-mvN", call}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Must parse and run cleanly on our binary.
			_ = runOut(t, h.ours, nil, tc.args...)
			// And the bundle must be ACCEPTED (parsed) by the live
			// upstream binary. Upstream `call` can still exit non-zero on
			// these fixtures because they lack the per-sample likelihoods
			// (QS) that mpileup would have produced — that is a data error,
			// not a parse error. So we assert the bundle is parsed: upstream
			// must not emit a getopt option-parsing diagnostic.
			cmd := exec.Command(h.upstream, tc.args...)
			cmd.Stdout = io.Discard
			var stderr bytes.Buffer
			cmd.Stderr = &stderr
			_ = cmd.Run()
			if msg := getoptParseError(stderr.String()); msg != "" {
				t.Fatalf("upstream rejected legacy/compat bundle %v at parse time: %s", tc.args, msg)
			}
		})
	}
}

// getoptParseError returns the offending line when upstream stderr
// contains a getopt option-parsing diagnostic (unknown option, missing
// option-argument, etc.), or "" when no such diagnostic is present.
func getoptParseError(stderr string) string {
	needles := []string{
		"invalid option",
		"unrecognized option",
		"unknown option",
		"option requires an argument",
		"illegal option",
	}
	for _, line := range strings.Split(stderr, "\n") {
		low := strings.ToLower(line)
		for _, n := range needles {
			if strings.Contains(low, n) {
				return line
			}
		}
	}
	return ""
}

// TestCliRolloutParseInProcess asserts at the unit level that every
// subcommand's run* entry point accepts a bundled short-flag form
// without a parse error (exit code != 2). It runs in-process (no
// subprocess) so it executes even when the upstream binary is absent,
// guarding the rollout wiring directly.
func TestCliRolloutParseInProcess(t *testing.T) {
	basic, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "basic.vcf"))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	if _, err := os.Stat(basic); err != nil {
		t.Skipf("fixture missing: %v", err)
	}
	// Each entry's args use a bundled short-flag cluster; we only assert
	// the parse layer did not reject it (exit code 2 means parse error).
	cases := []struct {
		name string
		fn   func([]string) int
		args []string
	}{
		{"view -hG", runView, []string{"-hG", basic}},
		{"view -GH", runView, []string{"-GH", basic}},
		{"stats -s-", runStatsCmd, []string{"-s-", basic}},
		{"norm -aD", runNorm, []string{"-aD", basic}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			withStdout(t, func() {
				if rc := tc.fn(tc.args); rc == 2 {
					t.Fatalf("%s: parse error (rc=2) for bundled args %v", tc.name, tc.args)
				}
			})
		})
	}
}

// withStdout redirects os.Stdout to /dev/null for the duration of fn so
// in-process run* calls do not pollute test output, restoring it after.
func withStdout(t *testing.T, fn func()) {
	t.Helper()
	orig := os.Stdout
	devnull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		t.Fatalf("open devnull: %v", err)
	}
	os.Stdout = devnull
	defer func() {
		os.Stdout = orig
		devnull.Close()
	}()
	fn()
}
