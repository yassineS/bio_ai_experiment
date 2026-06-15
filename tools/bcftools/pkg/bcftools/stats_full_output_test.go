package bcftools

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// stripStatsProvenance removes the lines that legitimately differ between the
// upstream C binary and the Go port: the two provenance lines (which embed the
// upstream version string and the literal argv) and every `##` header meta
// line (which the Go port does not echo). Everything that remains must match
// byte-for-byte.
func stripStatsProvenance(s string) string {
	lines := strings.Split(s, "\n")
	out := lines[:0]
	for _, l := range lines {
		switch {
		case strings.HasPrefix(l, "# This file was produced by"):
			continue
		case strings.HasPrefix(l, "# The command line was"):
			continue
		case strings.HasPrefix(l, "##"):
			continue
		}
		out = append(out, l)
	}
	return strings.Join(out, "\n")
}

// TestStatsFullOutputParity asserts that the Go port's complete `bcftools
// stats` output — every section, comment block, legend, and data row — matches
// the live upstream binary byte-for-byte on representative fixtures, once the
// provenance/`##` lines are stripped. It exercises both the per-sample mode
// (`-s -`, which enables PSC/PSI/HWE) and the sites-only mode (no `-s`).
func TestStatsFullOutputParity(t *testing.T) {
	bin := upstreamBcftools(t)
	root := mustRepoRoot()
	parity := filepath.Join(root, "tools", "bcftools", "testdata", "parity")

	fixtures := []string{"basic.vcf", "multi.vcf", "biallelic.vcf", "atom.vcf"}
	modes := []struct {
		name string
		args []string
	}{
		{name: "samples", args: []string{"stats", "-s", "-"}},
		{name: "sites", args: []string{"stats"}},
	}

	for _, fx := range fixtures {
		path := filepath.Join(parity, fx)
		for _, mode := range modes {
			t.Run(fx+"/"+mode.name, func(t *testing.T) {
				upArgs := append(append([]string{}, mode.args...), path)
				upOut, err := exec.Command(bin, upArgs...).Output()
				if err != nil {
					t.Fatalf("upstream bcftools %v: %v", upArgs, err)
				}

				var goBuf bytes.Buffer
				opts := StatsOptions{InputFile: path}
				if mode.name == "samples" {
					opts.Samples = []string{"-"}
				}
				if _, err := StatsFile(path, &goBuf, opts); err != nil {
					t.Fatalf("StatsFile(%s): %v", path, err)
				}

				want := stripStatsProvenance(string(upOut))
				got := stripStatsProvenance(goBuf.String())
				if want != got {
					t.Fatalf("full stats output mismatch for %s (%s):\n--- upstream ---\n%s\n--- go ---\n%s",
						fx, mode.name, want, got)
				}
			})
		}
	}
}
