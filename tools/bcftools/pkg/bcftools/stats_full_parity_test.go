package bcftools

// Whole-output parity tests for `bcftools stats`. Unlike the per-section
// tests (which assert one section at a time), these run the *entire*
// `bcftools stats -s -` report from the live upstream binary and compare it,
// byte-for-byte, against our port's full output — after stripping only the
// non-deterministic provenance lines (the version-stamped "# This file was
// produced by:" / "# The command line was:" header pair and any "##" meta
// lines). They are the regression guard for the TSTV, SiS, VAF and verbose
// comment-block work that closed the full-output gap.

import (
	"bytes"
	"os/exec"
	"strings"
	"testing"
)

// stripStatsProvenance drops the lines that legitimately differ between the
// upstream binary and our port: the two-line provenance header and any "##"
// meta lines. Everything else must match exactly.
func stripStatsProvenance(b []byte) string {
	var keep []string
	for _, line := range strings.Split(string(b), "\n") {
		switch {
		case strings.HasPrefix(line, "# This file was produced by"):
			continue
		case strings.HasPrefix(line, "# The command line was"):
			continue
		case strings.HasPrefix(line, "##"):
			continue
		}
		keep = append(keep, line)
	}
	return strings.Join(keep, "\n")
}

// TestStatsFullOutputParity asserts that `bcftools stats -s - FIXTURE`
// produces output identical to the live upstream binary (modulo provenance)
// across a SNP/indel/multi-allelic site fixture, a FORMAT/AD fixture that
// exercises the VAF section, and a multi-sample fixture.
func TestStatsFullOutputParity(t *testing.T) {
	upstream := upstreamBcftools(t)
	fixtures := []string{
		"basic.vcf",      // SNPs, an indel, a multi-allelic site, no FORMAT/AD
		"gt_plugins.vcf", // 4 samples with FORMAT/AD -> exercises the VAF section
		"multi.vcf",      // multi-sample
	}
	for _, fx := range fixtures {
		fx := fx
		t.Run(fx, func(t *testing.T) {
			path := parityPath(t, fx)

			out, err := exec.Command(upstream, "stats", "-s", "-", path).Output()
			if err != nil {
				t.Fatalf("upstream bcftools stats %s: %v", fx, err)
			}
			want := stripStatsProvenance(out)

			var got bytes.Buffer
			opts := StatsOptions{Samples: []string{"-"}, EnableSamples: true, InputFile: path}
			if _, err := StatsFile(path, &got, opts); err != nil {
				t.Fatalf("StatsFile %s: %v", fx, err)
			}
			gotStr := stripStatsProvenance(got.Bytes())

			if gotStr != want {
				t.Fatalf("full stats output mismatch for %s:\n--- upstream ---\n%s\n--- ours ---\n%s",
					fx, want, gotStr)
			}
		})
	}
}
