package main

// Live-upstream parity tests for BAM binary OUTPUT.
//
// Upstream `bedtools intersect` with a BAM query (-a/-abam) and no -bed writes
// the intersecting ALIGNMENTS back out as BAM by default. These tests drive our
// binary and the live upstream `bedtools intersect` over the same BAM fixtures,
// then compare by decoding BOTH BAM outputs to SAM text (header + records) with
// our own pkg/htsgo/sam reader and diffing byte for byte. BAM is BGZF-compressed
// so a raw byte diff of the .bam files is not meaningful; the SAM projection of
// the decoded stream is the stable comparison surface.
//
// Tests use t.Fatalf (never t.Skip): a missing/un-buildable upstream binary or
// any divergence is a hard failure. The upstream and our binaries are each built
// once via the shared upstreamBedtools / buildOurs helpers.

import (
	"bytes"
	"io"
	"os/exec"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// decodeBAMToSAM decodes a BGZF-wrapped BAM byte stream to SAM text (the header
// followed by one tab-delimited line per record), the stable surface for
// comparing two BAM outputs. An empty input (no bytes at all) decodes to the
// empty string; any decode error is a hard test failure.
func decodeBAMToSAM(t *testing.T, raw []byte) string {
	t.Helper()
	if len(raw) == 0 {
		return ""
	}
	r, err := sam.NewBAMReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("open BAM stream: %v", err)
	}
	var out bytes.Buffer
	w := sam.NewSAMWriter(&out)
	if err := w.WriteHeader(r.Header()); err != nil {
		t.Fatalf("write SAM header: %v", err)
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("read BAM record: %v", err)
		}
		if err := w.Write(rec); err != nil {
			t.Fatalf("write SAM record: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("flush SAM: %v", err)
	}
	return out.String()
}

// TestUpstreamBAMOutput_Intersect asserts that a BAM query without -bed produces
// BAM output identical (after decode to SAM) to upstream's, across the default
// mode and the alignment-level flags that BAM output honours (-u/-v/-wa), plus
// strand/fraction filters and the split/unmapped fixtures.
func TestUpstreamBAMOutput_Intersect(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []fixtureCase{
		// Default BAM output: each A alignment with >=1 overlap, once.
		{"default", "a.bam", []string{"a.bed"}, nil},
		{"default_via_abam", "a.bam", []string{"a.bed"}, nil}, // -abam alias exercised below
		{"unique", "a.bam", []string{"a.bed"}, []string{"-u"}},
		{"invert", "a.bam", []string{"a.bed"}, []string{"-v"}},
		{"writeA", "a.bam", []string{"a.bed"}, []string{"-wa"}},

		// -C is NOT gated for BAM-A: it falls through to default BAM output.
		{"countEach", "a.bam", []string{"a.bed"}, []string{"-C"}},

		// -wb / -loj are ignored (warn) and still produce default BAM output.
		{"writeB_ignored", "a.bam", []string{"a.bed"}, []string{"-wb"}},
		{"leftJoin_ignored", "a.bam", []string{"a.bed"}, []string{"-loj"}},

		// Strand and fraction filters apply to the overlap test.
		{"strand", "a.bam", []string{"a.bed"}, []string{"-s"}},
		{"opp_strand", "a.bam", []string{"a.bed"}, []string{"-S"}},
		{"fraction", "a.bam", []string{"a.bed"}, []string{"-f", "0.5"}},

		// Unmapped reads: absent under default, reported under -v.
		{"unmapped_default", "a_with_bothUnmapped.bam", []string{"a.bed"}, nil},
		{"unmapped_invert", "a_with_bothUnmapped.bam", []string{"a.bed"}, []string{"-v"}},
		{"oneUnmapped_default", "oneUnmapped.bam", []string{"j1.bed"}, nil},
		{"oneUnmapped_invert", "oneUnmapped.bam", []string{"j1.bed"}, []string{"-v"}},

		// Split (spliced) alignments: -split block-aware overlap, BAM output.
		{"split_match", "three_blocks_match.bam", []string{"three_blocks_match.bed"}, []string{"-split"}},
		{"split_match_s", "three_blocks_match.bam", []string{"three_blocks_match.bed"}, []string{"-split", "-s"}},

		// Fraction-boundary fixtures over the x/y design (BAM query, -wa).
		{"xy_f_F_wa", "x.bam", []string{"y.bed"}, []string{"-f", "0.21", "-F", "0.21", "-wa"}},
		{"xy_f_F2_wa", "x.bam", []string{"y.bed"}, []string{"-f", "0.19", "-F", "0.50", "-wa"}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			aFlag := "-a"
			if strings.HasSuffix(tc.name, "_via_abam") {
				aFlag = "-abam"
			}
			upArgs := append([]string{"intersect", aFlag, tc.a, "-b"}, tc.b...)
			upArgs = append(upArgs, tc.args...)
			ourArgs := append([]string{aFlag, tc.a, "-b"}, tc.b...)
			ourArgs = append(ourArgs, tc.args...)

			upOut, _ := runOut(t, dir, bt, upArgs...)
			ourOut, _ := runOut(t, dir, ours, ourArgs...)

			upSAM := decodeBAMToSAM(t, upOut)
			ourSAM := decodeBAMToSAM(t, ourOut)
			if upSAM != ourSAM {
				t.Fatalf("BAM-output (decoded SAM) mismatch for %s (args %v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.name, tc.args, upSAM, ourSAM)
			}
		})
	}
}

// TestUpstreamBAMOutput_RequiresBedErrors asserts that the flags upstream rejects
// for a BAM query without -bed (-c, -wo, -wao) produce the same ERROR banner on
// stderr and a non-zero exit in both binaries. Upstream prints its full help
// ahead of the banner; we compare the final stderr line (the banner) and the
// exit status, which is the load-bearing, stable part of the diagnostic.
func TestUpstreamBAMOutput_RequiresBedErrors(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []struct {
		name string
		args []string
	}{
		{"count", []string{"-c"}},
		{"writeOverlap", []string{"-wo"}},
		{"writeAllOverlap", []string{"-wao"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			upArgs := append([]string{"intersect", "-abam", "a.bam", "-b", "a.bed"}, tc.args...)
			ourArgs := append([]string{"-abam", "a.bam", "-b", "a.bed"}, tc.args...)

			upStdout, upErr := runOutCode(t, dir, bt, upArgs...)
			ourStdout, ourErr := runOutCode(t, dir, ours, ourArgs...)

			if upStdout.code == 0 || ourStdout.code == 0 {
				t.Fatalf("%s: expected non-zero exit (upstream=%d ours=%d)",
					tc.name, upStdout.code, ourStdout.code)
			}
			upBanner := lastLine(upErr)
			ourBanner := lastLine(ourErr)
			if upBanner != ourBanner {
				t.Fatalf("%s: ERROR banner mismatch\n--- upstream ---\n%q\n--- ours ---\n%q",
					tc.name, upBanner, ourBanner)
			}
			if !strings.Contains(ourBanner, "is not valid with BAM query input, unless bed output is specified with -bed option.") {
				t.Fatalf("%s: unexpected banner %q", tc.name, ourBanner)
			}
		})
	}
}

// TestUpstreamBAMOutput_IgnoredWarnings asserts the warn-and-ignore stderr
// diagnostics (-wb/-loj and -header) match upstream byte for byte for a BAM
// query without -bed.
func TestUpstreamBAMOutput_IgnoredWarnings(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []struct {
		name string
		args []string
	}{
		{"writeB", []string{"-wb"}},
		{"leftJoin", []string{"-loj"}},
		{"header", []string{"-header"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			upArgs := append([]string{"intersect", "-abam", "a.bam", "-b", "a.bed"}, tc.args...)
			ourArgs := append([]string{"-abam", "a.bam", "-b", "a.bed"}, tc.args...)
			_, upErr := runOut(t, dir, bt, upArgs...)
			_, ourErr := runOut(t, dir, ours, ourArgs...)
			if !bytes.Equal(upErr, ourErr) {
				t.Fatalf("%s: stderr mismatch\n--- upstream ---\n%q\n--- ours ---\n%q",
					tc.name, upErr, ourErr)
			}
		})
	}
}

// runResult bundles a captured stdout/stderr with the process exit code.
type runResult struct {
	out  []byte
	code int
}

// runOutCode runs a command in dir, returning its stdout (with exit code) and
// stderr separately. Unlike runOut it captures the exit code so error-path tests
// can assert a non-zero status.
func runOutCode(t *testing.T, dir, name string, args ...string) (stdout runResult, stderr []byte) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	var so, se bytes.Buffer
	cmd.Stdout = &so
	cmd.Stderr = &se
	code := 0
	if err := cmd.Run(); err != nil {
		if ee, ok := err.(*exec.ExitError); ok {
			code = ee.ExitCode()
		} else {
			t.Fatalf("run %s %v: %v", name, args, err)
		}
	}
	return runResult{out: so.Bytes(), code: code}, se.Bytes()
}

// lastLine returns the final non-empty line of b (trailing whitespace trimmed),
// or "" when there is none. Upstream prints its help ahead of the ERROR banner,
// so the banner is the last line.
func lastLine(b []byte) string {
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	for i := len(lines) - 1; i >= 0; i-- {
		if strings.TrimSpace(lines[i]) != "" {
			return strings.TrimSpace(lines[i])
		}
	}
	return ""
}
