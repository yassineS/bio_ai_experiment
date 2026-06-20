package main

// Live-upstream parity tests for CRAM binary OUTPUT.
//
// Upstream `bedtools intersect` with a CRAM query (-a) and no -bed writes the
// intersecting ALIGNMENTS back out as CRAM — but ONLY when a CRAM reference is
// available (via the global --cram-ref flag or the CRAM_REFERENCE environment
// variable). Its writer opens htslib mode "wc" when a reference is set and "wb"
// (BAM) otherwise (src/utils/BamTools/include/BamWriter.hpp). So:
//
//   - CRAM query, reference set  -> CRAM output
//   - CRAM query, no reference   -> BAM output
//   - CRAM query, -ubam          -> BAM output (uncompressed; we emit compressed)
//
// These tests drive our binary and the live upstream `bedtools intersect` over
// the same CRAM fixture, then compare by decoding BOTH CRAM outputs back to SAM
// (with the same reference) using our pkg/htsgo/cram reader. CRAM bytes are not
// identical across encoders (block layout / codec choices differ), so the SAM
// projection of the decoded stream is the stable comparison surface. Record
// lines (RNAME/POS/CIGAR/FLAG/MAPQ/SEQ/QUAL/aux) are compared byte for byte; the
// header is compared as a line set (our CRAM header round-trip preserves every
// @-line but may reorder @SQ/@RG/@PG relative to upstream, which is cosmetic).
//
// Tests use t.Fatalf (never t.Skip): a missing/un-buildable upstream binary or
// any record divergence is a hard failure.

import (
	"bytes"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// runOutEnv runs a command in dir with an explicit environment and returns its
// stdout, never failing on a non-zero exit. It is runOut with a custom env so
// the CRAM_REFERENCE-driven cases can be exercised.
func runOutEnv(t *testing.T, dir string, env []string, name string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	cmd.Env = env
	var so, se bytes.Buffer
	cmd.Stdout = &so
	cmd.Stderr = &se
	_ = cmd.Run()
	return so.Bytes()
}

// cramFixtureRef is the reference FASTA bedtools ships for its CRAM fixtures
// (a.cram / b.cram), relative to the intersect fixtures directory. The CRAM
// reads are reference-compressed against it, so it is required both to decode
// the query and to verify the re-emitted CRAM output.
const cramFixtureRef = "test_ref.fa"

// decodeCRAMToSAM decodes a CRAM byte stream (reconstructed against refPath) to
// a SAM projection: the header line set and one tab-delimited line per record.
// It is the stable surface for comparing two CRAM outputs whose raw bytes differ
// by encoder choice. An empty input decodes to two empty strings.
func decodeCRAMToSAM(t *testing.T, raw []byte, refPath string) (header []string, records string) {
	t.Helper()
	if len(raw) == 0 {
		return nil, ""
	}
	rr, err := cram.NewRecordReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("open CRAM stream: %v", err)
	}
	if refPath != "" {
		if err := rr.SetReferenceFASTA(refPath); err != nil {
			t.Fatalf("attach CRAM reference %s: %v", refPath, err)
		}
	}
	rr.UseRefCacheFromEnv()

	var full bytes.Buffer
	w := sam.NewSAMWriter(&full)
	if err := w.WriteHeader(rr.Header()); err != nil {
		t.Fatalf("write SAM header: %v", err)
	}
	for {
		rec, err := rr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("read CRAM record: %v", err)
		}
		if err := w.Write(rec); err != nil {
			t.Fatalf("write SAM record: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("flush SAM: %v", err)
	}
	return splitSAM(full.String())
}

// splitSAM partitions SAM text into its sorted header (@-prefixed) line set and
// the joined record block. The header is sorted so a cosmetic reordering of
// @SQ/@RG/@PG between two equivalent files does not register as a difference.
func splitSAM(samText string) (header []string, records string) {
	var recLines []string
	for _, line := range strings.Split(strings.TrimRight(samText, "\n"), "\n") {
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, "@") {
			header = append(header, line)
		} else {
			recLines = append(recLines, line)
		}
	}
	sort.Strings(header)
	return header, strings.Join(recLines, "\n")
}

// firstByteMagic returns up to the first four bytes of b, the file-type magic
// (CRAM = "CRAM", BAM = BGZF "\x1f\x8b").
func firstByteMagic(b []byte) []byte {
	if len(b) < 4 {
		return b
	}
	return b[:4]
}

// TestUpstreamCRAMOutput_Intersect asserts that a CRAM query without -bed, with
// a CRAM reference supplied (via --cram-ref and via the CRAM_REFERENCE env var),
// produces CRAM output whose decoded records match upstream's across the default
// mode and the alignment-level flags CRAM output honours (-u/-v/-wa).
func TestUpstreamCRAMOutput_Intersect(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)
	refAbs := filepath.Join(dir, cramFixtureRef)
	if _, err := os.Stat(refAbs); err != nil {
		t.Skipf("CRAM fixture reference missing (%v); run: git submodule update --init reference_code/bedtools", err)
	}

	cases := []struct {
		name string
		args []string
	}{
		{"default", nil},
		{"unique", []string{"-u"}},
		{"invert", []string{"-v"}},
		{"writeA", []string{"-wa"}},
		// -C is not gated for a CRAM/BAM query: it falls through to default output.
		{"countEach", []string{"-C"}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			upArgs := append([]string{"--cram-ref", cramFixtureRef, "intersect", "-a", "a.cram", "-b", "b.cram"}, tc.args...)
			ourArgs := append([]string{"--cram-ref", cramFixtureRef, "-a", "a.cram", "-b", "b.cram"}, tc.args...)

			upOut, _ := runOut(t, dir, bt, upArgs...)
			ourOut, _ := runOut(t, dir, ours, ourArgs...)

			// Both must be CRAM-framed.
			if got := string(firstByteMagic(upOut)); got != "CRAM" {
				t.Fatalf("upstream output is not CRAM (magic %q)", got)
			}
			if got := string(firstByteMagic(ourOut)); got != "CRAM" {
				t.Fatalf("our output is not CRAM (magic %q)", got)
			}

			upHdr, upRecs := decodeCRAMToSAM(t, upOut, refAbs)
			ourHdr, ourRecs := decodeCRAMToSAM(t, ourOut, refAbs)
			if upRecs != ourRecs {
				t.Fatalf("CRAM-output (decoded SAM) record mismatch for %s (args %v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.name, tc.args, upRecs, ourRecs)
			}
			if strings.Join(upHdr, "\n") != strings.Join(ourHdr, "\n") {
				t.Fatalf("CRAM-output header line-set mismatch for %s\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.name, strings.Join(upHdr, "\n"), strings.Join(ourHdr, "\n"))
			}
		})
	}
}

// TestUpstreamCRAMOutput_EnvReference asserts the CRAM_REFERENCE environment
// variable selects CRAM output exactly like the --cram-ref flag, matching
// upstream (which reads the same env var).
func TestUpstreamCRAMOutput_EnvReference(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)
	refAbs := filepath.Join(dir, cramFixtureRef)

	env := append(os.Environ(), "CRAM_REFERENCE="+cramFixtureRef)
	upOut := runOutEnv(t, dir, env, bt, "intersect", "-a", "a.cram", "-b", "b.cram")
	ourOut := runOutEnv(t, dir, env, ours, "-a", "a.cram", "-b", "b.cram")

	if got := string(firstByteMagic(upOut)); got != "CRAM" {
		t.Fatalf("upstream output is not CRAM under CRAM_REFERENCE (magic %q)", got)
	}
	if got := string(firstByteMagic(ourOut)); got != "CRAM" {
		t.Fatalf("our output is not CRAM under CRAM_REFERENCE (magic %q)", got)
	}
	_, upRecs := decodeCRAMToSAM(t, upOut, refAbs)
	_, ourRecs := decodeCRAMToSAM(t, ourOut, refAbs)
	if upRecs != ourRecs {
		t.Fatalf("CRAM_REFERENCE record mismatch\n--- upstream ---\n%s\n--- ours ---\n%s", upRecs, ourRecs)
	}
}

// TestUpstreamCRAMOutput_FramingSelection asserts the output FRAMING matches
// upstream's gating: a CRAM query with no reference writes BAM, while a CRAM
// query with a reference writes CRAM — even under -ubam, because upstream's
// format choice follows the reference alone (the -ubam compression hook is a
// no-op and does not change the format). The first-byte magic is the surface.
func TestUpstreamCRAMOutput_FramingSelection(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []struct {
		name      string
		args      []string
		wantMagic string // "CRAM" or BGZF magic for BAM
	}{
		{"no_reference_is_bam", []string{"intersect", "-a", "a.cram", "-b", "b.cram"}, "\x1f\x8b\x08\x04"},
		{"ubam_with_ref_still_cram", []string{"--cram-ref", cramFixtureRef, "intersect", "-a", "a.cram", "-b", "b.cram", "-ubam"}, "CRAM"},
		{"reference_is_cram", []string{"--cram-ref", cramFixtureRef, "intersect", "-a", "a.cram", "-b", "b.cram"}, "CRAM"},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			// Strip the leading "intersect" subcommand for our binary's args.
			ourArgs := stripIntersect(tc.args)
			upOut, _ := runOut(t, dir, bt, tc.args...)
			ourOut, _ := runOut(t, dir, ours, ourArgs...)

			upMagic := string(firstByteMagic(upOut))
			ourMagic := string(firstByteMagic(ourOut))
			if upMagic != tc.wantMagic {
				t.Fatalf("%s: upstream magic %q, want %q", tc.name, upMagic, tc.wantMagic)
			}
			if ourMagic != tc.wantMagic {
				t.Fatalf("%s: our magic %q, want %q", tc.name, ourMagic, tc.wantMagic)
			}
		})
	}
}

// stripIntersect removes a leading "--cram-ref VALUE" pair (kept) and the
// "intersect" subcommand token from an upstream argv so the same case can drive
// our subcommand-less binary. The global --cram-ref flag precedes "intersect"
// for upstream but is an ordinary flag for us, so it is retained in place.
func stripIntersect(args []string) []string {
	out := make([]string, 0, len(args))
	for _, a := range args {
		if a == "intersect" {
			continue
		}
		out = append(out, a)
	}
	return out
}
