package runner

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
)

// stripProvenance removes non-reproducible provenance lines from a text output
// stream so ByteExact comparisons are robust against tool-version stamps. It
// drops:
//
//   - SAM/BAM "@PG" and "@CO" header lines (program/comment provenance),
//   - VCF "##<tool>_..." command lines (e.g. ##bcftools_viewCommand=...,
//     ##samtoolsCommand=...) and "##source=" / "##fileDate=" headers.
//
// Data lines and structural headers (@SQ, ##contig, ##INFO, the column header)
// are preserved, so a genuine output difference still fails the comparison.
func stripProvenance(b []byte) []byte {
	lines := bytes.Split(b, []byte("\n"))
	out := lines[:0]
	for _, ln := range lines {
		switch {
		case bytes.HasPrefix(ln, []byte("@PG")), bytes.HasPrefix(ln, []byte("@CO")):
			continue
		case bytes.HasPrefix(ln, []byte("##source=")),
			bytes.HasPrefix(ln, []byte("##fileDate=")),
			bytes.HasPrefix(ln, []byte("##reference=")):
			continue
		case bytes.HasPrefix(ln, []byte("##")) && isCommandHeader(ln):
			continue
		}
		out = append(out, ln)
	}
	return bytes.Join(out, []byte("\n"))
}

// isCommandHeader matches VCF provenance command lines like
// "##bcftools_viewCommand=..." or "##samtoolsCommand=...".
func isCommandHeader(ln []byte) bool {
	s := string(ln)
	if !strings.HasPrefix(s, "##") {
		return false
	}
	eq := strings.IndexByte(s, '=')
	if eq < 0 {
		return false
	}
	key := s[2:eq]
	return strings.Contains(key, "Command") || strings.Contains(key, "_command") ||
		strings.Contains(strings.ToLower(key), "version")
}

// CompareResult holds the outcome of comparing two output streams.
type CompareResult struct {
	Equal        bool
	MaxDeviation float64 // for Similarity: largest numeric field deviation seen
	Detail       string  // human-readable explanation on mismatch
}

// CompareByteExact compares provenance-stripped streams for exact equality.
func CompareByteExact(ours, upstream []byte) CompareResult {
	a := stripProvenance(ours)
	b := stripProvenance(upstream)
	if bytes.Equal(a, b) {
		return CompareResult{Equal: true}
	}
	return CompareResult{Equal: false, Detail: firstDiff(a, b)}
}

// similarityEpsilon is the relative tolerance for numeric field comparison.
const similarityEpsilon = 1e-6

// CompareSimilarity compares streams structurally: identical non-numeric tokens
// and numeric tokens within a relative epsilon. It records the maximum relative
// deviation observed. Used for heuristic / float-scored / RNG paths where
// byte-exact equality is not expected but structural+numeric agreement is.
func CompareSimilarity(ours, upstream []byte) CompareResult {
	al := splitLines(stripProvenance(ours))
	bl := splitLines(stripProvenance(upstream))
	if len(al) != len(bl) {
		return CompareResult{Equal: false, Detail: fmt.Sprintf("line count differs: ours=%d upstream=%d", len(al), len(bl))}
	}
	var maxDev float64
	for i := range al {
		at := strings.Fields(al[i])
		bt := strings.Fields(bl[i])
		if len(at) != len(bt) {
			return CompareResult{Equal: false, MaxDeviation: maxDev,
				Detail: fmt.Sprintf("line %d field count differs", i+1)}
		}
		for j := range at {
			af, aok := parseNum(at[j])
			bf, bok := parseNum(bt[j])
			if aok && bok {
				dev := relDev(af, bf)
				if dev > maxDev {
					maxDev = dev
				}
				if dev > similarityEpsilon {
					return CompareResult{Equal: false, MaxDeviation: dev,
						Detail: fmt.Sprintf("line %d field %d numeric deviation %.3g (%v vs %v)", i+1, j+1, dev, af, bf)}
				}
				continue
			}
			if at[j] != bt[j] {
				return CompareResult{Equal: false, MaxDeviation: maxDev,
					Detail: fmt.Sprintf("line %d field %d differs: %q vs %q", i+1, j+1, at[j], bt[j])}
			}
		}
	}
	return CompareResult{Equal: true, MaxDeviation: maxDev}
}

func parseNum(s string) (float64, bool) {
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, false
	}
	return f, true
}

func relDev(a, b float64) float64 {
	if a == b {
		return 0
	}
	denom := math.Max(math.Abs(a), math.Abs(b))
	if denom == 0 {
		return 0
	}
	return math.Abs(a-b) / denom
}

func splitLines(b []byte) []string {
	s := strings.TrimRight(string(b), "\n")
	if s == "" {
		return nil
	}
	return strings.Split(s, "\n")
}

// firstDiff returns a short description of the first differing line.
func firstDiff(a, b []byte) string {
	al := splitLines(a)
	bl := splitLines(b)
	n := len(al)
	if len(bl) < n {
		n = len(bl)
	}
	for i := 0; i < n; i++ {
		if al[i] != bl[i] {
			return fmt.Sprintf("first diff at line %d:\n  ours:     %s\n  upstream: %s", i+1, trunc(al[i]), trunc(bl[i]))
		}
	}
	if len(al) != len(bl) {
		return fmt.Sprintf("line count differs: ours=%d upstream=%d", len(al), len(bl))
	}
	return "streams differ"
}

func trunc(s string) string {
	if len(s) > 200 {
		return s[:200] + "..."
	}
	return s
}

// CompareOutputFiles compares the named output files written under two prefixes
// (e.g. "<ourdir>/out" vs "<updir>/out"). suffixes are the per-file extensions
// relative to the prefix (e.g. ".frq", ".mosdepth.summary.txt"); each pair is
// read (decompressing transparently when the suffix ends in ".gz") and compared
// with the chosen mode. The first mismatching file fails the whole entry; a
// missing-on-exactly-one-side file is a divergence. This is how the vcftools
// and mosdepth matrices verify multi-file output.
func CompareOutputFiles(ourPrefix, upPrefix string, suffixes []string, mode matrix.CompareMode) CompareResult {
	var maxDev float64
	for _, sfx := range suffixes {
		ourPath := ourPrefix + sfx
		upPath := upPrefix + sfx
		ourBytes, ourErr := readMaybeGzip(ourPath)
		upBytes, upErr := readMaybeGzip(upPath)
		// Presence mismatch is a real divergence (one side wrote a file the
		// other did not).
		if (ourErr == nil) != (upErr == nil) {
			return CompareResult{Equal: false, MaxDeviation: maxDev,
				Detail: fmt.Sprintf("output file %q presence differs: ours_err=%v upstream_err=%v", sfx, ourErr, upErr)}
		}
		if ourErr != nil && upErr != nil {
			// Neither side produced the file: nothing to compare for this suffix.
			continue
		}
		var cmp CompareResult
		if mode == matrix.Similarity {
			cmp = CompareSimilarity(ourBytes, upBytes)
		} else {
			cmp = CompareByteExact(ourBytes, upBytes)
		}
		if cmp.MaxDeviation > maxDev {
			maxDev = cmp.MaxDeviation
		}
		if !cmp.Equal {
			return CompareResult{Equal: false, MaxDeviation: maxDev,
				Detail: fmt.Sprintf("output file %q: %s", sfx, cmp.Detail)}
		}
	}
	return CompareResult{Equal: true, MaxDeviation: maxDev}
}

// readMaybeGzip reads a file, transparently decompressing it when the path ends
// in ".gz" (BGZF is gzip-compatible). Comparing decompressed payloads isolates
// the data content from BGZF block-framing differences between our deflate
// backend and upstream's.
func readMaybeGzip(path string) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	if !strings.HasSuffix(path, ".gz") {
		return io.ReadAll(f)
	}
	gr, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer gr.Close()
	gr.Multistream(true)
	return io.ReadAll(gr)
}
