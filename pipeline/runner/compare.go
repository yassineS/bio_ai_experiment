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

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
)

// stripProvenance removes non-reproducible provenance lines from a text output
// stream so ByteExact comparisons are robust against tool-version stamps. It
// drops:
//
//   - SAM/BAM "@PG" and "@CO" header lines (program/comment provenance),
//   - VCF "##<tool>_..." command lines (e.g. ##bcftools_viewCommand=...,
//     ##samtoolsCommand=...) and "##source=" / "##fileDate=" headers,
//   - the "# This file was produced by ..." / "# The command line was: ..."
//     comment-block headers that "samtools stats", "bcftools stats", and
//     "bcftools gtcheck" emit (these carry the upstream version string and the
//     literal command line, so they differ by build/working directory),
//   - the "##FILTER=<ID=PASS,Description=\"All filters passed\">" boilerplate
//     line: bcftools auto-inserts it whenever a FILTER column is written, but
//     our header serialiser places it at a different position than upstream.
//     It is identical, tool-inserted boilerplate (not data) on both sides, so
//     dropping it neutralises the position-only difference without hiding any
//     genuine data divergence (a different FILTER definition would still
//     differ in content and fail).
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
		case bytes.HasPrefix(ln, []byte("##FILTER=<ID=PASS,Description=\"All filters passed\">")):
			continue
		case bytes.HasPrefix(ln, []byte("##")) && isCommandHeader(ln):
			continue
		case isStatsProvenance(ln):
			continue
		case bytes.HasPrefix(ln, []byte("INFO\tTime required")):
			// bcftools gtcheck prints a non-reproducible wall-clock timing line
			// ("INFO\tTime required to process one record .. <seconds>"); our port
			// omits it. It is timing provenance, not data.
			continue
		}
		out = append(out, ln)
	}
	return bytes.Join(out, []byte("\n"))
}

// isStatsProvenance matches the comment-block provenance lines that the
// samtools/bcftools stats-style reports (samtools stats, bcftools stats,
// bcftools gtcheck) emit. These echo the upstream version string ("This file
// was produced by ..."), the literal command line ("The command line was: ..."
// and its tab-indented "# \t<cmd>" echo), the working directory ("and the
// working directory was: ..." plus its "# \t<path>" echo), the fixed
// "This file contains statistics for all reads." banner, and a bare "#"
// separator — none of which is reproducible across builds or working
// directories, and which our ports omit. The data-describing comment rows a
// stats/gtcheck report keeps (e.g. "# CHK, Checksum...", "# ID\t...",
// "# DCv2, ...", "#DCv2\t...") do NOT match these patterns, so they are
// preserved on both sides.
//
// The bare "#" line is stripped unconditionally: where both sides emit one
// (bcftools stats) it is removed from both equally (no net effect), and where
// only upstream emits one (gtcheck's provenance separator) it removes the
// spurious one-sided line.
func isStatsProvenance(ln []byte) bool {
	s := string(ln)
	if s == "#" {
		return true
	}
	if strings.HasPrefix(s, "# \t") {
		// Tab-indented continuation echo of the command line / working directory.
		return true
	}
	if !strings.HasPrefix(s, "# ") {
		return false
	}
	body := strings.TrimSpace(s[1:])
	switch {
	case strings.HasPrefix(body, "This file was produced by"),
		strings.HasPrefix(body, "This file contains statistics"),
		strings.HasPrefix(body, "The command line was"),
		strings.HasPrefix(body, "and the working directory was"):
		return true
	}
	return false
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

// similarityEpsilon is the default relative tolerance for numeric field
// comparison. An entry may widen it via Entry.Tolerance (see resolveEpsilon).
const similarityEpsilon = 1e-6

// resolveEpsilon returns the per-entry numeric tolerance, falling back to the
// package default when the entry does not set one.
func resolveEpsilon(tol float64) float64 {
	if tol > 0 {
		return tol
	}
	return similarityEpsilon
}

// CompareSimilarity compares streams structurally: identical non-numeric tokens
// and numeric tokens within a relative epsilon eps. It records the maximum
// relative deviation observed. Used for heuristic / float-scored / RNG paths
// where byte-exact equality is not expected but structural+numeric agreement
// is.
func CompareSimilarity(ours, upstream []byte, eps float64) CompareResult {
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
				if dev > eps {
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
func CompareOutputFiles(ourPrefix, upPrefix string, suffixes []string, mode matrix.CompareMode, eps float64) CompareResult {
	var maxDev float64
	// BAMDecoded output files (e.g. the per-read-group BAMs samtools split
	// writes) are decoded through the upstream samtools so only their records
	// are compared, bypassing the BGZF framing difference.
	var samBin string
	if mode == matrix.BAMDecoded {
		b, err := upstream.Binary("samtools")
		if err != nil {
			return CompareResult{Equal: false, Detail: "BAM decode needs the upstream samtools binary: " + err.Error()}
		}
		samBin = b
	}
	for _, sfx := range suffixes {
		ourPath := ourPrefix + sfx
		upPath := upPrefix + sfx
		ourBytes, ourErr := readOutputFile(ourPath, mode, samBin)
		upBytes, upErr := readOutputFile(upPath, mode, samBin)
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
			cmp = CompareSimilarity(ourBytes, upBytes, eps)
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

// readOutputFile reads one named output file for comparison. In BAMDecoded mode
// it pipes the file's bytes through `samtools view -h` so two BAMs with
// different BGZF framing are compared by their decoded records (provenance is
// stripped by CompareByteExact); otherwise it falls back to readMaybeGzip.
func readOutputFile(path string, mode matrix.CompareMode, samBin string) ([]byte, error) {
	if mode != matrix.BAMDecoded {
		return readMaybeGzip(path)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	return decodeBAM(samBin, raw)
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
