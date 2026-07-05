package runner

import (
	"bytes"
	"compress/gzip"
	"crypto/md5"
	"fmt"
	"hash"
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
		if isProvenanceLine(ln) {
			continue
		}
		out = append(out, normaliseLine(ln))
	}
	return bytes.Join(out, []byte("\n"))
}

// normaliseLine rewrites a kept line to strip non-reproducible FIELDS that live
// inside an otherwise-data-bearing line (as opposed to isProvenanceLine, which
// drops whole lines). It is the SINGLE source of truth for the per-line rewrite,
// shared by the batch stripProvenance and the streaming provenanceFilter so the
// two paths can never drift.
//
// Today it strips the "UR:" tag from SAM "@SQ" header lines (as emitted by
// "samtools dict" and "samtools view -H"). UR is the reference URI —
// "UR:file://<absolute-path>" — which is inherently machine- and
// working-directory-dependent (and its file:// encoding can differ between our
// port and upstream for the same logical path), so it is provenance, not data.
// Every other @SQ field (SN, LN, M5, AN, AS, SP) is preserved, so a genuine
// dict/header divergence (a wrong checksum, length, or field order) still fails
// the comparison. Non-@SQ lines are returned unchanged.
func normaliseLine(ln []byte) []byte {
	if !bytes.HasPrefix(ln, []byte("@SQ\t")) {
		return ln
	}
	fields := bytes.Split(ln, []byte("\t"))
	kept := fields[:0]
	dropped := false
	for _, f := range fields {
		if bytes.HasPrefix(f, []byte("UR:")) {
			dropped = true
			continue
		}
		kept = append(kept, f)
	}
	if !dropped {
		return ln
	}
	return bytes.Join(kept, []byte("\t"))
}

// isProvenanceLine reports whether a single line is a non-reproducible
// provenance line that stripProvenance drops. It is the SINGLE source of truth
// for the per-line keep/drop decision, shared by the batch stripProvenance and
// the streaming provenanceFilter so the two paths can never drift: any line for
// which this returns true is removed identically by both. The patterns it
// matches are documented on stripProvenance.
func isProvenanceLine(ln []byte) bool {
	switch {
	case bytes.HasPrefix(ln, []byte("@PG")), bytes.HasPrefix(ln, []byte("@CO")):
		return true
	case bytes.HasPrefix(ln, []byte("##source=")),
		bytes.HasPrefix(ln, []byte("##fileDate=")),
		bytes.HasPrefix(ln, []byte("##reference=")):
		return true
	case bytes.HasPrefix(ln, []byte("##FILTER=<ID=PASS,Description=\"All filters passed\">")):
		return true
	case bytes.HasPrefix(ln, []byte("##")) && isCommandHeader(ln):
		return true
	case isStatsProvenance(ln):
		return true
	case bytes.HasPrefix(ln, []byte("INFO\tTime required")):
		// bcftools gtcheck prints a non-reproducible wall-clock timing line
		// ("INFO\tTime required to process one record .. <seconds>"); our port
		// omits it. It is timing provenance, not data.
		return true
	}
	return false
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

// StripProvenance removes non-reproducible provenance lines (tool-version
// stamps, @PG/@CO headers, ##<tool>_*Command= lines, stats banner comments,
// timing lines) from a text output stream, returning the normalized bytes.
//
// It is the exact normalization CompareByteExact and CompareSimilarity apply
// before comparing, exported so other harnesses (e.g. the differential fuzzer
// in pipeline/difffuzz) can reuse the SAME notion of "benign provenance
// difference" rather than re-deriving it. A divergence on StripProvenance'd
// bytes is therefore a genuine behavioral difference, not a known-benign stamp.
func StripProvenance(b []byte) []byte { return stripProvenance(b) }

// provenanceFilter is the STREAMING equivalent of stripProvenance: it forwards
// provenance-stripped bytes to an inner io.Writer while buffering at most one
// partial (not-yet-terminated) line, so its memory is O(longest line) rather
// than O(total output). The bytes it emits are byte-for-byte identical to
// bytes.Join(keptLines, "\n") where keptLines are the kept tokens of
// bytes.Split(input, "\n") — the exact normalization CompareByteExact applies.
//
// The trailing-newline / final-empty-token semantics of bytes.Split are
// reproduced precisely: a separator '\n' between two kept tokens is emitted, the
// final token (the bytes after the last '\n', possibly empty) is emitted with no
// trailing '\n', and a dropped token consumes its trailing '\n' too. To match
// bytes.Join, the '\n' separator is written BEFORE a kept token rather than
// after it (a "pending separator" model), which also yields the empty result
// for empty input and for input that is only provenance.
type provenanceFilter struct {
	inner   io.Writer
	line    []byte // bytes of the current (not-yet-terminated) line
	emitted bool   // whether any kept token has already been written
	needSep bool   // a '\n' separator is owed before the next kept token
	err     error
}

// newProvenanceFilter returns a provenanceFilter forwarding kept bytes to inner.
func newProvenanceFilter(inner io.Writer) *provenanceFilter {
	return &provenanceFilter{inner: inner}
}

// Write consumes p, splitting it on '\n'. Each complete line (terminated by a
// '\n' within the stream) is run through isProvenanceLine and, if kept, written
// to the inner writer preceded by any owed separator. Bytes after the final
// '\n' are retained as the partial current line for the next Write/Close. It
// always reports len(p) consumed unless the inner writer errors.
func (f *provenanceFilter) Write(p []byte) (int, error) {
	if f.err != nil {
		return 0, f.err
	}
	total := len(p)
	for {
		i := bytes.IndexByte(p, '\n')
		if i < 0 {
			f.line = append(f.line, p...)
			return total, nil
		}
		f.line = append(f.line, p[:i]...)
		if err := f.flushLine(true); err != nil {
			return total - len(p), err
		}
		p = p[i+1:]
	}
}

// flushLine emits the buffered current line as one completed token. terminated
// indicates the token was followed by a '\n' in the input (so a separator is
// owed before the NEXT kept token). A kept token is written preceded by any
// pending separator; a dropped token writes nothing but still owes a separator
// when it was terminated, exactly mirroring bytes.Join over the kept subset.
func (f *provenanceFilter) flushLine(terminated bool) error {
	keep := !isProvenanceLine(f.line)
	if keep {
		if f.needSep {
			if err := f.writeInner([]byte("\n")); err != nil {
				return err
			}
		}
		if err := f.writeInner(normaliseLine(f.line)); err != nil {
			return err
		}
		f.emitted = true
		f.needSep = false
	}
	if terminated {
		// Every '\n' in the original input is a separator between two tokens.
		// The next kept token (if any) must be preceded by exactly one '\n' per
		// original separator that followed a kept token; a separator that
		// followed a dropped token still owes a '\n' only if a kept token was
		// already emitted (bytes.Join puts separators BETWEEN kept tokens).
		if f.emitted {
			f.needSep = true
		}
	}
	f.line = f.line[:0]
	return nil
}

// writeInner forwards b to the inner writer, recording the first error.
func (f *provenanceFilter) writeInner(b []byte) error {
	if f.err != nil {
		return f.err
	}
	if _, err := f.inner.Write(b); err != nil {
		f.err = err
		return err
	}
	return nil
}

// Close flushes the final (unterminated) line — the token after the last '\n',
// which bytes.Split keeps as the final element. It must be called exactly once
// after the last Write to reproduce bytes.Join's trailing-token semantics.
func (f *provenanceFilter) Close() error {
	if f.err != nil {
		return f.err
	}
	if err := f.flushLine(false); err != nil {
		return err
	}
	return f.err
}

// streamHeadCap is the size of the provenance-stripped "head" window
// StreamDigest captures for diff snippets: 64 KiB is plenty to locate the first
// divergence in a header/early data region while keeping memory bounded.
const streamHeadCap = 64 << 10

// headWriter forwards every byte to an inner writer (an md5 hash) while
// retaining the first streamHeadCap bytes for diff snippets, so the captured
// head is itself O(64 KiB).
type headWriter struct {
	inner io.Writer
	head  []byte
}

func (h *headWriter) Write(p []byte) (int, error) {
	if len(h.head) < streamHeadCap {
		room := streamHeadCap - len(h.head)
		if room > len(p) {
			room = len(p)
		}
		h.head = append(h.head, p[:room]...)
	}
	return h.inner.Write(p)
}

// StreamDigester is the WRITER-side streaming digester: bytes written to it are
// provenance-stripped, fed into an md5 hash, and the first ~64 KiB of stripped
// output is retained as the head for diff snippets. A child process can write
// its stdout straight into a StreamDigester (as an exec.Cmd.Stdout sink) so its
// entire output is normalized and hashed without ever buffering more than one
// partial line plus the 64 KiB head — the memory-safe core of the realparity
// comparison path. After all writes, call Close once, then Sum/Head.
//
// It is the writer-facing twin of StreamDigest: for any byte slice b, writing b
// (in any chunking) to a StreamDigester and Close-ing it yields Sum() ==
// md5.Sum(stripProvenance(b)), the exact normalization CompareByteExact applies.
type StreamDigester struct {
	hash hash.Hash
	head *headWriter
	pf   *provenanceFilter
}

// NewStreamDigester returns a ready StreamDigester. Write the stream into it,
// Close it once, then read Sum()/Head().
func NewStreamDigester() *StreamDigester {
	h := md5.New()
	hw := &headWriter{inner: h}
	return &StreamDigester{hash: h, head: hw, pf: newProvenanceFilter(hw)}
}

// Write feeds bytes through the provenance filter into the hash and head window.
func (d *StreamDigester) Write(p []byte) (int, error) { return d.pf.Write(p) }

// Close flushes the final partial line. It must be called once after the last
// Write, before Sum/Head.
func (d *StreamDigester) Close() error { return d.pf.Close() }

// Sum returns the md5 of the provenance-stripped stream. Call after Close.
func (d *StreamDigester) Sum() [md5.Size]byte {
	var sum [md5.Size]byte
	copy(sum[:], d.hash.Sum(nil))
	return sum
}

// Head returns the first ~64 KiB of provenance-stripped output, for diff
// snippets on divergence. Memory is O(64 KiB).
func (d *StreamDigester) Head() []byte { return d.head.head }

// StreamDigest reads r to EOF, forwarding it through the provenance filter into
// an md5 hash, and returns the digest of the provenance-stripped stream along
// with the first ~64 KiB of that stripped output (the "head") for diff
// snippets. It is the streaming equivalent of md5(stripProvenance(all-of-r)):
// for any byte slice b, StreamDigest(bytes.NewReader(b)).sum equals
// md5.Sum(stripProvenance(b)) regardless of how b is chunked across Reads.
// Memory is O(64 KiB) plus one partial line, never O(len(r)).
func StreamDigest(r io.Reader) (sum [md5.Size]byte, head []byte, err error) {
	d := NewStreamDigester()
	if _, err = io.Copy(d, r); err != nil {
		return sum, d.Head(), err
	}
	if err = d.Close(); err != nil {
		return sum, d.Head(), err
	}
	return d.Sum(), d.Head(), nil
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
