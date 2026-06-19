// Package bedwindow ports `bedtools window` (aka windowBed): for each feature
// in A it examines a "window" — A's interval expanded by -w (or -l/-r) base
// pairs on each side — and reports every feature in B that overlaps that
// window. For each overlap the entire A and B records are reported (default
// mode), or one of the alternate writer modes (-u, -c, -v) is applied.
//
// Faithful-port notes:
//
//   - The window is added to A, not B (upstream AddWindow operates on the A
//     feature, then queries the B database with the fudged coordinates). This
//     matters for the asymmetric -l/-r and strand (-sw) cases: -l extends the
//     window upstream (lower coordinates) of A and -r downstream.
//
//   - B records are indexed in upstream's UCSC binning tree using their
//     ORIGINAL coordinates, and the per-A hit order is the bin-traversal order
//     (finest bin level first, then bin number ascending, then file order),
//     NOT plain file order and NOT B-start order. See binorder.go.
//
//   - Records are kept as their raw input text so every column (BED3..BED12 and
//     beyond) round-trips verbatim. The earlier implementation re-rendered B
//     from a typed record and truncated BED12 block columns; this preserves
//     them byte-for-byte.
package bedwindow

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Options configures Window.
type Options struct {
	// Left is the bp added upstream (to lower coordinates) of each A feature
	// when searching for overlaps in B (upstream -l, or -w which sets both).
	Left int
	// Right is the bp added downstream (to higher coordinates) of each A
	// feature (upstream -r, or -w which sets both).
	Right int

	// StrandWindows mirrors upstream -sw: define -l/-r relative to A's strand,
	// so for a negative-strand A the left/right slops swap. Disabled by default.
	StrandWindows bool

	// StrandSpec mirrors upstream -sm: only report B hits on the SAME strand as
	// A. Mutually exclusive with InverseStrand.
	StrandSpec bool
	// InverseStrand mirrors upstream -Sm: only report B hits on the OPPOSITE
	// strand to A. Mutually exclusive with StrandSpec.
	InverseStrand bool

	// WriteA emits the original A record only (upstream has no -wa for window;
	// retained for the CLI's -wa convenience alias).
	WriteA bool
	// WriteB emits the original B record only.
	WriteB bool
	// WriteAB emits `A<TAB>B` for each overlap. This is also the default when no
	// writer flag is set, matching upstream `bedtools window`.
	WriteAB bool
	// Count mirrors upstream -c: emit `A<TAB>count` (one row per A, count of B
	// overlaps within the window, 0 included).
	Count bool
	// AnyHit mirrors upstream -u: emit the original A record once if it has at
	// least one B overlap.
	AnyHit bool
	// Invert mirrors upstream -v: emit only A records with NO B overlap.
	Invert bool

	// MinOverlap is the minimum bp of overlap required. Upstream requires any
	// positive overlap (>0) regardless; this is retained for callers that want
	// a stricter threshold and defaults to 1.
	MinOverlap int
}

// rec is a parsed view of one raw BED line: the verbatim text plus the fields
// needed for window/overlap logic.
type rec struct {
	line   string
	chrom  string
	start  int
	end    int
	strand string
	order  int // in-chromosome insertion (file) order, for bin tie-breaks
}

// Window reads A from aR, B from bR, and writes results to w. It returns the
// number of output lines written.
func Window(aR, bR io.Reader, w io.Writer, opts Options) (int, error) {
	if opts.MinOverlap < 1 {
		opts.MinOverlap = 1
	}
	if opts.StrandSpec && opts.InverseStrand {
		return 0, fmt.Errorf("StrandSpec (-sm) and InverseStrand (-Sm) are mutually exclusive")
	}

	// We need B's strand only when a strand filter is active; skip interning it
	// otherwise. The chrom is always needed for the per-chromosome index.
	keepBStrand := opts.StrandSpec || opts.InverseStrand
	bRecs, err := readAll(bR, keepBStrand)
	if err != nil {
		return 0, fmt.Errorf("reading B: %w", err)
	}

	// Index B by chromosome in file order, stamping each record's in-chromosome
	// order so the bin tie-break can restore it. Upstream bins B by its
	// ORIGINAL coordinates.
	byChrom := make(map[string][]*rec)
	for i := range bRecs {
		b := &bRecs[i]
		b.order = len(byChrom[b.chrom])
		byChrom[b.chrom] = append(byChrom[b.chrom], b)
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	aReader := bufio.NewScanner(aR)
	aReader.Buffer(make([]byte, 64*1024), 16*1024*1024)
	written := 0
	// Reusable scratch for the A record (parsed from the scanner buffer), the
	// per-A hit buffer, and integer formatting, so the hot loop allocates
	// nothing per record beyond what output requires.
	var a rec
	var hits []*rec
	var aCI chromInterner
	var scratch []byte
	// A strand is needed only when a strand filter is active.
	keepAStrand := keepBStrand
	for aReader.Scan() {
		lineB := bytes.TrimRight(aReader.Bytes(), "\r")
		if isHeaderOrBlank(lineB) {
			continue
		}
		if perr := parseRecBytes(lineB, keepAStrand, &aCI, &a); perr != nil {
			return written, fmt.Errorf("reading A: %w", perr)
		}

		// Expand A's interval by the requested window to form the search range.
		winStart, winEnd := addWindow(&a, opts)

		hits = findHitsInto(&a, winStart, winEnd, byChrom[a.chrom], opts, hits[:0])

		switch {
		case opts.Invert:
			if len(hits) == 0 {
				if err := writeLine(bw, a.line); err != nil {
					return written, err
				}
				written++
			}
		case opts.Count:
			scratch = scratch[:0]
			scratch = append(scratch, a.line...)
			scratch = append(scratch, '\t')
			scratch = strconv.AppendInt(scratch, int64(len(hits)), 10)
			if err := writeLineBytes(bw, scratch); err != nil {
				return written, err
			}
			written++
		case opts.AnyHit:
			if len(hits) > 0 {
				if err := writeLine(bw, a.line); err != nil {
					return written, err
				}
				written++
			}
		default:
			for _, h := range hits {
				switch {
				case opts.WriteB:
					if err := writeLine(bw, h.line); err != nil {
						return written, err
					}
				case opts.WriteA:
					if err := writeLine(bw, a.line); err != nil {
						return written, err
					}
				default:
					// Default and -wa+-wb both emit A<TAB>B. Emit the two parts
					// directly to the writer to avoid concatenating a fresh
					// string per hit.
					if err := writeAB(bw, a.line, h.line); err != nil {
						return written, err
					}
				}
				written++
			}
		}
	}
	if err := aReader.Err(); err != nil {
		return written, fmt.Errorf("reading A: %w", err)
	}
	return written, nil
}

// isHeaderOrBlank reports whether a raw line is blank or a comment/track/browser
// header line that the reader skips, testing the byte buffer without allocating.
func isHeaderOrBlank(line []byte) bool {
	trimmed := bytes.TrimSpace(line)
	return len(trimmed) == 0 || trimmed[0] == '#' ||
		bytes.HasPrefix(trimmed, trackPrefix) || bytes.HasPrefix(trimmed, browserPrefix)
}

var (
	trackPrefix   = []byte("track")
	browserPrefix = []byte("browser")
)

// writeLineBytes writes a byte slice followed by a newline.
func writeLineBytes(bw *bufio.Writer, b []byte) error {
	if _, err := bw.Write(b); err != nil {
		return err
	}
	return bw.WriteByte('\n')
}

// writeAB writes `a<TAB>b\n` directly to the writer without allocating the
// concatenated string, the dominant default-mode allocation.
func writeAB(bw *bufio.Writer, a, b string) error {
	if _, err := bw.WriteString(a); err != nil {
		return err
	}
	if err := bw.WriteByte('\t'); err != nil {
		return err
	}
	if _, err := bw.WriteString(b); err != nil {
		return err
	}
	return bw.WriteByte('\n')
}

// addWindow expands A's interval by the requested slop, mirroring upstream
// BedWindow::AddWindow. The low end is clipped at 0. With StrandWindows the
// left/right slops swap for negative-strand A features.
func addWindow(a *rec, opts Options) (start, end int) {
	left, right := opts.Left, opts.Right
	if opts.StrandWindows && a.strand == "-" {
		left, right = right, left
	}
	// Upstream uses (int)(a.start - leftSlop) > 0 to decide; otherwise clamps to
	// 0. That clamps a result of exactly 0 to 0 too, which is the same value.
	start = a.start - left
	if start < 0 {
		start = 0
	}
	end = a.end + right
	return start, end
}

// findHitsInto returns the B records overlapping A's window [winStart, winEnd),
// applying the strand filter, in upstream's bin-traversal order. The overlap
// test matches upstream FindWindowOverlaps: a positive intersection of the B
// record with the window, and a positive overlap fraction relative to A's
// original length. Hits are appended into dst (pass hits[:0] to reuse the
// backing array across A-record queries), which is returned so the caller can
// recycle it.
func findHitsInto(a *rec, winStart, winEnd int, bRecs []*rec, opts Options, dst []*rec) []*rec {
	aLen := a.end - a.start
	hits := dst
	for _, b := range bRecs {
		s := winStart
		if b.start > s {
			s = b.start
		}
		e := winEnd
		if b.end < e {
			e = b.end
		}
		if s >= e {
			continue
		}
		overlapBases := e - s
		// Upstream requires (overlapBases / aLength) > 0, i.e. any positive
		// overlap. aLength == 0 would divide by zero in C++ float math
		// (yielding inf/nan, treated as >0); guard so any overlap counts.
		if aLen > 0 && overlapBases <= 0 {
			continue
		}
		if overlapBases < opts.MinOverlap {
			continue
		}
		if !strandOK(a.strand, b.strand, opts) {
			continue
		}
		hits = append(hits, b)
	}
	orderHitsByBin(hits)
	return hits
}

// strandOK applies the -sm/-Sm strand filter. With neither set, every overlap
// passes regardless of strand.
func strandOK(aStrand, bStrand string, opts Options) bool {
	if !opts.StrandSpec && !opts.InverseStrand {
		return true
	}
	same := aStrand == bStrand
	if opts.StrandSpec {
		return same
	}
	// InverseStrand (-Sm): opposite strand required.
	return !same
}

// readAll loads every BED record from r into a flat slice of recs, preserving
// the verbatim line text. The records are stored by value (one backing slice
// growth rather than a heap allocation per record); callers take &slice[i] to
// build the per-chromosome index. keepStrand controls whether column 6 is
// interned (only the strand-filtered modes need it).
func readAll(r io.Reader, keepStrand bool) ([]rec, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []rec
	var ci chromInterner
	for sc.Scan() {
		lineB := bytes.TrimRight(sc.Bytes(), "\r")
		if isHeaderOrBlank(lineB) {
			continue
		}
		var rc rec
		if err := parseRecBytes(lineB, keepStrand, &ci, &rc); err != nil {
			return nil, err
		}
		out = append(out, rc)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// parseRecBytes parses chrom/start/end (and strand, when keepStrand) directly
// from the scanner's byte buffer into dst, avoiding the strings.Split field
// slice and the per-coordinate Atoi string allocations of the old path. The
// verbatim line is copied (string(lineB) is a fresh allocation) so it does not
// alias the reused scanner buffer; the chromosome name is interned through ci so
// a run of records on one chromosome allocates the name once.
func parseRecBytes(lineB []byte, keepStrand bool, ci *chromInterner, dst *rec) error {
	// Locate the first up-to-6 tab-delimited columns without allocating a slice.
	var cols [6][]byte
	n := 0
	begin := 0
	for i := 0; i <= len(lineB); i++ {
		if i == len(lineB) || lineB[i] == '\t' {
			if n < len(cols) {
				cols[n] = lineB[begin:i]
			}
			n++
			begin = i + 1
		}
	}
	if n < 3 {
		return fmt.Errorf("BED record must have at least 3 fields, got %d: %q", n, string(lineB))
	}
	start, err := parsePosBytes(cols[1])
	if err != nil {
		return fmt.Errorf("invalid chromStart %q: %w", cols[1], err)
	}
	end, err := parsePosBytes(cols[2])
	if err != nil {
		return fmt.Errorf("invalid chromEnd %q: %w", cols[2], err)
	}
	dst.line = string(lineB)
	dst.chrom = ci.intern(cols[0])
	dst.start = start
	dst.end = end
	dst.strand = ""
	dst.order = 0
	if keepStrand && n >= 6 {
		dst.strand = internStrand(cols[5])
	}
	return nil
}

// chromInterner caches the most recently interned chromosome name so a run of
// records sharing a chromosome (the norm for sorted BED input) allocates the
// name once rather than per record. The `c.last == string(b)` compare is
// allocation-free: the compiler does not heap-allocate a []byte->string
// conversion used only as a comparison operand.
type chromInterner struct {
	last string
}

func (c *chromInterner) intern(b []byte) string {
	if c.last != "" && c.last == string(b) {
		return c.last
	}
	s := string(b)
	c.last = s
	return s
}

// internStrand maps the byte form of a strand column to a shared string,
// allocating nothing for the only values that occur in practice.
func internStrand(b []byte) string {
	switch string(b) {
	case "":
		return ""
	case "+":
		return "+"
	case "-":
		return "-"
	case ".":
		return "."
	}
	return string(b)
}

// parsePosBytes parses a BED coordinate from a byte slice. The fast path is
// plain optional-signed digits (no allocation); anything else — embedded
// whitespace, etc. — falls back to the trimming Atoi the old path used so
// behaviour and error text round-trip.
func parsePosBytes(b []byte) (int, error) {
	if len(b) > 0 {
		i := 0
		neg := false
		if b[0] == '+' || b[0] == '-' {
			neg = b[0] == '-'
			i++
		}
		if i < len(b) {
			val := 0
			ok := true
			for ; i < len(b); i++ {
				c := b[i]
				if c < '0' || c > '9' {
					ok = false
					break
				}
				val = val*10 + int(c-'0')
			}
			if ok {
				if neg {
					val = -val
				}
				return val, nil
			}
		}
	}
	return strconv.Atoi(strings.TrimSpace(string(b)))
}

func writeLine(bw *bufio.Writer, s string) error {
	if _, err := bw.WriteString(s); err != nil {
		return err
	}
	return bw.WriteByte('\n')
}
