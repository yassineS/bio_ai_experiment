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

	bRecs, err := readAll(bR)
	if err != nil {
		return 0, fmt.Errorf("reading B: %w", err)
	}

	// Index B by chromosome in file order, stamping each record's in-chromosome
	// order so the bin tie-break can restore it. Upstream bins B by its
	// ORIGINAL coordinates.
	byChrom := make(map[string][]*rec)
	for _, b := range bRecs {
		b.order = len(byChrom[b.chrom])
		byChrom[b.chrom] = append(byChrom[b.chrom], b)
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	aReader := bufio.NewScanner(aR)
	aReader.Buffer(make([]byte, 64*1024), 16*1024*1024)
	written := 0
	for aReader.Scan() {
		raw := aReader.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		a, perr := parseRec(raw)
		if perr != nil {
			return written, fmt.Errorf("reading A: %w", perr)
		}

		// Expand A's interval by the requested window to form the search range.
		winStart, winEnd := addWindow(a, opts)

		hits := findHits(a, winStart, winEnd, byChrom[a.chrom], opts)

		switch {
		case opts.Invert:
			if len(hits) == 0 {
				if err := writeLine(bw, a.line); err != nil {
					return written, err
				}
				written++
			}
		case opts.Count:
			if err := writeLine(bw, a.line+"\t"+strconv.Itoa(len(hits))); err != nil {
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
				var out string
				switch {
				case opts.WriteB:
					out = h.line
				case opts.WriteA:
					out = a.line
				default:
					// Default and -wa+-wb both emit A<TAB>B.
					out = a.line + "\t" + h.line
				}
				if err := writeLine(bw, out); err != nil {
					return written, err
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

// findHits returns the B records overlapping A's window [winStart, winEnd),
// applying the strand filter, in upstream's bin-traversal order. The overlap
// test matches upstream FindWindowOverlaps: a positive intersection of the B
// record with the window, and a positive overlap fraction relative to A's
// original length.
func findHits(a *rec, winStart, winEnd int, bRecs []*rec, opts Options) []*rec {
	aLen := a.end - a.start
	var hits []*rec
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

// readAll loads every BED record from r into raw recs, preserving the verbatim
// line text.
func readAll(r io.Reader) ([]*rec, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []*rec
	for sc.Scan() {
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		rc, err := parseRec(raw)
		if err != nil {
			return nil, err
		}
		out = append(out, rc)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// parseRec parses the chrom/start/end (and strand, if present at column 6) from
// a raw BED line while retaining the full text for verbatim output.
func parseRec(raw string) (*rec, error) {
	fields := strings.Split(raw, "\t")
	if len(fields) < 3 {
		return nil, fmt.Errorf("BED record must have at least 3 fields, got %d: %q", len(fields), raw)
	}
	start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromStart %q: %w", fields[1], err)
	}
	end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromEnd %q: %w", fields[2], err)
	}
	rc := &rec{line: raw, chrom: fields[0], start: start, end: end}
	if len(fields) >= 6 {
		rc.strand = strings.TrimSpace(fields[5])
	}
	return rc, nil
}

func writeLine(bw *bufio.Writer, s string) error {
	if _, err := bw.WriteString(s); err != nil {
		return err
	}
	return bw.WriteByte('\n')
}
