// Package bedshift shifts each feature of a BED/GFF/VCF file by a requested
// number of base pairs, clamping the result to the chromosome bounds. It
// mirrors the behaviour of `bedtools shift` (aka shiftBed).
//
// Features on the '+' strand (or with no strand) are shifted by ShiftPlus;
// features on the '-' strand are shifted by ShiftMinus. With Fractional set,
// the shift is interpreted as a fraction of the feature's length. Coordinates
// are clamped so the start stays within [0, chromSize-1] and the end within
// [1, chromSize], exactly as upstream does — including upstream's truncation of
// the floating-point shifted coordinate toward zero.
//
// Records round-trip as raw text so every input column is preserved.
package bedshift

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Options bundles the configuration for Shift.
type Options struct {
	// ShiftPlus is the number of base pairs to shift features on the '+'
	// strand (and features without a recognisable strand). It is stored as a
	// 32-bit float to match upstream, whose shift amounts are C floats.
	ShiftPlus float32
	// ShiftMinus is the number of base pairs to shift features on the '-'
	// strand. It too is a 32-bit float, matching upstream.
	ShiftMinus float32
	// Fractional, when true, treats ShiftPlus/ShiftMinus as fractions of each
	// feature's length rather than absolute base counts.
	Fractional bool
}

// ChromSizes maps chromosome name to its total length in bases.
type ChromSizes map[string]int64

// ReadChromSizes parses a chrom-sizes file (one `chrom\tsize` per line). It
// also accepts samtools-style .fai files (uses the first two whitespace-
// separated columns). Blank lines and comments (lines starting with '#') are
// skipped.
func ReadChromSizes(r io.Reader) (ChromSizes, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	out := make(ChromSizes)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("chrom-sizes line %q must have at least 2 fields", line)
		}
		size, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid size %q for chromosome %q: %v", fields[1], fields[0], err)
		}
		out[fields[0]] = size
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// chromSize returns the size upstream's GenomeFile::getChromSize would yield:
// the stored size for a known chromosome, or -1 for an unknown one (upstream
// returns (uint64_t)-1, which the caller reinterprets as the signed value -1).
func chromSize(sizes ChromSizes, chrom string) int64 {
	if s, ok := sizes[chrom]; ok {
		return s
	}
	return -1
}

// Shift reads BED records from in, shifts each according to opts using the
// chromosome sizes in sizes, and writes the results to out.
//
// Shift matches upstream's clamping semantics exactly: the shifted start is
// pinned to [0, chromSize-1] and the shifted end to [1, chromSize], with the
// floating-point shifted coordinate truncated toward zero when it lands inside
// the bounds. Comment, track, and browser lines are passed through only when
// printHeader is true (mirroring upstream's PrintHeader), otherwise dropped.
//
// Shift returns the number of records written.
func Shift(in io.Reader, out io.Writer, sizes ChromSizes, opts Options, printHeader bool) (int, error) {
	scanner := bufio.NewScanner(in)
	scanner.Buffer(make([]byte, 64*1024), 64*1024*1024)
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	written := 0
	headerDone := false
	for scanner.Scan() {
		raw := scanner.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			continue
		}
		// Header lines (comments, track, browser) are emitted before the first
		// data record when -header is requested, then suppressed thereafter.
		if isHeaderLine(trimmed) {
			if printHeader && !headerDone {
				if _, err := bw.WriteString(raw); err != nil {
					return written, err
				}
				if err := bw.WriteByte('\n'); err != nil {
					return written, err
				}
			}
			continue
		}
		headerDone = true

		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return written, fmt.Errorf("BED record must have at least 3 fields, got %d: %q", len(fields), raw)
		}
		start, err := strconv.ParseInt(strings.TrimSpace(fields[1]), 10, 64)
		if err != nil {
			return written, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.ParseInt(strings.TrimSpace(fields[2]), 10, 64)
		if err != nil {
			return written, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}

		strand := ""
		if len(fields) >= 6 {
			strand = fields[5]
		}

		newStart, newEnd := addShift(start, end, strand, chromSize(sizes, fields[0]), opts)
		fields[1] = strconv.FormatInt(newStart, 10)
		fields[2] = strconv.FormatInt(newEnd, 10)
		if _, err := bw.WriteString(strings.Join(fields, "\t")); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return written, err
		}
		written++
	}
	if err := scanner.Err(); err != nil {
		return written, err
	}
	return written, nil
}

// isHeaderLine reports whether a (trimmed) line is a header/comment line that
// upstream's BedFile treats as non-data: comments and track/browser lines.
func isHeaderLine(trimmed string) bool {
	return strings.HasPrefix(trimmed, "#") ||
		strings.HasPrefix(trimmed, "track") ||
		strings.HasPrefix(trimmed, "browser")
}

// addShift applies the shift to one feature's coordinates and returns the
// clamped result. It reproduces upstream BedShift::AddShift, including the
// double-precision arithmetic and truncation toward zero on assignment to the
// integer coordinate fields.
func addShift(start, end int64, strand string, size int64, opts Options) (int64, int64) {
	var shift float64
	if strand == "-" {
		shift = float64(opts.ShiftMinus)
	} else {
		shift = float64(opts.ShiftPlus)
	}
	if opts.Fractional {
		shift = shift * float64(end-start)
	}

	startF := float64(start) + shift
	var newStart int64
	switch {
	case startF < 0:
		newStart = 0
	case startF > float64(size-1):
		newStart = size - 1
	default:
		newStart = int64(startF) // C++ double->int64: truncation toward zero.
	}

	endF := float64(end) + shift
	var newEnd int64
	switch {
	case endF <= 0:
		newEnd = 1
	case endF > float64(size):
		newEnd = size
	default:
		newEnd = int64(endF) // truncation toward zero.
	}
	return newStart, newEnd
}
