// Package bedslop extends BED intervals by a fixed or percentage amount,
// clipping the result to chromosome boundaries. It mirrors the behaviour of
// `bedtools slop`.
//
// Records may be extended symmetrically (-b N), asymmetrically (-l L -r R), or
// as a fraction of the interval length (--pct). When -s/--strand is set,
// records on the '-' strand have their left/right semantics swapped before the
// extension is applied (i.e. left/right are interpreted relative to the
// transcribed strand). Intervals are clipped to [0, chromSize]; any interval
// whose extended length is non-positive is dropped with a stderr warning that
// includes the original input line.
//
// Records are kept as raw text lines so the full set of input columns
// round-trips through the slop unchanged.
package bedslop

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
)

// Options bundles the configuration for Slop.
type Options struct {
	// LeftAdd is the amount to extend on the left ("upstream") side. Negative
	// values shrink the interval. When Both is true LeftAdd is ignored.
	LeftAdd float64
	// RightAdd is the amount to extend on the right ("downstream") side.
	// Negative values shrink the interval. When Both is true RightAdd is
	// ignored.
	RightAdd float64
	// Both, when set, applies BothAdd to both ends and ignores LeftAdd/RightAdd.
	Both bool
	// BothAdd is the symmetric extension used when Both is true.
	BothAdd float64
	// StrandSpec, when true, swaps left/right semantics for entries on the '-'
	// strand. A BED entry must have at least 6 columns (col 6 = strand) for
	// strand to be recognised; rows without a recognisable strand are treated
	// as '+'.
	StrandSpec bool
	// Pct, when true, treats LeftAdd/RightAdd/BothAdd as fractions of the
	// interval length rather than absolute base counts.
	Pct bool
}

// ChromSizes maps chromosome name to its total length in bases.
type ChromSizes map[string]int

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
		size, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid size %q for chromosome %q: %v", fields[1], fields[0], err)
		}
		if size < 0 {
			return nil, fmt.Errorf("negative size %d for chromosome %q", size, fields[0])
		}
		out[fields[0]] = size
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// Slop reads BED records from in, extends them according to opts using the
// chromosome sizes in sizes, and writes the surviving records to out.
//
// A record is dropped (and a warning written to warn) when:
//
//   - its chromosome is not in sizes; or
//   - the extension shrinks it to a non-positive length; or
//   - its chromStart/chromEnd cannot be parsed.
//
// Slop returns the number of records written.
func Slop(in io.Reader, out io.Writer, warn io.Writer, sizes ChromSizes, opts Options) (int, error) {
	scanner := bufio.NewScanner(in)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	written := 0
	for scanner.Scan() {
		raw := scanner.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return written, fmt.Errorf("BED record must have at least 3 fields, got %d: %q", len(fields), raw)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return written, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return written, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		chrom := fields[0]
		chromSize, ok := sizes[chrom]
		if !ok {
			if warn != nil {
				fmt.Fprintf(warn, "warning: chromosome %q not in genome file; dropping record: %s\n", chrom, raw)
			}
			continue
		}

		strand := "+"
		if len(fields) >= 6 {
			strand = fields[5]
		}

		left, right := computeExtensions(start, end, strand, opts)

		newStart := start - left
		newEnd := end + right
		if newStart < 0 {
			newStart = 0
		}
		if newEnd > chromSize {
			newEnd = chromSize
		}
		if newStart >= newEnd {
			if warn != nil {
				fmt.Fprintf(warn, "warning: slop produced empty interval; dropping record: %s\n", raw)
			}
			continue
		}
		fields[1] = strconv.Itoa(newStart)
		fields[2] = strconv.Itoa(newEnd)
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

// computeExtensions returns the integer number of bases to add to the left
// and right of an interval given start/end/strand and opts. Negative results
// indicate the interval should shrink. With Pct=true the amounts are scaled
// by the interval length and rounded to the nearest integer.
func computeExtensions(start, end int, strand string, opts Options) (left, right int) {
	length := end - start
	var leftRaw, rightRaw float64
	if opts.Both {
		leftRaw = opts.BothAdd
		rightRaw = opts.BothAdd
	} else {
		leftRaw = opts.LeftAdd
		rightRaw = opts.RightAdd
	}
	if opts.StrandSpec && strand == "-" {
		leftRaw, rightRaw = rightRaw, leftRaw
	}
	if opts.Pct {
		leftRaw = leftRaw * float64(length)
		rightRaw = rightRaw * float64(length)
	}
	return int(math.Round(leftRaw)), int(math.Round(rightRaw))
}
