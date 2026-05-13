// Package bedflank emits the flanking regions of each BED interval, mirroring
// the behaviour of `bedtools flank`.
//
// For each interval [s, e) on a chromosome of size C, bedflank emits up to two
// new intervals:
//
//   - left flank:  [max(0, s-l), s)
//   - right flank: [e, min(C, e+r))
//
// Empty flanks are skipped. The original interval itself is NOT emitted (that
// would be the behaviour of `bedtools slop`). When -s/--strand is set, the
// flanks are interpreted relative to the transcribed strand: on '-' strand
// records, l and r are swapped before the flanks are computed.
//
// With Pct=true, the left and right values are treated as fractions of the
// interval length rather than absolute base counts. Records preserve all of
// their input columns (BED3 -> BED3, BED6 -> BED6, etc.).
package bedflank

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
)

// Options configures Flank.
type Options struct {
	// LeftAdd is the size of the left ("upstream") flank.
	LeftAdd float64
	// RightAdd is the size of the right ("downstream") flank.
	RightAdd float64
	// Both, when set, applies BothAdd to both flanks (ignores LeftAdd/RightAdd).
	Both bool
	// BothAdd is the symmetric flank size used when Both is true.
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

// ReadChromSizes parses a chrom-sizes file (one `chrom<TAB>size` per line).
// It also accepts samtools-style .fai files (uses the first two whitespace-
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

// Flank reads BED records from in, emits the surviving flanks to out, and
// writes warnings (e.g. for records on unknown chromosomes) to warn. It
// returns the number of records written.
func Flank(in io.Reader, out io.Writer, warn io.Writer, sizes ChromSizes, opts Options) (int, error) {
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
		if end < start {
			return written, fmt.Errorf("end < start (%d < %d)", end, start)
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

		leftAdd, rightAdd := computeExtensions(start, end, strand, opts)

		// Left flank: [max(0, start-left), start)
		if leftAdd > 0 {
			ls := start - leftAdd
			if ls < 0 {
				ls = 0
			}
			le := start
			if le > chromSize {
				le = chromSize
			}
			if le > ls {
				if err := writeFlank(bw, fields, ls, le); err != nil {
					return written, err
				}
				written++
			}
		}
		// Right flank: [end, min(chromSize, end+right))
		if rightAdd > 0 {
			rs := end
			if rs < 0 {
				rs = 0
			}
			re := end + rightAdd
			if re > chromSize {
				re = chromSize
			}
			if re > rs {
				if err := writeFlank(bw, fields, rs, re); err != nil {
					return written, err
				}
				written++
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return written, err
	}
	return written, nil
}

// writeFlank writes one flank record, replacing fields[1]/fields[2] with the
// requested coordinates and keeping all other columns untouched.
func writeFlank(bw *bufio.Writer, fields []string, start, end int) error {
	// Take a defensive copy of the [1] and [2] entries so we don't clobber
	// the per-record buffer between left/right flanks.
	out := make([]string, len(fields))
	copy(out, fields)
	out[1] = strconv.Itoa(start)
	out[2] = strconv.Itoa(end)
	if _, err := bw.WriteString(strings.Join(out, "\t")); err != nil {
		return err
	}
	return bw.WriteByte('\n')
}

// computeExtensions returns the (non-negative) integer flank sizes to apply
// on the left and right of an interval.
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
	l := int(math.Round(leftRaw))
	r := int(math.Round(rightRaw))
	if l < 0 {
		l = 0
	}
	if r < 0 {
		r = 0
	}
	return l, r
}
