// Package bed12tobed6 splits BED12 records into their constituent BED6 blocks,
// mirroring the behaviour of `bedtools bed12tobed6`.
//
// Each input record must have at least the 12 standard BED columns. For each
// of its blocks (BlockCount, BlockSizes, BlockStarts) one BED6 record is
// emitted with:
//
//   - chrom    = input chrom
//   - start    = inputStart + blockStart
//   - end      = inputStart + blockStart + blockSize
//   - name     = the input name (column 4)
//   - score    = the block index when -n is set, otherwise the parent
//     record's score (column 5), carried unchanged onto each block
//   - strand   = the input strand
//
// For any strand other than exactly "+" (i.e. "-", ".", or empty), `-n`
// reverses the per-block numbering so the first emitted block carries the
// highest index (matches upstream `bed12tobed6 -n`, which numbers i+1 only
// when strand == "+"; test case t5 covers the '-' strand).
//
// Records whose column count is not exactly 12 are treated by upstream as a
// single implicit block spanning the whole feature and are normalised to a
// 6-column BED6 record (chrom, start, end, name, score, strand), with empty
// strings substituted for any absent trailing columns. When -n is set such a
// record's score becomes the block number, which is always 1 for a single
// block. This mirrors bedtools' GetBedBlocks (fields.size() != 12 pushes the
// whole record) followed by ProcessBed always printing six fields.
package bed12tobed6

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// FieldCountError reports a record whose column count differs from the count
// established by the first data record. Upstream bedtools sets its expected
// BED field count from the first non-header line and aborts with
// "Differing number of BED fields encountered at line: N.  Exiting..." the
// moment a later line disagrees (bedFile.h parseBedLine). Error reproduces
// that message byte-for-byte so the CLI can print it verbatim.
type FieldCountError struct {
	// Line is the 1-based physical line number of the offending record.
	Line int
}

// Error returns the upstream bedtools diagnostic verbatim.
func (e *FieldCountError) Error() string {
	return fmt.Sprintf("Differing number of BED fields encountered at line: %d.  Exiting...", e.Line)
}

// Options configures the conversion.
type Options struct {
	// NumberBlocks, when true, sets the per-output `score` column to the
	// 1-based block index. For any strand other than exactly "+" the
	// numbering is reversed (last block becomes block 1) to match upstream
	// `bedtools bed12tobed6 -n` (which numbers i+1 only when strand == "+").
	NumberBlocks bool
}

// Convert reads BED12 records from in, splits each into BED6 records per
// block, and writes them to out. It returns the number of output records
// written.
//
// Lines starting with '#', 'track', or 'browser', plus blank lines, are
// silently skipped. Records whose column count is not exactly 12 are
// normalised to a 6-column BED6 record (see the package doc); 12-column
// records with zero blocks are emitted unchanged as a defensive pass-through.
func Convert(in io.Reader, out io.Writer, opts Options) (int, error) {
	scanner := bufio.NewScanner(in)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	written := 0
	lineNo := 0
	// expectFields mirrors bedtools' bedType: it is set from the first data
	// record's column count and every subsequent record must match it. -1
	// means "not yet established".
	expectFields := -1
	for scanner.Scan() {
		lineNo++
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
		// Enforce a uniform column count across records, exactly as upstream
		// bedtools does. The first data record fixes the expected count; any
		// later record with a different count aborts with the upstream
		// message. Blank/header lines above are skipped and never counted.
		if expectFields == -1 {
			expectFields = len(fields)
		} else if len(fields) != expectFields {
			return written, &FieldCountError{Line: lineNo}
		}
		if len(fields) != 12 {
			// Not a BED12 record. Upstream (GetBedBlocks) treats it as a
			// single whole-feature block and ProcessBed always prints six
			// columns, substituting empty strings for absent trailing
			// fields. Re-emit as a normalised BED6 record.
			if len(fields) < 3 {
				return written, fmt.Errorf("line %d: expected at least 3 columns, got %d", lineNo, len(fields))
			}
			startVal, err := strconv.Atoi(strings.TrimSpace(fields[1]))
			if err != nil {
				return written, fmt.Errorf("line %d: invalid chromStart %q: %v", lineNo, fields[1], err)
			}
			endVal, err := strconv.Atoi(strings.TrimSpace(fields[2]))
			if err != nil {
				return written, fmt.Errorf("line %d: invalid chromEnd %q: %v", lineNo, fields[2], err)
			}
			name := ""
			if len(fields) >= 4 {
				name = fields[3]
			}
			score := ""
			if len(fields) >= 5 {
				score = fields[4]
			}
			strand := ""
			if len(fields) >= 6 {
				strand = fields[5]
			}
			if opts.NumberBlocks {
				// A single implicit block: strand=="+" gives i+1==1, any
				// other strand gives size()-i==1, so the number is always 1.
				score = "1"
			}
			row := []string{
				fields[0],
				strconv.Itoa(startVal),
				strconv.Itoa(endVal),
				name,
				score,
				strand,
			}
			if _, err := bw.WriteString(strings.Join(row, "\t")); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			written++
			continue
		}

		chrom := fields[0]
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromStart %q: %v", lineNo, fields[1], err)
		}
		name := fields[3]
		parentScore := fields[4]
		strand := fields[5]
		blockCount, err := strconv.Atoi(strings.TrimSpace(fields[9]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid blockCount %q: %v", lineNo, fields[9], err)
		}
		if blockCount <= 0 {
			// No blocks → pass through unchanged.
			if _, err := bw.WriteString(raw); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			written++
			continue
		}

		sizes, err := parseIntList(fields[10])
		if err != nil {
			return written, fmt.Errorf("line %d: invalid blockSizes %q: %v", lineNo, fields[10], err)
		}
		starts, err := parseIntList(fields[11])
		if err != nil {
			return written, fmt.Errorf("line %d: invalid blockStarts %q: %v", lineNo, fields[11], err)
		}
		if len(sizes) != blockCount || len(starts) != blockCount {
			return written, fmt.Errorf("line %d: blockCount=%d but got %d sizes / %d starts", lineNo, blockCount, len(sizes), len(starts))
		}

		for i := 0; i < blockCount; i++ {
			bstart := start + starts[i]
			bend := bstart + sizes[i]
			score := parentScore
			if opts.NumberBlocks {
				// Upstream numbers blocks i+1 only when the strand is exactly
				// "+"; for every other strand value ("-", ".", or empty) it
				// reverses the numbering to blockCount-i (bed12ToBed6.cpp).
				idx := i + 1
				if strand != "+" {
					idx = blockCount - i
				}
				score = strconv.Itoa(idx)
			}
			row := []string{
				chrom,
				strconv.Itoa(bstart),
				strconv.Itoa(bend),
				name,
				score,
				strand,
			}
			if _, err := bw.WriteString(strings.Join(row, "\t")); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			written++
		}
	}
	if err := scanner.Err(); err != nil {
		return written, err
	}
	return written, nil
}

// parseIntList parses a trailing-comma-tolerant list of integers ("10,20,30,").
func parseIntList(s string) ([]int, error) {
	s = strings.TrimSpace(s)
	s = strings.TrimSuffix(s, ",")
	if s == "" {
		return nil, nil
	}
	parts := strings.Split(s, ",")
	out := make([]int, len(parts))
	for i, p := range parts {
		v, err := strconv.Atoi(strings.TrimSpace(p))
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}
