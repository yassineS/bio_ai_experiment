// Package bedspacing implements `bedtools spacing`: it walks a BED file (or
// any tab-delimited file with chrom/start/end in columns 1-3) and appends a
// column that reports the spacing between adjacent intervals on the same
// chromosome.
//
// Following upstream, the spacing token is one of:
//
//   - "."  — first interval on its chromosome.
//   - "-1" — overlaps the previous interval on the same chromosome.
//   - "0"  — exactly abuts the previous interval (prev.end == this.start).
//   - N    — otherwise: this.start - prev.end (positive gap in bases).
//
// The "previous" interval is tracked per-chromosome; ordering is whatever
// the input provides — bedspacing does not sort. To get the conventional
// genome-sorted spacing report, pipe a sorted BED in (`bedsort` or
// `sort -k1,1 -k2,2n`).
//
// The input columns are preserved verbatim; the spacing token is appended as
// a new trailing tab-separated column. Header lines (`#`, `track`, `browser`)
// and blank lines are passed through unchanged. Mirrors upstream's
// `bedtools spacing -i a.bed`.
//
// SAM/BAM input is also accepted: the stream is auto-detected and each mapped
// alignment is converted to its BED12 representation (matching upstream
// `bedtools spacing -i in.bam -bed`) before the spacing column is appended.
package bedspacing

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnbed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Spacing reads records from r, computes the spacing to the previous interval
// on the same chromosome, and writes each record back to w with the spacing
// appended as a new column. Returns the number of data records written.
//
// The input may be either a BED-like text stream or a SAM/BAM alignment
// stream; the two are auto-detected from the leading bytes (mirroring
// upstream `bedtools spacing -i <bed/gff/vcf/bam>`). For SAM/BAM input each
// mapped alignment is first converted to its BED12 representation — the same
// conversion bedtools applies under `-bed` — and the spacing column is then
// appended to that BED12 line. The feature span used for spacing is the whole
// reference span [start, end) of the alignment; there is no `-split` for
// spacing.
func Spacing(r io.Reader, w io.Writer) (int, error) {
	br := bufio.NewReader(r)
	head, _ := br.Peek(16)
	if alnbed.LooksLikeAlignment(head) {
		return spacingAlignment(br, w)
	}
	return spacingText(br, w)
}

// spacingTracker computes the spacing token for an interval against the most
// recent interval seen on the same chromosome, then records the new end.
type spacingTracker struct {
	prevEnd map[string]int
}

// newSpacingTracker returns an empty per-chromosome spacing tracker.
func newSpacingTracker() *spacingTracker {
	return &spacingTracker{prevEnd: make(map[string]int)}
}

// token returns the spacing token for [start, end) on chrom and advances the
// per-chromosome "previous end" pointer. The token is ".", "-1", "0", or the
// positive gap (start - prevEnd), matching upstream spacingFile.cpp.
func (s *spacingTracker) token(chrom string, start, end int) string {
	var tok string
	if pe, ok := s.prevEnd[chrom]; !ok {
		tok = "."
	} else if start < pe {
		tok = "-1"
	} else if start == pe {
		tok = "0"
	} else {
		tok = strconv.Itoa(start - pe)
	}
	// Upstream stores only the immediately preceding record per chrom
	// (reference_code/bedtools/src/spacingFile/spacingFile.cpp); the previous
	// pointer advances to this record's end after each comparison.
	s.prevEnd[chrom] = end
	return tok
}

// spacingAlignment computes spacing for a SAM/BAM stream. Each mapped
// alignment is converted to its BED12 line (the upstream `-bed` rendering)
// and the spacing token is appended as a trailing column. Spacing is computed
// on the whole reference span of each alignment.
func spacingAlignment(r io.Reader, w io.Writer) (int, error) {
	ar, err := alnbed.NewReader(r)
	if err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	track := newSpacingTracker()
	written := 0
	for {
		rec, err := ar.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return written, err
		}
		tok := track.token(rec.Chrom, rec.ChromStart, rec.ChromEnd)
		if _, err := bw.WriteString(formatBED12(rec)); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\t'); err != nil {
			return written, err
		}
		if _, err := bw.WriteString(tok); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return written, err
		}
		written++
	}
	return written, nil
}

// formatBED12 renders an alignment-derived BED12 record the way upstream
// bedtools prints BAM-as-BED under `-bed`: 12 tab-separated columns with the
// blockSizes and blockStarts lists keeping their trailing commas (e.g.
// "10,10," and "0,20,"). This matches `bedtools spacing -i in.bam -bed`
// byte-for-byte.
func formatBED12(rec *bed.Record) string {
	var b strings.Builder
	b.WriteString(rec.Chrom)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(rec.ChromStart))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(rec.ChromEnd))
	b.WriteByte('\t')
	b.WriteString(rec.Name)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(rec.Score))
	b.WriteByte('\t')
	b.WriteString(rec.Strand)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(rec.ThickStart))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(rec.ThickEnd))
	b.WriteByte('\t')
	b.WriteString(rec.ItemRGB)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(rec.BlockCount))
	b.WriteByte('\t')
	for _, sz := range rec.BlockSizes {
		b.WriteString(strconv.Itoa(sz))
		b.WriteByte(',')
	}
	b.WriteByte('\t')
	for _, st := range rec.BlockStarts {
		b.WriteString(strconv.Itoa(st))
		b.WriteByte(',')
	}
	return b.String()
}

// spacingText computes spacing for a BED-like text stream, echoing each input
// line verbatim with the spacing token appended as a trailing column.
func spacingText(r io.Reader, w io.Writer) (int, error) {
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	// Tracks the most recent end seen per chromosome.
	track := newSpacingTracker()
	written := 0
	lineNo := 0
	for sc.Scan() {
		lineNo++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			if _, err := bw.WriteString(raw); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			continue
		}
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			if _, err := bw.WriteString(raw); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return written, fmt.Errorf("line %d: BED record must have at least 3 columns: %q", lineNo, raw)
		}
		chrom := fields[0]
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromStart %q: %v", lineNo, fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromEnd %q: %v", lineNo, fields[2], err)
		}

		spacing := track.token(chrom, start, end)

		if _, err := bw.WriteString(raw); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\t'); err != nil {
			return written, err
		}
		if _, err := bw.WriteString(spacing); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return written, err
		}
		written++
	}
	if err := sc.Err(); err != nil {
		return written, err
	}
	return written, nil
}
