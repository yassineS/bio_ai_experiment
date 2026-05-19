// BEDPE (BED for Paired-End features) parser and writer.
//
// BEDPE is a 10-column tab-separated format used to represent two genomic
// intervals per record (typical use cases: structural-variant breakpoints,
// paired-end read pairs, Hi-C contacts). Columns:
//
//	chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 [extra...]
//
// Coordinates are 0-based, half-open — the same convention as BED. The two
// ends are independent intervals; chrom1 may differ from chrom2.
//
// A sentinel chrom value of "." means the corresponding end is unaligned
// (this is what bedtools writes when one end of a paired-end read has no
// mapping). Unaligned ends have start/end of -1.

package bed

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// BEDPE represents one BEDPE record. The two ends are independent intervals.
type BEDPE struct {
	Chrom1  string
	Start1  int
	End1    int
	Chrom2  string
	Start2  int
	End2    int
	Name    string
	Score   string // Score is kept as string to match bedtools (it allows non-numeric).
	Strand1 string
	Strand2 string
	// Extra holds any trailing tab-separated columns beyond the first 10.
	Extra []string
}

// BEDPEReader provides sequential access to BEDPE records.
type BEDPEReader struct {
	scanner *bufio.Scanner
	err     error
	line    int
}

// NewBEDPEReader returns a new BEDPE reader wrapping r.
func NewBEDPEReader(r io.Reader) *BEDPEReader {
	sc := bufio.NewScanner(r)
	// Allow long lines (some BEDPE records carry many extras / large CIGARs).
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	return &BEDPEReader{scanner: sc}
}

// Read returns the next BEDPE record, or io.EOF when the stream is exhausted.
// Empty/comment/track/browser lines are skipped. Records with fewer than 10
// fields return an error.
func (r *BEDPEReader) Read() (*BEDPE, error) {
	if r.err != nil {
		return nil, r.err
	}
	for r.scanner.Scan() {
		r.line++
		raw := r.scanner.Text()
		line := strings.TrimRight(raw, "\r\n")
		if line == "" || strings.HasPrefix(line, "#") ||
			strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			return nil, fmt.Errorf("BEDPE line %d: need >=10 tab-separated fields, got %d", r.line, len(fields))
		}
		rec := &BEDPE{
			Chrom1:  fields[0],
			Chrom2:  fields[3],
			Name:    fields[6],
			Score:   fields[7],
			Strand1: fields[8],
			Strand2: fields[9],
		}
		var err error
		if rec.Start1, err = strconv.Atoi(fields[1]); err != nil {
			return nil, fmt.Errorf("BEDPE line %d: invalid start1 %q: %v", r.line, fields[1], err)
		}
		if rec.End1, err = strconv.Atoi(fields[2]); err != nil {
			return nil, fmt.Errorf("BEDPE line %d: invalid end1 %q: %v", r.line, fields[2], err)
		}
		if rec.Start2, err = strconv.Atoi(fields[4]); err != nil {
			return nil, fmt.Errorf("BEDPE line %d: invalid start2 %q: %v", r.line, fields[4], err)
		}
		if rec.End2, err = strconv.Atoi(fields[5]); err != nil {
			return nil, fmt.Errorf("BEDPE line %d: invalid end2 %q: %v", r.line, fields[5], err)
		}
		if len(fields) > 10 {
			rec.Extra = append([]string(nil), fields[10:]...)
		}
		return rec, nil
	}
	if err := r.scanner.Err(); err != nil {
		r.err = err
		return nil, err
	}
	r.err = io.EOF
	return nil, io.EOF
}

// ReadAll reads every BEDPE record from the underlying stream.
func (r *BEDPEReader) ReadAll() ([]*BEDPE, error) {
	var out []*BEDPE
	for {
		rec, err := r.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, rec)
	}
}

// BEDPEWriter writes BEDPE records as 10-column tab-separated lines.
type BEDPEWriter struct {
	w *bufio.Writer
}

// NewBEDPEWriter returns a buffered writer that emits BEDPE records.
func NewBEDPEWriter(w io.Writer) *BEDPEWriter {
	return &BEDPEWriter{w: bufio.NewWriter(w)}
}

// Write writes a single BEDPE record followed by a newline.
func (w *BEDPEWriter) Write(rec *BEDPE) error {
	if _, err := w.w.WriteString(rec.String()); err != nil {
		return err
	}
	return w.w.WriteByte('\n')
}

// WriteRaw appends an arbitrary already-tab-joined line followed by a newline.
// Useful when callers want to emit a BEDPE record concatenated with a partner
// record (pairtopair output format).
func (w *BEDPEWriter) WriteRaw(line string) error {
	if _, err := w.w.WriteString(line); err != nil {
		return err
	}
	return w.w.WriteByte('\n')
}

// Flush flushes the underlying buffer.
func (w *BEDPEWriter) Flush() error { return w.w.Flush() }

// String returns the tab-joined 10-column representation of the record,
// including any Extra columns. No trailing newline.
func (b *BEDPE) String() string {
	var sb strings.Builder
	sb.Grow(64)
	sb.WriteString(b.Chrom1)
	sb.WriteByte('\t')
	sb.WriteString(strconv.Itoa(b.Start1))
	sb.WriteByte('\t')
	sb.WriteString(strconv.Itoa(b.End1))
	sb.WriteByte('\t')
	sb.WriteString(b.Chrom2)
	sb.WriteByte('\t')
	sb.WriteString(strconv.Itoa(b.Start2))
	sb.WriteByte('\t')
	sb.WriteString(strconv.Itoa(b.End2))
	sb.WriteByte('\t')
	sb.WriteString(b.Name)
	sb.WriteByte('\t')
	sb.WriteString(b.Score)
	sb.WriteByte('\t')
	sb.WriteString(b.Strand1)
	sb.WriteByte('\t')
	sb.WriteString(b.Strand2)
	for _, e := range b.Extra {
		sb.WriteByte('\t')
		sb.WriteString(e)
	}
	return sb.String()
}

// End1Record returns the first end as a plain BED Record.
// Unaligned ends (Chrom == ".") are returned with the original ("." chrom and
// possibly negative coordinates) — callers should skip them.
func (b *BEDPE) End1Record() *Record {
	return &Record{
		Chrom:      b.Chrom1,
		ChromStart: b.Start1,
		ChromEnd:   b.End1,
		Name:       b.Name,
		Strand:     b.Strand1,
	}
}

// End2Record returns the second end as a plain BED Record.
func (b *BEDPE) End2Record() *Record {
	return &Record{
		Chrom:      b.Chrom2,
		ChromStart: b.Start2,
		ChromEnd:   b.End2,
		Name:       b.Name,
		Strand:     b.Strand2,
	}
}
