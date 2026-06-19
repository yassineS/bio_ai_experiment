// Fast BED input parsing for the genomecov BED path.
//
// The shared pkg/htsgo/bed reader allocates a *Record, a scanner string
// (bufio.Scanner.Text) and a []string field split per line; for genomecov's
// histogram/bedGraph paths that parse cost dominates the allocation profile
// (~96% of allocations on the default histogram). genomecov only needs a
// handful of columns per record — chrom, start, end, strand, and (under
// -split) the BED12 block columns — so this file provides a reader that parses
// those directly from the scanner's byte buffer, copying only the chromosome
// (interned across a run of one chromosome) and strand. It is byte-for-byte
// equivalent to the shared reader for the columns genomecov consumes.
package bedgenomecov

import (
	"bufio"
	"fmt"
	"io"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// fastBEDSource reads BED records directly from a byte-oriented scanner,
// reusing a single *bed.Record across calls and avoiding the per-line string
// and field-slice allocations of the shared reader. It satisfies recordSource.
//
// keepBlocks controls whether the optional BED12 block columns (fields 10–12)
// are parsed; they are only needed under -split, so the common path skips them.
type fastBEDSource struct {
	sc         *bufio.Scanner
	rec        bed.Record
	ci         chromInterner
	keepBlocks bool
	err        error
}

// newFastBEDSource builds a fastBEDSource over r. keepBlocks must be true when
// the caller needs BED12 block columns (i.e. opts.Split).
func newFastBEDSource(r io.Reader, keepBlocks bool) *fastBEDSource {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 64*1024*1024)
	return &fastBEDSource{sc: sc, keepBlocks: keepBlocks}
}

// chromInterner caches the most recently interned chromosome name so a run of
// records sharing a chromosome (the norm for sorted BED input) allocates the
// name string once rather than per record. The `c.last == string(b)` compare is
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
// allocating nothing for the values that occur in practice.
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

// Read returns the next BED record, mirroring the shared bed.Reader's line
// handling: it skips blank, '#'-comment, "track" and "browser" lines (after a
// surrounding-whitespace trim) and requires at least three columns. The
// returned *bed.Record is reused between calls and is only valid until the next
// Read; the coverage accumulator consumes it immediately, so this is safe.
func (s *fastBEDSource) Read() (*bed.Record, error) {
	if s.err != nil {
		return nil, s.err
	}
	for s.sc.Scan() {
		line := trimSpaceBytes(s.sc.Bytes())
		if len(line) == 0 || line[0] == '#' ||
			hasBytePrefix(line, "track") || hasBytePrefix(line, "browser") {
			continue
		}
		if err := s.parse(line); err != nil {
			s.err = err
			return nil, err
		}
		return &s.rec, nil
	}
	if err := s.sc.Err(); err != nil {
		s.err = err
		return nil, err
	}
	s.err = io.EOF
	return nil, io.EOF
}

// parse fills s.rec from a single trimmed data line.
func (s *fastBEDSource) parse(line []byte) error {
	// Reset the reused record's variable fields. Block columns are cleared so a
	// non-BED12 record following a BED12 one does not retain stale blocks.
	r := &s.rec
	r.Strand = ""
	r.BlockCount = 0
	r.BlockSizes = nil
	r.BlockStarts = nil

	// Tokenize up to the columns we care about. Fields 0–5 are always needed;
	// fields 9–11 only under keepBlocks.
	var cols [12][]byte
	n := 0
	begin := 0
	for i := 0; i <= len(line); i++ {
		if i == len(line) || line[i] == '\t' {
			if n < len(cols) {
				cols[n] = line[begin:i]
			}
			n++
			begin = i + 1
		}
	}
	if n < 3 {
		return fmt.Errorf("BED record must have at least 3 fields, got %d", n)
	}
	start, err := parseInt(cols[1])
	if err != nil {
		return fmt.Errorf("invalid chromStart %s: %v", cols[1], err)
	}
	end, err := parseInt(cols[2])
	if err != nil {
		return fmt.Errorf("invalid chromEnd %s: %v", cols[2], err)
	}
	r.Chrom = s.ci.intern(cols[0])
	r.ChromStart = start
	r.ChromEnd = end
	if n > 5 {
		r.Strand = internStrand(cols[5])
	}
	if s.keepBlocks && n > 11 {
		bc, err := parseInt(cols[9])
		if err != nil {
			return fmt.Errorf("invalid blockCount %s: %v", cols[9], err)
		}
		r.BlockCount = bc
		sizes, err := parseIntCSV(cols[10])
		if err != nil {
			return fmt.Errorf("invalid block size: %v", err)
		}
		r.BlockSizes = sizes
		starts, err := parseIntCSV(cols[11])
		if err != nil {
			return fmt.Errorf("invalid block start: %v", err)
		}
		r.BlockStarts = starts
	}
	return nil
}

// parseInt parses a base-10 integer from a byte slice without allocating on the
// common plain-digit fast path. It rejects anything the shared reader's
// strconv.Atoi would reject.
func parseInt(b []byte) (int, error) {
	if len(b) == 0 {
		return 0, fmt.Errorf("invalid syntax")
	}
	i := 0
	neg := false
	if b[0] == '+' || b[0] == '-' {
		neg = b[0] == '-'
		i++
	}
	if i == len(b) {
		return 0, fmt.Errorf("invalid syntax")
	}
	val := 0
	for ; i < len(b); i++ {
		c := b[i]
		if c < '0' || c > '9' {
			// Fall back to strconv for exotic forms so behaviour matches Atoi.
			return strconv.Atoi(string(b))
		}
		val = val*10 + int(c-'0')
	}
	if neg {
		val = -val
	}
	return val, nil
}

// parseIntCSV parses a comma-separated list of integers (with an optional
// trailing comma), mirroring the shared reader's block-column handling.
func parseIntCSV(b []byte) ([]int, error) {
	// Drop a single trailing comma, matching strings.TrimSuffix(...,",").
	if len(b) > 0 && b[len(b)-1] == ',' {
		b = b[:len(b)-1]
	}
	// Count fields.
	count := 1
	for _, c := range b {
		if c == ',' {
			count++
		}
	}
	out := make([]int, 0, count)
	begin := 0
	for i := 0; i <= len(b); i++ {
		if i == len(b) || b[i] == ',' {
			v, err := parseInt(b[begin:i])
			if err != nil {
				return nil, fmt.Errorf("invalid number %s: %v", b[begin:i], err)
			}
			out = append(out, v)
			begin = i + 1
		}
	}
	return out, nil
}

// trimSpaceBytes trims leading/trailing ASCII whitespace from b, returning a
// subslice of b (no allocation). It matches the set strings.TrimSpace treats as
// space for the ASCII inputs BED files use.
func trimSpaceBytes(b []byte) []byte {
	start := 0
	for start < len(b) && asciiSpace(b[start]) {
		start++
	}
	end := len(b)
	for end > start && asciiSpace(b[end-1]) {
		end--
	}
	return b[start:end]
}

func asciiSpace(c byte) bool {
	return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r'
}

// hasBytePrefix reports whether b begins with prefix.
func hasBytePrefix(b []byte, prefix string) bool {
	if len(b) < len(prefix) {
		return false
	}
	return string(b[:len(prefix)]) == prefix
}
