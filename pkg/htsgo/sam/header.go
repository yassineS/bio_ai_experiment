package sam

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// HeaderLine is a single parsed header line. Tag is the two-letter line type
// without the leading '@' (HD, SQ, RG, PG, or any user-defined record). Fields
// preserves the order of tag/value pairs as they appeared in the input so that
// callers can round-trip headers byte-for-byte where possible. Comment lines
// ('@CO') have a single field with key "" and the comment text as value.
type HeaderLine struct {
	Tag    string
	Fields []HeaderField
}

// HeaderField is one TAG:VALUE pair on a header line. For comment lines the
// Tag is empty and Value carries the raw text after "@CO\t".
type HeaderField struct {
	Tag   string
	Value string
}

// Get returns the value for the first field with the given two-letter tag,
// and reports whether the tag was present.
func (h HeaderLine) Get(tag string) (string, bool) {
	for _, f := range h.Fields {
		if f.Tag == tag {
			return f.Value, true
		}
	}
	return "", false
}

// Reference describes one @SQ entry: a reference sequence name and length.
type Reference struct {
	Name   string
	Length int32
	// Extra holds any additional @SQ TAG:VALUE pairs beyond SN/LN, preserved
	// for header round-tripping.
	Extra []HeaderField
}

// ReadGroup describes one @RG entry. ID is the mandatory ID field; the rest
// of the line is kept verbatim in Extra so callers can re-emit unchanged.
type ReadGroup struct {
	ID    string
	Extra []HeaderField
}

// Program describes one @PG entry.
type Program struct {
	ID    string
	Extra []HeaderField
}

// Header is the parsed header of a SAM or BAM file.
//
// Lines holds every header line in input order so SAM round-trips preserve
// byte-for-byte ordering. The convenience slices (Refs, ReadGroups, Programs,
// Comments) are derived views populated alongside Lines for fast lookup.
type Header struct {
	// Lines is the raw, ordered sequence of header lines.
	Lines []HeaderLine
	// Refs is the ordered list of @SQ entries.
	Refs []Reference
	// ReadGroups is the ordered list of @RG entries.
	ReadGroups []ReadGroup
	// Programs is the ordered list of @PG entries.
	Programs []Program
	// Comments is the ordered list of @CO comment texts.
	Comments []string
	// HDFields holds the @HD line fields (if any) for direct access.
	HDFields []HeaderField
}

// ErrInvalidHeader is returned when a header line cannot be parsed.
var ErrInvalidHeader = errors.New("sam: invalid header line")

// ParseHeader reads SAM header lines from r until a non-header line (or EOF).
// It consumes the trailing newline of each parsed header line but does not
// consume past the first body record; callers using a bufio.Reader can peek
// to detect when the body starts. The returned firstBodyLine, if non-empty,
// is the first record line that ParseHeader has already read past — callers
// must process it before reading further from r.
func ParseHeader(r *bufio.Reader) (*Header, string, error) {
	h := &Header{}
	for {
		// Peek the leading byte without consuming, so a record line stays in
		// the buffer for the body parser.
		b, err := r.Peek(1)
		if err != nil {
			if err == io.EOF {
				return h, "", nil
			}
			return nil, "", err
		}
		if b[0] != '@' {
			return h, "", nil
		}
		line, err := r.ReadString('\n')
		if err != nil && err != io.EOF {
			return nil, "", err
		}
		line = strings.TrimRight(line, "\r\n")
		if line == "" {
			continue
		}
		hl, perr := parseHeaderLine(line)
		if perr != nil {
			return nil, "", perr
		}
		h.appendLine(hl)
	}
}

// appendLine records a parsed header line in both Lines and the typed slices.
func (h *Header) appendLine(hl HeaderLine) {
	h.Lines = append(h.Lines, hl)
	switch hl.Tag {
	case "HD":
		h.HDFields = hl.Fields
	case "SQ":
		ref := Reference{}
		for _, f := range hl.Fields {
			switch f.Tag {
			case "SN":
				ref.Name = f.Value
			case "LN":
				n, _ := strconv.ParseInt(f.Value, 10, 32)
				ref.Length = int32(n)
			default:
				ref.Extra = append(ref.Extra, f)
			}
		}
		h.Refs = append(h.Refs, ref)
	case "RG":
		rg := ReadGroup{}
		for _, f := range hl.Fields {
			if f.Tag == "ID" {
				rg.ID = f.Value
			} else {
				rg.Extra = append(rg.Extra, f)
			}
		}
		h.ReadGroups = append(h.ReadGroups, rg)
	case "PG":
		pg := Program{}
		for _, f := range hl.Fields {
			if f.Tag == "ID" {
				pg.ID = f.Value
			} else {
				pg.Extra = append(pg.Extra, f)
			}
		}
		h.Programs = append(h.Programs, pg)
	case "CO":
		if len(hl.Fields) > 0 {
			h.Comments = append(h.Comments, hl.Fields[0].Value)
		}
	}
}

// parseHeaderLine parses one "@XY\t..." text header line.
func parseHeaderLine(line string) (HeaderLine, error) {
	if len(line) < 3 || line[0] != '@' {
		return HeaderLine{}, fmt.Errorf("%w: %q", ErrInvalidHeader, line)
	}
	tag := line[1:3]
	rest := ""
	if len(line) > 3 {
		if line[3] != '\t' {
			return HeaderLine{}, fmt.Errorf("%w: missing tab after @%s", ErrInvalidHeader, tag)
		}
		rest = line[4:]
	}
	hl := HeaderLine{Tag: tag}
	if tag == "CO" {
		// Comment line: the whole remainder is the comment, no TAG:VALUE split.
		hl.Fields = []HeaderField{{Tag: "", Value: rest}}
		return hl, nil
	}
	if rest == "" {
		return hl, nil
	}
	for _, part := range strings.Split(rest, "\t") {
		if len(part) < 3 || part[2] != ':' {
			return HeaderLine{}, fmt.Errorf("%w: bad field %q", ErrInvalidHeader, part)
		}
		hl.Fields = append(hl.Fields, HeaderField{
			Tag:   part[:2],
			Value: part[3:],
		})
	}
	return hl, nil
}

// WriteTo serialises the header as text SAM. Each line is emitted in the order
// it was stored, with TAG:VALUE pairs joined by tabs and terminated by '\n'.
func (h *Header) WriteTo(w io.Writer) (int64, error) {
	bw := bufio.NewWriter(w)
	var total int64
	for _, line := range h.Lines {
		n, err := bw.WriteString("@" + line.Tag)
		total += int64(n)
		if err != nil {
			return total, err
		}
		if line.Tag == "CO" {
			if len(line.Fields) > 0 {
				n, err := bw.WriteString("\t" + line.Fields[0].Value)
				total += int64(n)
				if err != nil {
					return total, err
				}
			}
		} else {
			for _, f := range line.Fields {
				n, err := bw.WriteString("\t" + f.Tag + ":" + f.Value)
				total += int64(n)
				if err != nil {
					return total, err
				}
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return total, err
		}
		total++
	}
	return total, bw.Flush()
}

// Text returns the SAM-encoded header as a string, preserving the verbatim
// input order of the lines. Use TextCanonical for htslib's grouped emission
// order.
func (h *Header) Text() string {
	var sb strings.Builder
	for _, line := range h.Lines {
		writeHeaderLine(&sb, line)
	}
	return sb.String()
}

// TextCanonical returns the SAM-encoded header with @-lines grouped into
// htslib's canonical emission order rather than verbatim input order: the @HD
// line(s) first, then @CO comment lines, then @PG lines, then @RG lines, and
// finally @SQ lines, with every other (user-defined) line type appended last.
// Within each group the original input order is preserved, so @SQ reference
// order and @PG / @RG order are unchanged.
//
// This mirrors htslib's header rebuild (sam_hrecs_rebuild_text), which is the
// order samtools emits into BAM and CRAM headers. Text() keeps the verbatim
// input order for byte-faithful SAM round-tripping; TextCanonical() is used
// where a container format (CRAM) needs to byte-match upstream's reordered
// header.
func (h *Header) TextCanonical() string {
	// htslib's canonical grouping order. Any line type not listed here is
	// emitted after these groups, in input order, so nothing is dropped.
	order := []string{"HD", "CO", "PG", "RG", "SQ"}
	rank := make(map[string]int, len(order))
	for i, t := range order {
		rank[t] = i
	}
	// Stable-bucket the lines by type while preserving per-group input order.
	groups := make([][]HeaderLine, len(order)+1)
	for _, line := range h.Lines {
		if r, ok := rank[line.Tag]; ok {
			groups[r] = append(groups[r], line)
		} else {
			groups[len(order)] = append(groups[len(order)], line)
		}
	}
	var sb strings.Builder
	for _, group := range groups {
		for _, line := range group {
			writeHeaderLine(&sb, line)
		}
	}
	return sb.String()
}

// writeHeaderLine appends one serialised header line (including the trailing
// newline) to sb. It is shared by Text and TextCanonical so the two stay in
// lock-step on the per-line encoding and differ only in line ordering.
func writeHeaderLine(sb *strings.Builder, line HeaderLine) {
	sb.WriteByte('@')
	sb.WriteString(line.Tag)
	if line.Tag == "CO" {
		if len(line.Fields) > 0 {
			sb.WriteByte('\t')
			sb.WriteString(line.Fields[0].Value)
		}
	} else {
		for _, f := range line.Fields {
			sb.WriteByte('\t')
			sb.WriteString(f.Tag)
			sb.WriteByte(':')
			sb.WriteString(f.Value)
		}
	}
	sb.WriteByte('\n')
}

// RefIndex returns the position of the named reference in h.Refs, or -1 if
// the reference is unknown.
func (h *Header) RefIndex(name string) int {
	for i, r := range h.Refs {
		if r.Name == name {
			return i
		}
	}
	return -1
}

// ParseHeaderText parses a header that has already been read into memory.
// It is a thin wrapper around ParseHeader and is convenient for BAM, which
// loads the header text in one shot.
func ParseHeaderText(text string) (*Header, error) {
	r := bufio.NewReader(strings.NewReader(text))
	h, _, err := ParseHeader(r)
	return h, err
}
