package bcf

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Magic is the five-byte BCF v2.2 signature.
var Magic = [5]byte{'B', 'C', 'F', 0x02, 0x02}

// Errors surfaced by the BCF header parser.
var (
	ErrBadMagic  = errors.New("bcf: bad magic (expected BCF\\2\\2)")
	ErrTruncated = errors.New("bcf: truncated input")
)

// DictKind identifies which of the three BCF dictionaries an entry lives in.
type DictKind int

const (
	// DictContig holds CHROM names from ##contig lines.
	DictContig DictKind = iota
	// DictTagInfo holds INFO+FILTER tag names from ##INFO/##FILTER lines.
	DictTagInfo
	// DictTagFmt holds FORMAT tag names from ##FORMAT lines.
	DictTagFmt
)

// DictEntry describes one entry in a BCF dictionary. Number and Type come
// from the VCF meta-information line and are exposed so that the record
// decoder can pretty-print values back into VCF text.
type DictEntry struct {
	ID     string // tag identifier (CHROM name or INFO/FORMAT key)
	Number string // VCF "Number=" attribute (".", "A", "G", "R", or an integer)
	Type   string // VCF "Type=" attribute (Integer / Float / String / Flag / Character)
}

// Header is the parsed BCF header. Text is the verbatim VCF-style text
// header (without the trailing NUL byte); the typed slices give the
// dictionaries in the order htslib uses them on the wire.
type Header struct {
	Text     string
	VCF      *vcf.Header // a vcf.Header constructed from Text for convenience
	Contigs  []DictEntry // index = the int referred to by record CHROM
	InfoTags []DictEntry // dictionary for INFO keys *and* FILTER names
	FmtTags  []DictEntry // dictionary for FORMAT keys
	Samples  []string    // sample names from the #CHROM line

	// tail carries the buffered reader positioned at the first record byte
	// so Reader.Read can resume from where ReadHeader stopped.
	tail *bufio.Reader
}

// ContigName returns the n-th CHROM dictionary entry's ID, or "" if n is out
// of range. -1 (a wire "no chrom" marker) maps to "".
func (h *Header) ContigName(n int32) string {
	if n < 0 || int(n) >= len(h.Contigs) {
		return ""
	}
	return h.Contigs[n].ID
}

// InfoTag returns the n-th INFO/FILTER dictionary entry, or nil if n is out
// of range.
func (h *Header) InfoTag(n int32) *DictEntry {
	if n < 0 || int(n) >= len(h.InfoTags) {
		return nil
	}
	return &h.InfoTags[n]
}

// FmtTag returns the n-th FORMAT dictionary entry, or nil if n is out of
// range.
func (h *Header) FmtTag(n int32) *DictEntry {
	if n < 0 || int(n) >= len(h.FmtTags) {
		return nil
	}
	return &h.FmtTags[n]
}

// ReadHeader decodes a BCF magic + text header from r. It returns a parsed
// Header (text, dictionaries, samples) and leaves r positioned at the start
// of the first record.
//
// Callers must pass an already-BGZF-decompressed stream — typically the
// io.ReadCloser returned by pkg/bioformats/iohelper.OpenReader.
func ReadHeader(r io.Reader) (*Header, error) {
	br := bufio.NewReader(r)

	var magic [5]byte
	if _, err := io.ReadFull(br, magic[:]); err != nil {
		if err == io.EOF || err == io.ErrUnexpectedEOF {
			return nil, ErrTruncated
		}
		return nil, err
	}
	if magic != Magic {
		return nil, fmt.Errorf("%w: got %v", ErrBadMagic, magic)
	}

	var lText uint32
	if err := binary.Read(br, binary.LittleEndian, &lText); err != nil {
		return nil, wrapEOF(err)
	}
	if lText == 0 {
		return nil, fmt.Errorf("%w: empty header", ErrTruncated)
	}
	textBuf := make([]byte, lText)
	if _, err := io.ReadFull(br, textBuf); err != nil {
		return nil, wrapEOF(err)
	}
	// The header is NUL-terminated; trim trailing NULs so downstream parsers
	// don't choke on them.
	text := strings.TrimRight(string(textBuf), "\x00")

	hdr, err := parseTextHeader(text)
	if err != nil {
		return nil, err
	}

	// We may have peeked further into the stream via bufio — splice the
	// remaining buffered bytes back onto r so the record reader sees them.
	hdr.tail = br
	return hdr, nil
}

// parseTextHeader walks the VCF-style header lines and builds the three
// dictionaries plus the sample list. It is intentionally tolerant of out-of-
// order ##contig / ##INFO / ##FORMAT / ##FILTER lines — htslib emits them
// in declaration order, but the BCF spec only requires the dictionaries to
// be addressable by index.
func parseTextHeader(text string) (*Header, error) {
	h := &Header{Text: text, VCF: &vcf.Header{}}

	// The implicit PASS filter is always entry 0 of the INFO/FILTER dict
	// in htslib — register it before scanning any ##FILTER lines.
	h.InfoTags = append(h.InfoTags, DictEntry{ID: "PASS", Type: "Flag", Number: "0"})

	lines := strings.Split(text, "\n")
	for _, line := range lines {
		if line == "" {
			continue
		}
		switch {
		case strings.HasPrefix(line, "##contig="):
			entry := parseStructured(line[len("##contig="):])
			if entry.ID == "" {
				continue
			}
			h.Contigs = append(h.Contigs, entry)
			h.VCF.MetaInfo = append(h.VCF.MetaInfo, stripIDXAnnotation(line))
		case strings.HasPrefix(line, "##INFO="):
			entry := parseStructured(line[len("##INFO="):])
			if entry.ID != "" {
				h.InfoTags = append(h.InfoTags, entry)
			}
			h.VCF.MetaInfo = append(h.VCF.MetaInfo, stripIDXAnnotation(line))
		case strings.HasPrefix(line, "##FILTER="):
			entry := parseStructured(line[len("##FILTER="):])
			if entry.ID != "" && entry.ID != "PASS" {
				h.InfoTags = append(h.InfoTags, entry)
			}
			h.VCF.MetaInfo = append(h.VCF.MetaInfo, stripIDXAnnotation(line))
		case strings.HasPrefix(line, "##FORMAT="):
			entry := parseStructured(line[len("##FORMAT="):])
			if entry.ID != "" {
				h.FmtTags = append(h.FmtTags, entry)
			}
			h.VCF.MetaInfo = append(h.VCF.MetaInfo, stripIDXAnnotation(line))
		case strings.HasPrefix(line, "#CHROM"):
			fields := strings.Split(line, "\t")
			if len(fields) > 9 {
				h.Samples = append(h.Samples, fields[9:]...)
				h.VCF.Samples = append(h.VCF.Samples, fields[9:]...)
			}
		default:
			if strings.HasPrefix(line, "##") {
				h.VCF.MetaInfo = append(h.VCF.MetaInfo, stripIDXAnnotation(line))
			}
		}
	}
	return h, nil
}

// stripIDXAnnotation removes the htslib-private `,IDX=N>` suffix that
// htslib's bcf_hdr_format adds to ##INFO / ##FILTER / ##FORMAT / ##contig
// lines when reading a BCF file. The BCF text view emitted by upstream
// bcftools strips this annotation before printing, so we do the same to
// match byte-for-byte. Lines without an IDX= attribute pass through
// unchanged.
func stripIDXAnnotation(line string) string {
	// We only need to handle structured lines closed by ">" with a
	// possible trailing `,IDX=NNN` before the close.
	end := strings.LastIndexByte(line, '>')
	if end < 0 {
		return line
	}
	// Find a `,IDX=` that ends at `end`.
	idx := strings.LastIndex(line[:end], ",IDX=")
	if idx < 0 {
		return line
	}
	// Verify everything after IDX= up to end is decimal digits.
	for i := idx + len(",IDX="); i < end; i++ {
		if line[i] < '0' || line[i] > '9' {
			return line
		}
	}
	return line[:idx] + line[end:]
}

// parseStructured extracts the ID, Number, and Type attributes from a VCF
// structured meta-information value of the form "<ID=foo,Number=1,Type=Integer,...>".
// Unknown / missing attributes default to the empty string.
func parseStructured(s string) DictEntry {
	s = strings.TrimSpace(s)
	s = strings.TrimPrefix(s, "<")
	s = strings.TrimSuffix(s, ">")
	var entry DictEntry
	for _, kv := range splitStructured(s) {
		eq := strings.IndexByte(kv, '=')
		if eq < 0 {
			continue
		}
		key := strings.TrimSpace(kv[:eq])
		val := strings.TrimSpace(kv[eq+1:])
		val = strings.Trim(val, `"`)
		switch key {
		case "ID":
			entry.ID = val
		case "Number":
			entry.Number = val
		case "Type":
			entry.Type = val
		}
	}
	return entry
}

// splitStructured splits the inner part of a "<...>" attribute list on commas
// that are not inside double quotes. The VCF spec allows quoted commas inside
// Description="..." so a naive strings.Split won't do.
func splitStructured(s string) []string {
	var out []string
	depth := 0
	inQuote := false
	start := 0
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch {
		case c == '"':
			inQuote = !inQuote
		case c == '<' && !inQuote:
			depth++
		case c == '>' && !inQuote:
			depth--
		case c == ',' && !inQuote && depth == 0:
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	out = append(out, s[start:])
	return out
}

// tailReader returns the buffered reader positioned at the first record
// byte. It is used by Reader to continue from where ReadHeader stopped.
func (h *Header) tailReader() *bufio.Reader { return h.tail }

// wrapEOF turns io.EOF / io.ErrUnexpectedEOF into the package's ErrTruncated
// sentinel so callers can distinguish "ran out of bytes" from real I/O errors.
func wrapEOF(err error) error {
	if err == io.EOF || err == io.ErrUnexpectedEOF {
		return ErrTruncated
	}
	return err
}
