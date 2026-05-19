// Package gff provides a minimal GFF3 parser used by tools like
// bcftools csq for protein-coding consequence prediction.
//
// V1 SCOPE: only the columns needed to walk SNPs through CDS regions
// are parsed. The implementation is intentionally tiny — it understands
// the GFF3 9-column tab format, the `key=value;...` attributes, and the
// feature types `gene`, `mRNA` / `transcript`, `CDS`, `exon`. Anything
// else is skipped. This matches the v1 csq simplification documented in
// docs/PARITY_ROADMAP.md.
//
// Format reference: http://gmod.org/wiki/GFF3
package gff

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Strand encodes the GFF3 strand column.
type Strand byte

const (
	// StrandUnknown is "." or any unrecognised value.
	StrandUnknown Strand = '.'
	// StrandForward is "+".
	StrandForward Strand = '+'
	// StrandReverse is "-".
	StrandReverse Strand = '-'
)

// Feature is one parsed GFF3 record.
type Feature struct {
	Seqid      string            // CHROM / contig
	Source     string            // tool that produced the entry
	Type       string            // feature type (gene, mRNA, CDS, exon, ...)
	Start      int               // 1-based inclusive start
	End        int               // 1-based inclusive end
	Score      string            // score column, kept as string ('.' for missing)
	Strand     Strand            // strand
	Phase      int               // 0/1/2 for CDS, -1 when missing
	Attributes map[string]string // key=value pairs from column 9
}

// ID returns the GFF3 `ID=` attribute, or "" if absent.
func (f *Feature) ID() string {
	return f.Attributes["ID"]
}

// Parent returns the GFF3 `Parent=` attribute, or "" if absent.
// Note: GFF3 allows Parent to be a comma-list; we return the raw
// string. Callers that care about multi-parent records should split.
func (f *Feature) Parent() string {
	return f.Attributes["Parent"]
}

// Reader iterates over GFF3 features.
type Reader struct {
	scanner *bufio.Scanner
	err     error
	line    int
}

// NewReader constructs a Reader from r.
func NewReader(r io.Reader) *Reader {
	scanner := bufio.NewScanner(r)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)
	return &Reader{scanner: scanner}
}

// Read returns the next Feature, or io.EOF when the stream is done.
// Comment lines (`#...`) and blank lines are skipped silently.
func (r *Reader) Read() (*Feature, error) {
	if r.err != nil {
		return nil, r.err
	}
	for r.scanner.Scan() {
		r.line++
		line := r.scanner.Text()
		if line == "" || line[0] == '#' {
			continue
		}
		f, err := parseLine(line)
		if err != nil {
			return nil, fmt.Errorf("gff: line %d: %w", r.line, err)
		}
		return f, nil
	}
	if err := r.scanner.Err(); err != nil {
		r.err = err
		return nil, err
	}
	r.err = io.EOF
	return nil, io.EOF
}

// ReadAll consumes the reader and returns every parsed Feature.
func (r *Reader) ReadAll() ([]*Feature, error) {
	out := []*Feature{}
	for {
		f, err := r.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, f)
	}
}

// parseLine parses a single GFF3 record (9 tab-separated columns).
func parseLine(line string) (*Feature, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 9 {
		return nil, fmt.Errorf("expected 9 tab-separated columns, got %d", len(fields))
	}
	start, err := strconv.Atoi(fields[3])
	if err != nil {
		return nil, fmt.Errorf("bad start %q: %v", fields[3], err)
	}
	end, err := strconv.Atoi(fields[4])
	if err != nil {
		return nil, fmt.Errorf("bad end %q: %v", fields[4], err)
	}
	strand := StrandUnknown
	if len(fields[6]) > 0 {
		switch fields[6][0] {
		case '+':
			strand = StrandForward
		case '-':
			strand = StrandReverse
		}
	}
	phase := -1
	if fields[7] != "." && fields[7] != "" {
		p, err := strconv.Atoi(fields[7])
		if err != nil {
			return nil, fmt.Errorf("bad phase %q: %v", fields[7], err)
		}
		phase = p
	}
	return &Feature{
		Seqid:      fields[0],
		Source:     fields[1],
		Type:       fields[2],
		Start:      start,
		End:        end,
		Score:      fields[5],
		Strand:     strand,
		Phase:      phase,
		Attributes: parseAttributes(fields[8]),
	}, nil
}

// parseAttributes parses GFF3's column-9 key=value;key=value layout.
// Values are URL-decoded only minimally; we keep them as-is and trim
// surrounding whitespace. Empty input returns a non-nil empty map.
func parseAttributes(s string) map[string]string {
	out := make(map[string]string)
	if s == "" || s == "." {
		return out
	}
	for _, kv := range strings.Split(s, ";") {
		kv = strings.TrimSpace(kv)
		if kv == "" {
			continue
		}
		eq := strings.IndexByte(kv, '=')
		if eq < 0 {
			// Tolerate the older GFF2/GTF-ish `key "value"` form by
			// turning whitespace into the separator.
			sp := strings.IndexByte(kv, ' ')
			if sp < 0 {
				out[kv] = ""
				continue
			}
			out[kv[:sp]] = strings.Trim(kv[sp+1:], "\" ")
			continue
		}
		out[kv[:eq]] = kv[eq+1:]
	}
	return out
}
