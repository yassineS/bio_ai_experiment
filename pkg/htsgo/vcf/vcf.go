// Package vcf provides utilities for reading and writing VCF (Variant Call Format) files.
// VCF is a text file format for storing gene sequence variations.
//
// Format specification follows VCF v4.2 and later:
//   - Meta-information lines start with '##'
//   - Header line starts with '#CHROM'
//   - Data lines contain tab-separated fields
package vcf

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
	"unsafe"
)

// Variant represents a single VCF variant record.
type Variant struct {
	Chrom     string            // Chromosome
	Pos       int               // Position (1-based)
	ID        string            // Variant ID (or '.')
	Ref       string            // Reference allele
	Alt       []string          // Alternate alleles
	Qual      float64           // Quality score (or -1 if missing)
	Filter    []string          // Filter status (PASS or filter names)
	Info      map[string]string // INFO field key-value pairs
	InfoOrder []string          // INFO key insertion order (preserved for byte-for-byte parity with upstream)
	Format    []string          // FORMAT field tags
	Samples   []Sample          // Sample genotype data
	// RawTail, when non-empty, is the verbatim FORMAT + sample columns (the
	// line from the start of the FORMAT field, tab-separated). A reader in
	// KeepRawSamples mode fills it instead of parsing Format/Samples, and Write
	// re-emits it verbatim — byte-identical to parsing the columns into per-
	// sample maps and re-serialising them, for a well-formed record. It lets a
	// caller that never inspects sample data (e.g. bcftools isec) skip the
	// expensive map round-trip. Empty for normally-parsed records.
	RawTail string
	// RawLine, when non-empty, is the verbatim data line (no trailing newline) a
	// reader in KeepRawLine mode captured. In that mode only the sort-key columns
	// (CHROM/POS/REF/ALT) are parsed; the INFO map, FILTER, QUAL and per-sample
	// data are NOT populated. A caller that only needs to re-order records and
	// re-emit them unchanged (e.g. bcftools sort -O v) keys the sort on the parsed
	// columns and writes RawLine back verbatim, avoiding the full parse→re-encode
	// round-trip and the per-record maps it allocates. Empty for normally-parsed
	// records.
	RawLine string
}

// Sample represents genotype data for a single sample.
type Sample struct {
	Name string            // Sample name
	Data map[string]string // FORMAT field values
}

// Header represents VCF metadata and header information.
type Header struct {
	MetaInfo []string // Meta-information lines (##)
	Samples  []string // Sample names from header line
}

// Reader provides sequential access to VCF records.
type Reader struct {
	scanner *bufio.Scanner
	header  *Header
	err     error
	// fieldsBuf and valuesBuf are reused tab/colon split scratch slices, so a
	// record parse does not allocate them anew. They hold substrings of the
	// current line and are overwritten on the next Read, so callers must not
	// retain them.
	fieldsBuf []string
	valuesBuf []string
	// lineBuf is a reused copy of the current line's raw bytes. bufio.Scanner's
	// Bytes() return is only valid until the next Scan, so the read path copies
	// it here once per line (amortized zero-alloc after warmup) and parses from
	// the copy. The ReadInto fast path aliases lineBuf directly as a string via
	// unsafe.String (no per-line alloc); the retaining Read path materializes a
	// fresh string(lineBuf) so the returned Variant owns its own backing.
	lineBuf []byte
	// keepRawTail, set by KeepRawSamples, makes parseLine keep the FORMAT +
	// sample columns verbatim in Variant.RawTail instead of building per-sample
	// maps. See KeepRawSamples.
	keepRawTail bool
	// keepRawLine, set by KeepRawLine, makes parseLine capture the whole verbatim
	// data line in Variant.RawLine and parse only the sort-key columns
	// (CHROM/POS/REF/ALT). See KeepRawLine.
	keepRawLine bool
	// intern is a small string-interning table used only on the ReadInto path.
	// ReadInto's Chrom/Ref/Alt (and Filter) strings alias the reused line buffer
	// and become invalid on the next read, so a caller that retains them across
	// reads (e.g.
	// the vcftools streaming-stats accumulators, which store v.Chrom / the allele
	// strings per site) would otherwise hold dangling aliases. Interning collapses
	// each distinct value to a single stable, owned string: the whole-file set of
	// chromosome names collapses to a handful of entries, and the common short
	// alleles (A/C/G/T, ".") to a fixed few, so retained copies stay valid without
	// per-site allocation. The Read (owned) path already materialises fresh strings
	// and does not use this. See internStr.
	intern map[string]string
}

// internStr returns an owned, stable copy of s, reusing a single backing string
// per distinct value via the reader's intern table. The argument s aliases the
// reused line buffer (ReadInto path); the returned string owns its own memory and
// is safe to retain across subsequent reads. The table is bounded in practice by
// the small cardinality of the interned fields (chromosome names and short
// alleles), so it does not grow with the number of records.
func (r *Reader) internStr(s string) string {
	if v, ok := r.intern[s]; ok {
		return v
	}
	// string(...) forces a fresh allocation that owns its bytes, breaking the
	// alias into r.lineBuf before we store it as both key and value.
	owned := string([]byte(s))
	if r.intern == nil {
		r.intern = make(map[string]string)
	}
	r.intern[owned] = owned
	return owned
}

// internSlice interns each element of s in place using the supplied intern
// function, so a slice of aliasing strings (e.g. v.Alt on the ReadInto path)
// becomes a slice of stable owned strings that survive the next read.
func internSlice(intern func(string) string, s []string) {
	for i := range s {
		s[i] = intern(s[i])
	}
}

// KeepRawSamples toggles "shallow sample" parsing: when on, Read/ReadInto parse
// CHROM..INFO normally but leave the FORMAT + sample columns unparsed, exposing
// them verbatim as Variant.RawTail (Format/Samples stay empty); Write re-emits
// RawTail unchanged. For a well-formed record this is byte-identical to parsing
// the columns into per-sample maps and re-serialising them, but skips that whole
// round-trip — a large win for callers that only need CHROM/POS/REF/ALT/ID and
// re-emit records unchanged (e.g. bcftools isec). Because RawTail aliases the
// parsed line, only the retaining Read path is safe to keep across reads; do not
// retain a ReadInto Variant's RawTail.
func (r *Reader) KeepRawSamples(on bool) { r.keepRawTail = on }

// KeepRawLine toggles "raw line" parsing: when on, Read/ReadInto capture the
// whole verbatim data line in Variant.RawLine and parse only the sort-key
// columns (CHROM, POS, REF, ALT). The INFO/FILTER/QUAL fields and the per-sample
// data are left unparsed (Info empty, Qual 0, Filter/Format/Samples empty). It is
// for callers that re-order records and re-emit them unchanged — bcftools sort
// -O v — where re-emitting the captured line is byte-identical to a full
// parse→re-encode for a well-formed record while avoiding that round-trip and its
// per-record allocations. Because RawLine (and the parsed key columns) alias the
// reused line buffer on the ReadInto path, only the retaining Read path is safe
// to keep across reads; a ReadInto caller that buffers records must copy RawLine.
func (r *Reader) KeepRawLine(on bool) { r.keepRawLine = on }

// NewReader creates a new VCF reader from an io.Reader.
func NewReader(r io.Reader) *Reader {
	scanner := bufio.NewScanner(r)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024) // 10MB max token size
	return &Reader{
		scanner: scanner,
		header:  &Header{},
	}
}

// ReadHeader reads and parses the VCF header.
// Must be called before reading variants.
func (r *Reader) ReadHeader() (*Header, error) {
	for r.scanner.Scan() {
		// Header lines are retained (in MetaInfo / Samples), so Text() — which
		// copies the scanner's transient bytes into an owned string — is the
		// right call here. The header is a handful of lines, so its allocation
		// is negligible; the hot per-record path below is the one optimized to
		// avoid per-line strings.
		line := r.scanner.Text()

		// Meta-information line
		if strings.HasPrefix(line, "##") {
			r.header.MetaInfo = append(r.header.MetaInfo, line)
			continue
		}

		// Header line
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) < 8 {
				return nil, fmt.Errorf("invalid VCF header line: %s", line)
			}
			// Sample names start at column 9 (index 9)
			if len(fields) > 9 {
				r.header.Samples = fields[9:]
			}
			return r.header, nil
		}

		// Unexpected line before header
		return nil, fmt.Errorf("unexpected line before VCF header: %s", line)
	}

	if err := r.scanner.Err(); err != nil {
		return nil, err
	}

	return nil, fmt.Errorf("no VCF header found")
}

// Read reads the next variant record.
// Returns io.EOF when no more records are available.
//
// The returned Variant owns all its string fields (they are copied out of the
// reader's scratch buffer), so the caller may retain it across subsequent
// reads — unlike ReadInto, whose Variant must not be retained.
func (r *Reader) Read() (*Variant, error) {
	v := &Variant{}
	if err := r.readInto(v, true); err != nil {
		return nil, err
	}
	return v, nil
}

// ReadInto parses the next variant into v, reusing v's Info/Sample maps and its
// Alt/Filter/Format/InfoOrder/Samples slices instead of allocating fresh ones.
// It is for consume-and-discard scans (bcftools view, query): the caller must
// not retain v — or any string it exposes — across calls, since the next
// ReadInto overwrites them. It returns io.EOF at end of input, like Read.
func (r *Reader) ReadInto(v *Variant) error {
	return r.readInto(v, false)
}

// readInto advances to the next data line (skipping blanks and comments) and
// parses it into v. It is the shared body of Read and ReadInto.
//
// It reads the line via scanner.Bytes() (no per-line string allocation) into
// the reused r.lineBuf. When owned is true (the Read path) the line is
// materialized as a fresh string so v's fields own their backing memory and v
// may be retained. When owned is false (the ReadInto path) the line string
// aliases r.lineBuf with no allocation — valid only until the next read, which
// matches ReadInto's documented contract.
func (r *Reader) readInto(v *Variant, owned bool) error {
	if r.err != nil {
		return r.err
	}
	if r.header == nil || len(r.header.MetaInfo) == 0 {
		return fmt.Errorf("header not read; call ReadHeader() first")
	}
	for {
		if !r.scanner.Scan() {
			if err := r.scanner.Err(); err != nil {
				r.err = err
				return err
			}
			r.err = io.EOF
			return io.EOF
		}
		b := r.scanner.Bytes()
		if len(b) == 0 || b[0] == '#' {
			continue
		}
		// Copy the transient scanner bytes into the reused buffer; Bytes() is
		// only valid until the next Scan.
		r.lineBuf = append(r.lineBuf[:0], b...)
		var line string
		if owned {
			// Read: caller may retain v, so the string must own its memory.
			line = string(r.lineBuf)
		} else {
			// ReadInto: alias the buffer with no allocation. Safe because the
			// caller must not retain v (or its strings) past the next read.
			line = unsafe.String(unsafe.SliceData(r.lineBuf), len(r.lineBuf))
		}
		return r.parseLine(line, v, owned)
	}
}

// parseLine fills v from a single VCF data line, reusing v's existing maps and
// slices. On the Read (owned) path line is a freshly-allocated owned string, so
// all of v's fields own their backing. On the ReadInto path line aliases the
// reused buffer; parseLine then interns the fields a caller is documented as
// allowed to retain across reads (CHROM, REF, and each ALT allele) so those
// specific strings become stable owned copies. ID is not interned: no consumer
// retains it across sites and rsIDs are ~unique, so interning it would grow the
// table without bound for no benefit. The remaining aliasing fields (ID, INFO
// values, FILTER, FORMAT, per-sample data) keep ReadInto's discard-before-
// next-read contract.
func (r *Reader) parseLine(line string, v *Variant, owned bool) error {
	r.fieldsBuf = splitInto(r.fieldsBuf, line, '\t')
	fields := r.fieldsBuf
	if len(fields) < 8 {
		return fmt.Errorf("invalid VCF record: insufficient fields: %s", line)
	}
	pos, err := strconv.Atoi(fields[1])
	if err != nil {
		return fmt.Errorf("invalid position %s: %v", fields[1], err)
	}
	if r.keepRawLine {
		// Raw-line mode: parse only the sort-key columns and keep the whole line
		// verbatim; leave everything else unparsed/cleared.
		v.Chrom = fields[0]
		v.Pos = pos
		v.Ref = fields[3]
		v.Alt = splitInto(v.Alt, fields[4], ',')
		if !owned {
			v.Chrom = r.internStr(v.Chrom)
			v.Ref = r.internStr(v.Ref)
			internSlice(r.internStr, v.Alt)
		}
		v.RawLine = line
		v.ID = ""
		v.Qual = 0
		v.Filter = v.Filter[:0]
		v.InfoOrder = v.InfoOrder[:0]
		if v.Info != nil {
			clear(v.Info)
		}
		v.Format = v.Format[:0]
		v.Samples = v.Samples[:0]
		v.RawTail = ""
		return nil
	}
	v.Chrom = fields[0]
	v.Pos = pos
	v.ID = fields[2]
	v.Ref = fields[3]
	v.Alt = splitInto(v.Alt, fields[4], ',')
	if !owned {
		// ReadInto path: intern the fields a caller may retain across reads so
		// they become stable owned strings rather than aliases into r.lineBuf.
		// ID is deliberately not interned: no consumer retains v.ID by reference
		// across sites, and rsIDs are ~unique per site, so interning it would grow
		// the table O(n_sites) for no benefit.
		v.Chrom = r.internStr(v.Chrom)
		v.Ref = r.internStr(v.Ref)
		internSlice(r.internStr, v.Alt)
	}
	if fields[5] == "." {
		v.Qual = -1
	} else {
		qual, err := strconv.ParseFloat(fields[5], 64)
		if err != nil {
			return fmt.Errorf("invalid quality %s: %v", fields[5], err)
		}
		v.Qual = qual
	}
	v.Filter = parseFilterInto(v.Filter, fields[6])
	if !owned {
		// FILTER is retained across reads by the --FILTER-summary accumulator
		// (which keys a map on the joined FILTER string; a single-tag join
		// returns the aliasing element verbatim), so intern its tags too. FILTER
		// tags have tiny cardinality (PASS/./named filters), so this collapses to
		// a handful of stable strings.
		internSlice(r.internStr, v.Filter)
	}
	v.InfoOrder, v.Info = parseInfoInto(v.InfoOrder, v.Info, fields[7])

	if len(fields) > 8 && r.keepRawTail {
		// Shallow mode: keep the FORMAT + sample columns verbatim (the line from
		// the start of fields[8]) and skip the per-sample map round-trip.
		v.RawTail = rawTailOf(line)
		v.Format = v.Format[:0]
		v.Samples = v.Samples[:0]
		return nil
	}
	v.RawTail = ""
	if len(fields) > 8 {
		v.Format = splitInto(v.Format, fields[8], ':')
		n := len(r.header.Samples)
		if cap(v.Samples) >= n {
			v.Samples = v.Samples[:n]
		} else {
			v.Samples = make([]Sample, n)
		}
		for i := 0; i < n; i++ {
			s := &v.Samples[i]
			if 9+i >= len(fields) {
				// Fewer sample columns than header samples: match Read's
				// fresh make, which leaves a zero-value Sample here.
				s.Name, s.Data = "", nil
				continue
			}
			s.Name = r.header.Samples[i]
			if s.Data == nil {
				s.Data = make(map[string]string, len(v.Format))
			} else {
				clear(s.Data)
			}
			r.valuesBuf = splitInto(r.valuesBuf, fields[9+i], ':')
			for j, format := range v.Format {
				if j < len(r.valuesBuf) {
					s.Data[format] = r.valuesBuf[j]
				}
			}
		}
	} else {
		v.Format = v.Format[:0]
		v.Samples = v.Samples[:0]
	}
	return nil
}

// rawTailOf returns the substring of a VCF data line starting at the FORMAT
// column — i.e. everything after the 8th tab (CHROM POS ID REF ALT QUAL FILTER
// INFO are the first 8 fields). Returns "" if the line has no FORMAT column. The
// result aliases line, so it lives only as long as line.
func rawTailOf(line string) string {
	n := 0
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {
			if n == 7 {
				return line[i+1:]
			}
			n++
		}
	}
	return ""
}

// splitInto splits s on the single byte sep, appending the substrings into
// dst[:0] (reusing its backing array). It matches strings.Split for a one-byte
// separator. The substrings reference s, so they live only as long as s.
func splitInto(dst []string, s string, sep byte) []string {
	dst = dst[:0]
	for {
		i := strings.IndexByte(s, sep)
		if i < 0 {
			return append(dst, s)
		}
		dst = append(dst, s[:i])
		s = s[i+1:]
	}
}

// parseFilterInto is parseFilter writing into the reused slice dst.
func parseFilterInto(dst []string, filter string) []string {
	dst = dst[:0]
	if filter == "." || filter == "PASS" {
		return append(dst, filter)
	}
	return splitInto(dst, filter, ';')
}

// parseInfoInto is parseInfoWithOrder reusing the order slice and m map (m is
// cleared, or allocated when nil). It reproduces parseInfoWithOrder exactly:
// keys in first-seen order, "=" splits on the first '=', a bare key maps to "".
func parseInfoInto(order []string, m map[string]string, info string) ([]string, map[string]string) {
	if m == nil {
		m = make(map[string]string)
	} else {
		clear(m)
	}
	order = order[:0]
	if info == "." || info == "" {
		return order, m
	}
	rest := info
	for len(rest) > 0 {
		pair := rest
		if i := strings.IndexByte(rest, ';'); i >= 0 {
			pair, rest = rest[:i], rest[i+1:]
		} else {
			rest = ""
		}
		k, val, hasVal := pair, "", false
		if e := strings.IndexByte(pair, '='); e >= 0 {
			k, val, hasVal = pair[:e], pair[e+1:], true
		}
		if k == "" {
			continue
		}
		if _, seen := m[k]; !seen {
			order = append(order, k)
		}
		if hasVal {
			m[k] = val
		} else {
			m[k] = ""
		}
	}
	return order, m
}

// ReadAll reads all variant records from the reader.
func (r *Reader) ReadAll() ([]*Variant, error) {
	var variants []*Variant
	for {
		variant, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return variants, err
		}
		variants = append(variants, variant)
	}
	return variants, nil
}

// Writer provides sequential writing of VCF records.
type Writer struct {
	w        *bufio.Writer
	header   *Header
	num      []byte          // scratch for integer/float formatting, reused per record
	infoSeen map[string]bool // reused INFO key set for writeInfo's extras pass
}

// NewWriter creates a new VCF writer.
func NewWriter(w io.Writer, header *Header) *Writer {
	return &Writer{
		w:      bufio.NewWriter(w),
		header: header,
	}
}

// WriteHeader writes the VCF header.
// Must be called before writing variants.
func (w *Writer) WriteHeader() error {
	// Write meta-information
	for _, meta := range w.header.MetaInfo {
		if _, err := fmt.Fprintln(w.w, meta); err != nil {
			return err
		}
	}

	// Write header line
	headerLine := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
	if len(w.header.Samples) > 0 {
		headerLine += "\tFORMAT\t" + strings.Join(w.header.Samples, "\t")
	}
	if _, err := fmt.Fprintln(w.w, headerLine); err != nil {
		return err
	}

	return nil
}

// Write writes a variant record. The line is streamed field by field straight
// into the buffered writer — no intermediate per-field strings, join slices, or
// Sprintf — which is byte-identical to the join-based form but avoids the dozens
// of allocations per record that dominated VCF output.
func (w *Writer) Write(variant *Variant) error {
	bw := w.w
	if variant.RawLine != "" {
		// Raw-line record (KeepRawLine): the whole data line was captured
		// verbatim. Re-emitting it is byte-identical to parsing every column and
		// re-serialising, for a well-formed record.
		if _, err := bw.WriteString(variant.RawLine); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	bw.WriteString(variant.Chrom)
	bw.WriteByte('\t')
	w.num = strconv.AppendInt(w.num[:0], int64(variant.Pos), 10)
	bw.Write(w.num)
	bw.WriteByte('\t')
	bw.WriteString(variant.ID)
	bw.WriteByte('\t')
	bw.WriteString(variant.Ref)
	bw.WriteByte('\t')
	for i, a := range variant.Alt {
		if i > 0 {
			bw.WriteByte(',')
		}
		bw.WriteString(a)
	}
	bw.WriteByte('\t')

	// Quality
	if variant.Qual < 0 {
		bw.WriteByte('.')
	} else {
		bw.WriteString(formatQual(variant.Qual))
	}
	bw.WriteByte('\t')

	// Filter
	if len(variant.Filter) == 0 {
		bw.WriteByte('.')
	} else {
		for i, f := range variant.Filter {
			if i > 0 {
				bw.WriteByte(';')
			}
			bw.WriteString(f)
		}
	}
	bw.WriteByte('\t')

	// Info
	if !w.writeInfo(variant.Info, variant.InfoOrder) {
		bw.WriteByte('.')
	}

	// Format and samples
	if variant.RawTail != "" {
		// Shallow-parsed record (KeepRawSamples): the FORMAT + sample columns
		// were kept verbatim. Re-emitting them is byte-identical to parsing them
		// into per-sample maps and re-serialising, for a well-formed record.
		bw.WriteByte('\t')
		bw.WriteString(variant.RawTail)
	} else if len(variant.Samples) > 0 && len(variant.Format) > 0 {
		bw.WriteByte('\t')
		for i, f := range variant.Format {
			if i > 0 {
				bw.WriteByte(':')
			}
			bw.WriteString(f)
		}
		for _, sample := range variant.Samples {
			bw.WriteByte('\t')
			for i, format := range variant.Format {
				if i > 0 {
					bw.WriteByte(':')
				}
				if val, ok := sample.Data[format]; ok {
					bw.WriteString(val)
				} else {
					bw.WriteByte('.')
				}
			}
		}
	}

	return bw.WriteByte('\n')
}

// writeInfo streams the INFO column directly into the buffered writer,
// reproducing formatInfo's ordering exactly: keys in InfoOrder first (in order),
// then any remaining keys sorted ascending. It returns whether it wrote
// anything (false for an empty INFO map, where the caller emits "."). The
// extras pass reuses w.infoSeen to avoid allocating a set per record.
func (w *Writer) writeInfo(info map[string]string, order []string) bool {
	if len(info) == 0 {
		return false
	}
	bw := w.w
	wrote := false
	emit := func(k, v string) {
		if wrote {
			bw.WriteByte(';')
		}
		bw.WriteString(k)
		if v != "" {
			bw.WriteByte('=')
			bw.WriteString(v)
		}
		wrote = true
	}
	if w.infoSeen == nil {
		w.infoSeen = make(map[string]bool, len(info))
	}
	seen := w.infoSeen
	clear(seen)
	for _, k := range order {
		if v, ok := info[k]; ok {
			emit(k, v)
			seen[k] = true
		}
	}
	if len(seen) < len(info) {
		extras := make([]string, 0, len(info)-len(seen))
		for k := range info {
			if !seen[k] {
				extras = append(extras, k)
			}
		}
		sortStringsAsc(extras)
		for _, k := range extras {
			emit(k, info[k])
		}
	}
	return wrote
}

// Flush writes any buffered data to the underlying writer.
func (w *Writer) Flush() error {
	return w.w.Flush()
}

// WriteAll writes all variant records and flushes.
func (w *Writer) WriteAll(variants []*Variant) error {
	for _, variant := range variants {
		if err := w.Write(variant); err != nil {
			return err
		}
	}
	return w.Flush()
}

// formatQual formats a QUAL value byte-for-byte as upstream htslib's
// vcf_format() does. htslib stores QUAL as a 32-bit C float (bcf1_t.qual) and
// prints it with kputd, which is equivalent to C printf("%g", ...): six
// significant digits, switching to scientific notation (e+NN) for magnitudes
// outside the [0.0001, 999999] window (e.g. 4294967296 -> "4.29497e+09",
// 1000000 -> "1e+06"). The value is first narrowed to float32 so that the
// rounding at the sixth significant figure matches upstream exactly.
func formatQual(q float64) string {
	return FormatVCFFloat32(q)
}

// FormatVCFFloat32 renders a float for a VCF text field (QUAL, or a Float-typed
// INFO/FORMAT value) byte-for-byte as upstream htslib's kputd does. htslib
// stores these as 32-bit C floats and prints them with C's "%g" conversion
// (six significant digits, trailing zeros stripped, scientific notation for
// large/small magnitudes). The argument is narrowed to float32 first to match
// upstream's storage precision and rounding.
//
// Special spellings match htslib: "-0" for negative zero, "inf"/"-inf" for the
// infinities, and "nan" for NaN.
func FormatVCFFloat32(v float64) string {
	return formatVCFFloatG(float64(float32(v)))
}

// FormatVCFFloat64 renders a float for a VCF text field as htslib's kputd does
// (C "%g": six significant digits), but WITHOUT first narrowing to float32. Use
// this for values that htslib keeps in double precision; most VCF float fields
// (QUAL and Float-typed INFO/FORMAT) are float32 and should use
// FormatVCFFloat32 instead.
func FormatVCFFloat64(v float64) string {
	return formatVCFFloatG(v)
}

// formatVCFFloatG is the shared core of the VCF float formatters: it reproduces
// C printf("%g", v) plus htslib's special spellings. Go's
// strconv.FormatFloat(v, 'g', 6, 64) is byte-identical to C "%g" across the
// full double range (validated against the libc implementation), with the sole
// exception of negative zero, which Go prints as "0" but C prints as "-0".
func formatVCFFloatG(v float64) string {
	switch {
	case math.IsNaN(v):
		return "nan"
	case math.IsInf(v, 1):
		return "inf"
	case math.IsInf(v, -1):
		return "-inf"
	case v == 0 && math.Signbit(v):
		return "-0"
	}
	return strconv.FormatFloat(v, 'g', 6, 64)
}

// parseFilter parses the FILTER field.
func parseFilter(filter string) []string {
	if filter == "." || filter == "PASS" {
		return []string{filter}
	}
	return strings.Split(filter, ";")
}

// parseInfo parses the INFO field into a map.
// Deprecated: prefer parseInfoWithOrder so the original key order can be
// preserved on write. Retained for backward compatibility with callers that
// construct Variants by hand.
func parseInfo(info string) map[string]string {
	m, _ := parseInfoWithOrder(info)
	return m
}

// parseInfoWithOrder parses the INFO field, returning the key->value map and
// the keys in insertion order. Empty values denote flags.
func parseInfoWithOrder(info string) (map[string]string, []string) {
	result := make(map[string]string)
	if info == "." || info == "" {
		return result, nil
	}
	pairs := strings.Split(info, ";")
	order := make([]string, 0, len(pairs))
	for _, pair := range pairs {
		parts := strings.SplitN(pair, "=", 2)
		k := parts[0]
		if k == "" {
			continue
		}
		if _, seen := result[k]; !seen {
			order = append(order, k)
		}
		if len(parts) == 2 {
			result[k] = parts[1]
		} else {
			result[k] = ""
		}
	}
	return result, order
}

// formatInfo formats the INFO map into a string preserving the order in
// `order` first, then appending any keys present in `info` but missing from
// `order` (alphabetical for determinism). When `order` is nil/empty the
// function falls back to alphabetical order — this happens for variants
// constructed by hand without InfoOrder being set.
func formatInfo(info map[string]string, order []string) string {
	if len(info) == 0 {
		return ""
	}
	var pairs []string
	seen := make(map[string]bool, len(info))
	for _, k := range order {
		v, ok := info[k]
		if !ok {
			continue
		}
		seen[k] = true
		if v == "" {
			pairs = append(pairs, k)
		} else {
			pairs = append(pairs, fmt.Sprintf("%s=%s", k, v))
		}
	}
	if len(seen) < len(info) {
		extras := make([]string, 0, len(info)-len(seen))
		for k := range info {
			if !seen[k] {
				extras = append(extras, k)
			}
		}
		sortStringsAsc(extras)
		for _, k := range extras {
			v := info[k]
			if v == "" {
				pairs = append(pairs, k)
			} else {
				pairs = append(pairs, fmt.Sprintf("%s=%s", k, v))
			}
		}
	}
	return strings.Join(pairs, ";")
}

// sortStringsAsc does an in-place insertion sort. Used only for the
// "left-over" keys in formatInfo, which is typically empty or tiny.
func sortStringsAsc(s []string) {
	for i := 1; i < len(s); i++ {
		for j := i; j > 0 && s[j-1] > s[j]; j-- {
			s[j-1], s[j] = s[j], s[j-1]
		}
	}
}

// GetInfoString returns an INFO field value as a string.
func (v *Variant) GetInfoString(key string) (string, bool) {
	val, ok := v.Info[key]
	return val, ok
}

// InfoString reconstructs the full INFO column ("KEY=VAL;FLAG;..."),
// preserving the original key order recorded in InfoOrder. It returns
// "." for an empty INFO, matching the VCF missing-value convention.
func (v *Variant) InfoString() string {
	s := formatInfo(v.Info, v.InfoOrder)
	if s == "" {
		return "."
	}
	return s
}

// GetInfoInt returns an INFO field value as an integer.
func (v *Variant) GetInfoInt(key string) (int, error) {
	val, ok := v.Info[key]
	if !ok {
		return 0, fmt.Errorf("INFO key %s not found", key)
	}
	return strconv.Atoi(val)
}

// GetInfoFloat returns an INFO field value as a float.
func (v *Variant) GetInfoFloat(key string) (float64, error) {
	val, ok := v.Info[key]
	if !ok {
		return 0, fmt.Errorf("INFO key %s not found", key)
	}
	return strconv.ParseFloat(val, 64)
}

// GetSampleGenotype returns the genotype (GT) for a sample.
func (v *Variant) GetSampleGenotype(sampleName string) (string, bool) {
	for _, sample := range v.Samples {
		if sample.Name == sampleName {
			gt, ok := sample.Data["GT"]
			return gt, ok
		}
	}
	return "", false
}

// IsHomozygousRef checks if a sample is homozygous reference.
func (v *Variant) IsHomozygousRef(sampleName string) bool {
	gt, ok := v.GetSampleGenotype(sampleName)
	if !ok {
		return false
	}
	// Common patterns: 0/0, 0|0
	return gt == "0/0" || gt == "0|0"
}

// IsHeterozygous checks if a sample is heterozygous.
func (v *Variant) IsHeterozygous(sampleName string) bool {
	gt, ok := v.GetSampleGenotype(sampleName)
	if !ok {
		return false
	}
	// Common patterns: 0/1, 1/0, 0|1, 1|0
	return gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0"
}

// IsHomozygousAlt checks if a sample is homozygous alternate.
func (v *Variant) IsHomozygousAlt(sampleName string) bool {
	gt, ok := v.GetSampleGenotype(sampleName)
	if !ok {
		return false
	}
	// Common patterns: 1/1, 1|1 (assumes single alt allele)
	return gt == "1/1" || gt == "1|1"
}
