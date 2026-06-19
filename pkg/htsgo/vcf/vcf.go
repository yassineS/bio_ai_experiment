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
	"strconv"
	"strings"
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
}

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
func (r *Reader) Read() (*Variant, error) {
	v := &Variant{}
	if err := r.readInto(v); err != nil {
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
	return r.readInto(v)
}

// readInto advances to the next data line (skipping blanks and comments) and
// parses it into v. It is the shared body of Read and ReadInto.
func (r *Reader) readInto(v *Variant) error {
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
		line := r.scanner.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		return r.parseLine(line, v)
	}
}

// parseLine fills v from a single VCF data line, reusing v's existing maps and
// slices. The string fields reference line, which is valid until the next read.
func (r *Reader) parseLine(line string, v *Variant) error {
	r.fieldsBuf = splitInto(r.fieldsBuf, line, '\t')
	fields := r.fieldsBuf
	if len(fields) < 8 {
		return fmt.Errorf("invalid VCF record: insufficient fields: %s", line)
	}
	pos, err := strconv.Atoi(fields[1])
	if err != nil {
		return fmt.Errorf("invalid position %s: %v", fields[1], err)
	}
	v.Chrom = fields[0]
	v.Pos = pos
	v.ID = fields[2]
	v.Ref = fields[3]
	v.Alt = splitInto(v.Alt, fields[4], ',')
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
	v.InfoOrder, v.Info = parseInfoInto(v.InfoOrder, v.Info, fields[7])

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
	if len(variant.Samples) > 0 && len(variant.Format) > 0 {
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

// formatQual formats a QUAL value using the same minimal-precision rules as
// upstream htslib's vcf_format(): integer values print as integers (no
// trailing ".00"), otherwise the shortest representation of the float is used
// (via strconv with %g semantics). This matches bcftools byte-for-byte for the
// common case where the QUAL was an integer in the source file.
func formatQual(q float64) string {
	if q == float64(int64(q)) {
		return strconv.FormatInt(int64(q), 10)
	}
	return strconv.FormatFloat(q, 'g', -1, 64)
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
