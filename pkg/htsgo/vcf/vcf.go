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
	if r.err != nil {
		return nil, r.err
	}

	if r.header == nil || len(r.header.MetaInfo) == 0 {
		return nil, fmt.Errorf("header not read; call ReadHeader() first")
	}

	if !r.scanner.Scan() {
		if err := r.scanner.Err(); err != nil {
			r.err = err
			return nil, err
		}
		r.err = io.EOF
		return nil, io.EOF
	}

	line := r.scanner.Text()
	// Skip empty lines and comments
	if line == "" || strings.HasPrefix(line, "#") {
		return r.Read()
	}

	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return nil, fmt.Errorf("invalid VCF record: insufficient fields: %s", line)
	}

	// Parse required fields
	info, infoOrder := parseInfoWithOrder(fields[7])
	variant := &Variant{
		Chrom:     fields[0],
		ID:        fields[2],
		Ref:       fields[3],
		Filter:    parseFilter(fields[6]),
		Info:      info,
		InfoOrder: infoOrder,
	}

	// Parse position
	pos, err := strconv.Atoi(fields[1])
	if err != nil {
		return nil, fmt.Errorf("invalid position %s: %v", fields[1], err)
	}
	variant.Pos = pos

	// Parse alternate alleles
	variant.Alt = strings.Split(fields[4], ",")

	// Parse quality
	if fields[5] == "." {
		variant.Qual = -1
	} else {
		qual, err := strconv.ParseFloat(fields[5], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid quality %s: %v", fields[5], err)
		}
		variant.Qual = qual
	}

	// Parse FORMAT and sample data if present
	if len(fields) > 8 {
		variant.Format = strings.Split(fields[8], ":")
		variant.Samples = make([]Sample, len(r.header.Samples))

		for i, sampleData := range fields[9:] {
			if i >= len(r.header.Samples) {
				break
			}
			sample := Sample{
				Name: r.header.Samples[i],
				Data: make(map[string]string),
			}

			values := strings.Split(sampleData, ":")
			for j, format := range variant.Format {
				if j < len(values) {
					sample.Data[format] = values[j]
				}
			}
			variant.Samples[i] = sample
		}
	}

	return variant, nil
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
	w      *bufio.Writer
	header *Header
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

// Write writes a variant record.
func (w *Writer) Write(variant *Variant) error {
	// Build the line
	var fields []string

	// Required fields
	fields = append(fields,
		variant.Chrom,
		strconv.Itoa(variant.Pos),
		variant.ID,
		variant.Ref,
		strings.Join(variant.Alt, ","),
	)

	// Quality
	if variant.Qual < 0 {
		fields = append(fields, ".")
	} else {
		fields = append(fields, formatQual(variant.Qual))
	}

	// Filter
	if len(variant.Filter) == 0 {
		fields = append(fields, ".")
	} else {
		fields = append(fields, strings.Join(variant.Filter, ";"))
	}

	// Info
	infoStr := formatInfo(variant.Info, variant.InfoOrder)
	if infoStr == "" {
		infoStr = "."
	}
	fields = append(fields, infoStr)

	// Format and samples
	if len(variant.Samples) > 0 && len(variant.Format) > 0 {
		fields = append(fields, strings.Join(variant.Format, ":"))
		for _, sample := range variant.Samples {
			values := make([]string, len(variant.Format))
			for i, format := range variant.Format {
				if val, ok := sample.Data[format]; ok {
					values[i] = val
				} else {
					values[i] = "."
				}
			}
			fields = append(fields, strings.Join(values, ":"))
		}
	}

	// Write the line
	if _, err := fmt.Fprintln(w.w, strings.Join(fields, "\t")); err != nil {
		return err
	}

	return nil
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
