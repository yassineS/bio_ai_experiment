package bcf

import (
	"encoding/binary"
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Record is the decoded form of one BCF record. Numeric fields are stored in
// their native types; everything else (INFO, FORMAT) is decoded into
// strings so the result can be emitted as VCF text without further work.
type Record struct {
	ChromID  int32    // dictionary index; -1 for "no contig"
	Pos      int32    // 0-based on the wire
	Rlen     int32    // length of REF
	Qual     float32  // missing carries the IEEE-754 bit pattern 0x7F800001
	NInfo    uint16   // count of INFO fields
	NAllele  uint16   // count of alleles including REF
	NSample  uint32   // count of samples (only 24 bits used)
	NFmt     uint8    // count of FORMAT keys
	ID       string   // VCF ID column ("." if missing)
	Alleles  []string // REF (index 0) followed by ALTs
	Filters  []int32  // dictionary indices into Header.InfoTags
	InfoKeys []int32  // dictionary indices into Header.InfoTags, length NInfo
	InfoVals []TypedValue
	FmtKeys  []int32      // dictionary indices into Header.FmtTags, length NFmt
	FmtVals  []TypedValue // each has Length == NSample (or NSample * dimension)
}

// ID returns Record.ID, or "." if it is empty (the canonical VCF missing).
func (r *Record) IDOrMissing() string {
	if r.ID == "" {
		return "."
	}
	return r.ID
}

// QualOrMissing returns Qual unless it is the BCF "missing" pattern, in
// which case it returns NaN so callers can detect it with math.IsNaN.
func (r *Record) QualOrMissing() float32 {
	if math.Float32bits(r.Qual) == MissingFloat32 {
		return float32(math.NaN())
	}
	return r.Qual
}

// QualString returns the VCF text representation of Qual: "." for missing,
// the float formatted with up to 6 significant digits otherwise.
func (r *Record) QualString() string {
	if math.Float32bits(r.Qual) == MissingFloat32 {
		return "."
	}
	if r.Qual == float32(int32(r.Qual)) {
		// Print whole-number qualities without a decimal point to match
		// how bcftools renders them.
		return strconv.FormatInt(int64(r.Qual), 10)
	}
	return strconv.FormatFloat(float64(r.Qual), 'g', -1, 32)
}

// decodeShared parses the shared portion of a record (size = lShared) from
// buf into r. The header argument is accepted for future use (e.g. validating
// dictionary indices), but the decoder is dictionary-agnostic in this slice.
func decodeShared(buf []byte, _ *Header, r *Record) error {
	off := 0
	if len(buf) < 24 {
		return fmt.Errorf("bcf: shared portion truncated (%d bytes)", len(buf))
	}
	r.ChromID = int32(binary.LittleEndian.Uint32(buf[off : off+4]))
	off += 4
	r.Pos = int32(binary.LittleEndian.Uint32(buf[off : off+4]))
	off += 4
	r.Rlen = int32(binary.LittleEndian.Uint32(buf[off : off+4]))
	off += 4
	qualBits := binary.LittleEndian.Uint32(buf[off : off+4])
	r.Qual = math.Float32frombits(qualBits)
	off += 4
	r.NInfo = binary.LittleEndian.Uint16(buf[off : off+2])
	off += 2
	r.NAllele = binary.LittleEndian.Uint16(buf[off : off+2])
	off += 2
	packed := binary.LittleEndian.Uint32(buf[off : off+4])
	off += 4
	r.NSample = packed & 0x00FFFFFF
	r.NFmt = uint8(packed >> 24)

	// id
	id, err := DecodeTypedString(buf, &off)
	if err != nil {
		return fmt.Errorf("bcf: id: %w", err)
	}
	r.ID = id

	// REF + ALTs: n_allele typed strings in a row
	r.Alleles = make([]string, 0, r.NAllele)
	for i := uint16(0); i < r.NAllele; i++ {
		s, err := DecodeTypedString(buf, &off)
		if err != nil {
			return fmt.Errorf("bcf: allele %d: %w", i, err)
		}
		r.Alleles = append(r.Alleles, s)
	}

	// FILTER: typed int vector of dictionary indices
	filters, err := DecodeTypedInts(buf, &off)
	if err != nil {
		return fmt.Errorf("bcf: filter: %w", err)
	}
	r.Filters = filters

	// INFO: n_info pairs of (key, typed value)
	r.InfoKeys = make([]int32, 0, r.NInfo)
	r.InfoVals = make([]TypedValue, 0, r.NInfo)
	for i := uint16(0); i < r.NInfo; i++ {
		key, err := DecodeTypedInt(buf, &off)
		if err != nil {
			return fmt.Errorf("bcf: info key %d: %w", i, err)
		}
		val, err := DecodeTyped(buf, &off)
		if err != nil {
			return fmt.Errorf("bcf: info value %d: %w", i, err)
		}
		r.InfoKeys = append(r.InfoKeys, key)
		r.InfoVals = append(r.InfoVals, val)
	}
	if off != len(buf) {
		// Some encoders emit trailing padding; tolerate it but record the
		// number of trailing bytes for debug.
		// We intentionally do not return an error here — the BCF spec
		// allows additional shared-portion fields in future versions.
	}
	return nil
}

// decodeIndiv parses the per-sample portion of a record (size = lIndiv).
// It populates r.FmtKeys and r.FmtVals. The descriptor's size field is the
// PER-SAMPLE dimension; the payload spans nSample × size elements per
// FORMAT field. We use DecodeFormatTyped with r.NSample to read the right
// number of elements while keeping TypedValue.Length holding the
// per-sample dim for downstream stride-aware consumers.
func decodeIndiv(buf []byte, _ *Header, r *Record) error {
	if r.NFmt == 0 {
		return nil
	}
	off := 0
	r.FmtKeys = make([]int32, 0, r.NFmt)
	r.FmtVals = make([]TypedValue, 0, r.NFmt)
	nSample := int(r.NSample)
	if nSample == 0 {
		// Defensive: with no samples, treat each field as a singleton.
		nSample = 1
	}
	for i := uint8(0); i < r.NFmt; i++ {
		key, err := DecodeTypedInt(buf, &off)
		if err != nil {
			return fmt.Errorf("bcf: fmt key %d: %w", i, err)
		}
		val, err := DecodeFormatTyped(buf, &off, nSample)
		if err != nil {
			return fmt.Errorf("bcf: fmt value %d: %w", i, err)
		}
		r.FmtKeys = append(r.FmtKeys, key)
		r.FmtVals = append(r.FmtVals, val)
	}
	return nil
}

// ToVariant builds a vcf.Variant from r, resolving dictionary indices via
// hdr. The conversion produces the VCF-text view of the record (1-based POS,
// "." for missing strings, FORMAT/sample reconstruction).
func (r *Record) ToVariant(hdr *Header) *vcf.Variant {
	v := &vcf.Variant{
		Chrom: hdr.ContigName(r.ChromID),
		Pos:   int(r.Pos) + 1, // wire is 0-based, VCF text is 1-based
		ID:    r.IDOrMissing(),
		Info:  map[string]string{},
	}
	if len(r.Alleles) > 0 {
		v.Ref = r.Alleles[0]
		v.Alt = append([]string{}, r.Alleles[1:]...)
	}
	if math.Float32bits(r.Qual) == MissingFloat32 {
		v.Qual = -1
	} else {
		v.Qual = float64(r.Qual)
	}

	// Filters: empty → ".", a single entry with ID==PASS → "PASS".
	if len(r.Filters) == 0 {
		v.Filter = []string{"."}
	} else {
		for _, fi := range r.Filters {
			entry := hdr.InfoTag(fi)
			if entry == nil {
				v.Filter = append(v.Filter, ".")
				continue
			}
			v.Filter = append(v.Filter, entry.ID)
		}
	}

	// INFO — preserve key order from the BCF record for byte-for-byte
	// text-output parity with upstream bcftools.
	v.InfoOrder = make([]string, 0, len(r.InfoKeys))
	for i, key := range r.InfoKeys {
		entry := hdr.InfoTag(key)
		if entry == nil {
			continue
		}
		if _, seen := v.Info[entry.ID]; !seen {
			v.InfoOrder = append(v.InfoOrder, entry.ID)
		}
		v.Info[entry.ID] = formatTyped(r.InfoVals[i], entry)
	}

	// FORMAT / samples
	if len(r.FmtKeys) > 0 && len(hdr.Samples) > 0 {
		nSample := len(hdr.Samples)
		formatTags := make([]string, 0, len(r.FmtKeys))
		samples := make([]vcf.Sample, nSample)
		for i := range samples {
			samples[i] = vcf.Sample{Name: hdr.Samples[i], Data: map[string]string{}}
		}
		for i, key := range r.FmtKeys {
			entry := hdr.FmtTag(key)
			if entry == nil {
				continue
			}
			formatTags = append(formatTags, entry.ID)
			perSample := splitPerSample(r.FmtVals[i], nSample, entry, key, hdr)
			for s := 0; s < nSample; s++ {
				val := "."
				if s < len(perSample) {
					val = perSample[s]
				}
				samples[s].Data[entry.ID] = val
			}
		}
		v.Format = formatTags
		v.Samples = samples
	}

	return v
}

// formatTyped renders a single INFO TypedValue as the VCF text form. The
// dictionary entry is used to detect Type=Flag (no "=value" in VCF) and to
// pick the right numeric formatting.
func formatTyped(tv TypedValue, entry *DictEntry) string {
	if entry != nil && strings.EqualFold(entry.Type, "Flag") {
		// Flags carry no payload; the dictionary lookup told us this is
		// a flag, so VCF text is just the bare tag (handled by the caller
		// via the empty-string value convention).
		return ""
	}
	switch tv.Descriptor {
	case TypeInt8, TypeInt16, TypeInt32:
		parts := make([]string, 0, len(tv.Ints))
		for _, v := range tv.Ints {
			if v == MissingInt32 || v == EndOfVectorInt32 {
				parts = append(parts, ".")
				continue
			}
			parts = append(parts, strconv.FormatInt(int64(v), 10))
		}
		return strings.Join(parts, ",")
	case TypeFloat:
		parts := make([]string, 0, len(tv.Floats))
		for _, f := range tv.Floats {
			bits := math.Float32bits(f)
			if IsMissingFloat(bits) || IsEndOfVectorFloat(bits) {
				parts = append(parts, ".")
				continue
			}
			parts = append(parts, formatFloat(f))
		}
		return strings.Join(parts, ",")
	case TypeChar:
		return tv.String
	case TypeMissing:
		return "."
	}
	return ""
}

// splitPerSample turns a per-sample TypedValue into nSample VCF strings.
// FORMAT fields can be scalar (one value per sample) or vector (e.g. GT, AD)
// in which case `tv.Length` carries the per-sample dimension and the flat
// payload (tv.Ints / tv.Floats / tv.String) spans nSample × tv.Length
// elements. We re-pack with the correct separators.
func splitPerSample(tv TypedValue, nSample int, entry *DictEntry, _ int32, _ *Header) []string {
	if nSample == 0 {
		return nil
	}
	dim := tv.Length
	if dim < 1 {
		dim = 1
	}
	out := make([]string, nSample)
	isGT := entry != nil && entry.ID == "GT"

	switch tv.Descriptor {
	case TypeInt8, TypeInt16, TypeInt32:
		for s := 0; s < nSample; s++ {
			vals := make([]string, 0, dim)
			endOfVec := false
			for k := 0; k < dim; k++ {
				idx := s*dim + k
				if idx >= len(tv.Ints) {
					break
				}
				v := tv.Ints[idx]
				if v == EndOfVectorInt32 {
					endOfVec = true
					break
				}
				if isGT {
					vals = append(vals, formatGTAllele(v))
					continue
				}
				if v == MissingInt32 {
					vals = append(vals, ".")
					continue
				}
				vals = append(vals, strconv.FormatInt(int64(v), 10))
			}
			if isGT {
				out[s] = joinGTAlleles(vals, tv, s, dim)
				_ = endOfVec
				continue
			}
			if len(vals) == 0 {
				out[s] = "."
			} else {
				out[s] = strings.Join(vals, ",")
			}
		}
	case TypeFloat:
		for s := 0; s < nSample; s++ {
			vals := make([]string, 0, dim)
			for k := 0; k < dim; k++ {
				idx := s*dim + k
				if idx >= len(tv.Floats) {
					break
				}
				bits := math.Float32bits(tv.Floats[idx])
				if IsEndOfVectorFloat(bits) {
					break
				}
				if IsMissingFloat(bits) {
					vals = append(vals, ".")
					continue
				}
				vals = append(vals, formatFloat(tv.Floats[idx]))
			}
			if len(vals) == 0 {
				out[s] = "."
			} else {
				out[s] = strings.Join(vals, ",")
			}
		}
	case TypeChar:
		// FORMAT char vectors are encoded as a flat string with one
		// fixed-width per-sample slot on the wire. The slot width is the
		// descriptor's per-sample dimension (tv.Length), which htslib pads
		// with trailing NULs so every sample occupies the same number of
		// bytes. Strides derived from len(s)/nSample coincide with tv.Length
		// for well-formed input, but we use tv.Length directly so a final
		// sample whose payload happens to be all-NUL (and might be trimmed by
		// a non-conforming encoder) still lands in the right slot.
		s := tv.String
		stride := dim
		if stride*nSample > len(s) {
			// Fall back to an even split when the payload is shorter than the
			// declared width would imply (defensive; should not happen with
			// htslib output).
			stride = len(s) / nSample
			if stride == 0 {
				stride = len(s)
			}
		}
		for i := 0; i < nSample && i*stride < len(s); i++ {
			end := (i + 1) * stride
			if end > len(s) {
				end = len(s)
			}
			out[i] = strings.TrimRight(s[i*stride:end], "\x00")
			if out[i] == "" {
				out[i] = "."
			}
		}
	default:
		for s := 0; s < nSample; s++ {
			out[s] = "."
		}
	}
	return out
}

// formatGTAllele renders a single GT typed integer. The on-wire form is the
// allele index left-shifted by 1, with the low bit indicating "phased". A
// value of 0 (i.e. (-1 << 1)|0 == 0) is the conventional ".", and the
// missing sentinel maps to "." as well.
func formatGTAllele(v int32) string {
	if v == EndOfVectorInt32 {
		return ""
	}
	if v == MissingInt32 || v == 0 {
		return "."
	}
	allele := (v >> 1) - 1
	if allele < 0 {
		return "."
	}
	return strconv.FormatInt(int64(allele), 10)
}

// joinGTAlleles reconstructs the GT field. We need to look back at the
// raw int values to learn the "phased" bit on each non-first entry: phased
// alleles get '|' as the separator, unphased get '/'.
func joinGTAlleles(alleleStrs []string, tv TypedValue, s, dim int) string {
	if len(alleleStrs) == 0 {
		return "."
	}
	var b strings.Builder
	b.WriteString(alleleStrs[0])
	for k := 1; k < len(alleleStrs); k++ {
		idx := s*dim + k
		sep := "/"
		if idx < len(tv.Ints) && tv.Ints[idx]&1 == 1 {
			sep = "|"
		}
		b.WriteString(sep)
		b.WriteString(alleleStrs[k])
	}
	return b.String()
}

// formatFloat renders a float in the compact form bcftools uses.
func formatFloat(f float32) string {
	if f == float32(int32(f)) {
		return strconv.FormatInt(int64(f), 10)
	}
	return strconv.FormatFloat(float64(f), 'g', -1, 32)
}
