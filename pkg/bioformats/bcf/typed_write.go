package bcf

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// encodeRecord turns a VCF variant into the (shared, indiv) byte pair the
// outer (l_shared, l_indiv, body) frame wraps. The wire form follows the BCF
// 2.2 spec verbatim: chrom (int32), pos (int32, 0-based), rlen (int32),
// qual (float32), n_info (uint16), n_allele (uint16), n_sample+n_fmt (uint32
// packed 24/8), then typed: ID, alleles, filters, info pairs.
func (w *Writer) encodeRecord(v *vcf.Variant) ([]byte, []byte, error) {
	chromID, ok := w.chromIndex[v.Chrom]
	if !ok {
		// Unknown contig: write -1 (htslib's "no chrom") and let the consumer
		// surface the mismatch.
		chromID = -1
	}

	rlen := int32(len(v.Ref))
	if rlen <= 0 {
		rlen = 1
	}

	// Filters: convert names to dictionary indices. Missing names are silently
	// dropped — htslib does the same.
	filterIdx := make([]int32, 0, len(v.Filter))
	for _, name := range v.Filter {
		if name == "" || name == "." {
			continue
		}
		if idx, ok := w.infoIndex[name]; ok {
			filterIdx = append(filterIdx, idx)
		}
	}

	// INFO: collect (key index, dict entry, raw value text) tuples.
	type infoPair struct {
		key   int32
		entry DictEntry
		raw   string
	}
	pairs := make([]infoPair, 0, len(v.Info))
	for name, raw := range v.Info {
		idx, ok := w.infoIndex[name]
		if !ok {
			continue
		}
		var entry DictEntry
		if int(idx) < len(w.header.InfoTags) {
			entry = w.header.InfoTags[idx]
		}
		pairs = append(pairs, infoPair{key: idx, entry: entry, raw: raw})
	}

	// Shared portion.
	var shared bytes.Buffer
	binary.Write(&shared, binary.LittleEndian, uint32(chromID))
	binary.Write(&shared, binary.LittleEndian, uint32(v.Pos-1)) // 1-based → 0-based
	binary.Write(&shared, binary.LittleEndian, uint32(rlen))
	if v.Qual < 0 {
		// VCF Variant uses -1 to mean missing; emit the BCF sentinel.
		binary.Write(&shared, binary.LittleEndian, MissingFloat32)
	} else {
		binary.Write(&shared, binary.LittleEndian, math.Float32bits(float32(v.Qual)))
	}
	binary.Write(&shared, binary.LittleEndian, uint16(len(pairs)))

	nAllele := uint16(1 + len(v.Alt))
	if v.Ref == "" {
		nAllele = uint16(len(v.Alt))
	}
	binary.Write(&shared, binary.LittleEndian, nAllele)

	nSample := uint32(len(w.header.Samples))
	nFmt := uint8(len(v.Format))
	packed := (nSample & 0x00FFFFFF) | (uint32(nFmt) << 24)
	binary.Write(&shared, binary.LittleEndian, packed)

	// id
	if v.ID == "" || v.ID == "." {
		shared.Write(EncodeMissing())
	} else {
		shared.Write(EncodeTypedString(v.ID))
	}
	// REF + ALTs
	if v.Ref != "" {
		shared.Write(EncodeTypedString(v.Ref))
	}
	for _, alt := range v.Alt {
		shared.Write(EncodeTypedString(alt))
	}
	// Filters: missing (empty list) emitted as TypeMissing; one entry emitted
	// as a single int (smallest width); multiple entries as a typed vector.
	shared.Write(encodeInts(filterIdx))

	// INFO pairs
	for _, p := range pairs {
		shared.Write(encodeInts([]int32{p.key}))
		shared.Write(encodeInfoValue(p.entry, p.raw))
	}

	// Per-sample portion.
	var indiv bytes.Buffer
	if nFmt > 0 && nSample > 0 {
		for _, fmtTag := range v.Format {
			fIdx, ok := w.fmtIndex[fmtTag]
			if !ok {
				return nil, nil, fmt.Errorf("bcf: FORMAT tag %q not in header", fmtTag)
			}
			indiv.Write(encodeInts([]int32{fIdx}))
			var entry DictEntry
			if int(fIdx) < len(w.header.FmtTags) {
				entry = w.header.FmtTags[fIdx]
			}
			payload, err := encodeFormatField(entry, fmtTag, v.Samples, int(nSample))
			if err != nil {
				return nil, nil, err
			}
			indiv.Write(payload)
		}
	}

	return shared.Bytes(), indiv.Bytes(), nil
}

// encodeInfoValue returns the typed bytes (descriptor + payload) for one INFO
// value. The entry's declared Type and the textual raw value together drive
// the choice of TypeInt8/16/32, TypeFloat, or TypeChar. Flags emit a single
// missing-equivalent: htslib treats them as a typed integer "1" so consumers
// see them as present.
func encodeInfoValue(entry DictEntry, raw string) []byte {
	t := strings.ToLower(entry.Type)
	switch t {
	case "flag":
		// Flag = present. htslib encodes flags as a single int8 of value 1.
		return EncodeTypedInt8(1)
	case "integer":
		return encodeIntsFromText(raw)
	case "float":
		return encodeFloatsFromText(raw)
	default:
		// String / Character / unknown → typed char vector.
		if raw == "" || raw == "." {
			return EncodeMissing()
		}
		return EncodeTypedString(raw)
	}
}

// encodeInts picks the narrowest integer width that fits every value (and
// keeps room for the missing/end-of-vector sentinels at the chosen width).
func encodeInts(vs []int32) []byte {
	if len(vs) == 0 {
		return EncodeMissing()
	}
	width := pickIntWidth(vs)
	switch width {
	case 1:
		payload := make([]byte, len(vs))
		for i, v := range vs {
			payload[i] = byte(int8(v))
		}
		return encodeTypedRaw(TypeInt8, len(vs), payload)
	case 2:
		payload := make([]byte, len(vs)*2)
		for i, v := range vs {
			binary.LittleEndian.PutUint16(payload[i*2:], uint16(int16(v)))
		}
		return encodeTypedRaw(TypeInt16, len(vs), payload)
	default:
		payload := make([]byte, len(vs)*4)
		for i, v := range vs {
			binary.LittleEndian.PutUint32(payload[i*4:], uint32(v))
		}
		return encodeTypedRaw(TypeInt32, len(vs), payload)
	}
}

// pickIntWidth returns the smallest byte-width (1, 2, or 4) that can hold every
// value plus the missing/end-of-vector sentinels at that width. A value equal
// to the smaller-width's missing sentinel must escalate to the next width so
// it isn't mistaken for "missing" on decode.
func pickIntWidth(vs []int32) int {
	maxWidth := 1
	for _, v := range vs {
		w := intWidth(v)
		if w > maxWidth {
			maxWidth = w
		}
	}
	return maxWidth
}

// intWidth picks the narrowest width for a single value. Missing and
// end-of-vector inputs are reported at their own native widths.
func intWidth(v int32) int {
	switch v {
	case MissingInt32:
		return 1 // 0x80 fits as int8 missing
	case EndOfVectorInt32:
		return 1 // 0x81 fits as int8 vector-end
	}
	if v >= -120 && v <= 127 {
		return 1
	}
	if v >= -32760 && v <= 32767 {
		return 2
	}
	return 4
}

// encodeIntsFromText parses a comma-separated integer text and emits the
// typed encoding. "." entries become MissingInt32; empty input is a missing
// scalar.
func encodeIntsFromText(raw string) []byte {
	if raw == "" || raw == "." {
		return EncodeMissing()
	}
	parts := strings.Split(raw, ",")
	ints := make([]int32, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "." || p == "" {
			ints = append(ints, MissingInt32)
			continue
		}
		n, err := strconv.ParseInt(p, 10, 32)
		if err != nil {
			// Fall back to a string encoding if the value isn't a number.
			return EncodeTypedString(raw)
		}
		ints = append(ints, int32(n))
	}
	return encodeInts(ints)
}

// encodeFloatsFromText parses a comma-separated float text and emits the
// typed encoding. "." entries become MissingFloat32.
func encodeFloatsFromText(raw string) []byte {
	if raw == "" || raw == "." {
		return EncodeMissing()
	}
	parts := strings.Split(raw, ",")
	floats := make([]float32, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "." || p == "" {
			floats = append(floats, math.Float32frombits(MissingFloat32))
			continue
		}
		f, err := strconv.ParseFloat(p, 32)
		if err != nil {
			return EncodeTypedString(raw)
		}
		floats = append(floats, float32(f))
	}
	return EncodeTypedFloatVec(floats)
}

// encodeFormatField returns the typed bytes for one FORMAT field across all
// samples. Vectors are zero-padded with the end-of-vector sentinel so every
// sample has the same on-wire dimension, exactly as htslib does.
func encodeFormatField(entry DictEntry, key string, samples []vcf.Sample, nSample int) ([]byte, error) {
	t := strings.ToLower(entry.Type)
	if key == "GT" || (entry.ID == "GT" && t == "string") {
		return encodeFormatGT(samples, nSample), nil
	}
	switch t {
	case "integer":
		return encodeFormatInts(samples, nSample, key), nil
	case "float":
		return encodeFormatFloats(samples, nSample, key), nil
	default:
		return encodeFormatChars(samples, nSample, key), nil
	}
}

// encodeFormatGT writes a FORMAT/GT field. Allele indices are stored as
// (allele+1)<<1 | phased; missing is 0; the second-and-later alleles per
// sample carry the phased bit on themselves.
func encodeFormatGT(samples []vcf.Sample, nSample int) []byte {
	// Determine the max ploidy across samples.
	maxPloidy := 1
	parsed := make([][]int32, nSample)
	for i := 0; i < nSample; i++ {
		var gt string
		if i < len(samples) {
			gt = samples[i].Data["GT"]
		}
		if gt == "" {
			parsed[i] = []int32{MissingInt32}
			continue
		}
		parsed[i] = parseGT(gt)
		if len(parsed[i]) > maxPloidy {
			maxPloidy = len(parsed[i])
		}
	}
	// Build the flat int32 vector with EndOfVector padding.
	flat := make([]int32, nSample*maxPloidy)
	for i := 0; i < nSample; i++ {
		row := parsed[i]
		for k := 0; k < maxPloidy; k++ {
			if k < len(row) {
				flat[i*maxPloidy+k] = row[k]
			} else {
				flat[i*maxPloidy+k] = EndOfVectorInt32
			}
		}
	}
	return encodeIntsVec(flat)
}

// parseGT turns "0/1" / "1|0" / "./." into the on-wire int32 encoding.
// The first element does not carry a phased bit (htslib convention).
func parseGT(gt string) []int32 {
	if gt == "" || gt == "." {
		return []int32{MissingInt32}
	}
	// Replace pipes with slashes so we can split — and remember per-position
	// whether the original separator was a pipe.
	out := make([]int32, 0, 4)
	start := 0
	phased := false
	for i := 0; i <= len(gt); i++ {
		if i == len(gt) || gt[i] == '/' || gt[i] == '|' {
			tok := gt[start:i]
			var v int32
			if tok == "." || tok == "" {
				// missing genotype slot
				v = MissingInt32
			} else {
				n, err := strconv.Atoi(tok)
				if err != nil || n < 0 {
					v = MissingInt32
				} else {
					// allele+1 << 1, low bit = phased
					v = int32((n + 1) << 1)
					if phased && len(out) > 0 {
						v |= 1
					}
				}
			}
			out = append(out, v)
			if i < len(gt) {
				phased = gt[i] == '|'
				start = i + 1
			}
		}
	}
	return out
}

// encodeIntsVec is the vectorized counterpart of encodeInts — same width
// selection logic, but exposed for callers that already have flat int32 data.
func encodeIntsVec(vs []int32) []byte { return encodeInts(vs) }

// encodeFormatInts writes one FORMAT integer field across nSample samples.
// Per-sample multi-value fields are flat-padded with EndOfVectorInt32 to a
// uniform dimension.
func encodeFormatInts(samples []vcf.Sample, nSample int, key string) []byte {
	rows := make([][]int32, nSample)
	maxDim := 1
	for i := 0; i < nSample; i++ {
		var raw string
		if i < len(samples) {
			raw = samples[i].Data[key]
		}
		if raw == "" || raw == "." {
			rows[i] = []int32{MissingInt32}
			continue
		}
		parts := strings.Split(raw, ",")
		row := make([]int32, len(parts))
		for j, p := range parts {
			p = strings.TrimSpace(p)
			if p == "." || p == "" {
				row[j] = MissingInt32
				continue
			}
			n, err := strconv.ParseInt(p, 10, 32)
			if err != nil {
				row[j] = MissingInt32
				continue
			}
			row[j] = int32(n)
		}
		rows[i] = row
		if len(row) > maxDim {
			maxDim = len(row)
		}
	}
	flat := make([]int32, nSample*maxDim)
	for i := 0; i < nSample; i++ {
		for k := 0; k < maxDim; k++ {
			if k < len(rows[i]) {
				flat[i*maxDim+k] = rows[i][k]
			} else {
				flat[i*maxDim+k] = EndOfVectorInt32
			}
		}
	}
	return encodeInts(flat)
}

// encodeFormatFloats is encodeFormatInts for floats.
func encodeFormatFloats(samples []vcf.Sample, nSample int, key string) []byte {
	rows := make([][]float32, nSample)
	maxDim := 1
	for i := 0; i < nSample; i++ {
		var raw string
		if i < len(samples) {
			raw = samples[i].Data[key]
		}
		if raw == "" || raw == "." {
			rows[i] = []float32{math.Float32frombits(MissingFloat32)}
			continue
		}
		parts := strings.Split(raw, ",")
		row := make([]float32, len(parts))
		for j, p := range parts {
			p = strings.TrimSpace(p)
			if p == "." || p == "" {
				row[j] = math.Float32frombits(MissingFloat32)
				continue
			}
			f, err := strconv.ParseFloat(p, 32)
			if err != nil {
				row[j] = math.Float32frombits(MissingFloat32)
				continue
			}
			row[j] = float32(f)
		}
		rows[i] = row
		if len(row) > maxDim {
			maxDim = len(row)
		}
	}
	flat := make([]float32, nSample*maxDim)
	endBits := math.Float32frombits(EndOfVectorFloat)
	for i := 0; i < nSample; i++ {
		for k := 0; k < maxDim; k++ {
			if k < len(rows[i]) {
				flat[i*maxDim+k] = rows[i][k]
			} else {
				flat[i*maxDim+k] = endBits
			}
		}
	}
	return EncodeTypedFloatVec(flat)
}

// encodeFormatChars writes FORMAT char fields. The on-wire form is a flat
// char vector padded to a per-sample stride so consumers can carve it up
// by nSample.
func encodeFormatChars(samples []vcf.Sample, nSample int, key string) []byte {
	values := make([]string, nSample)
	maxLen := 1
	for i := 0; i < nSample; i++ {
		if i < len(samples) {
			values[i] = samples[i].Data[key]
		}
		if l := len(values[i]); l > maxLen {
			maxLen = l
		}
	}
	if maxLen == 0 {
		return EncodeMissing()
	}
	buf := make([]byte, nSample*maxLen)
	for i := 0; i < nSample; i++ {
		v := values[i]
		copy(buf[i*maxLen:], v)
		// Remaining bytes (after copy) are already zero — htslib uses NUL as
		// the trailing pad for char vectors.
	}
	return encodeTypedRaw(TypeChar, len(buf), buf)
}

// encodeRecordRaw re-emits an already-decoded Record without going through
// the vcf.Variant translation. It is used by Writer.WriteRecord and by tests
// that want to verify lossless round-trips.
func encodeRecordRaw(r *Record) ([]byte, []byte, error) {
	var shared bytes.Buffer
	binary.Write(&shared, binary.LittleEndian, uint32(r.ChromID))
	binary.Write(&shared, binary.LittleEndian, uint32(r.Pos))
	binary.Write(&shared, binary.LittleEndian, uint32(r.Rlen))
	binary.Write(&shared, binary.LittleEndian, math.Float32bits(r.Qual))
	binary.Write(&shared, binary.LittleEndian, uint16(len(r.InfoKeys)))
	binary.Write(&shared, binary.LittleEndian, uint16(len(r.Alleles)))
	packed := (r.NSample & 0x00FFFFFF) | (uint32(r.NFmt) << 24)
	binary.Write(&shared, binary.LittleEndian, packed)

	if r.ID == "" {
		shared.Write(EncodeMissing())
	} else {
		shared.Write(EncodeTypedString(r.ID))
	}
	for _, al := range r.Alleles {
		shared.Write(EncodeTypedString(al))
	}
	shared.Write(encodeInts(r.Filters))
	for i, k := range r.InfoKeys {
		shared.Write(encodeInts([]int32{k}))
		shared.Write(encodeTypedValue(r.InfoVals[i]))
	}

	var indiv bytes.Buffer
	for i, k := range r.FmtKeys {
		indiv.Write(encodeInts([]int32{k}))
		indiv.Write(encodeTypedValue(r.FmtVals[i]))
	}
	return shared.Bytes(), indiv.Bytes(), nil
}

// encodeTypedValue serialises a TypedValue back to its on-wire bytes. We
// keep the descriptor verbatim so a Reader→Writer round-trip is byte-stable.
func encodeTypedValue(tv TypedValue) []byte {
	switch tv.Descriptor {
	case TypeMissing:
		return EncodeMissing()
	case TypeInt8:
		payload := make([]byte, tv.Length)
		for i, v := range tv.Ints {
			if i >= tv.Length {
				break
			}
			var b int8
			switch v {
			case MissingInt32:
				b = MissingInt8
			case EndOfVectorInt32:
				b = EndOfVectorInt8
			default:
				b = int8(v)
			}
			payload[i] = byte(b)
		}
		return encodeTypedRaw(TypeInt8, tv.Length, payload)
	case TypeInt16:
		payload := make([]byte, tv.Length*2)
		for i, v := range tv.Ints {
			if i >= tv.Length {
				break
			}
			var s int16
			switch v {
			case MissingInt32:
				s = MissingInt16
			case EndOfVectorInt32:
				s = EndOfVectorInt16
			default:
				s = int16(v)
			}
			binary.LittleEndian.PutUint16(payload[i*2:], uint16(s))
		}
		return encodeTypedRaw(TypeInt16, tv.Length, payload)
	case TypeInt32:
		payload := make([]byte, tv.Length*4)
		for i, v := range tv.Ints {
			if i >= tv.Length {
				break
			}
			binary.LittleEndian.PutUint32(payload[i*4:], uint32(v))
		}
		return encodeTypedRaw(TypeInt32, tv.Length, payload)
	case TypeFloat:
		payload := make([]byte, tv.Length*4)
		for i, f := range tv.Floats {
			if i >= tv.Length {
				break
			}
			binary.LittleEndian.PutUint32(payload[i*4:], math.Float32bits(f))
		}
		return encodeTypedRaw(TypeFloat, tv.Length, payload)
	case TypeChar:
		return encodeTypedRaw(TypeChar, tv.Length, []byte(tv.String))
	}
	return EncodeMissing()
}
