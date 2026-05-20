package cram

import (
	"fmt"
)

// dataSeriesKey is a two-character CRAM data-series identifier (for
// example "BF" for the bit flags, "RL" for read lengths). The CRAM
// compression header maps each key to the Encoding that governs that
// series' on-disk representation.
type dataSeriesKey [2]byte

// String returns the two-character key as a string.
func (k dataSeriesKey) String() string { return string(k[:]) }

// tagKey is a three-byte CRAM tag identifier: the two-character SAM tag
// name followed by its one-character value type (for example the bytes
// 'M','D','Z' for an "MD" tag stored as a NUL-terminated string).
type tagKey [3]byte

// String returns the three-byte tag key as a string.
func (k tagKey) String() string { return string(k[:]) }

// PreservationMap holds the parsed CRAM preservation map: the boolean
// and raw-bytes settings that describe how the container's records were
// encoded. Only the entries the writer emitted are populated; the
// boolean fields carry their CRAM default when an entry is absent.
type PreservationMap struct {
	// ReadNamesIncluded is the RN entry: whether read names are stored
	// (true) or must be auto-generated on decode. CRAM default: true.
	ReadNamesIncluded bool
	// APDelta is the AP entry: whether the alignment-position data
	// series is stored as a delta from the previous record (true) or as
	// an absolute coordinate (false). CRAM default: true.
	APDelta bool
	// ReferenceRequired is the RR entry: whether decoding needs the
	// external reference sequence. CRAM default: true.
	ReferenceRequired bool
	// SubstitutionMatrix is the SM entry: the five raw bytes of the
	// base-substitution code matrix. It is nil when no SM entry was
	// written.
	SubstitutionMatrix []byte
	// TagDictionary is the TD entry: the raw tag-combination dictionary
	// bytes. It is nil when no TD entry was written. The dictionary is a
	// run of NUL-separated tag-key lists; C4b interprets it.
	TagDictionary []byte
	// hasRN, hasAP and hasRR record whether the corresponding entry was
	// present, so a caller can distinguish "absent" from "false".
	hasRN, hasAP, hasRR bool
}

// CompressionHeader is the fully-parsed data of a container's
// compression-header block. It carries the preservation map, the
// data-series encoding map (which encoding governs each two-character
// data series) and the tag-encoding map (which encoding governs each
// three-byte tag).
type CompressionHeader struct {
	// Preservation holds the parsed preservation map.
	Preservation PreservationMap
	// DataSeries maps each two-character data-series key to its
	// Encoding. A series absent from this map has no encoding declared.
	DataSeries map[dataSeriesKey]*Encoding
	// Tags maps each three-byte tag key to its Encoding.
	Tags map[tagKey]*Encoding
}

// Encoding returns the encoding declared for the named two-character
// data series, or nil if the series is absent from the data-series
// encoding map.
func (h *CompressionHeader) Encoding(key string) *Encoding {
	if len(key) != 2 {
		return nil
	}
	return h.DataSeries[dataSeriesKey{key[0], key[1]}]
}

// parseCompressionHeader parses the data of a CRAM compression-header
// block. The block payload, once decompressed, holds three consecutive
// maps: the preservation map, the data-series encoding map and the
// tag-encoding map. Each map is framed as an ITF-8 byte size followed
// by an ITF-8 entry count.
func parseCompressionHeader(p []byte) (*CompressionHeader, error) {
	h := &CompressionHeader{
		DataSeries: make(map[dataSeriesKey]*Encoding),
		Tags:       make(map[tagKey]*Encoding),
	}
	off := 0
	var err error
	if off, err = parsePreservationMap(p, off, &h.Preservation); err != nil {
		return nil, err
	}
	if off, err = parseDataSeriesMap(p, off, h.DataSeries); err != nil {
		return nil, err
	}
	if _, err = parseTagEncodingMap(p, off, h.Tags); err != nil {
		return nil, err
	}
	return h, nil
}

// mapFrame reads the ITF-8 byte-size and ITF-8 entry-count that prefix
// every CRAM compression-header map. It returns the offset of the first
// entry, the offset just past the map (size bytes past the entry
// count), and the entry count. The size field measures the bytes from
// the entry count onward.
func mapFrame(p []byte, off int, what string) (entryOff, mapEnd int, count int32, err error) {
	size, n, err := itf8At(p, off)
	if err != nil {
		return 0, 0, 0, fmt.Errorf("cram: %s map size: %w", what, err)
	}
	off += n
	if size < 0 {
		return 0, 0, 0, fmt.Errorf("cram: %s map declares negative size %d", what, size)
	}
	mapEnd = off + int(size)
	if mapEnd > len(p) || mapEnd < off {
		return 0, 0, 0, fmt.Errorf("cram: %s map (%d bytes) overruns the compression header", what, size)
	}
	count, n, err = itf8At(p, off)
	if err != nil {
		return 0, 0, 0, fmt.Errorf("cram: %s map entry count: %w", what, err)
	}
	off += n
	if count < 0 {
		return 0, 0, 0, fmt.Errorf("cram: %s map declares negative entry count %d", what, count)
	}
	return off, mapEnd, count, nil
}

// parsePreservationMap parses the CRAM preservation map starting at off.
// Each entry is a two-byte key followed by a value whose shape depends
// on the key: RN, AP and RR carry a single boolean byte; SM carries
// five raw bytes; TD carries an ITF-8-prefixed byte run. Unknown keys
// cannot be skipped (their value length is not self-describing), so an
// unknown key is reported as an error.
func parsePreservationMap(p []byte, off int, pm *PreservationMap) (int, error) {
	// CRAM defaults apply when an entry is absent.
	pm.ReadNamesIncluded = true
	pm.APDelta = true
	pm.ReferenceRequired = true

	entryOff, mapEnd, count, err := mapFrame(p, off, "preservation")
	if err != nil {
		return off, err
	}
	cur := entryOff
	for i := int32(0); i < count; i++ {
		if cur+2 > mapEnd {
			return off, fmt.Errorf("cram: preservation map entry %d: truncated key", i)
		}
		key := [2]byte{p[cur], p[cur+1]}
		cur += 2
		switch string(key[:]) {
		case "RN", "AP", "RR":
			if cur >= mapEnd {
				return off, fmt.Errorf("cram: preservation map %s entry: truncated boolean", key)
			}
			b := p[cur] != 0
			cur++
			switch string(key[:]) {
			case "RN":
				pm.ReadNamesIncluded, pm.hasRN = b, true
			case "AP":
				pm.APDelta, pm.hasAP = b, true
			case "RR":
				pm.ReferenceRequired, pm.hasRR = b, true
			}
		case "SM":
			// The substitution matrix is exactly five bytes.
			if cur+5 > mapEnd {
				return off, fmt.Errorf("cram: preservation map SM entry: truncated substitution matrix")
			}
			pm.SubstitutionMatrix = append([]byte(nil), p[cur:cur+5]...)
			cur += 5
		case "TD":
			// The tag dictionary is an ITF-8 byte count then that many
			// raw bytes.
			n, m, terr := itf8At(p, cur)
			if terr != nil {
				return off, fmt.Errorf("cram: preservation map TD entry size: %w", terr)
			}
			cur += m
			if n < 0 || cur+int(n) > mapEnd {
				return off, fmt.Errorf("cram: preservation map TD entry (%d bytes) overruns the map", n)
			}
			pm.TagDictionary = append([]byte(nil), p[cur:cur+int(n)]...)
			cur += int(n)
		default:
			return off, fmt.Errorf("cram: preservation map has unknown key %q (cannot skip a non-self-describing value)", key)
		}
	}
	if cur > mapEnd {
		return off, fmt.Errorf("cram: preservation map entries overran their declared size")
	}
	return mapEnd, nil
}

// parseDataSeriesMap parses the CRAM data-series encoding map starting
// at off. Each entry is a two-byte data-series key followed by one
// parsed Encoding.
func parseDataSeriesMap(p []byte, off int, out map[dataSeriesKey]*Encoding) (int, error) {
	entryOff, mapEnd, count, err := mapFrame(p, off, "data-series")
	if err != nil {
		return off, err
	}
	cur := entryOff
	for i := int32(0); i < count; i++ {
		if cur+2 > mapEnd {
			return off, fmt.Errorf("cram: data-series map entry %d: truncated key", i)
		}
		key := dataSeriesKey{p[cur], p[cur+1]}
		cur += 2
		enc, next, perr := parseEncoding(p[:mapEnd], cur)
		if perr != nil {
			return off, fmt.Errorf("cram: data-series %s: %w", key, perr)
		}
		cur = next
		out[key] = enc
	}
	if cur > mapEnd {
		return off, fmt.Errorf("cram: data-series map entries overran their declared size")
	}
	return mapEnd, nil
}

// parseTagEncodingMap parses the CRAM tag-encoding map starting at off.
// Each entry's key is an ITF-8 integer whose three low bytes are the
// tag's two-character name and one-character value type; the key is
// followed by one parsed Encoding.
func parseTagEncodingMap(p []byte, off int, out map[tagKey]*Encoding) (int, error) {
	entryOff, mapEnd, count, err := mapFrame(p, off, "tag-encoding")
	if err != nil {
		return off, err
	}
	cur := entryOff
	for i := int32(0); i < count; i++ {
		id, n, terr := itf8At(p[:mapEnd], cur)
		if terr != nil {
			return off, fmt.Errorf("cram: tag-encoding map entry %d key: %w", i, terr)
		}
		cur += n
		u := uint32(id)
		key := tagKey{byte(u >> 16), byte(u >> 8), byte(u)}
		enc, next, perr := parseEncoding(p[:mapEnd], cur)
		if perr != nil {
			return off, fmt.Errorf("cram: tag %s: %w", key, perr)
		}
		cur = next
		out[key] = enc
	}
	if cur > mapEnd {
		return off, fmt.Errorf("cram: tag-encoding map entries overran their declared size")
	}
	return mapEnd, nil
}
