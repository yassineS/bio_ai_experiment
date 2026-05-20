package cram

// SeriesKind classifies the value type of a CRAM data series, which
// determines how its EXTERNAL-encoded block is read: as a run of ITF-8
// integers or as raw bytes.
type SeriesKind int

// The CRAM data-series value kinds.
const (
	// SeriesInt is a series of integers. An EXTERNAL block carrying such
	// a series stores each value as an ITF-8 integer.
	SeriesInt SeriesKind = iota
	// SeriesByte is a series of single bytes or byte strings. An
	// EXTERNAL block carrying such a series stores the bytes verbatim,
	// one byte per value.
	SeriesByte
)

// seriesCatalogue records the value kind of every CRAM v3.0 data series,
// per the data-series table of the CRAM specification. A series absent
// from the catalogue defaults to SeriesInt, the common case.
//
// The byte-valued series are the ones whose EXTERNAL block holds raw
// bytes rather than ITF-8 integers: read bases (BA, BB), quality scores
// (QS, QQ), inserted/soft-clipped/stretch bytes and the base-call /
// mapping detail bytes. Reading one of these as ITF-8 would mis-parse
// any byte with the high bit set (a quality score above 0x80, for
// example), so the distinction is load-bearing.
var seriesCatalogue = map[string]SeriesKind{
	// Integer series (listed for completeness / documentation).
	"BF": SeriesInt, // bit flags
	"CF": SeriesInt, // CRAM flags
	"RI": SeriesInt, // reference id
	"RL": SeriesInt, // read lengths
	"AP": SeriesInt, // alignment positions
	"RG": SeriesInt, // read groups
	"MF": SeriesInt, // mate flags
	"NS": SeriesInt, // mate reference id
	"NP": SeriesInt, // mate alignment position
	"TS": SeriesInt, // template size
	"NF": SeriesInt, // distance to next fragment
	"TL": SeriesInt, // tag line (tag-combination index)
	"FN": SeriesInt, // number of read features
	"FP": SeriesInt, // read-feature position
	"DL": SeriesInt, // deletion length
	"RS": SeriesInt, // reference skip length
	"PD": SeriesInt, // padding length
	"HC": SeriesInt, // hard-clip length
	"MQ": SeriesInt, // mapping qualities
	"TC": SeriesInt, // tag counts (CRAM v1/v2 legacy)
	"TN": SeriesInt, // tag names (CRAM v1/v2 legacy)
	"TM": SeriesInt, // test marker
	"TV": SeriesInt, // test marker
	// Byte-valued series — EXTERNAL blocks hold raw bytes.
	"BA": SeriesByte, // single base
	"QS": SeriesByte, // quality score (one byte per base)
	"BB": SeriesByte, // stretches of bases
	"QQ": SeriesByte, // stretches of quality scores
	"IN": SeriesByte, // inserted bases
	"SC": SeriesByte, // soft-clipped bases
	"BS": SeriesByte, // base substitution codes
	"FC": SeriesByte, // read-feature codes (one ASCII byte each)
	"RN": SeriesByte, // read-name data series (distinct from the RN preservation-map key)
}

// SeriesValueKind returns the value kind of the named CRAM data series.
// A series not in the catalogue is treated as SeriesInt.
func SeriesValueKind(key string) SeriesKind {
	if k, ok := seriesCatalogue[key]; ok {
		return k
	}
	return SeriesInt
}
