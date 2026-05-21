package cram

// Quality-score binning is a lossy CRAM-encode transform: each Phred
// quality value is mapped through a small fixed lookup table to one of a
// handful of representative values before the writer stores it. Collapsing
// the quality alphabet from ~64 symbols to a few makes the QS data series
// far more compressible, which is a large file-size win at the cost of
// quality precision. Binning is always opt-in: the default RecordWriter
// performs no binning and is byte-for-byte lossless.
//
// The schemes here are the standard Illumina quality-recalibration tables
// published in Illumina's "Reducing Whole-Genome Data Storage Footprint"
// technical note. The 8-level table is the canonical scheme HiSeq 2500
// and later instruments apply; the 4-level and 2-level tables are the
// coarser variants the same note describes (the 2-level table matches the
// NovaSeq-style two-state binning).

// QualityBinning selects the lossy quality-binning scheme a RecordWriter
// applies to each record's QUAL before encoding. The zero value,
// BinningNone, performs no binning — the writer stays losslessly exact.
type QualityBinning int

const (
	// BinningNone disables quality binning. It is the default: quality
	// scores are stored verbatim and the CRAM file round-trips exactly.
	BinningNone QualityBinning = iota
	// BinningIllumina8 applies Illumina's canonical 8-level quality
	// binning, the scheme HiSeq 2500+ instruments use. Quality values are
	// collapsed to eight representatives (0, 6, 15, 22, 27, 33, 37, 40).
	BinningIllumina8
	// BinningIllumina4 applies Illumina's 4-level quality binning, a
	// coarser variant collapsing quality to four representatives
	// (0, 15, 25, 37).
	BinningIllumina4
	// BinningIllumina2 applies Illumina's 2-level quality binning, the
	// NovaSeq-style two-state scheme collapsing quality to two
	// representatives (6, 37).
	BinningIllumina2
)

// String returns a short human-readable name for the binning scheme.
func (b QualityBinning) String() string {
	switch b {
	case BinningIllumina8:
		return "illumina-8"
	case BinningIllumina4:
		return "illumina-4"
	case BinningIllumina2:
		return "illumina-2"
	default:
		return "none"
	}
}

// valid reports whether b is a known binning scheme.
func (b QualityBinning) valid() bool {
	return b >= BinningNone && b <= BinningIllumina2
}

// binBoundary describes one bin of a binning scheme: every input quality
// in [lo, hi] maps to the representative value rep.
type binBoundary struct {
	lo, hi, rep byte
}

// illumina8Bins is Illumina's canonical 8-level quality-binning table.
// The boundaries and representatives are taken directly from the Illumina
// "Reducing Whole-Genome Data Storage Footprint" technical note.
var illumina8Bins = []binBoundary{
	{0, 2, 0},
	{3, 9, 6},
	{10, 19, 15},
	{20, 24, 22},
	{25, 29, 27},
	{30, 34, 33},
	{35, 39, 37},
	{40, 255, 40},
}

// illumina4Bins is Illumina's 4-level quality-binning table.
var illumina4Bins = []binBoundary{
	{0, 9, 0},
	{10, 19, 15},
	{20, 29, 25},
	{30, 255, 37},
}

// illumina2Bins is Illumina's 2-level (NovaSeq-style) quality-binning
// table.
var illumina2Bins = []binBoundary{
	{0, 14, 6},
	{15, 255, 37},
}

// buildBinTable expands a bin-boundary list into a dense 256-entry lookup
// table. The boundaries must cover the whole 0-255 range with no gaps.
func buildBinTable(bins []binBoundary) [256]byte {
	var t [256]byte
	for _, b := range bins {
		for q := int(b.lo); q <= int(b.hi); q++ {
			t[q] = b.rep
		}
	}
	return t
}

// The dense 256-entry lookup tables, built once at package init from the
// boundary lists above. Each maps every possible quality byte to its
// scheme representative.
var (
	illumina8Table = buildBinTable(illumina8Bins)
	illumina4Table = buildBinTable(illumina4Bins)
	illumina2Table = buildBinTable(illumina2Bins)
)

// table returns the 256-entry lookup table for the binning scheme. For
// BinningNone (and any unknown scheme) it returns an identity table, so
// callers can map unconditionally without a special case.
func (b QualityBinning) table() [256]byte {
	switch b {
	case BinningIllumina8:
		return illumina8Table
	case BinningIllumina4:
		return illumina4Table
	case BinningIllumina2:
		return illumina2Table
	default:
		var identity [256]byte
		for i := range identity {
			identity[i] = byte(i)
		}
		return identity
	}
}

// BinQuality maps every byte of qual through the scheme's lookup table and
// returns a new slice holding the binned values. The input slice is never
// modified, so a caller's sam.Record quality stays untouched. For
// BinningNone the returned slice is an exact copy of the input. The SAM
// "no quality" sentinel (0xff) is passed through unchanged by every
// scheme, so an absent-quality record stays absent.
func (b QualityBinning) BinQuality(qual []byte) []byte {
	out := make([]byte, len(qual))
	if len(qual) == 0 {
		return out
	}
	tbl := b.table()
	for i, q := range qual {
		if q == 0xff {
			// Preserve the SAM no-quality sentinel rather than binning it
			// into a real representative; this keeps round-tripping of
			// records that carry no QUAL exact.
			out[i] = 0xff
			continue
		}
		out[i] = tbl[q]
	}
	return out
}
