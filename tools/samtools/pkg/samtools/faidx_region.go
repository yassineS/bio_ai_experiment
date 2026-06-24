package samtools

// Region parsing for faidx/fqidx. This reproduces htslib's hts_parse_region
// (the reference-id-aware variant used by fai_parse_region) closely enough to
// match `samtools faidx` byte-for-byte on the region forms it accepts:
//
//	chr            whole contig
//	chr:           whole contig
//	chr:0          whole contig (0 == "from the start")
//	chr:N          chr:N-<end of contig>
//	chr:N-         chr:N-<end of contig>
//	chr:-M         chr:1-M
//	chr:N-M        N..M inclusive (1-based)
//
// Coordinates accept thousands separators (commas) and k/m/g suffixes, like
// hts_parse_decimal. The whole string is tried as a contig name first (so a
// reference literally named "chr1:100" resolves to itself), matching htslib's
// "check whole name first" behaviour. Past-the-end positions are clamped here
// the same way fai_get_val does, and the unclamped request bounds are retained
// so the caller can reproduce the "Truncated sequence" warning.

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// htsPosMax mirrors htslib's HTS_POS_MAX sentinel for an open-ended region.
const htsPosMax = int64(^uint64(0) >> 1) // INT64_MAX

// parsedRegion is the resolved, clamped form of a region string ready for a
// FetchRaw call, plus the pre-clamp request bounds for truncation reporting.
type parsedRegion struct {
	name      string // resolved contig name
	beg0      int64  // clamped, 0-based inclusive start
	end0      int64  // clamped, 0-based exclusive end
	lineBases int64  // contig line length (for the "same as input" wrap)

	// reqBeg0 / reqEnd0 are the pre-clamp 0-based half-open bounds; used to
	// detect truncation. hasExplicitEnd is true when an explicit end coordinate
	// was supplied (i.e. end != HTS_POS_MAX), matching upstream's
	// `end < HTS_POS_MAX` guard.
	reqBeg0        int64
	reqEnd0        int64
	hasExplicitEnd bool
}

// parseFaidxRegion parses region against the supplied index. It returns
// ok == false when the contig cannot be found or the region is malformed
// (htslib returns NULL / tid < 0 in those cases), which the caller turns into
// the "not found / Failed to fetch" diagnostics.
func parseFaidxRegion(idx *fasta.Index, region string) (parsedRegion, bool) {
	tid, beg, end, ok := htsParseRegion(idx, region)
	if !ok {
		return parsedRegion{}, false
	}
	entry := idx.Entries()[tid]

	// Clamp exactly as fai_get_val: beg/end >= len snap to len, then beg>end
	// snaps beg to end.
	cbeg, cend := beg, end
	if cbeg >= entry.Length {
		cbeg = entry.Length
	}
	if cend >= entry.Length {
		cend = entry.Length
	}
	if cbeg > cend {
		cbeg = cend
	}

	return parsedRegion{
		name:           entry.Name,
		beg0:           cbeg,
		end0:           cend,
		lineBases:      entry.LineBases,
		reqBeg0:        beg,
		reqEnd0:        end,
		hasExplicitEnd: end < htsPosMax,
	}, true
}

// htsParseRegion is the core port of htslib's hts_parse_region for the FASTA
// index. It returns the contig index into idx.Entries(), the 0-based
// half-open [beg, end) bounds (end == htsPosMax for open-ended), and ok.
func htsParseRegion(idx *fasta.Index, s string) (tid int, beg, end int64, ok bool) {
	if s == "" {
		return -1, 0, 0, false
	}

	// Braced quoting "{name}" / "{name}:reg" disambiguation.
	if s[0] == '{' {
		closeIdx := indexByte(s, '}')
		if closeIdx < 0 {
			return -1, 0, 0, false
		}
		name := s[1:closeIdx]
		if closeIdx+1 < len(s) && s[closeIdx+1] == ':' {
			// Quoted with coordinates.
			id := lookupName(idx, name)
			if id < 0 {
				return -1, 0, 0, false
			}
			b, e, perr := parseCoords(s[closeIdx+2:])
			if perr {
				return -1, 0, 0, false
			}
			return id, b, e, true
		}
		// "{name}" — whole contig.
		id := lookupName(idx, name)
		if id < 0 {
			return -1, 0, 0, false
		}
		return id, 0, htsPosMax, true
	}

	colon := lastIndexByte(s, ':')
	if colon < 0 {
		// No colon: whole-contig name lookup.
		id := lookupName(idx, s)
		if id < 0 {
			return -1, 0, 0, false
		}
		return id, 0, htsPosMax, true
	}

	// Check the whole string as a name first.
	if id := lookupName(idx, s); id >= 0 {
		// Whole name matches; ambiguity is possible only if the pre-colon part
		// is ALSO a name (htslib errors out). Reproduce that by failing.
		if lookupName(idx, s[:colon]) >= 0 {
			return -1, 0, 0, false
		}
		return id, 0, htsPosMax, true
	}

	// Pre-colon part must be a valid name.
	id := lookupName(idx, s[:colon])
	if id < 0 {
		return -1, 0, 0, false
	}
	b, e, perr := parseCoords(s[colon+1:])
	if perr {
		return -1, 0, 0, false
	}
	return id, b, e, true
}

// parseCoords parses the post-colon coordinate string, returning 0-based
// half-open [beg, end). It mirrors the coordinate-defaulting branches of
// hts_parse_region exactly. perr is true on a malformed/empty-result region
// (beg >= end, or trailing garbage), which htslib maps to a NULL return.
func parseCoords(rest string) (beg, end int64, perr bool) {
	// hts_parse_decimal on the start; beg = value - 1.
	startVal, after, ok := parseDecimal(rest)
	if !ok {
		// No digits: hts_parse_decimal returns 0 with digits==0; beg = 0-1 =
		// -1, hyphen == rest (unmoved).
		startVal = 0
		after = rest
	}
	beg = startVal - 1

	if beg < 0 {
		// beg < 0 branch in hts_parse_region.
		afterCh := byte(0)
		if len(after) > 0 {
			afterCh = after[0]
		}
		if beg != -1 && afterCh == '-' && len(rest) > 0 {
			// "User specified zero, but we're 1-based" → error.
			return 0, 0, true
		}
		if afterCh == 0 || isDigit(afterCh) || afterCh == ',' {
			// interpret chr:-100 as chr:1-100; chr:0 (beg==-1) → whole contig.
			if beg == -1 {
				end = htsPosMax
			} else {
				end = -(beg + 1)
			}
			beg = 0
			return beg, end, false
		}
		if beg < -1 {
			// Unexpected string after region.
			return 0, 0, true
		}
	}

	// beg >= 0: inspect what follows the start number.
	if len(after) == 0 {
		// chr:N → chr:N-<end>
		end = htsPosMax
	} else if after[0] == '-' {
		// chr:N-M or chr:N-
		endVal, after2, ok2 := parseDecimal(after[1:])
		if !ok2 {
			// chr:N- (empty end) → end stays 0 → set to POS_MAX below.
			endVal = 0
			after2 = after[1:]
		}
		if len(after2) > 0 {
			return 0, 0, true // trailing garbage after end
		}
		end = endVal
	} else {
		return 0, 0, true // unexpected string after region
	}

	if end == 0 {
		end = htsPosMax // chr:N- → chr:N-<end>
	}
	if beg >= end {
		return 0, 0, true
	}
	return beg, end, false
}

// parseDecimal is a faidx-scoped port of hts_parse_decimal (with thousands
// separators enabled, the faidx default). It returns the parsed value, the
// remainder of the string after the number, and ok == (digits > 0). Commas
// embedded in the digits are skipped; a leading sign is honoured; k/m/g and
// scientific suffixes are supported.
func parseDecimal(s string) (val int64, rest string, ok bool) {
	i := 0
	for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
		i++
	}
	sign := int64(1)
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		if s[i] == '-' {
			sign = -1
		}
		i++
	}
	var n uint64
	digits := 0
	for i < len(s) {
		c := s[i]
		if c >= '0' && c <= '9' {
			n = n*10 + uint64(c-'0')
			digits++
			i++
		} else if c == ',' {
			// Thousands separator: skip.
			i++
		} else {
			break
		}
	}
	decimals := 0
	if i < len(s) && s[i] == '.' {
		i++
		for i < len(s) && s[i] >= '0' && s[i] <= '9' {
			n = n*10 + uint64(s[i]-'0')
			decimals++
			digits++
			i++
		}
	}
	e := 0
	if i < len(s) {
		switch s[i] {
		case 'e', 'E':
			i++
			esign := 1
			if i < len(s) && (s[i] == '+' || s[i] == '-') {
				if s[i] == '-' {
					esign = -1
				}
				i++
			}
			for i < len(s) && s[i] >= '0' && s[i] <= '9' {
				e = e*10 + int(s[i]-'0')
				i++
			}
			e *= esign
		case 'k', 'K':
			e += 3
			i++
		case 'm', 'M':
			e += 6
			i++
		case 'g', 'G':
			e += 9
			i++
		}
	}
	e -= decimals
	for e > 0 {
		n *= 10
		e--
	}
	for e < 0 {
		n /= 10
		e++
	}
	if digits == 0 {
		return 0, s, false
	}
	return sign * int64(n), s[i:], true
}

// lookupName returns the contig index for name, or -1 if absent.
func lookupName(idx *fasta.Index, name string) int {
	return idx.Pos(name)
}

func indexByte(s string, c byte) int {
	for i := 0; i < len(s); i++ {
		if s[i] == c {
			return i
		}
	}
	return -1
}

func lastIndexByte(s string, c byte) int {
	for i := len(s) - 1; i >= 0; i-- {
		if s[i] == c {
			return i
		}
	}
	return -1
}
