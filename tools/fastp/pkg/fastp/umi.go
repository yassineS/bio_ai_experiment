// UMI (Unique Molecular Identifier) extraction. Supports the same set of
// locations as upstream fastp:
//
//   - read1   : UMI is the prefix of the read 1 sequence
//   - read2   : UMI is the prefix of the read 2 sequence
//   - per_read: UMI is the prefix of BOTH read 1 and read 2 (PE only); the
//     R1 and R2 UMIs are joined with "_" before being added to the name
//   - index1  : UMI comes from the i7 index in the Illumina header
//     (everything before the first '+' in the last colon-separated field)
//   - index2  : UMI comes from the i5 index (everything after the '+')
//   - per_index: UMI is "i7_i5" assembled from both index fields
//
// In read/per_read modes the UMI bases (plus --umi_skip trailing bases)
// are removed from the sequence and quality of the affected record. In
// index/per_index modes the read sequence is left alone.
//
// The UMI is appended to the read name as ":<umi>" (no prefix set) or
// ":<prefix>_<umi>" (when --umi_prefix is provided), matching upstream
// fastp's umiprocessor.cpp.

package fastp

import (
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// UMI location identifiers (string values for the --umi_loc flag).
const (
	UMILocRead1    = "read1"
	UMILocRead2    = "read2"
	UMILocPerRead  = "per_read"
	UMILocIndex1   = "index1"
	UMILocIndex2   = "index2"
	UMILocPerIndex = "per_index"
)

// ValidUMILocation reports whether loc is one of the supported UMI
// location strings.
func ValidUMILocation(loc string) bool {
	switch loc {
	case UMILocRead1, UMILocRead2, UMILocPerRead, UMILocIndex1, UMILocIndex2, UMILocPerIndex:
		return true
	}
	return false
}

// applyUMI extracts the UMI from r1 (and r2 for paired-end) according to
// opts, appends it to the affected record names, and returns the updated
// records. The stats counter UMIProcessed is incremented once per
// successful extraction (counted per record, so PE per_read mode adds 2).
//
// When the read is shorter than UMILen+UMISkip, no UMI is extracted from
// that record and the record is returned unchanged; the stats counter is
// not incremented for it.
//
// For index-based modes (index1, index2, per_index) the original record
// sequence/quality is preserved and the UMI is parsed from the description.
func applyUMI(r1, r2 *fastq.Record, opts ProcessOptions, stats *ProcessStats) (*fastq.Record, *fastq.Record) {
	if !opts.UMI {
		return r1, r2
	}
	loc := opts.UMILoc
	if loc == "" {
		if r2 != nil {
			loc = UMILocPerRead
		} else {
			loc = UMILocRead1
		}
	}
	switch loc {
	case UMILocRead1:
		r1, ok := extractUMIFromSequence(r1, opts)
		if ok {
			stats.UMIProcessed++
			stats.UMIExtracted++
		}
		return r1, r2
	case UMILocRead2:
		r2, ok := extractUMIFromSequence(r2, opts)
		if ok {
			stats.UMIProcessed++
			stats.UMIExtracted++
		}
		return r1, r2
	case UMILocPerRead:
		// In per_read mode the R1 and R2 UMIs are joined into a single
		// UMI label so both records carry the same name suffix.
		umi1, r1New, ok1 := sliceUMI(r1, opts)
		umi2, r2New, ok2 := sliceUMI(r2, opts)
		if ok1 || ok2 {
			combined := umi1
			if ok1 && ok2 {
				combined = umi1 + "_" + umi2
			} else if !ok1 && ok2 {
				combined = umi2
			}
			r1New = appendUMIName(r1New, combined, opts.UMIPrefix)
			r2New = appendUMIName(r2New, combined, opts.UMIPrefix)
			if ok1 {
				stats.UMIProcessed++
			}
			if ok2 {
				stats.UMIProcessed++
			}
			stats.UMIExtracted += boolToInt(ok1) + boolToInt(ok2)
		}
		return r1New, r2New
	case UMILocIndex1, UMILocIndex2, UMILocPerIndex:
		i1, i2 := parseIndices(r1)
		var umi string
		switch loc {
		case UMILocIndex1:
			umi = i1
		case UMILocIndex2:
			umi = i2
		case UMILocPerIndex:
			switch {
			case i1 != "" && i2 != "":
				umi = i1 + "_" + i2
			case i1 != "":
				umi = i1
			case i2 != "":
				umi = i2
			}
		}
		if umi == "" {
			return r1, r2
		}
		r1New := appendUMIName(r1, umi, opts.UMIPrefix)
		var r2New *fastq.Record
		if r2 != nil {
			r2New = appendUMIName(r2, umi, opts.UMIPrefix)
			stats.UMIProcessed += 2
			stats.UMIExtracted += 2
		} else {
			stats.UMIProcessed++
			stats.UMIExtracted++
		}
		return r1New, r2New
	}
	return r1, r2
}

// sliceUMI returns the UMI bases extracted from the front of record's
// sequence (after honouring opts.UMISkip) along with a NEW record whose
// sequence/quality has those bases removed. If the record is too short
// or the UMI length is non-positive, the record is returned unchanged
// and ok is false.
//
// sliceUMI does NOT append the UMI to the name; callers do that so they
// can combine UMIs (per_read mode).
func sliceUMI(record *fastq.Record, opts ProcessOptions) (umi string, out *fastq.Record, ok bool) {
	if record == nil || opts.UMILen <= 0 {
		return "", record, false
	}
	end := opts.UMILen
	if end > len(record.Sequence) {
		return "", record, false
	}
	umi = string(record.Sequence[:end])

	skip := opts.UMISkip
	cut := end + skip
	if cut > len(record.Sequence) {
		cut = len(record.Sequence)
	}
	newSeq := make([]byte, len(record.Sequence)-cut)
	copy(newSeq, record.Sequence[cut:])
	newQual := make([]byte, len(record.Quality)-cut)
	copy(newQual, record.Quality[cut:])
	out = &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    newSeq,
		Quality:     newQual,
	}
	return umi, out, true
}

// extractUMIFromSequence is the convenience single-record wrapper around
// sliceUMI plus appendUMIName, used by the read1/read2 modes.
func extractUMIFromSequence(record *fastq.Record, opts ProcessOptions) (*fastq.Record, bool) {
	umi, out, ok := sliceUMI(record, opts)
	if !ok {
		return record, false
	}
	return appendUMIName(out, umi, opts.UMIPrefix), true
}

// appendUMIName appends the upstream-fastp-style UMI tag to record's ID
// and Description. Upstream's tag format is:
//
//	prefix == ""    -> ":<umi>"
//	prefix != ""    -> ":<prefix>_<umi>"
//
// (the leading colon is upstream's default delimiter — see
// umiprocessor.cpp `addUmiToName`). If umi is empty the record is
// returned unchanged.
func appendUMIName(record *fastq.Record, umi, prefix string) *fastq.Record {
	if record == nil || umi == "" {
		return record
	}
	var tag string
	if prefix == "" {
		tag = ":" + umi
	} else {
		tag = ":" + prefix + "_" + umi
	}
	out := *record
	out.ID = record.ID + tag
	if record.Description != "" {
		// Description is the full header line minus '@'; insert the tag
		// after the ID portion so downstream readers still see the
		// original metadata.
		if strings.HasPrefix(record.Description, record.ID) {
			rest := record.Description[len(record.ID):]
			out.Description = record.ID + tag + rest
		} else {
			out.Description = record.Description + tag
		}
	}
	return &out
}

// parseIndices extracts the i7 and i5 index strings from the Illumina-
// style header carried in record.Description. The standard header looks
// like "READID 1:N:0:ATCACG+CGATGT"; we look for the last colon-
// separated field of the second whitespace-separated token and split it
// on '+'. Returns ("", "") if no index is present.
func parseIndices(record *fastq.Record) (string, string) {
	if record == nil {
		return "", ""
	}
	header := record.Description
	if header == "" {
		header = record.ID
	}
	// Take the part after the first whitespace; that's the Illumina
	// metadata field (e.g. "1:N:0:ATCACG+CGATGT"). If no whitespace, fall
	// back to the whole string.
	if sp := strings.IndexAny(header, " \t"); sp >= 0 {
		header = header[sp+1:]
	}
	// The index is the last colon-separated field.
	colon := strings.LastIndex(header, ":")
	if colon < 0 {
		return "", ""
	}
	idx := header[colon+1:]
	if plus := strings.Index(idx, "+"); plus >= 0 {
		return idx[:plus], idx[plus+1:]
	}
	return idx, ""
}

func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}
