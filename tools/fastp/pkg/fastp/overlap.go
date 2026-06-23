// Overlap analysis and overlap-based base correction for paired-end reads.
//
// This file is a verbatim Go port of upstream fastp's
// OverlapAnalysis::analyze (reference_code/fastp/src/overlapanalysis.cpp)
// and BaseCorrector::correctByOverlapAnalysis
// (reference_code/fastp/src/basecorrector.cpp). It powers the
// --correction flag (overlap-based base correction, PE only) and the
// overlap-length/diff knobs (--overlap_len_require, --overlap_diff_limit,
// --overlap_diff_percent_limit).

package fastp

import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"

// OverlapAnalysisResult mirrors upstream fastp's OverlapResult struct.
// Overlapped reports whether a confident overlap was found. Offset is the
// number of bases R2 (reverse-complemented) is shifted relative to R1; a
// positive offset means R1 has a 5' overhang, a negative offset means the
// insert is shorter than the read length (adapter present). OverlapLen is
// the number of overlapping bases and Diff is the number of mismatches in
// that overlap.
type OverlapAnalysisResult struct {
	Overlapped bool
	Offset     int
	OverlapLen int
	Diff       int
	HasGap     bool
}

// analyzeOverlapPair is a verbatim port of upstream fastp's
// OverlapAnalysis::analyze (overlapanalysis.cpp:16-150). It compares R1
// against the reverse complement of R2 over a sliding offset (forward then
// reverse), returning the first confident overlap.
//
// diffLimit is the absolute mismatch cap (--overlap_diff_limit, default 5),
// overlapRequire is the minimum overlap length (--overlap_len_require,
// default 30), and diffPercentLimit is the per-overlap mismatch fraction
// (--overlap_diff_percent_limit / 100, default 0.2). The allowGap path of
// upstream is not implemented here (upstream defaults it off and the
// correction path explicitly disables it), so this port covers the
// no-gap forward and reverse passes that drive --correction and merge.
func analyzeOverlapPair(seq1, rcSeq2 string, diffLimit, overlapRequire int, diffPercentLimit float64) OverlapAnalysisResult {
	len1 := len(seq1)
	len2 := len(rcSeq2)

	const completeCompareRequire = 50

	// forward with no gap: slide R2 to the right by offset.
	for offset := 0; offset < len1-overlapRequire; offset++ {
		overlapLen := len1 - offset
		if len2 < overlapLen {
			overlapLen = len2
		}
		overlapDiffLimit := diffLimit
		if l := int(float64(overlapLen) * diffPercentLimit); l < overlapDiffLimit {
			overlapDiffLimit = l
		}

		diff := 0
		i := 0
		for ; i < overlapLen; i++ {
			if seq1[offset+i] != rcSeq2[i] {
				diff++
				if diff > overlapDiffLimit && i < completeCompareRequire {
					break
				}
			}
		}

		if diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > completeCompareRequire) {
			return OverlapAnalysisResult{Overlapped: true, Offset: offset, OverlapLen: overlapLen, Diff: diff}
		}
	}

	// reverse with no gap: slide R2 to the left (negative offset). This is
	// the adapter-present case where the insert is shorter than the read.
	for offset := 0; offset > -(len2 - overlapRequire); offset-- {
		absOffset := -offset
		overlapLen := len1
		if v := len2 - absOffset; v < overlapLen {
			overlapLen = v
		}
		overlapDiffLimit := diffLimit
		if l := int(float64(overlapLen) * diffPercentLimit); l < overlapDiffLimit {
			overlapDiffLimit = l
		}

		diff := 0
		i := 0
		for ; i < overlapLen; i++ {
			if seq1[i] != rcSeq2[absOffset+i] {
				diff++
				if diff > overlapDiffLimit && i < completeCompareRequire {
					break
				}
			}
		}

		if diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > completeCompareRequire) {
			return OverlapAnalysisResult{Overlapped: true, Offset: offset, OverlapLen: overlapLen, Diff: diff}
		}
	}

	return OverlapAnalysisResult{Overlapped: false}
}

// trimByOverlapAnalysis ports AdapterTrimmer::trimByOverlapAnalysis
// (adaptertrimmer.cpp:16). When the PE overlap shows the insert is shorter than
// the reads (Offset < 0, i.e. read-through), it trims both mates down to the
// insert length, removing the read-through adapter from BOTH reads regardless of
// the adapter sequence. This is the primary paired-end adapter-removal mechanism
// and works even when R1 and R2 adapters differ (e.g. Nextera), which a single
// shared 3' adapter sequence cannot. It mutates record1/record2 in place and
// returns the trimmed adapter tails plus whether a trim occurred.
func trimByOverlapAnalysis(record1, record2 *fastq.Record, ov OverlapAnalysisResult, frontTrimmed1, frontTrimmed2 int) (adapter1, adapter2 string, trimmed bool) {
	if !ov.Overlapped || ov.Offset >= 0 {
		return "", "", false
	}
	ol := ov.OverlapLen
	len1 := len(record1.Sequence)
	if v := ol + frontTrimmed2; v < len1 {
		len1 = v
	}
	len2 := len(record2.Sequence)
	if v := ol + frontTrimmed1; v < len2 {
		len2 = v
	}
	adapter1 = string(record1.Sequence[len1:])
	adapter2 = string(record2.Sequence[len2:])
	record1.Sequence = record1.Sequence[:len1]
	record1.Quality = record1.Quality[:len1]
	record2.Sequence = record2.Sequence[:len2]
	record2.Quality = record2.Quality[:len2]
	return adapter1, adapter2, true
}

// correctByOverlapAnalysis performs overlap-based base correction on a
// paired-end read pair. It is a verbatim port of
// BaseCorrector::correctByOverlapAnalysis (basecorrector.cpp:16-83): for
// each mismatched position inside the overlap, if one mate's base has high
// quality (>= Q30) and the other's is low (<= Q14), the low-quality base
// (and its quality) is overwritten with the complement of the high-quality
// base. It mutates record1/record2 in place and returns the number of
// corrected bases.
//
// ov must come from analyzeOverlapPair using R1's sequence and the reverse
// complement of R2's sequence, with allowGap disabled. Per upstream, no
// correction is performed when ov.Diff == 0 or the pair is not overlapped.
//
// It returns the number of corrected bases plus the number of corrected
// reads (1 if only one mate was edited, 2 if both were), matching upstream's
// FilterResult::incCorrectedReads accounting.
func correctByOverlapAnalysis(record1, record2 *fastq.Record, ov OverlapAnalysisResult, encoding fastq.QualityEncoding) (correctedBases, correctedReads int) {
	if ov.Diff == 0 || !ov.Overlapped {
		return 0, 0
	}
	offset := phredOffset(encoding)
	goodQual := byte(30 + offset)
	badQual := byte(14 + offset)

	ol := ov.OverlapLen
	start1 := 0
	if ov.Offset > 0 {
		start1 = ov.Offset
	}
	negOffset := 0
	if ov.Offset < 0 {
		negOffset = -ov.Offset
	}
	start2 := len(record2.Sequence) - negOffset - 1

	seq1 := record1.Sequence
	seq2 := record2.Sequence
	qual1 := record1.Quality
	qual2 := record2.Quality

	corrected := 0
	r1Corrected := false
	r2Corrected := false
	for i := 0; i < ol; i++ {
		p1 := start1 + i
		p2 := start2 - i
		if p1 < 0 || p1 >= len(seq1) || p2 < 0 || p2 >= len(seq2) {
			continue
		}
		if seq1[p1] != complementBase(seq2[p2]) {
			if qual1[p1] >= goodQual && qual2[p2] <= badQual {
				// use R1: rewrite R2's base/quality.
				seq2[p2] = complementBase(seq1[p1])
				qual2[p2] = qual1[p1]
				corrected++
				r2Corrected = true
			} else if qual2[p2] >= goodQual && qual1[p1] <= badQual {
				// use R2: rewrite R1's base/quality.
				seq1[p1] = complementBase(seq2[p2])
				qual1[p1] = qual2[p2]
				corrected++
				r1Corrected = true
			}
		}
	}
	if corrected > 0 {
		if r1Corrected && r2Corrected {
			correctedReads = 2
		} else {
			correctedReads = 1
		}
	}
	return corrected, correctedReads
}

// complementBase returns the Watson-Crick complement of a single base,
// matching upstream fastp's complement() (util.h). Non-ACGT bases
// (including N) map to N.
func complementBase(b byte) byte {
	switch b {
	case 'A', 'a':
		return 'T'
	case 'T', 't':
		return 'A'
	case 'C', 'c':
		return 'G'
	case 'G', 'g':
		return 'C'
	default:
		return 'N'
	}
}
