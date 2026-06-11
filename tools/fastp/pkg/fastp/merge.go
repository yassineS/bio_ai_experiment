// Overlap-driven merge writer, FASTA-list adapter trimming, and the
// upstream poly-X trimmer.
//
// This file ports the remaining fastp parity surface that builds on the
// existing OverlapAnalysis (overlap.go):
//
//   - mergeOverlappedPair: a verbatim Go port of upstream fastp's
//     OverlapAnalysis::merge (reference_code/fastp/src/overlapanalysis.cpp:
//     152-183), powering -m/--merge / --merged_out / --include_unmerged.
//   - trimByMultiSequences / trimBySequenceUpstream: verbatim ports of
//     AdapterTrimmer::trimByMultiSequences and AdapterTrimmer::trimBySequence
//     (reference_code/fastp/src/adaptertrimmer.cpp:47-170), powering
//     --adapter_fasta.
//   - trimPolyXUpstream: a verbatim port of PolyX::trimPolyX
//     (reference_code/fastp/src/polyx.cpp:49-116), powering --trim-poly-x
//     with its own --poly_x_min_len knob.
//   - matchWithOneInsertion: a verbatim port of
//     Matcher::matchWithOneInsertion (reference_code/fastp/src/matcher.cpp:
//     10-54), used by the adapter trimmer's gap-tolerant fallback.

package fastp

import (
	"io"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// atcgBases mirrors upstream's ATCG_BASES[] (common.h): index 0=A, 1=T,
// 2=C, 3=G. The order matches PolyX::trimPolyX's atcgNumbers[] buckets.
var atcgBases = [4]byte{'A', 'T', 'C', 'G'}

// mergeOverlappedPair builds a merged read from an overlapping pair, given
// the OverlapAnalysisResult produced by analyzeOverlapPair on R1's sequence
// and the reverse complement of R2's sequence. It is a verbatim Go port of
// upstream fastp's OverlapAnalysis::merge (overlapanalysis.cpp:152-183).
//
// The merged read takes R1's name with a " merged_<len1>_<len2>" suffix
// appended to its description, matching upstream byte-for-byte. Returns nil
// when ov is not overlapped (upstream returns NULL).
func mergeOverlappedPair(record1, record2 *fastq.Record, ov OverlapAnalysisResult) *fastq.Record {
	if !ov.Overlapped {
		return nil
	}
	ol := ov.OverlapLen

	// len1 = ol + max(0, offset); len2 = (offset>0) ? r2.len - ol : 0.
	len1 := ol
	if ov.Offset > 0 {
		len1 += ov.Offset
	}
	len2 := 0
	if ov.Offset > 0 {
		len2 = len(record2.Sequence) - ol
	}

	rcSeq2 := reverseComplement(string(record2.Sequence))
	rcQual2 := reverseSlice(record2.Quality)

	mergedSeq := make([]byte, 0, len1+len2)
	mergedQual := make([]byte, 0, len1+len2)
	mergedSeq = append(mergedSeq, record1.Sequence[:len1]...)
	mergedQual = append(mergedQual, record1.Quality[:len1]...)
	if ov.Offset > 0 {
		mergedSeq = append(mergedSeq, rcSeq2[ol:ol+len2]...)
		mergedQual = append(mergedQual, rcQual2[ol:ol+len2]...)
	}

	// Upstream: name = r1.name + " merged_" + len1 + "_" + len2.
	suffix := " merged_" + strconv.Itoa(len1) + "_" + strconv.Itoa(len2)
	return &fastq.Record{
		ID:          record1.ID,
		Description: record1.Description + suffix,
		Sequence:    mergedSeq,
		Quality:     mergedQual,
	}
}

// LoadAdapterFasta reads adapter sequences from a FASTA stream for
// --adapter_fasta. It mirrors upstream's Options::loadFastaAdapters
// (options.cpp:52-79): sequences shorter than 6 bp are skipped (the caller
// warns), and the returned slice is ordered by sorted full-header name
// because upstream iterates a std::map<string,string> in key order. It
// reuses the shared pkg/htsgo/fasta reader for parsing. Returns an empty
// slice when the file yields no usable adapters.
func LoadAdapterFasta(r io.Reader) (adapters []string, skipped []string) {
	records, err := fasta.NewReader(r).ReadAll()
	if err != nil {
		// On a malformed stream we still return whatever parsed so far;
		// upstream's FastaReader is similarly best-effort.
		return adapters, skipped
	}

	// Upstream keys contigs on the full header line, so sort by Description
	// to reproduce the std::map iteration order. Duplicate headers collapse
	// to the last value, matching std::map assignment.
	contigs := make(map[string]string, len(records))
	names := make([]string, 0, len(records))
	for _, rec := range records {
		if _, seen := contigs[rec.Description]; !seen {
			names = append(names, rec.Description)
		}
		contigs[rec.Description] = string(rec.Sequence)
	}
	sort.Strings(names)

	for _, n := range names {
		s := contigs[n]
		if len(s) >= 6 {
			adapters = append(adapters, s)
		} else {
			skipped = append(skipped, s)
		}
	}
	return adapters, skipped
}

// trimByMultiSequences applies every adapter in adapterList to a read,
// trimming the read in place. It is a verbatim port of upstream's
// AdapterTrimmer::trimByMultiSequences (adaptertrimmer.cpp:47-69).
//
// The minimum match length (matchReq) scales with the adapter-list size:
// 4 by default, 5 for >16 adapters, 6 for >256. Returns the number of
// adapter-trimmed bases and whether anything was trimmed.
func trimByMultiSequences(record *fastq.Record, adapterList []string, isR2 bool) (trimmedBases int, trimmed bool) {
	matchReq := 4
	if len(adapterList) > 16 {
		matchReq = 5
	}
	if len(adapterList) > 256 {
		matchReq = 6
	}
	originalLen := len(record.Sequence)
	for _, adapter := range adapterList {
		if trimBySequenceUpstream(record, adapter, isR2, matchReq) {
			trimmed = true
		}
	}
	if trimmed {
		trimmedBases = originalLen - len(record.Sequence)
	}
	return trimmedBases, trimmed
}

// trimBySequenceUpstream trims a single adapter sequence from a read,
// resizing record.Sequence/record.Quality in place. It is a verbatim port
// of upstream's AdapterTrimmer::trimBySequence (adaptertrimmer.cpp:71-170),
// including the negative-start A-tailing offset and the one-gap
// insertion/deletion fallbacks. Returns true when a trim was applied.
func trimBySequenceUpstream(record *fastq.Record, adapterseq string, isR2 bool, matchReq int) bool {
	const allowOneMismatchForEach = 8

	rlen := len(record.Sequence)
	alen := len(adapterseq)
	adata := adapterseq
	rdata := string(record.Sequence)

	if alen < matchReq {
		return false
	}

	pos := 0
	found := false
	start := 0
	if alen >= 16 {
		start = -4
	} else if alen >= 12 {
		start = -3
	} else if alen >= 8 {
		start = -2
	}

	// Exact match with hamming distance (no insertion or deletion).
	for pos = start; pos < rlen-matchReq; pos++ {
		cmplen := rlen - pos
		if alen < cmplen {
			cmplen = alen
		}
		allowedMismatch := cmplen / allowOneMismatchForEach
		mismatch := 0
		matched := true
		iStart := 0
		if -pos > iStart {
			iStart = -pos
		}
		for i := iStart; i < cmplen; i++ {
			if adata[i] != rdata[i+pos] {
				mismatch++
				if mismatch > allowedMismatch {
					matched = false
					break
				}
			}
		}
		if matched {
			found = true
			break
		}
	}

	// Fallback: one insertion in the sequence. Upstream passes the read
	// pointer rdata at offset 0 (NOT rdata+pos); the loop variable pos only
	// shrinks cmplen, it does not advance the compared window
	// (adaptertrimmer.cpp:120-130).
	if !found {
		for pos = 0; pos < rlen-matchReq-1; pos++ {
			cmplen := rlen - pos - 1
			if alen < cmplen {
				cmplen = alen
			}
			allowedMismatch := cmplen/allowOneMismatchForEach - 1
			if matchWithOneInsertion(rdata, adata, cmplen, allowedMismatch) {
				found = true
				break
			}
		}
	}

	// Fallback: one deletion in the sequence. As above, upstream passes
	// rdata at offset 0 as normalData (adaptertrimmer.cpp:136-147).
	if !found {
		for pos = 0; pos < rlen-matchReq; pos++ {
			cmplen := rlen - pos
			if alen-1 < cmplen {
				cmplen = alen - 1
			}
			allowedMismatch := cmplen/allowOneMismatchForEach - 1
			if matchWithOneInsertion(adata, rdata, cmplen, allowedMismatch) {
				found = true
				break
			}
		}
	}

	if found {
		if pos < 0 {
			// Whole read is adapter dimer; clear it.
			record.Sequence = record.Sequence[:0]
			record.Quality = record.Quality[:0]
		} else {
			record.Sequence = record.Sequence[:pos]
			record.Quality = record.Quality[:pos]
		}
		return true
	}
	return false
}

// matchWithOneInsertion is a verbatim port of upstream's
// Matcher::matchWithOneInsertion (matcher.cpp:10-54). It reports whether
// insData matches normalData over cmplen bases allowing a single insertion
// in insData, within diffLimit mismatches. insData must have at least
// cmplen+1 bytes available and normalData at least cmplen.
func matchWithOneInsertion(insData, normalData string, cmplen, diffLimit int) bool {
	if cmplen <= 0 {
		return false
	}
	if len(insData) < cmplen+1 || len(normalData) < cmplen {
		return false
	}
	accMismatchFromLeft := make([]int, cmplen)
	accMismatchFromRight := make([]int, cmplen)

	if insData[0] == normalData[0] {
		accMismatchFromLeft[0] = 0
	} else {
		accMismatchFromLeft[0] = 1
	}
	if insData[cmplen] == normalData[cmplen-1] {
		accMismatchFromRight[cmplen-1] = 0
	} else {
		accMismatchFromRight[cmplen-1] = 1
	}
	for i := 1; i < cmplen; i++ {
		if insData[i] != normalData[i] {
			accMismatchFromLeft[i] = accMismatchFromLeft[i-1] + 1
		} else {
			accMismatchFromLeft[i] = accMismatchFromLeft[i-1]
		}
		if accMismatchFromLeft[i]+accMismatchFromRight[cmplen-1] > diffLimit {
			break
		}
	}
	for i := cmplen - 2; i >= 0; i-- {
		if insData[i+1] != normalData[i] {
			accMismatchFromRight[i] = accMismatchFromRight[i+1] + 1
		} else {
			accMismatchFromRight[i] = accMismatchFromRight[i+1]
		}
		if accMismatchFromRight[i]+accMismatchFromLeft[0] > diffLimit {
			for p := 0; p < i; p++ {
				accMismatchFromRight[p] = diffLimit + 1
			}
			break
		}
	}
	for i := 1; i < cmplen; i++ {
		if accMismatchFromLeft[i-1]+accMismatchFromRight[cmplen-1] > diffLimit {
			return false
		}
		if accMismatchFromLeft[i-1]+accMismatchFromRight[i] <= diffLimit {
			return true
		}
	}
	return false
}

// trimPolyXUpstream trims a 3' poly-X run from seq and returns the new
// length (the index at which to truncate). It is a verbatim port of
// upstream fastp's PolyX::trimPolyX (reference_code/fastp/src/polyx.cpp:
// 49-116). Unlike the older naive consecutive-base counter, it tolerates 1
// mismatch per 8 bases scanned (capped at 5), counts N as every base, and
// picks the most frequent of A/T/C/G as the poly base.
//
// Returns len(seq) when no trim should be applied, and the count of trimmed
// bases via the second return value.
func trimPolyXUpstream(seq string, compareReq int) (newLen int, trimmedBases int) {
	const allowOneMismatchForEach = 8
	const maxMismatch = 5

	rlen := len(seq)
	var atcgNumbers [4]int
	pos := 0
	for pos = 0; pos < rlen; pos++ {
		switch seq[rlen-pos-1] {
		case 'A':
			atcgNumbers[0]++
		case 'T':
			atcgNumbers[1]++
		case 'C':
			atcgNumbers[2]++
		case 'G':
			atcgNumbers[3]++
		case 'N':
			atcgNumbers[0]++
			atcgNumbers[1]++
			atcgNumbers[2]++
			atcgNumbers[3]++
		}

		cmp := pos + 1
		allowedMismatch := cmp / allowOneMismatchForEach
		if maxMismatch < allowedMismatch {
			allowedMismatch = maxMismatch
		}

		needToBreak := true
		for b := 0; b < 4; b++ {
			if cmp-atcgNumbers[b] <= allowedMismatch {
				needToBreak = false
			}
		}
		if needToBreak && (pos >= allowOneMismatchForEach || pos+1 >= compareReq-1) {
			break
		}
	}

	if pos+1 >= compareReq {
		poly := 0
		maxCount := -1
		for b := 0; b < 4; b++ {
			if atcgNumbers[b] > maxCount {
				maxCount = atcgNumbers[b]
				poly = b
			}
		}
		polyBase := atcgBases[poly]
		// Upstream walks back to the last polyBase: `while(data[rlen-pos-1]
		// != polyBase && pos>=0) pos--`. Two index hazards exist in C that
		// are out-of-range in Go and must be handled to avoid a panic while
		// reproducing the observable result:
		//   - When the scan loop above never broke (e.g. an all-N tail),
		//     pos == rlen, so data[rlen-pos-1] == data[-1] (a byte before the
		//     buffer in C). Any non-polyBase byte there just keeps the
		//     walkback going.
		//   - At pos == -1, data[rlen] is C's NUL terminator (never
		//     polyBase), which then fails the pos>=0 test and ends the loop.
		// We model both by treating any index outside [0, rlen) as "not
		// polyBase". The loop still terminates because pos strictly
		// decreases and the pos>=0 guard ends it at pos == -1, yielding
		// newLen == rlen (no trim) for the no-literal-match case.
		for pos >= 0 {
			idx := rlen - pos - 1
			if idx >= 0 && idx < rlen && seq[idx] == polyBase {
				break
			}
			pos--
		}
		newLen = rlen - pos - 1
		return newLen, rlen - newLen
	}
	return rlen, 0
}
