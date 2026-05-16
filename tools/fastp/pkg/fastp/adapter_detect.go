// Adapter detection helpers.
//
// Two strategies are provided:
//
//   - DetectAdapterSE: single-end. Delegates to detectAdapterSEUpstream
//     in adapter_autodetect.go, which is a verbatim port of upstream
//     fastp's Evaluator::evalAdapterAndReadNum + checkKnownAdapters +
//     getAdapterWithSeed + NucleotideTree (reference_code/fastp/src/
//     evaluator.cpp, nucleotidetree.cpp). Returns "" if fewer than
//     10000 records are available (the upstream gate at
//     evaluator.cpp:344) — this matches upstream's "No adapter
//     detected" exit path byte-for-byte.
//   - DetectAdapterPE / DetectAdaptersFromPairs: paired-end, overlap based.
//     Aligns R1's 3' end against the reverse-complement of R2's 5' end with
//     up to a few mismatches. If a confident overlap is found and R1 has
//     unaligned 3' tail beyond the overlap, that tail is the candidate
//     adapter. The most common candidate across all pairs wins.
//
// Both helpers return "" when no clear candidate is found.

package fastp

import (
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// adapterDetectSampleSize is the maximum number of reads (or pairs) we
// buffer from the input stream to drive adapter detection. Upstream fastp
// uses 256k (evaluator.cpp:301); we use the same.
const adapterDetectSampleSize = adapterDetectReadLimit

// peCandidateMinCount is the minimum count a PE-overlap candidate must
// hit before we'll accept it as the consensus adapter for that pair.
// This threshold is local to DetectAdaptersFromPairs; the SE path uses
// upstream's evaluator gates instead.
const peCandidateMinCount = 10

// DetectAdapterSE detects an adapter sequence from a sample of
// single-end reads. This is a thin wrapper around
// detectAdapterSEUpstream which is the verbatim port of upstream
// fastp's SE adapter auto-detect algorithm (evaluator.cpp:295-526).
//
// Returns the detected adapter sequence, or "" if no clear candidate
// emerges. The upstream gate at evaluator.cpp:344 short-circuits when
// fewer than 10000 records are supplied — this is intentional and
// required for parity (upstream prints "No adapter detected for
// read1" and proceeds without adapter trimming).
func DetectAdapterSE(reads []*fastq.Record) string {
	if len(reads) == 0 {
		return ""
	}
	return detectAdapterSEUpstream(reads)
}

// DetectAdapterPE detects an adapter sequence from a single paired-end
// read pair using overlap analysis.
//
// In paired-end read-through sequencing, an insert shorter than the read
// length appears in both reads. After reverse-complementing R2 we
// expect:
//
//	R1            = insert + R1_adapter
//	revcomp(R2)   = revcomp(R2_adapter) + insert
//
// So the insert in R1 starts at index 0; the insert in revcomp(R2)
// starts at some offset k (where k = |R2_adapter|). When we find that
// offset k, the unaligned 3' tail of R1 is the R1 adapter.
//
// Returns "" if no overlap is found or R1's 3' tail is empty.
func DetectAdapterPE(r1, r2 *fastq.Record) string {
	if r1 == nil || r2 == nil {
		return ""
	}
	const minOverlap = 30
	const maxMismatchesPE = 5

	s1 := string(r1.Sequence)
	s2 := reverseComplement(string(r2.Sequence))
	if len(s1) < minOverlap || len(s2) < minOverlap {
		return ""
	}

	// Scan offsets k in s2 (the R2-adapter length we'd be trimming).
	// Larger k = shorter insert = longer R1 adapter tail.
	bestK := -1
	bestMismatches := maxMismatchesPE + 1
	bestOverlap := 0
	for k := 0; k <= len(s2)-minOverlap; k++ {
		overlapLen := len(s2) - k
		if overlapLen > len(s1) {
			overlapLen = len(s1)
		}
		if overlapLen < minOverlap {
			continue
		}
		mm := 0
		for i := 0; i < overlapLen; i++ {
			if s1[i] != s2[k+i] {
				mm++
				if mm > bestMismatches {
					break
				}
			}
		}
		if mm < bestMismatches || (mm == bestMismatches && overlapLen > bestOverlap) {
			bestMismatches = mm
			bestK = k
			bestOverlap = overlapLen
		}
	}

	if bestK < 0 || bestMismatches > maxMismatchesPE {
		return ""
	}

	// The R1 adapter tail is whatever in R1 extends beyond the overlap.
	tailStart := bestOverlap
	if tailStart >= len(s1) {
		// No tail - insert >= read length, no read-through.
		return ""
	}
	return s1[tailStart:]
}

// DetectAdaptersFromPairs runs DetectAdapterPE on each pair, then
// returns the most-common detected R1 adapter and (mirrored) R2 adapter.
// The returned strings are "" if no clear consensus emerges.
//
// The R2 adapter is derived by taking R2's 3' unaligned tail in the same
// way (i.e. running the algorithm symmetrically): for each pair we
// reverse-complement R1 and check against R2.
func DetectAdaptersFromPairs(pairs [][2]*fastq.Record) (r1Adapter, r2Adapter string) {
	r1Candidates := map[string]int{}
	r2Candidates := map[string]int{}
	for _, p := range pairs {
		if a := DetectAdapterPE(p[0], p[1]); a != "" {
			// Use a prefix-trim key so very-long unrelated tails don't
			// each look unique; we count by the first 20 bases.
			key := a
			if len(key) > 20 {
				key = key[:20]
			}
			r1Candidates[key]++
		}
		// Symmetric: swap roles to find R2 adapter.
		if a := DetectAdapterPE(p[1], p[0]); a != "" {
			key := a
			if len(key) > 20 {
				key = key[:20]
			}
			r2Candidates[key]++
		}
	}
	r1Adapter = pickConsensus(r1Candidates, len(pairs))
	r2Adapter = pickConsensus(r2Candidates, len(pairs))
	return r1Adapter, r2Adapter
}

// pickConsensus returns the most common candidate that appears in at
// least max(adapterMinExtendCount, 1% of total) pairs.
func pickConsensus(candidates map[string]int, total int) string {
	if len(candidates) == 0 {
		return ""
	}
	threshold := total / 100
	if threshold < peCandidateMinCount {
		threshold = peCandidateMinCount
	}
	// Pick top candidate.
	bestSeq := ""
	bestCount := 0
	for s, c := range candidates {
		if c > bestCount || (c == bestCount && s < bestSeq) {
			bestSeq = s
			bestCount = c
		}
	}
	if bestCount < threshold {
		return ""
	}
	return bestSeq
}
