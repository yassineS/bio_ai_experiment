// Adapter detection helpers.
//
// Two strategies are provided:
//
//   - DetectAdapterSE: single-end, k-mer-frequency based. Counts the
//     most common short k-mers near the 3' end of the first
//     ~adapterDetectSampleSize reads, then extends the leading k-mer
//     greedily as long as the next base is dominant (>=50% of cases).
//   - DetectAdapterPE / DetectAdaptersFromPairs: paired-end, overlap based.
//     Aligns R1's 3' end against the reverse-complement of R2's 5' end with
//     up to a few mismatches. If a confident overlap is found and R1 has
//     unaligned 3' tail beyond the overlap, that tail is the candidate
//     adapter. The most common candidate across all pairs wins.
//
// Both helpers return "" when no clear candidate is found.

package fastp

import (
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// adapterDetectSampleSize is the maximum number of reads (or pairs) we
// buffer from the input stream to drive adapter detection. Upstream fastp
// uses 256k; we use the same.
const adapterDetectSampleSize = 256 * 1024

// adapterKmerLen is the k-mer size used by DetectAdapterSE. fastp uses 10.
const adapterKmerLen = 10

// adapterMinExtendCount is the minimum number of supporting observations
// required to extend a candidate adapter by one more base in
// DetectAdapterSE.
const adapterMinExtendCount = 10

// adapterMaxAdapterLen is the maximum length to which a detected adapter
// will be extended, matching upstream fastp's behavior.
const adapterMaxAdapterLen = 60

// DetectAdapterSE detects an adapter sequence from a sample of
// single-end reads using k-mer frequencies. It scans the 3' portion of
// each read for k-mers, picks the most common k-mer as a seed, then
// greedily extends it base-by-base while a single base dominates the
// next position across all reads where the seed occurs.
//
// Returns the detected adapter sequence, or "" if no clear candidate
// emerges.
func DetectAdapterSE(reads []*fastq.Record) string {
	if len(reads) == 0 {
		return ""
	}

	// Step 1: count k-mers seen near the 3' end of reads. We look at
	// the last quarter of each read (a typical adapter is in the 3' end).
	type kmerCount struct {
		seq   string
		count int
	}
	counts := make(map[string]int, 1024)
	for _, r := range reads {
		if r == nil {
			continue
		}
		seq := string(r.Sequence)
		if len(seq) < adapterKmerLen {
			continue
		}
		// Scan the 3' half of the read.
		start := len(seq) / 2
		if start < 0 {
			start = 0
		}
		for i := start; i+adapterKmerLen <= len(seq); i++ {
			k := seq[i : i+adapterKmerLen]
			// Skip low-complexity / homopolymer k-mers.
			if isLowComplexityKmer(k) {
				continue
			}
			counts[k]++
		}
	}
	if len(counts) == 0 {
		return ""
	}

	// Step 2: pick the top k-mer. Require it to be at least 1% as common
	// as the dataset is large to avoid noise on tiny inputs.
	type kc struct {
		seq   string
		count int
	}
	all := make([]kc, 0, len(counts))
	for s, c := range counts {
		all = append(all, kc{s, c})
	}
	sort.Slice(all, func(i, j int) bool {
		if all[i].count != all[j].count {
			return all[i].count > all[j].count
		}
		return all[i].seq < all[j].seq
	})
	top := all[0]
	threshold := len(reads) / 100
	if threshold < adapterMinExtendCount {
		threshold = adapterMinExtendCount
	}
	if top.count < threshold {
		return ""
	}

	// Step 3a: greedy 3' extension. For each read containing the seed,
	// look at the base immediately following the seed; if one base
	// dominates with >= 50% support AND count >= adapterMinExtendCount,
	// append it and continue.
	adapter := top.seq
	for len(adapter) < adapterMaxAdapterLen {
		nextCounts := [5]int{}
		for _, r := range reads {
			if r == nil {
				continue
			}
			seq := string(r.Sequence)
			idx := indexOf(seq, adapter)
			if idx < 0 {
				continue
			}
			next := idx + len(adapter)
			if next >= len(seq) {
				continue
			}
			nextCounts[baseIndex(seq[next])]++
		}
		best, ok := dominantBase(nextCounts)
		if !ok {
			break
		}
		adapter += string("ACGT"[best])
	}

	// Step 3b: greedy 5' extension. The seed may have landed inside the
	// adapter rather than at its start; walk backward while a base
	// dominates the preceding position.
	for len(adapter) < adapterMaxAdapterLen {
		prevCounts := [5]int{}
		for _, r := range reads {
			if r == nil {
				continue
			}
			seq := string(r.Sequence)
			idx := indexOf(seq, adapter)
			if idx <= 0 {
				continue
			}
			prevCounts[baseIndex(seq[idx-1])]++
		}
		best, ok := dominantBase(prevCounts)
		if !ok {
			break
		}
		adapter = string("ACGT"[best]) + adapter
	}

	return adapter
}

// dominantBase returns the index in "ACGT" of the dominant base from a
// count vector (with index 4 reserved for N). It requires the dominant
// base to account for at least 50% of A+C+G+T support and to have at
// least adapterMinExtendCount observations.
func dominantBase(counts [5]int) (int, bool) {
	total := 0
	best := -1
	bestCount := 0
	for i := 0; i < 4; i++ {
		total += counts[i]
		if counts[i] > bestCount {
			bestCount = counts[i]
			best = i
		}
	}
	if total == 0 || bestCount < adapterMinExtendCount {
		return -1, false
	}
	if float64(bestCount)/float64(total) < 0.5 {
		return -1, false
	}
	return best, true
}

// isLowComplexityKmer returns true if a k-mer is degenerate (single
// base run or only two distinct bases dominating). Used to skip noisy
// seeds during DetectAdapterSE.
func isLowComplexityKmer(k string) bool {
	if len(k) == 0 {
		return true
	}
	freq := map[byte]int{}
	for i := 0; i < len(k); i++ {
		freq[k[i]]++
	}
	if len(freq) < 3 {
		return true
	}
	for _, c := range freq {
		if c*2 > len(k)*3/2 { // >= 75% one base
			return true
		}
	}
	return false
}

// indexOf returns the first index of needle in haystack, or -1.
func indexOf(haystack, needle string) int {
	n := len(needle)
	if n == 0 || n > len(haystack) {
		return -1
	}
	for i := 0; i+n <= len(haystack); i++ {
		if haystack[i:i+n] == needle {
			return i
		}
	}
	return -1
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
	if threshold < adapterMinExtendCount {
		threshold = adapterMinExtendCount
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
