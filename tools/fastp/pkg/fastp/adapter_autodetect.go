// Single-end adapter auto-detection, ported verbatim from upstream
// fastp's Evaluator (reference_code/fastp/src/evaluator.cpp:295-526)
// and NucleotideTree (reference_code/fastp/src/nucleotidetree.cpp).
//
// The algorithm has three stages, all preserved here byte-for-byte:
//
//  1. Sample up to READ_LIMIT records (no more than BASE_LIMIT bases) from
//     the FASTQ. If fewer than 10000 valid records are available, return
//     "" — same gate as upstream evaluator.cpp:344.
//  2. checkKnownAdapters: linear scan of every known adapter against every
//     sampled read, allowing 1 mismatch per 16 bases of overlap. If a
//     known adapter accumulates enough hits, that adapter wins.
//     (evaluator.cpp:207-293).
//  3. getAdapterWithSeed: bit-packed kmer histogram of the 3' portion of
//     reads (skipping 20 bp from the 5' end and shiftTail = max(1,
//     trim.tail1) bp from the 3' end). Pick the top-N highest-frequency
//     kmers that pass low-complexity / GC-bias / GGGG filters; for each
//     such seed, build forward + backward nucleotide trees from the
//     reads where the seed occurs, extract the dominant path
//     (children[i].count / total >= 0.95 with total >= 50) and
//     concatenate as: reverse(backward) + seed + forward. If the
//     resulting string matches a known adapter exactly we return that
//     known adapter; otherwise if the dominant walk reached a leaf we
//     return the de-novo adapter. (evaluator.cpp:362-470, nucleotidetree.cpp).
//
// The "verbatim" claim: every magic number, threshold, and loop bound
// here matches the upstream C++. Variable names mirror the upstream
// names where reasonable. See the per-function header comments for
// upstream line refs.

package fastp

import (
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Constants ported from evaluator.cpp.
const (
	// adapterDetectKeylen is the kmer length used when building the kmer
	// histogram. Upstream uses 10 (evaluator.cpp:367).
	adapterDetectKeylen = 10

	// adapterDetectTopN is the number of top-frequency kmers we try as
	// seeds. Upstream uses 10 (evaluator.cpp:387).
	adapterDetectTopN = 10

	// adapterDetectFoldThreshold gates whether a seed kmer is
	// sufficiently over-represented relative to the average kmer to be
	// worth growing into an adapter. Upstream uses 20 (evaluator.cpp:432).
	adapterDetectFoldThreshold = 20

	// adapterDetectMaxSearchLength caps how far into a read we look for
	// the seed kmer when building the nucleotide tree. Upstream uses 500
	// (evaluator.cpp:475).
	adapterDetectMaxSearchLength = 500

	// adapterDetectMaxLen caps the assembled adapter length; longer
	// strings are truncated. Upstream uses 60 (evaluator.cpp:510, 454).
	adapterDetectMaxLen = 60

	// adapterDetectMinRecords is the minimum number of records required
	// to attempt detection. Upstream uses 10000 (evaluator.cpp:344).
	adapterDetectMinRecords = 10000

	// adapterDetectReadLimit / BaseLimit cap how many reads we hold in
	// memory for detection. Upstream uses 256*1024 and 151*256*1024
	// (evaluator.cpp:301-302).
	adapterDetectReadLimit = 256 * 1024
	adapterDetectBaseLimit = int64(151) * 256 * 1024

	// shiftTailDefault is max(1, trim.tail1). Our port doesn't expose
	// trim_tail1 so we use the upstream default of 1 (evaluator.cpp:364).
	adapterDetectShiftTail = 1

	// checkKnownAdaptersMaxReads caps how many reads we scan against the
	// known-adapter table. Upstream uses 100000 (evaluator.cpp:213).
	checkKnownAdaptersMaxReads = 100000

	// checkKnownAdaptersMaxBases corresponds to MAX_CHECK_BASES.
	checkKnownAdaptersMaxBases = checkKnownAdaptersMaxReads * 1000

	// checkKnownAdaptersMaxHit short-circuits once any single adapter
	// passes this many positive hits. Upstream uses 1000
	// (evaluator.cpp:216).
	checkKnownAdaptersMaxHit = 1000

	// checkKnownAdaptersMatchReq is the minimum overlap length required
	// to consider a hit. Upstream uses 8 (evaluator.cpp:218).
	checkKnownAdaptersMatchReq = 8

	// checkKnownAdaptersOneMismatchPer is the denominator in
	// "allowedMismatch = cmplen / allowOneMismatchForEach". Upstream
	// uses 16 (evaluator.cpp:219).
	checkKnownAdaptersOneMismatchPer = 16
)

// nucleotideNode is the node type of the bit-encoded prefix tree used by
// upstream's NucleotideTree (nucleotidetree.cpp). The fan-out is 8 (not
// 4) so we can index children by `seq[i] & 0x07`, matching upstream's
// memory layout exactly. This is wasteful but keeps the byte-for-byte
// tie-breaking identical when several letters round-trip to the same
// hash bucket (i.e. it doesn't).
type nucleotideNode struct {
	count    int
	base     byte
	children [8]*nucleotideNode
}

// nucleotideTree is a prefix tree (trie) over A/T/C/G with bit-encoded
// child slots (children[base & 0x07]). Upstream source:
// nucleotidetree.cpp:32-88.
type nucleotideTree struct {
	root *nucleotideNode
}

// newNucleotideTree allocates an empty tree (nucleotidetree.cpp:32-35).
func newNucleotideTree() *nucleotideTree {
	return &nucleotideTree{root: &nucleotideNode{base: 'N'}}
}

// addSeq inserts a sequence into the tree, incrementing per-node
// counters along the path. Stops at the first 'N' or any other
// non-ACGT base. Upstream source: nucleotidetree.cpp:42-55.
func (t *nucleotideTree) addSeq(seq string) {
	cur := t.root
	for i := 0; i < len(seq); i++ {
		if seq[i] == 'N' {
			break
		}
		base := seq[i] & 0x07
		if cur.children[base] == nil {
			cur.children[base] = &nucleotideNode{base: seq[i]}
		}
		cur.children[base].count++
		cur = cur.children[base]
	}
}

// getDominantPath returns the longest path on which every step has a
// "dominant" child (>=95% of the local subtotal) and where the parent
// node has at least 50 total observations. reachedLeaf is set to false
// if we abort early because of no-dominant-child; otherwise true.
// Upstream source: nucleotidetree.cpp:57-87.
func (t *nucleotideTree) getDominantPath() (path string, reachedLeaf bool) {
	const ratioThreshold = 0.95
	const numThreshold = 50

	reachedLeaf = true
	var sb []byte
	cur := t.root
	for {
		total := 0
		for i := 0; i < 8; i++ {
			if cur.children[i] != nil {
				total += cur.children[i].count
			}
		}
		if total < numThreshold {
			break
		}
		hasDominant := false
		for i := 0; i < 8; i++ {
			if cur.children[i] == nil {
				continue
			}
			if float64(cur.children[i].count)/float64(total) >= ratioThreshold {
				hasDominant = true
				sb = append(sb, cur.children[i].base)
				cur = cur.children[i]
				break
			}
		}
		if !hasDominant {
			reachedLeaf = false
			break
		}
	}
	return string(sb), reachedLeaf
}

// evaluatorSeq2Int packs `keylen` bases starting at `pos` into a
// 2-bit-per-base integer. The base-to-code mapping is the upstream
// quirky one: A→0, T→1, C→2, G→3 (evaluator.cpp:564-613). Returns -1
// on any non-ACGT base.
//
// The lastVal optimisation rolls forward by one base: a previous
// 10-mer's value shifted left by 2 and OR'd with the new base's code
// gives the next 10-mer's value, with the high bits masked off.
//
// Upstream source: evaluator.cpp:564-613.
func evaluatorSeq2Int(seq string, pos, keylen, lastVal int) int {
	if lastVal >= 0 {
		mask := (1 << uint(keylen*2)) - 1
		key := (lastVal << 2) & mask
		base := seq[pos+keylen-1]
		switch base {
		case 'A':
			key += 0
		case 'T':
			key += 1
		case 'C':
			key += 2
		case 'G':
			key += 3
		default:
			return -1
		}
		return key
	}
	key := 0
	for i := pos; i < keylen+pos; i++ {
		key <<= 2
		base := seq[i]
		switch base {
		case 'A':
			key += 0
		case 'T':
			key += 1
		case 'C':
			key += 2
		case 'G':
			key += 3
		default:
			return -1
		}
	}
	return key
}

// evaluatorInt2Seq is the inverse of evaluatorSeq2Int: it decodes a
// 2-bit-packed integer back into a `seqlen`-character A/T/C/G string.
// Upstream source: evaluator.cpp:548-558.
func evaluatorInt2Seq(val uint32, seqlen int) string {
	bases := [4]byte{'A', 'T', 'C', 'G'}
	out := make([]byte, seqlen)
	for i := 0; i < seqlen; i++ {
		out[i] = 'N'
	}
	done := 0
	for done < seqlen {
		out[seqlen-done-1] = bases[val&0x03]
		val >>= 2
		done++
	}
	return string(out)
}

// reverseString returns s reversed (used in the backward-tree
// concatenation step). Upstream uses util.h:reverse(string).
func reverseString(s string) string {
	r := []byte(s)
	for i, j := 0, len(r)-1; i < j; i, j = i+1, j-1 {
		r[i], r[j] = r[j], r[i]
	}
	return string(r)
}

// checkKnownAdapters scans `reads` against every known adapter,
// allowing 1 mismatch per 16 bases of overlap, and returns the
// best-supported adapter if it meets the upstream confidence gate:
//
//	maxCount > checkedReads/50  OR
//	(maxCount > checkedReads/200 AND mismatches[adapter] < checkedReads).
//
// Returns "" if no adapter passes the gate. Upstream source:
// evaluator.cpp:207-293.
func checkKnownAdaptersSE(reads []*fastq.Record) string {
	known := knownAdapters()
	possibleCounts := make(map[string]int, len(known))
	mismatchesByAdapter := make(map[string]int, len(known))
	for adapter := range known {
		possibleCounts[adapter] = 0
		mismatchesByAdapter[adapter] = 0
	}

	checkedReads := 0
	checkedBases := 0
	curMaxCount := 0
	for _, r := range reads {
		if r == nil {
			continue
		}
		rdata := r.Sequence
		rlen := len(rdata)

		checkedReads++
		checkedBases += rlen
		if checkedReads > checkKnownAdaptersMaxReads || checkedBases > checkKnownAdaptersMaxBases {
			break
		}
		if curMaxCount > checkKnownAdaptersMaxHit {
			break
		}
		for adapter := range known {
			alen := len(adapter)
			if alen >= rlen {
				continue
			}
			// Not a candidate; skip for speedup. Mirrors evaluator.cpp:250-252.
			if curMaxCount > 20 && possibleCounts[adapter] < curMaxCount/10 {
				continue
			}
			adata := adapter
			for pos := 0; pos < rlen-checkKnownAdaptersMatchReq; pos++ {
				cmplen := rlen - pos
				if alen < cmplen {
					cmplen = alen
				}
				allowedMismatch := cmplen / checkKnownAdaptersOneMismatchPer
				mismatch := 0
				matched := true
				for i := 0; i < cmplen; i++ {
					if adata[i] != rdata[i+pos] {
						mismatch++
						if mismatch > allowedMismatch {
							matched = false
							break
						}
					}
				}
				if matched {
					possibleCounts[adapter]++
					if curMaxCount < possibleCounts[adapter] {
						curMaxCount = possibleCounts[adapter]
					}
					mismatchesByAdapter[adapter] += mismatch
					break
				}
			}
		}
	}

	bestAdapter := ""
	maxCount := 0
	// Note: Go map iteration order is randomised, but upstream's
	// std::map<string,int> iterates in lexicographic key order, which
	// matters when two adapters tie. To match upstream's tie-break we
	// iterate the keys in lexicographic order.
	for _, adapter := range sortedKeys(possibleCounts) {
		c := possibleCounts[adapter]
		if c > maxCount {
			bestAdapter = adapter
			maxCount = c
		}
	}
	if maxCount > checkedReads/50 ||
		(maxCount > checkedReads/200 && mismatchesByAdapter[bestAdapter] < checkedReads) {
		return bestAdapter
	}
	return ""
}

// sortedKeys returns the keys of m in lexicographic order. Used to
// make iteration deterministic and match upstream's std::map ordering.
func sortedKeys(m map[string]int) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	// Simple insertion sort: keys count is small (~200), no allocator
	// noise and no stdlib dep beyond what we already use.
	for i := 1; i < len(keys); i++ {
		for j := i; j > 0 && keys[j] < keys[j-1]; j-- {
			keys[j], keys[j-1] = keys[j-1], keys[j]
		}
	}
	return keys
}

// matchKnownAdapter checks whether `seq` exactly matches any known
// adapter as a prefix (zero mismatches). If yes, returns the matching
// adapter; otherwise "". Upstream source: evaluator.cpp:528-546.
func matchKnownAdapter(seq string) string {
	known := knownAdapters()
	for _, adapter := range sortedKeysFromStrMap(known) {
		if len(seq) < len(adapter) {
			continue
		}
		diff := 0
		end := len(adapter)
		if end > len(seq) {
			end = len(seq)
		}
		for i := 0; i < end; i++ {
			if adapter[i] != seq[i] {
				diff++
			}
		}
		if diff == 0 {
			return adapter
		}
	}
	return ""
}

// sortedKeysFromStrMap returns the keys of a map[string]string in
// lexicographic order. Mirrors std::map iteration order. Identical in
// intent to sortedKeys, kept separate so the call sites read naturally.
func sortedKeysFromStrMap(m map[string]string) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	for i := 1; i < len(keys); i++ {
		for j := i; j > 0 && keys[j] < keys[j-1]; j-- {
			keys[j], keys[j-1] = keys[j-1], keys[j]
		}
	}
	return keys
}

// getAdapterWithSeed builds forward + backward NucleotideTrees from
// reads whose 3'-shifted kmer window matches `seed`, derives the
// dominant path in each direction, and returns the concatenated
// adapter. Upstream source: evaluator.cpp:472-526.
func getAdapterWithSeed(seed int, reads []*fastq.Record, keylen int) string {
	const shiftTail = adapterDetectShiftTail
	maxSearchLen := adapterDetectMaxSearchLength

	forwardTree := newNucleotideTree()
	for _, r := range reads {
		if r == nil {
			continue
		}
		seq := r.Sequence
		key := -1
		for pos := 20; pos <= len(seq)-keylen-shiftTail && pos < maxSearchLen; pos++ {
			key = evaluatorSeq2Int(string(seq), pos, keylen, key)
			if key == seed {
				// Forward tree gets the suffix after the seed.
				tailLen := len(seq) - keylen - shiftTail - pos
				if tailLen < 0 {
					tailLen = 0
				}
				forwardTree.addSeq(string(seq[pos+keylen : pos+keylen+tailLen]))
			}
		}
	}
	forwardPath, forwardReachedLeaf := forwardTree.getDominantPath()

	backwardTree := newNucleotideTree()
	for _, r := range reads {
		if r == nil {
			continue
		}
		seq := r.Sequence
		key := -1
		for pos := 20; pos <= len(seq)-keylen-shiftTail && pos < maxSearchLen; pos++ {
			key = evaluatorSeq2Int(string(seq), pos, keylen, key)
			if key == seed {
				// Backward tree gets the prefix before the seed, reversed.
				prefix := string(seq[:pos])
				backwardTree.addSeq(reverseString(prefix))
			}
		}
	}
	backwardPath, backwardReachedLeaf := backwardTree.getDominantPath()

	// Upstream evaluator.cpp:489-507 declares `bool reachedLeaf = true;`
	// once and passes it by reference to BOTH getDominantPath calls,
	// which only ever AND it down. The final value is therefore
	// `true && !forwardFailed && !backwardFailed`. Reviewer-caught
	// regression on PR #122: discarding the forward result let an
	// adapter slip through when only the backward tree succeeded.
	reachedLeaf := forwardReachedLeaf && backwardReachedLeaf

	adapter := reverseString(backwardPath) + evaluatorInt2Seq(uint32(seed), keylen) + forwardPath
	if len(adapter) > adapterDetectMaxLen {
		adapter = adapter[:adapterDetectMaxLen]
	}

	matched := matchKnownAdapter(adapter)
	if matched != "" {
		return matched
	}
	if reachedLeaf {
		return adapter
	}
	return ""
}

// detectAdapterSEUpstream is the SE half of upstream's
// Evaluator::evalAdapterAndReadNum (evaluator.cpp:295-470), shorn of
// the read-count estimation (which we don't need; we already buffer
// the sample in `reads`). It returns the best-effort detected adapter,
// or "" if none can be derived. The caller is expected to feed it the
// raw sampled records; this function applies all the upstream gates
// internally, including the 10000-record minimum.
func detectAdapterSEUpstream(reads []*fastq.Record) string {
	// records = number of non-nil reads we received.
	records := 0
	for _, r := range reads {
		if r != nil {
			records++
		}
	}
	// evaluator.cpp:344 — we need at least 10000 valid records to
	// evaluate.
	if records < adapterDetectMinRecords {
		return ""
	}

	// Stage 1: try known-adapter substring match. evaluator.cpp:353.
	knownAdapter := checkKnownAdaptersSE(reads)
	if len(knownAdapter) > 8 {
		return knownAdapter
	}

	// Stage 2: kmer histogram + nucleotide-tree dominant-path
	// extraction. evaluator.cpp:362-460.
	const shiftTail = adapterDetectShiftTail
	keylen := adapterDetectKeylen
	size := 1 << uint(keylen*2)
	counts := make([]uint32, size)
	for _, r := range reads {
		if r == nil {
			continue
		}
		seq := r.Sequence
		key := -1
		for pos := 20; pos <= len(seq)-keylen-shiftTail; pos++ {
			key = evaluatorSeq2Int(string(seq), pos, keylen, key)
			if key >= 0 {
				counts[key]++
			}
		}
	}

	// AAAAAAAAAA = 0 is explicitly forced to 0 (evaluator.cpp:384).
	counts[0] = 0

	// Find the top-N kmers. evaluator.cpp:387-430. We replicate the
	// upstream insertion-sort exactly so tie ordering matches.
	topnum := adapterDetectTopN
	topkeys := make([]int, topnum)
	var total int64
	for k := 0; k < size; k++ {
		atcg := [4]int{}
		for i := 0; i < keylen; i++ {
			baseOfBit := (k >> uint(i*2)) & 0x03
			atcg[baseOfBit]++
		}
		// Low complexity: any single base accounts for >= keylen-4.
		lowComplexity := false
		for b := 0; b < 4; b++ {
			if atcg[b] >= keylen-4 {
				lowComplexity = true
				break
			}
		}
		if lowComplexity {
			continue
		}
		// Too many GC: C+G >= keylen-2.
		if atcg[2]+atcg[3] >= keylen-2 {
			continue
		}
		// Starts with GGGG: high bits == 0xff (0b11111111 = four G's
		// in the bits 12..19 position). evaluator.cpp:407-409.
		if k>>12 == 0xff {
			continue
		}

		val := counts[k]
		total += int64(val)
		// Upstream's "find slot in topkeys" loop (insertion-sort
		// style). evaluator.cpp:413-430.
		for t := topnum - 1; t >= 0; t-- {
			if val < counts[topkeys[t]] {
				if t < topnum-1 {
					for m := topnum - 1; m > t+1; m-- {
						topkeys[m] = topkeys[m-1]
					}
					topkeys[t+1] = k
				}
				break
			} else if t == 0 {
				for m := topnum - 1; m > t; m-- {
					topkeys[m] = topkeys[m-1]
				}
				topkeys[t] = k
			}
		}
	}

	// Try each top kmer as a seed in order. evaluator.cpp:432-460.
	foldThreshold := adapterDetectFoldThreshold
	for t := 0; t < topnum; t++ {
		key := topkeys[t]
		if key == 0 {
			continue
		}
		count := int64(counts[key])
		// Equivalent to upstream's: count < 10 || count*size < total * FOLD_THRESHOLD.
		if count < 10 || count*int64(size) < total*int64(foldThreshold) {
			break
		}
		seq := evaluatorInt2Seq(uint32(key), keylen)
		// Skip low-complexity-by-runs (less than 3 transitions).
		diff := 0
		for s := 0; s < len(seq)-1; s++ {
			if seq[s] != seq[s+1] {
				diff++
			}
		}
		if diff < 3 {
			continue
		}
		adapter := getAdapterWithSeed(key, reads, keylen)
		if adapter != "" {
			return adapter
		}
	}
	return ""
}
