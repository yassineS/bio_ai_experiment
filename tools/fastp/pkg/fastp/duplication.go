// Duplication evaluation. Approximates the upstream fastp algorithm with a
// fixed-size hash table sized by --dup_calc_accuracy. Each cell stores a
// uint16 occurrence counter for whichever k-mer hash maps to that cell
// first; subsequent reads whose hash collides with an existing non-empty
// cell are recorded as duplicates regardless of true sequence identity.
// This trades a small accuracy hit for predictable memory usage and
// streaming-friendly behaviour (matching upstream).
//
// Output statistics:
//
//   - Rate: fraction of input reads that map to an already-seen hash cell.
//   - Hist: histogram of duplicate-group sizes (count -> number of reads
//     observed at that occurrence count). hist[1] is the number of reads
//     that were seen exactly once (i.e. unique); hist[2] is the number of
//     reads that map to a hash seen twice in total, etc.
//
// Use case beyond reporting: when --dedup is set the same tracker is used
// to drop the second and later occurrences from the output FASTQ stream.

package fastp

import "hash/fnv"

// dupAccuracyMin and dupAccuracyMax are the inclusive bounds upstream
// fastp enforces on --dup_calc_accuracy. The default is 3.
const (
	dupAccuracyMin     = 1
	dupAccuracyMax     = 6
	dupAccuracyDefault = 3
	dupKeyLen          = 16
)

// DupTracker accumulates per-read duplication evidence using a fixed-size
// hash-collision table. It is intentionally cheap and approximate.
//
// The table has 1<<(17+accuracy) entries each holding a uint16 counter, so
// memory usage is roughly 2 * 2^(17+accuracy) bytes (≈ 1MB at accuracy=3).
//
// A DupTracker is not safe for concurrent use; callers running parallel
// workers should serialize access (e.g. compute hashes in parallel but
// update the table from a single goroutine).
type DupTracker struct {
	accuracy int
	mask     uint64
	cells    []uint16
	total    int64
	dup      int64
}

// NewDupTracker constructs a DupTracker with the given accuracy
// (clamped to [dupAccuracyMin, dupAccuracyMax]). Larger accuracy uses more
// memory but reduces spurious collisions on diverse libraries.
func NewDupTracker(accuracy int) *DupTracker {
	if accuracy < dupAccuracyMin {
		accuracy = dupAccuracyMin
	}
	if accuracy > dupAccuracyMax {
		accuracy = dupAccuracyMax
	}
	size := uint64(1) << uint(17+accuracy)
	return &DupTracker{
		accuracy: accuracy,
		mask:     size - 1,
		cells:    make([]uint16, size),
	}
}

// Observe records one read sequence and returns true if this read is a
// duplicate (i.e. the cell for its hash was already non-zero).
//
// Sequences shorter than dupKeyLen are hashed as a whole. Empty
// sequences are silently ignored and reported as non-duplicate.
func (t *DupTracker) Observe(seq []byte) bool {
	if t == nil || len(seq) == 0 {
		return false
	}
	key := seq
	if len(key) > dupKeyLen {
		key = key[:dupKeyLen]
	}
	h := fnv.New64a()
	_, _ = h.Write(key)
	idx := h.Sum64() & t.mask
	t.total++
	if t.cells[idx] == 0 {
		t.cells[idx] = 1
		return false
	}
	// Avoid uint16 overflow; cap at 65535 (sufficient for histogram).
	if t.cells[idx] < ^uint16(0) {
		t.cells[idx]++
	}
	t.dup++
	return true
}

// Rate returns the duplication rate as a fraction in [0, 1]. It is the
// number of observations that hit an already-occupied cell divided by the
// total number of observations.
func (t *DupTracker) Rate() float64 {
	if t == nil || t.total == 0 {
		return 0
	}
	return float64(t.dup) / float64(t.total)
}

// Histogram returns a map from occurrence-count -> number of reads
// observed at that count. hist[1] counts reads whose hash cell was
// unique; hist[k] for k>1 counts reads whose hash cell ended up with k
// total observations (so each such cell contributes k reads to hist[k]).
//
// The histogram is suitable for plotting "duplication levels" similar to
// upstream fastp's report.
func (t *DupTracker) Histogram() map[int]int64 {
	if t == nil {
		return nil
	}
	hist := make(map[int]int64)
	for _, c := range t.cells {
		if c == 0 {
			continue
		}
		k := int(c)
		hist[k] += int64(k)
	}
	return hist
}

// Total returns the total number of observations recorded.
func (t *DupTracker) Total() int64 {
	if t == nil {
		return 0
	}
	return t.total
}

// Accuracy returns the configured accuracy bucket (1..6).
func (t *DupTracker) Accuracy() int {
	if t == nil {
		return 0
	}
	return t.accuracy
}
