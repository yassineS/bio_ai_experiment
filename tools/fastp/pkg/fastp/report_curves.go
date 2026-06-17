// Per-read base-resolved curve and k-mer accumulation for the JSON report.
//
// Upstream fastp's Stats class (reference_code/fastp/src/stats.cpp) collects,
// for every read group (read1/read2, before/after filtering), a set of
// per-cycle, per-base counters that drive the JSON report's
// quality_curves / content_curves / kmer_count / q40_bases sub-fields. The
// Go port previously tracked only an aggregate mean-quality curve and an
// A/C/G/T/N base count for the *before* stream, which was enough for the
// summary counters but not for byte/structurally-exact per-read sub-fields.
//
// readCurves reproduces upstream's accumulation exactly:
//
//   - Per cycle c and base b in {A,T,C,G,N}: baseContent[b][c] (count) and
//     baseQualSum[b][c] (sum of phred scores). These feed both content_curves
//     (count/totalBase) and quality_curves (qualSum/baseContent, or the mean
//     when the base never occurs at that cycle), matching stats.cpp:182-201.
//   - totalBase[c] and totalQual[c]: the per-cycle totals used for the "mean"
//     quality curve and the content-curve denominators (stats.cpp:174-179).
//   - q40 base count: the number of bases with phred >= 40, matching the
//     mBaseQualHistogram[c>=40+33] sum in stats.cpp:169-171.
//   - kmer[1024]: 5-mer histogram over A/T/C/G runs (N resets the window),
//     base2val A=0,T=1,C=2,G=3, index ((kmer<<2)&0x3FC)|val, matching
//     statRead's kmer loop (stats.cpp:243-309).
//
// The base ordering matters: upstream's content_curves emit A,T,C,G,N,GC and
// quality_curves emit A,T,C,G,mean, so readCurves stores bases in the index
// order [A,T,C,G,N] and the emitters reproduce upstream's key order.

package fastp

import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"

// curveBaseOrder is upstream's base ordering for the per-base curves
// (stats.cpp qualNames/contentNames): A, T, C, G, then N.
var curveBaseOrder = [5]byte{'A', 'T', 'C', 'G', 'N'}

// curveBaseIndex maps a base byte to its slot in readCurves, using upstream's
// A,T,C,G,N ordering (note this differs from baseIndex, which is A,C,G,T,N).
// Anything not A/C/G/T is treated as N.
func curveBaseIndex(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'T', 't':
		return 1
	case 'C', 'c':
		return 2
	case 'G', 'g':
		return 3
	default:
		return 4
	}
}

// kmerBase2Val maps a base to upstream's base2val (A=0,T=1,C=2,G=3); any other
// base (including N) returns -1 to break the 5-mer window (stats.cpp:334-347).
func kmerBase2Val(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'T', 't':
		return 1
	case 'C', 'c':
		return 2
	case 'G', 'g':
		return 3
	default:
		return -1
	}
}

// readCurves accumulates the per-cycle, per-base counters and the 5-mer
// histogram for one read group, reproducing upstream Stats.
type readCurves struct {
	// baseContent[b][c] is the count of base b at cycle c; baseQualSum[b][c]
	// is the sum of phred scores of base b at cycle c. b indexes curveBaseOrder.
	baseContent [5][]int64
	baseQualSum [5][]int64
	// totalBase[c] / totalQual[c] are the per-cycle totals across all bases.
	totalBase []int64
	totalQual []int64
	// q40 is the number of bases with phred >= 40.
	q40 int64
	// kmer is the 5-mer histogram, indexed exactly as upstream's mKmer.
	kmer [1024]int64
	// reads is the number of reads observed; bases is the total base count.
	reads int64
	bases int64
}

// grow ensures the per-cycle slices can hold at least n cycles.
func (rc *readCurves) grow(n int) {
	if len(rc.totalBase) >= n {
		return
	}
	nb := make([]int64, n)
	copy(nb, rc.totalBase)
	rc.totalBase = nb
	nq := make([]int64, n)
	copy(nq, rc.totalQual)
	rc.totalQual = nq
	for b := 0; b < 5; b++ {
		c1 := make([]int64, n)
		copy(c1, rc.baseContent[b])
		rc.baseContent[b] = c1
		c2 := make([]int64, n)
		copy(c2, rc.baseQualSum[b])
		rc.baseQualSum[b] = c2
	}
}

// stat accumulates one record into the curves, mirroring Stats::statRead
// (stats.cpp:232-332): per-cycle base/quality tallies, the q40 histogram bump,
// and the N-resetting 5-mer window.
func (rc *readCurves) stat(record *fastq.Record, offset int) {
	if record == nil {
		return
	}
	n := len(record.Sequence)
	rc.grow(n)
	rc.reads++
	rc.bases += int64(n)

	kmer := 0
	needFullCompute := true
	for i := 0; i < n; i++ {
		base := record.Sequence[i]
		q := int(record.Quality[i]) - offset
		bi := curveBaseIndex(base)
		rc.baseContent[bi][i]++
		rc.baseQualSum[bi][i] += int64(q)
		rc.totalBase[i]++
		rc.totalQual[i] += int64(q)
		if q >= 40 {
			rc.q40++
		}
		// 5-mer accumulation (A/T/C/G only; N resets the window).
		if base == 'N' || base == 'n' {
			needFullCompute = true
			continue
		}
		if i < 4 {
			continue
		}
		if !needFullCompute {
			val := kmerBase2Val(base)
			if val < 0 {
				needFullCompute = true
				continue
			}
			kmer = ((kmer << 2) & 0x3FC) | val
			rc.kmer[kmer]++
		} else {
			valid := true
			kmer = 0
			for k := 0; k < 5; k++ {
				val := kmerBase2Val(record.Sequence[i-4+k])
				if val < 0 {
					valid = false
					break
				}
				kmer = ((kmer << 2) & 0x3FC) | val
			}
			if !valid {
				needFullCompute = true
				continue
			}
			rc.kmer[kmer]++
			needFullCompute = false
		}
	}
}

// cycles returns the number of cycles (the longest read length observed),
// matching upstream's mCycles = max read length.
func (rc *readCurves) cycles() int {
	return len(rc.totalBase)
}

// qualityCurves returns the per-base mean-quality curves keyed by upstream's
// quality_curves order (A,T,C,G,mean). For a base that never occurs at a cycle
// upstream substitutes the overall mean quality there (stats.cpp:193-194).
func (rc *readCurves) qualityCurves() map[string][]float64 {
	n := rc.cycles()
	if n == 0 {
		return nil
	}
	mean := make([]float64, n)
	for c := 0; c < n; c++ {
		if rc.totalBase[c] > 0 {
			mean[c] = float64(rc.totalQual[c]) / float64(rc.totalBase[c])
		}
	}
	out := make(map[string][]float64, 5)
	for bi := 0; bi < 4; bi++ { // A,T,C,G
		curve := make([]float64, n)
		for c := 0; c < n; c++ {
			if rc.baseContent[bi][c] == 0 {
				curve[c] = mean[c]
			} else {
				curve[c] = float64(rc.baseQualSum[bi][c]) / float64(rc.baseContent[bi][c])
			}
		}
		out[string(curveBaseOrder[bi])] = curve
	}
	out["mean"] = mean
	return out
}

// contentCurves returns the per-base content fractions keyed by upstream's
// content_curves order (A,T,C,G,N,GC), each value content/totalBase per cycle.
func (rc *readCurves) contentCurves() map[string][]float64 {
	n := rc.cycles()
	if n == 0 {
		return nil
	}
	out := make(map[string][]float64, 6)
	for bi := 0; bi < 5; bi++ { // A,T,C,G,N
		curve := make([]float64, n)
		for c := 0; c < n; c++ {
			if rc.totalBase[c] > 0 {
				curve[c] = float64(rc.baseContent[bi][c]) / float64(rc.totalBase[c])
			}
		}
		out[string(curveBaseOrder[bi])] = curve
	}
	gc := make([]float64, n)
	for c := 0; c < n; c++ {
		if rc.totalBase[c] > 0 {
			// C is index 2, G is index 3 in curveBaseOrder.
			gc[c] = float64(rc.baseContent[2][c]+rc.baseContent[3][c]) / float64(rc.totalBase[c])
		}
	}
	out["GC"] = gc
	return out
}

// kmerCounts returns the 5-mer histogram keyed by the kmer string in upstream's
// emission order (i in 0..63 outer, j in 0..15 inner; key = kmer3(i)+kmer2(j)).
func (rc *readCurves) kmerCounts() map[string]int64 {
	if rc.bases == 0 {
		return nil
	}
	out := make(map[string]int64, 1024)
	for i := 0; i < 64; i++ {
		first := kmer3(i)
		for j := 0; j < 16; j++ {
			target := (i << 4) + j
			out[first+kmer2(j)] = rc.kmer[target]
		}
	}
	return out
}

// kmer3 / kmer2 reproduce upstream's kmer-index-to-string mappings
// (stats.cpp:723-738), bases ordered {A,T,C,G}.
func kmer3(val int) string {
	bases := [4]byte{'A', 'T', 'C', 'G'}
	return string([]byte{
		bases[(val&0x30)>>4],
		bases[(val&0x0C)>>2],
		bases[val&0x03],
	})
}

func kmer2(val int) string {
	bases := [4]byte{'A', 'T', 'C', 'G'}
	return string([]byte{
		bases[(val&0x0C)>>2],
		bases[val&0x03],
	})
}
