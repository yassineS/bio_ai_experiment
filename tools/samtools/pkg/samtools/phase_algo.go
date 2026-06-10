package samtools

// Faithful Go ports of the algorithmic functions in
// reference_code/samtools/phase.c. Each function is annotated with its
// upstream line range. The bit-level semantics (uint32 masks, signed
// counts, etc.) are preserved exactly so byte-parity vs. upstream is
// achievable on the supported fixtures.

// count1 — phase.c:94. Tally per-haplotype combinations for one
// fragment's window of l hets (the l-window slides over the fragment
// in count_all). `seq` is a window of l values from frag.seq, with
// 0 = ambiguous and 1/2 = allele code. cnt is the (1<<l)-element
// histogram for that het's column.
func count1(l int, seq []int8, cnt []int32) {
	if seq[l-1] == 0 {
		return // last base ambiguous → contribute nothing
	}
	nAmbi := 0
	for i := 0; i < l; i++ {
		if seq[i] == 0 {
			nAmbi++
		}
	}
	if l-nAmbi <= 1 {
		return // only one SNP
	}
	// For each assignment of the ambiguous bases (2^nAmbi enumerations),
	// build the local haplotype index z and increment cnt[z].
	for x := uint32(0); x < (1 << uint(nAmbi)); x++ {
		var i, j int
		var z uint32
		for i = 0; i < l; i++ {
			var c uint32
			if seq[i] != 0 {
				c = uint32(seq[i]) - 1
			} else {
				c = (x >> uint(j)) & 1
				j++
			}
			z = z<<1 | c
		}
		cnt[z]++
	}
}

// countAll — phase.c:116. Builds the per-het histogram tables used by
// dynaprog (Viterbi). For each fragment with vlen >= 2, slides an
// l-wide window across its hets and calls count1 for the right-edge
// het's bucket. Returns a vpos-by-(1<<l) matrix.
func countAll(l int, vpos int, hash *fragKhash) [][]int32 {
	cntSz := 1 << uint(l)
	cnt := make([][]int32, vpos)
	for i := 0; i < vpos; i++ {
		cnt[i] = make([]int32, cntSz)
	}
	seq := make([]int8, l)
	for k := uint32(0); k < hash.end(); k++ {
		if !hash.exist(k) {
			continue
		}
		f := &hash.vals[k]
		if int(f.vpos) >= vpos || f.single != 0 {
			continue
		}
		if f.vlen == 1 {
			// "Such reads should be flagged as deleted previously if
			// everything is right." Mark as single now.
			f.single = 1
			continue
		}
		for j := 1; j < int(f.vlen); j++ {
			for i := 0; i < l; i++ {
				// seq[i] = j < l-1-i? 0 : f->seq[j - (l-1-i)]
				if j < l-1-i {
					seq[i] = 0
				} else {
					seq[i] = f.seq[j-(l-1-i)]
				}
			}
			count1(l, seq, cnt[int(f.vpos)+j])
		}
	}
	return cnt
}

// dynaprog — phase.c:163. Viterbi forward+backtrack over k-bit local
// haplotype states. Returns the per-het hap0 assignment (0 or 1) as
// `path[0..vpos)`.
func dynaprog(l int, vpos int, w [][]int32) []int8 {
	z := uint32(1) << uint(l-1)
	mask := (uint32(1) << uint(l)) - 1
	f0 := make([]int32, z)
	f1 := make([]int32, z)
	prev, curr := f0, f1
	b := make([][]int8, vpos)
	for i := 0; i < vpos; i++ {
		wi := w[i]
		bi := make([]int8, z)
		b[i] = bi
		// For each current state x (smaller haplotype), compute the two
		// predecessor candidates y0 = x>>1 and y1 = xc>>1 (xc = complement
		// of x within the l-bit window). Choose the higher.
		for x := uint32(0); x < z; x++ {
			xc := (^x) & mask
			y0 := x >> 1
			y1 := xc >> 1
			c0 := prev[y0] + wi[x] + wi[xc]
			c1 := prev[y1] + wi[x] + wi[xc]
			if c0 > c1 {
				bi[x] = 0
				curr[x] = c0
			} else {
				bi[x] = 1
				curr[x] = c1
			}
		}
		prev, curr = curr, prev
	}
	// Backtrack.
	var maxV int32 = 0
	var maxX uint32 = 0
	for x := uint32(0); x < z; x++ {
		if prev[x] > maxV {
			maxV = prev[x]
			maxX = x
		}
	}
	h := make([]int8, vpos)
	which := 0
	x := maxX
	for i := vpos - 1; i >= 0; i-- {
		if which != 0 {
			h[i] = int8((^x) & 1)
		} else {
			h[i] = int8(x & 1)
		}
		if b[i][x] != 0 {
			which = 1 - which
			x = ((^x) & mask) >> 1
		} else {
			x = x >> 1
		}
	}
	return h
}

// fragphase — phase.c:211. For each fragment, compute its haplotype
// assignment vs. the Viterbi path and, when `flip` is set, optionally
// flip one half of the fragment to repair chimeras. Returns the
// per-het pcnt[] vote tally used to drive M-line emission.
//
// The 64-bit pcnt[i] is packed as four 16-bit counters:
//
//	bits  0..15  supports0  (hap0 matches)
//	bits 16..31  errors0    (hap0 mismatches)
//	bits 32..47  supports1  (hap1 matches)
//	bits 48..63  errors1    (hap1 mismatches)
//
// The M-line uses each 16-bit field directly.
func fragphase(vpos int, path []int8, hash *fragKhash, flip bool) []uint64 {
	pcnt := make([]uint64, vpos)
	var left, right []uint32
	maxLen := 0
	for k := uint32(0); k < hash.end(); k++ {
		if !hash.exist(k) {
			continue
		}
		f := &hash.vals[k]
		if int(f.vpos) >= vpos {
			continue
		}
		var c [2]int
		for i := 0; i < int(f.vlen); i++ {
			if f.seq[i] == 0 {
				continue
			}
			// c[0] increments when seq matches path+1 (i.e. is on hap0
			// in the read's own frame). c[1] is the other side.
			if f.seq[i] == path[int(f.vpos)+i]+1 {
				c[0]++
			} else {
				c[1]++
			}
		}
		if c[0] > c[1] {
			f.phase = 0
		} else {
			f.phase = 1
		}
		// in = count on f.phase side (the majority), out = opposite.
		if f.phase == 0 {
			f.in = uint16(c[0])
			f.out = uint16(c[1])
		} else {
			f.in = uint16(c[1])
			f.out = uint16(c[0])
		}
		if f.in == f.out {
			f.phased = 0
		} else {
			f.phased = 1
		}
		if f.in != 0 && f.out != 0 && f.out < 3 && f.in <= f.out+1 {
			f.ambig = 1
		} else {
			f.ambig = 0
		}
		f.flip = 0
		if flip && c[0] >= 3 && c[1] >= 3 {
			// Grow left/right scratch arrays as needed.
			if int(f.vlen) > maxLen {
				maxLen = int(f.vlen)
				// kroundup-style growth — easiest in Go is just track
				// max len.
				left = make([]uint32, maxLen)
				right = make([]uint32, maxLen)
			} else {
				left = left[:int(f.vlen)]
				right = right[:int(f.vlen)]
			}
			var sum [2]int
			// left-counts: cumulative sums of (match, mismatch) of the
			// read's allele to the Viterbi path, walking left→right.
			sum[0], sum[1] = 0, 0
			for i := 0; i < int(f.vlen); i++ {
				if f.seq[i] != 0 {
					var c2 int
					if f.phase != 0 {
						c2 = 2 - int(f.seq[i])
					} else {
						c2 = int(f.seq[i]) - 1
					}
					if c2 == int(path[int(f.vpos)+i]) {
						sum[0]++
					} else {
						sum[1]++
					}
				}
				left[i] = uint32(sum[1])<<16 | uint32(sum[0])
			}
			// right-counts.
			sum[0], sum[1] = 0, 0
			for i := int(f.vlen) - 1; i >= 0; i-- {
				if f.seq[i] != 0 {
					var c2 int
					if f.phase != 0 {
						c2 = 2 - int(f.seq[i])
					} else {
						c2 = int(f.seq[i]) - 1
					}
					if c2 == int(path[int(f.vpos)+i]) {
						sum[0]++
					} else {
						sum[1]++
					}
				}
				right[i] = uint32(sum[1])<<16 | uint32(sum[0])
			}
			// find best flip point
			m := 0
			mi := -1
			md := -1
			for i := 0; i < int(f.vlen)-1; i++ {
				var a [2]int
				a[0] = int(left[i]&0xffff) + int(right[i+1]>>16&0xffff) - int(right[i+1]&0xffff)*flipPenalty
				a[1] = int(left[i]>>16&0xffff) + int(right[i+1]&0xffff) - int(right[i+1]>>16&0xffff)*flipPenalty
				if a[0] > a[1] {
					if a[0] > m {
						m = a[0]
						md = 0
						mi = i
					}
				} else {
					if a[1] > m {
						m = a[1]
						md = 1
						mi = i
					}
				}
			}
			if m-c[0] >= flipThreshold && m-c[1] >= flipThreshold {
				f.flip = 1
				if md == 0 {
					// flip the tail
					for i := mi + 1; i < int(f.vlen); i++ {
						if f.seq[i] == 1 {
							f.seq[i] = 2
						} else if f.seq[i] == 2 {
							f.seq[i] = 1
						}
					}
				} else {
					// flip the head
					for i := 0; i <= mi; i++ {
						if f.seq[i] == 1 {
							f.seq[i] = 2
						} else if f.seq[i] == 2 {
							f.seq[i] = 1
						}
					}
				}
			}
		}
		// update pcnt[]
		if f.single == 0 {
			for i := 0; i < int(f.vlen); i++ {
				if f.seq[i] == 0 {
					continue
				}
				var c2 int
				if f.phase != 0 {
					c2 = 2 - int(f.seq[i])
				} else {
					c2 = int(f.seq[i]) - 1
				}
				idx := int(f.vpos) + i
				if c2 == int(path[idx]) {
					if f.phase == 0 {
						pcnt[idx] += 1
					} else {
						pcnt[idx] += uint64(1) << 32
					}
				} else {
					if f.phase == 0 {
						pcnt[idx] += uint64(1) << 16
					} else {
						pcnt[idx] += uint64(1) << 48
					}
				}
			}
		}
	}
	return pcnt
}

// genmask — phase.c:302. Walks pcnt[] computing a sliding score for
// "this region looks like a chimeric block" and returns the list of
// (start, end) intervals that exceeded MASK_THRES.
func genmask(vpos int, pcnt []uint64, nOut *int) []uint64 {
	const maskThres = 3
	maxV := 0
	maxI := -1
	m := 0
	n := 0
	beg := 0
	score := 0
	var list []uint64
	for i := 0; i < vpos; i++ {
		x := pcnt[i]
		var c [4]int
		c[0] = int(x & 0xffff)
		c[1] = int(x >> 16 & 0xffff)
		c[2] = int(x >> 32 & 0xffff)
		c[3] = int(x >> 48 & 0xffff)
		pre := score
		var s int
		if c[1]+c[3] == 0 {
			s = -(c[0] + c[2])
		} else {
			s = c[1] + c[3] - 1
		}
		if c[3] > c[2] {
			s += c[3] - c[2]
		}
		if c[1] > c[0] {
			s += c[1] - c[0]
		}
		score += s
		if score < 0 {
			score = 0
		}
		if pre == 0 && score > 0 {
			beg = i
		}
		if (i == vpos-1 || score == 0) && maxV >= maskThres {
			if n == m {
				if m == 0 {
					m = 4
				} else {
					m <<= 1
				}
				newList := make([]uint64, m)
				copy(newList, list)
				list = newList
			}
			list[n] = uint64(beg)<<32 | uint64(maxI)
			n++
			i = maxI // reset i to max_i (upstream re-walks from there)
			score = 0
		} else if score > maxV {
			maxV = score
			maxI = i
		}
		if score == 0 {
			maxV = 0
		}
	}
	*nOut = n
	return list
}

// cleanSeqs — phase.c:332. Trims leading/trailing ambiguous bases from
// each frag's seq and deletes empty entries from the hash.
// Returns true iff any frag has vpos >= vpos (i.e. extends beyond the
// current block — used to compute min_pos for dump_aln).
func cleanSeqs(vpos int, hash *fragKhash) bool {
	ret := false
	for k := uint32(0); k < hash.end(); k++ {
		if !hash.exist(k) {
			continue
		}
		f := &hash.vals[k]
		if int(f.vpos) >= vpos {
			ret = true
			continue
		}
		beg := 0
		for ; beg < int(f.vlen); beg++ {
			if f.seq[beg] != 0 {
				break
			}
		}
		end := int(f.vlen) - 1
		for ; end >= 0; end-- {
			if f.seq[end] != 0 {
				break
			}
		}
		end++ // exclusive
		if end-beg <= 0 {
			hash.del(k)
		} else {
			if beg != 0 {
				copy(f.seq[:end-beg], f.seq[beg:end])
				// zero the tail to mimic memmove-ish behaviour
				for i := end - beg; i < int(f.vlen); i++ {
					f.seq[i] = 0
				}
			}
			f.vpos += int32(beg)
			f.vlen = uint16(end - beg)
			if f.vlen == 1 {
				f.single = 1
			} else {
				f.single = 0
			}
		}
	}
	return ret
}

// updateVpos — phase.c:488. Slides every frag's vpos down by `vpos`,
// deleting those that fall below zero.
func updateVpos(vpos int, hash *fragKhash) {
	for k := uint32(0); k < hash.end(); k++ {
		if !hash.exist(k) {
			continue
		}
		f := &hash.vals[k]
		if int(f.vpos) < vpos {
			hash.del(k)
		} else {
			f.vpos -= int32(vpos)
		}
	}
}

// gl2cns — phase.c:561. Reduces a 4×4 genotype-likelihood matrix to a
// consensus call: returns 0 for a homozygous best genotype, or
// 1<<18 | bigCode<<16 | smallCode | LOD<<2 for a heterozygous best.
//
// q is interpreted as a 16-element float32 matrix indexed by (i<<2|j)
// for i,j in 0..3.
func gl2cns(q []float32) uint32 {
	var minV, min2 float32 = 1e30, 1e30
	minIJ := -1
	for i := 0; i < 4; i++ {
		for j := i; j < 4; j++ {
			v := q[i<<2|j]
			if v < minV {
				minIJ = i<<2 | j
				min2 = minV
				minV = v
			} else if v < min2 {
				min2 = v
			}
		}
	}
	if minIJ < 0 {
		return 0
	}
	if (minIJ>>2)&3 == (minIJ & 3) {
		return 0
	}
	// LOD = (int)(min2 - min + .499)
	lod := int(min2 - minV + 0.499)
	return uint32(1)<<18 | uint32((minIJ>>2)&3)<<16 | uint32(minIJ&3) | uint32(lod)<<2
}

// kSize counts occupied buckets — kh_size(hash).
func kSize(hash *fragKhash) int {
	return int(hash.size)
}
