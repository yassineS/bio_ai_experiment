package samtools

// Port of htslib's ks_introsort_##name applied to fragPtr arrays
// (`ks_introsort_rseq` in phase.c). Faithful to ksort.h's KSORT_INIT_
// expansion so that the EV-line emit order matches upstream byte for
// byte on inputs that exercise the natural sort path.
//
// Only the introsort + insertsort + combsort + isort-stack pieces are
// ported; the heapsort/ksmall/shuffle paths are unused by phase.c on
// the natural code path.

// ksIsortStack mirrors the ks_isort_stack_t type in ksort.h.
type ksIsortStack struct {
	left, right int
	depth       int
}

// ksortRseq is a faithful port of ks_introsort_rseq(n, a). It sorts
// a[0..n) in place by frag.vpos (rseq_lt).
//
// The algorithm is the standard introsort: median-of-3 quicksort
// partition; once the working sub-array drops to ≤16 elements stop
// recursing (those tail-ranges get cleaned up by a final pass of
// `__ks_insertsort` over the whole array); if the recursion depth
// would exceed 2*ceil(log2 n) the current sub-array is delegated to
// combsort instead.
func ksortRseq(a []fragPtr) {
	n := len(a)
	if n < 1 {
		return
	}
	if n == 2 {
		if fragRseqLt(a[1], a[0]) {
			a[0], a[1] = a[1], a[0]
		}
		return
	}
	// d = 2; while (1 << d) < n: d++  → d = ceil(log2(n)). Match the
	// upstream loop's post-increment so d matches phase.c semantics.
	d := 2
	for (1 << d) < n {
		d++
	}
	d <<= 1
	stack := make([]ksIsortStack, 0, 64)
	s, t := 0, n-1
	for {
		if s < t {
			d--
			if d == 0 {
				ksCombsortRseq(a[s : t+1])
				t = s
				continue
			}
			// Median-of-three pivot selection.
			i, j := s, t
			k := i + ((j - i) >> 1) + 1
			if fragRseqLt(a[k], a[i]) {
				if fragRseqLt(a[k], a[j]) {
					k = j
				}
			} else if fragRseqLt(a[j], a[i]) {
				k = i
			} else {
				k = j
			}
			rp := a[k]
			if k != t {
				a[k], a[t] = a[t], a[k]
			}
			for {
				for {
					i++
					if !fragRseqLt(a[i], rp) {
						break
					}
				}
				for {
					j--
					if !(i <= j && fragRseqLt(rp, a[j])) {
						break
					}
				}
				if j <= i {
					break
				}
				a[i], a[j] = a[j], a[i]
			}
			a[i], a[t] = a[t], a[i]
			if i-s > t-i {
				if i-s > 16 {
					stack = append(stack, ksIsortStack{left: s, right: i - 1, depth: d})
				}
				if t-i > 16 {
					s = i + 1
				} else {
					s = t
				}
			} else {
				if t-i > 16 {
					stack = append(stack, ksIsortStack{left: i + 1, right: t, depth: d})
				}
				if i-s > 16 {
					t = i - 1
				} else {
					t = s
				}
			}
		} else {
			if len(stack) == 0 {
				ksInsertsortRseq(a)
				return
			}
			top := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			s = top.left
			t = top.right
			d = top.depth
		}
	}
}

// ksInsertsortRseq is the final cleanup pass — matches __ks_insertsort##name.
// Walks the array left-to-right, bubbling each element back to its
// place via swaps. Stable for equal keys (the algorithm only swaps on
// strict <).
func ksInsertsortRseq(a []fragPtr) {
	if len(a) < 2 {
		return
	}
	for i := 1; i < len(a); i++ {
		// upstream's __ks_insertsort##name walks a `t` (last) pointer
		// from a+1 to a+n-1, and for each position swaps backwards
		// while __sort_lt(*(j+1), *j). The two-pointer formulation
		// below matches that semantics.
		for j := i; j > 0 && fragRseqLt(a[j], a[j-1]); j-- {
			a[j], a[j-1] = a[j-1], a[j]
		}
	}
}

// ksCombsortRseq matches __ks_combsort##name — comb sort with a 1.3
// shrink factor, falling back to insertsort once gap == 1.
func ksCombsortRseq(a []fragPtr) {
	n := len(a)
	if n < 2 {
		return
	}
	const shrinkF = 1.2473309
	gap := n
	for swap := true; gap > 1 || swap; {
		// gap = (int)(gap / shrinkF + 0.5)
		gap = int(float64(gap)/shrinkF + 0.5)
		if gap == 9 || gap == 10 {
			gap = 11 // combsort 11 heuristic
		}
		if gap < 1 {
			gap = 1
		}
		swap = false
		for i := 0; i < n-gap; i++ {
			j := i + gap
			if fragRseqLt(a[j], a[i]) {
				a[i], a[j] = a[j], a[i]
				swap = true
			}
		}
	}
	if gap != 1 {
		ksInsertsortRseq(a)
	}
}
