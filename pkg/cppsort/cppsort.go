// Package cppsort provides a byte-faithful port of libstdc++'s std::sort
// (introsort) so Go tools can reproduce, exactly, the order C++ tools emit for
// elements that compare equal under a partial comparator.
//
// Why this exists: several bedtools subcommands (e.g. `sort`, `cluster`) call
// std::sort with a comparator that only compares part of the record (chromStart
// only). std::sort is NOT stable, so records that tie on that key come out in
// an order that is an artifact of the introsort algorithm — neither input order
// nor any simple secondary key (empirically the tie order is a mix of
// ascending and descending by the untested fields). Go's sort.Sort /
// sort.Slice use a different algorithm (pattern-defeating quicksort), so they
// produce a different — though equally "valid" — tie order, which diverges
// byte-for-byte from upstream.
//
// This package mirrors libstdc++'s <bits/stl_algo.h> introsort exactly:
//   - introsort_loop with depth limit 2*floor(log2(n)), median-of-3 pivot,
//     Hoare-style unguarded partition, heapsort fallback at depth 0;
//   - a final insertion-sort pass (guarded for the first _S_threshold=16
//     elements, unguarded after).
//
// The _S_threshold constant (16) and the median-of-3 pivot selection have been
// stable in libstdc++ for many years; the port is validated against the live
// upstream bedtools binary in the parity pipeline.
package cppsort

// threshold is libstdc++'s _S_threshold: ranges of this size or smaller are
// left for the final insertion-sort pass rather than recursed into.
const threshold = 16

// Sort orders s in place to match libstdc++ std::sort(first, last, less) for the
// same comparator. less(a, b) reports whether a must come before b (strict
// weak ordering), exactly like a C++ comparator.
func Sort[T any](s []T, less func(a, b T) bool) {
	if len(s) > 1 {
		introsortLoop(s, 0, len(s), 2*lg(len(s)), less)
		finalInsertionSort(s, less)
	}
}

// lg is libstdc++'s std::__lg: floor(log2(n)) for n >= 1 (and 0 for n <= 0),
// which seeds the introsort depth limit.
func lg(n int) int {
	k := 0
	for n > 1 {
		n >>= 1
		k++
	}
	return k
}

// introsortLoop mirrors std::__introsort_loop on the half-open range [lo, hi).
// It recurses on the right partition and loops on the left (as libstdc++ does),
// falling back to heapsort once the depth limit is exhausted.
func introsortLoop[T any](s []T, lo, hi, depthLimit int, less func(a, b T) bool) {
	for hi-lo > threshold {
		if depthLimit == 0 {
			heapSort(s, lo, hi, less)
			return
		}
		depthLimit--
		cut := unguardedPartitionPivot(s, lo, hi, less)
		introsortLoop(s, cut, hi, depthLimit, less)
		hi = cut
	}
}

// unguardedPartitionPivot mirrors std::__unguarded_partition_pivot: it moves the
// median of (first, mid, last-1) to first, then partitions [first+1, last).
func unguardedPartitionPivot[T any](s []T, lo, hi int, less func(a, b T) bool) int {
	mid := lo + (hi-lo)/2
	moveMedianToFirst(s, lo, lo+1, mid, hi-1, less)
	return unguardedPartition(s, lo+1, hi, lo, less)
}

// moveMedianToFirst mirrors std::__move_median_to_first: it places the median of
// elements a, b, c into result (all are indices).
func moveMedianToFirst[T any](s []T, result, a, b, c int, less func(x, y T) bool) {
	if less(s[a], s[b]) {
		if less(s[b], s[c]) {
			s[result], s[b] = s[b], s[result]
		} else if less(s[a], s[c]) {
			s[result], s[c] = s[c], s[result]
		} else {
			s[result], s[a] = s[a], s[result]
		}
	} else if less(s[a], s[c]) {
		s[result], s[a] = s[a], s[result]
	} else if less(s[b], s[c]) {
		s[result], s[c] = s[c], s[result]
	} else {
		s[result], s[b] = s[b], s[result]
	}
}

// unguardedPartition mirrors std::__unguarded_partition over [first, last) with
// the pivot held at index pivot. It returns the partition point.
func unguardedPartition[T any](s []T, first, last, pivot int, less func(a, b T) bool) int {
	for {
		for less(s[first], s[pivot]) {
			first++
		}
		last--
		for less(s[pivot], s[last]) {
			last--
		}
		if !(first < last) {
			return first
		}
		s[first], s[last] = s[last], s[first]
		first++
	}
}

// finalInsertionSort mirrors std::__final_insertion_sort: a guarded insertion
// sort over the first threshold elements, then an unguarded one over the rest
// (the guard is unnecessary there because a smaller element is guaranteed to
// exist to the left after introsort_loop's partitioning).
func finalInsertionSort[T any](s []T, less func(a, b T) bool) {
	n := len(s)
	if n > threshold {
		insertionSort(s, 0, threshold, less)
		for i := threshold; i < n; i++ {
			unguardedLinearInsert(s, i, less)
		}
	} else {
		insertionSort(s, 0, n, less)
	}
}

// insertionSort mirrors std::__insertion_sort over [lo, hi).
func insertionSort[T any](s []T, lo, hi int, less func(a, b T) bool) {
	if lo == hi {
		return
	}
	for i := lo + 1; i < hi; i++ {
		if less(s[i], s[lo]) {
			// val belongs at the front: rotate [lo, i) right by one.
			val := s[i]
			copy(s[lo+1:i+1], s[lo:i])
			s[lo] = val
		} else {
			unguardedLinearInsert(s, i, less)
		}
	}
}

// unguardedLinearInsert mirrors std::__unguarded_linear_insert: slide s[i] left
// past every strictly-greater neighbour. The caller guarantees a not-greater
// element exists to the left so no bounds check is needed.
func unguardedLinearInsert[T any](s []T, i int, less func(a, b T) bool) {
	val := s[i]
	j := i - 1
	for less(val, s[j]) {
		s[j+1] = s[j]
		j--
	}
	s[j+1] = val
}

// heapSort mirrors libstdc++'s introsort depth-limit fallback,
// std::__partial_sort(first, last, last) = make_heap + sort_heap over [lo, hi).
func heapSort[T any](s []T, lo, hi int, less func(a, b T) bool) {
	makeHeap(s, lo, hi, less)
	// sort_heap: repeatedly pop the max to the end.
	for last := hi; last-lo > 1; last-- {
		s[lo], s[last-1] = s[last-1], s[lo]
		adjustHeap(s, lo, lo, last-1, less)
	}
}

// makeHeap mirrors std::__make_heap over [lo, hi).
func makeHeap[T any](s []T, lo, hi int, less func(a, b T) bool) {
	n := hi - lo
	if n < 2 {
		return
	}
	for parent := (n - 2) / 2; ; parent-- {
		adjustHeap(s, lo, lo+parent, hi, less)
		if parent == 0 {
			return
		}
	}
}

// adjustHeap mirrors std::__adjust_heap: sift the element at index pos down the
// max-heap rooted in [lo, hi). It reproduces libstdc++'s "sift the hole to a
// leaf then push the value up" formulation so the comparison sequence — and
// thus the final order of equal elements — matches byte-for-byte.
func adjustHeap[T any](s []T, lo, pos, hi int, less func(a, b T) bool) {
	n := hi - lo
	holeIndex := pos - lo
	value := s[pos]
	topIndex := holeIndex
	secondChild := holeIndex
	for secondChild < (n-1)/2 {
		secondChild = 2 * (secondChild + 1)
		if less(s[lo+secondChild], s[lo+secondChild-1]) {
			secondChild--
		}
		s[lo+holeIndex] = s[lo+secondChild]
		holeIndex = secondChild
	}
	if n&1 == 0 && secondChild == (n-2)/2 {
		secondChild = 2 * (secondChild + 1)
		s[lo+holeIndex] = s[lo+secondChild-1]
		holeIndex = secondChild - 1
	}
	// push_heap the value from holeIndex up to topIndex.
	pushHeap(s, lo, holeIndex, topIndex, value, less)
}

// pushHeap mirrors std::__push_heap: bubble value up from holeIndex until it
// reaches topIndex or a parent that is not less than it.
func pushHeap[T any](s []T, lo, holeIndex, topIndex int, value T, less func(a, b T) bool) {
	parent := (holeIndex - 1) / 2
	for holeIndex > topIndex && less(s[lo+parent], value) {
		s[lo+holeIndex] = s[lo+parent]
		holeIndex = parent
		parent = (holeIndex - 1) / 2
	}
	s[lo+holeIndex] = value
}
