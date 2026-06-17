package cppsort

import (
	"math/rand"
	"sort"
	"testing"
)

// TestSortMatchesTotalOrder checks that for a TOTAL-order comparator cppsort
// produces a fully sorted result identical to sort.Ints, across a range of
// sizes that exercise the insertion-sort-only path (n<=16), the introsort
// recursion (n>16), and large inputs that trigger the heapsort depth-limit
// fallback.
func TestSortMatchesTotalOrder(t *testing.T) {
	rng := rand.New(rand.NewSource(12345))
	for _, n := range []int{0, 1, 2, 5, 16, 17, 50, 257, 1000, 5000} {
		for trial := 0; trial < 5; trial++ {
			a := make([]int, n)
			for i := range a {
				a[i] = rng.Intn(n + 1) // duplicates on purpose
			}
			b := append([]int(nil), a...)
			Sort(a, func(x, y int) bool { return x < y })
			sort.Ints(b)
			for i := range a {
				if a[i] != b[i] {
					t.Fatalf("n=%d trial=%d: not sorted at %d: got %v want %v", n, trial, i, a[i], b[i])
				}
			}
		}
	}
}

// TestSortStrictDescending verifies a descending comparator also fully orders.
func TestSortStrictDescending(t *testing.T) {
	a := []int{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6}
	Sort(a, func(x, y int) bool { return x > y })
	for i := 1; i < len(a); i++ {
		if a[i-1] < a[i] {
			t.Fatalf("not descending at %d: %v", i, a)
		}
	}
}

// TestSortIsPermutation confirms no elements are lost or duplicated (a hazard
// when porting the in-place partition/heap index arithmetic).
func TestSortIsPermutation(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	for _, n := range []int{17, 64, 333, 2048} {
		a := make([]int, n)
		for i := range a {
			a[i] = rng.Int()
		}
		before := map[int]int{}
		for _, v := range a {
			before[v]++
		}
		Sort(a, func(x, y int) bool { return x < y })
		after := map[int]int{}
		for _, v := range a {
			after[v]++
		}
		if len(before) != len(after) {
			t.Fatalf("n=%d: distinct count changed %d -> %d", n, len(before), len(after))
		}
		for k, c := range before {
			if after[k] != c {
				t.Fatalf("n=%d: multiplicity of %d changed %d -> %d", n, k, c, after[k])
			}
		}
	}
}
