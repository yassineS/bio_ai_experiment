package codec

import (
	"bytes"
	"math/rand"
	"sync"
	"testing"
)

// TestRANS4x8O1PoolConcurrentRoundTrip stresses the pooled order-1 encoder
// scratch (ransO1ScratchPool) from many goroutines. The pool reuses the large
// additive frequency/cumulative tables across blocks, so a missing zero-reset
// or an aliased scratch would corrupt one goroutine's frequency model and break
// its round-trip. Encoding a distinct, statistically-varied payload per
// iteration and decoding it back verifies both that reset() clears every stale
// count and that no two concurrent encodes share live scratch. This is the
// data-race companion to `go test -race`.
func TestRANS4x8O1PoolConcurrentRoundTrip(t *testing.T) {
	const goroutines = 8
	const iters = 40

	var wg sync.WaitGroup
	errs := make(chan error, goroutines*iters)
	for g := 0; g < goroutines; g++ {
		wg.Add(1)
		go func(seed int64) {
			defer wg.Done()
			r := rand.New(rand.NewSource(seed))
			for i := 0; i < iters; i++ {
				// Vary the alphabet size and length so different contexts are
				// exercised — this maximises the chance a stale table entry
				// from a previous block would show up.
				n := 100 + r.Intn(20000)
				alpha := 1 + r.Intn(255)
				in := make([]byte, n)
				for j := range in {
					in[j] = byte(r.Intn(alpha))
				}
				comp, err := RANS4x8Encode(in, 1)
				if err != nil {
					errs <- err
					return
				}
				got, err := RANS4x8Decode(comp)
				if err != nil {
					errs <- err
					return
				}
				if !bytes.Equal(got, in) {
					errs <- errRoundTrip(len(in), len(got))
					return
				}
			}
		}(int64(g) + 1)
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		if err != nil {
			t.Fatalf("concurrent order-1 round-trip failed: %v", err)
		}
	}
}

// TestRANS4x8O1PoolReuseDeterministic verifies the pooled scratch produces
// byte-identical output when the same input is encoded twice in a row (so the
// second encode reuses a scratch the first returned to the pool). A stale entry
// surviving reset() would perturb the frequency model and change the bytes.
func TestRANS4x8O1PoolReuseDeterministic(t *testing.T) {
	inputs := [][]byte{
		bytes.Repeat([]byte("ACGT"), 5000),
		[]byte("the quick brown fox " + repeat("ACGTN", 400)),
		fullAlphabet(9000),
	}
	for _, in := range inputs {
		// Prime the pool with an unrelated, differently-shaped payload so the
		// scratch reused on the real encode carries foreign counts that reset()
		// must clear.
		prime := bytes.Repeat([]byte{7, 200, 42}, 3000)
		if _, err := RANS4x8Encode(prime, 1); err != nil {
			t.Fatalf("prime encode: %v", err)
		}
		first, err := RANS4x8Encode(in, 1)
		if err != nil {
			t.Fatalf("first encode: %v", err)
		}
		if _, err := RANS4x8Encode(prime, 1); err != nil {
			t.Fatalf("re-prime encode: %v", err)
		}
		second, err := RANS4x8Encode(in, 1)
		if err != nil {
			t.Fatalf("second encode: %v", err)
		}
		if !bytes.Equal(first, second) {
			t.Fatalf("pooled scratch not reset: two encodes of the same input differ (%d vs %d bytes)",
				len(first), len(second))
		}
		got, err := RANS4x8Decode(first)
		if err != nil {
			t.Fatalf("decode: %v", err)
		}
		if !bytes.Equal(got, in) {
			t.Fatalf("round-trip mismatch after pool reuse")
		}
	}
}

type roundTripErr struct{ want, got int }

func (e roundTripErr) Error() string {
	return "round-trip length mismatch: want " + itoa(e.want) + " got " + itoa(e.got)
}

func errRoundTrip(want, got int) error { return roundTripErr{want, got} }
