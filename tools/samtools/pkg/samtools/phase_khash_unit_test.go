package samtools

// Binary-free unit tests for the phase fragment khash port
// (phase_khash.go). These do NOT invoke the upstream samtools binary;
// they exercise the in-place Cuckoo-style kick-out rehash that makes the
// bucket layout — and therefore the EV-line emit order — match upstream
// after the fragment table grows past its initial 16 buckets.

import "testing"

// refKhashKickout reproduces htslib's kh_resize_##name kick-out rehash
// independently of fragKhash, using the same probe sequence and the same
// in-place swap-chain. It returns the final bucket layout (occupied keys
// in bucket order, with 0 for empty/deleted) for a table built by
// inserting `keys` in order. Comparing fragKhash's iteration order
// against this reference pins the port to htslib's exact relocation
// semantics without needing the C binary.
func refKhashKickout(keys []uint64) []uint64 {
	const hashUpper = 0.77
	var nb, size, nocc, upper uint32
	var flags []uint32 // 2 bits per bucket; 10=empty, 01=del, 00=occupied
	var kk []uint64
	isEmpty := func(f []uint32, i uint32) bool { return (f[i>>4]>>((i&0xf)<<1))&2 != 0 }
	isEither := func(f []uint32, i uint32) bool { return (f[i>>4]>>((i&0xf)<<1))&3 != 0 }
	setOcc := func(f []uint32, i uint32) { f[i>>4] &^= 3 << ((i & 0xf) << 1) }
	setEmptyFalse := func(f []uint32, i uint32) { f[i>>4] &^= 2 << ((i & 0xf) << 1) }
	setDelTrue := func(f []uint32, i uint32) { f[i>>4] |= 1 << ((i & 0xf) << 1) }
	fsize := func(m uint32) uint32 {
		if m < 16 {
			return 1
		}
		return m >> 4
	}
	resize := func(newBuckets uint32) {
		newBuckets = kroundup32(newBuckets)
		if newBuckets < 4 {
			newBuckets = 4
		}
		if float64(size) >= float64(newBuckets)*hashUpper+0.5 {
			return
		}
		nf := make([]uint32, fsize(newBuckets))
		for i := range nf {
			nf[i] = 0xaaaaaaaa
		}
		if nb < newBuckets {
			nkk := make([]uint64, newBuckets)
			copy(nkk, kk)
			kk = nkk
		}
		mask := newBuckets - 1
		for j := uint32(0); j < nb; j++ {
			if isEither(flags, j) {
				continue
			}
			key := kk[j]
			setDelTrue(flags, j)
			for {
				i := uint32(key) & mask
				step := uint32(0)
				for !isEmpty(nf, i) {
					step++
					i = (i + step) & mask
				}
				setEmptyFalse(nf, i)
				if i < nb && !isEither(flags, i) {
					key, kk[i] = kk[i], key
					setDelTrue(flags, i)
				} else {
					kk[i] = key
					break
				}
			}
		}
		if nb > newBuckets {
			kk = kk[:newBuckets]
		}
		flags = nf
		nb = newBuckets
		nocc = size
		upper = uint32(float64(newBuckets)*hashUpper + 0.5)
	}
	put := func(key uint64) {
		if nocc >= upper {
			if nb > size*2 {
				resize(nb - 1)
			} else {
				resize(nb + 1)
			}
		}
		mask := nb - 1
		x := uint32(key) & mask
		if !isEmpty(flags, x) {
			step := uint32(0)
			last := x
			site := nb
			for !isEmpty(flags, x) && ((flags[x>>4]>>((x&0xf)<<1))&1 != 0 || kk[x] != key) {
				if (flags[x>>4]>>((x&0xf)<<1))&1 != 0 {
					site = x
				}
				step++
				x = (x + step) & mask
				if x == last {
					x = site
					break
				}
			}
			if isEmpty(flags, x) && site != nb {
				x = site
			}
		}
		if isEither(flags, x) {
			setOcc(flags, x)
			kk[x] = key
			size++
			nocc++
		} else {
			kk[x] = key // already present (kept)
		}
	}
	for _, key := range keys {
		put(key)
	}
	out := make([]uint64, nb)
	for j := uint32(0); j < nb; j++ {
		if !isEither(flags, j) {
			out[j] = kk[j]
		}
	}
	return out
}

// TestUnitFragKhashKickoutLayout asserts that fragKhash's bucket layout
// (and hence iteration order) matches the independent reference kick-out
// rehash across a range of key counts that straddle the grow thresholds
// (4 -> 8 -> 16 -> 32 -> 64 buckets).
func TestUnitFragKhashKickoutLayout(t *testing.T) {
	for _, n := range []int{1, 2, 3, 5, 6, 7, 12, 13, 16, 24, 25, 40, 49, 50, 64, 100} {
		keys := make([]uint64, n)
		for i := 0; i < n; i++ {
			// X31-ish spread so keys collide in low bits like real QNAME
			// hashes do, exercising the probe chains.
			keys[i] = x31HashString("read_" + itoa(i))
		}
		want := refKhashKickout(keys)

		h := newFragKhash()
		for _, key := range keys {
			h.put(key)
		}
		if int(h.end()) != len(want) {
			t.Fatalf("n=%d: bucket count %d != reference %d", n, h.end(), len(want))
		}
		for j := uint32(0); j < h.end(); j++ {
			var got uint64
			if h.exist(j) {
				got = h.keys[j]
			}
			if got != want[j] {
				t.Errorf("n=%d bucket %d: got key %d, want %d", n, j, got, want[j])
			}
		}
	}
}

// TestUnitFragKhashLookupAfterGrow verifies every inserted key remains
// findable after the table has grown through several rehashes, i.e. the
// kick-out relocation never loses an entry.
func TestUnitFragKhashLookupAfterGrow(t *testing.T) {
	h := newFragKhash()
	keys := make([]uint64, 80)
	for i := range keys {
		keys[i] = x31HashString("q" + itoa(i))
		bk, isNew := h.put(keys[i])
		if !isNew {
			t.Fatalf("key %d unexpectedly already present", i)
		}
		h.vals[bk].vpos = int32(i) // tag so we can verify identity
	}
	for i, key := range keys {
		bk, ok := h.get(key)
		if !ok {
			t.Fatalf("key %d (%d) not found after grow", i, key)
		}
		if int(h.vals[bk].vpos) != i {
			t.Errorf("key %d resolved to wrong value: vpos=%d", i, h.vals[bk].vpos)
		}
	}
	if int(h.size) != len(keys) {
		t.Errorf("size=%d, want %d", h.size, len(keys))
	}
}
