// Package samtools — minimal port of htslib's khash KHASH_MAP_INIT_INT64
// hash table, just enough to reproduce upstream samtools `phase`
// iteration order over the per-read fragment map.
//
// This is NOT a general khash port. It mirrors khash.h's bucket layout
// and resize policy bit-for-bit so that the EV-line order emitted by
// our phase port matches upstream samtools byte-for-byte on inputs
// that exercise the natural insertion/iteration path. The relevant
// invariants:
//
//   - n_buckets is always a power of two, computed via kroundup32; the
//     initial size is 0 (lazy-allocated on first put).
//   - upper_bound = floor(n_buckets * 0.77 + 0.5).
//   - Every kh_put call (existing key or new) first checks
//     `n_occupied >= upper_bound` and triggers a grow-resize to
//     `kroundup32(n_buckets + 1)`. This is the rule that produces the
//     16-bucket layout for 6 distinct keys × 2 puts each on the
//     canonical phase fixture.
//   - Linear probing with mask = n_buckets - 1. The probe sequence is
//     `i = (i + (++step)) & mask`, i.e. quadratic-by-step (1, 2, 3, ...
//     — not the traditional doubling — exactly as in khash.h).
//   - Iteration walks buckets 0..n_buckets-1 and yields occupied,
//     non-deleted slots in bucket order.
//
// The package only ports the integer-keyed map flavour; the set64
// variant used by phase.c's `-l` site list isn't needed for the
// default code path and is not implemented here.
package samtools

// khashFlagEmpty and khashFlagDel mirror khash.h's two-bit flag layout.
// Each bucket has two flag bits packed into a 32-bit word, four buckets
// per word. The flag pattern `10` means "empty" (the default) and `01`
// means "deleted". `00` means "occupied".
const (
	khashFlagEmpty byte = 2
	khashFlagDel   byte = 1
)

// fragKhash is the int64→frag table used by the upstream-style phase
// emitter. Keys are X31_hash_string(QNAME) values, values are frag
// structs allocated inline.
type fragKhash struct {
	nBuckets   uint32
	size       uint32
	nOccupied  uint32
	upperBound uint32
	flags      []uint32 // 4 flag bits per byte ×4 ≈ 2 bits per bucket
	keys       []uint64
	vals       []frag
}

// newFragKhash returns an empty table; the first kh_put will allocate.
func newFragKhash() *fragKhash { return &fragKhash{} }

// khHash64 is khash.h's kh_int64_hash_func for 64-bit keys:
//
//	(khint32_t)((key)>>33 ^ (key) ^ (key)<<11)
//
// phase.c instantiates its fragment table with KHASH_MAP_INIT_INT64, so the
// initial bucket of a key is this mixed hash masked by (nBuckets-1), NOT the
// raw low 32 bits. Reproducing the exact mix is required for the bucket layout
// — and therefore the unstable ks_introsort_rseq permutation that orders the
// EV supporting-read lines — to match upstream samtools byte-for-byte.
func khHash64(key uint64) uint32 {
	return uint32((key >> 33) ^ key ^ (key << 11))
}

// kroundup32 mirrors the macro of the same name in htslib.
func kroundup32(x uint32) uint32 {
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	x++
	return x
}

// fsize returns the flag-array length (in 32-bit words) for an m-bucket
// table; matches __ac_fsize(m) = (m<16 ? 1 : m>>4).
func khFsize(m uint32) uint32 {
	if m < 16 {
		return 1
	}
	return m >> 4
}

// flagIsEmpty / flagIsDel / flagIsEither match the bit-packing in khash.h.
func khFlagIsEmpty(flags []uint32, i uint32) bool {
	return (flags[i>>4]>>((i&0xf)<<1))&2 != 0
}
func khFlagIsDel(flags []uint32, i uint32) bool {
	return (flags[i>>4]>>((i&0xf)<<1))&1 != 0
}
func khFlagIsEither(flags []uint32, i uint32) bool {
	return (flags[i>>4]>>((i&0xf)<<1))&3 != 0
}
func khFlagSetIsBothFalse(flags []uint32, i uint32) {
	flags[i>>4] &^= 3 << ((i & 0xf) << 1)
}

// khFlagSetIsEmptyFalse clears the "empty" bit (bit 2) for bucket i,
// matching khash.h's __ac_set_isempty_false. Used by the in-place
// kick-out rehash to claim a slot in the new flag array.
func khFlagSetIsEmptyFalse(flags []uint32, i uint32) {
	flags[i>>4] &^= 2 << ((i & 0xf) << 1)
}

// khFlagSetIsDelTrue sets the "deleted" bit (bit 1) for bucket i,
// matching khash.h's __ac_set_isdel_true. The kick-out rehash uses it
// to mark an old-table slot vacated during relocation.
func khFlagSetIsDelTrue(flags []uint32, i uint32) {
	flags[i>>4] |= 1 << ((i & 0xf) << 1)
}

// resize grows (or shrinks) the table to at least new_n_buckets
// (rounded up to a power of two). This is a faithful port of
// kh_resize_##name in khash.h, including the in-place Cuckoo-style
// "kick-out" rehash. Reproducing that exact relocation order — rather
// than a clean rehash into a fresh array — is what makes the resulting
// bucket layout (and therefore the EV-line emit order, which is the
// unstable ks_introsort_rseq permutation of equal-vpos frags) match
// upstream samtools byte-for-byte even after the table grows.
func (h *fragKhash) resize(newBuckets uint32) {
	const hashUpper = 0.77
	// kroundup32 first, then clamp to >=4 — matches khash.h order.
	newBuckets = kroundup32(newBuckets)
	if newBuckets < 4 {
		newBuckets = 4
	}
	// "requested size is too small" → do nothing (j = 0 in khash.h).
	if float64(h.size) >= float64(newBuckets)*hashUpper+0.5 {
		return
	}
	newFlags := make([]uint32, khFsize(newBuckets))
	for i := range newFlags {
		newFlags[i] = 0xaaaaaaaa // memset 0xaa → all-empty
	}
	// Expand: keep existing keys/vals in place, grow the backing arrays.
	// (Shrink reuses the same arrays and truncates at the end.)
	if h.nBuckets < newBuckets {
		newKeys := make([]uint64, newBuckets)
		copy(newKeys, h.keys)
		h.keys = newKeys
		newVals := make([]frag, newBuckets)
		copy(newVals, h.vals)
		h.vals = newVals
	}
	newMask := newBuckets - 1
	// Kick-out rehash over the OLD buckets (h.nBuckets), reading the OLD
	// flags (h.flags) and writing occupancy into newFlags. h.keys/h.vals
	// are the combined (extended) arrays; relocations swap through them.
	for j := uint32(0); j < h.nBuckets; j++ {
		if khFlagIsEither(h.flags, j) {
			continue
		}
		key := h.keys[j]
		val := h.vals[j]
		khFlagSetIsDelTrue(h.flags, j) // vacate old slot
		for {
			i := khHash64(key) & newMask
			step := uint32(0)
			for !khFlagIsEmpty(newFlags, i) {
				step++
				i = (i + step) & newMask
			}
			khFlagSetIsEmptyFalse(newFlags, i)
			if i < h.nBuckets && khFlagIsEither(h.flags, i) == false {
				// Kick out the element currently occupying slot i.
				key, h.keys[i] = h.keys[i], key
				val, h.vals[i] = h.vals[i], val
				khFlagSetIsDelTrue(h.flags, i) // mark old slot deleted
			} else {
				// Empty (or beyond old range): write and stop.
				h.keys[i] = key
				h.vals[i] = val
				break
			}
		}
	}
	// Shrink: truncate the backing arrays to newBuckets.
	if h.nBuckets > newBuckets {
		h.keys = h.keys[:newBuckets]
		h.vals = h.vals[:newBuckets]
	}
	h.flags = newFlags
	h.nBuckets = newBuckets
	h.nOccupied = h.size
	h.upperBound = uint32(float64(newBuckets)*hashUpper + 0.5)
}

// put inserts `key` and returns (bucket index, isNew). isNew is true
// when the key was not previously present (and the bucket is now
// occupied with default-initialised value). Callers then read/write
// h.vals[bucket].
//
// Matches kh_put_##name in khash.h.
func (h *fragKhash) put(key uint64) (uint32, bool) {
	// Resize check on every put.
	if h.nOccupied >= h.upperBound {
		if h.nBuckets > h.size*2 {
			h.resize(h.nBuckets - 1) // clear deleted
		} else {
			h.resize(h.nBuckets + 1)
		}
	}
	mask := h.nBuckets - 1
	x := khHash64(key) & mask
	step := uint32(0)
	site := h.nBuckets
	// Probe for either:
	//   - an empty slot (insert here),
	//   - a deleted slot (remember the first one, keep probing for an
	//     existing matching key), or
	//   - an occupied slot with our key (found, no-op).
	if khFlagIsEmpty(h.flags, x) {
		// fast path
	} else {
		last := x
		for !khFlagIsEmpty(h.flags, x) && (khFlagIsDel(h.flags, x) || h.keys[x] != key) {
			if khFlagIsDel(h.flags, x) {
				site = x
			}
			step++
			x = (x + step) & mask
			if x == last {
				x = site
				break
			}
		}
		if khFlagIsEmpty(h.flags, x) && site != h.nBuckets {
			x = site
		}
	}
	if khFlagIsEmpty(h.flags, x) {
		// brand new entry
		h.keys[x] = key
		h.vals[x] = frag{}
		khFlagSetIsBothFalse(h.flags, x)
		h.size++
		h.nOccupied++
		return x, true
	}
	if khFlagIsDel(h.flags, x) {
		// recycle deleted slot
		h.keys[x] = key
		h.vals[x] = frag{}
		khFlagSetIsBothFalse(h.flags, x)
		h.size++
		return x, true
	}
	// key already present
	return x, false
}

// get looks up key without inserting. Returns the bucket index and
// true when the key is present, or (nBuckets, false) when it isn't.
// Matches kh_get_##name in khash.h.
func (h *fragKhash) get(key uint64) (uint32, bool) {
	if h.nBuckets == 0 {
		return 0, false
	}
	mask := h.nBuckets - 1
	x := khHash64(key) & mask
	step := uint32(0)
	last := x
	for !khFlagIsEmpty(h.flags, x) && (khFlagIsDel(h.flags, x) || h.keys[x] != key) {
		step++
		x = (x + step) & mask
		if x == last {
			return h.nBuckets, false
		}
	}
	if khFlagIsEither(h.flags, x) {
		return h.nBuckets, false
	}
	return x, true
}

// del marks bucket k as deleted. Matches kh_del_##name in khash.h.
func (h *fragKhash) del(k uint32) {
	if k >= h.nBuckets || khFlagIsEither(h.flags, k) {
		return
	}
	h.flags[k>>4] |= uint32(khashFlagDel) << ((k & 0xf) << 1)
	h.size--
}

// exist reports whether bucket k currently holds an entry.
func (h *fragKhash) exist(k uint32) bool {
	if k >= h.nBuckets {
		return false
	}
	return !khFlagIsEither(h.flags, k)
}

// end returns kh_end(h) — the number of buckets, suitable as the
// upper bound for an iteration `for k := 0; k < h.end(); k++`.
func (h *fragKhash) end() uint32 { return h.nBuckets }
