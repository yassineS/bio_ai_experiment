package samtools

// Binary-free unit tests for classifyDumpAln — the pure per-read →
// bucket routing logic ported from phase.c::dump_aln (phase.c:374-388).
// These run with submodules unpopulated (no upstream binary).

import "testing"

// fakeRNG is a phaseRNG that replays a fixed sequence of Float64 values,
// so the routing decisions that consult the RNG (which=3 random path and
// the is_flip<0.5 shuffle) are exercised deterministically.
type fakeRNG struct {
	vals []float64
	i    int
}

func (f *fakeRNG) Float64() float64 {
	if f.i >= len(f.vals) {
		return 0.0
	}
	v := f.vals[f.i]
	f.i++
	return v
}

// putFrag inserts a frag with the given flags under a synthetic key and
// returns that key, so classifyDumpAln can look it up.
func putFrag(h *fragKhash, key uint64, phased, flip, ambig, phase uint8) {
	k, _ := h.put(key)
	f := &h.vals[k]
	f.phased = phased
	f.flip = flip
	f.ambig = ambig
	f.phase = phase
}

func TestUnitClassifyDumpAln(t *testing.T) {
	const (
		keyConfident0 = uint64(1001)
		keyConfident1 = uint64(1002)
		keyChimera    = uint64(1003)
		keyUnphased   = uint64(1004)
		keyAmbig      = uint64(1005)
		keyAbsent     = uint64(9999)
	)

	tests := []struct {
		name       string
		key        uint64
		isFlip     bool
		dropAmbi   bool
		rng        []float64 // consulted only by the which=3 random fallback
		wantBucket int
		wantZP     bool
	}{
		{
			// Confident hap0, no is_flip: which = phase = 0, ZP annotated.
			name: "confident_hap0_noflip", key: keyConfident0,
			wantBucket: 0, wantZP: true,
		},
		{
			// Confident hap1, no is_flip: which = phase = 1.
			name: "confident_hap1_noflip", key: keyConfident1,
			wantBucket: 1, wantZP: true,
		},
		{
			// Confident hap0 WITH is_flip: which < 2 so 1 - which = 1.
			name: "confident_hap0_flip", key: keyConfident0, isFlip: true,
			wantBucket: 1, wantZP: true,
		},
		{
			// Confident hap1 WITH is_flip: 1 - 1 = 0.
			name: "confident_hap1_flip", key: keyConfident1, isFlip: true,
			wantBucket: 0, wantZP: true,
		},
		{
			// phased && flip → which = 2 (chimera). is_flip does NOT touch
			// bucket 2 (the which<2 guard). No ZP.
			name: "chimera_flip", key: keyChimera, isFlip: true,
			wantBucket: 2, wantZP: false,
		},
		{
			// phased == 0 → which = 3 → random. Upstream phase.c:388
			// `which = (drand48() < 0.5)`: rng < 0.5 → bucket 1.
			name: "unphased_random_low", key: keyUnphased,
			rng: []float64{0.25}, wantBucket: 1, wantZP: false,
		},
		{
			// phased == 0 → which = 3 → random. rng >= 0.5 → bucket 0.
			name: "unphased_random_high", key: keyUnphased,
			rng: []float64{0.75}, wantBucket: 0, wantZP: false,
		},
		{
			// ambig, NOT dropAmbi → which = 3 → random. rng < 0.5 → 1.
			name: "ambig_keep", key: keyAmbig,
			rng: []float64{0.25}, wantBucket: 1, wantZP: false,
		},
		{
			// ambig WITH dropAmbi → which = 2 (chimera), no RNG consulted.
			name: "ambig_drop", key: keyAmbig, dropAmbi: true,
			wantBucket: 2, wantZP: false,
		},
		{
			// Absent key → which = 3 → random. rng >= 0.5 → bucket 0.
			name: "absent_random", key: keyAbsent,
			rng: []float64{0.75}, wantBucket: 0, wantZP: false,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			h := newFragKhash()
			putFrag(h, keyConfident0, 1, 0, 0, 0)
			putFrag(h, keyConfident1, 1, 0, 0, 1)
			putFrag(h, keyChimera, 1, 1, 0, 0)
			putFrag(h, keyUnphased, 0, 0, 0, 0)
			putFrag(h, keyAmbig, 1, 0, 1, 0)

			rng := &fakeRNG{vals: tc.rng}
			bucket, zp := classifyDumpAln(tc.key, h, tc.isFlip, tc.dropAmbi, rng)
			if bucket != tc.wantBucket {
				t.Errorf("bucket = %d, want %d", bucket, tc.wantBucket)
			}
			if zp != tc.wantZP {
				t.Errorf("addZP = %v, want %v", zp, tc.wantZP)
			}
		})
	}
}

// TestUnitClassifyDumpAln_AmbigPrecedence verifies the ambig branch is
// checked BEFORE the flip branch (phase.c:378 orders the switch so an
// ambiguous read never reaches the phased&&flip arm). A read that is both
// ambig and flip must route by the ambig rule.
func TestUnitClassifyDumpAln_AmbigPrecedence(t *testing.T) {
	const key = uint64(42)
	h := newFragKhash()
	putFrag(h, key, 1, 1, 1, 0) // phased && flip && ambig

	// dropAmbi=false → ambig path → which=3 → random. rng 0.1 < 0.5 → 1.
	rng := &fakeRNG{vals: []float64{0.1}}
	if b, _ := classifyDumpAln(key, h, false, false, rng); b != 1 {
		t.Errorf("ambig+flip, keep: bucket = %d, want 1 (random low)", b)
	}
	// dropAmbi=true → ambig path → which=2.
	if b, _ := classifyDumpAln(key, h, false, true, &fakeRNG{}); b != 2 {
		t.Errorf("ambig+flip, drop: bucket = %d, want 2 (chimera)", b)
	}
}
