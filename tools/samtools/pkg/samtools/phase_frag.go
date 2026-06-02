package samtools

// phase fragment-map types and helpers, mirroring frag_t and the
// X31_hash_string / rseq_lt declarations in reference_code/samtools/phase.c.

// maxVars matches phase.c's MAX_VARS — the in-place capacity of frag.seq.
const maxVars = 256

// frag is the per-read structure, a Go port of `frag_t` (phase.c:69).
// Field semantics:
//
//	seq[0..vlen)   per-het allele code: 0 = ambiguous (read didn't
//	               make a confident base call at this het), 1 = base
//	               matched cns&3 (j, the larger ACGT code), 2 = base
//	               matched cns>>16&3 (i, the smaller ACGT code).
//	vpos           het-array index of seq[0] (block-local; offset by
//	               g.vpos_shift to recover the global het index).
//	beg, end       0-based reference start and exclusive end of the
//	               read; from bam_get_pos / bam_endpos.
//	vlen           valid length of seq (1..maxVars).
//	single         set when this frag has a single het (no phase
//	               information).
//	flip, phase, phased, ambig — fragphase flags (see phase.c:228-232).
//	in, out        per-fragment hap-agreement counts.
type frag struct {
	seq    [maxVars]int8
	vpos   int32
	beg    int32
	end    int32
	vlen   uint16
	single uint8 // 0 or 1
	flip   uint8
	phase  uint8
	phased uint8
	ambig  uint8
	in     uint16
	out    uint16
}

// x31HashString matches phase.c's static inline X31_hash_string:
//
//	uint64 h = *s;
//	for (++s; *s; ++s) h = (h << 5) - h + *s;
//
// — i.e. h(0)=s[0], then h(n) = 31*h(n-1) + s[n].
func x31HashString(s string) uint64 {
	if len(s) == 0 {
		return 0
	}
	h := uint64(s[0])
	for i := 1; i < len(s); i++ {
		h = (h << 5) - h + uint64(s[i])
	}
	return h
}

// fragPtr wraps a *frag for use with the ks_introsort port. It also
// carries the original bucket index so EV-line output (which prints the
// frag's contents) is reproducible against upstream khash iteration
// order. fragPtr is the "type_t" in `KSORT_INIT(rseq, frag_p, rseq_lt)`.
type fragPtr struct {
	f      *frag
	bucket uint32
}

// fragRseqLt is the comparator used by ks_introsort_rseq. Mirrors the
// macro `rseq_lt(a,b) ((a)->vpos < (b)->vpos)` in phase.c.
func fragRseqLt(a, b fragPtr) bool { return a.f.vpos < b.f.vpos }
