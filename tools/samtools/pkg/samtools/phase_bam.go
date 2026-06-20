package samtools

import (
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// bamSplitWriter owns the three per-haplotype BAM files written when
// `phase -b STR` is in effect. Files are opened lazily on construction
// and closed by the caller via Close.
type bamSplitWriter struct {
	files   [3]*os.File
	writers [3]*sam.BAMWriter
	names   [3]string
}

// newBAMSplitWriter opens `<prefix>.0.bam`, `<prefix>.1.bam` and
// `<prefix>.chimera.bam`, writes the BAM header (a copy of hdr) into
// each, and returns the bundle. On any error it closes whatever files
// it has already opened.
func newBAMSplitWriter(prefix string, hdr *sam.Header) (*bamSplitWriter, error) {
	suffix := [3]string{"0", "1", "chimera"}
	w := &bamSplitWriter{}
	for i, s := range suffix {
		name := fmt.Sprintf("%s.%s.bam", prefix, s)
		f, err := os.Create(name)
		if err != nil {
			w.cleanupOnError(i)
			return nil, fmt.Errorf("samtools phase: create %s: %w", name, err)
		}
		bw := sam.NewBAMWriter(f)
		if err := bw.WriteHeader(hdr); err != nil {
			f.Close()
			w.cleanupOnError(i)
			return nil, fmt.Errorf("samtools phase: write header to %s: %w", name, err)
		}
		w.files[i] = f
		w.writers[i] = bw
		w.names[i] = name
	}
	return w, nil
}

func (w *bamSplitWriter) cleanupOnError(upTo int) {
	for i := 0; i < upTo; i++ {
		if w.writers[i] != nil {
			_ = w.writers[i].Close()
		}
		if w.files[i] != nil {
			_ = w.files[i].Close()
		}
	}
}

// Close flushes and closes the three BAM streams. The first error
// encountered is returned; remaining streams are still closed.
func (w *bamSplitWriter) Close() error {
	var firstErr error
	for i := 0; i < 3; i++ {
		if w.writers[i] != nil {
			if err := w.writers[i].Close(); err != nil && firstErr == nil {
				firstErr = err
			}
		}
		if w.files[i] != nil {
			if err := w.files[i].Close(); err != nil && firstErr == nil {
				firstErr = err
			}
		}
	}
	return firstErr
}

// assignAndWrite classifies every record in recs against the phased
// het list and writes the record to one of the three BAM outputs.
//
// The classification mirrors upstream `phase.c` `fragphase` followed
// by `dump_aln`:
//
//   - Each read's allele calls at the phased hets are translated into
//     a per-haplotype vote (`c[0]` = votes for hap0, `c[1]` = votes
//     for hap1) using the same convention as the greedy chaining: a
//     read whose previously-seen allele at the first phased het was
//     `allele0` is on hap0, otherwise on hap1.
//   - A read with no allele calls at any phased het is routed
//     randomly to bucket 0 or 1 (matches upstream's drand48 path).
//   - A read with strong support on one haplotype (`c[majority] > 0`,
//     `c[minority] == 0`) goes to that haplotype's bucket.
//   - A read with non-trivial support on both haplotypes is examined
//     for a chimera flip point (the `fragphase` flip-point search,
//     guarded by `c[0] >= 3 && c[1] >= 3`). If a flip point survives
//     the FLIP_THRES/FLIP_PENALTY check the read is classified as
//     chimeric (bucket 2). Otherwise the majority bucket wins,
//     unless DropAmbiguous is set, in which case "ambiguous" reads
//     (small support on both sides) are also sent to bucket 2.
//   - When chimera repair is disabled (NoFixChimera), the flip-point
//     search is skipped — the read is sent to its majority bucket
//     regardless.
func (w *bamSplitWriter) assignAndWrite(
	recs []*sam.Record,
	hets []het,
	hap0AtHet hetHapMapping,
	opts PhaseOptions,
	rng phaseRNG,
) error {
	nReads := len(recs)
	nHets := len(hets)
	// Build a (read × het) matrix of allele calls. -1 = no call at
	// this het, 0 or 1 = allele index. Most cells will be -1; storing
	// densely is still cheap because both dimensions are O(thousands)
	// at most for phase's intended workload.
	readAllele := make([][]int8, nReads)
	for i := range readAllele {
		readAllele[i] = make([]int8, nHets)
		for j := range readAllele[i] {
			readAllele[i][j] = -1
		}
	}
	for hi, h := range hets {
		for _, s := range h.support {
			readAllele[s.readIdx][hi] = int8(s.allele)
		}
	}

	// Classify each read.
	for ri, rec := range recs {
		c0, c1 := 0, 0
		// First a quick majority pass.
		for hi := 0; hi < nHets; hi++ {
			a := readAllele[ri][hi]
			if a < 0 {
				continue
			}
			h0 := hap0AtHet[hi]
			if h0 < 0 {
				continue
			}
			// Which hap does this read's allele at this het support?
			// allele a == hap0 if a == h0, else hap1.
			if int8(a) == h0 {
				c0++
			} else {
				c1++
			}
		}
		var bucket int
		switch {
		case c0 == 0 && c1 == 0:
			// No phased-het evidence — random hap.
			if rng.Float64() < 0.5 {
				bucket = 0
			} else {
				bucket = 1
			}
		case c0 > 0 && c1 == 0:
			bucket = 0
		case c1 > 0 && c0 == 0:
			bucket = 1
		default:
			// Both haplotypes have support; chimera candidate. Apply
			// the chimera-repair flip-point search (upstream fragphase
			// `flip=1` branch) when enabled and counts are large
			// enough to be worth examining.
			noFix := opts.NoFixChimera || opts.FullRead
			chimera := false
			if !noFix && c0 >= 3 && c1 >= 3 {
				chimera = hasChimeraFlipPoint(readAllele[ri], hap0AtHet)
			} else if !noFix {
				// Light-weight ambiguity criterion mirroring
				// upstream's `f->ambig` flag: small support on the
				// minority side (≤ 2) and the majority barely
				// dominates. Send these to chimera when -A is set,
				// else to the majority bucket.
				minor := c0
				major := c1
				if c0 > c1 {
					minor, major = c1, c0
				}
				if opts.DropAmbiguous && minor > 0 && minor < 3 && major <= minor+1 {
					chimera = true
				}
			}
			if chimera {
				bucket = 2
			} else if c0 >= c1 {
				bucket = 0
			} else {
				bucket = 1
			}
		}
		if err := w.writers[bucket].Write(rec); err != nil {
			return fmt.Errorf("samtools phase: write %s: %w", w.names[bucket], err)
		}
	}
	return nil
}

// dumpAln drains the front of `queue` (starting at `cursor`) for
// reads whose alignment-end is at or before `minPos` (0-based ref
// position of the next not-yet-emitted het). For each drained read
// it looks up the qname's frag in `hash` and routes the record to
// one of the three BAM outputs per upstream `dump_aln` (phase.c:361):
//
//   - frag absent → which=3 → random hap (drand48 in upstream;
//     math/rand here, see PhaseOptions.RNGSeed).
//   - frag.ambig: drop_ambi ? bucket=2 : bucket=3 (random hap).
//   - frag.phased && frag.flip: bucket=2 (chimera).
//   - frag.phased==0: bucket=3 (random hap).
//   - else: bucket=frag.phase. Annotates ZP:A:Y (matches upstream
//     phase.c:384 `bam_aux_append(b, "ZP", 'A', 1, ...)`).
//
// Upstream additionally re-flips a non-chimeric read with 50/50
// probability (phase.c:386 `which = 1 - which` if is_flip). That
// re-flip is also reproduced — it's the upstream call's evidence-
// agnostic shuffle; with the same seed it lands the same way each
// run, and across two upstream invocations it may diverge by RNG
// state (rejection parity).
//
// Returns the new cursor (the number of reads dumped from the head).
//
// dropAmbi mirrors the FLAG_DROP_AMBI bit (upstream's -A).
func dumpAln(
	queue []*sam.Record,
	cursor int,
	minPos int32,
	hash *fragKhash,
	w *bamSplitWriter,
	rng phaseRNG,
	dropAmbi bool,
) (int, error) {
	if w == nil {
		return cursor, nil
	}
	// Upstream picks is_flip ONCE per dump_aln call (phase.c:365).
	isFlip := rng.Float64() < 0.5
	for cursor < len(queue) {
		rec := queue[cursor]
		// bam_endpos in 0-based-exclusive numerically equals our
		// Record.EndPosition() in 1-based-inclusive. The upstream
		// condition is `end > min_pos: break` so we drain while
		// EndPosition() <= minPos.
		if int32(rec.EndPosition()) > minPos {
			break
		}
		key := x31HashString(rec.QName)
		bucket, addZP := classifyDumpAln(key, hash, isFlip, dropAmbi, rng)
		if addZP {
			rec.Aux = append(rec.Aux, sam.Aux{Tag: "ZP", Type: 'A', Value: "Y"})
		}
		if err := w.writers[bucket].Write(rec); err != nil {
			return cursor, fmt.Errorf("samtools phase: write %s: %w", w.names[bucket], err)
		}
		cursor++
	}
	return cursor, nil
}

// classifyDumpAln runs upstream phase.c::dump_aln's per-read
// switch (phase.c:374-388) on the frag identified by `key`. Returns
// the destination bucket (0/1/2) and whether the record should be
// annotated with ZP:A:Y (true only for confidently phased reads).
func classifyDumpAln(
	key uint64,
	hash *fragKhash,
	isFlip bool,
	dropAmbi bool,
	rng phaseRNG,
) (int, bool) {
	bucket := 3
	addZP := false
	if k, ok := hash.get(key); ok {
		f := &hash.vals[k]
		switch {
		case f.ambig != 0:
			if dropAmbi {
				bucket = 2
			} else {
				bucket = 3
			}
		case f.phased != 0 && f.flip != 0:
			bucket = 2
		case f.phased == 0:
			bucket = 3
		default:
			// phased and not flipped — confident haplotype call.
			bucket = int(f.phase)
			addZP = true
		}
		// Upstream phase.c:386: `if (which < 2 && is_flip) which = 1 - which`.
		if bucket < 2 && isFlip {
			bucket = 1 - bucket
		}
	}
	if bucket == 3 {
		// Evidence-less read: route 50/50 via the RNG.
		if rng.Float64() < 0.5 {
			bucket = 1
		} else {
			bucket = 0
		}
	}
	return bucket, addZP
}

// hasChimeraFlipPoint searches for a per-read flip point that would
// turn a chimeric read into two haplotype-consistent halves. It is
// the Go port of the FLIP_THRES/FLIP_PENALTY scoring in upstream
// `fragphase` (samtools phase.c).
//
// For each candidate split index i (0 ≤ i < nHets-1), it computes
//
//	score_md0 = (#hap0-agreements in [0..i])      +
//	            (#hap1-agreements in [i+1..end])  -
//	            (#hap0-agreements in [i+1..end]) * FLIP_PENALTY
//	score_md1 = (#hap1-agreements in [0..i])      +
//	            (#hap0-agreements in [i+1..end])  -
//	            (#hap1-agreements in [i+1..end]) * FLIP_PENALTY
//
// and tracks the best score across both flip directions. The read is
// chimeric if the best score exceeds c[majority] + FLIP_THRES AND
// c[minority] + FLIP_THRES — i.e. the flip improves on simply
// assigning the whole read to either haplotype by at least
// FLIP_THRES votes on both sides.
func hasChimeraFlipPoint(alleles []int8, hap0AtHet []int8) bool {
	nHets := len(alleles)
	// Compute per-het hap0/hap1 vote contributions in the read's own
	// frame ("hap0" = phased hap0, "hap1" = phased hap1).
	v := make([]int8, nHets) // -1 ignore, 0 supports hap0, 1 supports hap1
	c0, c1 := 0, 0
	for i := 0; i < nHets; i++ {
		a := alleles[i]
		h0 := hap0AtHet[i]
		if a < 0 || h0 < 0 {
			v[i] = -1
			continue
		}
		if a == h0 {
			v[i] = 0
			c0++
		} else {
			v[i] = 1
			c1++
		}
	}
	// Prefix sums.
	type prefix struct{ p0, p1 int }
	left := make([]prefix, nHets)
	right := make([]prefix, nHets+1)
	{
		var p0, p1 int
		for i := 0; i < nHets; i++ {
			switch v[i] {
			case 0:
				p0++
			case 1:
				p1++
			}
			left[i] = prefix{p0, p1}
		}
	}
	{
		var p0, p1 int
		for i := nHets - 1; i >= 0; i-- {
			switch v[i] {
			case 0:
				p0++
			case 1:
				p1++
			}
			right[i] = prefix{p0, p1}
		}
	}
	bestScore := 0
	for i := 0; i < nHets-1; i++ {
		// Candidate 0: left half = hap0, right half = hap1.
		s0 := left[i].p0 + right[i+1].p1 - right[i+1].p0*flipPenalty
		// Candidate 1: left half = hap1, right half = hap0.
		s1 := left[i].p1 + right[i+1].p0 - right[i+1].p1*flipPenalty
		if s0 > bestScore {
			bestScore = s0
		}
		if s1 > bestScore {
			bestScore = s1
		}
	}
	// Matches upstream samtools phase.c fragphase() chimera-detection
	// criterion (`m - c[0] >= FLIP_THRES && m - c[1] >= FLIP_THRES`,
	// phase.c line 268, c. samtools 1.x).
	return bestScore-c0 >= flipThreshold && bestScore-c1 >= flipThreshold
}
