package samtools

import (
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FixmateOptions configures Fixmate. The defaults match `samtools fixmate`'s
// DEFAULT behaviour: fix RNEXT/PNEXT/TLEN/0x8 flags AND add the MC (mate CIGAR)
// and MQ (mate MAPQ) aux tags. Upstream (bam_mate.c: sync_mate -> sync_mq_mc)
// writes MC/MQ unconditionally on every pair, so this port does too. Optional
// knobs add the `ms` mate-score tag and remove unmapped reads.
type FixmateOptions struct {
	// AddMateScore (-m) writes the `ms` aux tag (sum of base qualities of
	// the mate's bases >= Q15). MQ (mate MAPQ) is written by default (upstream
	// adds MQ unconditionally in sync_mq_mc, not only under -m).
	AddMateScore bool
	// AddMateCigar (-c) is retained for CLI compatibility. Upstream's real -c
	// adds the CT template-cigar tag; the MC (mate CIGAR) tag it is commonly
	// associated with is now emitted BY DEFAULT (upstream writes MC in
	// sync_mq_mc regardless of any flag), so this field no longer gates MC.
	AddMateCigar bool
	// RemoveUnmapped (-r) drops records where both this read and its
	// mate are unmapped (and unpaired entirely-unmapped singletons).
	RemoveUnmapped bool
	// NoPG suppresses @PG injection. v1 never injects @PG so this is a
	// no-op kept for flag-compat.
	NoPG bool
	// Threads is accepted for upstream-CLI compatibility; ignored.
	Threads int
}

// Fixmate walks the input as consecutive pairs (the upstream contract:
// the input must be name-sorted or name-collated so the two mates of a
// pair sit adjacent). For each pair we compute and patch:
//   - RNext: mate's RNAME (or "=" for same chromosome)
//   - PNext: mate's POS (1-based)
//   - TLen: signed insert size, leftmost mate negative, rightmost positive
//   - 0x8 flag (mate unmapped) toggled to match reality
//   - 0x20 flag (mate reverse) toggled to match reality
//
// When AddMateCigar/AddMateScore are set, the corresponding MQ/MC/ms aux
// tags are added too. Singleton (unpaired-by-flag) records pass through
// unchanged. The implementation is single-pass and streams.
func Fixmate(in io.Reader, out io.Writer, opts FixmateOptions) error {
	br, err := sam.NewReader(in)
	if err != nil {
		return err
	}
	hdr := br.Header()
	bw := sam.NewBAMWriter(out)
	if err := bw.WriteHeader(hdr); err != nil {
		return err
	}

	var prev *sam.Record
	emit := func(r *sam.Record) error {
		if r == nil {
			return nil
		}
		if opts.RemoveUnmapped && fixmateShouldDrop(r) {
			return nil
		}
		return bw.Write(r)
	}

	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if prev == nil {
			prev = rec
			continue
		}
		if prev.QName == rec.QName && prev.IsPaired() && rec.IsPaired() {
			fixPair(prev, rec, opts)
			if err := emit(prev); err != nil {
				return err
			}
			if err := emit(rec); err != nil {
				return err
			}
			prev = nil
			continue
		}
		// prev was a singleton.
		if err := emit(prev); err != nil {
			return err
		}
		prev = rec
	}
	if err := emit(prev); err != nil {
		return err
	}
	return bw.Close()
}

// fixmateShouldDrop reports whether a record should be removed under -r.
// Upstream's policy: drop unmapped reads (we implement the standard
// "FLAG & 0x4 set" form).
func fixmateShouldDrop(r *sam.Record) bool {
	return r.IsUnmapped()
}

// fixPair patches the two records of a name-pair to make their
// mate-related fields consistent. Both records' QName already match.
func fixPair(a, b *sam.Record, opts FixmateOptions) {
	// Mate-related flag fixups.
	syncMateFlags(a, b)
	syncMateFlags(b, a)

	// RNext / PNext / TLen.
	a.RNext = mateRName(a.RName, b.RName)
	a.PNext = b.Pos
	b.RNext = mateRName(b.RName, a.RName)
	b.PNext = a.Pos
	tlen := computeTLen(a, b)
	a.TLen = tlen
	b.TLen = -tlen

	// MQ (mate MAPQ) and MC (mate CIGAR) are written BY DEFAULT, matching
	// upstream sync_mq_mc (called unconditionally from sync_mate for every
	// pair). The gating mirrors bam_mate.c exactly:
	//   - MQ is appended to dest only when the SOURCE (mate) is mapped.
	//   - MC is appended to dest when either the source OR dest is mapped.
	syncMQMC(a, b)
	syncMQMC(b, a)

	if opts.AddMateScore {
		setAuxInt(a, "ms", int64(mateScore(b)))
		setAuxInt(b, "ms", int64(mateScore(a)))
	}
}

// syncMQMC writes dest's MQ (mate MAPQ) and MC (mate CIGAR) aux tags from src,
// matching upstream bam_mate.c sync_mq_mc EXACTLY:
//   - MQ: copied from src->core.qual, but only when src is mapped.
//   - MC: src's CIGAR string, added when either src OR dest is mapped
//     (an all-unmapped pair gets no MC, as upstream leaves it out).
//
// Existing MQ/MC tags on dest are replaced (setAuxInt/setAuxString overwrite).
// An unmapped record's CIGAR renders as "*" (Cigar.String()), matching how
// bam_format_cigar emits a 0-length CIGAR for an unmapped read.
func syncMQMC(src, dest *sam.Record) {
	srcMapped := !src.IsUnmapped()
	destMapped := !dest.IsUnmapped()
	if srcMapped {
		setAuxInt(dest, "MQ", int64(src.MapQ))
	}
	if srcMapped || destMapped {
		setAuxString(dest, "MC", src.Cigar.String())
	}
}

// syncMateFlags copies the mate-state bits from src into dst's mate
// flags. Specifically: dst's 0x8 mirrors src's 0x4, and dst's 0x20
// mirrors src's 0x10.
func syncMateFlags(dst, src *sam.Record) {
	if src.IsUnmapped() {
		dst.Flag |= sam.FlagMateUnmapped
	} else {
		dst.Flag &^= sam.FlagMateUnmapped
	}
	if src.Flag&sam.FlagReverse != 0 {
		dst.Flag |= sam.FlagMateReverse
	} else {
		dst.Flag &^= sam.FlagMateReverse
	}
}

// mateRName returns "=" when both mates align to the same reference,
// otherwise it returns the mate's reference name verbatim.
func mateRName(self, mate string) string {
	if mate == "" || mate == "*" {
		return "*"
	}
	if mate == self && self != "" && self != "*" {
		return "="
	}
	return mate
}

// computeTLen returns the signed template length of a as part of a
// pair (a is the leftmost mate when its 5' position is smaller).
//
// We follow the SAM spec and upstream samtools: TLEN is unsigned distance
// between the leftmost mapped 5'-end and the rightmost mapped 5'-end,
// inclusive of both endpoints; sign is positive on the leftmost record
// and negative on the rightmost. Both records must be mapped on the same
// reference for a non-zero TLen — otherwise TLen is 0.
func computeTLen(a, b *sam.Record) int64 {
	if a.IsUnmapped() || b.IsUnmapped() {
		return 0
	}
	if a.RName == "" || a.RName != b.RName {
		return 0
	}
	aBeg, aEnd := a.Pos, a.EndPosition()
	bBeg, bEnd := b.Pos, b.EndPosition()
	left, right := aBeg, bEnd
	if bBeg < aBeg {
		left = bBeg
		right = aEnd
	}
	tlen := right - left + 1
	// Sign: a is leftmost ⇒ positive on a.
	if aBeg <= bBeg {
		return tlen
	}
	return -tlen
}

// mateScore is the sum of every base quality >= 15 in the read,
// matching upstream's `bam_aux_get_str("ms")` definition.
func mateScore(r *sam.Record) int {
	score := 0
	for _, q := range r.Qual {
		if q >= 15 {
			score += int(q)
		}
	}
	return score
}

// setAuxString sets or replaces a 'Z'-typed aux tag.
func setAuxString(r *sam.Record, tag, value string) {
	for i, a := range r.Aux {
		if a.Tag == tag {
			r.Aux[i].Type = 'Z'
			r.Aux[i].Value = value
			return
		}
	}
	r.Aux = append(r.Aux, sam.Aux{Tag: tag, Type: 'Z', Value: value})
}

// setAuxInt sets or replaces an integer-typed aux tag.
func setAuxInt(r *sam.Record, tag string, value int64) {
	for i, a := range r.Aux {
		if a.Tag == tag {
			r.Aux[i].Type = 'i'
			r.Aux[i].Value = value
			return
		}
	}
	r.Aux = append(r.Aux, sam.Aux{Tag: tag, Type: 'i', Value: value})
}
