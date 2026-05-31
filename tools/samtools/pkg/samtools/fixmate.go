package samtools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FixmateOptions configures Fixmate. The defaults match `samtools fixmate`'s
// minimal mode: just fix RNEXT/PNEXT/TLEN/0x8 flags. Optional knobs add
// MQ/MC/ms tags and remove unmapped reads.
type FixmateOptions struct {
	// AddMateScore (-m) writes the `ms` aux tag (sum of base qualities of
	// the mate's bases >= Q15).
	AddMateScore bool
	// AddMateCigar (-c) writes the `MC` aux tag (mate's CIGAR) and the
	// `MQ` aux tag (mate's MAPQ).
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

	// Upstream's `sync_mq_mc` always copies MQ (mate MAPQ) and MC
	// (mate CIGAR) onto both reads of a pair when at least one is
	// mapped, regardless of CLI flags. We mirror that here so the
	// default fixmate invocation produces upstream-compatible
	// output.
	if !a.IsUnmapped() {
		setAuxInt(b, "MQ", int64(a.MapQ))
	}
	if !b.IsUnmapped() {
		setAuxInt(a, "MQ", int64(b.MapQ))
	}
	if !a.IsUnmapped() || !b.IsUnmapped() {
		setAuxString(a, "MC", b.Cigar.String())
		setAuxString(b, "MC", a.Cigar.String())
	}
	if opts.AddMateScore {
		setAuxInt(a, "ms", int64(mateScore(b)))
		setAuxInt(b, "ms", int64(mateScore(a)))
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
func computeTLen(a, b *sam.Record) int32 {
	if a.IsUnmapped() || b.IsUnmapped() {
		return 0
	}
	if a.RName == "" || a.RName != b.RName {
		return 0
	}
	aBeg, aEnd := int32(a.Pos), a.EndPosition()
	bBeg, bEnd := int32(b.Pos), b.EndPosition()
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

// describe is unused but kept here as documentation of the BAM record
// fields a fixmate run touches: it's a single-line cheat sheet for
// readers of this file.
//
//nolint:unused
func describe() string {
	return fmt.Sprintf("Fixmate touches: Flag(0x8,0x20), RNext, PNext, TLen, +Aux(MQ,MC,ms)")
}
