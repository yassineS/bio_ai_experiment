package cram

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// reconstructMapped turns a mapped record's read-feature list into its
// SEQ bytes, per-base QUAL scores and CIGAR. The walk advances two
// cursors in lockstep — a position in the read and a position on the
// reference — applying each feature at its in-read coordinate and
// emitting CIGAR operations for the runs between and within features.
//
// In a reference-free CRAM every base of the read is supplied by a
// feature (a "b" stretch, an "i"/"B" single base, or a soft clip), so no
// external reference is consulted. A run of read positions not covered
// by any feature is treated as a reference match: in a reference-backed
// file those bases would be copied from the reference, and this decoder
// leaves a placeholder so the caller can see which bases need one.
func (rd *recordDecoder) reconstructMapped(feats []readFeature, readLen int32) (seq, qual []byte, cigar sam.Cigar, err error) {
	if readLen < 0 {
		return nil, nil, nil, errFormat("record declares a negative read length %d", readLen)
	}
	seq = make([]byte, readLen)
	qual = make([]byte, readLen)
	for i := range qual {
		qual[i] = 0xff // 0xff is the SAM "no quality" sentinel.
	}
	covered := make([]bool, readLen)

	var ops cigarBuilder
	readPos := int32(0) // 0-based cursor within the read.

	for fi := range feats {
		f := &feats[fi]
		// f.pos is 1-based. A feature that consumes read bases starts a
		// new read position; the gap before it is an implicit run of
		// matches. A non-consuming feature (a quality score, a deletion,
		// a reference skip, padding or a hard clip) carries no read
		// bases and may sit at a position the read cursor has already
		// passed — for example a Q feature setting the quality of a base
		// an earlier insertion supplied — so it is applied where it is
		// asked without disturbing the cursor.
		featStart := f.pos - 1
		if featStart < 0 || featStart > readLen {
			return nil, nil, nil, errFormat("read feature %d at in-read position %d is out of range (read length %d)",
				fi, f.pos, readLen)
		}
		if featConsumesRead(f.code) {
			if featStart < readPos {
				return nil, nil, nil, errFormat("read feature %d at in-read position %d is behind the read cursor %d",
					fi, f.pos, readPos)
			}
			if featStart > readPos {
				ops.add(sam.CigarMatch, featStart-readPos)
				readPos = featStart
			}
		}
		consumed, ferr := rd.applyFeature(f, seq, qual, covered, featStart, &ops)
		if ferr != nil {
			return nil, nil, nil, wrapf(ferr, "read feature %d (%c)", fi, f.code)
		}
		readPos += consumed
	}
	// Read positions after the last feature are reference matches.
	if readPos < readLen {
		ops.add(sam.CigarMatch, readLen-readPos)
		readPos = readLen
	}
	if readPos != readLen {
		return nil, nil, nil, errFormat("reconstructed read consumed %d of %d bases", readPos, readLen)
	}
	// Every base must have been supplied by a feature for a reference-free
	// file. A base no feature covered would, in a reference-backed file,
	// be copied from the external reference; this decoder has no
	// reference, so it fills such bases with 'N' and reports — through
	// the caller's needsReference accounting — that the record is only
	// partially reconstructed.
	for i := int32(0); i < readLen; i++ {
		if !covered[i] {
			seq[i] = 'N'
			rd.needsReference = true
		}
	}
	return seq, qual, ops.cigar(), nil
}

// applyFeature applies one read feature to the reconstruction buffers,
// writing any bases/qualities it carries at readPos and appending the
// CIGAR operation(s) it implies. It returns how many read bases the
// feature consumed (so the caller can advance its read cursor).
func (rd *recordDecoder) applyFeature(f *readFeature, seq, qual []byte, covered []bool, readPos int32, ops *cigarBuilder) (int32, error) {
	switch f.code {
	case featBases:
		n := int32(len(f.bases))
		if err := writeBases(seq, covered, readPos, f.bases); err != nil {
			return 0, err
		}
		ops.add(sam.CigarMatch, n)
		return n, nil
	case featScores:
		// A "q" feature carries quality scores for a stretch already
		// covered by bases; it consumes no read position itself.
		if int(readPos)+len(f.bases) > len(qual) {
			return 0, errFormat("quality-score stretch overruns the read")
		}
		copy(qual[readPos:], f.bases)
		return 0, nil
	case featBase:
		if err := writeBases(seq, covered, readPos, []byte{f.base}); err != nil {
			return 0, err
		}
		qual[readPos] = f.quality
		ops.add(sam.CigarMatch, 1)
		return 1, nil
	case featQualityScore:
		if int(readPos) >= len(qual) {
			return 0, errFormat("quality score overruns the read")
		}
		qual[readPos] = f.quality
		return 0, nil
	case featSubst:
		// A substitution names a single read base relative to the
		// reference base via the substitution matrix. Without an external
		// reference the read base cannot be resolved, so the position is
		// filled with 'N' and the record is flagged as reference-needing.
		if err := writeBases(seq, covered, readPos, []byte{'N'}); err != nil {
			return 0, err
		}
		rd.needsReference = true
		ops.add(sam.CigarMatch, 1)
		return 1, nil
	case featInsertion:
		n := int32(len(f.bases))
		if err := writeBases(seq, covered, readPos, f.bases); err != nil {
			return 0, err
		}
		ops.add(sam.CigarInsertion, n)
		return n, nil
	case featInsertBase:
		if err := writeBases(seq, covered, readPos, []byte{f.base}); err != nil {
			return 0, err
		}
		ops.add(sam.CigarInsertion, 1)
		return 1, nil
	case featSoftClip:
		n := int32(len(f.bases))
		if err := writeBases(seq, covered, readPos, f.bases); err != nil {
			return 0, err
		}
		ops.add(sam.CigarSoftClip, n)
		return n, nil
	case featDeletion:
		if f.length < 0 {
			return 0, errFormat("deletion feature has negative length %d", f.length)
		}
		ops.add(sam.CigarDeletion, f.length)
		return 0, nil
	case featRefSkip:
		if f.length < 0 {
			return 0, errFormat("reference-skip feature has negative length %d", f.length)
		}
		ops.add(sam.CigarSkipped, f.length)
		return 0, nil
	case featPadding:
		if f.length < 0 {
			return 0, errFormat("padding feature has negative length %d", f.length)
		}
		ops.add(sam.CigarPadding, f.length)
		return 0, nil
	case featHardClip:
		if f.length < 0 {
			return 0, errFormat("hard-clip feature has negative length %d", f.length)
		}
		ops.add(sam.CigarHardClip, f.length)
		return 0, nil
	default:
		return 0, errFormat("unknown read-feature code %#02x", f.code)
	}
}

// writeBases copies a feature's bases into the read buffer at readPos,
// marking those positions covered. It errors if the bases overrun the
// read, which only a malformed file produces.
func writeBases(seq []byte, covered []bool, readPos int32, bases []byte) error {
	if int(readPos)+len(bases) > len(seq) {
		return errFormat("feature bases overrun the read (%d bases at offset %d, read length %d)",
			len(bases), readPos, len(seq))
	}
	copy(seq[readPos:], bases)
	for i := range bases {
		covered[int(readPos)+i] = true
	}
	return nil
}

// featConsumesRead reports whether a read feature occupies read
// positions of its own — true for the base-carrying features (a base
// stretch, a single base, a substitution, an insertion or a soft clip)
// and false for features that annotate or skip without contributing
// read bases (quality scores, deletions, reference skips, padding and
// hard clips).
func featConsumesRead(code byte) bool {
	switch code {
	case featBases, featBase, featSubst, featInsertion, featInsertBase, featSoftClip:
		return true
	default:
		return false
	}
}

// cigarBuilder accumulates CIGAR operations, coalescing a run of the
// same operation into one CigarOp so the emitted CIGAR is canonical.
type cigarBuilder struct {
	ops sam.Cigar
}

// add appends one CIGAR operation of the given op code and length,
// merging it with the previous operation when the codes match. A
// zero-length operation is dropped.
func (b *cigarBuilder) add(op uint32, length int32) {
	if length <= 0 {
		return
	}
	if n := len(b.ops); n > 0 && b.ops[n-1].Op() == op {
		merged := b.ops[n-1].Length() + uint32(length)
		b.ops[n-1] = sam.CigarOp(merged<<4 | op)
		return
	}
	b.ops = append(b.ops, sam.CigarOp(uint32(length)<<4|op))
}

// cigar returns the accumulated CIGAR.
func (b *cigarBuilder) cigar() sam.Cigar { return b.ops }
