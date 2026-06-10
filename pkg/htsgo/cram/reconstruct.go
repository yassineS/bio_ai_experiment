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
// recPos is the record's 1-based reference alignment start (POS). A run
// of read positions not covered by any feature is a reference match:
// when an external reference is attached (rd.refBases non-nil) those
// bases are copied from the reference, and a substitution feature's read
// base is resolved through the substitution matrix; otherwise — the C4b
// reference-free path — such bases are filled with 'N' and the record is
// flagged as reference-needing.
func (rd *recordDecoder) reconstructMapped(feats []readFeature, readLen, recPos int32) (seq, qual []byte, cigar sam.Cigar, err error) {
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
	// refOffset is the count of reference bases consumed since recPos.
	// A read-and-reference operation (a match run or a substitution)
	// advances both cursors; a reference-only operation (a deletion or
	// reference skip) advances only refOffset; a read-only operation
	// (an insertion or soft clip) advances only readPos.
	refOffset := int32(0)

	for fi := range feats {
		f := &feats[fi]
		featStart := f.pos - 1
		if featStart < 0 || featStart > readLen {
			return nil, nil, nil, errFormat("read feature %d at in-read position %d is out of range (read length %d)",
				fi, f.pos, readLen)
		}
		// Every feature is positioned by its in-read coordinate. The run of
		// read positions between the read cursor and the feature is an
		// implicit reference match — emitted here, before the feature
		// itself, for EVERY feature kind (htslib's cram_decode.c does this
		// for deletions and reference skips too, not just base-consuming
		// features; emitting the gap only before consuming features
		// mis-orders e.g. "4M1D" into "1D" preceding the 4M). A feature
		// whose position equals the cursor (a quality annotation, or a
		// deletion right at the cursor) contributes no gap.
		if featStart > readPos {
			gap := featStart - readPos
			if cerr := rd.fillReferenceMatch(seq, covered, readPos, recPos+refOffset, gap); cerr != nil {
				return nil, nil, nil, wrapf(cerr, "read feature %d match run", fi)
			}
			ops.add(sam.CigarMatch, gap)
			readPos = featStart
			refOffset += gap
		} else if featStart < readPos && featConsumesRead(f.code) {
			return nil, nil, nil, errFormat("read feature %d at in-read position %d is behind the read cursor %d",
				fi, f.pos, readPos)
		}
		// A read-consuming feature is applied at the cursor (now equal to
		// featStart). A non-consuming positional feature — a quality
		// annotation — is applied at its own declared coordinate, which may
		// lie within bases already written, without moving the cursor.
		writePos := readPos
		if !featConsumesRead(f.code) {
			writePos = featStart
		}
		consumedRead, consumedRef, ferr := rd.applyFeature(f, seq, qual, covered, writePos, recPos+refOffset, &ops)
		if ferr != nil {
			return nil, nil, nil, wrapf(ferr, "read feature %d (%c)", fi, f.code)
		}
		readPos += consumedRead
		refOffset += consumedRef
	}
	// Read positions after the last feature are reference matches.
	if readPos < readLen {
		gap := readLen - readPos
		if cerr := rd.fillReferenceMatch(seq, covered, readPos, recPos+refOffset, gap); cerr != nil {
			return nil, nil, nil, wrapf(cerr, "trailing match run")
		}
		ops.add(sam.CigarMatch, gap)
		readPos = readLen
	}
	if readPos != readLen {
		return nil, nil, nil, errFormat("reconstructed read consumed %d of %d bases", readPos, readLen)
	}
	// Any base still uncovered would, in a reference-backed file, have
	// come from the external reference — fillReferenceMatch handled the
	// match runs, so a leftover here only occurs in the reference-free
	// fallback. Fill it with 'N' and report through needsReference.
	for i := int32(0); i < readLen; i++ {
		if !covered[i] {
			seq[i] = 'N'
			rd.needsReference = true
		}
	}
	return seq, qual, ops.cigar(), nil
}

// fillReferenceMatch writes a run of n reference-match bases into the
// read buffer starting at readPos, taking each base from the attached
// reference at the corresponding reference coordinate (1-based refPos
// for the first base). When no reference is attached it fills 'N' and
// sets needsReference, preserving the C4b fallback. It errors only when
// the reference span is too short for the run, which is a hard error:
// the reference does not cover the alignment the slice claims.
func (rd *recordDecoder) fillReferenceMatch(seq []byte, covered []bool, readPos, refPos, n int32) error {
	if n <= 0 {
		return nil
	}
	if int(readPos)+int(n) > len(seq) {
		return errFormat("match run of %d bases at read offset %d overruns the read (length %d)",
			n, readPos, len(seq))
	}
	if rd.refBases == nil {
		for i := int32(0); i < n; i++ {
			seq[readPos+i] = 'N'
			covered[readPos+i] = true
		}
		rd.needsReference = true
		return nil
	}
	idx := refPos - rd.refStart
	if idx < 0 || int(idx)+int(n) > len(rd.refBases) {
		return errFormat("match run needs reference bases %d-%d but the slice reference span covers %d-%d",
			refPos, refPos+n-1, rd.refStart, rd.refStart+int32(len(rd.refBases))-1)
	}
	copy(seq[readPos:readPos+n], rd.refBases[idx:idx+n])
	for i := int32(0); i < n; i++ {
		covered[readPos+i] = true
	}
	return nil
}

// referenceBaseAt returns the reference base at the 1-based reference
// coordinate refPos, or ok=false when no reference is attached or the
// coordinate falls outside the slice's resolved span.
func (rd *recordDecoder) referenceBaseAt(refPos int32) (byte, bool) {
	if rd.refBases == nil {
		return 0, false
	}
	idx := refPos - rd.refStart
	if idx < 0 || int(idx) >= len(rd.refBases) {
		return 0, false
	}
	return rd.refBases[idx], true
}

// applyFeature applies one read feature to the reconstruction buffers,
// writing any bases/qualities it carries at readPos and appending the
// CIGAR operation(s) it implies. refPos is the 1-based reference
// coordinate aligned with readPos. It returns how many read bases and
// how many reference bases the feature consumed, so the caller can
// advance its two cursors.
func (rd *recordDecoder) applyFeature(f *readFeature, seq, qual []byte, covered []bool, readPos, refPos int32, ops *cigarBuilder) (consumedRead, consumedRef int32, err error) {
	switch f.code {
	case featBases:
		n := int32(len(f.bases))
		if werr := writeBases(seq, covered, readPos, f.bases); werr != nil {
			return 0, 0, werr
		}
		ops.add(sam.CigarMatch, n)
		return n, n, nil
	case featScores:
		// A "q" feature carries quality scores for a stretch already
		// covered by bases; it consumes no read position itself.
		if int(readPos)+len(f.bases) > len(qual) {
			return 0, 0, errFormat("quality-score stretch overruns the read")
		}
		copy(qual[readPos:], f.bases)
		return 0, 0, nil
	case featBase:
		if werr := writeBases(seq, covered, readPos, []byte{f.base}); werr != nil {
			return 0, 0, werr
		}
		qual[readPos] = f.quality
		ops.add(sam.CigarMatch, 1)
		return 1, 1, nil
	case featQualityScore:
		if int(readPos) >= len(qual) {
			return 0, 0, errFormat("quality score overruns the read")
		}
		qual[readPos] = f.quality
		return 0, 0, nil
	case featSubst:
		// A substitution names the read base relative to the reference
		// base at this position via the substitution matrix. With an
		// external reference the read base is resolved exactly; without
		// one it is filled with 'N' and the record is flagged.
		base := byte('N')
		if refBase, ok := rd.referenceBaseAt(refPos); ok {
			base = rd.substMatrix.lookup(refBase, f.substCode)
		} else {
			rd.needsReference = true
		}
		if werr := writeBases(seq, covered, readPos, []byte{base}); werr != nil {
			return 0, 0, werr
		}
		ops.add(sam.CigarMatch, 1)
		return 1, 1, nil
	case featInsertion:
		n := int32(len(f.bases))
		if werr := writeBases(seq, covered, readPos, f.bases); werr != nil {
			return 0, 0, werr
		}
		ops.add(sam.CigarInsertion, n)
		return n, 0, nil
	case featInsertBase:
		if werr := writeBases(seq, covered, readPos, []byte{f.base}); werr != nil {
			return 0, 0, werr
		}
		ops.add(sam.CigarInsertion, 1)
		return 1, 0, nil
	case featSoftClip:
		n := int32(len(f.bases))
		if werr := writeBases(seq, covered, readPos, f.bases); werr != nil {
			return 0, 0, werr
		}
		ops.add(sam.CigarSoftClip, n)
		return n, 0, nil
	case featDeletion:
		if f.length < 0 {
			return 0, 0, errFormat("deletion feature has negative length %d", f.length)
		}
		ops.add(sam.CigarDeletion, f.length)
		return 0, f.length, nil
	case featRefSkip:
		if f.length < 0 {
			return 0, 0, errFormat("reference-skip feature has negative length %d", f.length)
		}
		ops.add(sam.CigarSkipped, f.length)
		return 0, f.length, nil
	case featPadding:
		if f.length < 0 {
			return 0, 0, errFormat("padding feature has negative length %d", f.length)
		}
		ops.add(sam.CigarPadding, f.length)
		return 0, 0, nil
	case featHardClip:
		if f.length < 0 {
			return 0, 0, errFormat("hard-clip feature has negative length %d", f.length)
		}
		ops.add(sam.CigarHardClip, f.length)
		return 0, 0, nil
	default:
		return 0, 0, errFormat("unknown read-feature code %#02x", f.code)
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
