package cram

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// encodeFeatures encodes a mapped record's read features. Reference-free
// CRAM carries the whole sequence in the features, so every CIGAR
// operation becomes one feature:
//
//   - M / = / X  → a base-stretch feature ('b'), carrying the literal
//     run of read bases (the BB data series);
//   - I          → an insertion feature ('I'), carrying the inserted
//     bases (the IN data series);
//   - S          → a soft-clip feature ('S'), carrying the clipped bases
//     (the SC data series);
//   - D          → a deletion feature ('D'), carrying the length (DL);
//   - N          → a reference-skip feature ('N'), carrying the length (RS);
//   - P          → a padding feature ('P'), carrying the length (PD);
//   - H          → a hard-clip feature ('H'), carrying the length (HC).
//
// This mirrors reconstructMapped exactly: the base-carrying features
// cover every read position so no reference-derived 'N' fill is needed,
// and the CIGAR is rebuilt from the features' implied operations.
//
// M, = and X are encoded identically — each becomes a base-stretch
// feature carrying the literal read bases. This matches htslib's CRAM
// encoder, whose cigar loop (cram_encode.c) folds BAM_CMATCH,
// BAM_CBASE_MATCH ('=') and BAM_CBASE_MISMATCH ('X') into a single case
// and deliberately "doesn't trust = and X to be correct", comparing every
// base against the reference instead. In the reference-free encoding this
// writer produces that per-base comparison reduces to "carry the bases
// verbatim", so an aligner run with --eqx (which emits =/X instead of M)
// round-trips to the same SAM record samtools decodes from its own
// reference-based CRAM (where =/X collapse to M plus MD/NM tags).
//
// The CIGAR back-step op (B, sam.CigarBack, op 9) is deliberately not
// supported: htslib's own CRAM encoder rejects it with "Unknown CIGAR op
// code" (the cigar switch in cram_encode.c has no BAM_CBACK case, and the
// consensus builder lists BACK in its cig_skip table), and htslib's SAM
// parser will not even construct a record carrying B without a
// length-mismatch error. Matching upstream, the writer returns a clear
// error rather than inventing a feature encoding htslib could never read
// back. See docs/CRAM_ROADMAP.md §6 (the lossy-names + =/X/B entry).
// referenceSpan returns the number of reference bases a CIGAR consumes — the
// sum of the M/=/X (match), D (deletion) and N (skip) operation lengths. It is
// used to bounds-check a read against its contig reference before choosing the
// reference-based encoding path.
func referenceSpan(cigar sam.Cigar) int {
	span := 0
	for _, op := range cigar {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch,
			sam.CigarDeletion, sam.CigarSkipped:
			span += int(op.Length())
		}
	}
	return span
}

func (e *recordEncoder) encodeFeatures(rec *sam.Record, readLen int) error {
	b := e.buffers
	seq := seqBytes(rec.Seq)

	// Reference-based encoding: when the reference for this record's contig is
	// available and the read's reference span lies within it, a match run is
	// diffed against the reference and only its mismatches are stored, as
	// substitution features (the matched bases are reconstructed from the
	// reference on decode — reconstructMapped fills inter-feature gaps as
	// reference matches). This is what makes CRAM small. ref is nil — keeping
	// the self-contained base-stretch encoding — when no reference was supplied
	// or the read falls outside it.
	var ref []byte
	if e.reference != nil && rec.Pos > 0 {
		if r, ok := e.reference[rec.RName]; ok {
			start := int(rec.Pos) - 1
			if start >= 0 && start+referenceSpan(rec.Cigar) <= len(r) {
				ref = r
				// Record that this container used the external reference, so
				// its compression header omits RR (reference required) and the
				// decoder loads the reference to fill the implicit match runs.
				e.usedReference = true
			}
		}
	}

	type feature struct {
		code      byte
		pos       int32 // 1-based in-read position.
		bases     []byte
		length    int32
		substCode byte // for featSubst: the BS substitution code.
	}
	var feats []feature

	readPos := int32(0) // 0-based cursor within the read.
	refPos := int32(0)  // 0-based cursor within ref (only used when ref != nil).
	if ref != nil {
		refPos = rec.Pos - 1
	}
	for _, op := range rec.Cigar {
		n := int32(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if int(readPos)+int(n) > len(seq) {
				return fmt.Errorf("CIGAR match run overruns SEQ")
			}
			if ref != nil {
				// Emit a substitution feature only where the read base differs
				// from the reference; equal positions are left implicit and
				// reconstructed from the reference. Equality is an exact byte
				// compare so the decoder's verbatim reference copy always
				// reproduces the read base (a soft-masked/lower-case reference
				// base never counts as a match).
				for k := int32(0); k < n; k++ {
					rb := ref[refPos+k]
					qb := seq[readPos+k]
					if qb != rb {
						feats = append(feats, feature{
							code:      featSubst,
							pos:       readPos + k + 1,
							substCode: substCodeFor(rb, qb),
						})
					}
				}
				refPos += n
			} else {
				feats = append(feats, feature{
					code:  featBases,
					pos:   readPos + 1,
					bases: append([]byte(nil), seq[readPos:readPos+n]...),
				})
			}
			readPos += n
		case sam.CigarInsertion:
			if int(readPos)+int(n) > len(seq) {
				return fmt.Errorf("CIGAR insertion overruns SEQ")
			}
			feats = append(feats, feature{
				code:  featInsertion,
				pos:   readPos + 1,
				bases: append([]byte(nil), seq[readPos:readPos+n]...),
			})
			readPos += n
		case sam.CigarSoftClip:
			if int(readPos)+int(n) > len(seq) {
				return fmt.Errorf("CIGAR soft clip overruns SEQ")
			}
			feats = append(feats, feature{
				code:  featSoftClip,
				pos:   readPos + 1,
				bases: append([]byte(nil), seq[readPos:readPos+n]...),
			})
			readPos += n
		case sam.CigarDeletion:
			feats = append(feats, feature{code: featDeletion, pos: readPos + 1, length: n})
			refPos += n
		case sam.CigarSkipped:
			feats = append(feats, feature{code: featRefSkip, pos: readPos + 1, length: n})
			refPos += n
		case sam.CigarPadding:
			feats = append(feats, feature{code: featPadding, pos: readPos + 1, length: n})
		case sam.CigarHardClip:
			feats = append(feats, feature{code: featHardClip, pos: readPos + 1, length: n})
		case sam.CigarBack:
			// htslib's CRAM encoder has no back-step case and rejects B with
			// "Unknown CIGAR op code"; its SAM parser cannot even build such a
			// record. We match that and reject B rather than emit a feature
			// stream upstream could never decode. See docs/CRAM_ROADMAP.md §6.
			return fmt.Errorf("CIGAR back-step op B is not supported (htslib's CRAM encoder rejects it as well)")
		default:
			return fmt.Errorf("CIGAR operation %c is not supported by the simple CRAM writer", op.Char())
		}
	}

	// A mapped record with no CIGAR but a non-empty sequence is encoded
	// as a single base-stretch covering the whole read; the reader's
	// reconstruction yields a plain match CIGAR for it.
	if len(rec.Cigar) == 0 && len(seq) > 0 {
		feats = append(feats, feature{
			code:  featBases,
			pos:   1,
			bases: append([]byte(nil), seq...),
		})
		readPos = int32(len(seq))
	}

	if readPos != int32(readLen) {
		return fmt.Errorf("CIGAR consumed %d read bases but SEQ has %d", readPos, readLen)
	}

	b.fn = e.putU(b.fn, int32(len(feats)))
	var prevPos int32
	for _, f := range feats {
		b.fc = append(b.fc, f.code)
		b.fp = e.putU(b.fp, f.pos-prevPos)
		prevPos = f.pos
		switch f.code {
		case featSubst:
			b.bs = append(b.bs, f.substCode)
		case featBases:
			b.bbLen = e.putU(b.bbLen, int32(len(f.bases)))
			b.bb = append(b.bb, f.bases...)
		case featInsertion:
			b.inLen = e.putU(b.inLen, int32(len(f.bases)))
			b.in = append(b.in, f.bases...)
		case featSoftClip:
			b.scLen = e.putU(b.scLen, int32(len(f.bases)))
			b.sc = append(b.sc, f.bases...)
		case featDeletion:
			b.dl = e.putU(b.dl, f.length)
		case featRefSkip:
			b.rs = e.putU(b.rs, f.length)
		case featPadding:
			b.pd = e.putU(b.pd, f.length)
		case featHardClip:
			b.hc = e.putU(b.hc, f.length)
		}
	}
	return nil
}
