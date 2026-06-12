// Package mdnm computes the SAM MD:Z and NM:i auxiliary tags by walking a
// record's CIGAR against a reference sequence. It is the single, shared
// implementation of the MD/NM algorithm used by both the samtools calmd
// subcommand and the CRAM decoder's reference-derived tag regeneration, so
// the two paths cannot drift apart.
//
// The walk mirrors upstream samtools' bam_fillmd1_core in
// reference_code/samtools/bam_md.c and htslib's equivalent in
// reference_code/htslib/cram/cram_decode.c:
//
//   - For every reference-consuming CIGAR op (M / = / X) the read base is
//     compared against the reference base. Matches extend a running
//     match-run counter. A mismatch flushes the counter as a decimal run
//     length, appends the reference base (uppercased), resets the counter,
//     and increments NM.
//   - A deletion ('D') flushes the counter, then appends '^' followed by the
//     deleted reference bases and adds the deletion length to NM.
//   - An insertion ('I') adds its length to NM and consumes only the read.
//     A soft clip ('S') consumes only the read. A reference skip ('N')
//     advances only the reference. Hard clip ('H') and padding ('P')
//     consume neither for MD/NM purposes.
//   - The final match-run counter is flushed at the end of the CIGAR, so the
//     MD string always ends with a (possibly zero) match-run length.
package mdnm

import (
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Compute walks rec's CIGAR against ref and returns the MD:Z string and the
// NM:i count. ref holds reference bases and refOffset is the 0-based
// reference coordinate of ref[0], so a slice-local reference span can be
// passed without materialising a whole chromosome: pass refOffset=0 with a
// whole-contig ref (rec.Pos indexes it directly), or refOffset=refStart-1
// with a span beginning at the 1-based coordinate refStart.
//
// When the CIGAR reaches a position the reference span does not cover the
// walk stops early, flushing the match-run as upstream's break semantics do;
// the MD and NM computed up to that point are returned. Callers that require
// the reference to fully cover the alignment should ensure ref spans
// [rec.Pos, rec.Pos+referenceLength) before calling.
func Compute(rec *sam.Record, ref []byte, refOffset int) (md string, nm int) {
	var (
		mdBuf   []byte
		matched int
		qpos    int
		rpos    = int(rec.Pos) - 1 - refOffset // index into ref of the alignment start.
	)
	seq := []byte(rec.Seq)

	flushRun := func() {
		mdBuf = strconv.AppendInt(mdBuf, int64(matched), 10)
		matched = 0
	}

	for _, op := range rec.Cigar {
		oplen := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			truncated := false
			for j := 0; j < oplen; j++ {
				if rpos+j < 0 || rpos+j >= len(ref) || qpos+j >= len(seq) {
					// Out of bounds — upstream "break" semantics: leave the
					// rest of the read alone. The MD up to here is still
					// emitted in the final flush.
					truncated = true
					break
				}
				if baseMatches(seq[qpos+j], ref[rpos+j]) {
					matched++
				} else {
					flushRun()
					mdBuf = append(mdBuf, upperByte(ref[rpos+j]))
					nm++
				}
			}
			if truncated {
				flushRun()
				return string(mdBuf), nm
			}
			rpos += oplen
			qpos += oplen
		case sam.CigarDeletion:
			flushRun()
			mdBuf = append(mdBuf, '^')
			actual := 0
			for j := 0; j < oplen; j++ {
				if rpos+j < 0 || rpos+j >= len(ref) {
					break
				}
				mdBuf = append(mdBuf, upperByte(ref[rpos+j]))
				actual++
			}
			nm += actual
			rpos += actual
			if actual < oplen {
				// The deletion could not be fully consumed from the
				// reference; upstream breaks out of the CIGAR loop, and the
				// trailing flush still appends the "0" terminator.
				flushRun()
				return string(mdBuf), nm
			}
		case sam.CigarInsertion:
			nm += oplen
			qpos += oplen
		case sam.CigarSoftClip:
			qpos += oplen
		case sam.CigarSkipped:
			rpos += oplen
		case sam.CigarHardClip, sam.CigarPadding:
			// Neither consumes query nor reference for MD/NM purposes.
		}
	}
	flushRun()
	return string(mdBuf), nm
}

// upperByte folds an ASCII byte to uppercase. Non-letters pass through.
func upperByte(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - ('a' - 'A')
	}
	return b
}

// baseMatches reports whether the read base readByte matches the reference
// base refByte for MD/NM purposes. It mirrors upstream bam_md.c's test
// `(c1==c2 && c1!=15 && c2!=15) || c1==0`: equal non-N bases match, an N on
// either side is a mismatch, and a read base of '=' (BAM 4-bit code 0)
// matches anything.
func baseMatches(readByte, refByte byte) bool {
	readBase := upperByte(readByte)
	if readBase == '=' {
		return true
	}
	refBase := upperByte(refByte)
	return readBase == refBase && readBase != 'N' && refBase != 'N'
}
