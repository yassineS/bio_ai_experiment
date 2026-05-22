// Package samtools — calmd subcommand.
//
// Calmd computes the MD and NM auxiliary tags by walking each record's CIGAR
// against the reference. The implementation mirrors upstream samtools
// `bam_fillmd1_core` in `reference_code/samtools/bam_md.c`:
//
//   - For every reference-consuming CIGAR op (M / = / X) we compare the
//     read base against the reference base. Matches add to the running
//     match-run counter. Mismatches flush the counter, append the reference
//     base (uppercase), reset the counter, and increment NM.
//   - Deletions ('D') flush the counter, then emit ^ followed by the deleted
//     reference bases and add the deletion length to NM.
//   - Insertions ('I') and soft-clip ('S') only consume query; insertions
//     add their length to NM. Reference-skip ('N') and pad ('P') do nothing
//     to MD/NM beyond moving the reference pointer.
//   - The final match-run counter is flushed at the end of the CIGAR.
//
// Unmapped records pass through unchanged. Records with seq "*" are passed
// through unchanged (matching upstream which prints a warning and skips
// them). Existing NM / MD aux tags are overwritten (with a stderr warning
// if the new value differs, mirroring upstream's "different NM/MD" line).
package samtools

import (
	"fmt"
	"io"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// CalmdOptions configures Calmd. The defaults (zero value) produce the same
// behaviour as a plain `samtools calmd in.bam ref.fa` invocation: update both
// MD and NM, emit SAM text on stdout, leave matched bases intact.
type CalmdOptions struct {
	// UseEqual rewrites match bases in the read sequence to '=' (the -e
	// flag). When set, NM and MD are still computed against the original
	// read base, only the emitted SEQ is rewritten.
	UseEqual bool
	// OutputBAM emits BAM (the -b flag). Implies BGZF output.
	OutputBAM bool
	// Uncompressed forces compress-level 0 on BAM output (the -u flag,
	// which also implies -b in upstream samtools).
	Uncompressed bool
	// ExtendedBAQ controls upstream's -E "extended-BAQ" mode: BAQ values may
	// exceed the input base qualities, trading specificity for sensitivity.
	// Only takes effect together with RealignBAQ.
	ExtendedBAQ bool
	// AdjustCapQ holds upstream's -A flag. With RealignBAQ it switches the
	// BAQ pass from "compute the BQ:Z: tag" to "cap base qualities by BAQ"
	// (writing a ZQ:Z: tag instead). Without RealignBAQ it has no effect.
	AdjustCapQ bool
	// RealignBAQ holds upstream's -r flag. When set, BAQ realignment runs
	// before the MD/NM fill: each record gets a BQ:Z: tag (or, with
	// AdjustCapQ, capped qualities plus a ZQ:Z: tag).
	RealignBAQ bool
	// CapMapQ holds upstream's -C flag. When greater than 10, each record's
	// MAPQ is capped by baq.SamCapMapq using CapMapQ as the threshold,
	// matching bam_md.c's `if (capQ > 10)` gate.
	CapMapQ int
	// Quiet suppresses the per-record "different NM/MD" stderr line when an
	// existing tag is overwritten with a different value.
	Quiet bool
	// DropTags drops every aux tag except RG (the -d flag). Upstream applies
	// this after the freshly-computed NM/MD have been written, so the
	// recomputed NM/MD are dropped too; only RG survives.
	DropTags bool
	// BinQual reduces base-quality resolution (the -q flag): each quality
	// value >= 3 is mapped to qual/10*10 + 7.
	BinQual bool
	// MaxNM, when > 0, masks the matching bases of any read whose computed
	// NM is >= MaxNM (the -n flag). Matching SEQ bases become '=' and their
	// qualities become 0. The emitted NM/MD are unaffected.
	MaxNM int
}

// Calmd reads SAM/BAM records from in, fills in MD + NM aux tags by
// re-walking the CIGAR against the reference FASTA at refPath, and writes
// the resulting stream to out. The output format is BAM when
// opts.OutputBAM or opts.Uncompressed is true, otherwise SAM text.
//
// warnW receives the per-record "different MD/NM" warnings; passing nil
// silences them entirely.
func Calmd(in io.Reader, out io.Writer, refPath string, opts CalmdOptions, warnW io.Writer) error {
	ra, err := fasta.OpenRandomAccess(refPath)
	if err != nil {
		return fmt.Errorf("samtools calmd: open reference %s: %w", refPath, err)
	}
	defer ra.Close()

	r, err := sam.NewReader(in)
	if err != nil {
		return fmt.Errorf("samtools calmd: open input: %w", err)
	}
	hdr := r.Header()

	var w sam.Writer
	if opts.OutputBAM || opts.Uncompressed {
		w = sam.NewBAMWriter(out)
	} else {
		w = sam.NewSAMWriter(out)
	}
	if err := w.WriteHeader(hdr); err != nil {
		return err
	}

	// Cache the reference sequence for the currently-active contig. Records
	// in SAM/BAM coming out of `samtools sort` (or any other coordinate-
	// sorted producer) are clustered by contig, so this single-slot cache
	// hits the common case without bringing the whole reference into memory.
	var curRef string
	var curSeq []byte

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if rec.IsUnmapped() || len(rec.Cigar) == 0 || rec.Seq == "" || rec.Seq == "*" {
			if err := w.Write(rec); err != nil {
				return err
			}
			continue
		}
		if rec.RName != curRef {
			seq, ferr := ra.Fetch(rec.RName, 0, ra.Length(rec.RName))
			if ferr != nil {
				return fmt.Errorf("samtools calmd: fetch %s: %w", rec.RName, ferr)
			}
			curRef = rec.RName
			curSeq = seq
		}
		// BAQ realignment and MAPQ capping run before the MD/NM fill,
		// mirroring upstream bam_md.c (sam_prob_realn -> sam_cap_mapq ->
		// bam_fillmd1_core). sam_prob_realn return codes -1/-3 are benign
		// "nothing to do" results; only -4 (alignment failure) is an error.
		if opts.RealignBAQ {
			flag := 0
			if opts.AdjustCapQ {
				flag |= baq.FlagApply
			}
			if opts.ExtendedBAQ {
				flag |= baq.FlagExtend
			}
			if r := baq.SamProbRealn(rec, curSeq, flag); r < -3 {
				return fmt.Errorf("samtools calmd: BAQ alignment failed for read %q", rec.QName)
			}
		}
		if opts.CapMapQ > 10 {
			if q := baq.SamCapMapq(rec, curSeq, opts.CapMapQ); q >= 0 && int(rec.MapQ) > q {
				rec.MapQ = uint8(q)
			}
		}
		md, nm, edited, eerr := fillMDNM(rec, curSeq, opts.UseEqual)
		if eerr != nil {
			return eerr
		}
		if opts.UseEqual && edited != "" {
			rec.Seq = edited
		}
		// Upstream ordering (bam_md.c): max-NM masking → write NM → write
		// MD → DROP_TAG → BIN_QUAL.
		if opts.MaxNM > 0 && nm >= opts.MaxNM {
			maskMatches(rec, curSeq)
		}
		updateNMAux(rec, nm, opts.Quiet, warnW)
		updateMDAux(rec, md, opts.Quiet, warnW)
		if opts.DropTags {
			dropOtherTags(rec)
		}
		if opts.BinQual {
			binQual(rec)
		}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	return w.Close()
}

// fillMDNM walks rec's CIGAR against ref (0-based contig bases). It returns
// the MD string, the NM count, and — when useEqual is true — a copy of the
// read sequence with M-op match bases rewritten to '='. The edited sequence
// is the empty string when useEqual is false (caller should not assign).
func fillMDNM(rec *sam.Record, ref []byte, useEqual bool) (string, int, string, error) {
	var (
		mdBuf   []byte
		nm      int
		matched int
		qpos    int
		rpos    = int(rec.Pos) - 1 // BAM/0-based on the reference
	)
	if rpos < 0 {
		rpos = 0
	}
	seq := []byte(rec.Seq)
	var edited []byte
	if useEqual {
		edited = make([]byte, len(seq))
		copy(edited, seq)
	}

	flushRun := func() {
		mdBuf = strconv.AppendInt(mdBuf, int64(matched), 10)
		matched = 0
	}

	for _, op := range rec.Cigar {
		opCode := op.Op()
		oplen := int(op.Length())
		switch opCode {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			truncated := false
			for j := 0; j < oplen; j++ {
				if rpos+j >= len(ref) || qpos+j >= len(seq) {
					// Out of bounds — upstream "break" semantics: leave
					// the rest of the read alone. The MD up to here is
					// still emitted in the final flush.
					truncated = true
					break
				}
				if baseMatches(seq[qpos+j], ref[rpos+j]) {
					matched++
					if useEqual {
						edited[qpos+j] = '='
					}
				} else {
					flushRun()
					mdBuf = append(mdBuf, upperByte(ref[rpos+j]))
					nm++
				}
			}
			if truncated {
				// Stop processing further CIGAR ops — flush match-run
				// and return.
				flushRun()
				return string(mdBuf), nm, string(editedOrEmpty(edited, useEqual)), nil
			}
			rpos += oplen
			qpos += oplen
		case sam.CigarDeletion:
			flushRun()
			mdBuf = append(mdBuf, '^')
			actual := 0
			for j := 0; j < oplen; j++ {
				if rpos+j >= len(ref) {
					break
				}
				mdBuf = append(mdBuf, upperByte(ref[rpos+j]))
				actual++
			}
			nm += actual
			rpos += actual
			if actual < oplen {
				// Upstream bam_md.c:121 breaks out of the whole CIGAR loop
				// when the deletion can't be fully consumed from the
				// reference; the trailing kputw(matched=0) at line 129
				// still appends the "0" terminator.
				flushRun()
				return string(mdBuf), nm, string(editedOrEmpty(edited, useEqual)), nil
			}
		case sam.CigarInsertion:
			nm += oplen
			qpos += oplen
		case sam.CigarSoftClip:
			qpos += oplen
		case sam.CigarSkipped:
			rpos += oplen
		case sam.CigarHardClip, sam.CigarPadding:
			// neither consumes query nor reference for our purposes
		}
	}
	flushRun()
	return string(mdBuf), nm, string(editedOrEmpty(edited, useEqual)), nil
}

// editedOrEmpty returns the edited buffer when useEqual is true, otherwise
// a nil slice so the caller's string conversion is the cheap empty string.
func editedOrEmpty(buf []byte, useEqual bool) []byte {
	if !useEqual {
		return nil
	}
	return buf
}

// upperByte folds an ASCII byte to uppercase. Non-letters pass through.
// (Named to avoid collision with the mpileup pileup column formatter's
// own `upper` helper.)
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

// maskMatches rewrites every matching M/=/X-op SEQ base of rec to '=' and
// sets its quality to 0, mirroring upstream's max-NM masking pass
// (bam_md.c:135-155). The CIGAR walk and base-match test reuse the same
// out-of-bounds break semantics and baseMatches comparison as fillMDNM.
func maskMatches(rec *sam.Record, ref []byte) {
	seq := []byte(rec.Seq)
	rpos := int(rec.Pos) - 1
	if rpos < 0 {
		rpos = 0
	}
	qpos := 0
	for _, op := range rec.Cigar {
		oplen := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			truncated := false
			for j := 0; j < oplen; j++ {
				if rpos+j >= len(ref) || qpos+j >= len(seq) {
					truncated = true
					break
				}
				if baseMatches(seq[qpos+j], ref[rpos+j]) {
					seq[qpos+j] = '='
					if qpos+j < len(rec.Qual) {
						rec.Qual[qpos+j] = 0
					}
				}
			}
			if truncated {
				rec.Seq = string(seq)
				return
			}
			rpos += oplen
			qpos += oplen
		case sam.CigarDeletion, sam.CigarSkipped:
			rpos += oplen
		case sam.CigarInsertion, sam.CigarSoftClip:
			qpos += oplen
		}
	}
	rec.Seq = string(seq)
}

// dropOtherTags removes every aux field except RG, mirroring upstream's
// DROP_TAG handling (bam_md.c:199-202). When the record has no RG tag all
// aux fields are dropped.
func dropOtherTags(rec *sam.Record) {
	kept := rec.Aux[:0]
	for _, a := range rec.Aux {
		if a.Tag == "RG" {
			kept = append(kept, a)
		}
	}
	rec.Aux = kept
	rebuildAuxIndex(rec)
}

// binQual reduces the resolution of rec's base qualities: each value >= 3
// becomes qual/10*10 + 7 (integer division), mirroring upstream's BIN_QUAL
// handling (bam_md.c:204-208). Values below 3 are left untouched.
func binQual(rec *sam.Record) {
	for i, q := range rec.Qual {
		if q >= 3 {
			rec.Qual[i] = q/10*10 + 7
		}
	}
}

// updateNMAux sets / replaces the NM:i: aux. Mirrors upstream's "different NM"
// stderr line on mismatch.
func updateNMAux(rec *sam.Record, nm int, quiet bool, warnW io.Writer) {
	for i, a := range rec.Aux {
		if a.Tag == "NM" {
			if old, ok := a.Int(); ok && old != int64(nm) {
				if !quiet && warnW != nil {
					fmt.Fprintf(warnW, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", rec.QName, old, nm)
				}
			}
			rec.Aux[i] = sam.Aux{Tag: "NM", Type: 'i', Value: int64(nm)}
			rebuildAuxIndex(rec)
			return
		}
	}
	rec.Aux = append(rec.Aux, sam.Aux{Tag: "NM", Type: 'i', Value: int64(nm)})
	rebuildAuxIndex(rec)
}

// updateMDAux sets / replaces the MD:Z: aux. Mirrors upstream's "different MD"
// stderr line on mismatch.
func updateMDAux(rec *sam.Record, md string, quiet bool, warnW io.Writer) {
	for i, a := range rec.Aux {
		if a.Tag == "MD" {
			if old, ok := a.String(); ok && old != md {
				if !quiet && warnW != nil {
					fmt.Fprintf(warnW, "[bam_fillmd1] different MD for read '%s': '%s' -> '%s'\n", rec.QName, old, md)
				}
			}
			rec.Aux[i] = sam.Aux{Tag: "MD", Type: 'Z', Value: md}
			rebuildAuxIndex(rec)
			return
		}
	}
	rec.Aux = append(rec.Aux, sam.Aux{Tag: "MD", Type: 'Z', Value: md})
	rebuildAuxIndex(rec)
}

// rebuildAuxIndex drops the record's cached tag→position lookup after the
// NM/MD aux fields have been mutated, so a subsequent GetAux rebuilds it
// from the current Aux slice.
func rebuildAuxIndex(rec *sam.Record) {
	rec.InvalidateAuxIndex()
}

// CalmdFile is a thin path-based wrapper that opens the input file via
// iohelper and dispatches to Calmd. It exists so the CLI driver can pass an
// input path through without having to manage stdin/gzip/BGZF detection itself.
func CalmdFile(inPath string, out io.Writer, refPath string, opts CalmdOptions, warnW io.Writer) error {
	if warnW == nil {
		warnW = os.Stderr
	}
	in, err := iohelper.OpenReader(inPath)
	if err != nil {
		return err
	}
	defer in.Close()
	return Calmd(in, out, refPath, opts, warnW)
}
