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

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/mdnm"
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
	// Threads is upstream's -@/--threads worker count. When > 1 it drives
	// block-parallel BGZF inflate on the input and, for compressed BAM
	// output, block-parallel BGZF deflate on the output. The MD/NM compute
	// itself is single-threaded (as in upstream bam_md.c); only the I/O
	// (de)compression is parallelised, so the emitted records are identical
	// regardless of the worker count.
	Threads int
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

	r, err := alnio.NewReaderThreaded(in, "", ReadDecodeThreads(opts.Threads))
	if err != nil {
		return fmt.Errorf("samtools calmd: open input: %w", err)
	}
	if rc, ok := r.(io.Closer); ok {
		defer rc.Close()
	}
	hdr := r.Header()

	var w sam.Writer
	if opts.OutputBAM || opts.Uncompressed {
		// Compressed BAM output honours -@ via the block-parallel BGZF
		// deflate writer; uncompressed BAM (-u) has no compression to
		// parallelise, so it keeps the plain writer.
		if opts.OutputBAM && !opts.Uncompressed && opts.Threads > 1 {
			w = sam.NewBAMWriterThreads(out, opts.Threads)
		} else {
			bw, err := sam.NewBAMWriterOptions(out, sam.BAMWriterOptions{Uncompressed: opts.Uncompressed})
			if err != nil {
				return err
			}
			w = bw
		}
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

// fillMDNM computes rec's MD string and NM count against ref (0-based
// contig bases) and — when useEqual is true — returns a copy of the read
// sequence with M-op match bases rewritten to '='. The edited sequence is
// the empty string when useEqual is false (caller should not assign).
//
// The MD/NM computation is delegated to the shared mdnm.Compute so calmd
// and the CRAM decoder share one implementation; ref is the whole contig
// here, so refOffset is 0. The '=' sequence rewrite (the calmd-only -e
// behaviour) is applied here, reusing the same base-match test mdnm uses.
func fillMDNM(rec *sam.Record, ref []byte, useEqual bool) (string, int, string, error) {
	md, nm := mdnm.Compute(rec, ref, 0)
	if !useEqual {
		return md, nm, "", nil
	}
	return md, nm, equalRewrite(rec, ref), nil
}

// equalRewrite returns a copy of rec's read sequence with every M/=/X-op base
// that matches the reference rewritten to '=', mirroring upstream calmd's -e
// flag. The walk reuses the same reference-coordinate bookkeeping and
// out-of-bounds break semantics as mdnm.Compute.
func equalRewrite(rec *sam.Record, ref []byte) string {
	seq := []byte(rec.Seq)
	edited := make([]byte, len(seq))
	copy(edited, seq)
	rpos := int(rec.Pos) - 1
	if rpos < 0 {
		rpos = 0
	}
	qpos := 0
	for _, op := range rec.Cigar {
		oplen := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for j := 0; j < oplen; j++ {
				if rpos+j >= len(ref) || qpos+j >= len(seq) {
					return string(edited)
				}
				if baseMatches(seq[qpos+j], ref[rpos+j]) {
					edited[qpos+j] = '='
				}
			}
			rpos += oplen
			qpos += oplen
		case sam.CigarDeletion, sam.CigarSkipped:
			rpos += oplen
		case sam.CigarInsertion, sam.CigarSoftClip:
			qpos += oplen
		}
	}
	return string(edited)
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

// updateNMAux sets / replaces the NM:i: aux, matching upstream bam_md.c's
// exact aux ordering: when an existing NM already equals the computed value
// it is left in place; when it differs (or is absent) the old entry is
// removed and the new one is appended at the END of the aux list
// (bam_aux_del + bam_aux_append). The "different NM" stderr line fires only
// on the mismatch path.
func updateNMAux(rec *sam.Record, nm int, quiet bool, warnW io.Writer) {
	for i, a := range rec.Aux {
		if a.Tag == "NM" {
			if old, ok := a.Int(); ok && old == int64(nm) {
				return // unchanged: keep the existing tag in place.
			}
			if old, ok := a.Int(); ok && !quiet && warnW != nil {
				fmt.Fprintf(warnW, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", rec.QName, old, nm)
			}
			rec.Aux = append(rec.Aux[:i], rec.Aux[i+1:]...)
			break
		}
	}
	rec.Aux = append(rec.Aux, sam.Aux{Tag: "NM", Type: 'i', Value: int64(nm)})
	rebuildAuxIndex(rec)
}

// updateMDAux sets / replaces the MD:Z: aux, matching upstream bam_md.c's
// exact aux ordering: an existing MD equal to the computed value is left in
// place; a differing (or absent) MD is removed and the new value appended at
// the END of the aux list. The "different MD" stderr line fires only on the
// mismatch path.
func updateMDAux(rec *sam.Record, md string, quiet bool, warnW io.Writer) {
	for i, a := range rec.Aux {
		if a.Tag == "MD" {
			if old, ok := a.String(); ok && old == md {
				return // unchanged: keep the existing tag in place.
			}
			if old, ok := a.String(); ok && !quiet && warnW != nil {
				fmt.Fprintf(warnW, "[bam_fillmd1] different MD for read '%s': '%s' -> '%s'\n", rec.QName, old, md)
			}
			rec.Aux = append(rec.Aux[:i], rec.Aux[i+1:]...)
			break
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
	// With -@ >= 2 open the raw (still-BGZF-framed) bytes so Calmd's
	// NewReaderThreaded can inflate the blocks in parallel; the standard
	// decompressing opener would hand it an already-inflated stream, so the
	// parallel input decode would never engage. This mirrors samtools sort /
	// stats (openStatsInput). The decoded records are identical either way.
	var in io.ReadCloser
	var err error
	if opts.Threads >= 2 {
		in, err = iohelper.OpenRaw(inPath)
	} else {
		in, err = iohelper.OpenReader(inPath)
	}
	if err != nil {
		return err
	}
	defer in.Close()
	return Calmd(in, out, refPath, opts, warnW)
}
