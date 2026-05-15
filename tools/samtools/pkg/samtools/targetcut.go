// Package samtools — targetcut subcommand.
//
// Targetcut walks a SAM/BAM stream and, for each mapped primary record,
// emits a FASTA entry containing the slice of the read sequence that
// aligns to the reference. Soft-clipped bases at either end are
// dropped; matches/mismatches and insertions contribute their query
// bases; deletions/skips/hard-clips/padding contribute nothing because
// they consume no query bases.
//
// The output is one record per line-pair:
//
//	>QNAME
//	<cut-sequence>
//
// Unmapped, secondary, supplementary, duplicate, and QC-fail records
// are skipped. Records whose SEQ is "*" or whose alignment has no
// CIGAR are skipped (there is no sequence to emit).
//
// Per-base minimum quality (-Q) filters bases out of the emitted
// sequence: a query base whose Phred quality is below MinBaseQ is
// excluded from the cut sequence. When SEQ is present but QUAL is
// missing ("*"), every base is treated as having quality 255 (i.e.
// unknown / unfiltered), matching SAM convention.
//
// This is a deliberately small, focused tool intended for slicing read
// fragments out of an aligned BAM. Upstream samtools ships a
// like-named subcommand (cut_target.c) that does something different —
// HMM consensus calling over fosmid pools. The Go reimplementation
// here follows the more commonly-requested "cut the aligned slice from
// each read" behaviour; the upstream HMM mode is tracked as deferred
// in docs/PARITY_ROADMAP.md.
package samtools

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// TargetcutOptions configures the Targetcut walker. The zero value
// matches the upstream defaults (min base quality 13, no extra
// filters).
type TargetcutOptions struct {
	// MinBaseQ drops query bases whose Phred quality is strictly less
	// than this. The upstream default is 13.
	MinBaseQ uint8
}

// DefaultTargetcutMinBaseQ is the default minimum base quality used by
// the CLI when -Q is not supplied.
const DefaultTargetcutMinBaseQ uint8 = 13

// Targetcut reads SAM/BAM records from in and writes one FASTA entry
// per kept record to out. It returns the number of records emitted
// and the first error encountered (if any).
func Targetcut(in io.Reader, out io.Writer, opts TargetcutOptions) (int, error) {
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, fmt.Errorf("samtools targetcut: open input: %w", err)
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	emitted := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return emitted, fmt.Errorf("samtools targetcut: read: %w", err)
		}
		if skipForTargetcut(rec) {
			continue
		}
		seq := cutSequence(rec, opts.MinBaseQ)
		if len(seq) == 0 {
			continue
		}
		if _, err := bw.WriteString(">"); err != nil {
			return emitted, err
		}
		if _, err := bw.WriteString(rec.QName); err != nil {
			return emitted, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return emitted, err
		}
		if _, err := bw.Write(seq); err != nil {
			return emitted, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return emitted, err
		}
		emitted++
	}
	return emitted, nil
}

// skipForTargetcut returns true if the record should be excluded from
// the FASTA output. We mirror the upstream cut_target.c read filter:
// drop unmapped / secondary / QC-fail / duplicate. Supplementary reads
// are also skipped (they would emit a partial slice that overlaps the
// primary alignment).
func skipForTargetcut(rec *sam.Record) bool {
	if rec.IsUnmapped() {
		return true
	}
	if rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
		return true
	}
	if rec.Seq == "" || rec.Seq == "*" {
		return true
	}
	if len(rec.Cigar) == 0 {
		return true
	}
	return false
}

// cutSequence returns the slice of rec.Seq corresponding to the
// aligned portion of the read. Soft-clip and hard-clip bases at
// either end are excluded; insertions are included (they are query
// bases that do not advance the reference); deletions / refskip /
// padding contribute nothing because they consume no query bases.
//
// Bases whose Phred quality is below minBaseQ are omitted. A QUAL
// field of "*" (length zero) is treated as "quality unknown" and no
// base is filtered.
func cutSequence(rec *sam.Record, minBaseQ uint8) []byte {
	qual := rec.Qual
	hasQual := len(qual) == len(rec.Seq)
	out := make([]byte, 0, len(rec.Seq))
	qpos := 0
	for _, op := range rec.Cigar {
		o := op.Op()
		n := int(op.Length())
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarInsertion:
			// Query-consuming + part of the aligned slice (insertions
			// included; soft-clip is handled in the other branch).
			for k := 0; k < n; k++ {
				idx := qpos + k
				if idx >= len(rec.Seq) {
					break
				}
				if hasQual && qual[idx] < minBaseQ {
					continue
				}
				out = append(out, rec.Seq[idx])
			}
			qpos += n
		case sam.CigarSoftClip:
			// Query-consuming but NOT part of the aligned slice.
			qpos += n
		case sam.CigarDeletion, sam.CigarSkipped, sam.CigarPadding, sam.CigarHardClip:
			// No query bases consumed; nothing to emit.
		default:
			// CigarBack and unknown ops: ignore.
		}
	}
	return out
}
