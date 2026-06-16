// Package bedbamtobed is a pure-Go reimplementation of `bedtools bamtobed`
// (aka bamToBed): it converts BAM/SAM alignments into BED6, blocked BED12, or
// BEDPE records.
//
// The conversion mirrors upstream bedtools v2.31.1 byte-for-byte, including its
// CIGAR→blocks splitting (on N for -split, on N and D for -splitD), the
// QNAME/1 and /2 mate suffixing, MAPQ-as-score default, and the two-mate
// joining performed for -bedpe. Where upstream has long-standing quirks (the
// extra column emitted by `-tag` combined with `-split`), this package
// reproduces them so existing pipelines keep working; those quirks are
// documented in docs/UPSTREAM_BUGS.md.
package bedbamtobed

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Options configures a BAM→BED conversion. The field set mirrors the upstream
// bamtobed command-line flags.
type Options struct {
	// WriteBedPE selects BEDPE output (one line per read pair). Mutually
	// exclusive with WriteBed12.
	WriteBedPE bool
	// Mate1First, valid only with WriteBedPE, always reports mate one as the
	// first BEDPE block.
	Mate1First bool
	// WriteBed12 selects blocked BED12 output. Forces splitting.
	WriteBed12 bool
	// ObeySplits reports each N-separated alignment block as its own BED line
	// (the -split flag). WriteBed12 implies this.
	ObeySplits bool
	// SplitOnDeletions also breaks blocks on D (deletion) CIGAR ops in addition
	// to N (the -splitD flag); it implies ObeySplits.
	SplitOnDeletions bool
	// UseEditDistance uses the NM tag as the BED score (the -ed flag). It is
	// shorthand for Tag = "NM" in BED mode and, in BEDPE mode, reports the
	// summed edit distance of the two mates.
	UseEditDistance bool
	// Tag, when non-empty, names a numeric BAM aux tag whose value is used as
	// the BED score in place of MAPQ.
	Tag string
	// Color is the R,G,B itemRgb string used for the 9th BED12 column.
	Color string
	// UseCigar appends the CIGAR string as a trailing column (BED6 only).
	UseCigar bool
}

// Validate checks an Options for the mutually exclusive / dependent flag
// combinations that upstream rejects, returning a non-nil error describing the
// first violation. The error text names the offending flag pairing.
func (o *Options) Validate() error {
	haveTag := o.Tag != "" && !o.UseEditDistance
	if o.UseEditDistance && o.ObeySplits {
		return fmt.Errorf("Cannot use -ed with -splits. Edit distances cannot be computed for each 'chunk'.")
	}
	if o.UseEditDistance && o.UseCigar {
		return fmt.Errorf("Cannot use -cigar with -splits.  Not yet supported.")
	}
	if o.UseEditDistance && haveTag {
		return fmt.Errorf("Cannot use -ed with -tag.  Choose one or the other.")
	}
	if o.WriteBedPE && haveTag {
		return fmt.Errorf("Cannot use -tag with -bedpe.")
	}
	if !o.WriteBedPE && o.Mate1First {
		return fmt.Errorf("Must use -mate1 with -bedpe.")
	}
	return nil
}

// scoreTag returns the aux tag to read for the BED score, accounting for -ed
// being shorthand for the NM tag.
func (o *Options) scoreTag() string {
	if o.UseEditDistance {
		return "NM"
	}
	return o.Tag
}

// block is a half-open [start,end) reference interval produced from a CIGAR.
type block struct {
	start, end int
}

// Run reads SAM or BAM alignments from r and writes the converted BED, BED12,
// or BEDPE lines to w according to opts. It returns the number of input
// alignment records consumed (matching upstream, which reads every record even
// when some produce no output).
func Run(r io.Reader, w io.Writer, opts Options) (int, error) {
	if err := opts.Validate(); err != nil {
		return 0, err
	}
	sr, err := sam.NewReader(r)
	if err != nil {
		return 0, err
	}
	if opts.WriteBedPE {
		return convertBedPE(sr, w, opts)
	}
	return convertBed(sr, w, opts)
}

// convertBed implements the non-BEDPE conversion (BED6 / BED12).
func convertBed(sr sam.Reader, w io.Writer, opts Options) (int, error) {
	n := 0
	for {
		rec, err := sr.Read()
		if err == io.EOF {
			return n, nil
		}
		if err != nil {
			return n, err
		}
		n++
		if rec.IsUnmapped() {
			continue
		}
		if opts.WriteBed12 {
			if err := printBed12(w, rec, opts); err != nil {
				return n, err
			}
		} else {
			if err := printBed(w, rec, opts); err != nil {
				return n, err
			}
		}
	}
}

// recordName builds the BED name: QNAME plus a /1 or /2 mate suffix when the
// record is flagged first or second in pair. (Upstream appends both if both
// bits are set, but in practice only one is.)
func recordName(rec *sam.Record) string {
	name := rec.QName
	if rec.IsRead1() {
		name += "/1"
	}
	if rec.IsRead2() {
		name += "/2"
	}
	return name
}

// strandOf returns "+" or "-" for the record's strand.
func strandOf(rec *sam.Record) string {
	if rec.Flag&sam.FlagReverse != 0 {
		return "-"
	}
	return "+"
}

// printTag reads the named numeric aux tag from rec and returns its decimal
// string. It returns an error naming the tag if it is absent, matching upstream
// (which exits 1 after the requested tag cannot be found).
func printTag(rec *sam.Record, tag string) (string, error) {
	a, ok := rec.GetAux(tag)
	if !ok {
		return "", fmt.Errorf("The requested tag (%s) was not found in the BAM file.  Exiting", tag)
	}
	if v, ok := a.Int(); ok {
		return strconv.FormatInt(v, 10), nil
	}
	// Non-integer tag: upstream's GetTag<uint32>/GetTag<int32> would fail and
	// it would report the tag as missing.
	return "", fmt.Errorf("The requested tag (%s) was not found in the BAM file.  Exiting", tag)
}

// printBed writes a single record as BED6 (optionally split into blocks, with
// CIGAR or a tag-based score).
func printBed(w io.Writer, rec *sam.Record, opts Options) error {
	strand := strandOf(rec)
	name := recordName(rec)
	chrom := rec.RName
	tag := opts.scoreTag()

	if !opts.ObeySplits {
		// Whole footprint as a single BED entry: [pos, pos+refLen).
		start := int(rec.Pos) - 1
		end := start + rec.Cigar.ReferenceLength()
		var score string
		if tag == "" {
			score = strconv.Itoa(int(rec.MapQ))
		} else {
			s, err := printTag(rec, tag)
			if err != nil {
				// Upstream prints the partial line up to the score column,
				// then the error to stderr and exits. We surface the partial
				// prefix on w too for parity, then return the error.
				fmt.Fprintf(w, "%s\t%d\t%d\t%s\t", chrom, start, end, name)
				return err
			}
			score = s
		}
		if !opts.UseCigar {
			_, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%s\t%s\n", chrom, start, end, name, score, strand)
			return err
		}
		_, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%s\t%s\t%s\n", chrom, start, end, name, score, strand, buildCigar(rec.Cigar))
		return err
	}

	// Split mode: one BED line per block.
	blocks := getBamBlocks(rec, opts.SplitOnDeletions, true)
	for _, b := range blocks {
		if tag == "" {
			if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%d\t%s\n",
				chrom, b.start, b.end, name, int(rec.MapQ), strand); err != nil {
				return err
			}
		} else {
			// Upstream quirk: in the tag+split branch it prints an extra
			// column (bam.Position) before the block start/end. Reproduced for
			// byte-parity (see docs/UPSTREAM_BUGS.md).
			s, err := printTag(rec, tag)
			if err != nil {
				return err
			}
			if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%s\t%s\t%s\n",
				chrom, int(rec.Pos)-1, b.start, b.end, name, s, strand); err != nil {
				return err
			}
		}
	}
	return nil
}

// printBed12 writes a single record as blocked BED12.
func printBed12(w io.Writer, rec *sam.Record, opts Options) error {
	strand := strandOf(rec)
	name := recordName(rec)
	chrom := rec.RName
	pos := int(rec.Pos) - 1
	alignmentEnd := pos + rec.Cigar.ReferenceLength()
	tag := opts.scoreTag()

	blocks := getBamBlocks(rec, opts.SplitOnDeletions, true)

	var score string
	if tag == "" {
		score = strconv.Itoa(int(rec.MapQ))
	} else {
		s, err := printTag(rec, tag)
		if err != nil {
			fmt.Fprintf(w, "%s\t%d\t%d\t%s\t", chrom, pos, alignmentEnd, name)
			return err
		}
		score = s
	}

	color := opts.Color
	if color == "" {
		color = "255,0,0"
	}

	var sb strings.Builder
	fmt.Fprintf(&sb, "%s\t%d\t%d\t%s\t%s\t%s\t", chrom, pos, alignmentEnd, name, score, strand)
	fmt.Fprintf(&sb, "%d\t%d\t%s\t%d\t", pos, alignmentEnd, color, len(blocks))

	// block sizes
	for i, b := range blocks {
		if i > 0 {
			sb.WriteByte(',')
		}
		sb.WriteString(strconv.Itoa(b.end - b.start))
	}
	sb.WriteByte('\t')
	// block starts (relative to pos)
	for i, b := range blocks {
		if i > 0 {
			sb.WriteByte(',')
		}
		sb.WriteString(strconv.Itoa(b.start - pos))
	}
	sb.WriteByte('\n')

	_, err := io.WriteString(w, sb.String())
	return err
}

// getBamBlocks mirrors upstream GetBamBlocks: it walks the CIGAR and produces
// reference-coordinate blocks, breaking on N (skip) ops always when
// breakOnSkipOps is set, and additionally on D (deletion) ops when
// breakOnDeletionOps is set. M/=/X consume reference into the current block;
// I/S/P/H do not advance the reference. A final block is always emitted even
// when zero-length, matching upstream.
func getBamBlocks(rec *sam.Record, breakOnDeletionOps, breakOnSkipOps bool) []block {
	var blocks []block
	currPosition := int(rec.Pos) - 1
	blockLength := 0
	for _, op := range rec.Cigar {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			blockLength += int(op.Length())
		case sam.CigarDeletion:
			if !breakOnDeletionOps {
				blockLength += int(op.Length())
			} else {
				if blockLength > 0 {
					blocks = append(blocks, block{currPosition, currPosition + blockLength})
				}
				currPosition += int(op.Length()) + blockLength
				blockLength = 0
			}
		case sam.CigarSkipped:
			if !breakOnSkipOps {
				blockLength += int(op.Length())
			} else {
				if blockLength > 0 {
					blocks = append(blocks, block{currPosition, currPosition + blockLength})
				}
				currPosition += int(op.Length()) + blockLength
				blockLength = 0
			}
		case sam.CigarInsertion, sam.CigarSoftClip, sam.CigarPadding, sam.CigarHardClip:
			// no reference advance, no block contribution
		}
	}
	blocks = append(blocks, block{currPosition, currPosition + blockLength})
	return blocks
}

// buildCigar renders a CIGAR using only the ops upstream's BuildCigarString
// emits (M/=/X/I/D/N/S/H/P), which is every op our parser produces. It matches
// Cigar.String() in practice but is kept explicit for parity with upstream.
func buildCigar(c sam.Cigar) string {
	var sb strings.Builder
	for _, op := range c {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch,
			sam.CigarInsertion, sam.CigarDeletion, sam.CigarSkipped,
			sam.CigarSoftClip, sam.CigarHardClip, sam.CigarPadding:
			sb.WriteString(strconv.FormatUint(uint64(op.Length()), 10))
			sb.WriteByte(op.Char())
		}
	}
	return sb.String()
}
