// Package alnbed converts aligned SAM/BAM records into BED intervals, the
// way upstream bedtools tools treat BAM input (the BamRecord::print /
// bam2bed conversion). Each mapped alignment becomes a BED12 record whose
// blocks are the reference-consuming runs of the CIGAR broken on N
// (skip/intron) operations, so a consumer that honours BED12 blocks gets
// "-split" (exon-aware) behaviour for free, and one that ignores blocks
// sees the whole-read span.
//
// The package is the shared foundation for BAM input across the bed* tools
// (genomecov, coverage, jaccard, …): a Reader yields *bed.Record values
// through the same Read() contract as bed.Reader, and References() exposes
// the @SQ chromosome ordering and lengths so tools like genomecov can take
// their genome from the BAM header.
package alnbed

import (
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Reader adapts a SAM/BAM stream to the bed.Reader Read() contract: each
// call returns the next mapped alignment as a BED12 *bed.Record. Unmapped
// reads (and reads with no reference name) are skipped, matching upstream's
// BAM-to-BED conversion which only emits mapped alignments.
type Reader struct {
	src sam.Reader
}

// NewReader wraps r (auto-detecting SAM text or BAM) and returns a Reader.
func NewReader(r io.Reader) (*Reader, error) {
	sr, err := sam.NewReader(r)
	if err != nil {
		return nil, err
	}
	return &Reader{src: sr}, nil
}

// References returns the @SQ entries (chromosome name + length) from the
// alignment header, in header order. Callers such as genomecov use this to
// build the genome when the input is a BAM (`-ibam`), where no separate
// genome file is supplied.
func (r *Reader) References() []sam.Reference {
	h := r.src.Header()
	if h == nil {
		return nil
	}
	return h.Refs
}

// Read returns the next mapped alignment as a BED12 *bed.Record, or io.EOF
// when the stream is exhausted. Unmapped reads are skipped.
func (r *Reader) Read() (*bed.Record, error) {
	for {
		rec, err := r.src.Read()
		if err != nil {
			return nil, err
		}
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" || rec.Pos <= 0 {
			continue
		}
		return ToBED12(rec), nil
	}
}

// LooksLikeAlignment reports whether head (the first bytes of an input
// stream) is a SAM or BAM alignment rather than a BED/text interval file. It
// recognises a BGZF-wrapped BAM, a raw "BAM\1" magic, and SAM text (which
// always begins with an '@' header line). A BED file never starts with '@'
// or the BAM magic, so this is a safe sniff for auto-routing tool input.
func LooksLikeAlignment(head []byte) bool {
	if len(head) >= 4 && head[0] == 'B' && head[1] == 'A' && head[2] == 'M' && head[3] == 0x01 {
		return true
	}
	if looksLikeBGZF(head) {
		return true
	}
	// SAM text starts with a header line ('@HD', '@SQ', ...).
	return len(head) >= 1 && head[0] == '@'
}

// looksLikeBGZF reports whether b begins with a BGZF gzip header (gzip magic
// with the BC subfield), the container BAM is stored in.
func looksLikeBGZF(b []byte) bool {
	if len(b) < 16 || b[0] != 0x1f || b[1] != 0x8b || b[2] != 0x08 || b[3]&0x04 == 0 {
		return false
	}
	if xlen := uint16(b[10]) | uint16(b[11])<<8; xlen < 6 {
		return false
	}
	return b[12] == 'B' && b[13] == 'C' && b[14] == 0x02 && b[15] == 0x00
}

// block is one reference-consuming run of an alignment.
type block struct{ start, end int }

// ToBED12 converts a single mapped SAM/BAM alignment into a BED12
// *bed.Record. Coordinates are 0-based half-open (BAM Pos is 1-based). The
// blocks are the CIGAR's reference-consuming runs, broken on N ops; Name is
// the read name, Score the MAPQ, Strand from the reverse flag. The whole-read
// span is [ChromStart, ChromEnd); a consumer that ignores blocks therefore
// sees the full span, while a -split-aware consumer sees the exon blocks.
func ToBED12(rec *sam.Record) *bed.Record {
	start := int(rec.Pos) - 1
	blocks := cigarBlocks(rec, start)
	end := start
	if len(blocks) > 0 {
		end = blocks[len(blocks)-1].end
	}
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	sizes := make([]int, len(blocks))
	starts := make([]int, len(blocks))
	for i, b := range blocks {
		sizes[i] = b.end - b.start
		starts[i] = b.start - start
	}
	// Build the raw BED12 column slice so callers that echo the record (or
	// re-read it through bed.NewReader) see a faithful 12-column line.
	var ssz, sst strings.Builder
	for _, s := range sizes {
		ssz.WriteString(strconv.Itoa(s))
		ssz.WriteByte(',')
	}
	for _, s := range starts {
		sst.WriteString(strconv.Itoa(s))
		sst.WriteByte(',')
	}
	return &bed.Record{
		Chrom:       rec.RName,
		ChromStart:  start,
		ChromEnd:    end,
		Name:        rec.QName,
		Score:       int(rec.MapQ),
		Strand:      strand,
		ThickStart:  start,
		ThickEnd:    end,
		ItemRGB:     "0,0,0",
		BlockCount:  len(blocks),
		BlockSizes:  sizes,
		BlockStarts: starts,
		ExtraFields: nil,
	}
}

// cigarBlocks returns the reference-consuming blocks of an alignment starting
// at 0-based refStart, breaking on each N (skip) op. M/=/X/D consume the
// reference; D extends the current block (a deletion does not split it,
// matching upstream's breakOnDeletionOps=false); N closes the current block
// and starts a new one after the gap; I/S/H/P do not consume reference.
func cigarBlocks(rec *sam.Record, refStart int) []block {
	var blocks []block
	cur := block{start: refStart, end: refStart}
	open := false
	pos := refStart
	for _, op := range rec.Cigar {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarDeletion:
			if !open {
				cur = block{start: pos, end: pos}
				open = true
			}
			pos += int(op.Length())
			cur.end = pos
		case sam.CigarSkipped:
			if open {
				blocks = append(blocks, cur)
				open = false
			}
			pos += int(op.Length())
		default:
			// I, S, H, P: no reference advance.
		}
	}
	if open {
		blocks = append(blocks, cur)
	}
	if len(blocks) == 0 {
		blocks = append(blocks, block{start: refStart, end: refStart})
	}
	return blocks
}
