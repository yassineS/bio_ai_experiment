package bamtobed

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FromBAMBED12 wraps a BGZF-wrapped BAM byte stream and emits BED12 text
// lines: `chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\
// titemRgb\tblockCount\tblockSizes\tblockStarts`. One row per primary
// alignment. Block layout comes from the CIGAR: each contiguous
// reference-consuming run (M/=/X) becomes one block; N (skip) and D
// (deletion) break blocks. Soft-clip / hard-clip / insertion / padding
// ops consume no reference.
//
// Used by `bedtools coverage -abam`: A flows verbatim to the output, so
// BAM records must be rendered as BED12 (which is what upstream does).
// Score is MAPQ; thickStart/thickEnd equal start/end; itemRgb is "0,0,0".
func FromBAMBED12(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		br, err := sam.NewBAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open BAM: %w", err))
			return
		}
		for {
			rec, err := br.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read BAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			if err := writeBED12(bw, rec); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// FromSAMBED12 is the SAM-text counterpart of FromBAMBED12 (same
// output format, same filter, same CIGAR-walk block layout).
func FromSAMBED12(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		sr, err := sam.NewSAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open SAM: %w", err))
			return
		}
		for {
			rec, err := sr.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read SAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			if err := writeBED12(bw, rec); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// writeBED12 renders one alignment record as a tab-separated BED12 line
// to bw. Shared between FromBAMBED12 and FromSAMBED12.
func writeBED12(bw *bufio.Writer, rec *sam.Record) error {
	refLen := rec.Cigar.ReferenceLength()
	if refLen <= 0 {
		return nil
	}
	start := int(rec.Pos) - 1
	if start < 0 {
		return nil
	}
	end := start + refLen
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	name := rec.QName
	if name == "" {
		name = "."
	}
	var sizes, starts []int
	cur := start
	runLen := 0
	runStart := cur
	flushBlock := func() {
		if runLen <= 0 {
			return
		}
		sizes = append(sizes, runLen)
		starts = append(starts, runStart-start)
		runLen = 0
	}
	for _, op := range rec.Cigar {
		oc := op.Op()
		l := int(op.Length())
		switch oc {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if runLen == 0 {
				runStart = cur
			}
			runLen += l
			cur += l
		case sam.CigarSkipped, sam.CigarDeletion:
			flushBlock()
			cur += l
			runStart = cur
		default:
		}
	}
	flushBlock()
	if len(sizes) == 0 {
		sizes = append(sizes, end-start)
		starts = append(starts, 0)
	}
	_, err := fmt.Fprintf(bw,
		"%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,0\t%d\t%s\t%s\n",
		rec.RName, start, end, name, rec.MapQ, strand,
		start, end, len(sizes), joinIntsComma(sizes), joinIntsComma(starts))
	return err
}

// joinIntsComma renders a slice of ints as "a,b,c," (trailing comma is
// what UCSC BED12 expects).
func joinIntsComma(xs []int) string {
	if len(xs) == 0 {
		return ""
	}
	var sb strings.Builder
	for _, x := range xs {
		sb.WriteString(strconv.Itoa(x))
		sb.WriteByte(',')
	}
	return sb.String()
}
