package bamtobed

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FromBAMSAMText wraps a BGZF-wrapped BAM byte stream and emits one
// tab-separated row per primary alignment with columns matching how
// upstream `bedtools groupby -i x.bam` ingests BAM input. Specifically:
//
//	col 1: QNAME           col 7: RNEXT
//	col 2: FLAG            col 8: PNEXT (0-based, as stored in BAM)
//	col 3: RNAME           col 9: TLEN
//	col 4: POS (0-based)   col 10: SEQ
//	col 5: MAPQ            col 11: QUAL
//	col 6: CIGAR           col 12+: AUX...
//
// POS is emitted 0-based (matching bedtools' internal `BAM.Position`
// field) rather than the 1-based SAM convention — this is what makes
// `bedtools groupby -i x.bam -g 1,3 -c 4 -o mean` produce upstream's
// reference values for tests like groupby.t17. The header is consumed
// but dropped; records that fail the standard alignment filter are kept
// here unlike FromBAM because groupby has no concept of alignment
// filtering.
func FromBAMSAMText(r io.Reader) io.Reader {
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
			if werr := writeBAMSAMTextRow(bw, rec); werr != nil {
				_ = pw.CloseWithError(werr)
				return
			}
		}
	}()
	return pr
}

// writeBAMSAMTextRow renders one alignment as the BAM-as-TSV line
// described on FromBAMSAMText. POS and PNEXT are zero-based.
func writeBAMSAMTextRow(bw *bufio.Writer, rec *sam.Record) error {
	rname := rec.RName
	if rname == "" {
		rname = "*"
	}
	rnext := rec.RNext
	if rnext == "" {
		rnext = "*"
	}
	cigar := rec.Cigar.String()
	seq := rec.Seq
	if seq == "" {
		seq = "*"
	}
	qual := "*"
	if len(rec.Qual) > 0 {
		all0xff := true
		for _, b := range rec.Qual {
			if b != 0xff {
				all0xff = false
				break
			}
		}
		if !all0xff {
			buf := make([]byte, len(rec.Qual))
			for i, b := range rec.Qual {
				buf[i] = b + 33
			}
			qual = string(buf)
		}
	}
	pos := int64(rec.Pos) - 1
	pnext := int64(rec.PNext) - 1
	if _, err := fmt.Fprintf(bw,
		"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		rec.QName, rec.Flag, rname, pos, rec.MapQ, cigar,
		rnext, pnext, rec.TLen, seq, qual); err != nil {
		return err
	}
	for _, a := range rec.Aux {
		if err := bw.WriteByte('\t'); err != nil {
			return err
		}
		if _, err := bw.WriteString(a.FormatSAM()); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}
