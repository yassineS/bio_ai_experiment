package bedgroupby

import (
	"bufio"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnbed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// lineSource yields successive input lines for the grouping engine. It hides
// whether the underlying input is a text/BED stream (scanned line by line) or
// a SAM/BAM alignment stream (each mapped record rendered to its groupby
// column line). next reports ok=false at end of input.
type lineSource struct {
	scanner *bufio.Scanner
	bam     *samColumnReader
}

// newLineSource peeks the first bytes of r to decide whether it is a SAM/BAM
// alignment stream or a text stream, and returns a lineSource over the right
// reader. The peek is non-destructive: bytes are buffered and replayed.
func newLineSource(r io.Reader) (*lineSource, error) {
	br := bufio.NewReader(r)
	head, _ := br.Peek(16)
	if alnbed.LooksLikeAlignment(head) {
		bam, err := newSAMColumnReader(br)
		if err != nil {
			return nil, err
		}
		return &lineSource{bam: bam}, nil
	}
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 8*1024*1024)
	return &lineSource{scanner: sc}, nil
}

// next returns the next input line. ok is false once the input is exhausted.
func (s *lineSource) next() (line string, ok bool, err error) {
	if s.bam != nil {
		l, err := s.bam.nextLine()
		if err == io.EOF {
			return "", false, nil
		}
		if err != nil {
			return "", false, err
		}
		return l, true, nil
	}
	if s.scanner.Scan() {
		return s.scanner.Text(), true, nil
	}
	return "", false, s.scanner.Err()
}

// samColumns renders a mapped SAM/BAM alignment into the tab-delimited column
// layout that upstream `bedtools groupby` groups over when its input is a BAM
// file. Unlike `bedtools bamtobed` (which emits BED12), groupby reads the BAM
// through bedtools' BamRecord::getField, whose columns are the SAM fields with
// two adjustments: column 4 is the 0-based start (POS-1, the BED start) and
// column 6 is the CIGAR rendered op-char-then-length (e.g. "5M" -> "M5").
//
// The 11 columns, in order, are:
//
//	1  QNAME  query name
//	2  FLAG   bitwise flag (upstream errors if grouped/aggregated on; emitted
//	          here so non-FLAG columns line up positionally)
//	3  RNAME  reference (chrom) name
//	4  start  0-based left-most position (POS-1)
//	5  MAPQ   mapping quality
//	6  CIGAR  CIGAR string, op char before length ("M5")
//	7  RNEXT  mate reference name ("*"->".", "="->RNAME)
//	8  PNEXT  0-based mate position (PNEXT-1)
//	9  TLEN   observed template length
//	10 SEQ    segment sequence
//	11 QUAL   ASCII Phred+33 base qualities
//
// This mirrors reference_code/bedtools/src/utils/FileRecordTools/Records/
// BamRecord.cpp (getField + buildCigarStr).
func samColumns(rec *sam.Record) string {
	start := int(rec.Pos) - 1

	rnext := rec.RNext
	switch rnext {
	case "", "*":
		rnext = "."
	case "=":
		rnext = rec.RName
	}

	seq := rec.Seq
	if seq == "" {
		seq = "*"
	}

	qual := "*"
	if len(rec.Qual) > 0 {
		hasQual := false
		for _, q := range rec.Qual {
			if q != 0xff {
				hasQual = true
				break
			}
		}
		if hasQual {
			var sb strings.Builder
			sb.Grow(len(rec.Qual))
			for _, q := range rec.Qual {
				sb.WriteByte(q + 33)
			}
			qual = sb.String()
		}
	}

	fields := []string{
		rec.QName,                                // 1 QNAME
		strconv.FormatUint(uint64(rec.Flag), 10), // 2 FLAG
		rec.RName,                                // 3 RNAME
		strconv.Itoa(start),                      // 4 start (0-based)
		strconv.FormatUint(uint64(rec.MapQ), 10), // 5 MAPQ
		bamCigarString(rec.Cigar),                // 6 CIGAR
		rnext,                                    // 7 RNEXT
		strconv.Itoa(int(rec.PNext) - 1),         // 8 PNEXT (0-based)
		strconv.FormatInt(int64(rec.TLen), 10),   // 9 TLEN
		seq,                                      // 10 SEQ
		qual,                                     // 11 QUAL
	}
	return strings.Join(fields, "\t")
}

// bamCigarString renders a CIGAR the way bedtools' BamRecord does: each
// operation is written as its op character followed by its length, in
// alignment order (e.g. 5M -> "M5", 3M2I5M -> "M3I2M5"). An empty CIGAR
// yields "*".
func bamCigarString(c sam.Cigar) string {
	if len(c) == 0 {
		return "*"
	}
	var sb strings.Builder
	for _, op := range c {
		sb.WriteByte(op.Char())
		sb.WriteString(strconv.FormatUint(uint64(op.Length()), 10))
	}
	return sb.String()
}

// samColumnReader streams a SAM/BAM source, yielding for each mapped
// alignment the tab-delimited groupby column line (see samColumns). Unmapped
// reads are skipped, matching upstream's BAM-to-record conversion which only
// groups mapped alignments.
type samColumnReader struct {
	src sam.Reader
}

// newSAMColumnReader wraps r (auto-detecting SAM text or BAM) for column-line
// streaming.
func newSAMColumnReader(r io.Reader) (*samColumnReader, error) {
	sr, err := sam.NewReader(r)
	if err != nil {
		return nil, err
	}
	return &samColumnReader{src: sr}, nil
}

// nextLine returns the next mapped alignment rendered as a groupby column
// line, or io.EOF when the stream is exhausted.
func (r *samColumnReader) nextLine() (string, error) {
	for {
		rec, err := r.src.Read()
		if err != nil {
			return "", err
		}
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" || rec.Pos <= 0 {
			continue
		}
		return samColumns(rec), nil
	}
}
