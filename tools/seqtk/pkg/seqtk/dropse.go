// dropse.go: port of upstream "seqtk dropse" — read an interleaved
// FASTA/FASTQ stream and emit only records whose immediate neighbour
// in the stream has the same name modulo a trailing "/<digit>" suffix
// (any single digit 0-9, matching upstream's `isdigit` check at
// seqtk.c:1659). Singletons (records whose neighbour does not match)
// are silently dropped. This is the inverse of "mergepe" minus any
// records that have no partner.
//
// Algorithm (1:1 port of reference_code/seqtk/seqtk.c::stk_dropse):
//
//	last := empty
//	for each record r in stream:
//	    if last is set:
//	        if same_pair(last, r):
//	            emit last, emit r
//	            last := empty
//	        else:
//	            last := r        // r becomes the new candidate
//	    else:
//	        last := r
//
// same_pair(p, q) is true iff len(p) == len(q) AND the first len-2
// bytes are equal AND, when p/q both end with a '/<digit>', that
// suffix is ignored. Upstream only honours the '/' separator (not
// '_1'/'_2'); we match upstream byte-for-byte.
//
// Note: the trailing dangling singleton (a record left in "last" at
// EOF with no following match) is dropped, matching upstream.

package seqtk

import (
	"bufio"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// samePairName reports whether the two record names p and q name two
// mates of the same fragment, using upstream's exact comparison from
// stk_dropse: identical length, identical bytes except possibly the
// final two bytes when both records have a '/<digit>' tail.
func samePairName(p, q string) bool {
	if len(p) != len(q) {
		return false
	}
	l := len(p)
	if l > 2 &&
		p[l-2] == '/' && q[l-2] == '/' &&
		isASCIIDigit(p[l-1]) && isASCIIDigit(q[l-1]) {
		l -= 2
	}
	return p[:l] == q[:l]
}

// isASCIIDigit reports whether b is one of '0'..'9'.
func isASCIIDigit(b byte) bool { return b >= '0' && b <= '9' }

// Dropse reads an interleaved FASTA or FASTQ stream from in (format
// auto-detected from the first non-whitespace byte) and writes back
// only the paired records to w. A record at position i is emitted iff
// its name and the name of the record at position i-1 (or i+1) compare
// equal under samePairName. Unpaired records ("singletons") are
// silently dropped, matching upstream "seqtk dropse".
//
// Output preserves the input format (FASTQ stays FASTQ, FASTA stays
// FASTA) and the original record headers, including comments.
func Dropse(in io.Reader, w io.Writer) error {
	br, isFastq := peekIsFastq(in)
	if isFastq {
		return dropseFastq(br, w)
	}
	return dropseFasta(br, w)
}

func dropseFastq(in io.Reader, w io.Writer) error {
	r := fastq.NewReader(in, fastq.Phred33)
	wr := fastq.NewWriter(w, fastq.Phred33)
	var last *fastq.Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if last != nil {
			if samePairName(last.ID, rec.ID) {
				if err := wr.Write(last); err != nil {
					return err
				}
				if err := wr.Write(rec); err != nil {
					return err
				}
				last = nil
				continue
			}
			// Replace candidate with current record (matches upstream's
			// cpy_kseq(&last, seq) branch).
			last = rec
			continue
		}
		last = rec
	}
	return wr.Flush()
}

func dropseFasta(in io.Reader, w io.Writer) error {
	r := fasta.NewReader(in)
	bw := bufio.NewWriter(w)
	var last *fasta.Record
	emit := func(rec *fasta.Record) error {
		// Write a single-line (un-wrapped) FASTA record, matching
		// upstream's stk_printseq with line_len == 0.
		if err := bw.WriteByte('>'); err != nil {
			return err
		}
		if _, err := bw.WriteString(rec.Description); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Sequence); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if last != nil {
			if samePairName(last.ID, rec.ID) {
				if err := emit(last); err != nil {
					return err
				}
				if err := emit(rec); err != nil {
					return err
				}
				last = nil
				continue
			}
			last = rec
			continue
		}
		last = rec
	}
	return bw.Flush()
}
