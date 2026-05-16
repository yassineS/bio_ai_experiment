// pair.go: project-extension "seqtk pair" — split an interleaved
// FASTA/FASTQ stream into two parallel mate streams (the inverse of
// "seqtk mergepe"). Upstream seqtk v1.5 has no "pair" subcommand
// (verified against reference_code/seqtk/seqtk.c v1.5-r133); this is
// a project-original convenience that does not extend the upstream
// flag surface — only positional arguments are accepted.
//
// Algorithm: read pairs of consecutive records from the input. The
// (2k)-th record is written verbatim to out1, the (2k+1)-th to out2.
// If the input contains an odd number of records the trailing
// singleton is treated as an error (the file is not actually paired).
//
// Mate-name parity (e.g. matching '/1' vs '/2') is NOT enforced: the
// records' positions in the stream are taken to define the pairing,
// in line with how "mergepe" interleaves them. Callers that need
// strict name parity should run "dropse" first.

package seqtk

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Pair reads an interleaved FASTA or FASTQ stream from in and writes
// every odd-indexed record (1st, 3rd, 5th, ...) to out1 and every
// even-indexed record (2nd, 4th, 6th, ...) to out2, preserving the
// input format on both outputs. An error is returned if the input
// contains an odd total number of records, indicating an unpaired
// trailing read.
func Pair(in io.Reader, out1, out2 io.Writer) error {
	br, isFastq := peekIsFastq(in)
	if isFastq {
		return pairFastq(br, out1, out2)
	}
	return pairFasta(br, out1, out2)
}

func pairFastq(in io.Reader, out1, out2 io.Writer) error {
	r := fastq.NewReader(in, fastq.Phred33)
	w1 := fastq.NewWriter(out1, fastq.Phred33)
	w2 := fastq.NewWriter(out2, fastq.Phred33)
	idx := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("pair: read record %d: %w", idx+1, err)
		}
		if idx%2 == 0 {
			if err := w1.Write(rec); err != nil {
				return err
			}
		} else {
			if err := w2.Write(rec); err != nil {
				return err
			}
		}
		idx++
	}
	if idx%2 != 0 {
		return fmt.Errorf("pair: input has odd record count (%d): last record is unpaired", idx)
	}
	if err := w1.Flush(); err != nil {
		return err
	}
	return w2.Flush()
}

func pairFasta(in io.Reader, out1, out2 io.Writer) error {
	r := fasta.NewReader(in)
	w1 := fasta.NewWriter(out1, 0)
	w2 := fasta.NewWriter(out2, 0)
	idx := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("pair: read record %d: %w", idx+1, err)
		}
		if idx%2 == 0 {
			if err := w1.Write(rec); err != nil {
				return err
			}
		} else {
			if err := w2.Write(rec); err != nil {
				return err
			}
		}
		idx++
	}
	if idx%2 != 0 {
		return fmt.Errorf("pair: input has odd record count (%d): last record is unpaired", idx)
	}
	if err := w1.Flush(); err != nil {
		return err
	}
	return w2.Flush()
}
