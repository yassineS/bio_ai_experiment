// compression.go: subcommands that compress sequence content.
//
// Currently this file holds the HPC (homopolymer-compression) implementation
// for the seqtk hpc subcommand.

package seqtk

import (
	"bufio"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// HPC reads a FASTA (or FASTQ — sequence-only) stream from in and writes
// homopolymer-compressed FASTA records to w. Every maximal run of identical
// bytes in a sequence is collapsed to a single byte; runs of length 1 are
// preserved as one byte. The record name and (where present) the
// description are preserved on output; the compressed sequence is emitted on
// a single line with no wrapping, matching the upstream "seqtk hpc" output.
//
// Empty input sequences are skipped (no record emitted), matching upstream.
func HPC(in io.Reader, w io.Writer) error {
	br, isFastq := peekIsFastq(in)
	bw := bufio.NewWriter(w)

	emit := func(header string, seq []byte) error {
		if len(seq) == 0 {
			return nil
		}
		compressed := CompressHomopolymers(seq)
		if _, err := bw.WriteString(">"); err != nil {
			return err
		}
		if _, err := bw.WriteString(header); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if _, err := bw.Write(compressed); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}

	if isFastq {
		// We still support FASTQ input by reading sequences only; the
		// output is always FASTA.
		// Re-use the FASTQ reader via existing helpers.
		return hpcFromFastq(br, bw, emit)
	}

	reader := fasta.NewReader(br)
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if err := emit(rec.Description, rec.Sequence); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// hpcFromFastq reads FASTQ records from br and emits homopolymer-compressed
// FASTA records via emit, then flushes bw.
func hpcFromFastq(br io.Reader, bw *bufio.Writer, emit func(header string, seq []byte) error) error {
	r := fastq.NewReader(br, fastq.Phred33)
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if err := emit(rec.Description, rec.Sequence); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// CompressHomopolymers returns a new byte slice in which every maximal run of
// identical bytes in seq has been collapsed to a single byte. The first byte
// of each run is the one kept (so the result preserves the case present at
// the start of each run). An empty input yields an empty result.
//
// Examples:
//
//	CompressHomopolymers([]byte("AAACCGT"))   == []byte("ACGT")
//	CompressHomopolymers([]byte("aaAACCt"))   == []byte("aACt")
//	CompressHomopolymers([]byte(""))          == []byte{}
func CompressHomopolymers(seq []byte) []byte {
	if len(seq) == 0 {
		return []byte{}
	}
	// Worst case (no runs) the result is the same length as the input.
	out := make([]byte, 0, len(seq))
	out = append(out, seq[0])
	for i := 1; i < len(seq); i++ {
		if seq[i] != seq[i-1] {
			out = append(out, seq[i])
		}
	}
	return out
}
