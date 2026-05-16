// split.go: implementation of `seqtk split`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_split (v1.5-r133).
//
// Behaviour: round-robin every record from the input across N output
// files named "<prefix>.<5-digit, 1-based>.fa" (note the literal ".fa"
// suffix -- upstream uses it even for FASTQ input). Within each output
// file the input format is preserved (FASTA stays FASTA, FASTQ stays
// FASTQ) and sequence/quality lines are wrapped at lineLen characters
// when lineLen > 0; lineLen <= 0 keeps everything on a single line,
// matching upstream's stk_printseq(out[i%n], seq, len) call.
//
// Upstream flags (verified against the getopt loop at seqtk.c:1039):
//
//	-n INT    number of output files (default 10)
//	-l INT    line length for sequence/quality wrapping (default 0)
//
// Positional arguments: <prefix> <in.fa>  (the second may be "-" for
// stdin). The CLI layer handles input opening / decompression via
// seqtk.OpenInput; the output files themselves are plain (uncompressed)
// to match upstream byte-for-byte.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// SplitOptions configures Split.
type SplitOptions struct {
	// N is the number of output files. Upstream default is 10.
	N int
	// LineLen wraps sequence (and FASTQ quality) lines at LineLen
	// characters when > 0. 0 (the upstream default) emits each
	// sequence/quality on a single line.
	LineLen int
	// Prefix is the output-file name prefix; files are named
	// "<Prefix>.<5-digit 1-based index>.fa".
	Prefix string
}

// DefaultSplitN is the upstream default for `seqtk split -n` (see
// seqtk.c:1036).
const DefaultSplitN = 10

// Split reads every record from r and writes them round-robin across N
// output files (created on disk via os.Create). The first record goes
// to <prefix>.00001.fa, the second to <prefix>.00002.fa, ..., the
// (N+1)-th wraps back to <prefix>.00001.fa, and so on, matching
// upstream "seqtk split" byte-for-byte.
//
// The function takes care of creating, buffering and closing each
// output file. It returns the first I/O error it encounters; any
// partially-written files are still closed (but not removed -- this
// mirrors upstream behaviour, which also leaves them on disk).
func Split(r io.Reader, opts SplitOptions) error {
	if opts.N <= 0 {
		return fmt.Errorf("seqtk split: -n must be > 0, got %d", opts.N)
	}
	if opts.Prefix == "" {
		return fmt.Errorf("seqtk split: prefix is required")
	}

	// Open all N output files up front, matching upstream which
	// opens them all (and aborts on the first failure) before
	// reading any records.
	outs := make([]*bufio.Writer, opts.N)
	files := make([]*os.File, opts.N)
	closeAll := func() {
		for i := range outs {
			if outs[i] != nil {
				_ = outs[i].Flush()
			}
			if files[i] != nil {
				_ = files[i].Close()
			}
		}
	}
	for i := 0; i < opts.N; i++ {
		name := fmt.Sprintf("%s.%05d.fa", opts.Prefix, i+1)
		f, err := os.Create(name)
		if err != nil {
			closeAll()
			return fmt.Errorf("seqtk split: failed to create %s: %w", name, err)
		}
		files[i] = f
		outs[i] = bufio.NewWriter(f)
	}

	br, isFastq := peekIsFastq(r)
	var i int
	err := func() error {
		if isFastq {
			fr := fastq.NewReader(br, fastq.Phred33)
			for {
				rec, err := fr.Read()
				if err == io.EOF {
					return nil
				}
				if err != nil {
					return err
				}
				if err := writeSplitFastq(outs[i%opts.N], rec, opts.LineLen); err != nil {
					return err
				}
				i++
			}
		}
		fr := fasta.NewReader(br)
		for {
			rec, err := fr.Read()
			if err == io.EOF {
				return nil
			}
			if err != nil {
				return err
			}
			if err := writeSplitFasta(outs[i%opts.N], rec, opts.LineLen); err != nil {
				return err
			}
			i++
		}
	}()
	closeAll()
	return err
}

// writeSplitFasta emits one FASTA record to bw, wrapping the sequence
// at lineLen characters when lineLen > 0 (1:1 port of upstream
// stk_printstr's wrap behaviour).
func writeSplitFasta(bw *bufio.Writer, rec *fasta.Record, lineLen int) error {
	if err := bw.WriteByte('>'); err != nil {
		return err
	}
	if _, err := bw.WriteString(rec.Description); err != nil {
		return err
	}
	if err := bw.WriteByte('\n'); err != nil {
		return err
	}
	return writeWrappedBytes(bw, rec.Sequence, lineLen)
}

// writeSplitFastq emits one FASTQ record to bw, wrapping both the
// sequence and the quality lines at lineLen characters when lineLen
// > 0 (1:1 port of upstream stk_printseq, which wraps both via
// stk_printstr).
func writeSplitFastq(bw *bufio.Writer, rec *fastq.Record, lineLen int) error {
	if err := bw.WriteByte('@'); err != nil {
		return err
	}
	if _, err := bw.WriteString(rec.Description); err != nil {
		return err
	}
	if err := bw.WriteByte('\n'); err != nil {
		return err
	}
	if err := writeWrappedBytes(bw, rec.Sequence, lineLen); err != nil {
		return err
	}
	if _, err := bw.WriteString("+\n"); err != nil {
		return err
	}
	return writeWrappedBytes(bw, rec.Quality, lineLen)
}

// writeWrappedBytes writes b to bw followed by a trailing '\n'. When
// lineLen > 0 the bytes are split into successive lineLen-byte lines
// (one '\n' per chunk); when lineLen <= 0 b is written as a single
// line. Matches upstream stk_printstr (seqtk.c:237) byte-for-byte:
// our callers always emit the '\n' that follows the header themselves
// (whereas upstream emits it as the first byte inside stk_printstr),
// so for empty input we add no extra '\n' here in either branch --
// upstream's empty-string output is exactly one trailing '\n' in both
// branches, which our callers have already written.
func writeWrappedBytes(bw *bufio.Writer, b []byte, lineLen int) error {
	if lineLen <= 0 {
		if _, err := bw.Write(b); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	for off := 0; off < len(b); off += lineLen {
		end := off + lineLen
		if end > len(b) {
			end = len(b)
		}
		if _, err := bw.Write(b[off:end]); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	return nil
}
