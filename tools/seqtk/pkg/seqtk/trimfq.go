// trimfq.go: an upstream-faithful port of `seqtk trimfq`, which trims FASTQ
// reads with a modified Mott algorithm (stk_trimfq in
// reference_code/seqtk/seqtk.c) or, when -b/-e/-L are given, by fixed offsets.
//
// The pre-existing TrimQuality helper in seqtk.go does a plain Phred-threshold
// trim and never matched upstream (its parity test is skipped). TrimFQ replaces
// that path for the CLI and reproduces stk_trimfq byte-for-byte, including:
//
//	-q FLOAT  error-rate threshold for the Mott algorithm [0.05]
//	-l INT    maximally trim down to INT bp (Mott path) [30]
//	-b INT    trim INT bp from the left  (non-zero disables -q/-l)
//	-e INT    trim INT bp from the right (non-zero disables -q/-l)
//	-L INT    retain at most INT bp from the 5'-end (disables -q/-l)

package seqtk

import (
	"bufio"
	"io"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TrimFQOptions holds the option tail for TrimFQ, mirroring stk_trimfq.
type TrimFQOptions struct {
	ErrorRate float64 // -q: error-rate threshold for the Mott algorithm
	MinLen    int     // -l: maximally trim down to this many bp (Mott path)
	Left      int     // -b: trim this many bp from the left
	Right     int     // -e: trim this many bp from the right
	FixedLen  int     // -L: retain at most this many bp from the 5'-end (-1 = off)
}

// DefaultTrimFQOptions returns TrimFQ's upstream defaults.
func DefaultTrimFQOptions() TrimFQOptions {
	return TrimFQOptions{ErrorRate: 0.05, MinLen: 30, FixedLen: -1}
}

// TrimFQ reads a FASTA/FASTQ stream from in (auto-detected) and writes the
// trimmed records to w, reproducing upstream's stk_trimfq. FASTA input (or a
// zero-length quality record) takes the fixed-offset path only.
func TrimFQ(in io.Reader, w io.Writer, opts TrimFQOptions) error {
	// q_int2real[i] = 10^(-(i-33)/10): the error probability for Phred i.
	var qInt2Real [128]float64
	for i := 0; i < 128; i++ {
		qInt2Real[i] = math.Pow(10, -float64(i-33)/10)
	}

	br, isFastq := peekIsFastq(in)
	bw := bufio.NewWriter(w)

	trim := func(rec *seqRecord) error {
		ql := len(rec.qual)
		var beg, end int
		if opts.Left != 0 || opts.Right != 0 || opts.FixedLen > 0 {
			// Fixed-offset path (disables the Mott algorithm).
			beg = opts.Left
			end = len(rec.seq) - opts.Right
			if beg >= end {
				beg, end = 0, 0
			}
			if opts.FixedLen > 0 && end-beg > opts.FixedLen {
				end = beg + opts.FixedLen
			}
		} else if ql > opts.MinLen {
			// Modified Mott algorithm.
			var s, max float64
			tmp := 0
			beg, end = 0, ql
			for i := 0; i < ql; i++ {
				q := int(rec.qual[i])
				if q < 36 {
					q = 36
				}
				if q > 127 {
					q = 127
				}
				s += opts.ErrorRate - qInt2Real[q]
				if s > max {
					max = s
					beg = tmp
					end = i + 1
				}
				if s < 0 {
					s = 0
					tmp = i + 1
				}
			}
			if max == 0 {
				// All low quality: keep the first MinLen bp.
				beg, end = 0, opts.MinLen
			}
			if end-beg < opts.MinLen {
				// Window-based fallback: pick the highest-sum MinLen window.
				is := 0
				for i := 0; i < opts.MinLen; i++ {
					is += int(rec.qual[i]) - 33
				}
				imax := is
				beg = 0
				for i := opts.MinLen; i < ql; i++ {
					is += int(rec.qual[i]) - int(rec.qual[i-opts.MinLen])
					if imax < is {
						imax = is
						beg = i - opts.MinLen + 1
					}
				}
				end = beg + opts.MinLen
			}
		} else {
			beg, end = 0, len(rec.seq)
		}

		if end > len(rec.seq) {
			end = len(rec.seq)
		}
		if beg > end {
			beg = end
		}

		marker := byte('>')
		if isFastq {
			marker = '@'
		}
		if err := bw.WriteByte(marker); err != nil {
			return err
		}
		if _, err := bw.WriteString(rec.name); err != nil {
			return err
		}
		if rec.comment != "" {
			if err := bw.WriteByte(' '); err != nil {
				return err
			}
			if _, err := bw.WriteString(rec.comment); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if _, err := bw.Write(rec.seq[beg:end]); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if isFastq {
			if _, err := bw.WriteString("+\n"); err != nil {
				return err
			}
			if _, err := bw.Write(rec.qual[beg:end]); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
		}
		return nil
	}

	if isFastq {
		r := fastq.NewReader(br, fastq.Phred33)
		for {
			fr, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			name, comment := splitNameComment(fr.Description)
			if err := trim(&seqRecord{name: name, comment: comment, seq: fr.Sequence, qual: fr.Quality}); err != nil {
				return err
			}
		}
	} else {
		r := fasta.NewReader(br)
		for {
			fr, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			name, comment := splitNameComment(fr.Description)
			if err := trim(&seqRecord{name: name, comment: comment, seq: fr.Sequence}); err != nil {
				return err
			}
		}
	}
	return bw.Flush()
}
