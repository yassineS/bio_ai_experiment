// sample.go: an upstream-faithful port of `seqtk sample`, including the krand
// 64-bit Mersenne-Twister RNG it uses, so the seeded subsampling is
// byte-for-byte identical to stk_sample in reference_code/seqtk/seqtk.c.
//
// The high-level Sample helper in seqtk.go uses a deterministic every-Nth
// selection that never matched upstream (its parity test is skipped). SampleN
// reproduces upstream's reservoir sampler (and, with -2, its two-pass mode).

package seqtk

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// krand is a 64-bit Mersenne-Twister (MT19937-64) PRNG, a direct port of
// upstream seqtk's krand_t. It is used so that SampleN's seeded output matches
// upstream's stk_sample byte-for-byte.
type krand struct {
	mti int
	mt  [krNN]uint64
}

const (
	krNN = 312
	krMM = 156
	krUM = 0xFFFFFFFF80000000 // most significant 33 bits
	krLM = 0x7FFFFFFF         // least significant 31 bits
)

// newKrand seeds a krand exactly as upstream's kr_srand / kr_srand0.
func newKrand(seed uint64) *krand {
	kr := &krand{}
	kr.mt[0] = seed
	for kr.mti = 1; kr.mti < krNN; kr.mti++ {
		kr.mt[kr.mti] = 6364136223846793005*(kr.mt[kr.mti-1]^(kr.mt[kr.mti-1]>>62)) + uint64(kr.mti)
	}
	return kr
}

// rand returns the next 64-bit value, mirroring upstream's kr_rand.
func (kr *krand) rand() uint64 {
	mag01 := [2]uint64{0, 0xB5026F5AA96619E9}
	var x uint64
	if kr.mti >= krNN {
		var i int
		if kr.mti == krNN+1 {
			*kr = *newKrand(5489)
		}
		for i = 0; i < krNN-krMM; i++ {
			x = (kr.mt[i] & krUM) | (kr.mt[i+1] & krLM)
			kr.mt[i] = kr.mt[i+krMM] ^ (x >> 1) ^ mag01[x&1]
		}
		for ; i < krNN-1; i++ {
			x = (kr.mt[i] & krUM) | (kr.mt[i+1] & krLM)
			kr.mt[i] = kr.mt[i+(krMM-krNN)] ^ (x >> 1) ^ mag01[x&1]
		}
		x = (kr.mt[krNN-1] & krUM) | (kr.mt[0] & krLM)
		kr.mt[krNN-1] = kr.mt[krMM-1] ^ (x >> 1) ^ mag01[x&1]
		kr.mti = 0
	}
	x = kr.mt[kr.mti]
	kr.mti++
	x ^= (x >> 29) & 0x5555555555555555
	x ^= (x << 17) & 0x71D67FFFEDA60000
	x ^= (x << 37) & 0xFFF7EEE000000000
	x ^= x >> 43
	return x
}

// drand returns a float64 in [0,1), mirroring upstream's kr_drand macro.
func (kr *krand) drand() float64 {
	return float64(kr.rand()>>11) * (1.0 / 9007199254740992.0)
}

// SampleN draws a fixed NUMBER of records from in using upstream's seeded
// reservoir sampler and writes them to w, matching `seqtk sample [-2] -s SEED
// <in> <num>` byte-for-byte. When twoPass is true the second (memory-frugal)
// pass over reader2 is used and output is in input order; otherwise the
// streaming reservoir is used and output is in reservoir-slot order. The
// reader2 factory is only consulted when twoPass is true (it must re-open the
// same input from the start).
func SampleN(in io.Reader, w io.Writer, num uint64, seed int64, twoPass bool, reopen func() (io.ReadCloser, error)) error {
	kr := newKrand(uint64(seed))
	bw := bufio.NewWriter(w)

	if !twoPass {
		// Streaming reservoir: buffer up to num records, replacing slots.
		buf := make([]*seqRecord, num)
		var nSeqs uint64
		if err := forEachSeqRecord(in, func(rec *seqRecord) error {
			r := kr.drand()
			nSeqs++
			var y uint64
			if nSeqs-1 < num {
				y = nSeqs - 1
			} else {
				y = uint64(r * float64(nSeqs))
			}
			if y < num {
				buf[y] = cloneSeqRecord(rec)
			}
			return nil
		}); err != nil {
			return err
		}
		for _, rec := range buf {
			if rec != nil {
				if err := writeSampleRecord(bw, rec); err != nil {
					return err
				}
			}
		}
		return bw.Flush()
	}

	// Two-pass mode: 1st pass picks the reservoir of record indices, 2nd pass
	// re-reads and emits the chosen records in input order.
	selected := make(map[uint64]bool, num)
	{
		buf := make([]uint64, num)
		for i := range buf {
			buf[i] = ^uint64(0)
		}
		var nSeqs uint64
		if err := forEachSeqRecord(in, func(_ *seqRecord) error {
			r := kr.drand()
			nSeqs++
			var y uint64
			if nSeqs-1 < num {
				y = nSeqs - 1
			} else {
				y = uint64(r * float64(nSeqs))
			}
			if y < num {
				buf[y] = nSeqs
			}
			return nil
		}); err != nil {
			return err
		}
		for _, v := range buf {
			selected[v] = true
		}
	}

	rc, err := reopen()
	if err != nil {
		return err
	}
	defer rc.Close()
	var nSeqs uint64
	if err := forEachSeqRecord(rc, func(rec *seqRecord) error {
		nSeqs++
		if selected[nSeqs] {
			return writeSampleRecord(bw, rec)
		}
		return nil
	}); err != nil {
		return err
	}
	return bw.Flush()
}

// SampleFraction draws each record independently with probability frac using
// upstream's seeded RNG, matching `seqtk sample -s SEED <in> <frac>` for
// 0 < frac < 1. Output preserves input order (upstream streams it directly).
func SampleFraction(in io.Reader, w io.Writer, frac float64, seed int64) error {
	kr := newKrand(uint64(seed))
	bw := bufio.NewWriter(w)
	if err := forEachSeqRecord(in, func(rec *seqRecord) error {
		if kr.drand() < frac {
			return writeSampleRecord(bw, rec)
		}
		return nil
	}); err != nil {
		return err
	}
	return bw.Flush()
}

// cloneSeqRecord deep-copies a seqRecord so it survives the reservoir.
func cloneSeqRecord(rec *seqRecord) *seqRecord {
	c := &seqRecord{name: rec.name, comment: rec.comment}
	c.seq = append([]byte(nil), rec.seq...)
	if rec.qual != nil {
		c.qual = append([]byte(nil), rec.qual...)
	}
	return c
}

// writeSampleRecord writes a record the way upstream's stk_printseq does for
// sample: no line wrapping (UINT_MAX line length).
func writeSampleRecord(bw *bufio.Writer, rec *seqRecord) error {
	return writeSeqRecord(bw, rec, 0)
}

// forEachSeqRecord reads a FASTA/FASTQ stream (auto-detected) and invokes fn
// for every record, splitting the header into name/comment like upstream.
func forEachSeqRecord(in io.Reader, fn func(*seqRecord) error) error {
	br, isFastq := peekIsFastq(in)
	if isFastq {
		r := fastq.NewReader(br, fastq.Phred33)
		for {
			fr, err := r.Read()
			if err == io.EOF {
				return nil
			}
			if err != nil {
				return err
			}
			name, comment := splitNameComment(fr.Description)
			if err := fn(&seqRecord{name: name, comment: comment, seq: fr.Sequence, qual: fr.Quality}); err != nil {
				return err
			}
		}
	}
	r := fasta.NewReader(br)
	for {
		fr, err := r.Read()
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		name, comment := splitNameComment(fr.Description)
		if err := fn(&seqRecord{name: name, comment: comment, seq: fr.Sequence}); err != nil {
			return err
		}
	}
}

// ErrTwoPassStdin is returned when -2 is requested with stdin input, which
// upstream rejects because it cannot rewind a stream.
var ErrTwoPassStdin = fmt.Errorf("sample: in the 2-pass mode, the input cannot be stdin")
