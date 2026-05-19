// mergefa.go: implementation of `seqtk mergefa`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_mergefa (v1.5-r133).
// Behaviour: walks two FASTA/FASTQ streams in parallel and emits, for
// every paired position, a single base derived from the pair via
// upstream's seq_nt16 / seq_nt16_rev tables. The exact merge rule
// depends on the four boolean modes (`-i`, `-h`, `-m`, `-r`) and the
// quality threshold (`-q`). All five flags are real upstream surface
// (see the getopt string "himrq:" at seqtk.c:774).
//
// Output is FASTA, wrapped at 60 bases per line (upstream `l%60==0`).
// Case encodes confidence: an output base is uppercase only when both
// input bases are uppercase (or in the OR-modes when either is); below-
// quality bases are lowercased before the merge.
//
// Counters (same, diff, hom-het, het-hom, het-het) are emitted to
// stderr in the same format as upstream's fprintf at seqtk.c:868.
//
// Flag semantics, from the upstream getopt loop / inner switch:
//
//   -q INT  quality threshold; bases with PHRED+33 quality < q are
//           lowercased before merging (default 0).
//   -i      "intersect" mode: merged value is c[0] & c[1]; if it
//           collapses to 0, the result becomes 'x' (lowercase X) since
//           seq_nt16_rev_table[0] = 'X' and the case rule lowercases.
//   -m      "mask" mode: if either input is N (seq_nt16 == 15) the
//           result is lowercased; otherwise behave like -i.
//   -r      "randhet" mode: resolve hets using upstream's lrand48
//           coin-flips. Reproduces the upstream algorithm but uses
//           Go's math/rand (RNG-byte-parity is documented as
//           non-goal in docs/PARITY_ROADMAP.md#rng-policy).
//   -h      "haploid" mode: heterozygous merges become lowercase
//           (their case bit is cleared via is_upper = 0).
//
// -i and -m are mutually exclusive, matching upstream's early-exit
// check at seqtk.c:783.

package seqtk

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// MergefaOptions configures Mergefa. All zero-values match upstream's
// default behaviour (`seqtk mergefa <a.fa> <b.fa>` with no flags).
type MergefaOptions struct {
	// Quality is the PHRED+33 threshold below which a FASTQ base is
	// lowercased before merging. 0 disables the threshold check and
	// matches upstream's default. The upstream code reads
	// `seq->qual.s[l] - 33` regardless of the configured encoding;
	// non-FASTQ inputs have an empty qual buffer and the test is
	// skipped per-base.
	Quality int
	// Intersect enables `-i` (intersect mode).
	Intersect bool
	// Haploid enables `-h` (suppress hets via case lowering).
	Haploid bool
	// Mask enables `-m` (lowercase when one input is N).
	Mask bool
	// RandHet enables `-r` (random allele picking on hets).
	RandHet bool
	// Seed lets tests deterministically seed the RandHet RNG. 0
	// means "use upstream's srand48(11) constant" — chosen for
	// reproducibility within our tool (see PARITY_ROADMAP.md).
	Seed int64
}

// Mergefa reads two FASTA/FASTQ streams from r1 and r2 and writes the
// merged FASTA to w. Quality data from FASTQ inputs participates in
// the Quality threshold check; FASTQ output is never produced — the
// merge result is always FASTA, matching upstream.
//
// Mergefa returns an error if -i and -m are both set; otherwise it
// surfaces I/O / parsing errors and returns nil at clean EOF.
func Mergefa(r1, r2 io.Reader, w io.Writer, opts MergefaOptions) error {
	return mergefaImpl(r1, r2, w, os.Stderr, opts)
}

// mergefaRec is the small union used by the merge loop so FASTA and
// FASTQ inputs share a code path. qual is nil for FASTA.
type mergefaRec struct {
	name string
	seq  []byte
	qual []byte
}

// mergefaIter returns a closure that yields the next mergefaRec from
// the underlying FASTA / FASTQ stream, surfacing io.EOF unchanged.
func mergefaIter(r io.Reader) func() (*mergefaRec, error) {
	br, isFq := peekIsFastq(r)
	if isFq {
		// Upstream subtracts the literal "- 33" regardless of the
		// configured encoding, so we always read as Phred33.
		fr := fastq.NewReader(br, fastq.Phred33)
		return func() (*mergefaRec, error) {
			rec, err := fr.Read()
			if err != nil {
				return nil, err
			}
			return &mergefaRec{name: rec.ID, seq: rec.Sequence, qual: rec.Quality}, nil
		}
	}
	fr := fasta.NewReader(br)
	return func() (*mergefaRec, error) {
		rec, err := fr.Read()
		if err != nil {
			return nil, err
		}
		return &mergefaRec{name: rec.ID, seq: rec.Sequence}, nil
	}
}

// mergefaImpl is the testable core: lets tests capture the upstream
// stderr-style warnings and counter line on a custom writer.
func mergefaImpl(r1, r2 io.Reader, w, warn io.Writer, opts MergefaOptions) error {
	if opts.Mask && opts.Intersect {
		return errors.New("`-i' and `-m' cannot be applied at the same time")
	}
	it1 := mergefaIter(r1)
	it2 := mergefaIter(r2)

	seed := opts.Seed
	if seed == 0 {
		// Upstream uses srand48(11); mirror with a fixed Go seed
		// for reproducibility within our tool. See
		// docs/PARITY_ROADMAP.md#rng-policy.
		seed = 11
	}
	rng := rand.New(rand.NewSource(seed))

	bw := bufio.NewWriter(w)
	var cnt [5]uint64

	for {
		s1, err := it1()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		s2, err := it2()
		if err == io.EOF {
			// Upstream calls kseq_read on seq[1] without checking
			// the return value; an early EOF leaves seq[1] with
			// stale data. We mirror by zeroing out the mate so
			// min_l drops to 0 and the per-base loop is a no-op.
			fmt.Fprintf(warn, "[mergefa] second stream ended before first: %s has no mate record\n", s1.name)
			s2 = &mergefaRec{}
		} else if err != nil {
			return err
		}
		if s1.name != s2.name {
			fmt.Fprintf(warn, "[mergefa] Different sequence names: %s != %s\n", s1.name, s2.name)
		}
		if len(s1.seq) != len(s2.seq) {
			fmt.Fprintf(warn, "[mergefa] Unequal sequence length: %d != %d\n", len(s1.seq), len(s2.seq))
		}
		minL := len(s1.seq)
		if len(s2.seq) < minL {
			minL = len(s2.seq)
		}
		if _, err := fmt.Fprintf(bw, ">%s", s1.name); err != nil {
			return err
		}
		for l := 0; l < minL; l++ {
			c0 := s1.seq[l]
			c1 := s2.seq[l]
			// Quality lowering, matching upstream's per-side check
			// `if (seq[i]->qual.l && seq[i]->qual.s[l] - 33 < qual)`.
			if l < len(s1.qual) && int(s1.qual[l])-33 < opts.Quality {
				c0 = toLowerByte(c0)
			}
			if l < len(s2.qual) && int(s2.qual[l])-33 < opts.Quality {
				c1 = toLowerByte(c1)
			}
			var isUpper bool
			switch {
			case opts.Intersect, opts.Mask:
				isUpper = isASCIIUpper(c0) || isASCIIUpper(c1)
			default:
				isUpper = isASCIIUpper(c0) && isASCIIUpper(c1)
			}
			n0 := int(seqNT16Table[c0])
			n1 := int(seqNT16Table[c1])
			// Upstream: `if (c[0] == 0) c[0] = 15;` — X (seq_nt16
			// == 0) is widened to N for the merge math.
			if n0 == 0 {
				n0 = 15
			}
			if n1 == 0 {
				n1 = 15
			}
			b0 := int(bitCntTable[n0])
			b1 := int(bitCntTable[n1])
			if isUpper {
				switch {
				case b0 == 1 && b1 == 1:
					if n0 == n1 {
						cnt[0]++
					} else {
						cnt[1]++
					}
				case b0 == 1 && b1 == 2:
					cnt[2]++
				case b0 == 2 && b1 == 1:
					cnt[3]++
				case b0 == 2 && b1 == 2:
					cnt[4]++
				}
			}
			if opts.Haploid && (b0 > 1 || b1 > 1) {
				isUpper = false
			}
			var merged int
			switch {
			case opts.Intersect:
				merged = n0 & n1
				if merged == 0 {
					isUpper = false
				}
			case opts.Mask:
				if n0 == 15 || n1 == 15 {
					isUpper = false
				}
				merged = n0 & n1
				if merged == 0 {
					isUpper = false
				}
			case opts.RandHet:
				switch {
				case b0 == 1 && b1 == 1: // two homs
					merged = n0 | n1
				case (b0 == 1 && b1 == 2 || b0 == 2 && b1 == 1) && (n0&n1) != 0:
					// one hom, one het
					if rng.Int63()&1 == 1 {
						merged = n0 & n1
					} else {
						merged = n0 | n1
					}
				case b0 == 2 && b1 == 2 && n0 == n1:
					// double hets
					merged = n0
					if rng.Int63()&1 == 1 {
						if rng.Int63()&1 == 1 {
							// pick the "larger" allele
							for i := 8; i >= 1; i >>= 1 {
								if merged&i != 0 {
									merged &= i
								}
							}
						} else {
							// pick the "smaller" allele
							for i := 1; i <= 8; i <<= 1 {
								if merged&i != 0 {
									merged &= i
								}
							}
						}
					}
				default:
					merged = n0
					isUpper = false
				}
			default:
				merged = n0 | n1
			}
			out := seqNT16RevTable[merged&15]
			if !isUpper {
				out = toLowerByte(out)
			}
			if l%60 == 0 {
				if err := bw.WriteByte('\n'); err != nil {
					return err
				}
			}
			if err := bw.WriteByte(out); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	// Emit the upstream counter summary on the warn stream
	// (upstream uses stderr at seqtk.c:868).
	fmt.Fprintf(warn, "[mergefa] (same,diff,hom-het,het-hom,het-het)=(%d,%d,%d,%d,%d)\n",
		cnt[0], cnt[1], cnt[2], cnt[3], cnt[4])
	return bw.Flush()
}

// isASCIIUpper reports whether b is an ASCII uppercase letter (A-Z),
// matching C's `isupper` semantics for the ASCII subset used by
// upstream seqtk.
func isASCIIUpper(b byte) bool {
	return b >= 'A' && b <= 'Z'
}

// seqNT16RevTable is upstream's seq_nt16_rev_table: the inverse of
// seqNT16Table for the 16-state IUPAC alphabet. Indexing is by 4-bit
// IUPAC code (0..15); position 0 is 'X' (no base), position 15 is 'N'.
var seqNT16RevTable = [16]byte{
	'X', 'A', 'C', 'M',
	'G', 'R', 'S', 'V',
	'T', 'W', 'Y', 'H',
	'K', 'D', 'B', 'N',
}
