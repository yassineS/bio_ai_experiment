// hrun.go: implementation of `seqtk hrun`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_hrun (v1.5-r133,
// lines 1174-1204). For every FASTA record, walk the sequence and emit
// one BED4 row per maximal run of identical bytes whose length is at
// least minLen:
//
//	name\t<0-based start>\t<0-based end>\t<byte>\n
//
// where `end` = start + run-length (so the interval is half-open and
// `end - start` is the run length). The comparison is byte-exact —
// upstream does NOT case-fold and does NOT normalise IUPAC codes, so
// `AAaa` is two runs of length 2 (the lowercase `a` differs from the
// uppercase `A`).
//
// Upstream getopt surface: NONE. The first positional argument is the
// input file; an optional second positional argument overrides the
// minimum-length default (`if (argc == optind + 2) min_len =
// atoi(argv[optind+1]);`). The Go cmd/ layer exposes this knob as a
// flag named `-l/--min-len` for consistency with the rest of the
// project; the positional form is still accepted to match the upstream
// invocation pattern verbatim.
//
// Upstream quirk reproduced byte-for-byte: the "flush the open run on
// end-of-stream" `if (l >= min_len) printf(...);` lives OUTSIDE the
// kseq_read loop (line 1200). That means it fires AT MOST ONCE for the
// entire input, using the last record's name and the run state left
// over from that record (or, if the last record is empty, with the
// previous record's name and a degenerate l == 1 carried from
// `ks->seq.s[0]` UB — the upstream behaviour we observed empirically
// is that an empty trailing record swallows the would-be flush because
// the empty kseq leaves `l == 1 < min_len`). We mirror that quirk
// here: see the "trailing flush" comment below.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
)

// HrunOptions configures Hrun. MinLen is the minimum run length to
// report (upstream's optional positional argument, default 7).
type HrunOptions struct {
	MinLen int
}

// DefaultHrunMinLen is the upstream default for `seqtk hrun`
// (`int min_len = 7;` at seqtk.c:1178).
const DefaultHrunMinLen = 7

// Hrun walks every FASTA record in r and writes one BED4 row per
// maximal run of identical bytes of length >= opts.MinLen. See the
// package-level comment for the exact output format and the
// reproduced upstream "trailing flush" quirk.
//
// opts.MinLen <= 0 is silently treated like `>= 1`: every run is
// reported. (Upstream's `if (l >= min_len)` admits MinLen <= 0
// equivalently — every run has l >= 1.)
func Hrun(r io.Reader, w io.Writer, opts HrunOptions) error {
	bw := bufio.NewWriter(w)

	// Carry the loop-trailing run state across records, mirroring
	// upstream's `c, l, beg` (which are loop-local in C but their
	// values persist after the while-loop because they are declared
	// outside it — see seqtk.c:1179). After processing the last
	// non-empty record these hold the final open run; the post-loop
	// `if (l >= min_len)` reads them.
	var (
		lastName string
		c        byte
		l        int64 = 0
		beg      int64 = 0
	)
	// haveRecord tracks whether we've seen at least one record so we
	// can mirror upstream's behaviour at the trailing flush: when no
	// record was ever read, `ks->name.s` is NULL and upstream would
	// crash; here we simply skip the flush, which is the same
	// observable behaviour (no output) as an upstream invocation
	// with no records — empty input is a no-op upstream too because
	// the post-loop flush's `l` is 0 in that case.
	haveRecord := false

	err := scanFASTA(r, func(name string, seq []byte) error {
		lastName = name
		haveRecord = true
		if len(seq) == 0 {
			// Empty record: upstream reads `ks->seq.s[0]` (UB,
			// usually the NUL terminator) and sets l = 1, beg = 0.
			// The inner for loop doesn't execute. We replicate
			// that state to keep the trailing-flush quirk
			// byte-identical: c becomes 0 (NUL) and l becomes 1,
			// so the post-loop `if (l >= min_len)` is l(1) >=
			// min_len which is false for any min_len >= 2 (the
			// upstream default is 7) — meaning empty trailing
			// records silently swallow the trailing flush. This
			// matches what we observed by running the upstream
			// binary on a fixture ending with an empty record.
			c = 0
			l = 1
			beg = 0
			return nil
		}
		c = seq[0]
		l = 1
		beg = 0
		for i := 1; i < len(seq); i++ {
			if seq[i] != c {
				if l >= int64(opts.MinLen) {
					if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%c\n", name, beg, beg+l, c); err != nil {
						return err
					}
				}
				c = seq[i]
				l = 1
				beg = int64(i)
			} else {
				l++
			}
		}
		return nil
	})
	if err != nil {
		return err
	}

	// Trailing flush: upstream `if (l >= min_len) printf(...);` at
	// seqtk.c:1200, executed exactly once after the kseq_read loop
	// exits. We replicate the single-shot flush using the carried
	// state from the last record processed.
	if haveRecord && l >= int64(opts.MinLen) {
		if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%c\n", lastName, beg, beg+l, c); err != nil {
			return err
		}
	}
	return bw.Flush()
}
