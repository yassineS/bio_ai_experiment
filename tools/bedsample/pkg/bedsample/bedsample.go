// Package bedsample implements `bedtools sample`: it draws N random records
// without replacement from a BED-like input and writes them to the output in
// the relative order they appeared in the file.
//
// The sampling is reservoir-style: a single linear pass over the input,
// O(N) memory. The seed is configurable (`-seed`) for deterministic
// replays; when seed=0 we fall back to a time-based seed.
//
// Mirrors upstream `bedtools sample -i FILE -n N [-seed SEED] [-header]`
// byte-for-byte for a given seed: the reservoir replacement uses a Go port
// of the same std::mt19937_64 engine and rejection-sampling bound upstream
// uses (see mt19937.go), and the kept records are emitted in reservoir-slot
// order — exactly as upstream's `giveFinalReport` does for BED output (it
// only re-sorts when the output type is BAM). Notes vs. upstream:
//
//   - The reservoir guarantees uniform sampling without replacement when N
//     <= total records. When N > total records we emit every record and
//     return an error from the library so the CLI wrapper can mirror
//     upstream's "Input file has fewer records than the requested number"
//     diagnostic.
package bedsample

import (
	"bufio"
	"fmt"
	"io"
	"strings"
	"time"
)

// Options configures Sample.
type Options struct {
	// N is the number of records to draw (without replacement).
	N int
	// Seed seeds the deterministic PRNG. When 0, Sample uses time-based seed.
	Seed int64
	// Header, when true, forwards header lines ('#', 'track', 'browser')
	// from the input to the output verbatim *before* the sampled records.
	// When false (default), header lines are dropped. Mirrors upstream's
	// `-header` flag.
	Header bool
}

// ErrTooFewRecords is returned when the input has fewer data records than
// the requested N. Mirrors upstream's "Input file has fewer records than the
// requested number of output records" error.
type ErrTooFewRecords struct {
	Have int
	Want int
}

func (e *ErrTooFewRecords) Error() string {
	return fmt.Sprintf("Input file has fewer records than the requested number of output records (have=%d want=%d)", e.Have, e.Want)
}

// Sample reads BED-like records from r, draws opts.N at random (using
// reservoir sampling), and writes the result to w. Returns the number of
// records written.
func Sample(r io.Reader, w io.Writer, opts Options) (int, error) {
	if opts.N <= 0 {
		return 0, fmt.Errorf("-n must be > 0 (got %d)", opts.N)
	}

	seed := opts.Seed
	if seed == 0 {
		seed = time.Now().UnixNano()
	}
	// Upstream seeds std::mt19937_64 with the (int) seed value; mirror that
	// by reinterpreting the seed's low bits as the unsigned engine seed.
	rng := newMT19937_64(uint64(seed))

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	// reservoir holds up to N kept lines in upstream's slot order: a new
	// record either fills the next empty slot (during the fill phase) or
	// replaces an existing slot chosen by the RNG. The final output is this
	// slice in slot order — no re-sort, matching upstream's BED path.
	reservoir := make([]string, 0, opts.N)

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	// total mirrors upstream's _currRecordNum: it counts every data record
	// seen so far, and is incremented BEFORE the replacement draw so the
	// draw range is [0, total) for the total-th record (1-based).
	total := 0
	for sc.Scan() {
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			if opts.Header {
				if _, err := bw.WriteString(raw); err != nil {
					return 0, err
				}
				if err := bw.WriteByte('\n'); err != nil {
					return 0, err
				}
			}
			continue
		}
		if trimmed == "" {
			// Drop blank lines silently — matches upstream behaviour.
			continue
		}

		total++
		if len(reservoir) < opts.N {
			// Fill phase: no RNG draw, exactly as upstream's keepRecord.
			reservoir = append(reservoir, raw)
			continue
		}
		// Replacement phase: draw an index in [0, total) and replace that
		// slot when it lands inside the reservoir.
		idx := rng.randRange(uint64(total))
		if idx < uint64(opts.N) {
			reservoir[idx] = raw
		}
	}
	if err := sc.Err(); err != nil {
		return 0, err
	}

	if total < opts.N {
		return 0, &ErrTooFewRecords{Have: total, Want: opts.N}
	}

	// Emit in reservoir-slot order (upstream does not re-sort BED output).
	for _, line := range reservoir {
		if _, err := bw.WriteString(line); err != nil {
			return 0, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return 0, err
		}
	}
	return len(reservoir), nil
}
