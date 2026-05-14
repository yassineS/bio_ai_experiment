// Package bedsample implements `bedtools sample`: it draws N random records
// without replacement from a BED-like input and writes them to the output in
// the relative order they appeared in the file.
//
// The sampling is reservoir-style: a single linear pass over the input,
// O(N) memory. The seed is configurable (`-seed`) for deterministic
// replays; when seed=0 we fall back to a time-based seed.
//
// Mirrors upstream `bedtools sample -i FILE -n N [-seed SEED] [-header]`.
// Notes vs. upstream:
//
//   - Output order is the input file order (preserves the first-appearance
//     ordering of sampled records). Upstream behaves the same way in
//     practice — it never re-sorts.
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
	"math/rand"
	"sort"
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
	rng := rand.New(rand.NewSource(seed))

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	// reservoirSlot tracks the original 0-based index of each sampled
	// record so we can emit them in input order.
	type slot struct {
		index int
		line  string
	}
	reservoir := make([]slot, 0, opts.N)

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
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

		if len(reservoir) < opts.N {
			reservoir = append(reservoir, slot{index: total, line: raw})
		} else {
			// Pick a random index in [0, total]. If it falls inside the
			// reservoir, replace that slot.
			j := rng.Intn(total + 1)
			if j < opts.N {
				reservoir[j] = slot{index: total, line: raw}
			}
		}
		total++
	}
	if err := sc.Err(); err != nil {
		return 0, err
	}

	if total < opts.N {
		return 0, &ErrTooFewRecords{Have: total, Want: opts.N}
	}

	// Emit in input order (sorted by `index`).
	sort.Slice(reservoir, func(i, j int) bool {
		return reservoir[i].index < reservoir[j].index
	})
	for _, s := range reservoir {
		if _, err := bw.WriteString(s.line); err != nil {
			return 0, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return 0, err
		}
	}
	return len(reservoir), nil
}
