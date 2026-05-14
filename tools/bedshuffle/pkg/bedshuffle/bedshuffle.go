// Package bedshuffle randomly relocates BED intervals across a genome,
// mirroring `bedtools shuffle`.
//
// For each input interval the package draws a new (chrom, start) coordinate
// uniformly at random and emits a record with the same length, keeping name /
// score / strand / extras intact. Various flags constrain the draw:
//
//   - Include: restrict placements to regions listed in an "include" BED.
//   - Exclude: forbid placements that overlap any region listed in an
//     "exclude" BED.
//   - Chrom: keep each shuffled interval on its original chromosome.
//
// Placement uses a deterministic *math/rand.Rand seeded by Seed. Each draw
// is retried up to MaxRetries (default 1000) before the interval is reported
// as unplaceable and an error is returned with a clear message.
package bedshuffle

import (
	"bufio"
	"fmt"
	"io"
	"math/rand"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options controls Shuffle.
type Options struct {
	// Genome is a tab-delimited chromosome-size map; key is the chrom name,
	// value is the chromosome length in bp.
	Genome map[string]int

	// Include limits placements to regions listed here (per chrom).
	// nil means "place anywhere on the chromosome".
	Include []*bed.Record
	// Exclude rejects placements that overlap any region here.
	Exclude []*bed.Record

	// Chrom: when true, each shuffled interval stays on its original
	// chromosome.
	Chrom bool

	// Seed for the RNG. Same seed + same inputs ⇒ same output.
	Seed int64

	// MaxRetries caps the number of placement attempts per interval. 0 ⇒
	// default of 1000 (mirrors upstream).
	MaxRetries int
}

// Shuffle reads BED records from r, draws new coordinates for each according
// to opts, and writes the relocated records to w. Returns the number of
// records placed; the first interval that cannot be placed within
// MaxRetries triggers an error after the records placed so far have been
// flushed.
func Shuffle(r io.Reader, w io.Writer, opts Options) (int, error) {
	if opts.Genome == nil || len(opts.Genome) == 0 {
		return 0, fmt.Errorf("shuffle: -g (genome file) is required")
	}
	maxRetries := opts.MaxRetries
	if maxRetries == 0 {
		maxRetries = 1000
	}

	rng := rand.New(rand.NewSource(opts.Seed))

	// Index include/exclude into per-chrom sorted slices so we can sample.
	include := indexByChrom(opts.Include)
	exclude := indexByChrom(opts.Exclude)

	// Precompute include offsets per chrom for uniform sampling: a list of
	// (start, end) pairs and the cumulative bp count.
	includeRegions := map[string][]region{}
	includeCum := map[string]int{}
	for chrom, recs := range include {
		var regions []region
		total := 0
		for _, rec := range recs {
			start := rec.ChromStart
			end := rec.ChromEnd
			if size, ok := opts.Genome[chrom]; ok && end > size {
				end = size
			}
			if start < 0 {
				start = 0
			}
			if end > start {
				regions = append(regions, region{start: start, end: end, cum: total})
				total += end - start
			}
		}
		includeRegions[chrom] = regions
		includeCum[chrom] = total
	}

	// Cached sorted chrom list for "any-chrom" weighted sampling.
	weightedChroms, chromCum, totalGenome := genomeWeighting(opts.Genome)
	weightedIncl, inclCum, totalIncl := includeWeighting(includeCum)

	// Sort exclude records for fast overlap checks (linear scan with early
	// break works for small lists; if perf matters we can swap in an
	// interval tree later).
	for chrom := range exclude {
		sort.SliceStable(exclude[chrom], func(i, j int) bool {
			return exclude[chrom][i].ChromStart < exclude[chrom][j].ChromStart
		})
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	count := 0
	for scanner.Scan() {
		line := scanner.Text()
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return count, fmt.Errorf("BED record must have at least 3 columns, got %d", len(fields))
		}
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return count, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return count, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		length := end - start
		if length <= 0 {
			return count, fmt.Errorf("interval %s:%d-%d has non-positive length", fields[0], start, end)
		}

		newChrom, newStart, ok := drawPlacement(rng, fields[0], length, opts,
			includeRegions, weightedChroms, chromCum, totalGenome,
			weightedIncl, inclCum, totalIncl,
			exclude, maxRetries)
		if !ok {
			return count, fmt.Errorf("Error, line %d: tried %d potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on.",
				count+1, maxRetries)
		}
		// Write the new record: replace chrom/start/end columns, keep the
		// rest verbatim.
		out := append([]string(nil), fields...)
		out[0] = newChrom
		out[1] = strconv.Itoa(newStart)
		out[2] = strconv.Itoa(newStart + length)
		if _, err := fmt.Fprintln(bw, strings.Join(out, "\t")); err != nil {
			return count, err
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return count, err
	}
	return count, nil
}

// region holds a half-open [start, end) range and its cumulative bp offset
// inside an include list.
type region struct {
	start, end, cum int
}

// drawPlacement returns a (chrom, start, ok) tuple for the next placement.
// It honours Chrom, Include, and Exclude. Returns ok=false after maxRetries.
func drawPlacement(
	rng *rand.Rand,
	origChrom string,
	length int,
	opts Options,
	includeRegions map[string][]region,
	weightedChroms []string, chromCum []int, totalGenome int,
	weightedIncl []string, inclCum []int, totalIncl int,
	exclude map[string][]*bed.Record,
	maxRetries int,
) (string, int, bool) {
	for try := 0; try < maxRetries; try++ {
		chrom, start, ok := drawOne(rng, origChrom, length, opts,
			includeRegions, weightedChroms, chromCum, totalGenome,
			weightedIncl, inclCum, totalIncl)
		if !ok {
			continue
		}
		// Honour Exclude.
		if hasOverlap(exclude[chrom], start, start+length) {
			continue
		}
		return chrom, start, true
	}
	return "", 0, false
}

// drawOne performs a single weighted draw. Returns ok=false if the chrom has
// no room or the include list is empty for that chrom.
func drawOne(
	rng *rand.Rand,
	origChrom string,
	length int,
	opts Options,
	includeRegions map[string][]region,
	weightedChroms []string, chromCum []int, totalGenome int,
	weightedIncl []string, inclCum []int, totalIncl int,
) (string, int, bool) {
	if len(opts.Include) > 0 {
		// Include mode: sample uniformly over the include intervals.
		var chrom string
		if opts.Chrom {
			chrom = origChrom
		} else {
			if totalIncl <= 0 {
				return "", 0, false
			}
			chrom = pickByWeight(rng, weightedIncl, inclCum, totalIncl)
		}
		regs := includeRegions[chrom]
		if len(regs) == 0 {
			return "", 0, false
		}
		// Filter regions that fit the interval length.
		var fit []region
		for _, r := range regs {
			if r.end-r.start >= length {
				fit = append(fit, r)
			}
		}
		if len(fit) == 0 {
			return "", 0, false
		}
		// Uniform over bp inside the fitting regions.
		totalFit := 0
		offsets := make([]int, len(fit))
		for i, r := range fit {
			offsets[i] = totalFit
			totalFit += (r.end - r.start - length + 1)
		}
		if totalFit <= 0 {
			return "", 0, false
		}
		bp := rng.Intn(totalFit)
		// Find which fit region this falls into via linear scan; len(fit) is
		// small per chrom so this is fine.
		idx := sort.Search(len(offsets), func(i int) bool { return offsets[i] > bp })
		if idx == 0 {
			idx = 1
		}
		idx--
		r := fit[idx]
		start := r.start + (bp - offsets[idx])
		if start+length > r.end {
			return "", 0, false
		}
		return chrom, start, true
	}

	// No include: sample uniformly across the genome (per-chrom weighted).
	var chrom string
	if opts.Chrom {
		chrom = origChrom
	} else {
		if totalGenome <= 0 {
			return "", 0, false
		}
		chrom = pickByWeight(rng, weightedChroms, chromCum, totalGenome)
	}
	size, ok := opts.Genome[chrom]
	if !ok || size < length {
		return "", 0, false
	}
	start := rng.Intn(size - length + 1)
	return chrom, start, true
}

// pickByWeight selects a chromosome with probability proportional to its
// genome size. `keys` is the chrom-name list, `cum` is the cumulative-bp
// table aligned with keys, total is the sum of all weights.
func pickByWeight(rng *rand.Rand, keys []string, cum []int, total int) string {
	v := rng.Intn(total)
	// Find first cum[i] > v.
	idx := sort.SearchInts(cum, v+1)
	if idx >= len(keys) {
		idx = len(keys) - 1
	}
	return keys[idx]
}

// indexByChrom buckets records by chromosome.
func indexByChrom(recs []*bed.Record) map[string][]*bed.Record {
	out := map[string][]*bed.Record{}
	for _, r := range recs {
		out[r.Chrom] = append(out[r.Chrom], r)
	}
	return out
}

// hasOverlap reports whether [start, end) overlaps any of recs (linear scan).
func hasOverlap(recs []*bed.Record, start, end int) bool {
	for _, r := range recs {
		if r.ChromStart < end && r.ChromEnd > start {
			return true
		}
	}
	return false
}

// genomeWeighting returns (chrom-name list sorted, cumulative bp at each,
// total genome size) for weighted chrom selection.
func genomeWeighting(g map[string]int) ([]string, []int, int) {
	keys := make([]string, 0, len(g))
	for k := range g {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	cum := make([]int, len(keys))
	total := 0
	for i, k := range keys {
		total += g[k]
		cum[i] = total
	}
	return keys, cum, total
}

// includeWeighting is the analogous helper for the include list: a chrom is
// weighted by the total bp listed for it in the include BED.
func includeWeighting(includeCum map[string]int) ([]string, []int, int) {
	keys := make([]string, 0, len(includeCum))
	for k, n := range includeCum {
		if n > 0 {
			keys = append(keys, k)
		}
	}
	sort.Strings(keys)
	cum := make([]int, len(keys))
	total := 0
	for i, k := range keys {
		total += includeCum[k]
		cum[i] = total
	}
	return keys, cum, total
}

// ParseGenome reads a tab-separated chrom-size file ("chr\tlen\n") and
// returns the chrom-to-size map. Comment lines starting with '#' are
// ignored.
func ParseGenome(r io.Reader) (map[string]int, error) {
	out := map[string]int{}
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("genome line %q: expected 'chrom\\tsize'", line)
		}
		n, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid size %q in genome line: %v", fields[1], err)
		}
		out[fields[0]] = n
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// ParseBED is a small helper that reads BED records from r using the shared
// bed.Reader and returns them as a slice. Useful for loading include /
// exclude BED files in the CLI.
func ParseBED(r io.Reader) ([]*bed.Record, error) {
	return bed.NewReader(r).ReadAll()
}
