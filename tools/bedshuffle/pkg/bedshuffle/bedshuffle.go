// Package bedshuffle randomly relocates BED intervals across a genome,
// mirroring `bedtools shuffle`.
//
// For each input interval the package draws a new (chrom, start) coordinate
// and emits a record with the same length, keeping name / score / strand /
// extras intact. Various flags constrain the draw:
//
//   - Include: restrict placements to regions listed in an "include" BED,
//     weighted by interval size (upstream's -incl).
//   - Exclude: forbid placements that overlap any region listed in an
//     "exclude" BED (upstream's -excl).
//   - Chrom: keep each shuffled interval on its original chromosome
//     (upstream's -chrom; implies "choose chrom, then position").
//   - ChromFirst: pick a chromosome uniformly first, then a position within it
//     (upstream's -chromFirst).
//
// # Byte-for-byte parity with upstream
//
// Placement uses a faithful Go port of std::mt19937_64 (mt19937.go) — the exact
// 64-bit Mersenne Twister that upstream bedtools' Random.cpp uses (the default,
// non-USE_RAND build, which the shipped bedtools binary is). The draw order,
// the rand_range rejection bound, the genome-projection layout (chromosomes in
// genome-file order), and the per-mode retry loops all mirror
// reference_code/bedtools/src/shuffleBed/shuffleBed.cpp exactly, so
// `bedshuffle -seed N` reproduces `bedtools shuffle -seed N` byte-for-byte.
package bedshuffle

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Options controls Shuffle.
type Options struct {
	// Genome is a tab-delimited chromosome-size map; key is the chrom name,
	// value is the chromosome length in bp.
	Genome map[string]int

	// GenomeOrder lists the chromosome names in genome-file order. Upstream
	// bedtools builds its start-offset table and "project a genome-wide
	// position onto a chromosome" logic in file order, so byte-for-byte parity
	// requires that exact order. When nil, Shuffle falls back to the map keys
	// sorted lexicographically (deterministic, but not parity-exact for the
	// default genome-projection mode). ParseGenomeOrdered populates this.
	GenomeOrder []string

	// Include limits placements to regions listed here, weighted by interval
	// size (upstream -incl). nil means "place anywhere on the chromosome".
	Include []*bed.Record
	// Exclude rejects placements that overlap any region here (upstream -excl).
	Exclude []*bed.Record

	// Chrom (`-chrom`): keep each shuffled interval on its original
	// chromosome. Upstream sets both chooseChrom and sameChrom.
	Chrom bool

	// ChromFirst (`-chromFirst`): choose the destination chromosome uniformly
	// at random first (each chromosome weighted equally), then a position
	// within it — instead of the default of projecting a genome-wide position
	// (weighting a chromosome by its size). Note: upstream's -incl path ignores
	// -chromFirst entirely (the include selection is always size-weighted), so
	// ChromFirst only affects the no-include case.
	ChromFirst bool

	// AllowBeyondChromEnd (`-allowBeyondChromEnd`): permit a shuffled interval
	// to be clamped to the chromosome end instead of being redrawn when it
	// would exceed the end. Upstream's preventExceedingChromEnd is the negation
	// of this; the default (false) redraws until the interval fits.
	AllowBeyondChromEnd bool

	// Seed for the RNG. Same seed + same inputs ⇒ same output. Mirrors
	// upstream `-seed`.
	Seed int64

	// MaxRetries caps the number of placement attempts per interval for the
	// exclude / include retry loops. 0 ⇒ default of 1000 (mirrors upstream
	// -maxTries).
	MaxRetries int
}

// orderedGenome holds the genome in upstream's file order with the cumulative
// start-offset table used by projectOnGenome.
type orderedGenome struct {
	names        []string
	startOffsets []uint64
	sizes        map[string]uint64
	genomeSize   uint64
}

// buildOrderedGenome constructs the file-order start-offset table that mirrors
// GenomeFile::loadGenomeFileIntoMap. When order is empty it falls back to the
// lexicographically sorted map keys.
func buildOrderedGenome(g map[string]int, order []string) *orderedGenome {
	names := order
	if len(names) == 0 {
		names = make([]string, 0, len(g))
		for k := range g {
			names = append(names, k)
		}
		sort.Strings(names)
	}
	og := &orderedGenome{
		names:        names,
		startOffsets: make([]uint64, 0, len(names)),
		sizes:        make(map[string]uint64, len(names)),
	}
	for _, name := range names {
		size := uint64(g[name])
		og.startOffsets = append(og.startOffsets, og.genomeSize)
		og.sizes[name] = size
		og.genomeSize += size
	}
	return og
}

// projectOnGenome maps a genome-wide position to a (chrom, start) pair,
// matching GenomeFile::projectOnGenome: lower_bound(startOffsets, pos+1) then
// chrom = names[i-1], start = pos - startOffsets[i-1].
func (og *orderedGenome) projectOnGenome(pos uint64) (string, uint64) {
	// lower_bound: first index whose offset >= pos+1.
	i := sort.Search(len(og.startOffsets), func(k int) bool {
		return og.startOffsets[k] >= pos+1
	})
	chrom := og.names[i-1]
	start := pos - og.startOffsets[i-1]
	return chrom, start
}

// Shuffle reads BED records from r, draws new coordinates for each according
// to opts, and writes the relocated records to w. Returns the number of
// records placed; the first interval that cannot be placed within MaxRetries
// triggers an error.
func Shuffle(r io.Reader, w io.Writer, opts Options) (int, error) {
	if len(opts.Genome) == 0 {
		return 0, fmt.Errorf("shuffle: -g (genome file) is required")
	}
	maxRetries := opts.MaxRetries
	if maxRetries == 0 {
		maxRetries = 1000
	}

	mt := newMT19937_64(uint64(opts.Seed))
	genome := buildOrderedGenome(opts.Genome, opts.GenomeOrder)

	// Build the size-weighted include table exactly as
	// BedFile::assignWeightsBasedOnSize does: sort the include intervals by
	// size ascending, then assign each a cumulative weight equal to the
	// running fraction of total include bp.
	var inclWeighted []weightedInterval
	if len(opts.Include) > 0 {
		inclWeighted = buildIncludeWeights(opts.Include)
	}

	// Index exclude regions by chrom for the overlap check, sorted by start.
	exclude := indexByChrom(opts.Exclude)
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
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
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

		newChrom, newStart, newEnd, ok := chooseLocus(mt, genome, inclWeighted,
			exclude, opts, fields[0], uint64(length), maxRetries)
		if !ok {
			return count, fmt.Errorf("Error, line %d: tried %d potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on.",
				count+1, maxRetries)
		}

		out := append([]string(nil), fields...)
		out[0] = newChrom
		out[1] = strconv.FormatUint(newStart, 10)
		out[2] = strconv.FormatUint(newEnd, 10)
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

// chooseLocus dispatches to the include / exclude / plain placement loop that
// matches the corresponding upstream BedShuffle method, returning the placed
// (chrom, start, end). ok is false only when an exclude/include retry loop is
// exhausted (mirrors upstream's "tried N loci" error path).
func chooseLocus(
	mt *mt19937_64,
	genome *orderedGenome,
	inclWeighted []weightedInterval,
	exclude map[string][]*bed.Record,
	opts Options,
	origChrom string,
	length uint64,
	maxRetries int,
) (string, uint64, uint64, bool) {
	hasInclude := len(inclWeighted) > 0
	hasExclude := len(exclude) > 0

	switch {
	case !hasInclude && !hasExclude:
		// BedShuffle::Shuffle → ChooseLocus.
		return chooseLocusPlain(mt, genome, opts, origChrom, length, maxRetries)

	case hasExclude && !hasInclude:
		// BedShuffle::ShuffleWithExclusions: redraw (via ChooseLocus) while the
		// placement overlaps an exclude region, up to maxTries.
		for tries := 0; tries <= maxRetries; tries++ {
			c, s, e, ok := chooseLocusPlain(mt, genome, opts, origChrom, length, maxRetries)
			if ok && !hasOverlap(exclude[c], int(s), int(e)) {
				return c, s, e, true
			}
		}
		return "", 0, 0, false

	case hasInclude && !hasExclude:
		// BedShuffle::ShuffleWithInclusions: redraw (via
		// ChooseLocusFromInclusionFile) while end > chromSize, up to maxTries.
		for tries := 0; tries <= maxRetries; tries++ {
			c, s, e := chooseLocusFromInclude(mt, inclWeighted, length)
			if e <= genome.sizes[c] {
				return c, s, e, true
			}
		}
		return "", 0, 0, false

	default:
		// BedShuffle::ShuffleWithInclusionsAndExclusions: redraw (via
		// ChooseLocusFromInclusionFile) while the placement overlaps an exclude
		// region, up to maxTries.
		for tries := 0; tries <= maxRetries; tries++ {
			c, s, e := chooseLocusFromInclude(mt, inclWeighted, length)
			if !hasOverlap(exclude[c], int(s), int(e)) {
				return c, s, e, true
			}
		}
		return "", 0, 0, false
	}
}

// chooseLocusPlain mirrors BedShuffle::ChooseLocus for the non-include modes.
// It honours Chrom (sameChrom), ChromFirst (chooseChrom without sameChrom), and
// the default genome-projection mode, plus the preventExceedingChromEnd retry.
func chooseLocusPlain(
	mt *mt19937_64,
	genome *orderedGenome,
	opts Options,
	origChrom string,
	length uint64,
	maxRetries int,
) (string, uint64, uint64, bool) {
	prevent := !opts.AllowBeyondChromEnd

	if !opts.Chrom && !opts.ChromFirst {
		// Default: project a uniform genome-wide position onto a chromosome.
		// Mirrors ChooseLocus's bounded do-while: retry up to _maxTries while
		// the interval exceeds the chrom end. Upstream warns and "ignores the
		// entry" if it never fits; this port surfaces that as ok=false so the
		// caller skips the line (fix-on-port: upstream actually drops it too).
		for tries := 0; tries <= maxRetries; tries++ {
			pos := mt.randRange(genome.genomeSize)
			chrom, start := genome.projectOnGenome(pos)
			end := start + length
			chromSize := genome.sizes[chrom]
			if end > chromSize && !prevent {
				return chrom, start, chromSize, true
			}
			if end <= chromSize {
				return chrom, start, end, true
			}
		}
		return "", 0, 0, false
	}

	// chooseChrom modes: -chrom (sameChrom) keeps origChrom; -chromFirst picks
	// a chromosome uniformly then a position within it. Upstream's else branch
	// loops `while (end > chromSize)` with no try bound; we add the same
	// maxTries guard the default path uses so a too-large interval cannot hang.
	if opts.Chrom {
		if _, ok := genome.sizes[origChrom]; !ok {
			// Upstream's getChromSize returns -1 here, producing garbage
			// coordinates; we treat the missing chromosome as unplaceable.
			return "", 0, 0, false
		}
	}
	for tries := 0; tries <= maxRetries; tries++ {
		var chrom string
		var chromSize uint64
		if !opts.Chrom {
			// -chromFirst: pick a chromosome uniformly by index.
			idx := mt.randRange(uint64(len(genome.names)))
			chrom = genome.names[idx]
			chromSize = genome.sizes[chrom]
		} else {
			chrom = origChrom
			chromSize = genome.sizes[chrom]
		}
		start := mt.randRange(chromSize)
		end := start + length
		if end > chromSize && !prevent {
			return chrom, start, chromSize, true
		}
		if end <= chromSize {
			return chrom, start, end, true
		}
	}
	return "", 0, 0, false
}

// weightedInterval is an include interval with its cumulative size weight,
// matching BED.weight after assignWeightsBasedOnSize.
type weightedInterval struct {
	chrom  string
	start  uint64
	size   uint64
	weight float64
}

// buildIncludeWeights replicates BedFile::assignWeightsBasedOnSize: sort the
// include intervals by size ascending (stable, matching std::sort's behaviour
// for the parity fixtures), then assign each a cumulative weight equal to the
// running fraction of total include bp.
func buildIncludeWeights(include []*bed.Record) []weightedInterval {
	ivs := make([]weightedInterval, len(include))
	for i, r := range include {
		size := r.ChromEnd - r.ChromStart
		if size < 0 {
			size = 0
		}
		ivs[i] = weightedInterval{chrom: r.Chrom, start: uint64(r.ChromStart), size: uint64(size)}
	}
	sort.SliceStable(ivs, func(i, j int) bool { return ivs[i].size < ivs[j].size })

	var totalSize uint64
	for i := range ivs {
		totalSize += ivs[i].size
	}
	var totalWeight float64
	for i := range ivs {
		if totalSize > 0 {
			totalWeight += float64(ivs[i].size) / float64(totalSize)
		}
		ivs[i].weight = totalWeight
	}
	return ivs
}

// chooseLocusFromInclude mirrors BedShuffle::ChooseLocusFromInclusionFile:
// draw a proportion, size-weighted-search for the include interval, then a
// uniform start within it.
func chooseLocusFromInclude(
	mt *mt19937_64,
	inclWeighted []weightedInterval,
	length uint64,
) (string, uint64, uint64) {
	runif := mt.randProportion()
	iv := sizeWeightedSearch(inclWeighted, runif)
	start := iv.start + mt.randRange(iv.size)
	return iv.chrom, start, start + length
}

// sizeWeightedSearch mirrors BedFile::sizeWeightedSearch: upper_bound with
// CompareByWeight, i.e. the first interval whose cumulative weight is strictly
// greater than val.
func sizeWeightedSearch(ivs []weightedInterval, val float64) weightedInterval {
	idx := sort.Search(len(ivs), func(i int) bool { return ivs[i].weight > val })
	if idx >= len(ivs) {
		idx = len(ivs) - 1
	}
	return ivs[idx]
}

// indexByChrom buckets records by chromosome.
func indexByChrom(recs []*bed.Record) map[string][]*bed.Record {
	out := map[string][]*bed.Record{}
	for _, r := range recs {
		out[r.Chrom] = append(out[r.Chrom], r)
	}
	return out
}

// hasOverlap reports whether [start, end) overlaps any of recs. This mirrors
// upstream's default overlap fraction of 1E-9 (any ≥1bp overlap).
func hasOverlap(recs []*bed.Record, start, end int) bool {
	for _, r := range recs {
		if r.ChromStart < end && r.ChromEnd > start {
			return true
		}
	}
	return false
}

// ParseGenome reads a tab-separated chrom-size file ("chrom\tsize\n") and
// returns the chrom-to-size map. Comment lines starting with '#' are ignored.
// For byte-for-byte upstream parity use ParseGenomeOrdered, which also returns
// the chromosome file order.
func ParseGenome(r io.Reader) (map[string]int, error) {
	g, _, err := ParseGenomeOrdered(r)
	return g, err
}

// ParseGenomeOrdered reads a tab-separated chrom-size file and returns both the
// chrom-to-size map and the chromosome names in file order. The order is what
// upstream bedtools uses to build its genome-projection table, so it is
// required for byte-for-byte parity in the default placement mode.
func ParseGenomeOrdered(r io.Reader) (map[string]int, []string, error) {
	out := map[string]int{}
	var order []string
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, nil, fmt.Errorf("genome line %q: expected 'chrom\\tsize'", line)
		}
		n, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, nil, fmt.Errorf("invalid size %q in genome line: %v", fields[1], err)
		}
		if _, seen := out[fields[0]]; !seen {
			order = append(order, fields[0])
		}
		out[fields[0]] = n
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, err
	}
	return out, order, nil
}

// ParseBED is a small helper that reads BED records from r using the shared
// bed.Reader and returns them as a slice. Useful for loading include /
// exclude BED files in the CLI.
func ParseBED(r io.Reader) ([]*bed.Record, error) {
	return bed.NewReader(r).ReadAll()
}
