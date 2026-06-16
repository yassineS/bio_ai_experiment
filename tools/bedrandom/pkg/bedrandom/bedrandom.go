// Package bedrandom is a Go re-implementation of `bedtools random` (the
// upstream randomBed tool). It generates a requested number of fixed-length
// intervals placed uniformly at random across a genome and writes them as
// BED6 records.
//
// The placement algorithm mirrors upstream randomBed.cpp byte-for-byte: it
// draws a random offset into the concatenated genome with rand_range, projects
// that offset back onto a chromosome (file order, cumulative start offsets),
// rejects any interval whose end runs past the chromosome length and redraws,
// then draws strand with a second rand_range(2). The random-number generator is
// the same std::mt19937_64 engine and rejection-sampling rand_range bound that
// upstream's default (non-USE_RAND) build uses, so a seeded run reproduces
// `bedtools random -seed N` exactly.
package bedrandom

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// DefaultLength is the default interval length (`-l`), matching upstream's
// default of 100.
const DefaultLength = 100

// DefaultNum is the default number of intervals to generate (`-n`), matching
// upstream's default of 1,000,000.
const DefaultNum = 1000000

// Options controls a single bedrandom run.
type Options struct {
	// Length is the length of each generated interval (upstream `-l`). Must be
	// positive.
	Length int

	// Num is the number of intervals to generate (upstream `-n`).
	Num int

	// Seed is the RNG seed (upstream `-seed`). Generate always seeds the engine
	// with this value; the CLI is responsible for deriving a time+pid seed when
	// the user did not supply one (see HaveSeed).
	Seed int

	// HaveSeed reports whether Seed was explicitly supplied by the user. The
	// core library does not branch on it (Generate always uses Seed); it is
	// carried here so the CLI can decide whether to substitute a time+pid seed,
	// matching upstream behaviour.
	HaveSeed bool
}

// Genome holds chromosome sizes in the order they appear in the genome file,
// together with the cumulative start offset of each chromosome within the
// concatenated genome. This ordering and the cumulative offsets are what make
// projectOnGenome byte-exact against upstream's GenomeFile.
type Genome struct {
	names   []string
	sizes   []int64
	offsets []int64 // offsets[i] is the genome offset where names[i] begins.
	total   int64
}

// Size returns the total concatenated genome size (sum of all chromosome
// sizes), matching upstream GenomeFile::getGenomeSize.
func (g *Genome) Size() int64 { return g.total }

// NumChroms returns the number of chromosomes in the genome.
func (g *Genome) NumChroms() int { return len(g.names) }

// projectOnGenome maps a 0-based genome offset back to a (chrom, localStart)
// pair, reproducing upstream GenomeFile::projectOnGenome exactly: it does a
// lower_bound (first start offset >= pos+1) over the start-offset vector, takes
// the index just before it, and subtracts that chromosome's start offset.
func (g *Genome) projectOnGenome(pos int64) (string, int64) {
	target := pos + 1
	i := sort.Search(len(g.offsets), func(i int) bool {
		return g.offsets[i] >= target
	})
	chrom := g.names[i-1]
	start := pos - g.offsets[i-1]
	return chrom, start
}

// chromSize returns the size of the named chromosome. The name always comes
// from projectOnGenome, so it is guaranteed present.
func (g *Genome) chromSize(name string) int64 {
	for i, n := range g.names {
		if n == name {
			return g.sizes[i]
		}
	}
	return -1
}

// ParseGenome reads a genome (chrom-sizes) file from r, preserving file order
// and building the cumulative start-offset table. The format is one
// `chrom<TAB>size` (or whitespace-separated) record per line; lines whose first
// field begins with `#` and blank lines are skipped, matching upstream
// GenomeFile parsing. A line with fewer than two fields is an error; a
// non-numeric size is silently skipped (matching upstream's strtol guard).
func ParseGenome(r io.Reader) (*Genome, error) {
	g := &Genome{}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	lineNum := 0
	for sc.Scan() {
		lineNum++
		line := sc.Text()
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		// Upstream tests genomeFields[0].find("#") != 0, i.e. it skips the
		// line only when '#' is the first character of the first field.
		if strings.HasPrefix(fields[0], "#") {
			continue
		}
		if len(fields) < 2 {
			return nil, fmt.Errorf("genome file line %d: expected at least 2 fields, got %d", lineNum, len(fields))
		}
		size, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			continue
		}
		g.names = append(g.names, fields[0])
		g.sizes = append(g.sizes, size)
		g.offsets = append(g.offsets, g.total)
		g.total += size
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if len(g.names) == 0 {
		return nil, errors.New("genome file contains no chromosomes")
	}
	return g, nil
}

// Generate produces opts.Num random intervals over genome g and writes them as
// BED6 records to w. It seeds a fresh std::mt19937_64 engine with opts.Seed and
// reproduces upstream's draw order exactly:
//
//	for each interval:
//	    do { off = rand_range(genomeSize); project; } while (end > chromSize)
//	    strand = rand_range(2) >= 1 ? '+' : '-'
//	    print chrom start end index length strand
//
// where index counts 1..Num and the score column holds the interval length.
// Generate returns the number of records written. Length must be positive.
func Generate(g *Genome, w io.Writer, opts Options) (int, error) {
	if opts.Length <= 0 {
		return 0, fmt.Errorf("interval length must be positive, got %d", opts.Length)
	}
	if opts.Num < 0 {
		return 0, fmt.Errorf("number of intervals must be non-negative, got %d", opts.Num)
	}
	if g.total <= 0 {
		return 0, errors.New("genome size is zero")
	}
	length := int64(opts.Length)

	// std::mt19937_64::seed takes the int seed; an int converts to the engine's
	// 64-bit result_type. Match that conversion (sign-extend negatives) here.
	mt := newMT19937_64(uint64(int64(opts.Seed)))
	bw := bufio.NewWriter(w)

	genomeSize := uint64(g.total)
	written := 0
	for n := 0; n < opts.Num; n++ {
		var chrom string
		var start, end int64
		for {
			off := int64(mt.randRange(genomeSize))
			chrom, start = g.projectOnGenome(off)
			end = start + length
			if end <= g.chromSize(chrom) {
				break
			}
		}
		index := n + 1
		strand := byte('-')
		if mt.randRange(2) >= 1 {
			strand = '+'
		}
		// %s\t%d\t%d\t%d\t%d\t%c\n
		bw.WriteString(chrom)
		bw.WriteByte('\t')
		bw.WriteString(strconv.FormatInt(start, 10))
		bw.WriteByte('\t')
		bw.WriteString(strconv.FormatInt(end, 10))
		bw.WriteByte('\t')
		bw.WriteString(strconv.Itoa(index))
		bw.WriteByte('\t')
		bw.WriteString(strconv.FormatInt(end-start, 10))
		bw.WriteByte('\t')
		bw.WriteByte(strand)
		bw.WriteByte('\n')
		written++
	}
	if err := bw.Flush(); err != nil {
		return written, err
	}
	return written, nil
}
