// Package bedgenomecov computes per-base or summarised coverage of BED/bedGraph
// intervals against a chromosome-size genome file, mirroring the behaviour of
// `bedtools genomecov`.
//
// Three output modes are supported:
//
//   - Histogram (default): for every (chrom, depth) pair, report the number of
//     bases at that depth, the chromosome size and the fraction of bases at
//     that depth. A `genome` row aggregates across all chromosomes.
//   - bedGraph: report contiguous runs of constant depth as
//     `chrom\tstart\tend\tdepth`. Zero-depth runs are emitted only when
//     ReportZeroBedGraph is true (the `-bga` flag).
//   - Per-base: one line per base, `chrom\tposition\tdepth` (1-based). When
//     SkipZero is true (the `-dz` flag) zero-depth bases are skipped.
//
// Coverage is computed per chromosome with a single int slice sized to the
// chromosome length declared in the genome file. This is simple, fast and
// matches bedtools' typical memory footprint for vertebrate-scale inputs.
package bedgenomecov

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Mode selects the output style for Run.
type Mode int

const (
	// ModeHistogram is the default `bedtools genomecov` output.
	ModeHistogram Mode = iota
	// ModeBedGraph emits non-zero runs of constant depth (the `-bg` flag).
	ModeBedGraph
	// ModeBedGraphAll emits every run of constant depth, including zero (the `-bga` flag).
	ModeBedGraphAll
	// ModePerBase emits one record per base position (the `-d` flag).
	ModePerBase
	// ModePerBaseNonZero emits one record per non-zero base position (the `-dz` flag).
	ModePerBaseNonZero
)

// Options configures the coverage computation.
type Options struct {
	Mode       Mode
	Strand     string  // "+", "-" or "" for both
	MaxDepth   int     // cap depth in histograms (0 = uncapped)
	Scale      float64 // multiplier applied to every depth (default 1.0)
	FivePrime  bool    // count only the 5'-most base of each interval
	ThreePrime bool    // count only the 3'-most base of each interval
	TrackLine  bool    // emit a UCSC trackline header in bedGraph modes
	TrackOpts  string  // optional `key=value` string appended after `track`
}

// GenomeSize holds chromosome ordering and lengths parsed from a genome file.
type GenomeSize struct {
	Order  []string       // chromosomes in the order they appeared
	Length map[string]int // chrom -> length
}

// ReadGenome parses a chromosome-sizes file (`chrom\tsize` per line). Lines
// starting with `#` or empty lines are skipped. Returns an error if a size is
// missing or non-positive.
func ReadGenome(r io.Reader) (*GenomeSize, error) {
	g := &GenomeSize{Length: map[string]int{}}
	scanner := bufio.NewScanner(r)
	// Allow long genome lines.
	buf := make([]byte, 0, 1<<16)
	scanner.Buffer(buf, 1<<24)
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		raw := strings.TrimRight(scanner.Text(), "\r")
		line := strings.TrimSpace(raw)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("genome file line %d: need chrom and size, got %q", lineNo, raw)
		}
		size, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("genome file line %d: invalid size %q: %w", lineNo, fields[1], err)
		}
		if size <= 0 {
			return nil, fmt.Errorf("genome file line %d: size must be > 0, got %d", lineNo, size)
		}
		if _, dup := g.Length[fields[0]]; dup {
			return nil, fmt.Errorf("genome file line %d: chromosome %q listed twice", lineNo, fields[0])
		}
		g.Order = append(g.Order, fields[0])
		g.Length[fields[0]] = size
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading genome file: %w", err)
	}
	if len(g.Order) == 0 {
		return nil, errors.New("genome file is empty")
	}
	return g, nil
}

// Run reads BED intervals from input, computes coverage against the genome, and
// writes the configured output to writer. It is the only public entry point.
func Run(input io.Reader, genome *GenomeSize, writer io.Writer, opts Options) error {
	if genome == nil || len(genome.Order) == 0 {
		return errors.New("genome size info is required")
	}
	if opts.FivePrime && opts.ThreePrime {
		return errors.New("cannot combine 5' and 3' end-only counting")
	}
	if opts.Strand != "" && opts.Strand != "+" && opts.Strand != "-" {
		return fmt.Errorf("strand must be \"+\" or \"-\", got %q", opts.Strand)
	}
	if opts.Scale == 0 {
		opts.Scale = 1.0
	}

	// Allocate per-chromosome depth arrays in genome-declared order.
	depth := make(map[string][]int, len(genome.Order))
	for _, chrom := range genome.Order {
		depth[chrom] = make([]int, genome.Length[chrom])
	}

	bw := bufio.NewWriter(writer)
	defer bw.Flush()

	// Ingest BED records and bump counts.
	br := bed.NewReader(input)
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("reading input: %w", err)
		}
		arr, ok := depth[rec.Chrom]
		if !ok {
			// Skip silently: chromosome not in genome file (matches bedtools).
			continue
		}
		if opts.Strand != "" && rec.Strand != opts.Strand {
			continue
		}
		start, end := clamp(rec.ChromStart, rec.ChromEnd, len(arr))
		if start >= end {
			continue
		}
		switch {
		case opts.FivePrime:
			pos := start
			if rec.Strand == "-" {
				pos = end - 1
			}
			if pos >= 0 && pos < len(arr) {
				arr[pos]++
			}
		case opts.ThreePrime:
			pos := end - 1
			if rec.Strand == "-" {
				pos = start
			}
			if pos >= 0 && pos < len(arr) {
				arr[pos]++
			}
		default:
			for i := start; i < end; i++ {
				arr[i]++
			}
		}
	}

	// Trackline for bedGraph modes.
	if opts.TrackLine && (opts.Mode == ModeBedGraph || opts.Mode == ModeBedGraphAll) {
		line := "track type=bedGraph"
		if opts.TrackOpts != "" {
			line += " " + opts.TrackOpts
		}
		if _, err := fmt.Fprintln(bw, line); err != nil {
			return err
		}
	}

	switch opts.Mode {
	case ModeBedGraph, ModeBedGraphAll:
		return writeBedGraph(bw, depth, genome, opts)
	case ModePerBase, ModePerBaseNonZero:
		return writePerBase(bw, depth, genome, opts)
	default:
		return writeHistogram(bw, depth, genome, opts)
	}
}

// clamp limits a BED interval to [0, length).
func clamp(start, end, length int) (int, int) {
	if start < 0 {
		start = 0
	}
	if end > length {
		end = length
	}
	return start, end
}

// scaledDepth applies the scale factor.
func scaledDepth(d int, scale float64) float64 {
	return float64(d) * scale
}

// histogramKey converts a (possibly scaled) depth to the integer bucket used
// in histogram output. The cap is applied after scaling so that `-max` always
// limits the reported depth.
func histogramKey(d int, opts Options) int {
	v := int(scaledDepth(d, opts.Scale))
	if v < 0 {
		v = 0
	}
	if opts.MaxDepth > 0 && v > opts.MaxDepth {
		v = opts.MaxDepth
	}
	return v
}

// formatDepth renders a (possibly scaled) depth, emitting an integer when the
// scale is 1.0 (the common case) and a compact float otherwise. Matches
// bedtools' `%g`-style formatting under `-scale`.
func formatDepth(d float64, scale float64) string {
	if scale == 1.0 {
		return strconv.Itoa(int(d))
	}
	return strconv.FormatFloat(d, 'g', 10, 64)
}

func writeBedGraph(w *bufio.Writer, depth map[string][]int, g *GenomeSize, opts Options) error {
	includeZero := opts.Mode == ModeBedGraphAll
	for _, chrom := range g.Order {
		arr := depth[chrom]
		i := 0
		for i < len(arr) {
			d := scaledDepth(arr[i], opts.Scale)
			j := i + 1
			for j < len(arr) && scaledDepth(arr[j], opts.Scale) == d {
				j++
			}
			if d != 0 || includeZero {
				if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\n", chrom, i, j, formatDepth(d, opts.Scale)); err != nil {
					return err
				}
			}
			i = j
		}
	}
	return nil
}

func writePerBase(w *bufio.Writer, depth map[string][]int, g *GenomeSize, opts Options) error {
	skipZero := opts.Mode == ModePerBaseNonZero
	for _, chrom := range g.Order {
		arr := depth[chrom]
		for i, raw := range arr {
			d := scaledDepth(raw, opts.Scale)
			if skipZero && d == 0 {
				continue
			}
			if _, err := fmt.Fprintf(w, "%s\t%d\t%s\n", chrom, i+1, formatDepth(d, opts.Scale)); err != nil {
				return err
			}
		}
	}
	return nil
}

// writeHistogram emits per-chromosome histograms followed by a `genome` row,
// matching the exact column layout of `bedtools genomecov`.
func writeHistogram(w *bufio.Writer, depth map[string][]int, g *GenomeSize, opts Options) error {
	genome := map[int]int{}
	genomeSize := 0
	for _, chrom := range g.Order {
		arr := depth[chrom]
		counts := map[int]int{}
		for _, raw := range arr {
			d := histogramKey(raw, opts)
			counts[d]++
			genome[d]++
		}
		chromSize := len(arr)
		genomeSize += chromSize
		for _, d := range sortedKeys(counts) {
			n := counts[d]
			frac := 0.0
			if chromSize > 0 {
				frac = float64(n) / float64(chromSize)
			}
			if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%s\n", chrom, d, n, chromSize, formatFraction(frac)); err != nil {
				return err
			}
		}
	}
	for _, d := range sortedKeys(genome) {
		n := genome[d]
		frac := 0.0
		if genomeSize > 0 {
			frac = float64(n) / float64(genomeSize)
		}
		if _, err := fmt.Fprintf(w, "genome\t%d\t%d\t%d\t%s\n", d, n, genomeSize, formatFraction(frac)); err != nil {
			return err
		}
	}
	return nil
}

// sortedKeys returns the keys of m in ascending order.
func sortedKeys(m map[int]int) []int {
	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	return keys
}

// formatFraction renders a fraction with up to 10 significant digits and
// without trailing-zero noise, matching the bedtools output style.
func formatFraction(f float64) string {
	s := strconv.FormatFloat(f, 'g', 10, 64)
	return s
}
