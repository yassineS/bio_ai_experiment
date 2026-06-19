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

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
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
	// Split ("-split") treats a BED12 record as its blocks: coverage is
	// counted over each block interval rather than the whole record span,
	// matching upstream bedtools genomecov -split.
	Split bool
	// PairedCoverage ("-pc", BAM only) counts coverage of paired-end
	// fragments: each properly-paired read with a positive insert size (TLEN)
	// contributes the fragment span [pos, pos+TLEN); the negative-TLEN mate is
	// skipped so each fragment is counted once.
	PairedCoverage bool
	// FragmentSize ("-fs N", BAM only), when > 0, forces every read to cover a
	// fixed N-base fragment anchored at its 5' end instead of its alignment
	// length: forward reads cover [pos, pos+N), reverse reads cover [end-N, end).
	FragmentSize int
	// CRAMReference names a FASTA file used to reconstruct reference-backed
	// reads when the alignment input (-ibam) is a CRAM. It is ignored for SAM
	// and BAM input, which carry their sequence inline. When empty a
	// reference-backed CRAM still decodes (reference-derived bases fill with
	// 'N'), and the REF_CACHE/REF_PATH environment variables are honoured.
	CRAMReference string
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

// recordSource is the minimal interface the coverage accumulator consumes;
// both *bed.Reader (BED input) and *alnbed.Reader (BAM/SAM input) satisfy it.
type recordSource interface {
	Read() (*bed.Record, error)
}

// Run reads BED intervals from input, computes coverage against the genome, and
// writes the configured output to writer.
func Run(input io.Reader, genome *GenomeSize, writer io.Writer, opts Options) error {
	if genome == nil || len(genome.Order) == 0 {
		return errors.New("genome size info is required")
	}
	// The fast BED source parses directly from the scanner buffer, avoiding the
	// shared reader's per-line string/field-slice allocations. Block columns are
	// only consumed under -split, so they are parsed only then.
	return runCore(newFastBEDSource(input, opts.Split), genome, writer, opts)
}

// runCore is the BED-input coverage pipeline: per-record coverage intervals
// honour opts.Split (a BED12 record is split into its blocks only under
// -split).
func runCore(src recordSource, genome *GenomeSize, writer io.Writer, opts Options) error {
	return runCoreBlocks(src, genome, writer, opts, false)
}

// RunBAM reads a SAM/BAM/CRAM alignment stream from input, derives the genome
// from the alignment header's @SQ entries (mirroring upstream `genomecov
// -ibam`, where no separate genome file is supplied), and writes coverage to
// writer. Each mapped alignment contributes its reference blocks, split on a
// deletion (D) op and — under opts.Split — also on a skip (N) op, matching
// upstream getBamBlocks. The input format is auto-detected; for a
// reference-backed CRAM, opts.CRAMReference names the decode FASTA (the
// REF_CACHE/REF_PATH environment variables are also honoured).
func RunBAM(input io.Reader, writer io.Writer, opts Options) error {
	if opts.PairedCoverage && opts.FragmentSize > 0 {
		return errors.New("cannot combine -pc and -fs")
	}
	sr, err := alnio.NewReaderWithReference(input, opts.CRAMReference)
	if err != nil {
		return fmt.Errorf("reading alignment input: %w", err)
	}
	hdr := sr.Header()
	genome := &GenomeSize{Length: map[string]int{}}
	if hdr != nil {
		for _, ref := range hdr.Refs {
			genome.Order = append(genome.Order, ref.Name)
			genome.Length[ref.Name] = int(ref.Length)
		}
	}
	if len(genome.Order) == 0 {
		return errors.New("alignment header has no @SQ reference entries")
	}
	return runCoreBlocks(&alnFragmentSource{sr: sr, opts: opts}, genome, writer, opts, true)
}

// alnFragmentSource adapts a SAM/BAM reader to the recordSource contract,
// converting each mapped alignment into the *bed.Record genomecov should
// count. The shape depends on opts: by default the BED12 block record
// (alnbed.ToBED12, so -split works); under -pc the paired-fragment span; and
// under -fs the strand-anchored fixed-size span.
type alnFragmentSource struct {
	sr   sam.Reader
	opts Options
}

// Read returns the next alignment-derived record, skipping unmapped reads and
// (under -pc) reads that do not start a counted fragment.
func (s *alnFragmentSource) Read() (*bed.Record, error) {
	for {
		rec, err := s.sr.Read()
		if err != nil {
			return nil, err
		}
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" || rec.Pos <= 0 {
			continue
		}
		switch {
		case s.opts.PairedCoverage:
			// Count each proper pair once, from the read carrying the positive
			// insert size: fragment = [pos, pos+TLEN).
			if !rec.IsProperPair() || rec.TLen <= 0 {
				continue
			}
			start := int(rec.Pos) - 1
			return spanRecord(rec, start, start+int(rec.TLen)), nil
		case s.opts.FragmentSize > 0:
			start := int(rec.Pos) - 1
			end := int(rec.EndPosition())
			if rec.Flag&sam.FlagReverse != 0 {
				return spanRecord(rec, end-s.opts.FragmentSize, end), nil
			}
			return spanRecord(rec, start, start+s.opts.FragmentSize), nil
		default:
			// genomecov's BAM path computes coverage blocks itself rather than
			// reusing the generic alnbed.ToBED12 BED12 conversion. Upstream calls
			// GetBamBlocks(..., breakOnDeletionOps=!_ignoreD, obeySplits=_split):
			// it always breaks the alignment on a D (deletion) op — _ignoreD is
			// not exposed by genomecov's CLI, so it is always false — and breaks
			// on an N (skip) op only under -split. The shared alnbed.ToBED12 only
			// ever breaks on N (the common BAM-to-BED convention used by the
			// other bed* tools), so it would merge a deletion into one span; the
			// genomecov-specific helper below restores the upstream behaviour.
			return genomecovAlnRecord(rec, s.opts.Split), nil
		}
	}
}

// genomecovAlnRecord builds the *bed.Record genomecov counts for one mapped
// alignment, with reference-consuming CIGAR blocks split per the upstream
// GetBamBlocks rules: always break on a D (deletion) op, and additionally
// break on an N (skip) op when split is true. The resulting blocks are stored
// as BED12 block columns so the coverage accumulator (which always iterates
// the per-record blocks for alignment input) counts each block once.
func genomecovAlnRecord(rec *sam.Record, split bool) *bed.Record {
	start := int(rec.Pos) - 1
	blocks := genomecovBlocks(rec, start, split)
	end := start
	if len(blocks) > 0 {
		end = blocks[len(blocks)-1][1]
	}
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	sizes := make([]int, len(blocks))
	starts := make([]int, len(blocks))
	for i, b := range blocks {
		sizes[i] = b[1] - b[0]
		starts[i] = b[0] - start
	}
	return &bed.Record{
		Chrom:       rec.RName,
		ChromStart:  start,
		ChromEnd:    end,
		Name:        rec.QName,
		Score:       int(rec.MapQ),
		Strand:      strand,
		BlockCount:  len(blocks),
		BlockSizes:  sizes,
		BlockStarts: starts,
	}
}

// genomecovBlocks returns the 0-based half-open reference blocks of an
// alignment beginning at refStart, mirroring upstream getBamBlocks. M/=/X
// consume the reference and extend the current block; a D op always closes the
// current block and starts a new one after the deletion (breakOnDeletionOps is
// always true for genomecov); an N op closes the current block and, when
// breakOnN (-split) is set, starts a new one after the gap — without -split the
// N still advances the reference but the surrounding blocks are not separated
// here (upstream's obeySplits=false path never reaches getBamBlocks with N
// breaks, but N is exceedingly rare in non-split data). I/S/H/P do not consume
// reference.
func genomecovBlocks(rec *sam.Record, refStart int, breakOnN bool) [][2]int {
	var blocks [][2]int
	pos := refStart
	blockStart := refStart
	open := false
	for _, op := range rec.Cigar {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if !open {
				blockStart = pos
				open = true
			}
			pos += int(op.Length())
		case sam.CigarDeletion:
			// Always break on a deletion: close the current block (if any) and
			// resume after the deleted reference bases.
			if open {
				blocks = append(blocks, [2]int{blockStart, pos})
				open = false
			}
			pos += int(op.Length())
			blockStart = pos
		case sam.CigarSkipped:
			if breakOnN {
				if open {
					blocks = append(blocks, [2]int{blockStart, pos})
					open = false
				}
				pos += int(op.Length())
				blockStart = pos
			} else {
				// Without -split, an N op extends the span across the gap.
				if !open {
					blockStart = pos
					open = true
				}
				pos += int(op.Length())
			}
		default:
			// I, S, H, P: no reference advance.
		}
	}
	if open {
		blocks = append(blocks, [2]int{blockStart, pos})
	}
	if len(blocks) == 0 {
		blocks = append(blocks, [2]int{refStart, pos})
	}
	return blocks
}

// spanRecord builds a single-span (block-free) BED record on rec's chromosome
// and strand covering [start, end). Used for the -pc and -fs fragment spans.
func spanRecord(rec *sam.Record, start, end int) *bed.Record {
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	return &bed.Record{Chrom: rec.RName, ChromStart: start, ChromEnd: end, Strand: strand}
}

// runCoreBlocks is the shared coverage pipeline used by both the BED (Run) and
// the alignment (RunBAM) entry points. When alwaysBlocks is true the per-record
// coverage intervals are taken from the record's BED12 blocks regardless of
// opts.Split: the alignment source has already encoded genomecov's CIGAR
// D/N-split semantics into those blocks, so the accumulator must always honour
// them.
func runCoreBlocks(src recordSource, genome *GenomeSize, writer io.Writer, opts Options, alwaysBlocks bool) error {
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

	// Ingest records and bump counts.
	br := src
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
			// End-only counting works on the whole-record extent, not the
			// per-block split (upstream applies -5/-3 to the read, not blocks).
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
			// Under -split a BED12 record contributes coverage over each of
			// its blocks; otherwise the whole [start,end) span. For alignment
			// input (alwaysBlocks) the blocks already encode the genomecov
			// D/N-split semantics, so they are always honoured.
			for _, iv := range recordCoverIntervals(rec, opts.Split || alwaysBlocks) {
				s, e := clamp(iv[0], iv[1], len(arr))
				for i := s; i < e; i++ {
					arr[i]++
				}
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

// recordCoverIntervals returns the 0-based half-open intervals a record
// contributes coverage over. With split=true and a BED12 record carrying
// blocks, that is one interval per block (block start = ChromStart +
// BlockStarts[i], length BlockSizes[i]); otherwise the single [ChromStart,
// ChromEnd) span.
func recordCoverIntervals(rec *bed.Record, split bool) [][2]int {
	if split && rec.BlockCount > 0 && len(rec.BlockSizes) > 0 {
		out := make([][2]int, 0, len(rec.BlockSizes))
		for i := range rec.BlockSizes {
			s := rec.ChromStart
			if i < len(rec.BlockStarts) {
				s += rec.BlockStarts[i]
			}
			out = append(out, [2]int{s, s + rec.BlockSizes[i]})
		}
		return out
	}
	return [][2]int{{rec.ChromStart, rec.ChromEnd}}
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

func writeBedGraph(w *bufio.Writer, depth map[string][]int, g *GenomeSize, opts Options) error {
	includeZero := opts.Mode == ModeBedGraphAll
	intScale := opts.Scale == 1.0
	var buf []byte
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
				// chrom\tstart\tend\tdepth\n built into a reusable scratch buffer,
				// then written once. Avoids fmt.Fprintf's per-call interface
				// boxing and string allocations.
				buf = buf[:0]
				buf = append(buf, chrom...)
				buf = append(buf, '\t')
				buf = strconv.AppendInt(buf, int64(i), 10)
				buf = append(buf, '\t')
				buf = strconv.AppendInt(buf, int64(j), 10)
				buf = append(buf, '\t')
				buf = appendDepth(buf, d, intScale)
				buf = append(buf, '\n')
				if _, err := w.Write(buf); err != nil {
					return err
				}
			}
			i = j
		}
	}
	return nil
}

func writePerBase(w *bufio.Writer, depth map[string][]int, g *GenomeSize, opts Options) error {
	// Upstream `bedtools genomecov` reports per-base positions with an offset
	// of 0 in zero-based mode (-dz) and 1 otherwise (-d), and only -dz skips
	// zero-depth positions (offset = _eachBaseZeroBased ? 0 : 1; the print
	// guard is `depth>0 || !_eachBaseZeroBased`).
	zeroBased := opts.Mode == ModePerBaseNonZero
	offset := 1
	if zeroBased {
		offset = 0
	}
	intScale := opts.Scale == 1.0
	// One reusable scratch buffer for the whole pass: this path emits one line
	// per base (millions for vertebrate-scale genomes), so the fmt.Fprintf-per-
	// line allocations dominated the profile (~99% of allocations in -d/-dz).
	var buf []byte
	for _, chrom := range g.Order {
		arr := depth[chrom]
		for i, raw := range arr {
			d := scaledDepth(raw, opts.Scale)
			if zeroBased && d == 0 {
				continue
			}
			buf = buf[:0]
			buf = append(buf, chrom...)
			buf = append(buf, '\t')
			buf = strconv.AppendInt(buf, int64(i+offset), 10)
			buf = append(buf, '\t')
			buf = appendDepth(buf, d, intScale)
			buf = append(buf, '\n')
			if _, err := w.Write(buf); err != nil {
				return err
			}
		}
	}
	return nil
}

// appendDepth appends a (possibly scaled) depth to buf, emitting an integer
// when the scale is 1.0 (the common case) and a compact float otherwise,
// matching bedtools' `%g`-style formatting under `-scale`. It writes into a
// reusable buffer instead of allocating a string per call.
func appendDepth(buf []byte, d float64, intScale bool) []byte {
	if intScale {
		return strconv.AppendInt(buf, int64(d), 10)
	}
	return strconv.AppendFloat(buf, d, 'g', 10, 64)
}

// writeHistogram emits per-chromosome histograms followed by a `genome` row,
// matching the exact column layout of `bedtools genomecov`.
func writeHistogram(w *bufio.Writer, depth map[string][]int, g *GenomeSize, opts Options) error {
	genome := map[int]int{}
	genomeSize := 0
	var buf []byte
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
			buf = appendHistLine(buf[:0], chrom, d, n, chromSize, frac)
			if _, err := w.Write(buf); err != nil {
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
		buf = appendHistLine(buf[:0], "genome", d, n, genomeSize, frac)
		if _, err := w.Write(buf); err != nil {
			return err
		}
	}
	return nil
}

// appendHistLine appends a histogram row "chrom\tdepth\tcount\tsize\tfrac\n" to
// buf using strconv.Append* in place of fmt.Fprintf, matching the exact column
// layout of `bedtools genomecov`.
func appendHistLine(buf []byte, chrom string, depth, count, size int, frac float64) []byte {
	buf = append(buf, chrom...)
	buf = append(buf, '\t')
	buf = strconv.AppendInt(buf, int64(depth), 10)
	buf = append(buf, '\t')
	buf = strconv.AppendInt(buf, int64(count), 10)
	buf = append(buf, '\t')
	buf = strconv.AppendInt(buf, int64(size), 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, frac, 'g', 6, 64)
	buf = append(buf, '\n')
	return buf
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

// formatFraction renders a fraction with 6 significant digits and no
// trailing-zero noise. Upstream `bedtools genomecov` prints histogram
// fractions via C++ ostream's default precision (6), so we match that here.
func formatFraction(f float64) string {
	s := strconv.FormatFloat(f, 'g', 6, 64)
	return s
}
