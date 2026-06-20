// Package bedtobam is a pure-Go reimplementation of `bedtools bedtobam` (aka
// bedToBam): it converts BED (or BED12) records into BAM alignments against a
// header built from a genome (chrom sizes) file.
//
// The output is byte-for-byte equivalent to upstream bedtools v2.31.1 once both
// BAMs are decoded to SAM: the same @HD/@PG/@SQ header (including the literal
// "VN:V<version>" program tag and the genome-file path in each @SQ AS field),
// the same 0-based POS, MAPQ, FLAG (0, or 0x10 for "-" strand), empty SEQ/QUAL,
// and the same CIGAR — a single <len>M for plain BED, or an N/M spliced CIGAR
// derived from the BED12 blocks.
package bedtobam

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// UpstreamVersion is the bedtools version string embedded in the @PG header
// line, matching the binary this port targets.
const UpstreamVersion = "v2.31.1"

// Options configures a BED→BAM conversion. The fields mirror the upstream
// bedtobam command-line flags.
type Options struct {
	// MapQ is the MAPQ value stamped on every emitted record (0..255).
	MapQ int
	// BED12 interprets the input as BED12 and derives a spliced CIGAR from the
	// block sizes/starts.
	BED12 bool
	// Uncompressed writes an uncompressed BAM (the -ubam flag). The decoded
	// records are identical either way; only the BGZF compression level differs.
	Uncompressed bool
	// GenomeFileName is the path used verbatim in each @SQ AS: field, matching
	// upstream which records the genome file path it was given.
	GenomeFileName string
}

// Chrom is one (name, length) entry from a genome / chrom-sizes file.
type Chrom struct {
	Name   string
	Length int32
}

// Genome is the ordered list of reference sequences parsed from a genome file.
// Order is preserved exactly as in the file (upstream does not sort), since it
// drives both the @SQ order and the BAM reference IDs.
type Genome struct {
	Chroms []Chrom
	index  map[string]int
}

// ReadGenome parses a two-column (chrom<TAB or whitespace>size) genome file.
// Blank lines are skipped. The chromosome order in the file is preserved.
func ReadGenome(r io.Reader) (*Genome, error) {
	g := &Genome{index: map[string]int{}}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 1<<20)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		if strings.TrimSpace(line) == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("bedtobam: malformed genome line %q", line)
		}
		n, err := strconv.ParseInt(fields[1], 10, 32)
		if err != nil {
			return nil, fmt.Errorf("bedtobam: bad chrom size in %q: %w", line, err)
		}
		if _, dup := g.index[fields[0]]; dup {
			continue
		}
		g.index[fields[0]] = len(g.Chroms)
		g.Chroms = append(g.Chroms, Chrom{Name: fields[0], Length: int32(n)})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return g, nil
}

// Header builds the BAM/SAM header for the genome, matching upstream's
// MakeBamHeader: an @HD line, a @PG line naming BEDTools_bedToBam, and one @SQ
// line per chromosome carrying SN, AS (the genome file path) and LN.
func (g *Genome) Header(genomeFileName string) *sam.Header {
	h := &sam.Header{}
	hd := sam.HeaderLine{Tag: "HD", Fields: []sam.HeaderField{
		{Tag: "VN", Value: "1.0"},
		{Tag: "SO", Value: "unsorted"},
	}}
	pg := sam.HeaderLine{Tag: "PG", Fields: []sam.HeaderField{
		{Tag: "ID", Value: "BEDTools_bedToBam"},
		{Tag: "VN", Value: "V" + UpstreamVersion},
	}}
	h.Lines = append(h.Lines, hd, pg)
	h.HDFields = hd.Fields
	h.Programs = append(h.Programs, sam.Program{ID: "BEDTools_bedToBam", Extra: pg.Fields[1:]})
	for _, c := range g.Chroms {
		sq := sam.HeaderLine{Tag: "SQ", Fields: []sam.HeaderField{
			{Tag: "SN", Value: c.Name},
			{Tag: "AS", Value: genomeFileName},
			{Tag: "LN", Value: strconv.FormatInt(int64(c.Length), 10)},
		}}
		h.Lines = append(h.Lines, sq)
		h.Refs = append(h.Refs, sam.Reference{
			Name:   c.Name,
			Length: c.Length,
			Extra:  []sam.HeaderField{{Tag: "AS", Value: genomeFileName}},
		})
	}
	return h
}

// bedRecord is the minimal parsed view of an input BED line.
type bedRecord struct {
	chrom  string
	start  int
	end    int
	name   string
	strand string
	fields []string
}

// Run reads BED records from r, converts each to a BAM alignment against the
// genome, and writes a BAM stream to w. It returns the number of records
// written. A record whose chromosome is absent from the genome, or a line with
// fewer than four columns (no name), is a fatal error, matching upstream.
func Run(r io.Reader, w io.Writer, g *Genome, opts Options) (int, error) {
	if opts.MapQ < 0 || opts.MapQ > 255 {
		return 0, fmt.Errorf("bedtobam: MAPQ must be in range [0,255]")
	}
	bw := sam.NewBAMWriter(w)
	if err := bw.WriteHeader(g.Header(opts.GenomeFileName)); err != nil {
		return 0, err
	}

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 1<<20)
	n := 0
	lineNum := 0
	for sc.Scan() {
		lineNum++
		line := strings.TrimRight(sc.Text(), "\r")
		if line == "" || strings.HasPrefix(line, "#") ||
			strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 4 {
			return n, fmt.Errorf("Error: BED entry without name found at line: %d.  Exiting!", lineNum)
		}
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return n, fmt.Errorf("bedtobam: bad start at line %d: %w", lineNum, err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return n, fmt.Errorf("bedtobam: bad end at line %d: %w", lineNum, err)
		}
		bed := bedRecord{
			chrom:  fields[0],
			start:  start,
			end:    end,
			name:   fields[3],
			fields: fields,
		}
		if len(fields) >= 6 {
			bed.strand = fields[5]
		}
		rec, err := convert(bed, opts, lineNum)
		if err != nil {
			return n, err
		}
		// Upstream resolves the reference id via std::map::operator[], which
		// silently inserts a default value of 0 for a chromosome that is not in
		// the genome file — so an unknown chrom is written against the FIRST
		// @SQ reference rather than rejected. We reproduce that quirk for
		// byte-parity (see docs/UPSTREAM_BUGS.md); it requires at least one
		// reference in the genome.
		if _, ok := g.index[bed.chrom]; !ok {
			if len(g.Chroms) == 0 {
				return n, fmt.Errorf("bedtobam: empty genome file; cannot map chromosome %q (line %d)", bed.chrom, lineNum)
			}
			rec.RName = g.Chroms[0].Name
		}
		if err := bw.Write(rec); err != nil {
			return n, err
		}
		n++
	}
	if err := sc.Err(); err != nil {
		return n, err
	}
	if err := bw.Close(); err != nil {
		return n, err
	}
	return n, nil
}

// convert turns one BED record into a sam.Record, building either a single
// <len>M CIGAR (plain BED) or an N/M spliced CIGAR (BED12).
func convert(bed bedRecord, opts Options, lineNum int) (*sam.Record, error) {
	rec := &sam.Record{
		QName: bed.name,
		RName: bed.chrom,
		Pos:   int64(bed.start) + 1, // sam.Record uses 1-based POS.
		MapQ:  uint8(opts.MapQ),
		RNext: "*",
		PNext: 0,
		TLen:  0,
	}
	if bed.strand == "-" {
		rec.Flag |= sam.FlagReverse
	}

	if !opts.BED12 {
		rec.Cigar = sam.Cigar{sam.CigarOp(uint32(bed.end-bed.start)<<4 | sam.CigarMatch)}
		return rec, nil
	}

	if len(bed.fields) != 12 {
		return nil, fmt.Errorf("You've indicated that the input file is in BED12 format, yet the relevant fields cannot be found.  Exiting.")
	}
	blockCount, err := strconv.Atoi(bed.fields[9])
	if err != nil {
		return nil, fmt.Errorf("bedtobam: bad blockCount at line %d: %w", lineNum, err)
	}
	blockSizes, err := splitInts(bed.fields[10])
	if err != nil {
		return nil, fmt.Errorf("bedtobam: bad blockSizes at line %d: %w", lineNum, err)
	}
	blockStarts, err := splitInts(bed.fields[11])
	if err != nil {
		return nil, fmt.Errorf("bedtobam: bad blockStarts at line %d: %w", lineNum, err)
	}
	if len(blockSizes) != blockCount {
		return nil, fmt.Errorf("Error: Number of BED blocks does not match blockCount at line: %d.  Exiting!", lineNum)
	}

	var cig sam.Cigar
	addOp := func(length int, op uint32) {
		cig = append(cig, sam.CigarOp(uint32(length)<<4|op))
	}
	// Leading skip if the first block does not start at bed.start.
	if blockStarts[0] > 0 {
		addOp(blockStarts[0], sam.CigarSkipped)
	}
	for i := 0; i < blockCount-1; i++ {
		addOp(blockSizes[i], sam.CigarMatch)
		if blockStarts[i+1] > blockStarts[i]+blockSizes[i] {
			addOp(blockStarts[i+1]-(blockStarts[i]+blockSizes[i]), sam.CigarSkipped)
		}
	}
	addOp(blockSizes[blockCount-1], sam.CigarMatch)
	rec.Cigar = cig
	return rec, nil
}

// splitInts parses a comma-delimited list of ints, ignoring a trailing comma
// and empty trailing element (the BED convention).
func splitInts(s string) ([]int, error) {
	s = strings.TrimRight(s, ",")
	if s == "" {
		return nil, nil
	}
	parts := strings.Split(s, ",")
	out := make([]int, 0, len(parts))
	for _, p := range parts {
		if p == "" {
			continue
		}
		v, err := strconv.Atoi(p)
		if err != nil {
			return nil, err
		}
		out = append(out, v)
	}
	return out, nil
}
