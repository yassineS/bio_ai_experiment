// Package bedcluster ports `bedtools cluster`: it groups overlapping
// intervals into clusters and tags each input record with a cluster ID.
//
// Unlike `bedmerge`, no records are collapsed — every input record appears
// in the output, with an extra trailing column holding its cluster ID.
//
// Two intervals end up in the same cluster when they are within MaxDistance
// bp of each other (default 0 = require overlap, matching upstream's
// `-d 0` default). With StrandSpec, only same-strand intervals cluster.
package bedcluster

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/cppsort"
)

// Options configures Cluster.
type Options struct {
	// MaxDistance is the largest gap (in bp) between two intervals that still
	// puts them in the same cluster. Default 0 = require overlap (book-ends
	// included: an interval ending at 100 and one starting at 100 cluster
	// when MaxDistance=0, matching upstream).
	MaxDistance int
	// StrandSpec, when true, only clusters intervals with the same strand.
	StrandSpec bool
}

// row is one parsed input line. We keep the raw (CR/LF-trimmed) line bytes so
// the record can be re-emitted verbatim — the original column count and exact
// tab-delimited text — with only the cluster ID appended, avoiding a per-record
// split-then-join round trip.
type row struct {
	chrom  string
	start  int
	end    int
	strand string
	raw    []byte
}

// Cluster reads BED records from r, assigns each one a cluster ID, and writes
// each record back to w with the cluster ID appended as the final column.
//
// Records are emitted in genomic-sort order (chrom, then — for StrandSpec —
// strand, then start, then end), which matches upstream's behaviour: the
// upstream tool requires sorted input and groups by strand internally when
// `-s` is set, producing per-strand contiguous output.
//
// Returns the number of records written.
func Cluster(r io.Reader, w io.Writer, opts Options) (int, error) {
	rows, err := readAll(r)
	if err != nil {
		return 0, err
	}
	if len(rows) == 0 {
		return 0, nil
	}

	// Reproduce upstream clusterBed exactly (clusterBed.cpp): records are
	// bucketed by chromosome (std::map, lexicographic) and each bucket is
	// std::sort'd by start ALONE (introsort — equal-start records land in its
	// artifact order, not input order), via loadBedFileIntoMapNoBin. Then, with
	// -s, two passes emit the "+" then the "-" records of each chromosome in
	// that introsorted order (other strands, e.g. ".", are dropped — upstream's
	// strands vector is exactly {"+","-"}). Without -s the introsorted bucket is
	// emitted as-is. assignClusters then sweeps; it already opens a new cluster
	// on a strand change, matching upstream's per-strand `end = -1` reset.
	rows = clusterOrder(rows, opts.StrandSpec)

	clusterIDs := assignClusters(rows, opts)

	bw := bufio.NewWriter(w)
	defer bw.Flush()
	// Emit each record directly from its retained line bytes: original text,
	// a tab, the cluster ID (formatted into a reusable scratch buffer), and a
	// newline. This avoids the per-record split-then-Join and Itoa allocations
	// the previous []string round trip incurred.
	var scratch []byte
	for i := range rows {
		if _, err := bw.Write(rows[i].raw); err != nil {
			return 0, err
		}
		if err := bw.WriteByte('\t'); err != nil {
			return 0, err
		}
		scratch = strconv.AppendInt(scratch[:0], int64(clusterIDs[i]), 10)
		if _, err := bw.Write(scratch); err != nil {
			return 0, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return 0, err
		}
	}
	return len(rows), nil
}

// assignClusters returns a slice (parallel to the sorted `rows`) of cluster
// IDs. Cluster IDs are 1-based and assigned in walk order.
// clusterOrder reproduces upstream clusterBed's record ordering: per-chromosome
// (lexicographic) introsort by start, then — when strandSpec is set — a stable
// "+"-then-"-" partition within each chromosome, dropping records on any other
// strand (matching upstream's two fixed strand passes).
func clusterOrder(rows []row, strandSpec bool) []row {
	// Without -s, upstream uses the streaming GetNextBed(force_sorted) path,
	// which processes records in (chrom, start) order. We sort STABLY by
	// (chrom, start) — equal-start records keep input order, matching the
	// stream — which is a lenient fix-on-port (upstream errors on unsorted input
	// rather than sorting it). Only the -s path goes through
	// loadBedFileIntoMapNoBin's per-chromosome introsort, whose unstable
	// equal-start order we must reproduce.
	if !strandSpec {
		sort.SliceStable(rows, func(i, j int) bool {
			if rows[i].chrom != rows[j].chrom {
				return rows[i].chrom < rows[j].chrom
			}
			return rows[i].start < rows[j].start
		})
		return rows
	}
	buckets := map[string][]row{}
	var chroms []string
	for _, r := range rows {
		if _, ok := buckets[r.chrom]; !ok {
			chroms = append(chroms, r.chrom)
		}
		buckets[r.chrom] = append(buckets[r.chrom], r)
	}
	sort.Strings(chroms)
	out := make([]row, 0, len(rows))
	for _, chrom := range chroms {
		b := buckets[chrom]
		cppsort.Sort(b, func(x, y row) bool { return x.start < y.start })
		if strandSpec {
			for _, want := range [2]string{"+", "-"} {
				for _, r := range b {
					if r.strand == want {
						out = append(out, r)
					}
				}
			}
		} else {
			out = append(out, b...)
		}
	}
	return out
}

func assignClusters(rows []row, opts Options) []int {
	ids := make([]int, len(rows))
	if len(rows) == 0 {
		return ids
	}
	curID := 1
	curChrom := rows[0].chrom
	curStrand := rows[0].strand
	curEnd := rows[0].end
	ids[0] = curID

	for i := 1; i < len(rows); i++ {
		r := rows[i]
		newCluster := false
		switch {
		case r.chrom != curChrom:
			newCluster = true
		case opts.StrandSpec && r.strand != curStrand:
			newCluster = true
		case r.start > curEnd+opts.MaxDistance:
			// Strict gap: a gap > MaxDistance opens a new cluster.
			newCluster = true
		}
		if newCluster {
			curID++
			curChrom = r.chrom
			curStrand = r.strand
			curEnd = r.end
		} else if r.end > curEnd {
			curEnd = r.end
		}
		ids[i] = curID
	}
	return ids
}

// Header/track prefixes recognised on text input, kept as byte slices so the
// read loop can test the scanner buffer without allocating a string per line.
var (
	trackPrefix   = []byte("track")
	browserPrefix = []byte("browser")
)

func readAll(r io.Reader) ([]row, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var rows []row
	var ci chromInterner
	for sc.Scan() {
		// Trim trailing CR/LF from the scanner's buffer in place. The scanner
		// already strips the line's '\n'; this also drops a trailing '\r' so the
		// emitted text matches the previous TrimRight(line, "\r\n") behaviour.
		line := bytes.TrimRight(sc.Bytes(), "\r\n")
		trimmed := bytes.TrimSpace(line)
		if len(trimmed) == 0 || trimmed[0] == '#' ||
			bytes.HasPrefix(trimmed, trackPrefix) || bytes.HasPrefix(trimmed, browserPrefix) {
			continue
		}
		// Tokenize on tabs without allocating a []string: collect up to the
		// first six column spans, which is all the parse needs (chrom, start,
		// end, strand).
		var cols [6][]byte
		n := 0
		begin := 0
		for i := 0; i <= len(line); i++ {
			if i == len(line) || line[i] == '\t' {
				if n < len(cols) {
					cols[n] = line[begin:i]
				}
				n++
				begin = i + 1
			}
		}
		if n < 3 {
			return nil, fmt.Errorf("BED record needs >=3 columns, got %d in line: %q", n, line)
		}
		start, err := parseIntBytes(cols[1])
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %q: %v", cols[1], err)
		}
		end, err := parseIntBytes(cols[2])
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %q: %v", cols[2], err)
		}
		strand := ""
		if n >= 6 {
			strand = internStrand(cols[5])
		}
		// The scanner reuses its buffer, so retain a private copy of the line.
		raw := make([]byte, len(line))
		copy(raw, line)
		rows = append(rows, row{
			chrom:  ci.intern(cols[0]),
			start:  start,
			end:    end,
			strand: strand,
			raw:    raw,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("reading BED: %w", err)
	}
	return rows, nil
}

// chromInterner caches the most recently interned chromosome name so a run of
// records sharing a chromosome (the norm for sorted BED input) allocates the
// name string once rather than per record.
type chromInterner struct {
	last string
}

func (c *chromInterner) intern(b []byte) string {
	if c.last != "" && c.last == string(b) {
		return c.last
	}
	s := string(b)
	c.last = s
	return s
}

// internStrand maps the byte form of a strand column to a shared string,
// allocating nothing for the values that occur in practice.
func internStrand(b []byte) string {
	switch string(b) {
	case "":
		return ""
	case "+":
		return "+"
	case "-":
		return "-"
	case ".":
		return "."
	}
	return string(b)
}

// parseIntBytes parses a base-10 integer from a byte slice, accepting an
// optional leading sign, without allocating a string on the fast path. It
// trims surrounding ASCII whitespace first to match the previous
// strconv.Atoi(strings.TrimSpace(...)) behaviour.
func parseIntBytes(b []byte) (int, error) {
	b = bytes.TrimSpace(b)
	if len(b) == 0 {
		return 0, fmt.Errorf("empty number")
	}
	i := 0
	neg := false
	if b[0] == '+' || b[0] == '-' {
		neg = b[0] == '-'
		i++
	}
	if i >= len(b) {
		return 0, fmt.Errorf("invalid number %q", b)
	}
	val := 0
	for ; i < len(b); i++ {
		c := b[i]
		if c < '0' || c > '9' {
			return 0, fmt.Errorf("invalid number %q", b)
		}
		val = val*10 + int(c-'0')
	}
	if neg {
		val = -val
	}
	return val, nil
}
