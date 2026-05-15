// Package bedtag implements `bedtools tag`: it annotates each record in A
// with the names (or another column) of all overlapping records in B.
//
// The tool emits each A record verbatim, with a single appended TSV column
// containing the comma-separated list of tags from B records that overlap A.
// When no B record overlaps an A record, an empty column is appended.
//
// Options:
//
//   - Names: a comma-separated label list to use INSTEAD of B's column 4.
//     The i-th name labels the i-th input B file. Requires Files (multi-B).
//   - Labels: when set, each tag is prefixed with `<filename>=` so callers
//     can tell which B file contributed which tag.
//   - StrandSpec: if true, only same-strand overlaps tag.
//   - InverseStrand: if true, only opposite-strand overlaps tag.
//   - MinOverlap, FractionA, FractionB: same overlap-requirement semantics
//     as bedintersect.
package bedtag

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options configures Tag.
type Options struct {
	// TagColumn is the 1-based column number from each B record that
	// supplies the tag value. Default: 4 (the BED "name" column).
	TagColumn int

	// Names, when non-empty, REPLACES B-derived tags with one fixed name
	// per input B file. Length must equal len(BFiles) when both are set.
	// Mutually exclusive with TagColumn.
	Names []string

	// Labels, when true, prefixes each tag with `<source>=` where source
	// is the B-file label (the file's path, or the matching Names entry
	// when Names is set).
	Labels bool

	// StrandSpec: only same-strand B records contribute tags.
	StrandSpec bool
	// InverseStrand: only opposite-strand B records contribute tags.
	// Mutually exclusive with StrandSpec.
	InverseStrand bool

	// MinOverlap is the minimum bp overlap required (default: 1).
	MinOverlap int
	// FractionA is the minimum fraction of A that must overlap.
	FractionA float64
	// FractionB is the minimum fraction of B that must overlap.
	FractionB float64
}

// Source represents a single B input together with its label.
type Source struct {
	Name   string    // label (used when Labels or Names is set)
	Reader io.Reader // stream
}

// Tag reads A from aR and one or more B sources, then writes each A record
// to w with a tag column appended. Returns the number of A records written.
func Tag(aR io.Reader, bSources []Source, w io.Writer, opts Options) (int, error) {
	if opts.MinOverlap < 1 {
		opts.MinOverlap = 1
	}
	if opts.TagColumn == 0 {
		opts.TagColumn = 4
	}
	if len(opts.Names) > 0 && len(opts.Names) != len(bSources) {
		return 0, fmt.Errorf("len(Names)=%d must equal number of B sources (%d)",
			len(opts.Names), len(bSources))
	}
	if opts.StrandSpec && opts.InverseStrand {
		return 0, fmt.Errorf("StrandSpec and InverseStrand are mutually exclusive")
	}

	// Build per-source per-chromosome interval-tree indexes.
	type bIndex struct {
		label string
		// rec stores the loaded record together with its raw fields so we can
		// extract an arbitrary column for the tag value.
		byChrom map[string]*bed.IntervalTree
		records []*recordWithFields
		fixed   string // when Names is set, every record has this fixed tag.
	}
	indexes := make([]bIndex, len(bSources))
	for i, src := range bSources {
		recs, err := readWithFields(src.Reader)
		if err != nil {
			return 0, fmt.Errorf("reading B source %s: %w", src.Name, err)
		}
		// Sort by chrom/start for tree building.
		sort.SliceStable(recs, func(a, b int) bool {
			if recs[a].rec.Chrom != recs[b].rec.Chrom {
				return recs[a].rec.Chrom < recs[b].rec.Chrom
			}
			return recs[a].rec.ChromStart < recs[b].rec.ChromStart
		})
		byChrom := map[string][]*bed.Record{}
		for _, rwf := range recs {
			byChrom[rwf.rec.Chrom] = append(byChrom[rwf.rec.Chrom], rwf.rec)
		}
		trees := map[string]*bed.IntervalTree{}
		for chrom, list := range byChrom {
			trees[chrom] = bed.NewIntervalTree(list)
		}
		idx := bIndex{
			label:   src.Name,
			byChrom: trees,
			records: recs,
		}
		if len(opts.Names) > 0 {
			idx.fixed = opts.Names[i]
		}
		indexes[i] = idx
	}

	// Build a recordPtr -> raw fields map for each source so column extraction is O(1).
	type fieldMap = map[*bed.Record][]string
	fieldMaps := make([]fieldMap, len(indexes))
	for i, idx := range indexes {
		m := make(fieldMap, len(idx.records))
		for _, rwf := range idx.records {
			m[rwf.rec] = rwf.fields
		}
		fieldMaps[i] = m
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	aReader := bed.NewReader(aR)
	written := 0
	for {
		recA, err := aReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return written, fmt.Errorf("reading A: %w", err)
		}
		var tags []string
		for i, idx := range indexes {
			tree, ok := idx.byChrom[recA.Chrom]
			if !ok {
				continue
			}
			candidates := tree.Query(recA)
			for _, b := range candidates {
				if !overlapPasses(recA, b, opts) {
					continue
				}
				var tag string
				if idx.fixed != "" {
					tag = idx.fixed
				} else {
					tag = extractColumn(fieldMaps[i][b], opts.TagColumn)
				}
				if opts.Labels && idx.fixed == "" {
					// When fixed (Names) is set, the label IS the tag — no
					// further prefixing.
					tag = idx.label + "=" + tag
				}
				tags = append(tags, tag)
			}
		}
		// Re-emit A with the tag column appended.
		out := emitRecord(recA, strings.Join(tags, ","))
		if _, err := fmt.Fprintln(bw, out); err != nil {
			return written, err
		}
		written++
	}
	return written, nil
}

// emitRecord prints `recA` as a tab-separated line with `tagCol` appended
// as the last column. We rebuild the line from recA's normal columns up to
// strand, and any ExtraFields, then append the tag.
func emitRecord(r *bed.Record, tagCol string) string {
	cols := []string{r.Chrom, strconv.Itoa(r.ChromStart), strconv.Itoa(r.ChromEnd)}
	hasName := r.Name != ""
	hasScore := r.Score != 0 || hasName // emit zero-score when there are downstream cols
	hasStrand := r.Strand != ""
	if hasName || hasScore || hasStrand {
		if hasName {
			cols = append(cols, r.Name)
		} else if hasScore || hasStrand {
			cols = append(cols, ".")
		}
		if hasScore || hasStrand {
			cols = append(cols, strconv.Itoa(r.Score))
		}
		if hasStrand {
			cols = append(cols, r.Strand)
		}
	}
	cols = append(cols, r.ExtraFields...)
	cols = append(cols, tagCol)
	return strings.Join(cols, "\t")
}

// recordWithFields keeps a parsed *bed.Record together with its raw input
// fields for arbitrary-column extraction.
type recordWithFields struct {
	rec    *bed.Record
	fields []string
}

func readWithFields(r io.Reader) ([]*recordWithFields, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []*recordWithFields
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("BED record needs >=3 columns, got %d in line: %q", len(fields), line)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		rec := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) >= 4 {
			rec.Name = fields[3]
		}
		if len(fields) >= 6 {
			rec.Strand = fields[5]
		}
		out = append(out, &recordWithFields{rec: rec, fields: fields})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// extractColumn returns the 1-based column from fields, or "." when missing.
func extractColumn(fields []string, col int) string {
	if col < 1 || col > len(fields) {
		return "."
	}
	v := fields[col-1]
	if v == "" {
		return "."
	}
	return v
}

// overlapPasses returns true when a/b overlap satisfies all overlap criteria
// (min length, fraction A, fraction B, strand requirements).
func overlapPasses(a, b *bed.Record, opts Options) bool {
	overlap := overlapBP(a, b)
	if overlap < opts.MinOverlap {
		return false
	}
	if opts.FractionA > 0 {
		lenA := a.ChromEnd - a.ChromStart
		if lenA == 0 || float64(overlap)/float64(lenA) < opts.FractionA {
			return false
		}
	}
	if opts.FractionB > 0 {
		lenB := b.ChromEnd - b.ChromStart
		if lenB == 0 || float64(overlap)/float64(lenB) < opts.FractionB {
			return false
		}
	}
	if opts.StrandSpec {
		if a.Strand == "" || b.Strand == "" || a.Strand != b.Strand {
			return false
		}
	}
	if opts.InverseStrand {
		if a.Strand == "" || b.Strand == "" || a.Strand == b.Strand {
			return false
		}
	}
	return true
}

func overlapBP(a, b *bed.Record) int {
	s := a.ChromStart
	if b.ChromStart > s {
		s = b.ChromStart
	}
	e := a.ChromEnd
	if b.ChromEnd < e {
		e = b.ChromEnd
	}
	if e <= s {
		return 0
	}
	return e - s
}
