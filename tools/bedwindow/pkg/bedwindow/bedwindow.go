// Package bedwindow ports `bedtools window`: it computes the overlap of A
// against B after expanding each B interval by a window (-w, -l, -r). It
// otherwise behaves like bedintersect with the same set of writer modes
// (default = intersection coords, -wa = original A, -wb = original B,
// -c = count of B overlaps for each A, -v = invert).
//
// This implementation expands B intervals at load time (clipping at 0 on
// the low end) and then runs the same overlap finding logic as
// bedintersect via per-chromosome interval trees from the shared
// pkg/bioformats/bed package.
package bedwindow

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options configures Window.
type Options struct {
	// Left is the bp to extend each B interval to the left (lower coordinates).
	// Negative values shrink. Result is clipped at 0.
	Left int
	// Right is the bp to extend each B interval to the right.
	Right int

	// StrandSpec, when true, requires same-strand overlap (matches `-sm`).
	StrandSpec bool
	// InverseStrand, when true, requires opposite-strand overlap (matches `-sw`
	// in the upstream sense of "different strand" — naming follows our
	// project convention).
	InverseStrand bool

	// WriteA — emit the original A record.
	WriteA bool
	// WriteB — emit the original B record (overrides WriteA).
	WriteB bool
	// WriteAB — emit `A<TAB>B` for each overlap (matches `bedtools window` default
	// when neither -u/-c/-v is set: in upstream, the default emits A and B
	// concatenated for each overlap).
	WriteAB bool
	// Count — emit `A<TAB>count` instead of one row per overlap.
	Count bool
	// Invert — emit only A records with no B overlap.
	Invert bool
	// MinOverlap — minimum bp overlap required (default 1).
	MinOverlap int
}

// ExpandedB is a B record with its expanded coordinates, kept alongside its
// original coordinates for the -wb writer mode.
type ExpandedB struct {
	Orig     *bed.Record
	Expanded *bed.Record
}

// Window reads A from aR, B from bR, and writes results to w.
func Window(aR, bR io.Reader, w io.Writer, opts Options) (int, error) {
	if opts.MinOverlap < 1 {
		opts.MinOverlap = 1
	}
	if opts.StrandSpec && opts.InverseStrand {
		return 0, fmt.Errorf("StrandSpec (-sm) and InverseStrand (-sw) are mutually exclusive")
	}

	bRecs, err := readAll(bR)
	if err != nil {
		return 0, fmt.Errorf("reading B: %w", err)
	}

	// Expand each B record by [-Left, +Right].
	expanded := make([]*ExpandedB, 0, len(bRecs))
	for _, b := range bRecs {
		exStart := b.ChromStart - opts.Left
		if exStart < 0 {
			exStart = 0
		}
		exEnd := b.ChromEnd + opts.Right
		if exEnd <= exStart {
			// Window collapsed the interval to nothing; skip.
			continue
		}
		ex := &bed.Record{
			Chrom:      b.Chrom,
			ChromStart: exStart,
			ChromEnd:   exEnd,
			Strand:     b.Strand,
		}
		expanded = append(expanded, &ExpandedB{Orig: b, Expanded: ex})
	}

	// Sort and tree-index by chrom/start of the expanded coordinates.
	sort.SliceStable(expanded, func(i, j int) bool {
		if expanded[i].Expanded.Chrom != expanded[j].Expanded.Chrom {
			return expanded[i].Expanded.Chrom < expanded[j].Expanded.Chrom
		}
		return expanded[i].Expanded.ChromStart < expanded[j].Expanded.ChromStart
	})
	byChrom := map[string][]*bed.Record{}
	exToOrig := map[*bed.Record]*bed.Record{}
	for _, ex := range expanded {
		byChrom[ex.Expanded.Chrom] = append(byChrom[ex.Expanded.Chrom], ex.Expanded)
		exToOrig[ex.Expanded] = ex.Orig
	}
	trees := map[string]*bed.IntervalTree{}
	for chrom, list := range byChrom {
		trees[chrom] = bed.NewIntervalTree(list)
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

		var hits []*bed.Record
		if tree, ok := trees[recA.Chrom]; ok {
			candidates := tree.Query(recA)
			for _, ex := range candidates {
				orig := exToOrig[ex]
				if !overlapPasses(recA, ex, orig, opts) {
					continue
				}
				hits = append(hits, orig)
			}
		}

		if opts.Invert {
			if len(hits) == 0 {
				if _, err := fmt.Fprintln(bw, formatRecord(recA)); err != nil {
					return written, err
				}
				written++
			}
			continue
		}
		if opts.Count {
			out := formatRecord(recA) + "\t" + strconv.Itoa(len(hits))
			if _, err := fmt.Fprintln(bw, out); err != nil {
				return written, err
			}
			written++
			continue
		}
		if len(hits) == 0 {
			continue
		}
		for _, hit := range hits {
			var out string
			switch {
			case opts.WriteAB:
				out = formatRecord(recA) + "\t" + formatRecord(hit)
			case opts.WriteB:
				out = formatRecord(hit)
			case opts.WriteA:
				out = formatRecord(recA)
			default:
				// Default for `bedtools window` (no writer flag) is A<TAB>B
				// (matches upstream behaviour). We adopt the same default.
				out = formatRecord(recA) + "\t" + formatRecord(hit)
			}
			if _, err := fmt.Fprintln(bw, out); err != nil {
				return written, err
			}
			written++
		}
	}
	return written, nil
}

func overlapPasses(a, ex, _ *bed.Record, opts Options) bool {
	overlap := overlapBP(a, ex)
	if overlap < opts.MinOverlap {
		return false
	}
	if opts.StrandSpec {
		if a.Strand == "" || ex.Strand == "" || a.Strand != ex.Strand {
			return false
		}
	}
	if opts.InverseStrand {
		if a.Strand == "" || ex.Strand == "" || a.Strand == ex.Strand {
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

// readAll loads BED records from r using the shared bed.Reader.
func readAll(r io.Reader) ([]*bed.Record, error) {
	br := bed.NewReader(r)
	var out []*bed.Record
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		out = append(out, rec)
	}
	return out, nil
}

// formatRecord re-renders a *bed.Record back to a tab-separated line. The
// number of columns matches what was originally read (we infer that from
// which fields are populated; ExtraFields are appended verbatim).
func formatRecord(r *bed.Record) string {
	cols := []string{r.Chrom, strconv.Itoa(r.ChromStart), strconv.Itoa(r.ChromEnd)}
	if r.Name != "" {
		cols = append(cols, r.Name)
	}
	if r.Score != 0 || r.Strand != "" {
		// Score may be the literal 0; emit "." for empty name slot.
		if r.Name == "" {
			cols = append(cols, ".")
		}
		cols = append(cols, strconv.Itoa(r.Score))
	}
	if r.Strand != "" {
		cols = append(cols, r.Strand)
	}
	if len(r.ExtraFields) > 0 {
		cols = append(cols, r.ExtraFields...)
	}
	return strings.Join(cols, "\t")
}
