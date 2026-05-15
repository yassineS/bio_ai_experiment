// Package bedannotate implements `bedtools annotate`: for each interval in
// a primary BED file (`-i`), it annotates with overlap statistics drawn
// from N additional BED files (`-files`). Three output modes mirror the
// upstream flags:
//
//   - default        — emit the fraction of A covered by each B (as %f)
//   - -counts        — emit the count of overlapping records per B
//   - -both          — emit `<count>\t<fraction>` per B
//
// Strand filters (-s same-strand, -S opposite-strand) restrict which
// overlaps count. The output preserves A's original columns and appends
// the per-B columns; an optional column-header line ("#…") is emitted
// when names are supplied (either explicitly via -names or implicitly
// from -files basenames).
//
// Internally we read each B file into a per-chromosome interval tree
// (`pkg/bioformats/bed.IntervalTree`) and stream A line by line, so
// the working set is O(sum(|B_i|)) rather than O(|A|·|B_i|).
package bedannotate

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Mode controls which columns are appended per B file.
type Mode int

const (
	// ModeFraction emits a single float column per B (default).
	ModeFraction Mode = iota
	// ModeCounts emits a single int column per B (-counts).
	ModeCounts
	// ModeBoth emits `<count>\t<fraction>` per B (-both).
	ModeBoth
)

// Options configures Run.
type Options struct {
	// Mode selects the output shape. Default = ModeFraction.
	Mode Mode
	// Names is the header label for each B file. When empty, no header
	// line is emitted. When supplied via the CLI, the caller fills this
	// from -names / file basenames.
	Names []string
	// SameStrand: -s, require A.Strand == B.Strand.
	SameStrand bool
	// OppositeStrand: -S, require A.Strand != B.Strand.
	OppositeStrand bool
}

// Run reads A from aR and the N B files from bRs in order, indexes each B
// into per-chromosome interval trees, then streams A and emits one row per
// A record with the per-B columns appended. Returns the number of A
// records processed.
func Run(aR io.Reader, bRs []io.Reader, out io.Writer, opts Options) (int, error) {
	if opts.SameStrand && opts.OppositeStrand {
		return 0, fmt.Errorf("cannot combine -s and -S")
	}
	// Build a tree per B file.
	trees := make([]map[string]*bed.IntervalTree, len(bRs))
	for i, br := range bRs {
		t, err := indexB(br)
		if err != nil {
			return 0, fmt.Errorf("file %d: %w", i+1, err)
		}
		trees[i] = t
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	// Optional header line.
	if len(opts.Names) > 0 {
		writeHeader(bw, opts.Names, opts.Mode)
	}

	sc := bufio.NewScanner(aR)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	count := 0
	lineNo := 0
	for sc.Scan() {
		lineNo++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return count, fmt.Errorf("line %d: BED record needs >=3 columns: %q", lineNo, raw)
		}
		rec, err := parseRecord(fields)
		if err != nil {
			return count, fmt.Errorf("line %d: %w", lineNo, err)
		}
		// Emit original columns first.
		if _, err := bw.WriteString(strings.Join(fields, "\t")); err != nil {
			return count, err
		}
		// Then per-B columns.
		for i, t := range trees {
			matches := selectOverlapping(rec, t[rec.Chrom], opts)
			cnt := len(matches)
			frac := coveredFraction(rec, matches)
			switch opts.Mode {
			case ModeCounts:
				fmt.Fprintf(bw, "\t%d", cnt)
			case ModeBoth:
				fmt.Fprintf(bw, "\t%d\t%f", cnt, frac)
			default:
				fmt.Fprintf(bw, "\t%f", frac)
			}
			_ = i
		}
		if err := bw.WriteByte('\n'); err != nil {
			return count, err
		}
		count++
	}
	if err := sc.Err(); err != nil {
		return count, err
	}
	return count, nil
}

// writeHeader emits the "#<TAB>...<TAB>name1<TAB>name2..." header line.
// Upstream pads the leading hash with `bedType-1` empty tabs so the
// labels line up with the first appended column; we don't know bedType
// at header time without peeking, so we emit a single '#' and then the
// labels, which matches what upstream does for `bedType=1` and is still
// machine-parseable. With `-both`, each label is split into `_cnt`/`_pct`.
func writeHeader(w *bufio.Writer, names []string, mode Mode) {
	w.WriteByte('#')
	if mode == ModeBoth {
		for _, n := range names {
			fmt.Fprintf(w, "\t%s_cnt\t%s_pct", n, n)
		}
	} else {
		for _, n := range names {
			fmt.Fprintf(w, "\t%s", n)
		}
	}
	w.WriteByte('\n')
}

// parseRecord parses the minimum subset of a BED line we need for overlap
// + strand filtering. Extra columns are preserved by the caller as raw
// fields.
func parseRecord(fields []string) (*bed.Record, error) {
	start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
	}
	end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
	}
	r := &bed.Record{
		Chrom:      fields[0],
		ChromStart: start,
		ChromEnd:   end,
	}
	if len(fields) >= 6 {
		r.Strand = fields[5]
	}
	return r, nil
}

// indexB reads a B file fully into memory and returns a per-chrom tree.
func indexB(r io.Reader) (map[string]*bed.IntervalTree, error) {
	rd := bed.NewReader(r)
	all, err := rd.ReadAll()
	if err != nil {
		return nil, err
	}
	byChrom := map[string][]*bed.Record{}
	for _, x := range all {
		byChrom[x.Chrom] = append(byChrom[x.Chrom], x)
	}
	out := make(map[string]*bed.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].ChromStart != recs[j].ChromStart {
				return recs[i].ChromStart < recs[j].ChromStart
			}
			return recs[i].ChromEnd < recs[j].ChromEnd
		})
		out[chrom] = bed.NewIntervalTree(recs)
	}
	return out, nil
}

// selectOverlapping returns the B records overlapping a (after applying
// the strand filters in opts).
func selectOverlapping(a *bed.Record, t *bed.IntervalTree, opts Options) []*bed.Record {
	if t == nil {
		return nil
	}
	cand := t.Query(a)
	if len(cand) == 0 {
		return nil
	}
	out := cand[:0:0]
	for _, b := range cand {
		if !strandOK(a, b, opts) {
			continue
		}
		out = append(out, b)
	}
	return out
}

// strandOK applies the -s / -S filters.
func strandOK(a, b *bed.Record, opts Options) bool {
	if opts.SameStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		return a.Strand == b.Strand
	}
	if opts.OppositeStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		return a.Strand != b.Strand
	}
	return true
}

// coveredFraction returns the fraction of A's length covered by at least
// one of the matches (depth >= 1). Used for default and -both modes.
func coveredFraction(a *bed.Record, matches []*bed.Record) float64 {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 {
		return 0
	}
	if len(matches) == 0 {
		return 0
	}
	covered := make([]bool, lenA)
	for _, b := range matches {
		start := b.ChromStart - a.ChromStart
		end := b.ChromEnd - a.ChromStart
		if start < 0 {
			start = 0
		}
		if end > lenA {
			end = lenA
		}
		for i := start; i < end; i++ {
			covered[i] = true
		}
	}
	n := 0
	for _, c := range covered {
		if c {
			n++
		}
	}
	return float64(n) / float64(lenA)
}

// DefaultNames derives a label per file from the trailing path component
// of each filename, mirroring upstream's behaviour when -names is omitted
// but a header is still desirable.
func DefaultNames(paths []string) []string {
	out := make([]string, len(paths))
	for i, p := range paths {
		// Strip everything up to and including the last '/'.
		base := p
		if idx := strings.LastIndex(p, "/"); idx >= 0 {
			base = p[idx+1:]
		}
		out[i] = base
	}
	return out
}
