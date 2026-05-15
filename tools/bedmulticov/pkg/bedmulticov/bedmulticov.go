// Package bedmulticov implements `bedtools multicov`: for each interval in
// a primary BED file (`-bed`), it reports the count of overlapping records
// from each of N input files (`-files` / `-bams`). The output is the
// original A columns followed by N integer columns — one per input file —
// holding the overlap count.
//
// Upstream supports BED *and* indexed BAM inputs. This port now supports
// both. BAM inputs are decoded via pkg/bioformats/sam.NewBAMReader and the
// MAPQ filter (`-q`) and per-position depth cap (`-D`) are honoured. CRAM
// remains deferred (we don't have a CRAM reader yet — see
// docs/CRAM_DESIGN.md); the CLI surfaces a clear error in that case.
//
// Internally each input file is loaded into a per-chromosome interval
// tree (`pkg/bioformats/bed.IntervalTree`), and the A file is streamed
// line-by-line. Optional strand filters (-s same / -S opposite),
// fraction-of-A (-f), fraction-of-B (-F), and reciprocal (-r) thresholds
// mirror upstream's semantics.
//
// `-split` (`Options.Split`) toggles BAM CIGAR `N`-op block-aware coverage:
// instead of treating each alignment as a single span covering its full
// reference footprint, the alignment is decomposed into its contiguous
// M/=/X reference-consuming runs and an A interval is counted once iff
// *any* block has positive overlap. When `-f` is also set, upstream's
// semantics divide the total overlap by the sum of the BAM blocks'
// lengths (NOT by A's length) and use a strict `>` comparison — a quirk
// of bedtools 2.x that this port preserves byte-for-byte. See
// `reference_code/bedtools/src/multiBamCov/multiBamCov.cpp::FindBlockedOverlaps`.
package bedmulticov

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	bedpkg "github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// Options configures Run.
type Options struct {
	// FractionA is the minimum fraction of A that must be covered by a
	// single B record for it to count as overlapping (mirrors `-f`).
	// 0 = any positive overlap counts. Range: (0, 1].
	FractionA float64
	// FractionB is the analogous minimum fraction of B covered by A
	// (mirrors `-F`). 0 = unconstrained.
	FractionB float64
	// Reciprocal mirrors `-r`: when set together with FractionA, the
	// same threshold is also applied to B (equivalent to FractionB=FractionA).
	Reciprocal bool
	// SameStrand mirrors `-s`: only count B records on the same strand as A.
	SameStrand bool
	// OppositeStrand mirrors `-S`: only count B records on the opposite
	// strand from A.
	OppositeStrand bool
	// MinMAPQ mirrors `-q`: BAM alignments with MAPQ below this threshold
	// are skipped during indexing. Ignored for BED inputs.
	MinMAPQ int
	// MaxDepth mirrors `-D`: cap the reported count per A interval per BAM
	// input at this many overlapping alignments. 0 disables the cap; the
	// default in upstream (and in the CLI wrapper) is 64000. Ignored for
	// BED inputs.
	MaxDepth int
	// Split mirrors `-split`: when set, BAM alignments are decomposed by
	// CIGAR — only the contiguous M/=/X reference-consuming runs count as
	// coverage, and any spanning N-op gap is skipped. The alignment is
	// counted at most once per A interval even if multiple of its blocks
	// overlap. When combined with FractionA, the threshold is applied to
	// `total_block_overlap / sum_of_BAM_block_lengths` using a strict `>`
	// comparison (mirrors upstream's exact arithmetic). Ignored for BED
	// inputs.
	Split bool
}

// SourceKind tags an input as BED or BAM (CRAM is rejected at the CLI layer).
type SourceKind int

const (
	// SourceBED is a plain BED file (any number of columns ≥3).
	SourceBED SourceKind = iota
	// SourceBAM is a BGZF-wrapped BAM file (decoded via
	// pkg/bioformats/sam.NewBAMReader). Each primary alignment contributes
	// one interval over [Pos-1, Pos-1+ReferenceLength()) on its reference.
	SourceBAM
)

// Source pairs an io.Reader with its file-format tag so a single Run call
// can mix BED and BAM inputs. The order of Sources in a slice determines
// the order of count columns in the output.
type Source struct {
	Reader io.Reader
	Kind   SourceKind
}

// Run is the BED-only convenience entry point: every reader in bRs is
// treated as a BED file. Use RunSources to mix BAM inputs in.
func Run(aR io.Reader, bRs []io.Reader, out io.Writer, opts Options) (int, error) {
	srcs := make([]Source, len(bRs))
	for i, br := range bRs {
		srcs[i] = Source{Reader: br, Kind: SourceBED}
	}
	return RunSources(aR, srcs, out, opts)
}

// RunSources reads A from aR and the N B inputs from srcs in order. Each
// input is indexed per chromosome (BED records are read with the bed
// package; BAM records are decoded with pkg/bioformats/sam and -q MAPQ
// filtered up front). RunSources then streams A and emits one row per A
// record with one count column appended per source. -D, if set, caps the
// reported per-A count per BAM input. Returns the number of A records
// processed.
func RunSources(aR io.Reader, srcs []Source, out io.Writer, opts Options) (int, error) {
	if opts.SameStrand && opts.OppositeStrand {
		return 0, fmt.Errorf("cannot combine -s and -S")
	}
	if opts.FractionA < 0 || opts.FractionA > 1 {
		return 0, fmt.Errorf("-f must be in [0,1], got %g", opts.FractionA)
	}
	if opts.FractionB < 0 || opts.FractionB > 1 {
		return 0, fmt.Errorf("-F must be in [0,1], got %g", opts.FractionB)
	}
	if opts.Reciprocal && opts.FractionA <= 0 {
		return 0, fmt.Errorf("-r requires -f to be specified")
	}
	if opts.MinMAPQ < 0 || opts.MinMAPQ > 255 {
		return 0, fmt.Errorf("-q must be in [0,255], got %d", opts.MinMAPQ)
	}
	if opts.MaxDepth < 0 {
		return 0, fmt.Errorf("-D must be ≥0, got %d", opts.MaxDepth)
	}
	// Reciprocal: apply FractionA threshold to B as well.
	effFracB := opts.FractionB
	if opts.Reciprocal && effFracB < opts.FractionA {
		effFracB = opts.FractionA
	}

	trees := make([]map[string]*bedpkg.IntervalTree, len(srcs))
	splitMeta := make([][]splitAln, len(srcs)) // populated only when src is split-mode BAM
	kinds := make([]SourceKind, len(srcs))
	for i, s := range srcs {
		switch s.Kind {
		case SourceBED:
			t, err := indexBED(s.Reader)
			if err != nil {
				return 0, fmt.Errorf("file %d (BED): %w", i+1, err)
			}
			trees[i] = t
		case SourceBAM:
			if opts.Split {
				t, meta, err := indexBAMSplit(s.Reader, opts.MinMAPQ)
				if err != nil {
					return 0, fmt.Errorf("file %d (BAM): %w", i+1, err)
				}
				trees[i] = t
				splitMeta[i] = meta
			} else {
				t, err := indexBAM(s.Reader, opts.MinMAPQ)
				if err != nil {
					return 0, fmt.Errorf("file %d (BAM): %w", i+1, err)
				}
				trees[i] = t
			}
		default:
			return 0, fmt.Errorf("file %d: unknown source kind %d", i+1, s.Kind)
		}
		kinds[i] = s.Kind
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

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
		// Emit A's original columns verbatim, then one count per source.
		if _, err := bw.WriteString(strings.Join(fields, "\t")); err != nil {
			return count, err
		}
		for i, t := range trees {
			var n int
			if kinds[i] == SourceBAM && opts.Split {
				n = countOverlapsSplit(rec, t[rec.Chrom], splitMeta[i], opts)
			} else {
				n = countOverlaps(rec, t[rec.Chrom], opts, effFracB)
			}
			// MaxDepth caps the reported count for BAM inputs.
			if kinds[i] == SourceBAM && opts.MaxDepth > 0 && n > opts.MaxDepth {
				n = opts.MaxDepth
			}
			if _, err := fmt.Fprintf(bw, "\t%d", n); err != nil {
				return count, err
			}
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

// countOverlaps returns the number of B records that overlap a after
// applying strand and fraction filters.
func countOverlaps(a *bedpkg.Record, t *bedpkg.IntervalTree, opts Options, effFracB float64) int {
	if t == nil {
		return 0
	}
	cand := t.Query(a)
	if len(cand) == 0 {
		return 0
	}
	lenA := a.ChromEnd - a.ChromStart
	n := 0
	for _, b := range cand {
		if !strandOK(a, b, opts) {
			continue
		}
		overlap := overlapLen(a, b)
		if overlap <= 0 {
			continue
		}
		if opts.FractionA > 0 && lenA > 0 {
			if float64(overlap)/float64(lenA) < opts.FractionA {
				continue
			}
		}
		if effFracB > 0 {
			lenB := b.ChromEnd - b.ChromStart
			if lenB <= 0 {
				continue
			}
			if float64(overlap)/float64(lenB) < effFracB {
				continue
			}
		}
		n++
	}
	return n
}

// overlapLen returns the length of the intersection of a and b's spans.
// 0 if disjoint.
func overlapLen(a, b *bedpkg.Record) int {
	start := a.ChromStart
	if b.ChromStart > start {
		start = b.ChromStart
	}
	end := a.ChromEnd
	if b.ChromEnd < end {
		end = b.ChromEnd
	}
	if end <= start {
		return 0
	}
	return end - start
}

// strandOK applies the -s / -S filters. Missing strand on either side is
// treated as "no match" under a strand filter (matches upstream).
func strandOK(a, b *bedpkg.Record, opts Options) bool {
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

// parseRecord parses the minimum subset of a BED line we need for overlap
// + strand filtering. Extra columns are preserved by the caller as raw
// fields.
func parseRecord(fields []string) (*bedpkg.Record, error) {
	start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
	}
	end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
	}
	r := &bedpkg.Record{
		Chrom:      fields[0],
		ChromStart: start,
		ChromEnd:   end,
	}
	if len(fields) >= 6 {
		r.Strand = fields[5]
	}
	return r, nil
}

// indexBED reads a B file fully into memory and returns a per-chrom tree.
func indexBED(r io.Reader) (map[string]*bedpkg.IntervalTree, error) {
	rd := bedpkg.NewReader(r)
	all, err := rd.ReadAll()
	if err != nil {
		return nil, err
	}
	byChrom := map[string][]*bedpkg.Record{}
	for _, x := range all {
		byChrom[x.Chrom] = append(byChrom[x.Chrom], x)
	}
	return buildTrees(byChrom), nil
}

// indexBAM decodes every primary alignment from a BGZF-wrapped BAM stream,
// drops unmapped / secondary / supplementary / duplicate / QC-fail records
// (matching `bedtools multicov`'s default record filter), enforces the
// caller-supplied MAPQ floor, and returns a per-reference interval tree
// over the alignments' reference spans. The recorded interval is
// [Pos-1, Pos-1+ReferenceLength()) and the BAM strand is the FLAG-derived
// `-` / `+` so `-s` / `-S` keep working on BAM inputs.
//
// `-split` (per-block CIGAR coverage) is handled by indexBAMSplit; this
// function is only called when Options.Split is false. With Split off the
// recorded interval is the full reference footprint of each alignment,
// which is what upstream emits without `-split`.
func indexBAM(r io.Reader, minMAPQ int) (map[string]*bedpkg.IntervalTree, error) {
	br, err := sam.NewBAMReader(r)
	if err != nil {
		return nil, err
	}
	byChrom := map[string][]*bedpkg.Record{}
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
			rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		if int(rec.MapQ) < minMAPQ {
			continue
		}
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			continue
		}
		// BAM stores POS as 1-based; convert to 0-based half-open BED-style.
		start := int(rec.Pos) - 1
		if start < 0 {
			continue
		}
		end := start + refLen
		strand := "+"
		if rec.Flag&sam.FlagReverse != 0 {
			strand = "-"
		}
		b := &bedpkg.Record{
			Chrom:      rec.RName,
			ChromStart: start,
			ChromEnd:   end,
			Strand:     strand,
		}
		byChrom[b.Chrom] = append(byChrom[b.Chrom], b)
	}
	return buildTrees(byChrom), nil
}

// splitAln is per-alignment metadata for `-split` mode: the strand and
// the total reference footprint summed across all emitted blocks. Indexed
// by the value stored in `bed.Record.Score` of every block contributed by
// this alignment.
type splitAln struct {
	strand    string
	footprint int
}

// indexBAMSplit is the `-split` counterpart of indexBAM: it walks each
// alignment's CIGAR and emits one `bed.Record` per contiguous reference-
// consuming op-run (M/=/X). N ops break a block, advance the reference,
// and contribute zero coverage — the whole point of `-split`. D, I, S, H,
// P are handled to upstream's semantics:
//
//   - M, =, X advance both query and reference and extend the current block.
//   - D advances reference only and EXTENDS the current block (matches
//     upstream's `breakOnDeletionOps=false` for multicov).
//   - N advances reference only and BREAKS the current block (matches
//     upstream's `breakOnSkipOps=true`).
//   - I, S, P consume neither reference position — ignored for block math.
//   - H consumes neither — ignored.
//
// The per-block `bed.Record` is tagged with `Score = alignment index` so
// `countOverlapsSplit` can group hits back to their parent alignment for
// the once-per-alignment counting rule and the footprint-fraction filter.
// The returned `[]splitAln` is indexed by that same alignment index.
func indexBAMSplit(r io.Reader, minMAPQ int) (map[string]*bedpkg.IntervalTree, []splitAln, error) {
	br, err := sam.NewBAMReader(r)
	if err != nil {
		return nil, nil, err
	}
	byChrom := map[string][]*bedpkg.Record{}
	var alns []splitAln
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, err
		}
		if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
			rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		if int(rec.MapQ) < minMAPQ {
			continue
		}
		// BAM stores 1-based POS; bedtools internally uses 0-based.
		curPos := int(rec.Pos) - 1
		if curPos < 0 {
			continue
		}
		strand := "+"
		if rec.Flag&sam.FlagReverse != 0 {
			strand = "-"
		}

		// Walk the CIGAR, emitting one block per contiguous reference-
		// consuming run that is broken by an N op.
		blockLen := 0
		blockStart := curPos
		footprint := 0
		alnIdx := len(alns)
		// blocks emitted by this alignment, deferred so we don't append to
		// byChrom until we know the alignment actually produced ≥1 block.
		var pending []*bedpkg.Record
		emit := func() {
			if blockLen <= 0 {
				return
			}
			pending = append(pending, &bedpkg.Record{
				Chrom:      rec.RName,
				ChromStart: blockStart,
				ChromEnd:   blockStart + blockLen,
				Strand:     strand,
				Score:      alnIdx,
			})
			footprint += blockLen
		}
		for _, op := range rec.Cigar {
			switch op.Op() {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				if blockLen == 0 {
					blockStart = curPos
				}
				l := int(op.Length())
				blockLen += l
				curPos += l
			case sam.CigarDeletion:
				// Extend block (upstream: breakOnDeletionOps=false).
				if blockLen == 0 {
					blockStart = curPos
				}
				l := int(op.Length())
				blockLen += l
				curPos += l
			case sam.CigarSkipped: // N
				emit()
				curPos += int(op.Length())
				blockLen = 0
			case sam.CigarInsertion, sam.CigarSoftClip, sam.CigarHardClip, sam.CigarPadding:
				// No reference advance, no block contribution.
			}
		}
		emit()
		if len(pending) == 0 {
			continue
		}
		alns = append(alns, splitAln{strand: strand, footprint: footprint})
		for _, b := range pending {
			byChrom[b.Chrom] = append(byChrom[b.Chrom], b)
		}
	}
	return buildTrees(byChrom), alns, nil
}

// countOverlapsSplit counts the number of distinct BAM alignments whose
// emitted blocks overlap A, applying strand and fraction filters using
// upstream's `-split` semantics:
//
//   - The alignment counts once iff at least one block has positive
//     overlap with A.
//   - When FractionA > 0, the threshold is applied to
//     total_block_overlap / alignment.footprint with a strict `>`
//     comparison — exactly mirroring
//     `multiBamCov.cpp::FindBlockedOverlaps`. (Note: this differs from
//     the non-split path, which uses overlap/lenA ≥ frac. The asymmetry
//     is upstream's behaviour, not a port artifact.)
//   - Reciprocal additionally requires
//     total_block_overlap / lenA > FractionA.
//   - SameStrand / OppositeStrand filter on the alignment's strand vs A's.
func countOverlapsSplit(a *bedpkg.Record, t *bedpkg.IntervalTree, alns []splitAln, opts Options) int {
	if t == nil {
		return 0
	}
	cand := t.Query(a)
	if len(cand) == 0 {
		return 0
	}
	// Aggregate overlap per alignment index.
	perAln := make(map[int]int, len(cand))
	for _, b := range cand {
		ov := overlapLen(a, b)
		if ov <= 0 {
			continue
		}
		perAln[b.Score] += ov
	}
	lenA := a.ChromEnd - a.ChromStart
	n := 0
	for idx, totOverlap := range perAln {
		if idx < 0 || idx >= len(alns) {
			continue
		}
		al := alns[idx]
		// Strand filter on the parent alignment.
		if opts.SameStrand {
			if a.Strand == "" || al.strand == "" || a.Strand != al.strand {
				continue
			}
		}
		if opts.OppositeStrand {
			if a.Strand == "" || al.strand == "" || a.Strand == al.strand {
				continue
			}
		}
		if opts.FractionA > 0 {
			if al.footprint <= 0 {
				continue
			}
			// Upstream uses strict `>` here, not `>=`.
			if float64(totOverlap)/float64(al.footprint) <= opts.FractionA {
				continue
			}
			if opts.Reciprocal {
				if lenA <= 0 {
					continue
				}
				if float64(totOverlap)/float64(lenA) <= opts.FractionA {
					continue
				}
			}
		}
		n++
	}
	return n
}

// buildTrees turns a per-chrom record map into per-chrom interval trees,
// sorted by (start,end) so the tree is balanced.
func buildTrees(byChrom map[string][]*bedpkg.Record) map[string]*bedpkg.IntervalTree {
	out := make(map[string]*bedpkg.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].ChromStart != recs[j].ChromStart {
				return recs[i].ChromStart < recs[j].ChromStart
			}
			return recs[i].ChromEnd < recs[j].ChromEnd
		})
		out[chrom] = bedpkg.NewIntervalTree(recs)
	}
	return out
}
