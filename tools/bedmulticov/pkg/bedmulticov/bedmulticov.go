// Package bedmulticov implements `bedtools multicov`: for each interval in
// a primary BED file (`-bed`), it reports the count of overlapping records
// from each of N input files (`-files` / `-bams`). The output is the
// original A columns followed by N integer columns — one per input file —
// holding the overlap count.
//
// Upstream supports BED, BAM, *and* CRAM inputs. This port supports all
// three. BAM and CRAM inputs are decoded via pkg/htsgo/alnio (which routes
// CRAM to pkg/htsgo/cram and BAM to pkg/htsgo/sam); the MAPQ filter (`-q`)
// and per-interval depth cap (`-D`) are honoured for both. A CRAM input may
// supply a decode reference via Options.Reference (`-T` / the REF_CACHE
// environment variable), though multicov only reads alignment coordinates
// and so decodes CRAM correctly even without a reference.
//
// Internally the *A* file (`-bed`, typically small) is loaded into a
// per-chromosome interval tree (`pkg/htsgo/bed.IntervalTree`) and each B
// input is then streamed a record at a time and discarded after use — the
// data flow is inverted relative to the naive approach so peak memory is
// O(A intervals) plus one decoded read, not O(B alignments). Each streamed
// B record queries the A tree and increments a per-A counter for every A
// interval it overlaps (subject to the filters below); the A rows are then
// emitted verbatim in their original input order with one count column per
// source appended. Optional strand filters (-s same / -S opposite),
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

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	bedpkg "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
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
	// Reference is the path to a FASTA reference used to decode
	// reference-backed CRAM inputs (the analogue of `samtools view -T` /
	// upstream htslib's `CRAM_OPT_REFERENCE`). It is ignored for BED and
	// BAM inputs, which carry their coordinates inline. multicov only needs
	// each alignment's reference span (POS, CIGAR, FLAG, MAPQ), none of
	// which require base reconstruction, so a CRAM input decodes correctly
	// even when Reference is empty; the option exists for completeness. The
	// REF_CACHE environment variable is also honoured.
	Reference string
}

// SourceKind tags an input as BED, BAM, or CRAM.
type SourceKind int

const (
	// SourceBED is a plain BED file (any number of columns ≥3).
	SourceBED SourceKind = iota
	// SourceBAM is a BGZF-wrapped BAM file (decoded via
	// pkg/htsgo/sam.NewBAMReader). Each primary alignment contributes
	// one interval over [Pos-1, Pos-1+ReferenceLength()) on its reference.
	SourceBAM
	// SourceCRAM is a CRAM file (decoded via pkg/htsgo/alnio, which routes
	// to pkg/htsgo/cram). It is handled identically to SourceBAM once
	// decoded: each primary alignment contributes one reference-span
	// interval. Options.Reference supplies the optional decode reference.
	SourceCRAM
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

// RunSources reads A from aR and the N B inputs from srcs in order. The
// data flow is inverted for memory: the (typically small) A file is loaded
// into a per-chromosome interval tree while the B inputs are streamed one
// record at a time and discarded after use (BED records are read with the
// bed package; BAM/CRAM records are decoded with pkg/htsgo/alnio and -q
// MAPQ filtered on the fly). Each streamed B record queries the A tree and
// increments a per-A counter for every A interval it overlaps, so peak
// memory is O(A intervals) plus one decoded read rather than O(B records).
// RunSources then emits one row per A record — in original input order —
// with one count column appended per source. -D, if set, caps the reported
// per-A count per BAM/CRAM input at emit time. Returns the number of A
// records processed.
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

	// Load A into an ordered slice (for verbatim emit) and a per-chrom
	// interval tree. Each A record's Score field carries its 0-based index
	// so a tree hit can be mapped back to its counters.
	aRaw, aRecs, err := readARecords(aR)
	if err != nil {
		return 0, err
	}
	aTrees := buildATrees(aRecs)

	// counts[aIdx][srcIdx] accumulates overlaps for A record aIdx from
	// source srcIdx.
	counts := make([][]int, len(aRecs))
	for i := range counts {
		counts[i] = make([]int, len(srcs))
	}
	kinds := make([]SourceKind, len(srcs))

	// Stream each source once, discarding each record after it has been
	// scored against the A tree.
	for i, s := range srcs {
		kinds[i] = s.Kind
		switch s.Kind {
		case SourceBED:
			if err := streamBED(s.Reader, aTrees, counts, i, opts, effFracB); err != nil {
				return 0, fmt.Errorf("file %d (BED): %w", i+1, err)
			}
		case SourceBAM, SourceCRAM:
			label := "BAM"
			if s.Kind == SourceCRAM {
				label = "CRAM"
			}
			ar, err := alnio.NewReaderWithReference(s.Reader, opts.Reference)
			if err != nil {
				return 0, fmt.Errorf("file %d (%s): %w", i+1, label, err)
			}
			if opts.Split {
				if err := streamAlignmentsSplit(ar, aTrees, counts, i, opts); err != nil {
					return 0, fmt.Errorf("file %d (%s): %w", i+1, label, err)
				}
			} else {
				if err := streamAlignments(ar, aTrees, counts, i, opts, effFracB); err != nil {
					return 0, fmt.Errorf("file %d (%s): %w", i+1, label, err)
				}
			}
		default:
			return 0, fmt.Errorf("file %d: unknown source kind %d", i+1, s.Kind)
		}
	}

	// Emit A rows verbatim in original order, appending the per-source count.
	bw := bufio.NewWriter(out)
	defer bw.Flush()
	for aIdx, raw := range aRaw {
		if _, err := bw.WriteString(raw); err != nil {
			return aIdx, err
		}
		for srcIdx := range srcs {
			n := counts[aIdx][srcIdx]
			// MaxDepth caps the reported count for BAM/CRAM inputs.
			if isAlignment(kinds[srcIdx]) && opts.MaxDepth > 0 && n > opts.MaxDepth {
				n = opts.MaxDepth
			}
			if _, err := fmt.Fprintf(bw, "\t%d", n); err != nil {
				return aIdx, err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return aIdx, err
		}
	}
	return len(aRecs), nil
}

// readARecords reads the primary (-bed) A file fully into memory: it returns
// the raw output columns (tab-joined) in input order and the parsed records
// used for overlap/strand tests. Comment/track/browser lines are skipped and
// each parsed record's Score field is set to its 0-based index so a later
// interval-tree hit can be mapped back to its output row. The A file is the
// small input; keeping it resident is what makes the streaming B pass O(A).
func readARecords(aR io.Reader) ([]string, []*bedpkg.Record, error) {
	sc := bufio.NewScanner(aR)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var raws []string
	var recs []*bedpkg.Record
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
			return nil, nil, fmt.Errorf("line %d: BED record needs >=3 columns: %q", lineNo, raw)
		}
		rec, err := parseRecord(fields)
		if err != nil {
			return nil, nil, fmt.Errorf("line %d: %w", lineNo, err)
		}
		rec.Score = len(recs)
		recs = append(recs, rec)
		raws = append(raws, strings.Join(fields, "\t"))
	}
	if err := sc.Err(); err != nil {
		return nil, nil, err
	}
	return raws, recs, nil
}

// buildATrees groups the A records by chromosome and builds one balanced
// interval tree per chromosome. Records keep their Score-encoded index
// across the sort so tree hits still map back to the right output row.
func buildATrees(recs []*bedpkg.Record) map[string]*bedpkg.IntervalTree {
	byChrom := map[string][]*bedpkg.Record{}
	for _, r := range recs {
		byChrom[r.Chrom] = append(byChrom[r.Chrom], r)
	}
	return buildTrees(byChrom)
}

// streamBED reads a B BED file one record at a time and scores each against
// the A tree. Records are not retained, so peak memory stays O(A).
func streamBED(r io.Reader, aTrees map[string]*bedpkg.IntervalTree, counts [][]int, srcIdx int, opts Options, effFracB float64) error {
	rd := bedpkg.NewReader(r)
	for {
		b, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		countBIntoA(b, aTrees[b.Chrom], counts, srcIdx, opts, effFracB)
	}
	return nil
}

// streamAlignments is the BAM/CRAM counterpart of streamBED for the
// non-split path. It decodes each alignment, turns its reference span into a
// transient B record, scores it against the A tree, and discards it.
// bedtools multicov counts reads "regardless of the BAM FLAG field", so only
// duplicate and QC-fail records are dropped (plus the MAPQ floor); unmapped,
// secondary, and supplementary reads are all counted. An unmapped read (or one
// with no CIGAR) is given a 1bp span [POS-1, POS) — bam_endpos semantics — so
// it still overlaps its enclosing interval. Each alignment contributes at most
// +1 to every A interval it overlaps (subject to the strand/fraction filters).
func streamAlignments(br sam.Reader, aTrees map[string]*bedpkg.IntervalTree, counts [][]int, srcIdx int, opts Options, effFracB float64) error {
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		// bedtools multicov counts reads "regardless of the BAM FLAG field":
		// it drops only duplicate / QC-fail records (plus the MAPQ floor and
		// the strand/properOnly filters). Unmapped, secondary, and
		// supplementary reads are NOT dropped — matching multiBamCov's default.
		if rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		if int(rec.MapQ) < opts.MinMAPQ {
			continue
		}
		// BAM stores POS as 1-based; convert to 0-based half-open BED-style.
		start := int(rec.Pos) - 1
		if start < 0 {
			continue
		}
		// bam_endpos semantics: a mapped read spans [start, start+refLen),
		// but an unmapped read (or one with no CIGAR) is given a 1bp span
		// [start, start+1) so it still overlaps its enclosing interval — this
		// is what upstream counts and what our old code wrongly dropped.
		refLen := rec.Cigar.ReferenceLength()
		end := start + refLen
		if rec.IsUnmapped() || len(rec.Cigar) == 0 {
			end = start + 1
		}
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
		countBIntoA(b, aTrees[rec.RName], counts, srcIdx, opts, effFracB)
	}
	return nil
}

// countBIntoA scores a single B record against the A interval tree,
// incrementing counts[a.Score][srcIdx] for every A interval that survives
// the strand, fraction-of-A (-f), and fraction-of-B (-F/-r) filters. It is
// the inverted-flow analogue of the old countOverlaps: the same (A,B) pairs
// are evaluated against the same predicate, only the iteration is driven by
// B instead of A, so the resulting counts are identical.
func countBIntoA(b *bedpkg.Record, t *bedpkg.IntervalTree, counts [][]int, srcIdx int, opts Options, effFracB float64) {
	if t == nil {
		return
	}
	cand := t.Query(b)
	if len(cand) == 0 {
		return
	}
	lenB := b.ChromEnd - b.ChromStart
	for _, a := range cand {
		if !strandOK(a, b, opts) {
			continue
		}
		overlap := overlapLen(a, b)
		if overlap <= 0 {
			continue
		}
		if opts.FractionA > 0 {
			lenA := a.ChromEnd - a.ChromStart
			if lenA > 0 {
				if float64(overlap)/float64(lenA) < opts.FractionA {
					continue
				}
			}
		}
		if effFracB > 0 {
			if lenB <= 0 {
				continue
			}
			if float64(overlap)/float64(lenB) < effFracB {
				continue
			}
		}
		counts[a.Score][srcIdx]++
	}
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

// isAlignment reports whether a SourceKind is a SAM-family alignment input
// (BAM or CRAM), as opposed to a plain BED file. Alignment inputs share the
// `-split`, `-q`, and `-D` handling that BED inputs ignore.
func isAlignment(k SourceKind) bool {
	return k == SourceBAM || k == SourceCRAM
}

// streamAlignmentsSplit is the `-split` counterpart of streamAlignments. It
// folds the per-alignment CIGAR decomposition into the streaming loop: each
// primary alignment is decomposed locally into its contiguous reference-
// consuming op-runs (M/=/X), scored against the A tree, and discarded — no
// per-block tree over the whole BAM is ever built. D, I, S, H, P are handled
// to upstream's semantics:
//
//   - M, =, X advance both query and reference and extend the current block.
//   - D advances reference only and EXTENDS the current block (matches
//     upstream's `breakOnDeletionOps=false` for multicov).
//   - N advances reference only and BREAKS the current block (matches
//     upstream's `breakOnSkipOps=true`).
//   - I, S, P consume neither reference position — ignored for block math.
//   - H consumes neither — ignored.
//
// Counting mirrors the old countOverlapsSplit exactly, only inverted so the
// aggregation is keyed by A index instead of alignment index:
//
//   - The alignment counts once per A iff at least one of its blocks has
//     positive overlap with that A.
//   - When FractionA > 0, the threshold is applied to
//     total_block_overlap / footprint with a strict `>` comparison — exactly
//     mirroring `multiBamCov.cpp::FindBlockedOverlaps`. (This differs from
//     the non-split path, which uses overlap/lenA ≥ frac; the asymmetry is
//     upstream's behaviour, not a port artifact.)
//   - Reciprocal additionally requires total_block_overlap / lenA > FractionA.
//   - SameStrand / OppositeStrand filter on the alignment's strand vs A's.
//
// footprint is the sum of every positive-length block of the alignment,
// independent of which A intervals (if any) it overlaps.
func streamAlignmentsSplit(br sam.Reader, aTrees map[string]*bedpkg.IntervalTree, counts [][]int, srcIdx int, opts Options) error {
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		// Same FLAG policy as the non-split path: drop only duplicate /
		// QC-fail (plus the MAPQ floor). Unmapped / no-CIGAR reads emit no
		// blocks below and so contribute 0 naturally; secondary and
		// supplementary reads are counted consistently with the non-split path.
		if rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		if int(rec.MapQ) < opts.MinMAPQ {
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
		// consuming run that is broken by an N op. Blocks are transient and
		// local to this alignment; they are discarded once scored.
		blockLen := 0
		blockStart := curPos
		footprint := 0
		var blocks []*bedpkg.Record
		emit := func() {
			if blockLen <= 0 {
				return
			}
			blocks = append(blocks, &bedpkg.Record{
				ChromStart: blockStart,
				ChromEnd:   blockStart + blockLen,
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
		if len(blocks) == 0 {
			continue
		}
		t := aTrees[rec.RName]
		if t == nil {
			continue
		}

		// Aggregate this alignment's total overlap per A interval (keyed by
		// the A record's Score-encoded index) across all of its blocks.
		perA := make(map[int]int)
		aByIdx := make(map[int]*bedpkg.Record)
		for _, blk := range blocks {
			for _, a := range t.Query(blk) {
				ov := overlapLen(a, blk)
				if ov <= 0 {
					continue
				}
				perA[a.Score] += ov
				aByIdx[a.Score] = a
			}
		}
		for idx, totOverlap := range perA {
			a := aByIdx[idx]
			// Strand filter on the parent alignment.
			if opts.SameStrand {
				if a.Strand == "" || strand == "" || a.Strand != strand {
					continue
				}
			}
			if opts.OppositeStrand {
				if a.Strand == "" || strand == "" || a.Strand == strand {
					continue
				}
			}
			if opts.FractionA > 0 {
				if footprint <= 0 {
					continue
				}
				// Upstream uses strict `>` here, not `>=`.
				if float64(totOverlap)/float64(footprint) <= opts.FractionA {
					continue
				}
				if opts.Reciprocal {
					lenA := a.ChromEnd - a.ChromStart
					if lenA <= 0 {
						continue
					}
					if float64(totOverlap)/float64(lenA) <= opts.FractionA {
						continue
					}
				}
			}
			counts[idx][srcIdx]++
		}
	}
	return nil
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
