// Package samtools — phase subcommand.
//
// Phase walks a coordinate-sorted SAM/BAM stream against an in-memory
// pileup, identifies heterozygous SNP positions, assigns each read at
// each het to one of two haplotype clusters based on its base, then
// chains those clusters across overlapping reads to produce phased
// blocks.
//
// Upstream samtools (`reference_code/samtools/phase.c`) drives the
// chaining with a Markov-chain-Monte-Carlo solver; the v1 Go port
// implements only the common-case greedy chaining. For each pair of
// adjacent het sites we count the number of reads that span both. If
// the same-allele count outweighs the opposite-allele count we keep
// the labels aligned; if opposite outweighs same we flip the labels;
// if neither dominates we emit `0` (ambiguous) for the current het.
// The MCMC fallback that resolves chimeras and tied junctions is
// deliberately deferred — see docs/PARITY_ROADMAP.md.
//
// Output format is the tab-separated stream documented in the user
// spec:
//
//	PS<TAB>chrom<TAB>pos<TAB>{0|1|2}
//
// where 0 = ambiguous (no consistent cluster), 1 = hap1, 2 = hap2.
// One line per het SNP, in coordinate order. Het positions are
// 1-based to match SAM POS.
//
// A new phase block is implicitly started whenever the distance (in
// number of intervening het sites) to the previous successfully-phased
// het exceeds the block-merge window `k`. The CLI default for `k` is
// 13, matching upstream.
package samtools

import (
	"bufio"
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// PhaseOptions configures the Phase walker.
type PhaseOptions struct {
	// BlockWindow is the maximum number of intervening (non-phased)
	// het sites between two successfully-phased sites before a new
	// phase block starts. Matches upstream samtools phase -k.
	BlockWindow int
	// MinMAPQ drops records whose MAPQ is strictly less than this.
	MinMAPQ uint8
	// MinBaseQ drops query bases whose Phred quality is strictly less
	// than this.
	MinBaseQ uint8
	// MaxDepth caps the number of reads observed at any one het. The
	// upstream default is 256.
	MaxDepth int
	// FullRead, when set, mirrors upstream's -F flag: use the full
	// read regardless of soft-clipped bases. v1 always uses the
	// aligned slice; the flag is accepted on the CLI but its current
	// effect is a no-op (see PARITY_ROADMAP.md).
	FullRead bool
	// DropAmbiguous, when set, mirrors upstream's -A flag: indicate
	// "drop" in the chimera/dropped output. v1 has no per-read
	// chimera output, so this flag is accepted but informational.
	DropAmbiguous bool
	// OutputPrefix is upstream's -b STR option. v1 emits the phased-
	// TSV stream to the writer the caller passes; this option is
	// accepted on the CLI but not yet wired through to per-haplotype
	// BAM splitting (which is the upstream behaviour).
	OutputPrefix string
}

// Phase default constants matching upstream samtools phase.c.
const (
	DefaultPhaseBlockWindow  = 13
	DefaultPhaseMinMAPQ      = 13
	DefaultPhaseMinBaseQ     = 13
	DefaultPhaseMaxDepth     = 256
	DefaultPhaseOutputPrefix = ""
)

// Phase reads SAM/BAM records from in, identifies het SNPs, and writes
// the phased-position TSV to out. Returns the number of het sites
// emitted and the first error encountered.
func Phase(in io.Reader, out io.Writer, opts PhaseOptions) (int, error) {
	if opts.BlockWindow == 0 {
		opts.BlockWindow = DefaultPhaseBlockWindow
	}
	if opts.MaxDepth == 0 {
		opts.MaxDepth = DefaultPhaseMaxDepth
	}
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, fmt.Errorf("samtools phase: open input: %w", err)
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	// We need a pileup. The simplest correct approach for v1 is to
	// read every record, group them by reference, then for each
	// reference scan the read set and build per-position base
	// observations. This is O(reference-length × depth) memory but
	// works for the small-to-medium BAMs the original phase tool was
	// designed for (upstream caps depth at 256 too).
	byRef := make(map[string][]*sam.Record)
	refOrder := []string{}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("samtools phase: read: %w", err)
		}
		if rec.IsUnmapped() {
			continue
		}
		if rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
			continue
		}
		if rec.MapQ < opts.MinMAPQ {
			continue
		}
		if rec.RName == "" || rec.RName == "*" {
			continue
		}
		if _, seen := byRef[rec.RName]; !seen {
			refOrder = append(refOrder, rec.RName)
		}
		byRef[rec.RName] = append(byRef[rec.RName], rec)
	}

	emitted := 0
	for _, ref := range refOrder {
		recs := byRef[ref]
		sort.SliceStable(recs, func(i, j int) bool { return recs[i].Pos < recs[j].Pos })
		hets, err := callHets(recs, opts)
		if err != nil {
			return emitted, err
		}
		phased := phaseHets(hets, opts)
		for _, h := range phased {
			line := fmt.Sprintf("PS\t%s\t%d\t%d\n", ref, h.pos, h.label)
			if _, err := bw.WriteString(line); err != nil {
				return emitted, err
			}
			emitted++
		}
	}
	return emitted, nil
}

// het is one heterozygous SNP candidate. allele0 / allele1 are the
// two most-supported bases at this position; reads listed by index
// in recs along with the allele each supports (0 or 1).
type het struct {
	pos     int32
	allele0 byte
	allele1 byte
	support []readSupport
}

// readSupport pairs a record index (into the per-reference slice
// passed to callHets) with the allele assignment at one het site
// (0 = allele0, 1 = allele1).
type readSupport struct {
	readIdx int
	allele  int // 0 or 1
}

// callHets scans the per-reference records and returns one het entry
// per position where:
//   - at least two distinct query bases were observed, AND
//   - the two most-common bases each have ≥ 2 supporting reads (a
//     minimal "het call"), AND
//   - the supporting reads' query quality at that position passed the
//     opts.MinBaseQ filter.
//
// Het positions are 1-based (matching SAM POS) and returned in
// coordinate order.
func callHets(recs []*sam.Record, opts PhaseOptions) ([]het, error) {
	// per-position: map base → list of (readIdx).
	type bucket struct {
		// bases[b] holds read indices that called base b at this
		// reference position. Index by 'A'/'C'/'G'/'T' lowered to
		// 0..3 via baseIdx.
		bases [4][]int
	}
	buckets := make(map[int32]*bucket)

	for i, rec := range recs {
		walkAlignment(rec, func(refPos int32, queryBase byte, queryQ uint8) {
			if queryQ < opts.MinBaseQ {
				return
			}
			bi, ok := baseIdx(queryBase)
			if !ok {
				return
			}
			b := buckets[refPos]
			if b == nil {
				b = &bucket{}
				buckets[refPos] = b
			}
			if len(b.bases[0])+len(b.bases[1])+len(b.bases[2])+len(b.bases[3]) >= opts.MaxDepth {
				return
			}
			b.bases[bi] = append(b.bases[bi], i)
		})
	}

	// Collect het positions in order.
	positions := make([]int32, 0, len(buckets))
	for p := range buckets {
		positions = append(positions, p)
	}
	sort.Slice(positions, func(i, j int) bool { return positions[i] < positions[j] })

	out := make([]het, 0, len(positions))
	for _, p := range positions {
		b := buckets[p]
		// Find top two alleles.
		type ac struct {
			base  byte
			count int
			reads []int
		}
		var alleles []ac
		for bi := 0; bi < 4; bi++ {
			if len(b.bases[bi]) > 0 {
				alleles = append(alleles, ac{base: "ACGT"[bi], count: len(b.bases[bi]), reads: b.bases[bi]})
			}
		}
		if len(alleles) < 2 {
			continue
		}
		sort.SliceStable(alleles, func(i, j int) bool { return alleles[i].count > alleles[j].count })
		if alleles[0].count < 2 || alleles[1].count < 2 {
			continue
		}
		h := het{pos: p, allele0: alleles[0].base, allele1: alleles[1].base}
		// Only reads that pick allele0 or allele1 contribute to phasing;
		// reads picking a third / fourth allele are ignored for this site.
		for _, ri := range alleles[0].reads {
			h.support = append(h.support, readSupport{readIdx: ri, allele: 0})
		}
		for _, ri := range alleles[1].reads {
			h.support = append(h.support, readSupport{readIdx: ri, allele: 1})
		}
		sort.SliceStable(h.support, func(i, j int) bool { return h.support[i].readIdx < h.support[j].readIdx })
		out = append(out, h)
	}
	return out, nil
}

// phasedSite is one output row: pos and the {0,1,2} label.
type phasedSite struct {
	pos   int32
	label int // 0 ambig, 1 hap1, 2 hap2
}

// phaseHets walks the het list in coordinate order and produces one
// phasedSite per input het. The first het that picks up any support
// is labelled with allele0 -> hap1, allele1 -> hap2. Each subsequent
// het is compared with the previous successfully-phased site: we
// count the read overlaps that agree with the current labelling vs.
// the flipped labelling. The winner sets the label; a tie emits
// label 0 (ambiguous). Blocks reset when more than opts.BlockWindow
// hets in a row are unphased.
func phaseHets(hets []het, opts PhaseOptions) []phasedSite {
	out := make([]phasedSite, 0, len(hets))
	// Per-read: which allele did it pick at the previous phased het?
	// readAssign[readIdx] = 0 / 1 / -1 (no assignment).
	readAssign := map[int]int{}
	prevPhasedIdx := -1
	consecUnphased := 0

	for i, h := range hets {
		if prevPhasedIdx < 0 {
			// Block start: label the first het arbitrarily.
			out = append(out, phasedSite{pos: h.pos, label: 1})
			for _, s := range h.support {
				readAssign[s.readIdx] = s.allele
			}
			prevPhasedIdx = i
			consecUnphased = 0
			continue
		}
		// Count overlap between this het's supporting reads and the
		// previous-het assignment.
		same, opposite := 0, 0
		for _, s := range h.support {
			prevAllele, ok := readAssign[s.readIdx]
			if !ok {
				continue
			}
			if prevAllele == s.allele {
				same++
			} else {
				opposite++
			}
		}
		var label int
		switch {
		case same == 0 && opposite == 0:
			label = 0 // no overlap → ambiguous
		case same == opposite:
			label = 0 // tied → ambiguous
		case same > opposite:
			label = 1 // current allele0 aligns with hap1
		default:
			label = 2 // current allele0 aligns with hap2 (labels flipped)
		}
		out = append(out, phasedSite{pos: h.pos, label: label})
		if label == 0 {
			consecUnphased++
			if consecUnphased > opts.BlockWindow {
				// Block break — reset.
				prevPhasedIdx = -1
				readAssign = map[int]int{}
				consecUnphased = 0
			}
			continue
		}
		// Update readAssign so future hets chain off this one. If the
		// label flipped (label==2) we record the OPPOSITE allele as
		// the "hap1 marker" so the next het's comparison stays
		// consistent.
		flip := label == 2
		for _, s := range h.support {
			a := s.allele
			if flip {
				a ^= 1
			}
			readAssign[s.readIdx] = a
		}
		prevPhasedIdx = i
		consecUnphased = 0
	}
	return out
}

// baseIdx maps an upper- or lower-case ACGT base to a 0..3 index.
func baseIdx(b byte) (int, bool) {
	switch b {
	case 'A', 'a':
		return 0, true
	case 'C', 'c':
		return 1, true
	case 'G', 'g':
		return 2, true
	case 'T', 't':
		return 3, true
	}
	return 0, false
}

// walkAlignment iterates the aligned bases of rec, invoking fn for
// each (1-based refPos, queryBase, queryQual) tuple where the CIGAR
// op consumes both reference and query (M/=/X). Insertions, deletions,
// soft- and hard-clips, refskips and padding are walked but never
// produce a callback — they are not phaseable positions.
func walkAlignment(rec *sam.Record, fn func(refPos int32, queryBase byte, queryQ uint8)) {
	if rec.Pos <= 0 || len(rec.Cigar) == 0 || rec.Seq == "" || rec.Seq == "*" {
		return
	}
	refPos := rec.Pos // 1-based
	qpos := 0
	hasQual := len(rec.Qual) == len(rec.Seq)
	for _, op := range rec.Cigar {
		o := op.Op()
		n := int(op.Length())
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for k := 0; k < n; k++ {
				if qpos+k >= len(rec.Seq) {
					break
				}
				var q uint8 = 255
				if hasQual {
					q = rec.Qual[qpos+k]
				}
				fn(refPos+int32(k), rec.Seq[qpos+k], q)
			}
			refPos += int32(n)
			qpos += n
		case sam.CigarInsertion, sam.CigarSoftClip:
			qpos += n
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += int32(n)
		case sam.CigarHardClip, sam.CigarPadding:
			// no movement on either axis
		}
	}
}
