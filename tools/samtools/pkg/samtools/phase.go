// Package samtools — phase subcommand.
//
// Phase walks a coordinate-sorted SAM/BAM stream against an in-memory
// pileup, identifies heterozygous SNP positions, assigns each read at
// each het to one of two haplotype clusters based on its base, then
// chains those clusters across overlapping reads to produce phased
// blocks.
//
// The CLI default path (PhaseOptions.UpstreamSchema, set by the
// `samtools phase` command) is a faithful, fully DETERMINISTIC port of
// upstream phase.c — there is no MCMC anywhere in upstream. The
// per-site haplotype assignment comes from `dynaprog` (phase_algo.go),
// a Viterbi-style dynamic program over k-bit local-haplotype states
// that fills an int8 `path[]`; `fragphase` then assigns each fragment
// to a haplotype vs. `path` and, when chimera repair is enabled (the
// default; upstream's `FLAG_FIX_CHIMERA`, cleared by
// `-F`/`--no-fix-chimera`), finds the best per-read flip point via a
// forward/backward sum scan (FLIP_PENALTY/FLIP_THRES). `genmask` emits
// the FL masked regions. The whole text stream (PS/FL/M/EV) is
// byte-identical to upstream; see phase_emit.go and the live oracles.
//
// A separate LEGACY path (UpstreamSchema=false) drives a simpler greedy
// adjacent-het chainer (phaseHets, below) that emits a v1 PS-label TSV.
// For each pair of adjacent het sites it counts the reads spanning both
// and keeps/flips/ambiguates the label by majority vote. It is NOT the
// CLI emit; it survives only for the in-process v1 unit tests.
//
// The only random source in upstream phase is the `drand48()` that
// routes evidence-less reads (and applies the per-call is_flip shuffle)
// into the 0/1/chimera output buckets in `-b` mode — phasing itself is
// RNG-free. On the upstream-schema path that RNG is an in-tree
// byte-exact port of glibc's default-seeded `drand48()`
// (phase_drand48.go), so the `-b` split agrees with upstream
// record-for-record — not merely up to a 0<->1 relabelling. The legacy
// v1 path retains a `math/rand` source seeded by RNGSeed.
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
// When the caller supplies `OutputPrefix` (upstream's `-b STR`), three
// BAM files are written alongside the TSV stream:
//
//	<prefix>.0.bam       — reads assigned to haplotype 0
//	<prefix>.1.bam       — reads assigned to haplotype 1
//	<prefix>.chimera.bam — reads that span both haplotypes (chimeric)
//
// Reads with no allele evidence are routed randomly to .0 or .1.
// When `DropAmbiguous` is set (upstream's `-A`), reads with weak but
// non-zero evidence on both haplotypes are routed to the chimera
// bucket instead of being kept in their majority bucket.
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
	"math/rand"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
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
	// than this. Mirrors upstream samtools phase -Q/--min-BQ.
	MinBaseQ uint8
	// MinVarLOD is the minimum heterozygous Phred-scaled LOD a pileup
	// column must reach to be treated as a variant (het) site. Mirrors
	// upstream samtools phase -q (g.min_varLOD; default 37). Columns
	// whose gl2cns LOD is below this are skipped — they emit neither an
	// M line nor a phase block. A NEGATIVE value selects the upstream
	// default (37) so existing library callers keep upstream behaviour;
	// a zero or positive value is used verbatim (e.g. -q 0 accepts any
	// LOD, exactly like upstream).
	MinVarLOD int
	// MaxDepth caps the number of reads observed at any one het. The
	// upstream default is 256.
	MaxDepth int
	// FullRead, when set, mirrors upstream's -F flag: when CLEARED it
	// disables the deterministic chimera repair (upstream's
	// FLAG_FIX_CHIMERA — there is no MCMC). The
	// historical name FullRead is retained for backwards compatibility
	// with the v1 API; semantically the field now means "skip chimera
	// repair when true". See NoFixChimera for the canonical accessor.
	//
	// Deprecated: use NoFixChimera in new code.
	FullRead bool
	// NoFixChimera, when set, disables the chimera-repair pass. This
	// mirrors upstream `-F`/`--no-fix-chimera`. When unset (the
	// default) the chimera-repair pass runs after greedy phasing.
	NoFixChimera bool
	// DropAmbiguous, when set, mirrors upstream's -A flag: ambiguous
	// reads are routed to the chimera bucket in the `-b` BAM split
	// instead of being kept in their majority bucket.
	DropAmbiguous bool
	// OutputPrefix is upstream's -b STR option. When non-empty, three
	// BAM files are written: `<prefix>.0.bam`, `<prefix>.1.bam`, and
	// `<prefix>.chimera.bam`.
	OutputPrefix string
	// RNGSeed seeds the random number generator used to route
	// truly-ambiguous reads to the 0/1 buckets in `-b` mode. When
	// zero, the default seed (1) is used; tests should pin this for
	// determinism.
	RNGSeed int64
	// UpstreamSchema selects the byte-faithful upstream samtools phase
	// emit pipeline (CC banner + PS / FL / M / EV / //). When false,
	// the legacy v1 PS-label TSV is emitted instead — that path is
	// preserved for the in-process unit tests that exercise the greedy
	// chainer. The CLI sets this to true by default; library callers
	// must opt in.
	UpstreamSchema bool
	// SiteListPath is the path to a phase site list (CLI `-l FILE`).
	// Each line is `CHROM\tPOS` (1-based). Sites in this set are
	// always treated as candidate hets regardless of the LOD test,
	// matching phase.c::loadpos.
	SiteListPath string
	// ListExclusive mirrors upstream `-e`. When set together with
	// SiteListPath, ONLY positions in the list are phased; sites
	// outside the list are dropped even when their LOD exceeds the
	// threshold. Matches phase.c FLAG_LIST_EXCL.
	ListExclusive bool
}

// Phase default constants matching upstream samtools phase.c.
const (
	DefaultPhaseBlockWindow  = 13
	DefaultPhaseMinMAPQ      = 13
	DefaultPhaseMinBaseQ     = 13
	DefaultPhaseMinVarLOD    = 37
	DefaultPhaseMaxDepth     = 256
	DefaultPhaseOutputPrefix = ""

	// flipPenalty and flipThreshold mirror upstream phase.c's
	// FLIP_PENALTY and FLIP_THRES — the scoring constants used by
	// the chimera-repair flip-point search in `fragphase`.
	flipPenalty   = 2
	flipThreshold = 4
)

// loadPhaseSites parses a phase site list file of the upstream
// loadpos shape: each non-comment line is "CHROM<sep>POS" with
// 1-based POS. Whitespace separators are tolerated. Returned map is
// indexed [chrom][pos1].
func loadPhaseSites(path string) (map[string]map[int]struct{}, error) {
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	out := map[string]map[int]struct{}{}
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || line[0] == '#' {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}
		pos, perr := strconv.Atoi(fields[1])
		if perr != nil || pos <= 0 {
			continue
		}
		if _, ok := out[fields[0]]; !ok {
			out[fields[0]] = map[int]struct{}{}
		}
		out[fields[0]][pos] = struct{}{}
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// Phase reads SAM/BAM records from in, identifies het SNPs, and writes
// the phased-position TSV to out. Returns the number of het sites
// emitted and the first error encountered.
//
// When opts.OutputPrefix is non-empty, three BAM files are written
// alongside the TSV stream (see the package doc comment). The input
// must be a SAM stream when BAM output is requested; the BAM writer
// needs the full header from the reader.
//
// The TSV output is emitted in upstream samtools `phase` schema (CC
// banner + PS / FL / M / EV / //). The implementation is a faithful
// port of reference_code/samtools/phase.c, including the bit-exact
// khash/ksort iteration order so EV-line ordering matches upstream on
// the canonical fixtures.
//
// The -b BAM split is wired through upstream's dump_aln routing on this
// path: after each block's phaseEmit, the per-read frag flags
// (phase/phased/flip/ambig) computed by the dynaprog+fragphase pass
// drive dispatch into <prefix>.{0,1,chimera}.bam exactly as phase.c
// does, including the drand48()-seeded is_flip shuffle and the
// evidence-less 50/50 routing (reproduced bit-for-bit by the in-tree
// drand48 port — see phase_drand48.go). The -A (DropAmbiguous) and -F
// (NoFixChimera) flags feed straight into that routing. The split is
// byte-validated against the upstream binary in TestLivePhaseBamSplit.
func Phase(in io.Reader, out io.Writer, opts PhaseOptions) (int, error) {
	if opts.BlockWindow == 0 {
		opts.BlockWindow = DefaultPhaseBlockWindow
	}
	if opts.MaxDepth == 0 {
		opts.MaxDepth = DefaultPhaseMaxDepth
	}
	// When the legacy V1 emission is requested (used by the existing
	// TestPhase_TwoHetsConsistentChain / SecondHetLabelOne / etc. unit
	// tests that assert on the PS-label TSV), short-circuit to the
	// greedy chainer. When UpstreamSchema is set (the default for the
	// CLI when invoked as `samtools phase`), use the byte-faithful
	// upstream pipeline.
	if !opts.UpstreamSchema {
		return phaseLegacyTSV(in, out, opts)
	}
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, fmt.Errorf("samtools phase: open input: %w", err)
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()

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
		// phase.c filters: skip BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP.
		// Note: BAM_FSUPPLEMENTARY is NOT in this list — phase.c will
		// pick up supplementary alignments.
		if rec.IsUnmapped() {
			continue
		}
		if rec.Flag&(sam.FlagSecondary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
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

	// Open the per-haplotype BAM writers up front if requested.
	var bamSplit *bamSplitWriter
	if opts.OutputPrefix != "" {
		bs, err := newBAMSplitWriter(opts.OutputPrefix, r.Header())
		if err != nil {
			return 0, err
		}
		defer bs.Close()
		bamSplit = bs
	}

	// Upstream samtools phase routes -b reads via drand48() and never
	// calls srand48(), so it runs from glibc's default seed state
	// (Xi = 0). Our drand48 port reproduces that exact sequence, which
	// is what makes the .0/.1/.chimera split agree with upstream
	// record-for-record (not merely up to a 0<->1 relabelling). The
	// opts.RNGSeed knob is therefore ignored on the upstream-schema
	// path — it only influences the legacy math/rand path below.
	rng := newDrand48()

	// CC banner emitted once at the start.
	emitPhaseBanner(bw)

	runner := &upstreamPhaseRunner{
		k:          DefaultPhaseBlockWindow,
		minBaseQ:   DefaultPhaseMinBaseQ,
		minVarLOD:  DefaultPhaseMinVarLOD,
		maxDepth:   opts.MaxDepth,
		fixChimera: !(opts.NoFixChimera || opts.FullRead),
	}
	if opts.MinBaseQ != 0 {
		runner.minBaseQ = opts.MinBaseQ
	}
	if opts.MinVarLOD >= 0 {
		runner.minVarLOD = opts.MinVarLOD
	}
	if opts.BlockWindow != 0 {
		runner.k = opts.BlockWindow
	}
	if opts.SiteListPath != "" {
		set, lerr := loadPhaseSites(opts.SiteListPath)
		if lerr != nil {
			return 0, fmt.Errorf("samtools phase: %w", lerr)
		}
		runner.siteSet = set
		runner.listExcl = opts.ListExclusive
	}

	emitted := 0
	// vposShift (the global het index) is carried across references: upstream
	// runs one continuous pileup and only resets at a tid change, just before
	// flushing the previous chromosome's trailing buffer. runUpstreamPhase
	// applies that reset before its final flush (for all but the last ref), so
	// the counter must start at 0 once here, not per reference.
	runner.vposShift = 0
	// One fragment hash for the WHOLE input, carried across references. Upstream
	// runs a single continuous pileup with one `seqs` table; at a chromosome
	// change it flushes the previous chr's trailing block then calls
	// update_vpos(0x7fffffff) to delete every fragment (leaving the table full of
	// tombstones at its current n_buckets) before the new chr's reads arrive.
	// Re-creating the table per reference would reset n_buckets to 0 and grow it
	// afresh, diverging the bucket layout — and therefore the unstable EV-line
	// emit order — from the second chromosome on.
	hash := newFragKhash()
	for ri, ref := range refOrder {
		recs := byRef[ref]
		sort.SliceStable(recs, func(i, j int) bool { return recs[i].Pos < recs[j].Pos })
		// runUpstreamPhase interleaves phaseEmit with dump_aln when
		// bamSplit is non-nil, matching upstream phase.c byte-for-byte
		// on reads with confident haplotype evidence. Evidence-less
		// reads use math/rand (seeded RNG) rather than upstream's
		// drand48 — see PhaseOptions.RNGSeed.
		n, err := runUpstreamPhase(runner, hash, recs, ref, ri == len(refOrder)-1, bw, bamSplit, rng, opts)
		if err != nil {
			return emitted, err
		}
		emitted += n
	}
	return emitted, nil
}

// phaseLegacyTSV is the v1 PS-label TSV emitter, preserved for the
// in-process unit tests that assert on the simpler schema. New
// callers should set PhaseOptions.UpstreamSchema = true.
func phaseLegacyTSV(in io.Reader, out io.Writer, opts PhaseOptions) (int, error) {
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, fmt.Errorf("samtools phase: open input: %w", err)
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()

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

	var bamSplit *bamSplitWriter
	if opts.OutputPrefix != "" {
		bs, err := newBAMSplitWriter(opts.OutputPrefix, r.Header())
		if err != nil {
			return 0, err
		}
		defer bs.Close()
		bamSplit = bs
	}

	seed := opts.RNGSeed
	if seed == 0 {
		seed = 1
	}
	rng := rand.New(rand.NewSource(seed))

	emitted := 0
	for _, ref := range refOrder {
		recs := byRef[ref]
		sort.SliceStable(recs, func(i, j int) bool { return recs[i].Pos < recs[j].Pos })
		hets, err := callHets(recs, opts)
		if err != nil {
			return emitted, err
		}
		phased, mapping := phaseHetsWithMapping(hets, opts)
		for _, h := range phased {
			line := fmt.Sprintf("PS\t%s\t%d\t%d\n", ref, h.pos, h.label)
			if _, err := bw.WriteString(line); err != nil {
				return emitted, err
			}
			emitted++
		}
		if bamSplit != nil {
			if err := bamSplit.assignAndWrite(recs, hets, mapping, opts, rng); err != nil {
				return emitted, err
			}
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
//   - the two most-common bases each have >= 2 supporting reads (a
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
	sites, _ := phaseHetsWithMapping(hets, opts)
	return sites
}

// hetHapMapping records, for each het, what the per-block "hap0"
// allele index is (0 means allele0 of this het represents hap0; 1
// means allele1 of this het represents hap0; -1 means the het is
// not phased into the current block).
type hetHapMapping []int8

// phaseHetsWithMapping is the workhorse for phaseHets. It returns
// the per-het labels and, alongside, the per-het hap0-allele
// mapping consistent with the block's labelling. The mapping is
// what callers need to classify a read by haplotype (see the
// chimera-repair pass in phase_bam.go).
func phaseHetsWithMapping(hets []het, opts PhaseOptions) ([]phasedSite, hetHapMapping) {
	out := make([]phasedSite, 0, len(hets))
	mapping := make(hetHapMapping, len(hets))
	for i := range mapping {
		mapping[i] = -1
	}
	// Per-read: which block-frame haplotype (0 or 1) was this read on
	// at the previous phased het? The block frame is fixed at the
	// first het of the block; "blockHap" is allele-index–independent
	// and stays consistent across flips so that a non-chimeric read's
	// assignment never appears to flip from het to het.
	readBlockHap := map[int]int{}
	prevPhasedIdx := -1
	consecUnphased := 0
	// cumFlip is the cumulative number of label-2 (flip) events seen
	// since the block start, mod 2. mapping[i] = cumFlip (after the
	// flip at het i is applied) — i.e. cumFlip == 0 ⇒ allele0
	// represents block hap0; cumFlip == 1 ⇒ allele1 represents block
	// hap0.
	cumFlip := 0

	for i, h := range hets {
		if prevPhasedIdx < 0 {
			// Block start: label the first het arbitrarily; record
			// each read's block-frame hap = its allele index here.
			out = append(out, phasedSite{pos: h.pos, label: 1})
			cumFlip = 0
			mapping[i] = 0
			for _, s := range h.support {
				readBlockHap[s.readIdx] = s.allele
			}
			prevPhasedIdx = i
			consecUnphased = 0
			continue
		}
		// Translate each read's local allele at this het into the
		// current block frame: blockHap = s.allele XOR cumFlip.
		// Compare against the prior block-frame hap of the same read.
		// "same" supports keeping the labelling (label=1, no flip);
		// "opposite" supports flipping the labelling (label=2). Note
		// that for a non-chimeric read the block-frame hap is
		// invariant across hets, so the count is robust against a
		// chimera read's mid-block jump.
		same, opposite := 0, 0
		for _, s := range h.support {
			prev, ok := readBlockHap[s.readIdx]
			if !ok {
				continue
			}
			curr := s.allele ^ cumFlip
			if prev == curr {
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
				readBlockHap = map[int]int{}
				consecUnphased = 0
				cumFlip = 0
			}
			continue
		}
		// Apply the flip (if any) to the running cumulative flip
		// count. Then update readBlockHap with each supporting read's
		// block-frame hap at this het (s.allele XOR cumFlip).
		if label == 2 {
			cumFlip ^= 1
		}
		mapping[i] = int8(cumFlip)
		for _, s := range h.support {
			readBlockHap[s.readIdx] = s.allele ^ cumFlip
		}
		prevPhasedIdx = i
		consecUnphased = 0
	}
	return out, mapping
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
	refPos := int32(rec.Pos) // 1-based
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
