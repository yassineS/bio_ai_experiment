// Package vcftools — --hapcount BED: per-BED-bin haplotype-count summaries.
//
// Ported from upstream variant_file_output.cpp:1169-1401
// (variant_file::output_haplotype_count). For each interval in the BED
// file the routine emits the number of SNPs that fell into the bin
// (N_SNP), the number of unique kept-individual haplotypes observed
// across those SNPs (N_UNIQ_HAPS), and the histogram of haplotype
// multiplicities (N_GROUPS + {COUNT:MULT} pairs).
//
// # Upstream bugs FIXED in this port
//
// Per project policy (CLAUDE.md "don't replicate upstream bugs") and the
// wave-14 precedent (PR #138, fix-on-port for the .ifreqburden INDV-
// label bug), three upstream defects in `output_haplotype_count` are
// CORRECTED here rather than mirrored. See docs/UPSTREAM_BUGS.md for
// the corresponding writeup.
//
//  1. prev_bin_idx shift on bin change (upstream lines 1314-1315). The
//     unconditional `prev_bin_idx = bin_idx; bin_idx = ui;` at the start
//     of every successful bin-match leaves `prev_bin_idx` pointing at
//     the OLD bin AFTER a within-chromosome bin-transition flush has
//     already fired. The next flush then overwrites SNP_count /
//     haplotype_count for the OLD bin with the NEW bin's values and
//     ACCUMULATES the new bin's multiplicity histogram into the old
//     bin's `haplotype_frequencies`. Concretely, in our chr-1 fixture
//     bin (100,500] should report N_SNP=4 but the buggy upstream binary
//     reports N_SNP=1 (the value from bin (1000,2000]).
//
//     Fix: when a within-chromosome bin transition is detected, flush
//     the OLD bin's data into its OWN slot BEFORE reassigning bin_idx.
//
//  2. End-of-stream read-after-free (upstream lines 1370-1400). The
//     final-flush block reads `e->include_indv[ui]` AFTER `delete e;`
//     at line 1370. Observed behaviour with a glibc-built upstream
//     binary: when `have_data == true` the buffered last chromosome is
//     silently dropped from the output file; when `have_data == false`
//     the read-after-free skips the inner loop and emits all-zero rows.
//
//     Fix: run the chrom-transition flush + per-chrom emit
//     unconditionally at end-of-stream using the kept-sample list we
//     own (no freed-pointer access).
//
//  3. BED first-line silent skip (upstream line 1183:
//     `BED.ignore(numeric_limits<streamsize>::max(), '\n')`). The first
//     BED line is unconditionally discarded, even when it is real data
//     (no header). A user with a header-less BED silently loses one
//     bin.
//
//     Fix: auto-detect headers. Skip the first line ONLY if it starts
//     with `#`, `track`, `browser`, or is blank. Otherwise treat the
//     line as data.
package vcftools

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// hapcountBin represents one (POS1, POS2) entry from the BED file.
// Membership for a 1-based VCF position p is `p > POS1 && p <= POS2`
// (upstream's half-open-from-the-right convention at line 1311).
type hapcountBin struct {
	pos1 int
	pos2 int
}

// hapcountRunner accumulates the per-bin haplotype tables and emits
// `<prefix>.hapcount` at end-of-stream.
type hapcountRunner struct {
	prefix string

	// keptIndvCount is the number of kept individuals, recorded at
	// runner construction. The runner sizes the per-haplotype slice to
	// `2 * keptIndvCount`.
	keptIndvCount int

	// chrOrder + chrToIdx + binPositions track the BED file in
	// first-appearance chromosome order (mirroring upstream's
	// chr_to_idx assignment at lines 1196-1201). Per-chr bins are
	// sorted by (pos1, pos2) after parsing.
	chrOrder     []string
	chrToIdx     map[string]int
	binPositions [][]hapcountBin

	// Per-chromosome accumulators (sized to the chromosome's bin
	// count). One slot per (chrIdx, binIdx).
	haplotypeCount       [][]int         // unique haplotype count per bin.
	snpCount             [][]int         // N_SNP per bin.
	haplotypeFrequencies [][]map[int]int // multiplicity -> #groups histogram per bin.
	minUi                []int           // per-chr search-loop low-water mark.

	// haplotypes is a per-active-bin allele vector for each kept
	// haplotype (2 * keptIndvCount slots). Cleared inside flushBin.
	haplotypes [][]int

	// Running state mirroring upstream local variables.
	binIdx     int
	prevBinIdx int
	prevIdx    int
	prevCHROM  string
	haveData   bool

	// out is the destination stream; we keep it on the runner so the
	// last per-chrom transition writes incrementally rather than
	// buffering everything.
	out       io.Writer
	closer    io.Closer
	headerOut bool
}

// shouldSkipBEDHeader reports whether the supplied first BED line looks
// like header / commentary rather than data. Returns true for blank
// lines and for lines beginning with `#`, `track`, or `browser` (the
// three common BED-header conventions). Otherwise returns false so the
// caller treats the line as a data row.
//
// This is the FIX for upstream bug #3 — the unconditional first-line
// skip at variant_file_output.cpp:1183 silently dropped header-less BED
// data.
func shouldSkipBEDHeader(line string) bool {
	trimmed := strings.TrimSpace(line)
	if trimmed == "" {
		return true
	}
	if trimmed[0] == '#' {
		return true
	}
	if strings.HasPrefix(trimmed, "track") {
		return true
	}
	if strings.HasPrefix(trimmed, "browser") {
		return true
	}
	return false
}

// newHapcountRunner parses the BED file (auto-detecting headers — see
// shouldSkipBEDHeader), validates that no chromosome has overlapping
// bins (upstream lines 1208-1216 LOG.error on overlap), and prepares
// the per-chromosome accumulators. Returns an error on missing /
// unreadable BED, malformed BED rows, or overlapping bins.
func newHapcountRunner(prefix, bedPath string, keptIndvCount int) (*hapcountRunner, error) {
	f, err := iohelper.OpenReader(bedPath)
	if err != nil {
		return nil, fmt.Errorf("Could not open BED file: %s", bedPath)
	}
	defer f.Close()

	chrOrder := []string{}
	chrToIdx := map[string]int{}
	binPositions := [][]hapcountBin{}

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	lineNum := 0
	first := true
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()
		if first {
			first = false
			if shouldSkipBEDHeader(line) {
				continue
			}
			// fall through: treat as data.
		}
		if len(line) == 0 {
			continue
		}
		if line[0] == '#' {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 3 {
			return nil, fmt.Errorf("--hapcount: malformed BED line %d in %s: %q", lineNum, bedPath, line)
		}
		pos1, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("--hapcount: invalid POS1 on BED line %d in %s: %v", lineNum, bedPath, err)
		}
		pos2, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("--hapcount: invalid POS2 on BED line %d in %s: %v", lineNum, bedPath, err)
		}
		chrom := fields[0]
		idx, ok := chrToIdx[chrom]
		if !ok {
			idx = len(chrOrder)
			chrOrder = append(chrOrder, chrom)
			chrToIdx[chrom] = idx
			binPositions = append(binPositions, nil)
		}
		binPositions[idx] = append(binPositions[idx], hapcountBin{pos1: pos1, pos2: pos2})
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("--hapcount: reading BED file %s: %w", bedPath, err)
	}

	// Sort + non-overlap check (upstream lines 1208-1216).
	for i := range binPositions {
		sort.Slice(binPositions[i], func(a, b int) bool {
			if binPositions[i][a].pos1 != binPositions[i][b].pos1 {
				return binPositions[i][a].pos1 < binPositions[i][b].pos1
			}
			return binPositions[i][a].pos2 < binPositions[i][b].pos2
		})
		for j := 1; j < len(binPositions[i]); j++ {
			if binPositions[i][j-1].pos2 > binPositions[i][j].pos1 {
				return nil, fmt.Errorf("BED file must be non-overlapping.")
			}
		}
	}

	r := &hapcountRunner{
		prefix:               prefix,
		keptIndvCount:        keptIndvCount,
		chrOrder:             chrOrder,
		chrToIdx:             chrToIdx,
		binPositions:         binPositions,
		haplotypes:           make([][]int, 2*keptIndvCount),
		haplotypeCount:       make([][]int, len(chrOrder)),
		snpCount:             make([][]int, len(chrOrder)),
		haplotypeFrequencies: make([][]map[int]int, len(chrOrder)),
		minUi:                make([]int, len(chrOrder)),
		binIdx:               0,
		prevBinIdx:           -1,
		prevIdx:              -1,
	}
	for i := range r.haplotypeCount {
		r.haplotypeCount[i] = make([]int, len(binPositions[i]))
		r.snpCount[i] = make([]int, len(binPositions[i]))
		r.haplotypeFrequencies[i] = make([]map[int]int, len(binPositions[i]))
		for j := range r.haplotypeFrequencies[i] {
			r.haplotypeFrequencies[i][j] = make(map[int]int)
		}
	}
	for i := range r.haplotypes {
		r.haplotypes[i] = make([]int, 0, 16)
	}
	return r, nil
}

// open lazily creates the output file and writes the header row. Called
// from the first chromosome transition or from close().
func (r *hapcountRunner) open() error {
	if r.headerOut {
		return nil
	}
	f, err := iohelper.OpenWriter(r.prefix + ".hapcount")
	if err != nil {
		return fmt.Errorf("Could not open Haplotype Output File: %s", r.prefix+".hapcount")
	}
	r.closer = f
	r.out = bufio.NewWriter(f)
	if _, err := io.WriteString(r.out, "#CHROM\tBIN_START\tBIN_END\tN_SNP\tN_UNIQ_HAPS\tN_GROUPS\t{MULTIPLICITY:FREQ}\n"); err != nil {
		return err
	}
	r.headerOut = true
	return nil
}

// addVariant processes one already-filtered (phased, kept-sample,
// passes-all-filters) variant. Mirrors upstream's main per-site block
// at variant_file_output.cpp:1247-1369 but with the prev_bin_idx-shift
// FIX (upstream bug #1) applied: on a within-chromosome bin transition
// we flush the OLD bin's accumulators into the OLD bin's slot BEFORE
// reassigning bin_idx.
func (r *hapcountRunner) addVariant(v *vcf.Variant) error {
	if r == nil {
		return nil
	}
	chrom := v.Chrom
	idx, ok := r.chrToIdx[chrom]
	if !ok {
		// Upstream's `continue` at line 1262 — chrom not in BED, drop.
		return nil
	}
	pos := v.Pos

	if idx != r.prevIdx {
		// Chrom transition — flush the previous chromosome and emit
		// its rows. Mirrors lines 1265-1305.
		if err := r.flushChromTransition(idx); err != nil {
			return err
		}
	}

	// Per-site BED-bin search (lines 1307-1320). Upstream uses
	// `min_ui[idx]` as a low-water mark so subsequent sites don't
	// re-scan bins already passed.
	found := false
	var matchedUi int
	maxUi := len(r.binPositions[idx])
	for ui := r.minUi[idx]; ui < maxUi; ui++ {
		if pos > r.binPositions[idx][ui].pos1 && pos <= r.binPositions[idx][ui].pos2 {
			found = true
			matchedUi = ui
			break
		} else if pos > r.binPositions[idx][ui].pos2 {
			r.minUi[idx] = ui + 1
		}
	}

	// FIX for upstream bug #1: when a bin transition is detected,
	// flush the OLD bin's accumulators into the OLD bin's slot BEFORE
	// reassigning bin_idx. The reviewer-recommended formulation.
	if found {
		if r.haveData && r.binIdx != matchedUi {
			// Bin transition: flush the OLD bin first.
			r.flushBin(r.prevIdx, r.binIdx)
			r.haveData = false
		}
		r.prevBinIdx = r.binIdx
		r.binIdx = matchedUi
	} else if r.haveData {
		// Site fell outside every bin AND we had data in flight ->
		// flush the in-flight bin (bin_idx) into its slot.
		r.flushBin(r.prevIdx, r.binIdx)
		r.haveData = false
	}

	if found {
		// Diploid-only guard (upstream `is_diploid` check at line 1350).
		if !isFullyDiploid(v) {
			return nil
		}
		r.haveData = true
		// Push per-individual alleles into the haplotype slices.
		// Upstream iterates `for ui=0..N_indv` filtered by
		// include_indv; our `v.Samples` is already the kept-sample
		// list, so we iterate directly.
		for ui := 0; ui < r.keptIndvCount; ui++ {
			a1, a2 := -1, -1
			if ui < len(v.Samples) {
				if x, y, ok := diploidAlleles(getGT(v, ui)); ok {
					a1 = x
					a2 = y
				}
			}
			r.haplotypes[2*ui] = append(r.haplotypes[2*ui], a1)
			r.haplotypes[2*ui+1] = append(r.haplotypes[2*ui+1], a2)
		}
	}
	return nil
}

// flushBin computes per-bin haplotype_set + multiplicity histogram for
// (chrIdx, binIdx) using the current contents of `haplotypes`, then
// clears the haplotypes slices. Mirrors upstream's per-flush inner
// block at lines 1324-1341 — except SNP_count / haplotype_count for the
// flushed bin are SET, not overwritten with stale prev_bin_idx data.
func (r *hapcountRunner) flushBin(chrIdx, binIdx int) {
	if chrIdx < 0 || chrIdx >= len(r.snpCount) {
		return
	}
	if binIdx < 0 || binIdx >= len(r.snpCount[chrIdx]) {
		return
	}
	hapSet := make(map[string]int)
	snpLen := 0
	for ui := 0; ui < r.keptIndvCount; ui++ {
		left := r.haplotypes[2*ui]
		right := r.haplotypes[2*ui+1]
		hapSet[encodeHap(left)]++
		hapSet[encodeHap(right)]++
		// snpLen is the same for every ui within the same flush (every
		// individual gets one append at each accepted site), so we
		// capture it once.
		snpLen = len(left)
		r.haplotypes[2*ui] = r.haplotypes[2*ui][:0]
		r.haplotypes[2*ui+1] = r.haplotypes[2*ui+1][:0]
	}
	r.snpCount[chrIdx][binIdx] = snpLen
	r.haplotypeCount[chrIdx][binIdx] = len(hapSet)
	for _, mult := range hapSet {
		r.haplotypeFrequencies[chrIdx][binIdx][mult]++
	}
}

// encodeHap serialises a haplotype-allele vector to a string suitable
// for use as a map key. We use a comma-separated encoding so missing
// alleles (-1) are distinguishable from real allele 0.
func encodeHap(a []int) string {
	var buf bytes.Buffer
	for _, x := range a {
		_, _ = fmt.Fprintf(&buf, "%d,", x)
	}
	return buf.String()
}

// flushChromTransition handles the chrom-change branch at upstream
// lines 1265-1305: flush any in-flight bin into prev_chr's
// accumulators, emit ALL of prev_chr's bin rows, then reset the
// per-chr accumulators for the new chromosome.
func (r *hapcountRunner) flushChromTransition(newIdx int) error {
	if err := r.open(); err != nil {
		return err
	}
	if r.haveData {
		r.flushBin(r.prevIdx, r.binIdx)
	}
	r.haveData = false

	if r.prevIdx >= 0 && r.prevIdx < len(r.binPositions) {
		if err := r.writeChrom(r.prevIdx, r.prevCHROM); err != nil {
			return err
		}
	}

	// Reset for new chrom (upstream lines 1298-1304).
	r.binIdx = 0
	r.prevBinIdx = -1
	r.prevIdx = newIdx
	r.prevCHROM = r.chrOrder[newIdx]
	return nil
}

// writeChrom emits every bin of `chrIdx` to the output stream, using
// `chrom` as the leading column. Mirrors lines 1287-1295.
func (r *hapcountRunner) writeChrom(chrIdx int, chrom string) error {
	for ui := 0; ui < len(r.haplotypeCount[chrIdx]); ui++ {
		bin := r.binPositions[chrIdx][ui]
		var sb strings.Builder
		fmt.Fprintf(&sb, "%s\t%d\t%d\t%d\t%d\t%d", chrom, bin.pos1, bin.pos2,
			r.snpCount[chrIdx][ui], r.haplotypeCount[chrIdx][ui],
			len(r.haplotypeFrequencies[chrIdx][ui]))
		// Sort multiplicities ascending for std::map<int,int> parity.
		mults := make([]int, 0, len(r.haplotypeFrequencies[chrIdx][ui]))
		for m := range r.haplotypeFrequencies[chrIdx][ui] {
			mults = append(mults, m)
		}
		sort.Ints(mults)
		for _, m := range mults {
			fmt.Fprintf(&sb, "\t%d:%d", r.haplotypeFrequencies[chrIdx][ui][m], m)
		}
		sb.WriteByte('\n')
		if _, err := io.WriteString(r.out, sb.String()); err != nil {
			return err
		}
	}
	return nil
}

// close finalises output. FIX for upstream bug #2: the last seen
// chromosome's bins are emitted UNCONDITIONALLY (after a final in-flight
// bin flush if any), without ever reading from freed memory. The header
// row is always written so an empty / wholly-non-matching input still
// produces a parseable one-line output file.
func (r *hapcountRunner) close() error {
	if r == nil {
		return nil
	}
	if err := r.open(); err != nil {
		return err
	}
	if r.haveData {
		r.flushBin(r.prevIdx, r.binIdx)
		r.haveData = false
	}
	if r.prevIdx >= 0 && r.prevIdx < len(r.binPositions) {
		if err := r.writeChrom(r.prevIdx, r.prevCHROM); err != nil {
			return err
		}
	}
	if bw, ok := r.out.(*bufio.Writer); ok {
		if err := bw.Flush(); err != nil {
			return err
		}
	}
	if r.closer != nil {
		return r.closer.Close()
	}
	return nil
}
