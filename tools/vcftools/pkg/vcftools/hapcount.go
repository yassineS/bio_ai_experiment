// --hapcount BED: per-BED-bin haplotype-count summaries.
//
// Mirrors upstream variant_file_output.cpp:1169-1401
// (variant_file::output_haplotype_count). For each interval in the BED
// file the routine emits the number of SNPs that fell into the bin
// (N_SNP), the number of unique kept-individual haplotypes observed
// across those SNPs (N_UNIQ_HAPS), and the histogram of haplotype
// multiplicities (N_GROUPS + {COUNT:MULT} pairs).
//
// Upstream's BED is read with `BED.ignore(numeric_limits<streamsize>::max(),
// '\n')` (variant_file_output.cpp:1183) BEFORE the parse loop, so the
// FIRST line of the BED file is unconditionally skipped (treated as a
// mandatory header). Subsequent lines starting with `#` and blank lines
// are also dropped (lines 1188-1189). We replicate exactly so a BED
// file's first data row is silently dropped unless the user supplied a
// throwaway header row.
//
// Phased-only: upstream sets `phased_only=true` for --hapcount in
// parameters.cpp:248, so the global `passes phased site filter` gate in
// `Run` already drops `0/1`-style sites before they reach this runner.
//
// Diploid-only: upstream skips any site that fails `entry::is_diploid()`
// (line 1350-1354). We match by using the same `isFullyDiploid` helper
// shared with burden.go.
//
// Allele coordinates: upstream uses `e->get_indv_GENOTYPE_ids(ui, geno)`
// which returns the per-haplotype allele index (REF=0, ALT_i=i+1) or -1
// for missing (`.|.`). We use `diploidAlleles` which has the same
// convention.
//
// PRESERVED UPSTREAM BUG: lines 1314-1315 unconditionally assign
//
//	prev_bin_idx = bin_idx;
//	bin_idx      = ui;
//
// at the start of every per-site search loop's successful match. The
// flush-trigger predicate at line 1322 is then
//
//	if ((found == false) || (prev_bin_idx != bin_idx))
//
// — so after a real bin-change event (within a chromosome), the next
// per-site iteration leaves `prev_bin_idx` pointing at the OLD bin
// index, even though the data has moved on. The next time a flush
// fires (either at the next bin change or at the end-of-chromosome
// transition), `SNP_count[prev_bin_idx]` and
// `haplotype_count[prev_bin_idx]` get OVERWRITTEN with the new bin's
// values, AND `haplotype_frequencies[prev_bin_idx]` gets the new bin's
// histogram ADDED to the old bin's histogram. The N_GROUPS column
// reflects the union and the N_SNP / N_UNIQ_HAPS columns reflect the
// latest write. We replicate this byte-for-byte because that's what
// real users of upstream vcftools see in their `.hapcount` output.
//
// PRESERVED UPSTREAM TRUNCATION: lines 1370-1400 read `e->include_indv`
// AFTER `delete e;` (line 1370). In practice this read-after-free
// either segfaults or returns garbage on glibc systems; observed
// behaviour on the upstream binary in this repo is that the last
// chromosome's bins do NOT appear in the output file (the process
// exits with stdout/stderr partially flushed and the final-flush block
// produces no rows). We mirror this by SUPPRESSING the final-flush
// step entirely — only chromosomes followed by another chromosome in
// the input VCF are emitted. Test fixtures should end with a sentinel
// chromosome to force the final transition. Tracked in
// docs/UPSTREAM_BUGS.md alongside the indv-freq-burden label bug.
package vcftools

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// hapcountBin represents one (POS1, POS2) entry from the BED file. Membership
// for a 1-based VCF position p is `p > POS1 && p <= POS2` (upstream's
// half-open-from-the-right convention; line 1311).
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

	// binPositions holds the per-chromosome sorted bin list. Ordering
	// follows the BED-file first-appearance order so output chromosomes
	// match upstream's chr_to_idx assignment order (variant_file_output.cpp:1196-1201).
	chrOrder     []string
	chrToIdx     map[string]int
	binPositions [][]hapcountBin

	// Running state (mirrors upstream local variables).
	haplotypes           [][]int         // 2*keptIndvCount slots, each a per-bin allele sequence.
	haplotypeCount       [][]int         // per-chr per-bin int.
	snpCount             [][]int         // per-chr per-bin int.
	haplotypeFrequencies [][]map[int]int // per-chr per-bin multiplicity histogram.
	minUi                []int           // per-chr search-loop low water mark.
	binIdx               int
	prevBinIdx           int
	prevIdx              int
	prevCHROM            string
	haveData             bool

	// out is the destination stream; we keep it on the runner so the
	// last per-chrom transition writes incrementally rather than
	// buffering everything.
	out       io.Writer
	closer    io.Closer
	headerOut bool
}

// newHapcountRunner parses the BED file (with the upstream-mandatory
// first-line-skip), validates that no chromosome has overlapping bins
// (upstream lines 1208-1216 LOG.error on overlap), and prepares the
// per-chromosome accumulators. Returns an error on missing/unreadable
// BED file, malformed BED rows, or overlapping bins.
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
	// Upstream-mandatory first-line-skip (variant_file_output.cpp:1183).
	if scanner.Scan() {
		// header line discarded.
	}
	lineNum := 1
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}
		if line[0] == '#' {
			continue
		}
		// `ss >> CHROM >> POS1 >> POS2`: whitespace-separated, ignore
		// trailing fields.
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

// open lazily creates the output file (or stdout) and writes the header
// row. Called from the first transition or close().
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
// at variant_file_output.cpp:1247-1369. Returns an error only on output
// failure during a chrom-transition flush.
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
		// Chrom transition — flush the previous chromosome and emit its
		// rows. Mirrors lines 1265-1305.
		if err := r.flushChromTransition(idx); err != nil {
			return err
		}
	}

	// Per-site BED-bin search (lines 1307-1320). Upstream uses
	// `min_ui[idx]` as a low-water mark so subsequent sites don't
	// re-scan bins already passed.
	found := false
	maxUi := len(r.binPositions[idx])
	for ui := r.minUi[idx]; ui < maxUi; ui++ {
		if pos > r.binPositions[idx][ui].pos1 && pos <= r.binPositions[idx][ui].pos2 {
			found = true
			// PRESERVED UPSTREAM BUG (see file-level docs): assign
			// prev_bin_idx BEFORE bin_idx, unconditionally — the bug
			// shifts data into the wrong bin's accumulators when the
			// next flush fires.
			r.prevBinIdx = r.binIdx
			r.binIdx = ui
			break
		} else if pos > r.binPositions[idx][ui].pos2 {
			r.minUi[idx] = ui + 1
		}
	}

	// Bin-change flush (lines 1322-1343).
	if !found || r.prevBinIdx != r.binIdx {
		if r.haveData {
			r.flushBin(r.prevIdx, r.prevBinIdx)
		}
		r.haveData = false
	}

	if found {
		r.haveData = true
		// Diploid-only guard (upstream `is_diploid` check at line 1350).
		if !isFullyDiploid(v) {
			return nil
		}
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
// clears the haplotypes slices. Mirrors upstream's per-flush inner block
// at lines 1324-1341. Updates SNP_count, haplotype_count, and
// haplotype_frequencies in place (haplotype_frequencies is ACCUMULATED
// — see the file-level bug note).
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
		// Upstream takes the size of haplotypes[2*ui] (line 1334).
		// This is the same for every ui within the same flush (every
		// individual's haplotype was push_back'd at every accepted
		// site), so we capture it once.
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
// for use as a map key. We use a binary-safe encoding so missing
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
		r.flushBin(r.prevIdx, r.prevBinIdx)
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
		// First three cols + N_SNP + N_UNIQ_HAPS — upstream emits the
		// N_GROUPS column separately. Multiplicity:Frequency pairs are
		// emitted in ascending-multiplicity order (std::map iteration
		// order over int keys).
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

// close finalises output. Mirrors upstream's read-after-free
// final-flush block (lines 1370-1400) by emitting the last chromosome's
// bins ONLY when `have_data == false` at end-of-stream. When
// `have_data == true`, upstream's freed-pointer access truncates the
// output file silently — we replicate that by NOT emitting the
// last-chrom rows in that case. See the file-level "PRESERVED UPSTREAM
// TRUNCATION" note for the empirical evidence behind this.
//
// The header row is always written so an empty / wholly-non-matching
// input still produces a parseable one-line output file.
func (r *hapcountRunner) close() error {
	if r == nil {
		return nil
	}
	if err := r.open(); err != nil {
		return err
	}
	// Mirror upstream's `if (idx == prev_idx)` final-flush. The only
	// path that reaches output is have_data == false (otherwise
	// upstream's read-after-free truncates the file). When
	// have_data == false, `haplotype_set` stays empty, so
	// haplotype_count / haplotype_frequencies for the bin are not
	// touched — the existing per-bin defaults (zero) are emitted.
	if !r.haveData && r.prevIdx >= 0 && r.prevIdx < len(r.binPositions) {
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
