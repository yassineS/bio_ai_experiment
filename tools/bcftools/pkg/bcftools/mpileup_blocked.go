// bcftools mpileup — block-streaming driver for the indexed `-r` fast path.
//
// Background. The whole-buffer path (MpileupFile → writeMpileupVCF →
// emitChromMpileup) drains EVERY read of a chromosome into a
// map[chrom][]*sam.Record and then runs whole-contig passes over it.
// Two of those passes allocate a SECOND contig-sized structure:
// applyMpileupBAQ's `column := make([][]*mpileupReadBAQInfo, maxPos)`
// (sized to the contig) and the three full-Qual snapshot maps
// (snapshotRawQuals / snapshotFirstMateQuals / snapshotPostMergeQuals).
// On human chr20 a whole-contig `-r 20` peaked at ~6.6 GB versus
// upstream's ~106 MB.
//
// Fix. When every input is an indexed, coordinate-sorted BAM (the
// region-reader fast path applies) and the simple one-column-per-BAM
// sample model is in force, this driver processes each contig/region in
// fixed-grid coordinate BLOCKS. For each block [blockBeg, blockEnd) it:
//
//   - reads ONLY the records overlapping
//     [blockBeg-leftFlank, blockEnd+rightFlank) via a fresh indexed
//     region read. leftFlank is wide enough that every in-block read's
//     UPSTREAM mate co-resides and the depth-cap / overlap / BAQ /
//     snapshot state has converged by blockBeg; rightFlank is the
//     symmetric guarantee that every in-block read's DOWNSTREAM mate and
//     every read the indel-realignment window reaches co-reside by
//     blockEnd (so a column near the right edge sees the same overlap /
//     snapshot / indel state a whole-contig run would);
//   - runs the UNCHANGED per-position pipeline (depth cap, BAQ, overlap
//     merge, indel) over that block-local read set — so the BAQ column
//     index and the Qual snapshots are sized to the BLOCK, not the
//     contig;
//   - emits records only for columns in [blockBeg, blockEnd); the left
//     flank columns [blockBeg-leftFlank, blockBeg) and the right flank
//     columns [blockEnd, blockEnd+rightFlank) are processed for context
//     only (the adjacent blocks emit them);
//   - drops the block's reads before the next block.
//
// Byte-exactness. The per-position math, the events-window grid, and the
// output format are all untouched — only WHICH reads and WHICH columns
// are resident per block change. The left flank guarantees that, by the
// time iteration reaches blockBeg, every read overlapping blockBeg has
// had its mate paired, its BAQ applied at its true first eligible column,
// its overlap quals merged and its pre/post-merge snapshots taken exactly
// as a whole-contig run would — so the columns this block emits are
// byte-identical to the whole-buffer path. The shared *biasLeak threads
// the SNP→indel bias scalar across blocks and chromosomes just as the
// whole-buffer path threads it across positions.
package bcftools

import (
	"fmt"
	"io"
	"runtime/debug"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// mpileupBlockWidth is the coordinate width (reference bases) of one
// streaming block. The block-local read buffer, BAQ column index and Qual
// snapshots all scale with this rather than with the whole contig, so it
// is the dominant peak-RSS knob. 0.2 Mbp keeps the resident set a few
// hundred MB — well inside upstream's envelope (a whole human chr20
// `-r 20` peaks ~680 MB ≈ 6x upstream's 106 MB, down from ~6.6 GB) — while
// the left-flank re-read overhead (mpileupBlockFlank / mpileupBlockWidth)
// stays a couple of percent. It is a var (not a const) only so tests can
// shrink it to exercise block boundaries on small synthetic inputs;
// production never changes it.
var mpileupBlockWidth = 200_000

// mpileupBlockFlank is the left context span (reference bases) read before
// each block so every in-block read's mate co-resides and the stateful
// left-to-right passes (depth cap, overlap pairing, BAQ, snapshots) have
// converged to their whole-contig values by blockBeg. It must comfortably
// exceed maxReadLen + the largest proper-pair insert. Only proper-pair,
// close mates need co-residence: the overlap classifier already excludes
// far/discordant mates (isize >= 2*readLen && mpos >= end) from buffering,
// so a large |TLEN| does not require a wide flank. 5 kbp is ~50 Illumina
// read-lengths and many times the largest observed proper-pair insert,
// guaranteeing both mate co-residence and depth-cap convergence while
// costing only a couple of percent extra reads at a 0.2 Mbp block width.
// A var (not a const) only so tests can shrink it alongside
// mpileupBlockWidth.
var mpileupBlockFlank = 5_000

// mpileupBlockRightFlank is the RIGHT context span (reference bases) read
// AFTER each block so the per-column state at columns near blockEnd
// converges to its whole-contig (whole-buffer) value. It is the symmetric
// analogue of mpileupBlockFlank: just as the left flank guarantees that
// every in-block read's UPSTREAM mate co-resides (so overlap pairing,
// snapshots and the depth cap have converged by blockBeg), the right flank
// guarantees that every in-block read's DOWNSTREAM mate and every read the
// indel-realignment window reaches co-reside by blockEnd.
//
// Why a right flank is needed at all. A column C in [blockBeg, blockEnd) is
// emitted by this block. Two whole-contig effects reach DOWNSTREAM of C and
// therefore need reads that start after C to be resident:
//
//  1. Smart-overlap merge. A read R that piles at C can be the FIRST mate of
//     a proper pair whose SECOND mate starts downstream of C (up to one
//     insert length past C). When that second mate is resident,
//     classifyMatePairs + applySmartOverlaps merge the pair's overlapping
//     quals and the pre/post-merge snapshots differ from the no-mate case —
//     visibly shifting the indel QS/PL of a column just left of the pair's
//     overlap (proven on the GIAB indel 20:30142484, where the whole-buffer
//     result only stabilises once reads ~76 bp downstream — the second mates
//     — are resident). Without the right flank a block truncates those
//     second mates at blockEnd and the right-edge column's overlap state
//     diverges from the whole-buffer path.
//  2. Indel realignment window. bcfCallGapPrep aligns each read over the
//     reference window [pos-IndelWinSize, pos+IndelWinSize-types[0]); with
//     the default IndelWinSize (DefaultMpileupIndelSize = 110) and a
//     deletion this reaches at most ~110 bases + the max deletion past C,
//     i.e. < a read length, far inside one insert.
//
// Sizing. Both effects are bounded by the largest proper-pair insert plus
// one read length downstream of C — exactly the quantity mpileupBlockFlank
// already covers on the left. The indel window (110 + maxDel + a read
// length to align it) is the smaller of the two and is also covered. 5 kbp
// is ~50 Illumina read-lengths and many times the largest observed
// proper-pair insert, so it provably covers both the downstream-mate and
// indel-window reach for every emitted column while costing only a couple of
// percent extra reads at a 0.2 Mbp block width. A var (not a const) only so
// tests can shrink it alongside mpileupBlockWidth.
var mpileupBlockRightFlank = 5_000

// mpileupBlockInput names one indexed BAM column for the block driver:
// the file path (re-opened per block as a fresh region reader) and the
// output sample name. The block driver never holds a reader open across
// blocks — it re-seeks the index for each block's [flank, end) span.
type mpileupBlockInput struct {
	path   string
	sample string
}

// writeMpileupVCFBlocked is the block-streaming equivalent of
// writeMpileupVCF for the indexed one-column-per-BAM case. It shares the
// header / writer / errmod / biasLeak setup, then walks each chromosome's
// requested windows in fixed coordinate blocks, re-reading each block's
// reads (block + left flank) from the index and emitting only the block's
// own columns. regWindows is the resolved -r/-R/-t/-T post-filter (never
// nil here: the block path is only taken when an indexed seek region
// exists).
func writeMpileupVCFBlocked(out io.Writer, opts MpileupOptions, ref *fasta.RandomAccess,
	chromOrder []string, chromLen map[string]int,
	inputs []mpileupBlockInput, samples []string,
	regWindows map[string][][2]int) error {

	hdr := buildMpileupHeader(opts, chromOrder, chromLen, samples)
	w, finish, err := openMpileupOutput(out, opts, hdr)
	if err != nil {
		return err
	}
	defer finish()
	if len(opts.GVCFRange) > 0 {
		w = newGVCFBlocker(w, opts.GVCFRange)
	}
	if err := w.WriteHeader(); err != nil {
		return err
	}

	em := errmod.Init(1.0 - mpileupTheta)
	// One biasLeak threaded across every block and chromosome, exactly as
	// the whole-buffer path threads it across every position (see the
	// long comment in writeMpileupVCF).
	leak := biasLeak{bq: 0, mqs: 0, bqOK: true, mqsOK: true}

	nIn := len(inputs)
	for _, chrom := range chromOrder {
		windows, ok := regWindows[chrom]
		if !ok || len(windows) == 0 {
			continue
		}
		refLen := chromLen[chrom]
		if refLen <= 0 {
			continue
		}
		refSlab, err := ref.Fetch(chrom, 0, int64(refLen))
		if err != nil {
			return fmt.Errorf("bcftools mpileup: fetch %s: %w", chrom, err)
		}
		// windows are 1-based inclusive, pre-sorted and merged by
		// parseMpileupRegions. Tile each window into fixed-grid blocks.
		for _, win := range windows {
			// Convert to 0-based half-open and clamp to the contig.
			wBeg0 := win[0] - 1
			if wBeg0 < 0 {
				wBeg0 = 0
			}
			wEnd0 := win[1] // win[1] is inclusive 1-based == exclusive 0-based
			if wEnd0 > refLen {
				wEnd0 = refLen
			}
			if wBeg0 >= wEnd0 {
				continue
			}
			for blockBeg := wBeg0; blockBeg < wEnd0; blockBeg += mpileupBlockWidth {
				blockEnd := blockBeg + mpileupBlockWidth
				if blockEnd > wEnd0 {
					blockEnd = wEnd0
				}
				flankBeg := blockBeg - mpileupBlockFlank
				if flankBeg < 0 {
					flankBeg = 0
				}
				// The right flank is the symmetric analogue of the left
				// flank: it pulls in the downstream mates and indel-window
				// reads that columns near blockEnd need so their overlap /
				// snapshot / indel state converges to the whole-buffer value
				// (see mpileupBlockRightFlank). These columns are CONTEXT
				// ONLY — the emit gate below stays [blockBeg, blockEnd) so
				// they are never written and biasLeak is never updated for
				// them.
				//
				// Clamp the right flank to the REGION end (wEnd0), never just
				// the contig length. wEnd0 is the user's requested region
				// right boundary (0-based exclusive). For an INTERNAL block
				// edge (blockEnd < wEnd0) the flank still extends past
				// blockEnd — up to min(blockEnd+flank, wEnd0) — so a mid-region
				// indel at an internal boundary keeps its downstream context
				// (the attempt-2 fix). For the FINAL block of the region
				// (blockEnd == wEnd0) this clamps the flank to wEnd0, so the
				// last emitted column at the region boundary uses exactly the
				// reads a region-bounded query [X, wEnd0] sees rather than
				// reading downstream context PAST the user's requested region
				// end. A whole-contig query has wEnd0 == refLen, so nothing is
				// ever past it and the clamp is a no-op there.
				flankEnd := blockEnd + mpileupBlockRightFlank
				if flankEnd > wEnd0 {
					flankEnd = wEnd0
				}
				// Read this block's reads (block + left flank + right flank)
				// for every input via a fresh indexed region read. Region
				// strings are 1-based inclusive.
				regionStr := fmt.Sprintf("%s:%d-%d", chrom, flankBeg+1, flankEnd)
				perInputChromRecs := make([][]*sam.Record, nIn)
				anyHit := false
				for i := range inputs {
					recs, rerr := readMpileupBlockRecords(inputs[i].path, regionStr, opts)
					if rerr != nil {
						return fmt.Errorf("bcftools mpileup: %s: %w", inputs[i].path, rerr)
					}
					perInputChromRecs[i] = recs
					if len(recs) > 0 {
						anyHit = true
					}
				}
				if !anyHit {
					continue
				}
				// Emit only this block's own columns [blockBeg, blockEnd);
				// the flank columns are context. regWindows still post-
				// filters inside emitChromMpileup so partial-window edges
				// and multi-window chroms are handled exactly as before.
				if err := emitChromMpileup(w, em, chrom, refSlab, refLen,
					perInputChromRecs, opts, regWindows, &leak, blockBeg, blockEnd); err != nil {
					return err
				}
				// Drop the block's reads so the next block's allocation
				// reuses the freed heap rather than growing it.
				for i := range perInputChromRecs {
					perInputChromRecs[i] = nil
				}
				perInputChromRecs = nil
				// Return the freed block heap to the OS. Each block churns
				// hundreds of MB of *sam.Record (Seq/Qual/Cigar/Name slices)
				// plus a block-sized BAQ column index and Qual snapshots;
				// without this the Go runtime lets the heap high-water-mark
				// grow across blocks (peak RSS stayed multi-GB even though
				// only one block is live), which is exactly the whole-contig
				// blow-up this driver exists to avoid. FreeOSMemory runs a GC
				// and madvise-frees the span, bounding peak RSS to ~one
				// block's footprint at the cost of one GC per block (a few MB
				// of work amortised over millions of columns).
				debug.FreeOSMemory()
			}
		}
	}
	return w.Flush()
}

// readMpileupBlockRecords opens a fresh indexed region reader for one
// input over regionStr, applies the same record-level filters as
// mpileupReadBAM (flag bits, MAPQ via mpileupKeepRecord), and returns the
// kept records sorted by start position. It mirrors mpileupReadBAM's
// per-record handling exactly so a block's read set is identical to the
// corresponding slice of the whole-buffer path's bucketed reads.
func readMpileupBlockRecords(path, regionStr string, opts MpileupOptions) ([]*sam.Record, error) {
	rr, err := openBAMRegionReader(path, []string{regionStr})
	if err != nil {
		return nil, err
	}
	if rr == nil {
		// Should not happen: the block path is only entered when the
		// indexed reader was available. Be defensive and signal the
		// caller to fall back rather than silently emit nothing.
		return nil, fmt.Errorf("indexed region reader unavailable for %s", path)
	}
	defer rr.Close()
	var out []*sam.Record
	for {
		rec, rerr := rr.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return nil, rerr
		}
		if !mpileupKeepRecord(rec, opts) {
			continue
		}
		out = append(out, rec)
	}
	sort.SliceStable(out, func(i, j int) bool { return out[i].Pos < out[j].Pos })
	return out, nil
}
