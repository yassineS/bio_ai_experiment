// Package samtools — targetcut subcommand.
//
// Targetcut is a faithful Go port of upstream samtools'
// `cut_target.c` (Heng Li, 2011). The upstream tool consumes a
// SAM/BAM stream, builds a per-position pileup, calls a per-position
// consensus base + quality with the MAQ error model
// (`errmod_cal`), then runs a 2-state Viterbi HMM over the per-chrom
// consensus track to segment "covered, callable" regions away from
// "uncovered or low-information" regions. One consensus SAM record
// is emitted per identified region. The tool was originally designed
// for fosmid-pool sequencing where each region is a separately
// captured insert; in practice it works on any aligned BAM.
//
// Output: one tab-separated SAM line per region, with shape
//
//	<chrom>:<start>-<end>\t0\t<chrom>\t<start>\t60\t<len>M\t*\t0\t0\t<seq>\t<qual>
//
// where <seq> is the per-position consensus base ('A'/'C'/'G'/'T' or
// 'N' when no read provided usable evidence at that position) and
// <qual> is the per-position consensus quality (Phred+33-encoded;
// the field upstream right-shifts the packed score by >>2 before
// adding 33).
//
// Flags (matching upstream cut_target.c getopt):
//
//	-Q INT  Per-base quality cutoff (default 13).
//	-i INT  HMM entry penalty (the 0->1 transition penalty; default 14000).
//	-0 INT  HMM emission score for "no information at this position"
//	        in the "inside" state (default -4).
//	-1 INT  HMM emission score for "depth but no callable base" in
//	        "inside" state (default 1).
//	-2 INT  HMM emission score for "callable base" in "inside" state
//	        (default 6).
//	-f FILE Reference FASTA. When supplied, every per-record SEQ is
//	        run through pkg/htsgo/baq.SamProbRealn (apply+extend mode,
//	        flag = 1<<1|1) before its bases enter the pileup, matching
//	        upstream cut_target.c's read_aln.
//
// Compatibility note: prior to this port the Go re-implementation
// shipped a simpler "cut the aligned slice from each read to FASTA"
// mode. That mode is retained behind the SimpleMode option (CLI
// --simple) for backward compatibility; the default is now the
// upstream HMM consensus path.
package samtools

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// DefaultTargetcutMinBaseQ is the per-base quality cutoff used when -Q
// is not supplied (matches the upstream default of 13).
const DefaultTargetcutMinBaseQ uint8 = 13

// DefaultTargetcutEntryPenalty is the magnitude of the 0->1 transition
// penalty in the segmentation HMM. Upstream stores this as
// `g_param.p[0][1] = -atoi(optarg)` with a built-in default of 14000.
const DefaultTargetcutEntryPenalty = 14000

// Default emission scores for the "inside" state, indexed by the
// per-position class c (0 = no depth, 1 = depth but no callable base,
// 2 = callable base). These mirror the upstream default
// `g_param.e[1] = {-4,1,6}` literal.
const (
	DefaultTargetcutEmissionNoInfo   = -4
	DefaultTargetcutEmissionDepth    = 1
	DefaultTargetcutEmissionCallable = 6
)

// errModDepCorr is the depcorr value passed to errmod_init upstream:
// `errmod_init(1. - ERR_DEP)` with ERR_DEP = 0.83.
const errModDepCorr = 1.0 - 0.83

// TargetcutOptions configures the targetcut walker. The zero value
// matches the upstream defaults exactly.
type TargetcutOptions struct {
	// MinBaseQ drops bases with Phred quality strictly below this from
	// per-position consensus calling (upstream `-Q`, default 13).
	MinBaseQ uint8
	// EntryPenalty is the magnitude of the 0->1 transition penalty in
	// the segmentation HMM (upstream `-i`, default 14000). The value
	// is stored as a positive int and negated when applied.
	EntryPenalty int
	// EmissionNoInfo / EmissionDepth / EmissionCallable are the three
	// emission scores for the "inside" state. Upstream defaults:
	// {-4, 1, 6}. The "outside" state emissions are always {0,0,0}.
	EmissionNoInfo   int
	EmissionDepth    int
	EmissionCallable int

	// SimpleMode reverts to the legacy aligned-slice FASTA mode shipped
	// in v1 of this re-implementation. When false (the default), the
	// faithful upstream HMM-consensus path runs.
	SimpleMode bool

	// FastaRef is the optional reference FASTA path (upstream `-f`).
	// When non-empty, every record surviving the read filter is run
	// through pkg/htsgo/baq.SamProbRealn before its bases enter the
	// per-position pileup, matching upstream cut_target.c's read_aln
	// behaviour (`sam_prob_realn(b, g->ref, g->len, 1<<1|1)`). Has no
	// effect in SimpleMode.
	FastaRef string
}

// defaultedOptions returns opts with zero-valued fields filled in from
// the upstream defaults. This makes the zero TargetcutOptions{} behave
// like upstream's default invocation.
func defaultedOptions(opts TargetcutOptions) TargetcutOptions {
	if opts.MinBaseQ == 0 {
		opts.MinBaseQ = DefaultTargetcutMinBaseQ
	}
	if opts.EntryPenalty == 0 {
		opts.EntryPenalty = DefaultTargetcutEntryPenalty
	}
	if opts.EmissionNoInfo == 0 && opts.EmissionDepth == 0 && opts.EmissionCallable == 0 {
		opts.EmissionNoInfo = DefaultTargetcutEmissionNoInfo
		opts.EmissionDepth = DefaultTargetcutEmissionDepth
		opts.EmissionCallable = DefaultTargetcutEmissionCallable
	}
	return opts
}

// Targetcut reads SAM/BAM records from in and writes either one
// consensus SAM record per identified region (HMM mode, the default)
// or one FASTA entry per aligned read (SimpleMode). It returns the
// number of records emitted and the first error encountered.
func Targetcut(in io.Reader, out io.Writer, opts TargetcutOptions) (int, error) {
	if opts.SimpleMode {
		return targetcutSimple(in, out, opts)
	}
	return targetcutHMM(in, out, defaultedOptions(opts))
}

// ----------------------------------------------------------------------
// HMM consensus path — the faithful port of cut_target.c.
// ----------------------------------------------------------------------

// targetcutHMM streams the input once, buckets records by chromosome,
// then, in @SQ-header order, runs gencns + the 2-state Viterbi on each
// per-chrom consensus track and emits one SAM record per region.
func targetcutHMM(in io.Reader, out io.Writer, opts TargetcutOptions) (int, error) {
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, fmt.Errorf("samtools targetcut: open input: %w", err)
	}
	hdr := r.Header()

	// Bucket records by chrom (upstream pileup is per-chrom; records on
	// chroms not present in the header are silently dropped, matching
	// upstream's `tid >= 0` guard).
	type bucket struct {
		chrom string
		recs  []*sam.Record
	}
	chromIdx := map[string]int{}
	buckets := make([]bucket, 0, len(hdr.Refs))
	for _, ref := range hdr.Refs {
		chromIdx[ref.Name] = len(buckets)
		buckets = append(buckets, bucket{chrom: ref.Name})
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("samtools targetcut: read: %w", err)
		}
		if skipForTargetcutHMM(rec) {
			continue
		}
		idx, ok := chromIdx[rec.RName]
		if !ok {
			continue
		}
		buckets[idx].recs = append(buckets[idx].recs, rec)
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	em := errmod.Init(errModDepCorr)

	// Optional reference FASTA: when -f is supplied we run BAQ
	// realignment on every record before it enters the per-chrom
	// pileup, matching upstream cut_target.c::read_aln. The reference
	// is loaded lazily per chromosome (mirroring upstream's
	// per-tid fai_fetch64 cache).
	var ref *fasta.RandomAccess
	if opts.FastaRef != "" {
		ref, err = fasta.OpenRandomAccess(opts.FastaRef)
		if err != nil {
			return 0, fmt.Errorf("samtools targetcut: open reference %s: %w", opts.FastaRef, err)
		}
		defer ref.Close()
	}

	emitted := 0
	for i, b := range buckets {
		if len(b.recs) == 0 {
			continue
		}
		refLen := int(hdr.Refs[i].Length)
		if refLen <= 0 {
			continue
		}
		if ref != nil {
			refSeq, err := ref.Fetch(b.chrom, 0, int64(refLen))
			if err != nil {
				// Upstream tolerates a missing contig (fai_fetch64 may
				// return NULL) and simply skips BAQ for that chrom; we
				// do the same.
				refSeq = nil
			}
			if refSeq != nil {
				applyTargetcutBAQ(b.recs, refSeq)
			}
		}
		cns := buildConsensusTrack(b.recs, refLen, em, opts.MinBaseQ)
		n, err := emitTargetcutRegions(bw, b.chrom, cns, opts)
		if err != nil {
			return emitted, err
		}
		emitted += n
	}
	return emitted, nil
}

// applyTargetcutBAQ runs pkg/htsgo/baq.SamProbRealn on each record in
// apply+extend mode, mirroring upstream cut_target.c::read_aln which
// invokes `sam_prob_realn(b, g->ref, g->len, 1<<1|1)` once per record
// when a `-f` reference is supplied. The call mutates rec.Qual and
// adds a ZQ aux tag; we ignore SamProbRealn's return value because
// "no work needed" results (-1, -3) are not errors for our purposes
// — they leave the record's qualities untouched, exactly as upstream
// does.
func applyTargetcutBAQ(recs []*sam.Record, ref []byte) {
	flag := baq.FlagApply | baq.FlagExtend
	for _, rec := range recs {
		_ = baq.SamProbRealn(rec, ref, flag)
	}
}

// skipForTargetcutHMM mirrors the read filter inside cut_target.c::read_aln:
// drop unmapped / secondary / qcfail / dup. (Supplementary records pass
// upstream — they share the same flag treatment as primary.)
func skipForTargetcutHMM(rec *sam.Record) bool {
	if rec.Flag&(sam.FlagUnmapped|sam.FlagSecondary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
		return true
	}
	if rec.Pos <= 0 || rec.RName == "" || rec.RName == "*" {
		return true
	}
	return false
}

// consensusCell is the packed (cns_qual, base, depth) tuple gencns
// returns upstream. We store the upstream uint16 verbatim so the
// downstream Viterbi can use the exact same right-shift logic.
//
//	bits  0..7  : depth (capped at 255)
//	bits  8..9  : base ('A'/'C'/'G'/'T' = 0..3)
//	bits 10..15 : consensus quality, capped at 63
//
// A cell of 0 means "no read covered this position with usable info".
type consensusCell uint16

func (c consensusCell) depth() int    { return int(c) & 0xff }
func (c consensusCell) base() byte    { return byte((c >> 8) & 0x3) }
func (c consensusCell) cnsQual() byte { return byte((c >> 10) & 0x3f) }

// hasCallableBase reports whether the cell encodes a callable base
// (i.e. upstream's `cns[i]>>8 != 0` — the high byte is non-zero).
// Upstream uses this to distinguish "depth but no usable base" (c=1)
// from "callable base" (c=2). Note that when cnsQual is 0 AND base is
// 0 ('A') the high byte is 0 too — this is the exact ambiguity in
// upstream, preserved here.
func (c consensusCell) hasCallableBase() bool { return (c >> 8) != 0 }

// buildConsensusTrack walks every reference position 0..refLen and
// calls gencns over the per-position pileup. Returns a refLen-length
// vector of consensusCells.
func buildConsensusTrack(recs []*sam.Record, refLen int, em *errmod.Errmod, minBQ uint8) []consensusCell {
	// Build per-position contributions in a single CIGAR pass per read.
	// Each contribution is the packed `(q<<5 | strand<<4 | base)` byte
	// upstream feeds into errmod.
	contribs := make([][]uint16, refLen)
	for _, rec := range recs {
		appendTargetcutContribs(rec, contribs, refLen, minBQ)
	}

	out := make([]consensusCell, refLen)
	// Scratch buffers reused across positions; errmod takes a uint16
	// vector + writes a 4x4 float likelihood matrix.
	q := make([]float32, 16)
	for pos := 0; pos < refLen; pos++ {
		bases := contribs[pos]
		if len(bases) == 0 {
			continue
		}
		out[pos] = gencns(em, bases, q)
		// Free the per-position slice eagerly; targetcut is memory-bound
		// on large chroms and this keeps live working set close to one
		// chrom's worth of cells + the current position's contributions.
		contribs[pos] = nil
	}
	return out
}

// appendTargetcutContribs walks rec's CIGAR and appends one packed
// uint16 (q<<5 | strand<<4 | base) to contribs[pos] for every
// reference-aligned base that survives the upstream filters:
//   - skip is_refskip / is_del (D and N are not bases)
//   - skip baseQ < minBQ
//   - skip non-ACGT (`b > 3`)
//
// q is min(baseQ, mapQ) clamped to [4, 63].
func appendTargetcutContribs(rec *sam.Record, contribs [][]uint16, refLen int, minBQ uint8) {
	refPos := int(rec.Pos) - 1
	queryPos := 0
	isReverse := rec.Flag&sam.FlagReverse != 0
	hasQual := len(rec.Qual) == len(rec.Seq)
	mapQ := int(rec.MapQ)
	for _, op := range rec.Cigar {
		l := int(op.Length())
		o := op.Op()
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for k := 0; k < l; k++ {
				p := refPos + k
				q := queryPos + k
				if p < 0 || p >= refLen {
					continue
				}
				if q >= len(rec.Seq) {
					continue
				}
				b := nucBase4(rec.Seq[q])
				if b > 3 {
					continue
				}
				var baseQ int
				if hasQual {
					baseQ = int(rec.Qual[q])
				} else {
					// Upstream uses bam_get_qual which returns 0xff per
					// element when QUAL is unset; the >= minBQ check
					// then passes for any reasonable minBQ, so we mirror
					// that with a sentinel "high enough" value.
					baseQ = 0xff
				}
				if baseQ < int(minBQ) {
					continue
				}
				qq := baseQ
				if mapQ < qq {
					qq = mapQ
				}
				if qq < 4 {
					qq = 4
				}
				if qq > 63 {
					qq = 63
				}
				var strandBit uint16
				if isReverse {
					strandBit = 1
				}
				packed := uint16(qq)<<5 | strandBit<<4 | uint16(b)
				contribs[p] = append(contribs[p], packed)
			}
			refPos += l
			queryPos += l
		case sam.CigarInsertion, sam.CigarSoftClip:
			queryPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// No-op.
		}
	}
}

// nucBase4 maps an ASCII IUPAC base character to its 2-bit ACGT code.
// Anything outside ACGT (case-insensitive) returns 4 to signal "skip".
// Mirrors htslib's `seq_nt16_int[seq_nt16_table[base]]` chain.
func nucBase4(b byte) byte {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	}
	return 4
}

// gencns calls the consensus base at one pileup column. It is the
// upstream `gencns` function in Go: feeds the packed bases into the
// MAQ error model, picks the best base by score, and packs the result
// in the upstream uint16 shape.
func gencns(em *errmod.Errmod, bases []uint16, q []float32) consensusCell {
	em.Cal(bases, 4, q)
	// sum[i] holds (q[i*4+i] rounded) << 2 | i  — encoding base index
	// in the low 2 bits so the upcoming sort keeps them together.
	var sum [4]int
	for i := 0; i < 4; i++ {
		sum[i] = int(q[i*4+i]+0.499)<<2 | i
	}
	// Insertion sort, smallest first (upstream's identical sort).
	for i := 1; i < 4; i++ {
		for j := i; j > 0 && sum[j] < sum[j-1]; j-- {
			sum[j], sum[j-1] = sum[j-1], sum[j]
		}
	}
	qual := (sum[1] >> 2) - (sum[0] >> 2)
	if qual > 63 {
		qual = 63
	}
	if qual < 0 {
		qual = 0
	}
	depth := len(bases)
	if depth > 255 {
		depth = 255
	}
	best := sum[0] & 3
	// ret = (qual << 2 | base) << 8 | depth
	ret := uint16(qual)<<10 | uint16(best)<<8 | uint16(depth)
	return consensusCell(ret)
}

// emitTargetcutRegions runs the 2-state Viterbi over cns and emits one
// SAM-format consensus record per "inside" run. Returns the number of
// records emitted.
func emitTargetcutRegions(bw *bufio.Writer, chrom string, cns []consensusCell, opts TargetcutOptions) (int, error) {
	l := len(cns)
	if l == 0 {
		return 0, nil
	}

	e0 := [3]int{0, 0, 0}
	e1 := [3]int{opts.EmissionNoInfo, opts.EmissionDepth, opts.EmissionCallable}
	// p[from][to]; upstream: p={{0,-entryPenalty},{0,0}}
	pTransition := [2][2]int{
		{0, -opts.EntryPenalty},
		{0, 0},
	}

	b := make([]byte, l)
	var prev0, prev1 int // running f[prev][0], f[prev][1]
	for i := 0; i < l; i++ {
		c := classifyCell(cns[i])
		// curr[0] := max(prev[0]+e0[c]+p[0][0], prev[1]+e0[c]+p[1][0])
		tmp0 := prev0 + e0[c] + pTransition[0][0]
		tmp1 := prev1 + e0[c] + pTransition[1][0]
		var curr0, curr1 int
		if tmp0 > tmp1 {
			curr0 = tmp0
			b[i] = 0
		} else {
			curr0 = tmp1
			b[i] = 1
		}
		// curr[1] := max(prev[0]+e1[c]+p[0][1], prev[1]+e1[c]+p[1][1])
		tmp0 = prev0 + e1[c] + pTransition[0][1]
		tmp1 = prev1 + e1[c] + pTransition[1][1]
		if tmp0 > tmp1 {
			curr1 = tmp0
			// b[i] bit 1 = 0 (predecessor of state 1 was state 0)
		} else {
			curr1 = tmp1
			b[i] |= 1 << 1
		}
		prev0, prev1 = curr0, curr1
	}
	// Backtrack: bit 2 records the final state per position.
	var s byte
	if prev0 > prev1 {
		s = 0
	} else {
		s = 1
	}
	for i := l - 1; i > 0; i-- {
		b[i] |= s << 2
		s = (b[i] >> s) & 1
	}
	// Note: upstream cut_target.c's backtrack loop ranges (l-1, 0],
	// so b[0]'s bit-2 (state) flag is left at zero. We replicate that
	// quirk verbatim — position 0 never participates in a region,
	// even when the rest of the chrom is in state 1. This is byte-
	// equivalent to upstream `process_cns` behaviour.

	// Emit a SAM line per contiguous run of state-1 positions.
	emitted := 0
	start := -1 // 0-based start of the current run, -1 = no run
	// Upstream iterates i in [0, l]; at i==l we flush a trailing run.
	for i := 0; i <= l; i++ {
		var inside bool
		if i < l {
			inside = ((b[i] >> 2) & 3) != 0
		}
		switch {
		case !inside && start >= 0:
			if err := writeTargetcutRegion(bw, chrom, start, i, cns); err != nil {
				return emitted, err
			}
			emitted++
			start = -1
		case inside && start < 0:
			start = i
		}
	}
	return emitted, nil
}

// classifyCell returns the upstream `c` integer for the HMM emission
// table: 0 = empty cell, 1 = depth but no callable base, 2 = callable
// base. Mirrors upstream `(cns[i] == 0)? 0 : (cns[i]>>8 == 0)? 1 : 2`.
func classifyCell(cell consensusCell) int {
	if cell == 0 {
		return 0
	}
	if (cell >> 8) == 0 {
		return 1
	}
	return 2
}

// writeTargetcutRegion emits one consensus SAM record covering the
// half-open region [start0, end0) on chrom, mirroring upstream's
// process_cns printf format exactly:
//
//	%s:%d-%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t<seq>\t<qual>\n
//
// Seq/qual are derived from the per-position consensus cells: 'N' for
// empty cells, "ACGT"[base] otherwise; qual is `(cell>>8)>>2 + 33`,
// matching upstream's `putchar(33 + (cns[j]>>8>>2))`.
func writeTargetcutRegion(bw *bufio.Writer, chrom string, start0, end0 int, cns []consensusCell) error {
	s1 := start0 + 1
	e1 := end0
	width := end0 - start0
	// QNAME = "<chrom>:<s1>-<e1>"
	if _, err := bw.WriteString(chrom); err != nil {
		return err
	}
	if err := bw.WriteByte(':'); err != nil {
		return err
	}
	if _, err := bw.WriteString(strconv.Itoa(s1)); err != nil {
		return err
	}
	if err := bw.WriteByte('-'); err != nil {
		return err
	}
	if _, err := bw.WriteString(strconv.Itoa(e1)); err != nil {
		return err
	}
	// "\t0\t<chrom>\t<s1>\t60\t<width>M\t*\t0\t0\t"
	if _, err := bw.WriteString("\t0\t"); err != nil {
		return err
	}
	if _, err := bw.WriteString(chrom); err != nil {
		return err
	}
	if err := bw.WriteByte('\t'); err != nil {
		return err
	}
	if _, err := bw.WriteString(strconv.Itoa(s1)); err != nil {
		return err
	}
	if _, err := bw.WriteString("\t60\t"); err != nil {
		return err
	}
	if _, err := bw.WriteString(strconv.Itoa(width)); err != nil {
		return err
	}
	if _, err := bw.WriteString("M\t*\t0\t0\t"); err != nil {
		return err
	}
	// SEQ
	for j := start0; j < end0; j++ {
		if cns[j].hasCallableBase() {
			if err := bw.WriteByte("ACGT"[cns[j].base()]); err != nil {
				return err
			}
		} else {
			if err := bw.WriteByte('N'); err != nil {
				return err
			}
		}
	}
	if err := bw.WriteByte('\t'); err != nil {
		return err
	}
	// QUAL — upstream writes `33 + (cns[j]>>8>>2)`. cnsQual() does the
	// `>>2` shift on the already-shifted high byte.
	for j := start0; j < end0; j++ {
		if err := bw.WriteByte(33 + cns[j].cnsQual()); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}

// ----------------------------------------------------------------------
// Simple mode — the legacy aligned-slice FASTA cut, kept for users who
// were relying on the v1 behaviour. Opt-in via TargetcutOptions.SimpleMode
// (CLI --simple).
// ----------------------------------------------------------------------

// targetcutSimple is the original aligned-slice FASTA implementation,
// preserved here for backward compatibility. See the package doc for
// why the default switched to the HMM path.
func targetcutSimple(in io.Reader, out io.Writer, opts TargetcutOptions) (int, error) {
	if opts.MinBaseQ == 0 {
		opts.MinBaseQ = DefaultTargetcutMinBaseQ
	}
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, fmt.Errorf("samtools targetcut: open input: %w", err)
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	emitted := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return emitted, fmt.Errorf("samtools targetcut: read: %w", err)
		}
		if skipForTargetcutSimple(rec) {
			continue
		}
		seq := cutSequence(rec, opts.MinBaseQ)
		if len(seq) == 0 {
			continue
		}
		if _, err := bw.WriteString(">"); err != nil {
			return emitted, err
		}
		if _, err := bw.WriteString(rec.QName); err != nil {
			return emitted, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return emitted, err
		}
		if _, err := bw.Write(seq); err != nil {
			return emitted, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return emitted, err
		}
		emitted++
	}
	return emitted, nil
}

// skipForTargetcutSimple returns true for records that the simple FASTA
// mode should ignore — same filter as targetcutHMM plus supplementary,
// SEQ="*", and empty CIGAR.
func skipForTargetcutSimple(rec *sam.Record) bool {
	if rec.IsUnmapped() {
		return true
	}
	if rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
		return true
	}
	if rec.Seq == "" || rec.Seq == "*" {
		return true
	}
	if len(rec.Cigar) == 0 {
		return true
	}
	return false
}

// cutSequence returns the slice of rec.Seq corresponding to the
// aligned portion of the read. Soft-clip / hard-clip bases at either
// end are excluded; insertions are included (they are query bases
// that do not advance the reference); deletions / refskip / padding
// contribute nothing because they consume no query bases.
//
// Bases whose Phred quality is below minBaseQ are omitted. A QUAL
// field of "*" (length zero) is treated as "quality unknown" and no
// base is filtered.
func cutSequence(rec *sam.Record, minBaseQ uint8) []byte {
	qual := rec.Qual
	hasQual := len(qual) == len(rec.Seq)
	out := make([]byte, 0, len(rec.Seq))
	qpos := 0
	for _, op := range rec.Cigar {
		o := op.Op()
		n := int(op.Length())
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarInsertion:
			for k := 0; k < n; k++ {
				idx := qpos + k
				if idx >= len(rec.Seq) {
					break
				}
				if hasQual && qual[idx] < minBaseQ {
					continue
				}
				out = append(out, rec.Seq[idx])
			}
			qpos += n
		case sam.CigarSoftClip:
			qpos += n
		case sam.CigarDeletion, sam.CigarSkipped, sam.CigarPadding, sam.CigarHardClip:
			// No query bases consumed.
		default:
			// CigarBack and unknown ops: ignore.
		}
	}
	return out
}

// SortConsensusForTesting is exported for test packages that want to
// reproduce gencns's deterministic insertion-sort ordering without
// pulling in the unexported helpers.
func SortConsensusForTesting(s []int) {
	sort.Ints(s)
}
