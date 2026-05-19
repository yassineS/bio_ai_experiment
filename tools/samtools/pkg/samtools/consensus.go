package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// ConsensusFormat selects the output format for samtools consensus.
type ConsensusFormat int

const (
	// ConsensusFASTA emits one FASTA record per contig.
	ConsensusFASTA ConsensusFormat = iota
	// ConsensusFASTQ emits one FASTQ record per contig; quality bytes
	// carry the per-position call confidence (Phred+33, capped at 93).
	ConsensusFASTQ
	// ConsensusPileup emits the samtools consensus pileup row format:
	// "chrom\tpos\tnth\tdepth\tcall\tcq\tseq\tqual\n". nth is always
	// 0 in v1; upstream emits insertion-column rows with nth>0 when
	// --show-ins is yes (the default), and this is tracked as a
	// deferred item in docs/PARITY_ROADMAP.md#samtools.
	ConsensusPileup
)

// ConsensusMode selects the per-position calling algorithm.
type ConsensusMode int

const (
	// ConsensusModeSimple uses the per-position majority-vote / qual-sum
	// scoring identical to upstream `calculate_consensus_simple`.
	ConsensusModeSimple ConsensusMode = iota
	// ConsensusModeBayesian is upstream's Gap5-derived posterior caller.
	// Accepted but not implemented in v1; ConsensusFile emits a stderr
	// warning and falls back to ConsensusModeSimple. Upstream's default
	// mode is MODE_RECALL (a flavour of bayesian), so default invocation
	// of our binary takes this fallback path.
	ConsensusModeBayesian
)

// DefaultConsensusLineLen is the upstream samtools consensus line wrap
// width (`-l/--line-len`) for FASTA / FASTQ output.
const DefaultConsensusLineLen = 70

// DefaultConsensusMinBaseQ is the upstream samtools consensus `--min-BQ`
// default (0 — distinct from mpileup's 13).
const DefaultConsensusMinBaseQ = 0

// ConsensusOptions configures samtools consensus.
//
// Zero-valued defaults are filled in by applyConsensusDefaults to match
// upstream `bam_consensus.c::main_consensus`'s `consensus_opts` struct
// initialisers (`call_fract=0.75`, `het_fract=0.5`, `min_depth=1`,
// `line_len=70`, `show_ins=1`, `show_del=0`, `use_qual=0`, `ambig=0`).
type ConsensusOptions struct {
	// Input is the path to a BAM/SAM file. When empty, ConsensusReader
	// is used instead.
	Input string
	// Format selects FASTA / FASTQ / Pileup output (CLI -f/--format).
	Format ConsensusFormat
	// Mode picks the per-position algorithm (CLI -m/--mode). Upstream's
	// default is MODE_RECALL (a bayesian mode); ConsensusFile/Consensus
	// honour ConsensusModeBayesian by emitting a stderr warning and
	// falling back to ConsensusModeSimple, since v1 only implements
	// simple mode.
	Mode ConsensusMode
	// AllPositions emits zero-coverage positions as 'N' (CLI -a).
	// When false, contigs with no covered positions emit nothing.
	AllPositions bool
	// Regions is a list of "chr[:start-end]" specifiers (CLI -r).
	Regions []string
	// BEDPath is a BED file restricting emission. v1 reuses this for the
	// CLI's `-l/--positions` path; upstream doesn't expose -l on
	// consensus, but this is a useful test-driver hook and was on the
	// closed prior PR's wire.
	BEDPath string
	// MinDepth filters out positions with strictly less than this many
	// usable reads (CLI -d/--min-depth, default 1).
	MinDepth int
	// MinCallFraction is the only fraction gate. Upstream
	// (bam_consensus.c:1988-1994) calls a position 'N' (or '*' if call1
	// was the gap bucket) when usedScore < MinCallFraction * tscore.
	// The CLI flag is -c/--call-fract (default 0.75).
	MinCallFraction float64
	// MinHetFraction is the minimum proportion of the second-best score
	// over the best required to call a heterozygous IUPAC ambiguity
	// (CLI -H/--het-fract, default 0.5). Only meaningful with
	// AmbigCodes=true. Upstream: bam_consensus.c:1982.
	MinHetFraction float64
	// MinMAPQ skips reads with MAPQ below this value (CLI --min-MQ).
	// Note: there is NO short alias for this flag upstream — `-q`
	// upstream is `--use-qual`, not min MAPQ.
	MinMAPQ uint8
	// MinBaseQ skips bases with quality below this value (CLI --min-BQ).
	MinBaseQ uint8
	// AmbigCodes enables IUPAC ambiguity codes for heterozygous calls
	// (CLI -A/--ambig).
	AmbigCodes bool
	// UseQual, when true, weights each base contribution by its
	// per-base quality (the `q` multiplier in upstream's score sums).
	// When false (the upstream default), each base contributes the
	// unit weight from `seqi2{A,C,G,T}` — equivalent to a pure
	// frequency count. CLI: -q/--use-qual / --no-use-qual.
	UseQual bool
	// LineLen wraps FASTA / FASTQ lines after this many bases (CLI -l).
	// 0 means use DefaultConsensusLineLen.
	LineLen int
	// ShowDel renders deletion placeholder '*' in the consensus
	// sequence. Upstream default is `0` (no '*'). Pileup mode also
	// honours this — when ShowDel is false, pileup rows whose call is
	// '*' are suppressed (matching bam_consensus.c:2244).
	ShowDel bool
	// NoShowIns suppresses inserted bases in FASTA/FASTQ output. The
	// upstream default is "show insertions" (show_ins=1); the zero
	// value of NoShowIns (false) reproduces it.
	NoShowIns bool
	// MarkIns prepends '+' to every inserted base/qual emitted in
	// FASTA/FASTQ — upstream's --mark-ins. v1 wires the option but
	// only implements it when NoShowIns=false.
	MarkIns bool
	// HetOnly, when true, suppresses calls that are not heterozygous
	// (upstream --het-only). v1 accepts it but does not implement the
	// filter; tracked in docs/PARITY_ROADMAP.md.
	HetOnly bool
	// IgnoreOverlaps is accepted for CLI compatibility but not
	// implemented (v1 does not deduplicate mate-pair overlaps in the
	// consensus walker).
	IgnoreOverlaps bool
	// CountOrphans accepts reads with unmapped mates / anomalous pair
	// flags. Upstream consensus' default `excl_flags` only drops
	// UNMAP|SECONDARY|QCFAIL|DUP — paired-but-not-proper reads are
	// kept. We mirror that by always running with CountOrphans=true.
	CountOrphans bool
	// Output is the path to write to (default stdout). CLI -o.
	Output string
	// Threads is accepted but ignored in v1.
	Threads int
}

// ConsensusFile is the file-path entry point for `samtools consensus`.
// Opens the input file and delegates to Consensus. errOut, when
// non-nil, receives the bayesian-fallback stderr warning.
func ConsensusFile(opts ConsensusOptions, out io.Writer, errOut io.Writer) error {
	if opts.Input == "" {
		return fmt.Errorf("samtools consensus: no input file")
	}
	if opts.Mode == ConsensusModeBayesian {
		if errOut != nil {
			fmt.Fprintln(errOut, "samtools consensus: --mode bayesian not yet implemented; falling back to simple (tracked in docs/PARITY_ROADMAP.md#samtools)")
		}
		opts.Mode = ConsensusModeSimple
	}
	f, err := os.Open(opts.Input)
	if err != nil {
		return fmt.Errorf("samtools consensus: %w", err)
	}
	defer f.Close()
	return Consensus(f, out, opts)
}

// Consensus is the streaming entry point: read records from `in` and
// emit the consensus per the configured format to `out`.
//
// Note: Consensus does not emit the bayesian-fallback warning by itself
// (it has no stderr handle). The CLI path goes through ConsensusFile,
// which routes the warning to os.Stderr; library callers that need the
// warning should normalise opts.Mode themselves before invoking
// Consensus.
func Consensus(in io.Reader, out io.Writer, opts ConsensusOptions) error {
	applyConsensusDefaults(&opts)
	// Library-level safety: if a caller passes Mode=Bayesian directly
	// to Consensus (not through ConsensusFile), still fall back to
	// simple so we don't pretend to support a mode we don't implement.
	if opts.Mode == ConsensusModeBayesian {
		opts.Mode = ConsensusModeSimple
	}
	rd, err := sam.NewReader(in)
	if err != nil {
		return fmt.Errorf("samtools consensus: %w", err)
	}
	hdr := rd.Header()

	// Bucket records by chromosome up front. We reuse the same filter
	// strategy as the mpileup walker (drop unmapped/secondary/etc,
	// apply MAPQ floor), with the consensus-specific tweak that we
	// always accept orphan mates (matching upstream consensus, which
	// only excludes UNMAP|SECONDARY|QCFAIL|DUP).
	mpopts := MpileupOptions{
		MinMAPQ:      opts.MinMAPQ,
		MinBaseQ:     opts.MinBaseQ,
		CountOrphans: true,
	}
	byChrom, err := bucketByChrom(rd, mpopts, hdr)
	if err != nil {
		return fmt.Errorf("samtools consensus: %w", err)
	}

	// Resolve regions.
	resolved, _, err := ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if err != nil {
		return err
	}
	regionByChrom := map[string][][2]int{}
	for _, r := range resolved {
		end0 := r.End0
		if end0 > 1<<29 {
			end0 = int(refLengthForName(hdr, r.Region.Chrom))
		}
		regionByChrom[r.Region.Chrom] = append(regionByChrom[r.Region.Chrom], [2]int{r.Beg0, end0})
	}

	// Optional BED filter (-l/--positions on our CLI; upstream consensus
	// has no equivalent, this is a v1 testing convenience).
	var posFilter *positionFilter
	if opts.BEDPath != "" {
		pf, perr := loadPositionsFile(opts.BEDPath)
		if perr != nil {
			return fmt.Errorf("samtools consensus: %w", perr)
		}
		posFilter = pf
	}

	// Decide chromosome walk order.
	var chromsToWalk []string
	switch {
	case len(opts.Regions) > 0:
		seen := map[string]struct{}{}
		for _, r := range resolved {
			if _, ok := seen[r.Region.Chrom]; ok {
				continue
			}
			seen[r.Region.Chrom] = struct{}{}
			chromsToWalk = append(chromsToWalk, r.Region.Chrom)
		}
	case opts.AllPositions:
		// Walk every contig in the header so empty contigs emit N's.
		for _, ref := range hdr.Refs {
			chromsToWalk = append(chromsToWalk, ref.Name)
		}
	default:
		// Only contigs touched by at least one record, in header order.
		hit := map[string]struct{}{}
		for chrom := range byChrom {
			hit[chrom] = struct{}{}
		}
		for _, ref := range hdr.Refs {
			if _, ok := hit[ref.Name]; ok {
				chromsToWalk = append(chromsToWalk, ref.Name)
			}
		}
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	for _, chrom := range chromsToWalk {
		recs := byChrom[chrom]
		refLen := int(refLengthForName(hdr, chrom))
		if refLen <= 0 {
			continue
		}

		var windows [][2]int
		if iv, ok := regionByChrom[chrom]; ok {
			windows = mergeIntervals(iv)
		} else {
			windows = [][2]int{{0, refLen}}
		}

		if err := emitConsensusContig(bw, chrom, refLen, windows, recs, posFilter, opts); err != nil {
			return err
		}
	}
	return nil
}

// applyConsensusDefaults fills zero-valued options with their CLI
// defaults to match upstream `consensus_opts` initialisers (see
// bam_consensus.c:2981+).
//
//   - LineLen <= 0 -> 70.
//   - MinDepth <= 0 -> 1.
//   - MinCallFraction <= 0 -> 0.75.
//   - MinHetFraction <= 0 -> 0.5.
//
// We deliberately do NOT clobber opts.ShowDel for the pileup format:
// upstream honours --show-del in pileup too (bam_consensus.c:2244
// suppresses rows whose call is '*' when show_del is false), and the
// reviewer's correctness finding #5 requires we mirror that.
func applyConsensusDefaults(opts *ConsensusOptions) {
	if opts.LineLen <= 0 {
		opts.LineLen = DefaultConsensusLineLen
	}
	if opts.MinDepth <= 0 {
		opts.MinDepth = 1
	}
	if opts.MinCallFraction <= 0 {
		opts.MinCallFraction = 0.75
	}
	if opts.MinHetFraction <= 0 {
		opts.MinHetFraction = 0.5
	}
}

// consensusCall is a single per-position call.
type consensusCall struct {
	base byte // 'A','C','G','T','*','N' or IUPAC ambig code when AmbigCodes set
	// qual is the per-position confidence as an integer 0..100
	// (upstream's formula: `100 * used_score / tscore`). The pileup
	// format emits this number verbatim; FASTA/FASTQ encode it as a
	// Phred+33 byte capped at 93.
	qual  int
	depth int // usable depth (post-MinBaseQ)
}

// emitConsensusContig walks every window on chrom and emits per-format
// records (FASTA/FASTQ/Pileup).
func emitConsensusContig(bw *bufio.Writer, chrom string, refLen int, windows [][2]int,
	recs []*sam.Record, posFilter *positionFilter, opts ConsensusOptions) error {

	// For FASTA/FASTQ we accumulate one buffer per contig and emit at
	// the end. For pileup we stream line-by-line.
	var seqBuf, qualBuf []byte

	for _, w := range windows {
		beg0 := w[0]
		end0 := w[1]
		if end0 > refLen {
			end0 = refLen
		}
		if beg0 < 0 {
			beg0 = 0
		}
		if beg0 >= end0 {
			continue
		}
		// Build per-position event slices for this window by reusing
		// the mpileup accumulator. This guarantees the consensus is
		// byte-faithful with what `samtools mpileup` reports for the
		// same input.
		width := end0 - beg0
		events := make([][]pileupEvent, width)
		accOpts := MpileupOptions{
			MinBaseQ: opts.MinBaseQ,
			MinMAPQ:  opts.MinMAPQ,
		}
		for ridx, rec := range recs {
			accumulateRecordEvents(rec, ridx, beg0, end0, events, accOpts, nil, chrom)
		}

		// Walk positions.
		for pos0 := beg0; pos0 < end0; pos0++ {
			pos1 := pos0 + 1
			if posFilter != nil && !posFilter.contains(chrom, pos1) {
				continue
			}
			col := pos0 - beg0
			call, totalDepth := callConsensus(events[col], opts)

			switch opts.Format {
			case ConsensusPileup:
				if totalDepth == 0 && !opts.AllPositions {
					continue
				}
				// Honour --show-del in pileup mode too: when the
				// call is '*' and ShowDel is false, suppress the
				// row, matching bam_consensus.c:2244.
				if call.base == '*' && !opts.ShowDel {
					continue
				}
				if err := writeConsensusPileupRow(bw, chrom, pos1, totalDepth, call, events[col], opts); err != nil {
					return err
				}
			default:
				// FASTA / FASTQ accumulate.
				if call.base == 0 {
					// No coverage AND not -a -> skip.
					if !opts.AllPositions {
						continue
					}
					seqBuf = append(seqBuf, 'N')
					qualBuf = append(qualBuf, '!')
					continue
				}
				if call.base == '*' && !opts.ShowDel {
					// Suppress deletion placeholder.
					continue
				}
				seqBuf = append(seqBuf, call.base)
				qualBuf = append(qualBuf, phredByte(call.qual))
				if !opts.NoShowIns {
					if insSeq, insQuals := callConsensusInsertion(events[col], opts); len(insSeq) > 0 {
						if opts.MarkIns {
							// Upstream prepends a single '+' marker; we
							// keep that semantic and mirror it for qual.
							seqBuf = append(seqBuf, '+')
							qualBuf = append(qualBuf, '+')
						}
						seqBuf = append(seqBuf, insSeq...)
						qualBuf = append(qualBuf, insQuals...)
					}
				}
			}
		}
	}

	if opts.Format != ConsensusPileup {
		if len(seqBuf) == 0 && !opts.AllPositions {
			return nil
		}
		// AllPositions on an untouched contig: emit a record full of
		// N's, matching upstream's `-aa` semantics.
		if len(seqBuf) == 0 && opts.AllPositions {
			seqBuf = make([]byte, refLen)
			for i := range seqBuf {
				seqBuf[i] = 'N'
			}
			qualBuf = make([]byte, refLen)
			for i := range qualBuf {
				qualBuf[i] = '!'
			}
		}
		writeFastaFastqRecord(bw, chrom, seqBuf, qualBuf, opts)
	}
	return nil
}

// callConsensus runs the simple per-position consensus on the event
// slice. Returns a call (base=0 when there's no usable coverage AND
// the caller should treat as a skip) and the usable depth.
//
// The score model mirrors upstream's `calculate_consensus_simple`
// (bam_consensus.c:1900-2006) bit-for-bit:
//
//   - bases bucket into A/C/G/T/* (seqi codes 1,2,4,8,16);
//   - ambiguous IUPAC bases distribute weight across A/C/G/T per the
//     seqi2A/C/G/T tables;
//   - the score multiplier is the base quality when UseQual=true,
//     else 1 (frequency-only count) — upstream's `use_qual=0` default;
//   - call1/call2 are best/second-best by score;
//   - if AmbigCodes AND score2 >= het_fract*score1, the call is the
//     bitwise OR of call1|call2 (an IUPAC ambiguity);
//   - the single fraction gate is `used_score < call_fract * tscore`
//     (the reviewer's correctness finding #1: this is the ONLY
//     fraction gate upstream applies);
//   - depth gate: tot_depth < min_depth promotes to N (or '*' if
//     call1 was the gap bucket).
func callConsensus(evs []pileupEvent, opts ConsensusOptions) (consensusCall, int) {
	// Score buckets indexed by 4-bit "seqi" code:
	//   A=1, C=2, G=4, T=8, *=16 (out of band).
	var freq [17]int
	var score [17]uint64
	totDepth := 0
	for _, e := range evs {
		if e.dropped {
			continue
		}
		switch e.kind {
		case pileupEventBase:
			if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
				continue
			}
			b := upper(e.base)
			seqi := baseToSeqi(b)
			if seqi == 0 {
				// 'N' / unknown
				continue
			}
			// Frequency-only count by default (upstream use_qual=0).
			// With UseQual=true, weight by per-base quality. Upstream
			// bam_consensus.c:1937-1953 computes Q = wt * (use_qual?q:1)
			// and increments freq/score only when Q != 0; a qual-0 base
			// therefore contributes nothing under UseQual=true.
			var q uint64 = 1
			if opts.UseQual {
				q = uint64(e.qual)
				if q == 0 {
					continue
				}
			}
			// Map ambiguous IUPAC bases by upstream's seqi2A/C/G/T
			// matrices.
			for _, slot := range [4]struct {
				m   [16]int
				bit int
			}{
				{seqi2A, 1}, {seqi2C, 2}, {seqi2G, 4}, {seqi2T, 8},
			} {
				wt := slot.m[seqi&15]
				if wt == 0 {
					continue
				}
				freq[slot.bit]++
				score[slot.bit] += q * uint64(wt)
			}
			totDepth++
		case pileupEventDel:
			freq[16]++
			score[16] += 8
			totDepth++
		case pileupEventRefSkip:
			// Reference skips don't contribute to consensus.
		}
	}

	if totDepth == 0 {
		return consensusCall{}, 0
	}

	// Total score across A/C/G/T/*.
	var tscore uint64
	for _, c := range []int{1, 2, 4, 8, 16} {
		tscore += score[c]
	}

	// Best and second-best.
	call1, call2 := 15, 15
	var s1, s2 uint64
	for _, c := range []int{1, 2, 4, 8, 16} {
		if score[c] > s1 {
			s2 = s1
			call2 = call1
			s1 = score[c]
			call1 = c
		} else if score[c] > s2 {
			s2 = score[c]
			call2 = c
		}
	}

	used := call1
	usedScore := s1
	// Het condition mirrors upstream (bam_consensus.c:1982):
	//   `score2 >= het_fract * score1 && opts->ambig`.
	// No s1>0 guard: when s1 is 0, s2 is also 0 and the call ends up
	// at N anyway via the call_fract gate below — staying consistent
	// with upstream is cheaper than guarding.
	// Het condition mirrors upstream bam_consensus.c:1982 exactly:
	// score2 >= het_fract * score1 && ambig. No extra guards.
	if opts.AmbigCodes && float64(s2) >= opts.MinHetFraction*float64(s1) {
		used |= call2
		usedScore += s2
	}

	// Single fraction gate (upstream bam_consensus.c:1988-1994).
	depthOK := totDepth >= opts.MinDepth
	callOK := tscore > 0 && float64(usedScore) >= opts.MinCallFraction*float64(tscore)
	if !depthOK || !callOK {
		// "But note shallow gaps are still called gaps, not N, as
		//  we're still more confident there is no base than it is
		//  A, C, G or T." — upstream comment.
		if call1 == 16 {
			used = 16
		} else {
			used = 0 // 'N'
		}
	}

	// "N" + IUPAC + gap-row tables from upstream.
	const het = "NACMGRSVTWYHKDBN" + "*ac?g???t???????"
	base := het[used]
	if base == '?' {
		base = 'N'
	}

	// Confidence as a 0..100 integer (upstream's formula).
	qual := 0
	if used != 0 && tscore > 0 {
		qual = int(100 * float64(usedScore) / float64(tscore))
	}
	return consensusCall{base: base, qual: qual, depth: totDepth}, totDepth
}

// callConsensusInsertion derives a consensus for an inserted-sequence
// annotation attached to the current event slice. The simplest
// faithful behaviour (matching upstream's "show_ins yes" / simple
// mode) is: if a majority of reads carry an insertion at this position,
// emit its consensus sequence; otherwise nothing.
//
// The gate uses MinCallFraction (reviewer correctness finding #8) —
// there is no separate "min fraction" knob.
func callConsensusInsertion(evs []pileupEvent, opts ConsensusOptions) ([]byte, []byte) {
	withIns := 0
	insSeqs := [][]byte{}
	live := 0
	for _, e := range evs {
		if e.dropped || e.kind != pileupEventBase {
			continue
		}
		if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
			continue
		}
		live++
		if e.insAfter != "" {
			withIns++
			insSeqs = append(insSeqs, []byte(strings.ToUpper(e.insAfter)))
		}
	}
	if live == 0 || withIns == 0 {
		return nil, nil
	}
	if float64(withIns)/float64(live) < opts.MinCallFraction {
		return nil, nil
	}
	maxLen := 0
	for _, s := range insSeqs {
		if len(s) > maxLen {
			maxLen = len(s)
		}
	}
	if maxLen == 0 {
		return nil, nil
	}
	seq := make([]byte, maxLen)
	qual := make([]byte, maxLen)
	for c := 0; c < maxLen; c++ {
		var cnt [256]int
		for _, s := range insSeqs {
			if c >= len(s) {
				continue
			}
			cnt[s[c]]++
		}
		bestB := byte('N')
		bestC := 0
		total := 0
		for b, n := range cnt {
			if n > 0 {
				total += n
			}
			if n > bestC {
				bestC = n
				bestB = byte(b)
			}
		}
		if total == 0 || float64(bestC)/float64(total) < opts.MinCallFraction {
			seq[c] = 'N'
			qual[c] = '!'
			continue
		}
		seq[c] = bestB
		q := int(100 * float64(bestC) / float64(total))
		qual[c] = phredByte(q)
	}
	return seq, qual
}

// writeConsensusPileupRow emits one row of the samtools consensus
// pileup schema: chrom\tpos\tnth\tdepth\tcall\tcq\tseq\tqual\n.
//
// nth is always 0 in v1. Upstream's `--show-ins yes` default emits
// extra rows with nth>0 for the columns of each insertion; v1 does
// not yet emit those rows (deferred — docs/PARITY_ROADMAP.md).
func writeConsensusPileupRow(bw *bufio.Writer, chrom string, pos1, depth int,
	call consensusCall, evs []pileupEvent, opts ConsensusOptions) error {

	sortEvents(evs)
	var seq, qual []byte
	for _, e := range evs {
		if e.dropped {
			continue
		}
		var b byte
		var q byte
		switch e.kind {
		case pileupEventBase:
			if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
				continue
			}
			b = upper(e.base)
			if e.isReverse {
				b = lower(b)
			}
			q = e.qual
		case pileupEventDel:
			b = '*'
			if e.isReverse {
				b = '#'
			}
			q = e.qual
		case pileupEventRefSkip:
			continue
		}
		seq = append(seq, b)
		qual = append(qual, q+33)
	}

	cb := call.base
	if cb == 0 {
		cb = 'N'
	}
	if _, err := bw.WriteString(chrom); err != nil {
		return err
	}
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(pos1))
	bw.WriteByte('\t')
	bw.WriteByte('0') // nth
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(depth))
	bw.WriteByte('\t')
	bw.WriteByte(cb)
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(call.qual))
	bw.WriteByte('\t')
	if len(seq) == 0 {
		bw.WriteByte('*')
	} else {
		bw.Write(seq)
	}
	bw.WriteByte('\t')
	if len(qual) == 0 {
		bw.WriteByte('*')
	} else {
		bw.Write(qual)
	}
	bw.WriteByte('\n')
	return nil
}

// writeFastaFastqRecord emits one >/@ record with the accumulated seq
// and (for FASTQ) quality buffers wrapped at opts.LineLen.
func writeFastaFastqRecord(bw *bufio.Writer, chrom string, seq, qual []byte, opts ConsensusOptions) {
	if len(seq) == 0 {
		return
	}
	if opts.Format == ConsensusFASTQ {
		bw.WriteByte('@')
	} else {
		bw.WriteByte('>')
	}
	bw.WriteString(chrom)
	bw.WriteByte('\n')
	for i := 0; i < len(seq); i += opts.LineLen {
		end := i + opts.LineLen
		if end > len(seq) {
			end = len(seq)
		}
		bw.Write(seq[i:end])
		bw.WriteByte('\n')
	}
	if opts.Format == ConsensusFASTQ {
		bw.WriteString("+\n")
		for i := 0; i < len(qual); i += opts.LineLen {
			end := i + opts.LineLen
			if end > len(qual) {
				end = len(qual)
			}
			bw.Write(qual[i:end])
			bw.WriteByte('\n')
		}
	}
}

// phredByte converts a 0..100 Phred-ish confidence to a printable
// ASCII byte (Phred+33), clamped at 93 so output stays within
// printable ASCII (33..126).
func phredByte(q int) byte {
	if q < 0 {
		q = 0
	}
	if q > 93 {
		q = 93
	}
	return byte(q) + 33
}

// baseToSeqi maps an ASCII base byte to its BAM 4-bit "seqi" code.
//
//	'='->0 'A'->1 'C'->2 'M'->3 'G'->4 'R'->5 'S'->6 'V'->7
//	'T'->8 'W'->9 'Y'->10 'H'->11 'K'->12 'D'->13 'B'->14 'N'->15
//
// Unrecognised bases map to 15 ('N'), letting them be ignored as
// ambiguous (upstream behaviour).
func baseToSeqi(b byte) int {
	switch b {
	case '=':
		return 0
	case 'A':
		return 1
	case 'C':
		return 2
	case 'M':
		return 3
	case 'G':
		return 4
	case 'R':
		return 5
	case 'S':
		return 6
	case 'V':
		return 7
	case 'T', 'U':
		return 8
	case 'W':
		return 9
	case 'Y':
		return 10
	case 'H':
		return 11
	case 'K':
		return 12
	case 'D':
		return 13
	case 'B':
		return 14
	case 'N':
		return 15
	}
	return 15
}

// Pre-computed upstream mapping tables (seqi -> weight on a pure
// base). Mirrors `bam_consensus.c::calculate_consensus_simple::
// seqi2A/C/G/T`.
var (
	seqi2A = [16]int{0, 8, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 0, 2, 0, 1}
	seqi2C = [16]int{0, 0, 8, 4, 0, 0, 4, 2, 0, 0, 4, 2, 0, 0, 2, 1}
	seqi2G = [16]int{0, 0, 0, 0, 8, 4, 4, 1, 0, 0, 0, 0, 4, 2, 2, 1}
	seqi2T = [16]int{0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 4, 2, 8, 2, 2, 1}
)

// ParseConsensusFormat maps a CLI -f/--format value to a
// ConsensusFormat. The lookup is case-insensitive.
func ParseConsensusFormat(s string) (ConsensusFormat, error) {
	switch strings.ToLower(s) {
	case "", "fasta", "fa":
		return ConsensusFASTA, nil
	case "fastq", "fq":
		return ConsensusFASTQ, nil
	case "pileup":
		return ConsensusPileup, nil
	default:
		return 0, fmt.Errorf("samtools consensus: unknown format %q (want fasta|fastq|pileup)", s)
	}
}

// ParseConsensusMode maps a CLI -m/--mode value to a ConsensusMode.
// Upstream accepts "simple", "bayesian", "bayesian_r", "bayesian_m",
// "bayesian_p", "bayesian_116" — we collapse every bayesian variant
// onto ConsensusModeBayesian since v1 only implements simple.
func ParseConsensusMode(s string) (ConsensusMode, error) {
	switch strings.ToLower(s) {
	case "", "simple":
		return ConsensusModeSimple, nil
	case "bayesian", "bayesian_r", "bayesian_m", "bayesian_p", "bayesian_116":
		return ConsensusModeBayesian, nil
	default:
		return 0, fmt.Errorf("samtools consensus: unknown mode %q (want simple|bayesian)", s)
	}
}
