package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
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
	// "chrom\tpos\tnth\tdepth\tcall\tcq\tseq\tqual\n". The nth==0 row is
	// the reference column; insertion columns are emitted as additional
	// rows with nth>0 when --show-ins is yes (the default), matching
	// upstream for both simple and bayesian modes.
	ConsensusPileup
)

// ConsensusMode selects the per-position calling algorithm.
type ConsensusMode int

const (
	// ConsensusModeSimple uses the per-position majority-vote / qual-sum
	// scoring identical to upstream `calculate_consensus_simple`.
	ConsensusModeSimple ConsensusMode = iota
	// ConsensusModeBayesian is upstream's Gap5-derived posterior caller.
	// The specific bayesian flavour (BAYES_116 / RECALL / PRECISE / MIXED)
	// is selected by ConsensusOptions.BayesianSubMode. Upstream's default
	// mode is MODE_RECALL, so the bare `samtools consensus` invocation
	// runs ConsensusModeBayesian with the RECALL sub-mode.
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
	// default is MODE_RECALL (a bayesian mode). Both ConsensusModeSimple
	// and ConsensusModeBayesian (the Gap5 posterior caller, including its
	// RECALL/PRECISE/MIXED/116 sub-modes) are fully implemented; the bare
	// `samtools consensus` invocation runs ConsensusModeBayesian.
	Mode ConsensusMode
	// AllPositions emits zero-coverage positions as 'N' across the full
	// length of every contig that has at least one read (CLI -a). When
	// false, only the covered span of each touched contig is emitted.
	AllPositions bool
	// AllContigs additionally emits records/rows for contigs with no
	// reads at all (CLI -aa). Implies AllPositions.
	AllContigs bool
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
	// MarkIns prepends '_' to every inserted column emitted in
	// FASTA/FASTQ — upstream's --mark-ins (the marker byte is '_' in both
	// the seq and qual streams; bam_consensus.c:2409-2412). Only effective
	// when NoShowIns=false.
	MarkIns bool
	// HetOnly corresponds to upstream's --het-only flag (CLI --het-only).
	// When true, the consensus is restricted to HETEROZYGOUS-called
	// positions: homozygous and no-call positions are suppressed. In
	// FASTA/FASTQ they are rendered as 'N' (coordinates preserved); in
	// pileup the rows are omitted entirely. Het-ness is determined
	// independently of AmbigCodes (simple mode: score2 >= het_fract*score1
	// on a confident call; bayesian mode: positive het log-odds on a
	// confident call).
	//
	// NOTE: upstream samtools (through 1.22) parses --het-only into its
	// options struct but never reads it again, so the flag is a no-op
	// there (a dead-option bug). We implement the intended behaviour. See
	// docs/UPSTREAM_BUGS.md and docs/PARITY_ROADMAP.md.
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

	// Reference is the path to an indexed reference FASTA (-T/--reference).
	// When set, every no-coverage / gap position that would otherwise be
	// emitted as 'N' is filled from the reference instead, mirroring upstream
	// bam_consensus.c (update_ref + empty_pileup2 + the basic_fasta /
	// ref_or_Ns gap fills). The substitution applies ONLY to positions with
	// no usable coverage; genuinely-called columns (including low-depth 'N'
	// calls) keep their computed call. The .fai index must exist alongside the
	// FASTA (or be buildable). When empty, no reference is loaded and gaps
	// stay 'N'.
	Reference string
	// RefQual is the phred quality (0-93) assigned to reference-filled
	// positions in FASTQ output (--ref-qual, default 0). It is the offset
	// added to '!' for every gap base substituted from the reference; it has
	// no effect on FASTA output or on covered positions. Mirrors upstream
	// opts->ref_qual.
	RefQual int

	// --- Bayesian-mode knobs (only meaningful when Mode is
	// ConsensusModeBayesian). Zero values are filled in by
	// applyConsensusDefaults to match upstream main_consensus. ---

	// BayesianSubMode selects the bayesian parameter set. It is not set
	// directly by external callers; use SetBayesianMode with a CLI mode
	// string instead. Zero is treated as the RECALL set (the upstream
	// default) by applyConsensusDefaults.
	BayesianSubMode bayesianMode
	// ConsCutoff is the upstream -C/--cutoff: bayesian calls with a
	// confidence below this become 'N'. Default 10. Because 0 is a valid
	// explicit value, ConsCutoffSet records whether the caller set it.
	ConsCutoff int
	// ConsCutoffSet records whether ConsCutoff was explicitly set; when
	// false, applyConsensusDefaults fills the upstream default of 10.
	ConsCutoffSet bool
	// PHet is the heterozygous-site prior (--P-het, default 1e-3).
	PHet float64
	// PIndel is the indel prior (--P-indel, default 2e-4).
	PIndel float64
	// HetScale scales the heterozygous likelihood (--het-scale, default 1).
	HetScale float64
	// AdjQual enables the localised base-quality adjustment in the
	// NM-halo computation (--adj-qual / --no-adj-qual; default on).
	AdjQual bool
	// AdjQualSet records whether AdjQual was explicitly set by the caller;
	// when false, applyConsensusDefaults turns AdjQual on.
	AdjQualSet bool
	// UseMQual enables the MAPQ-based quality adjustment (--use-MQ /
	// --no-use-MQ; default on).
	UseMQual bool
	// UseMQualSet records whether UseMQual was explicitly set.
	UseMQualSet bool
	// NMAdjust enables the NM-halo localised-MAPQ adjustment (--adj-MQ /
	// --no-adj-MQ; default on).
	NMAdjust bool
	// NMAdjustSet records whether NMAdjust was explicitly set.
	NMAdjustSet bool
	// NMHalo is the window radius for the local NM count (--NM-halo,
	// default 50).
	NMHalo int
	// SCCost is the soft-clip penalty added to localNM (--SC-cost,
	// default 60).
	SCCost int
	// ScaleMQual scales the adjusted MAPQ (--scale-MQ, default 1.0).
	ScaleMQual float64
	// LowMQual / HighMQual clamp the adjusted MAPQ (--low-MQ default 1,
	// --high-MQ default 60).
	LowMQual  int
	HighMQual int
	// DefaultQual substitutes for missing per-base qualities
	// (--default-qual, default 10).
	DefaultQual int
	// HomopolyFix enables the homopolymer quality fix (-p/--homopoly-fix);
	// when non-zero it is also the poly_adj multiplier (--homopoly-score).
	HomopolyFix float64
	// HomopolyRedux is the poly-length quality reduction multiplier
	// (--homopoly-redux, default 0.01).
	HomopolyRedux float64
	// HomopolyReduxSet records whether HomopolyRedux was explicitly set.
	HomopolyReduxSet bool
}

// ConsensusFile is the file-path entry point for `samtools consensus`.
// Opens the input file and delegates to Consensus. errOut is retained for
// API compatibility; both bayesian and simple modes are fully implemented,
// so no fallback warning is emitted.
func ConsensusFile(opts ConsensusOptions, out io.Writer, errOut io.Writer) error {
	if opts.Input == "" {
		return fmt.Errorf("samtools consensus: no input file")
	}
	_ = errOut
	f, err := os.Open(opts.Input)
	if err != nil {
		return fmt.Errorf("samtools consensus: %w", err)
	}
	defer f.Close()
	return Consensus(f, out, opts)
}

// Consensus is the streaming entry point: read records from `in` and
// emit the consensus per the configured format to `out`. Both the simple
// frequency caller and the bayesian Gap5 caller (the upstream default) are
// implemented; opts.Mode selects between them.
func Consensus(in io.Reader, out io.Writer, opts ConsensusOptions) error {
	applyConsensusDefaults(&opts)
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
	resolved, _, err := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
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
	case opts.AllContigs:
		// -aa: walk every contig in the header so empty contigs emit N's.
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

	// -T/--reference: load the reference FASTA so no-coverage / gap positions
	// substitute the reference base for 'N' (bam_consensus.c update_ref).
	var ref *consensusRef
	if opts.Reference != "" {
		ref, err = loadConsensusRef(opts.Reference)
		if err != nil {
			return fmt.Errorf("samtools consensus: %w", err)
		}
		defer ref.close()
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

		if err := emitConsensusContig(bw, chrom, refLen, windows, recs, posFilter, ref, opts); err != nil {
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
	// Bayesian-mode defaults, mirroring upstream main_consensus's
	// consensus_opts initialisers.
	if opts.BayesianSubMode == 0 {
		opts.BayesianSubMode = modeRecall
	}
	if !opts.ConsCutoffSet && opts.ConsCutoff == 0 {
		opts.ConsCutoff = 10
	}
	if opts.PHet == 0 {
		opts.PHet = defaultPHet
	}
	if opts.PIndel == 0 {
		opts.PIndel = defaultPIndel
	}
	if opts.HetScale == 0 {
		opts.HetScale = defaultHetScale
	}
	if !opts.AdjQualSet {
		opts.AdjQual = true
	}
	if !opts.UseMQualSet {
		opts.UseMQual = true
	}
	if !opts.NMAdjustSet {
		opts.NMAdjust = true
	}
	if opts.NMHalo == 0 {
		opts.NMHalo = 50
	}
	if opts.SCCost == 0 {
		opts.SCCost = 60
	}
	if opts.ScaleMQual == 0 {
		opts.ScaleMQual = 1.0
	}
	if opts.LowMQual == 0 {
		opts.LowMQual = 1
	}
	if opts.HighMQual == 0 {
		opts.HighMQual = 60
	}
	if opts.DefaultQual == 0 {
		opts.DefaultQual = 10
	}
	if !opts.HomopolyReduxSet {
		opts.HomopolyRedux = 0.01
	}
}

// bayesOptionsFrom builds the internal bayesOptions from a fully-defaulted
// ConsensusOptions.
func bayesOptionsFrom(opts ConsensusOptions) bayesOptions {
	return bayesOptions{
		mode:        opts.BayesianSubMode,
		useMQual:    opts.UseMQual,
		adjQual:     opts.AdjQual,
		nmAdjust:    opts.NMAdjust,
		nmHalo:      opts.NMHalo,
		scCost:      opts.SCCost,
		scaleMQual:  opts.ScaleMQual,
		lowMQual:    opts.LowMQual,
		highMQual:   opts.HighMQual,
		defaultQual: opts.DefaultQual,
		minQual:     int(opts.MinBaseQ),
		minDepth:    opts.MinDepth,
		consCutoff:  opts.ConsCutoff,
		ambig:       opts.AmbigCodes,
		pHet:        opts.PHet,
		pIndel:      opts.PIndel,
		hetScale:    opts.HetScale,
		homopolyFix: opts.HomopolyFix,
		homopolyRed: opts.HomopolyRedux,
	}
}

// consensusNaNQual is the quality value upstream prints for a simple-mode
// pileup column that has reads but no scoring base (tscore == 0) yet whose
// call escaped the "downgrade to N" branch — typically an intron position
// covered only by ref-skip (CIGAR N) reads. Upstream computes the quality as
// `100.0 * used_score / tscore` == 100.0*0/0 == NaN and casts it to a 32-bit
// C int, which yields INT_MIN on every platform samtools targets. We hold the
// exact sentinel so the pileup quality column matches byte-for-byte.
const consensusNaNQual = -2147483648

// consensusCall is a single per-position call.
type consensusCall struct {
	base byte // 'A','C','G','T','*','N' or IUPAC ambig code when AmbigCodes set
	// qual is the per-position confidence as an integer 0..100
	// (upstream's formula: `100 * used_score / tscore`). The pileup
	// format emits this number verbatim; FASTA/FASTQ encode it as a
	// Phred+33 byte capped at 93.
	qual  int
	depth int // usable depth (post-MinBaseQ)
	// isHet reports whether this position was called heterozygous,
	// determined INDEPENDENTLY of AmbigCodes. In simple mode it is
	// `score2 >= het_fract*score1` on a confidently-called position; in
	// bayesian mode it is a positive het log-odds on a confident call.
	// Used to implement --het-only (HetOnly), which suppresses every
	// non-heterozygous position from the output.
	isHet bool
}

// consensusRef provides per-contig reference bases for the -T/--reference
// substitution. It mirrors upstream bam_consensus.c's update_ref: the full
// contig sequence is fetched once and cached (upstream caches one contig at a
// time; we keep them all since the walk is contig-major). Bases are returned
// 0-based, preserving the FASTA's original case (so soft-masked lowercase
// bases are emitted verbatim, exactly as upstream's c->ref[pos] does).
//
// For a plain FASTA the contig is read raw via the .fai line geometry so byte
// case is preserved. For a BGZF-compressed FASTA we fall back to the shared
// fasta.RandomAccess.Fetch path, which canonicalises to uppercase — a rare
// case for -T and one upstream itself does not soft-mask differently.
type consensusRef struct {
	ra    *fasta.RandomAccess
	idx   *fasta.Index
	path  string
	plain bool // true when path is a non-BGZF FASTA we can read raw
	cache map[string][]byte
}

// loadConsensusRef opens path (an indexed FASTA; the .fai is built on demand
// when absent) for reference-base substitution.
func loadConsensusRef(path string) (*consensusRef, error) {
	ra, err := fasta.OpenRandomAccess(path)
	if err != nil {
		return nil, fmt.Errorf("could not load reference %q: %w", path, err)
	}
	return &consensusRef{
		ra:    ra,
		idx:   ra.Index(),
		path:  path,
		plain: !fastaIsGzip(path),
		cache: make(map[string][]byte),
	}, nil
}

// fastaIsGzip reports whether path begins with the gzip/BGZF magic (0x1f 0x8b),
// in which case raw byte-offset reads are not possible and the consensus
// reference falls back to the (uppercasing) random-access fetch.
func fastaIsGzip(path string) bool {
	f, err := os.Open(path)
	if err != nil {
		return false
	}
	defer f.Close()
	var magic [2]byte
	if _, err := io.ReadFull(f, magic[:]); err != nil {
		return false
	}
	return magic[0] == 0x1f && magic[1] == 0x8b
}

// close releases the underlying FASTA handle.
func (r *consensusRef) close() {
	if r != nil && r.ra != nil {
		r.ra.Close()
	}
}

// contig returns the full reference sequence for name, preserving the FASTA's
// original byte case for plain FASTAs. It returns nil when the contig is not
// in the reference (upstream's update_ref returns <0 and the caller falls back
// to 'N').
func (r *consensusRef) contig(name string) []byte {
	if r == nil {
		return nil
	}
	if seq, ok := r.cache[name]; ok {
		return seq
	}
	entry, ok := r.idx.Get(name)
	if !ok {
		r.cache[name] = nil
		return nil
	}
	var seq []byte
	if r.plain {
		var err error
		if seq, err = readRawContig(r.path, entry); err != nil {
			seq = nil
		}
	}
	if seq == nil {
		// BGZF input or a raw-read failure: fall back to the uppercasing
		// random-access fetch.
		if b, err := r.ra.Fetch(name, 0, entry.Length); err == nil {
			seq = b
		}
	}
	r.cache[name] = seq
	return seq
}

// base returns the 0-based reference base at pos on name, or 'N' when the
// reference does not cover it. The substitution table is otherwise verbatim.
func (r *consensusRef) base(name string, pos0 int) byte {
	seq := r.contig(name)
	if pos0 < 0 || pos0 >= len(seq) {
		return 'N'
	}
	return seq[pos0]
}

// readRawContig reads the whole contig from a plain FASTA preserving byte
// case. fasta.RandomAccess.Fetch uppercases its output (it is built for
// case-insensitive left-alignment), but upstream consensus copies the
// reference bytes verbatim, so we re-read the raw bytes here using the .fai
// line geometry to keep soft-masked lowercase bases byte-faithful.
func readRawContig(path string, entry fasta.IndexEntry) ([]byte, error) {
	if entry.Length == 0 {
		return []byte{}, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	// Total raw byte span of the contig including line terminators, derived
	// from the .fai line geometry (LineBases bases per LineWidth bytes).
	fullLines := entry.Length / entry.LineBases
	rem := entry.Length % entry.LineBases
	rawLen := fullLines * entry.LineWidth
	if rem > 0 {
		rawLen += (entry.LineWidth - entry.LineBases) + rem
	}
	buf := make([]byte, rawLen)
	if _, err := f.ReadAt(buf, entry.Offset); err != nil && err != io.EOF {
		return nil, err
	}
	out := make([]byte, 0, entry.Length)
	for _, b := range buf {
		if b == '\n' || b == '\r' {
			continue
		}
		out = append(out, b)
	}
	if int64(len(out)) != entry.Length {
		return nil, fmt.Errorf("reference %s: read %d bases, expected %d", entry.Name, len(out), entry.Length)
	}
	return out, nil
}

// emitConsensusContig walks every window on chrom and emits per-format
// records (FASTA/FASTQ/Pileup). ref is non-nil when -T/--reference is in
// effect, in which case no-coverage / gap positions substitute the reference
// base for 'N'.
func emitConsensusContig(bw *bufio.Writer, chrom string, refLen int, windows [][2]int,
	recs []*sam.Record, posFilter *positionFilter, ref *consensusRef, opts ConsensusOptions) error {

	// For FASTA/FASTQ we accumulate one buffer per contig and emit at
	// the end. For pileup we stream line-by-line. firstCovIdx/lastCovIdx
	// bracket the covered span within seqBuf: every position is appended
	// (uncovered ones as 'N') so internal gaps are filled, then leading
	// and trailing N runs are trimmed unless -a is set.
	var seqBuf, qualBuf []byte
	firstCovIdx, lastCovIdx := -1, -1

	// pileupLastPos tracks upstream's c->last_pos for the -a/--all-positions
	// pileup placeholder mechanism (basic_pileup, bam_consensus.c:2202). It
	// is the 1-based reference position of the last row actually emitted; 0
	// means "nothing emitted yet on this contig" (upstream resets last_pos to
	// the region start — 0 for a whole-contig walk — at the first column of a
	// new tid). Zero-coverage and suppressed deletion-only positions are NOT
	// emitted inline; they are filled lazily as placeholder rows when the next
	// genuine pileup column is reached, and finally by a tail fill at the end
	// of the contig. Because a suppressed deletion column does not advance
	// last_pos, the same gap span is re-emitted at each following column —
	// reproducing upstream's quirky duplicate placeholder rows at deletion
	// sites.
	pileupAll := opts.Format == ConsensusPileup && opts.AllPositions
	// pileupLastPos tracks the 1-based reference position of the last pileup
	// row actually emitted within the current window; it is reset to the
	// window's start at the top of each window so placeholder fill is scoped
	// to the requested span and never bridges the void between two disjoint
	// regions. Upstream clamps the fill to the iterator's [beg, end) when a
	// region is given (bam_consensus.c:2832-2842); processing one window at a
	// time reproduces that per region.
	pileupLastPos := 0

	// For bayesian mode, build the per-read NM-halo state once per
	// contig (indexed by the record's position in recs) and the
	// parameter-set matrices once.
	var bayes bayesOptions
	var bayesProbs bayesProbSet
	var bayesReads []*bayesRead
	if opts.Mode == ConsensusModeBayesian {
		bayes = bayesOptionsFrom(opts)
		bayesProbs = buildBayesProbSet(bayes)
		bayesReads = make([]*bayesRead, len(recs))
		for i, rec := range recs {
			bayesReads[i] = nmInit(rec, bayes)
		}
	}

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
		// Reset the placeholder cursor to this window's start (1-based
		// last_pos == beg0 means "nothing emitted yet inside the window", so
		// the first fill covers beg0+1..). Whole-contig walks use beg0 == 0.
		pileupLastPos = beg0
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
			var call consensusCall
			var totalDepth int
			if opts.Mode == ConsensusModeBayesian {
				call, totalDepth = callConsensusBayesian(events[col], recs, bayesReads, bayes, bayesProbs)
			} else {
				call, totalDepth = callConsensus(events[col], opts)
			}

			switch opts.Format {
			case ConsensusPileup:
				column := hasPileupColumn(events[col])
				if !pileupAll {
					// Without -a, zero-coverage positions never produce a
					// row at all.
					if !column {
						continue
					}
				} else if !column {
					// With -a, a genuine zero-coverage position is NOT a
					// pileup callback in upstream; it is emitted only by the
					// placeholder gap mechanism when the next real column is
					// reached (or by the contig tail fill). Defer it.
					continue
				}
				// We are at a genuine pileup column. Upstream's basic_pileup
				// is invoked once per nth (nth==0 base column, then one call
				// per insertion column) and re-runs its gap-fill / suppression
				// / last_pos block on EVERY call (bam_consensus.c:2185-2298).
				// We model that here by looping nth over the base column and
				// each insertion column, performing the same per-nth steps so
				// the output matches upstream line-for-line — including its
				// quirk of re-emitting a leading gap at each nth while last_pos
				// stays unadvanced through a suppressed nth==0 deletion.
				var insCols []bayesInsertionColumn
				if !opts.NoShowIns {
					insCols = consensusInsertionColumns(events[col], recs, bayesReads, bayes, bayesProbs, pos1, opts)
				}
				// --het-only is an intentional divergence (upstream parses but
				// never acts on it — a dead option; see docs/UPSTREAM_BUGS.md).
				// A non-het position is dropped entirely: under -a we still
				// advance pileupLastPos past it so the het-suppressed column is
				// not resurrected as a zero-depth placeholder row by a later
				// column's gap fill, honouring --het-only's "drop entirely"
				// contract. This short-circuits the whole position (nth==0 and
				// all insertion columns), since the position is suppressed.
				if opts.HetOnly && !call.isHet {
					if pileupAll {
						pileupLastPos = pos1
					}
					continue
				}
				for nth := 0; nth <= len(insCols); nth++ {
					// Per-nth gap fill: upstream re-runs empty_pileup2 at the
					// top of every basic_pileup call (bam_consensus.c:2227), so
					// a gap preceding a position whose nth==0 deletion row is
					// suppressed (leaving last_pos unadvanced) is re-emitted at
					// the next nth — reproduced here by gap-filling inside the
					// nth loop rather than once before it.
					if pileupAll && pos1 > pileupLastPos+1 {
						if err := writeEmptyPileupRows(bw, chrom, pileupLastPos+1, pos1-1, posFilter, ref); err != nil {
							return err
						}
					}
					if nth == 0 {
						// Honour --show-del: when the call is '*' and ShowDel
						// is false, suppress the nth==0 row (bam_consensus.c:
						// 2244). This suppresses ONLY the reference row; the
						// nth>0 insertion columns follow on their own merits.
						// A suppressed deletion column does NOT advance
						// pileupLastPos (upstream returns early without updating
						// c->last_pos), so the gap before it is re-filled at the
						// next nth — upstream's duplicate placeholder rows.
						if call.base == '*' && !opts.ShowDel {
							continue
						}
						if err := writeConsensusPileupRow(bw, chrom, pos1, totalDepth, call, events[col], opts); err != nil {
							return err
						}
						// Upstream sets c->last_pos = pos when any row for the
						// position is emitted, including an insertion row; we set
						// it here once the first non-suppressed nth emits.
						pileupLastPos = pos1
						continue
					}
					// nth>0 insertion column. Upstream emits one pileup row per
					// inserted column when --show-ins is on (the default), for
					// both simple and bayesian modes.
					ic := insCols[nth-1]
					if ic.call.base == '*' && !opts.ShowDel {
						continue
					}
					if err := writeConsensusInsertionPileupRow(bw, chrom, pos1, nth, ic, opts); err != nil {
						return err
					}
					// Emitting an insertion row also advances last_pos
					// (bam_consensus.c:2295). When the nth==0 deletion was
					// suppressed this is the only place the position's last_pos
					// gets set, preventing a spurious re-emission downstream.
					pileupLastPos = pos1
				}
			default:
				// FASTA / FASTQ accumulate. Every position is appended
				// (uncovered ones as 'N', or the reference base under
				// -T/--reference) so internal gaps fill; the covered span is
				// bracketed by firstCovIdx/lastCovIdx and leading/trailing N is
				// trimmed unless -a. Upstream's basic_fasta gap fill
				// (bam_consensus.c:2423-2436) copies c->ref[pos] for the gap
				// bytes and sets their qual to ref_qual+'!' (else 'N'/'!').
				if call.base == 0 {
					if ref != nil {
						seqBuf = append(seqBuf, ref.base(chrom, pos0))
						qualBuf = append(qualBuf, byte(opts.RefQual)+'!')
					} else {
						seqBuf = append(seqBuf, 'N')
						qualBuf = append(qualBuf, '!')
					}
					continue
				}
				// --het-only: render every non-heterozygous position as
				// 'N' to preserve coordinates (homozygous and no-call
				// positions are masked, not deleted). We treat the 'N'
				// exactly like an uncovered position — it fills internal
				// gaps but does not extend the covered span, so leading
				// and trailing non-het runs trim away (unless -a). The
				// intended behaviour the flag implies; upstream samtools
				// parses --het-only but never acts on it (a dead option —
				// see docs/UPSTREAM_BUGS.md).
				if opts.HetOnly && !call.isHet {
					seqBuf = append(seqBuf, 'N')
					qualBuf = append(qualBuf, '!')
					continue
				}
				// Emit the reference (nth==0) base unless it is a
				// deletion placeholder suppressed by --show-del off. The
				// suppression is scoped to the nth==0 base only: upstream
				// invokes basic_fasta independently per nth, so a
				// deletion-called reference column still has its nth>0
				// insertion columns emitted below.
				if !(call.base == '*' && !opts.ShowDel) {
					if firstCovIdx < 0 {
						firstCovIdx = len(seqBuf)
					}
					seqBuf = append(seqBuf, call.base)
					qualBuf = append(qualBuf, phredByte(call.qual))
					lastCovIdx = len(seqBuf)
				}
				if !opts.NoShowIns {
					// Insertion columns (nth>0). Upstream's basic_fasta
					// emits the inserted column's call whenever cb != '*'
					// (bam_consensus.c:2439) — including 'N' calls — and,
					// under --mark-ins, prepends a single '_' to BOTH the
					// seq and qual stream once per inserted column
					// (bam_consensus.c:2409-2412). A '*' inserted call is
					// never emitted. This holds for simple and bayesian
					// modes alike.
					insCols := consensusInsertionColumns(events[col], recs, bayesReads, bayes, bayesProbs, pos1, opts)
					for _, ic := range insCols {
						cb := ic.call.base
						if cb == 0 {
							cb = 'N'
						}
						if cb == '*' {
							continue
						}
						if firstCovIdx < 0 {
							firstCovIdx = len(seqBuf)
						}
						if opts.MarkIns {
							seqBuf = append(seqBuf, '_')
							qualBuf = append(qualBuf, '_')
						}
						seqBuf = append(seqBuf, cb)
						qualBuf = append(qualBuf, phredByte(ic.call.qual))
						lastCovIdx = len(seqBuf)
					}
				}
			}
		}

		// Pileup -a tail fill: emit placeholder rows for the remainder of
		// this window after the last emitted row, through to the window end
		// (basic_pileup's end-of-loop empty_pileup2, bam_consensus.c:2832-
		// 2842). This covers a deletion-only or zero-coverage run trailing
		// the final genuine row, and the window's leading gap when it has no
		// coverage at all.
		if pileupAll && pileupLastPos < end0 {
			if err := writeEmptyPileupRows(bw, chrom, pileupLastPos+1, end0, posFilter, ref); err != nil {
				return err
			}
		}
	}

	if opts.Format != ConsensusPileup {
		// Trim leading/trailing uncovered N runs unless -a is set:
		// upstream emits only the covered span of each contig by
		// default, internal gaps included, and extends to the full
		// contig only with -a.
		if !opts.AllPositions {
			if firstCovIdx < 0 {
				// No coverage at all and not -a: emit nothing.
				return nil
			}
			seqBuf = seqBuf[firstCovIdx:lastCovIdx]
			qualBuf = qualBuf[firstCovIdx:lastCovIdx]
		}
		if len(seqBuf) == 0 && !opts.AllContigs {
			return nil
		}
		// -aa on an untouched contig: emit a record full of N's, or of the
		// reference bases under -T/--reference (upstream's empty-chr branch
		// fills the whole contig via append_cons -> ref_or_Ns,
		// bam_consensus.c:2636-2658).
		if len(seqBuf) == 0 && opts.AllContigs {
			seqBuf = make([]byte, refLen)
			qualBuf = make([]byte, refLen)
			for i := range seqBuf {
				if ref != nil {
					seqBuf[i] = ref.base(chrom, i)
					qualBuf[i] = byte(opts.RefQual) + '!'
				} else {
					seqBuf[i] = 'N'
					qualBuf[i] = '!'
				}
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
		if e.kind == pileupEventRefSkip {
			// Reference skips contribute no base/gap SCORE (upstream's
			// base4 is 0, so every seqi2X[0] weight is 0), but upstream
			// still counts them in tot_depth (bam_consensus.c:1955 runs
			// unconditionally after the b<16 branch). Upstream's pileup
			// engine assigns a ref-skip quality of 0 (consensus_pileup.c:
			// 216), so a non-zero --min-BQ excludes it via the same min-qual
			// gate that drops a base (q < min_qual → continue, before the
			// tot_depth++). We therefore gate on quality 0 — independent of
			// the event's carried (read-base) quality — and only then count
			// it. An intron-only column thus stays a confident 'N' (no
			// score) rather than an empty no-call.
			if opts.MinBaseQ == 0 {
				totDepth++
			}
			continue
		}
		// Upstream applies the min-qual gate to EVERY event (base and
		// gap alike) before bucketing, using the event's carried quality
		// (bam_consensus.c:1925-1931). A pad/deletion event below the
		// floor is therefore dropped just like a base.
		if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
			continue
		}
		switch e.kind {
		case pileupEventBase:
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
			// Gap bucket. Upstream weights the fixed gap score by the
			// event quality under --use-qual: score[16] += 8 * (use_qual
			// ? q : 1) (bam_consensus.c:1953).
			var q uint64 = 1
			if opts.UseQual {
				q = uint64(e.qual)
			}
			freq[16]++
			score[16] += 8 * q
			totDepth++
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
	//
	// hetSite records whether the two top alleles satisfy the
	// heterozygosity test INDEPENDENTLY of --ambig. It is the gate for
	// --het-only (HetOnly). We require s1>0 so that an empty/uncalled
	// column is not spuriously flagged het (when s1==0, s2==0 too and
	// 0 >= 0 would otherwise hold).
	hetSite := s1 > 0 && float64(s2) >= opts.MinHetFraction*float64(s1)
	if opts.AmbigCodes && hetSite {
		used |= call2
		usedScore += s2
	}

	// Single fraction gate, mirroring upstream bam_consensus.c:1988-1994
	// EXACTLY: `tot_depth < min_depth || used_score < call_fract*tscore`.
	// Crucially there is NO `tscore > 0` guard: when tscore is 0 (a column
	// with reads but no scoring base — e.g. an intron covered only by
	// ref-skip reads), the comparison `0 < 0` is false, so the downgrade is
	// NOT triggered and used_base retains call1, which stays at its initial
	// 15 ('N') because no score ever beat 0. A `tscore > 0` guard would
	// wrongly force the downgrade here and lose upstream's INT_MIN quality.
	notCall := totDepth < opts.MinDepth ||
		float64(usedScore) < opts.MinCallFraction*float64(tscore)
	// isHet for --het-only is computed independently of --ambig: a
	// position counts as heterozygous when the two top alleles pass the
	// het-fract test (hetSite) AND the het-inclusive call (call1|call2)
	// is itself confident — i.e. it clears the same depth/call-fraction
	// gates upstream uses to accept a call. We evaluate the call gate on
	// the het-inclusive score (s1+s2) so the determination does not
	// depend on whether --ambig happened to widen usedScore above. The
	// `tscore > 0` guard here is deliberate: --het-only is our own
	// extension (no upstream oracle), and a het call on a zero-score
	// column is meaningless.
	hetInclusiveScore := s1
	if hetSite {
		hetInclusiveScore += s2
	}
	hetCallOK := tscore > 0 && float64(hetInclusiveScore) >= opts.MinCallFraction*float64(tscore)
	isHet := hetSite && totDepth >= opts.MinDepth && hetCallOK
	if notCall {
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

	// Confidence, mirroring upstream's `*qual = used_base ? 100.0 *
	// used_score / tscore : 0` (bam_consensus.c:2003). When used_base is
	// non-zero but tscore is 0 (the intron / ref-skip-only column whose
	// call1 stayed 15), upstream computes 100.0*0/0 == NaN and casts it to
	// a C int, which on every platform samtools targets yields INT_MIN
	// (-2147483648). We reproduce that sentinel exactly so the pileup row's
	// quality column is byte-for-byte identical.
	qual := 0
	if used != 0 {
		if tscore == 0 {
			qual = consensusNaNQual
		} else {
			qual = int(100 * float64(usedScore) / float64(tscore))
		}
	}
	return consensusCall{base: base, qual: qual, depth: totDepth, isHet: isHet}, totDepth
}

// callConsensusBayesian runs the Gap5-derived bayesian caller on one
// column. It builds the per-base bayesPileupBase view from the pileup
// events (recovering each read's NM-halo state by readIdx), invokes
// calculateConsensusGap5m, and maps the result to a consensusCall.
//
// The returned base is 0 (and depth 0) when the column has no usable
// coverage so the FASTA/FASTQ caller treats it as a skip, matching the
// simple path.
func callConsensusBayesian(evs []pileupEvent, recs []*sam.Record,
	bayesReads []*bayesRead, o bayesOptions, ps bayesProbSet) (consensusCall, int) {

	bases := make([]bayesPileupBase, 0, len(evs))
	td := 0
	for _, e := range evs {
		if e.dropped {
			continue
		}
		var bp bayesPileupBase
		bp.mapQ = e.mapq
		if e.readIdx >= 0 && e.readIdx < len(bayesReads) {
			bp.read = bayesReads[e.readIdx]
		}
		if e.readIdx >= 0 && e.readIdx < len(recs) {
			bp.readPos0 = int(recs[e.readIdx].Pos) - 1
		}
		switch e.kind {
		case pileupEventBase:
			bp.base4 = byte(baseToSeqi(upper(e.base)))
			bp.qual = e.qual
			bp.seqOff = e.readBP - 1
			// An exon base abutting a ref-skip (CIGAR N) run carries
			// upstream's p->ref_skip flag, so the Gap5/bayesian caller
			// excludes it from the consensus depth (bam_consensus.c:1333)
			// even though it is a real, displayed base. Mirror that by
			// flagging the gap5 view as a ref-skip; the displayed pileup
			// column still renders e.base via writeConsensusPileupRow.
			if e.refSkipBoundary {
				bp.refSkip = true
			}
		case pileupEventDel:
			bp.base4 = 16
			bp.qual = e.qual
			bp.seqOff = e.readBP - 1
			// A deletion ('*') position directly abutting a ref-skip run
			// also carries upstream's p->ref_skip flag (consensus_pileup.c:
			// 240 tests `p->base != '.'`, which a deletion satisfies), so
			// the Gap5/bayesian caller excludes it from cons.depth just
			// like an exon base boundary. Mirror that here.
			if e.refSkipBoundary {
				bp.refSkip = true
			}
		case pileupEventRefSkip:
			bp.refSkip = true
			bp.seqOff = e.readBP - 1
		}
		bases = append(bases, bp)
		td++
	}

	cons := calculateConsensusGap5m(bases, o.useMQual, td, o, ps)
	cb, cq := bayesCallToBase(cons, o)

	// Depth 0: signal "skip" to the FASTA/FASTQ accumulator.
	if cons.depth == 0 && cons.call == 4 {
		// A genuine all-N / empty column. Upstream emits 'N' qual 0
		// in pileup; the FASTA path skips it unless -a is set, which
		// is handled by the base==0 sentinel below.
		if td == 0 {
			return consensusCall{}, 0
		}
	}
	// isHet for --het-only, determined independently of --ambig. A
	// positive het log-odds means the bayesian caller favours a
	// two-allele genotype; we additionally require that genotype to
	// survive the SAME confidence gates bayesCallToBase applies when it
	// emits a het IUPAC base (with --ambig assumed on): the depth floor,
	// and the cutoff-downgrade that would otherwise turn the call into
	// 'N'. Mirroring bayesCallToBase exactly (rather than a looser
	// hetLogOdd>=consCutoff test) keeps the het-only decision consistent
	// with the base the caller would actually render — including the
	// corner case of a heterozygote involving a gap allele, which
	// bayesCallToBase exempts from the cutoff downgrade.
	isHet := cons.hetLogOdd > 0 && cons.depth >= o.minDepth
	if isHet {
		// Replicate bayesCallToBase's downgrade-to-N test for the het
		// branch (cb is the het IUPAC code, cq = hetLogOdd). The base is
		// never '*' on the het branch, so only the cutoff + gap-allele
		// exemption matter.
		if cons.hetLogOdd < o.consCutoff &&
			cons.hetCall%5 != 4 && cons.hetCall/5 != 4 {
			isHet = false
		}
	}
	return consensusCall{base: cb, qual: cq, depth: cons.depth, isHet: isHet}, td
}

// bayesInsertionColumn is one nth>0 insertion column produced for a
// reference position: the bayesian call plus the per-read base/qual bytes
// for the pileup seq/qual columns.
type bayesInsertionColumn struct {
	call  consensusCall
	depth int    // raw read count in this column (all spanning reads)
	seq   []byte // per-read base bytes (upper/lower-cased), '*' for pads
	qual  []byte // per-read qual bytes (Phred+33)
}

// callConsensusBayesianInsertions builds the nth>0 insertion columns for a
// reference column and runs the bayesian caller on each. A read in the base
// column participates only when it survives into the insertion column (see
// spansInsertionColumn): a read with an nth inserted base contributes it,
// any other surviving read contributes a '*' pad — matching upstream's
// pileup engine, which pads non-inserting reads but drops reads that end at
// the reference position. pos1 is the 1-based reference position.
func callConsensusBayesianInsertions(evs []pileupEvent, recs []*sam.Record,
	bayesReads []*bayesRead, o bayesOptions, ps bayesProbSet, pos1 int) []bayesInsertionColumn {

	// Stable order so the per-read seq/qual columns match upstream.
	sortEvents(evs)

	maxIns := 0
	for _, e := range evs {
		if e.dropped || e.kind != pileupEventBase {
			continue
		}
		if len(e.insAfter) > maxIns {
			maxIns = len(e.insAfter)
		}
	}
	if maxIns == 0 {
		return nil
	}

	cols := make([]bayesInsertionColumn, maxIns)
	for nth := 1; nth <= maxIns; nth++ {
		bases := make([]bayesPileupBase, 0, len(evs))
		seq := make([]byte, 0, len(evs))
		qual := make([]byte, 0, len(evs))
		td := 0
		for _, e := range evs {
			if e.dropped || e.kind == pileupEventRefSkip {
				continue
			}
			if !spansInsertionColumn(e, recs, pos1) {
				continue
			}
			var bp bayesPileupBase
			bp.mapQ = e.mapq
			if e.readIdx >= 0 && e.readIdx < len(bayesReads) {
				bp.read = bayesReads[e.readIdx]
			}
			var b byte = '*'
			var q byte
			if e.kind == pileupEventBase && nth <= len(e.insAfter) {
				ib := upper(e.insAfter[nth-1])
				bp.base4 = byte(baseToSeqi(ib))
				// The nth inserted base sits at query offset
				// (readBP-1)+nth in the read's SEQ.
				bp.seqOff = e.readBP - 1 + nth
				var rec *sam.Record
				if e.readIdx >= 0 && e.readIdx < len(recs) {
					rec = recs[e.readIdx]
				}
				if rec != nil && bp.seqOff >= 0 && bp.seqOff < len(rec.Qual) {
					q = rec.Qual[bp.seqOff]
				}
				bp.qual = q
				b = ib
				if e.isReverse {
					b = lower(b)
				}
			} else {
				// Pad: '*' base. Upstream's pileup engine carries the
				// read's current base quality and position into the
				// insertion column for non-inserting reads.
				bp.base4 = 16
				bp.seqOff = e.readBP - 1
				q = e.qual
				bp.qual = q
				if e.isReverse {
					b = '#'
				}
			}
			bases = append(bases, bp)
			seq = append(seq, b)
			qual = append(qual, q+33)
			td++
		}
		cons := calculateConsensusGap5m(bases, o.useMQual, td, o, ps)
		cb, cq := bayesCallToBase(cons, o)
		cols[nth-1] = bayesInsertionColumn{
			call:  consensusCall{base: cb, qual: cq, depth: cons.depth},
			depth: td,
			seq:   seq,
			qual:  qual,
		}
	}
	return cols
}

// callConsensusSimpleInsertions builds the nth>0 insertion columns for a
// reference column and runs the simple frequency caller on each, mirroring
// upstream's pileup engine (consensus_pileup.c): for each inserted column
// every spanning read participates — a read with an nth inserted base
// contributes it, otherwise it contributes a '*' pad. The simple caller
// (callConsensus) is then run on the synthesised column exactly as
// consensus_base() dispatches to calculate_consensus_simple() for nth>0
// columns. This makes simple-mode insertion handling byte-faithful with
// upstream in both FASTA/FASTQ and pileup output.
//
// The returned columns carry the per-read seq/qual bytes (upper/lower-cased
// by strand, '*'/'#' for pads) needed by the pileup seq/qual fields. pos1 is
// the 1-based reference position whose insertion follows.
func callConsensusSimpleInsertions(evs []pileupEvent, recs []*sam.Record,
	pos1 int, opts ConsensusOptions) []bayesInsertionColumn {

	// Stable order so the per-read seq/qual columns match upstream.
	sortEvents(evs)

	maxIns := 0
	for _, e := range evs {
		if e.dropped || e.kind != pileupEventBase {
			continue
		}
		if len(e.insAfter) > maxIns {
			maxIns = len(e.insAfter)
		}
	}
	if maxIns == 0 {
		return nil
	}

	cols := make([]bayesInsertionColumn, maxIns)
	for nth := 1; nth <= maxIns; nth++ {
		colEvs := make([]pileupEvent, 0, len(evs))
		seq := make([]byte, 0, len(evs))
		qual := make([]byte, 0, len(evs))
		td := 0
		for _, e := range evs {
			if e.dropped || e.kind == pileupEventRefSkip {
				continue
			}
			if !spansInsertionColumn(e, recs, pos1) {
				continue
			}
			var ce pileupEvent
			ce.readIdx = e.readIdx
			ce.mapq = e.mapq
			ce.isReverse = e.isReverse
			var b byte = '*'
			var q byte
			if e.kind == pileupEventBase && nth <= len(e.insAfter) {
				ib := upper(e.insAfter[nth-1])
				ce.kind = pileupEventBase
				ce.base = ib
				// The nth inserted base sits at query offset
				// (readBP-1)+nth in the read's SEQ.
				seqOff := e.readBP - 1 + nth
				var rec *sam.Record
				if e.readIdx >= 0 && e.readIdx < len(recs) {
					rec = recs[e.readIdx]
				}
				if rec != nil && seqOff >= 0 && seqOff < len(rec.Qual) {
					q = rec.Qual[seqOff]
				}
				ce.qual = q
				b = ib
				if e.isReverse {
					b = lower(b)
				}
			} else {
				// Pad: '*' deletion placeholder. Upstream's pileup engine
				// carries the read's current base quality into the column.
				ce.kind = pileupEventDel
				q = e.qual
				ce.qual = q
				if e.isReverse {
					b = '#'
				}
			}
			colEvs = append(colEvs, ce)
			seq = append(seq, b)
			qual = append(qual, q+33)
			td++
		}
		call, _ := callConsensus(colEvs, opts)
		cols[nth-1] = bayesInsertionColumn{
			call:  call,
			depth: td,
			seq:   seq,
			qual:  qual,
		}
	}
	return cols
}

// consensusInsertionColumns dispatches to the simple or bayesian insertion
// column builder per opts.Mode, returning one bayesInsertionColumn per
// inserted nth-column. This is the single entry point both the FASTA/FASTQ
// and pileup emitters use so the two modes share identical insertion
// semantics with upstream. pos1 is the 1-based reference position of the
// base column whose insertion follows.
func consensusInsertionColumns(evs []pileupEvent, recs []*sam.Record,
	bayesReads []*bayesRead, bayes bayesOptions, bayesProbs bayesProbSet,
	pos1 int, opts ConsensusOptions) []bayesInsertionColumn {

	if opts.Mode == ConsensusModeBayesian {
		return callConsensusBayesianInsertions(evs, recs, bayesReads, bayes, bayesProbs, pos1)
	}
	return callConsensusSimpleInsertions(evs, recs, pos1, opts)
}

// spansInsertionColumn reports whether the read behind event e participates
// in the insertion column that follows reference position pos1. Upstream's
// pileup engine (consensus_pileup.c::get_next_base) only keeps a read in an
// nth>0 column while the read still has unconsumed bases: a read whose
// alignment terminates exactly at pos1 reaches EOF and is removed before the
// insertion column, so it neither inserts nor pads there. A read therefore
// participates iff it carries an insertion at pos1 (insAfter set) or its
// alignment extends to a later reference position.
func spansInsertionColumn(e pileupEvent, recs []*sam.Record, pos1 int) bool {
	if e.insAfter != "" {
		return true
	}
	if e.readIdx < 0 || e.readIdx >= len(recs) {
		return false
	}
	return int(recs[e.readIdx].EndPosition()) > pos1
}

// hasPileupColumn reports whether the position's event slice represents a
// genuine pileup column — i.e. at least one read overlaps it (via a base,
// deletion, or reference-skip CIGAR op). This mirrors the set of positions at
// which upstream's pileup engine invokes its per-column callback
// (basic_pileup); zero-coverage positions are never callbacks and are emitted
// only via the placeholder-row gap mechanism under -a/--all-positions.
func hasPileupColumn(evs []pileupEvent) bool {
	for _, e := range evs {
		if !e.dropped {
			return true
		}
	}
	return false
}

// writeEmptyPileupRows emits one placeholder pileup row per reference position
// in the inclusive 1-based range [start1, end1], reproducing upstream's
// empty_pileup2 (bam_consensus.c:2101). Each row is
// "<chrom>\t<pos>\t0\t0\t<C>\t0\t*\t*\n": nth 0, depth 0, the call <C>, quality
// 0, and '*' for both the seq and qual columns. The call is 'N' unless a
// -T/--reference FASTA is loaded (ref != nil), in which case it is the
// reference base at that position (upstream's `rseq ? rseq[i] : 'N'`). A
// reference position outside the loaded contig falls back to 'N'.
//
// When posFilter is non-nil (the -l/--positions convenience), positions it
// excludes are skipped, so the placeholder fill never emits rows for
// coordinates the caller restricted out — matching the inline row path, which
// also honours the filter.
func writeEmptyPileupRows(bw *bufio.Writer, chrom string, start1, end1 int, posFilter *positionFilter, ref *consensusRef) error {
	for pos1 := start1; pos1 <= end1; pos1++ {
		if posFilter != nil && !posFilter.contains(chrom, pos1) {
			continue
		}
		if _, err := bw.WriteString(chrom); err != nil {
			return err
		}
		bw.WriteByte('\t')
		bw.WriteString(strconv.Itoa(pos1))
		bw.WriteString("\t0\t0\t")
		cb := byte('N')
		if ref != nil {
			cb = ref.base(chrom, pos1-1)
		}
		bw.WriteByte(cb)
		bw.WriteString("\t0\t*\t*\n")
	}
	return nil
}

// writeConsensusPileupRow emits one row of the samtools consensus
// pileup schema: chrom\tpos\tnth\tdepth\tcall\tcq\tseq\tqual\n.
//
// This writes the nth==0 (reference) row. The nth>0 insertion-column
// rows that upstream emits under `--show-ins yes` (the default) are
// written separately by writeConsensusInsertionPileupRow.
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
			// The displayed seq/qual columns and the row depth are the
			// RAW pileup column — upstream's basic_pileup loop emits every
			// overlapping read's base/qual with no min-BQ filter
			// (bam_consensus.c:2274-2285); --min-BQ affects only the
			// consensus call, not the displayed column.
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
			// Upstream's consensus pileup engine keeps a ref-skip (CIGAR N)
			// read in the column, counts it in depth, and renders its base as
			// '.' with quality 0 (consensus_pileup.c:214-220 sets p->base='.',
			// p->qual=0; basic_pileup then emits it in the seq/qual columns and
			// counts it via depth++). It is NOT lowercased on the reverse strand
			// — '.' has no case. So an intron position covered only by spliced
			// reads prints those reads' '.' rather than collapsing to depth 0.
			b = '.'
			q = 0
		}
		if q > 93 {
			q = 93 // upstream MIN(qual,93) before +'!'
		}
		seq = append(seq, b)
		qual = append(qual, q+33)
	}
	// Row depth is the raw column read count (the number of displayed
	// reads), matching upstream's `depth` argument — distinct from the
	// min-BQ-filtered depth used by the consensus caller.
	depth = len(seq)

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

// writeConsensusInsertionPileupRow emits one nth>0 pileup row for an
// inserted column: chrom\tpos\tnth\tdepth\tcall\tcq\tseq\tqual\n.
func writeConsensusInsertionPileupRow(bw *bufio.Writer, chrom string, pos1, nth int,
	ic bayesInsertionColumn, opts ConsensusOptions) error {

	cb := ic.call.base
	if cb == 0 {
		cb = 'N'
	}
	if _, err := bw.WriteString(chrom); err != nil {
		return err
	}
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(pos1))
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(nth))
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(ic.depth))
	bw.WriteByte('\t')
	bw.WriteByte(cb)
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(ic.call.qual))
	bw.WriteByte('\t')
	if len(ic.seq) == 0 {
		bw.WriteByte('*')
	} else {
		bw.Write(ic.seq)
	}
	bw.WriteByte('\t')
	if len(ic.qual) == 0 {
		bw.WriteByte('*')
	} else {
		bw.Write(ic.qual)
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
// "bayesian_p", "bayesian_116". The bayesian sub-mode is reported by
// ParseConsensusBayesianMode.
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

// SetBayesianMode selects the bayesian parameter set from a CLI -m/--mode
// string. "bayesian" and "bayesian_r" select RECALL (the upstream
// default); "bayesian_p" PRECISE, "bayesian_m" MIXED, "bayesian_116" the
// samtools-1.16 parameter set. Any other value (including the simple-mode
// strings) leaves the RECALL set, which is unused unless Mode is
// ConsensusModeBayesian.
func (o *ConsensusOptions) SetBayesianMode(s string) {
	switch strings.ToLower(s) {
	case "bayesian_p":
		o.BayesianSubMode = modePrecise
	case "bayesian_m":
		o.BayesianSubMode = modeMixed
	case "bayesian_116":
		o.BayesianSubMode = modeBayes116
	default:
		o.BayesianSubMode = modeRecall
	}
}
