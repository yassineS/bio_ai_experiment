// bcftools mpileup — generate per-position genotype likelihoods from
// BAM input. This is the upstream input to `bcftools call`.
//
// The genotype-likelihood model is the MAQ error model ported from
// bcftools' bam2bcf.c (see errmod.go and bam2bcf.go). For every covered
// reference position mpileup emits one BCF/VCF record carrying:
//
//   - REF and the ALT alleles ordered by coverage-normalised quality
//     sum, always followed by the `<*>` "unseen" symbolic allele.
//   - QUAL fixed at 0 (upstream leaves QUAL for `bcftools call` to set).
//   - INFO/DP (raw read depth), INFO/I16 (the 16-slot calling aux tag),
//     INFO/QS (per-allele quality sums) and INFO/MQ0F.
//   - FORMAT/PL — the multi-allelic phred-scaled genotype likelihoods,
//     one upper-triangle grid of n_alleles*(n_alleles+1)/2 values per
//     sample.
//
// BAQ realignment (slice 3) is wired: mapped reads are run through
// `pkg/htsgo/baq.SamProbRealn` in apply+extend mode before their bases
// enter the pileup, matching upstream's `sam_prob_realn(b, ref, ref_len,
// 3)` call in mpileup.c. `-B/--no-BAQ` disables it and `-E/--redo-BAQ`
// forces recomputation (flag 7). By default upstream enables
// MPLP_REALN | MPLP_REALN_PARTIAL (mpileup.c:1389), so realignment is
// PARTIAL: the per-column has_indel/soft-clip heuristic and the per-read
// spanning check skip reads that do not need BAQ. `-D/--full-BAQ` clears
// MPLP_REALN_PARTIAL (mpileup.c:1567), forcing full BAQ — every read on
// the chromosome is realigned. The port mirrors both modes via
// opts.FullBAQ. For indel-free inputs the two paths coincide.
//
// One faithful-port caveat: upstream's per-column `p->indel` term (an
// indel event adjacent to the column, supplied by the pileup engine) is
// not available without indel detection (slice 4), so for indel-bearing
// inputs the partial heuristic is a slight underestimate.
//
// Deferred slices (see docs/PARITY_ROADMAP.md#bcftools):
//
//   - Slice 4: the bias annotations VDB / SGB / RPBZ / MQBZ / BQBZ /
//     MQSBZ / SCBZ. Records are emitted without those INFO tags.
//   - Indel calling (bam2bcf_indel.c) — every indel knob is accepted at
//     the CLI but inert. The MPLP_REALN_PARTIAL BAQ-skip heuristic is
//     tracked here too: it cannot be faithful without indel detection.
//
// Upstream reference: reference_code/bcftools/mpileup.c (the driver) and
// bam2bcf.c (bcf_call_glfgen / bcf_call_combine / bcf_call2bcf).
package bcftools

import (
	"compress/gzip"
	"container/heap"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Defaults that match upstream bcftools mpileup (mpileup.c:1381-1383).
const (
	// DefaultMpileupMaxDepth is upstream `-d` default.
	DefaultMpileupMaxDepth = 250
	// DefaultMpileupMinMQ is upstream `-q` default.
	DefaultMpileupMinMQ uint8 = 0
	// DefaultMpileupMinBQ is upstream `-Q` default (mpileup.c:1381). The
	// obsolete samtools default of 13 was wrong for bcftools mpileup.
	DefaultMpileupMinBQ uint8 = 1
	// DefaultMpileupMaxBQ is upstream `--max-BQ` default (mpileup.c:1382).
	DefaultMpileupMaxBQ uint8 = 60
	// DefaultMpileupDeltaBQ is upstream `--delta-BQ` default
	// (mpileup.c:1383): a base quality is capped at neighbour_qual+delta.
	DefaultMpileupDeltaBQ = 30
	// DefaultMpileupTandemQual is upstream `-h` default (indel-aware
	// homopolymer penalty). Accepted at CLI but unused.
	DefaultMpileupTandemQual = 500
	// DefaultMpileupExtProb is upstream `--ext-prob`. Unused.
	DefaultMpileupExtProb = 20
	// DefaultMpileupGapFrac is upstream `--gap-frac`. Unused.
	DefaultMpileupGapFrac = 0.05
	// DefaultMpileupOpenProb is upstream `--open-prob`. Unused.
	DefaultMpileupOpenProb = 40
	// DefaultMpileupIndelBias is upstream `--indel-bias`. Unused.
	DefaultMpileupIndelBias = 1.00
	// DefaultMpileupIndelSize is upstream `--indel-size`. Unused.
	DefaultMpileupIndelSize = 110
	// DefaultMpileupMinIReads is upstream `--min-ireads`. The
	// upstream default in mpileup.c:1387 is 2; the value 1 is only
	// applied via the `--indels 1.12` config preset (mpileup.c:1738),
	// which this port does not (yet) honour.
	DefaultMpileupMinIReads = 2
	// DefaultMpileupMaxIDepth is upstream `--max-idepth`. Unused.
	DefaultMpileupMaxIDepth = 250
	// DefaultMpileupARProb is upstream `--ar-prob`. Unused.
	DefaultMpileupARProb = 1e-4
	// mpileupTheta is upstream's CALL_DEFTHETA (bam2bcf.c:39); the errmod
	// depth-correlation parameter is 1 - theta.
	mpileupTheta = 0.83
)

// AmbigReadsMode selects how `--ambig-reads` compensates the per-allele
// AD counts for low-quality REF-looking reads at an indel site. Values
// mirror upstream's B2B_DROP / B2B_INC_AD / B2B_INC_AD0 (bam2bcf.h:81-83).
type AmbigReadsMode int

const (
	// AmbigReadsDrop is the upstream default (`drop`): discard
	// ambiguous reads entirely.
	AmbigReadsDrop AmbigReadsMode = 0
	// AmbigReadsIncAD is `incAD`: distribute ambiguous reads across the
	// per-allele AD slots in proportion to the existing allele counts.
	AmbigReadsIncAD AmbigReadsMode = 1
	// AmbigReadsIncAD0 is `incAD0`: claim every ambiguous read as a
	// REF support, adding to AD slot 0.
	AmbigReadsIncAD0 AmbigReadsMode = 2
)

// parseAmbigReads maps an `--ambig-reads STR` option value to its
// AmbigReadsMode. Comparison is case-insensitive, matching upstream
// (mpileup.c:1779-1782). An empty input maps to AmbigReadsDrop.
func parseAmbigReads(s string) (AmbigReadsMode, error) {
	switch {
	case s == "":
		return AmbigReadsDrop, nil
	case strings.EqualFold(s, "drop"):
		return AmbigReadsDrop, nil
	case strings.EqualFold(s, "incAD"):
		return AmbigReadsIncAD, nil
	case strings.EqualFold(s, "incAD0"):
		return AmbigReadsIncAD0, nil
	}
	return AmbigReadsDrop, fmt.Errorf("the option to --ambig-reads not recognised: %q", s)
}

// FmtFlag bit assignments, mirroring B2B_FMT_*/B2B_INFO_* in
// reference_code/bcftools/bam2bcf.h:46-75. Only the bits that this port
// actually consumes drive emission today; the rest are parsed and
// accepted so `-a` token lists match upstream verbatim.
const (
	B2BFmtDP     uint32 = 1 << 0
	B2BFmtSP     uint32 = 1 << 1
	B2BFmtDV     uint32 = 1 << 2
	B2BFmtDP4    uint32 = 1 << 3
	B2BFmtDPR    uint32 = 1 << 4
	B2BInfoDPR   uint32 = 1 << 5
	B2BFmtAD     uint32 = 1 << 6
	B2BFmtADF    uint32 = 1 << 7
	B2BFmtADR    uint32 = 1 << 8
	B2BInfoAD    uint32 = 1 << 9
	B2BInfoADF   uint32 = 1 << 10
	B2BInfoADR   uint32 = 1 << 11
	B2BInfoSCR   uint32 = 1 << 12
	B2BFmtSCR    uint32 = 1 << 13
	B2BInfoVDB   uint32 = 1 << 14
	B2BFmtQS     uint32 = 1 << 15
	B2BFmtNMBZ   uint32 = 1 << 16
	B2BInfoNMBZ  uint32 = 1 << 17
	B2BInfoBQBZ  uint32 = 1 << 18
	B2BInfoMQBZ  uint32 = 1 << 19
	B2BInfoMQSBZ uint32 = 1 << 20
	B2BInfoRPBZ  uint32 = 1 << 21
	B2BInfoSCBZ  uint32 = 1 << 22
	B2BInfoSGB   uint32 = 1 << 23
	B2BFmtQM     uint32 = 1 << 24
	B2BInfoNM    uint32 = 1 << 25
	B2BInfoMQ0F  uint32 = 1 << 26
	B2BInfoIDV   uint32 = 1 << 27
	B2BInfoIMF   uint32 = 1 << 28
	B2BInfoFS    uint32 = 1 << 29
)

// DefaultMpileupFmtFlag mirrors the bitset assigned in mpileup.c:1399:
// the bias INFO tags that mpileup always emits when a real ALT is
// present, plus the default FORMAT/AD per-sample tag.
const DefaultMpileupFmtFlag = B2BInfoBQBZ | B2BInfoIDV | B2BInfoIMF |
	B2BInfoMQ0F | B2BInfoMQBZ | B2BInfoMQSBZ | B2BInfoRPBZ | B2BInfoSCBZ |
	B2BInfoSGB | B2BInfoVDB | B2BFmtAD

// parseFormatFlag is the Go port of mpileup.c:1141 parse_format_flag. It
// updates *flag with the bits selected by the comma-separated `str`:
// each token is a FORMAT/INFO tag name, optionally prefixed with "-" to
// clear the bit. Both bare names ("AD") and prefixed names ("FORMAT/AD",
// "INFO/AD") are accepted (case-insensitive). An unknown token returns
// an error; callers wire that to a non-zero exit. The upstream warnings
// for deprecated tags (DP4/DPR/DV) are intentionally not emitted — they
// just clutter test output; the bits are still toggled correctly.
func parseFormatFlag(flag *uint32, str string) error {
	if str == "" {
		return nil
	}
	for _, raw := range strings.Split(str, ",") {
		tag := strings.TrimSpace(raw)
		if tag == "" {
			continue
		}
		exclude := false
		if tag[0] == '-' {
			exclude = true
			tag = tag[1:]
		}
		// upper-case keys are compared case-insensitively via strings.EqualFold
		var bit uint32
		switch {
		case mpileupTagMatch(tag, "AD"):
			bit = B2BFmtAD
		case mpileupTagMatch(tag, "ADF"):
			bit = B2BFmtADF
		case mpileupTagMatch(tag, "ADR"):
			bit = B2BFmtADR
		case mpileupTagMatch(tag, "DP"):
			bit = B2BFmtDP
		case mpileupTagMatch(tag, "DP4"):
			bit = B2BFmtDP4
		case mpileupTagMatch(tag, "DPR"):
			bit = B2BFmtDPR
		case mpileupTagMatch(tag, "DV"):
			bit = B2BFmtDV
		case mpileupTagMatch(tag, "NMBZ"):
			bit = B2BFmtNMBZ
		case mpileupTagMatch(tag, "QM"):
			bit = B2BFmtQM
		case mpileupTagMatch(tag, "QS"):
			bit = B2BFmtQS
		case mpileupTagMatch(tag, "SP"):
			bit = B2BFmtSP
		case mpileupTagMatch(tag, "SCR"):
			bit = B2BFmtSCR
		case strings.EqualFold(tag, "INFO/DPR"):
			bit = B2BInfoDPR
		case strings.EqualFold(tag, "INFO/AD"):
			bit = B2BInfoAD
		case strings.EqualFold(tag, "INFO/ADF"):
			bit = B2BInfoADF
		case strings.EqualFold(tag, "INFO/ADR"):
			bit = B2BInfoADR
		case strings.EqualFold(tag, "INFO/BQBZ"):
			bit = B2BInfoBQBZ
		case strings.EqualFold(tag, "INFO/FS"):
			bit = B2BInfoFS
		case strings.EqualFold(tag, "INFO/IDV"):
			bit = B2BInfoIDV
		case strings.EqualFold(tag, "INFO/IMF"):
			bit = B2BInfoIMF
		case strings.EqualFold(tag, "INFO/MQ0F"):
			bit = B2BInfoMQ0F
		case strings.EqualFold(tag, "INFO/MQBZ"):
			bit = B2BInfoMQBZ
		case strings.EqualFold(tag, "INFO/NM"):
			bit = B2BInfoNM
		case strings.EqualFold(tag, "INFO/NMBZ"):
			bit = B2BInfoNMBZ
		case strings.EqualFold(tag, "INFO/RPBZ"):
			bit = B2BInfoRPBZ
		case strings.EqualFold(tag, "INFO/SCBZ"):
			bit = B2BInfoSCBZ
		case strings.EqualFold(tag, "INFO/SCR"):
			bit = B2BInfoSCR
		case strings.EqualFold(tag, "INFO/SGB"):
			bit = B2BInfoSGB
		case strings.EqualFold(tag, "INFO/VDB"):
			bit = B2BInfoVDB
		default:
			return fmt.Errorf("could not parse tag %q in %q", tag, str)
		}
		if exclude {
			*flag &^= bit
		} else {
			*flag |= bit
		}
	}
	return nil
}

// mpileupTagMatch reports whether tag (the bare name supplied by the
// user) names the FORMAT tag fmtName. Bare names ("AD") and the
// "FORMAT/" or "FMT/" prefix ("FORMAT/AD", "FMT/AD") all succeed;
// comparison is case-insensitive, mirroring SET_FMT_FLAG in
// mpileup.c:1120-1122.
func mpileupTagMatch(tag, fmtName string) bool {
	if strings.EqualFold(tag, fmtName) {
		return true
	}
	if strings.EqualFold(tag, "FORMAT/"+fmtName) {
		return true
	}
	if strings.EqualFold(tag, "FMT/"+fmtName) {
		return true
	}
	return false
}

// MpileupOptions configures bcftools mpileup. Fields are 1:1 with the
// upstream getopt_long table in `mpileup.c`. Knobs the model does not
// consume are tagged "accepted; unused" and tracked in PARITY_ROADMAP.
type MpileupOptions struct {
	// Inputs is the list of BAM/SAM paths to pile up. Multi-BAM input
	// yields one sample column per BAM (sample name comes from the @RG
	// SM tag if uniform, otherwise the basename of the input).
	Inputs []string
	// FastaRef is upstream's -f/--fasta-ref. Required: every emitted
	// record needs the REF base.
	FastaRef string
	// BamList is upstream's -b/--bam-list. Files listed one per line
	// (lines starting with '#' and blank lines are ignored) are
	// appended to Inputs.
	BamList string

	// Regions is upstream's -r/--regions (`chr:beg-end[,...]`).
	Regions []string
	// RegionsFile is upstream's -R/--regions-file (BED-like).
	RegionsFile string
	// Targets is upstream's -t/--targets (`chr:beg-end[,...]`).
	Targets []string
	// TargetsFile is upstream's -T/--targets-file (BED-like).
	TargetsFile string

	// Samples is upstream's -s/--samples (comma list).
	Samples []string
	// SamplesFile is upstream's -S/--samples-file.
	SamplesFile string

	// MaxDepth is upstream's -d/--max-depth (default 250).
	MaxDepth int
	// MinMQ is upstream's -q/--min-MQ (default 0).
	MinMQ uint8
	// MinBQ is upstream's -Q/--min-BQ (default 1).
	MinBQ uint8
	// MaxBQ is upstream's --max-BQ cap (default 60).
	MaxBQ uint8
	// DeltaBQ is upstream's --delta-BQ (default 30): a base quality is
	// capped at neighbour_qual+DeltaBQ.
	DeltaBQ int

	// CountOrphans is upstream's -A/--count-orphans.
	CountOrphans bool
	// IgnoreOverlaps is upstream's -x/--ignore-overlaps.
	IgnoreOverlaps bool
	// NoBAQ is upstream's -B/--no-BAQ. When set, BAQ realignment is
	// skipped and raw (delta_baseQ-capped) base qualities are used.
	NoBAQ bool
	// RedoBAQ is upstream's -E/--redo-BAQ. When set, BAQ is recomputed
	// from scratch, discarding any pre-existing BQ tag (baq.FlagRedo).
	RedoBAQ bool
	// FullBAQ is upstream's -D/--full-BAQ. By default mpileup does
	// PARTIAL realignment (upstream's MPLP_REALN_PARTIAL, on by
	// default): the per-column indel/soft-clip heuristic and the
	// per-read spanning check skip reads that do not need BAQ. When
	// FullBAQ is set, -D clears MPLP_REALN_PARTIAL so every read on the
	// chromosome is BAQ-realigned ("full BAQ").
	FullBAQ bool
	// AdjustMQ is upstream's -C/--adjust-mq. Accepted; ignored.
	AdjustMQ int

	// Annotate is upstream's -a/--annotate list (FORMAT/INFO tags to
	// include). Tokens are parsed into FmtFlag by validateMpileupOptions
	// (see parseFormatFlag, the Go port of mpileup.c:1141 parse_format_flag);
	// a "-TAG" token clears the bit, mirroring upstream. The default set
	// (mpileup.c:1399) is BQBZ/IDV/IMF/MQ0F/MQBZ/MQSBZ/RPBZ/SCBZ/SGB/VDB
	// + FORMAT/AD.
	Annotate string
	// FmtFlag is the resolved bitset selected by Annotate; the
	// B2BFmtAD/B2BInfoAD/... constants are the bit assignments. It is
	// populated by validateMpileupOptions from Annotate (and the upstream
	// default). External callers can also set it directly.
	FmtFlag uint32

	// ReadGroups is upstream's -G/--read-groups. Accepted; ignored.
	ReadGroups string
	// IgnoreRG is upstream's --ignore-RG (long-only). Accepted; ignored.
	IgnoreRG bool

	// Platforms is upstream's -P/--platforms. Accepted; ignored.
	Platforms string

	// Config is upstream's -X/--config (predefined indel-model preset).
	// Accepted; ignored (no indel realigner).
	Config string

	// PerSampleMF is upstream's -p/--per-sample-mF. Accepted; ignored.
	PerSampleMF bool

	// Seed is upstream's --seed (random seed for subsampling).
	// Accepted; ignored (no subsampling).
	Seed int64

	// TandemQual is upstream's -h/--tandem-qual. Accepted; ignored.
	TandemQual int
	// ExtProb is upstream's --ext-prob. Accepted; ignored.
	ExtProb int
	// GapFrac is upstream's --gap-frac. Accepted; ignored.
	GapFrac float64
	// OpenProb is upstream's --open-prob. Accepted; ignored.
	OpenProb int
	// IndelBias is upstream's --indel-bias. Accepted; ignored.
	IndelBias float64
	// IndelSize is upstream's --indel-size. Accepted; ignored.
	IndelSize int
	// MinIReads is upstream's --min-ireads. Accepted; ignored.
	MinIReads int
	// MaxIDepth is upstream's --max-idepth. Accepted; ignored.
	MaxIDepth int
	// ARProb is upstream's --ar-prob. Accepted; ignored.
	ARProb float64
	// AmbigReads is upstream's --ambig-reads / --ar string form.
	// validateMpileupOptions parses it into AmbigReadsMode.
	AmbigReads string
	// AmbigReadsMode is the parsed --ambig-reads selection that drives
	// the indel-branch glfgen's ADR/ADF compensation. Callers can also
	// set it directly when the string form is left empty.
	AmbigReadsMode AmbigReadsMode
	// MaxReadLen is upstream's -M/--max-read-len. Accepted; ignored.
	MaxReadLen int

	// DelBias is upstream's --del-bias (hidden). Accepted; ignored.
	DelBias float64
	// PolyMQual is upstream's --poly-mqual. Accepted; ignored.
	PolyMQual bool
	// ScoreVsRef is upstream's --score-vs-ref. Accepted; ignored.
	ScoreVsRef float64
	// SeqQOffset is upstream's --seqq-offset. Accepted; ignored.
	SeqQOffset int

	// SkipIndels is upstream's -I/--skip-indels. mpileup never emits
	// indel records yet so the flag is effectively the default.
	SkipIndels bool
	// IndelsCNS enables upstream's --indels-cns / --indels-2.0
	// consensus-based indel caller (Go port at bam2bcf_indelcns.go).
	// When false the legacy probabilistic caller
	// (bam2bcf_indel_align.go) is used; this is the default.
	IndelsCNS bool
	// NoIndelsCNS is upstream's --no-indels-cns. It is the inverse of
	// IndelsCNS; the CLI wiring resolves which takes precedence.
	NoIndelsCNS bool

	// GVCFBlock is upstream's -g/--gvcf. Accepted; one record per
	// covered position is always emitted (no gVCF blocking yet).
	GVCFBlock string

	// NoReference is upstream's --no-reference (skip the FASTA REF
	// check). Accepted; the FASTA REF is always used.
	NoReference bool

	// OutputFormat is upstream's -O/--output-type (v|z|u|b).
	OutputFormat OutputFormat
	// Output is upstream's -o/--output (default stdout).
	Output string
	// CompressLevel is upstream's --compression-level (gzip level for -O z).
	CompressLevel int

	// Threads is upstream's --threads (accepted; single-threaded).
	Threads int
	// NoVersion is upstream's --no-version (omit the version line).
	NoVersion bool

	// Verbosity is upstream's -v/--verbosity (accepted; ignored).
	Verbosity int

	// FlagIncl / FlagExcl / FlagAny / FlagLS are the raw user-supplied
	// string forms of upstream's --rf/--ff/--ls/etc. They are parsed by
	// validateMpileupOptions into the typed RflagSkip* masks below.
	//
	//   - FlagIncl  is upstream's --rf / --lu / --skip-all-unset
	//     (mpileup.c:1413+1418): skip reads with all of these bits unset.
	//   - FlagExcl  is upstream's --ff / --ns / --skip-any-set
	//     (mpileup.c:1415+1419): skip reads with any of these bits set.
	//   - FlagAny   is upstream's --nu / --skip-any-unset
	//     (mpileup.c:1416): skip reads with any of these bits unset.
	//   - FlagLS    is upstream's --ls / --skip-all-set
	//     (mpileup.c:1419): skip reads with all of these bits set.
	FlagIncl string
	FlagExcl string
	FlagAny  string
	FlagLS   string

	// RflagSkipAnyUnset / RflagSkipAllUnset / RflagSkipAnySet /
	// RflagSkipAllSet are the parsed BAM-flag masks driving the
	// upstream `mplp_func` per-read filters (mpileup.c:208-211).
	// validateMpileupOptions populates them from the corresponding
	// Flag* string fields; callers can also set them directly. Both
	// `0` and "no bits" are no-op (matching upstream's `if (mask)` gates).
	RflagSkipAnyUnset uint16
	RflagSkipAllUnset uint16
	RflagSkipAnySet   uint16
	RflagSkipAllSet   uint16
}

// mpileupBAQFlag derives the realn flag passed to baq.SamProbRealn from
// the mpileup options, mirroring mpileup.c:548
// `sam_prob_realn(b, ref, ref_len, (flag & MPLP_REDO_BAQ) ? 7 : 3)`.
// mpileup always realigns in apply+extend mode (3); -E adds FlagRedo (7).
func mpileupBAQFlag(opts MpileupOptions) int {
	flag := baq.FlagApply | baq.FlagExtend
	if opts.RedoBAQ {
		flag |= baq.FlagRedo
	}
	return flag
}

// MpileupFile is the file-path entry point. It opens every input BAM,
// the FASTA reference, and writes BCF or VCF to out.
func MpileupFile(opts MpileupOptions, out io.Writer) error {
	if err := validateMpileupOptions(&opts); err != nil {
		return err
	}
	inputs, err := resolveMpileupInputs(opts)
	if err != nil {
		return err
	}
	if len(inputs) == 0 {
		return fmt.Errorf("bcftools mpileup: no input BAM files")
	}
	if opts.FastaRef == "" {
		return fmt.Errorf("bcftools mpileup: -f/--fasta-ref is required")
	}
	ref, err := fasta.OpenRandomAccess(opts.FastaRef)
	if err != nil {
		return fmt.Errorf("bcftools mpileup: open reference: %w", err)
	}
	defer ref.Close()

	// Open every BAM input and bind to a sam.Reader.
	type input struct {
		path   string
		file   *os.File
		reader sam.Reader
		sample string
	}
	in := make([]input, 0, len(inputs))
	defer func() {
		for _, x := range in {
			if x.file != nil {
				_ = x.file.Close()
			}
		}
	}()
	for _, p := range inputs {
		f, err := os.Open(p)
		if err != nil {
			return fmt.Errorf("bcftools mpileup: %w", err)
		}
		rd, err := sam.NewReader(f)
		if err != nil {
			_ = f.Close()
			return fmt.Errorf("bcftools mpileup: %s: %w", p, err)
		}
		in = append(in, input{path: p, file: f, reader: rd, sample: deriveSample(rd, p)})
	}

	// Pull records bucketed by chrom for every input.
	perInputRecs := make([]map[string][]*sam.Record, len(in))
	for i, x := range in {
		recs, err := mpileupReadBAM(x.reader, opts)
		if err != nil {
			return fmt.Errorf("bcftools mpileup: %s: %w", x.path, err)
		}
		perInputRecs[i] = recs
	}

	// Resolve the chromosome iteration order: prefer the first input's header.
	hdr0 := in[0].reader.Header()
	chromOrder := make([]string, 0, len(hdr0.Refs))
	chromLen := make(map[string]int, len(hdr0.Refs))
	for _, r := range hdr0.Refs {
		chromOrder = append(chromOrder, r.Name)
		chromLen[r.Name] = int(r.Length)
	}

	// Parse region/target windows. -r and -t are both treated as
	// post-filters (no BAI seek path).
	regWindows, err := parseMpileupRegions(opts, chromLen)
	if err != nil {
		return err
	}

	// Sample names for the #CHROM line and FORMAT column.
	samples := make([]string, len(in))
	for i, x := range in {
		samples[i] = x.sample
	}
	if len(opts.Samples) > 0 || opts.SamplesFile != "" {
		want := map[string]struct{}{}
		for _, s := range opts.Samples {
			want[s] = struct{}{}
		}
		if opts.SamplesFile != "" {
			names, err := LoadSamplesFile(opts.SamplesFile)
			if err != nil {
				return fmt.Errorf("bcftools mpileup: %w", err)
			}
			for _, s := range names {
				want[s] = struct{}{}
			}
		}
		keep := in[:0]
		keepRecs := perInputRecs[:0]
		keepSamp := samples[:0]
		for i, x := range in {
			if _, ok := want[x.sample]; !ok {
				continue
			}
			keep = append(keep, x)
			keepRecs = append(keepRecs, perInputRecs[i])
			keepSamp = append(keepSamp, x.sample)
		}
		in = keep
		perInputRecs = keepRecs
		samples = keepSamp
	}

	return writeMpileupVCF(out, opts, ref, chromOrder, chromLen, perInputRecs, samples, regWindows)
}

// validateMpileupOptions applies upstream's defaults (mpileup.c:1381-1383).
// It also parses opts.Annotate into opts.FmtFlag, layering the
// user-specified `-a` tokens on top of DefaultMpileupFmtFlag (the
// "* "-marked default annotation set in mpileup.c:1399). When FmtFlag
// has already been set non-zero by the caller, it is taken as the
// starting bitset and only the Annotate tokens further modify it.
func validateMpileupOptions(opts *MpileupOptions) error {
	if opts.MaxDepth == 0 {
		opts.MaxDepth = DefaultMpileupMaxDepth
	}
	if opts.MinBQ == 0 {
		opts.MinBQ = DefaultMpileupMinBQ
	}
	if opts.MaxBQ == 0 {
		opts.MaxBQ = DefaultMpileupMaxBQ
	}
	if opts.DeltaBQ == 0 {
		opts.DeltaBQ = DefaultMpileupDeltaBQ
	}
	if opts.FmtFlag == 0 {
		opts.FmtFlag = DefaultMpileupFmtFlag
	}
	if err := parseFormatFlag(&opts.FmtFlag, opts.Annotate); err != nil {
		return fmt.Errorf("bcftools mpileup: -a/--annotate: %w", err)
	}
	// --ambig-reads: parse the string form when AmbigReadsMode hasn't
	// been set directly. An empty string means "stay with default
	// AmbigReadsDrop".
	if opts.AmbigReads != "" && opts.AmbigReadsMode == AmbigReadsDrop {
		mode, err := parseAmbigReads(opts.AmbigReads)
		if err != nil {
			return fmt.Errorf("bcftools mpileup: %w", err)
		}
		opts.AmbigReadsMode = mode
	}
	// --skip-* BAM-flag masks. The string form is parsed via
	// parseBAMFlagString; either bare integers (e.g. "0x14", "20") or
	// comma-separated names ("PAIRED,DUP") are accepted, mirroring
	// htslib's bam_str2flag (sam.c:5290).
	if v, err := parseBAMFlagString(opts.FlagAny); err != nil {
		return fmt.Errorf("bcftools mpileup: --skip-any-unset/--nu: %w", err)
	} else if opts.RflagSkipAnyUnset == 0 {
		opts.RflagSkipAnyUnset = v
	}
	if v, err := parseBAMFlagString(opts.FlagIncl); err != nil {
		return fmt.Errorf("bcftools mpileup: --skip-all-unset/--rf/--lu: %w", err)
	} else if opts.RflagSkipAllUnset == 0 {
		opts.RflagSkipAllUnset = v
	}
	if v, err := parseBAMFlagString(opts.FlagExcl); err != nil {
		return fmt.Errorf("bcftools mpileup: --skip-any-set/--ff/--ns: %w", err)
	} else if opts.RflagSkipAnySet == 0 {
		opts.RflagSkipAnySet = v
	}
	if v, err := parseBAMFlagString(opts.FlagLS); err != nil {
		return fmt.Errorf("bcftools mpileup: --skip-all-set/--ls: %w", err)
	} else if opts.RflagSkipAllSet == 0 {
		opts.RflagSkipAllSet = v
	}
	return nil
}

// parseBAMFlagString is the Go port of htslib's bam_str2flag
// (sam.c:5290). An empty string returns 0 (no-op). A bare integer
// (decimal, hex with 0x prefix, octal with 0 prefix) is parsed
// directly. Otherwise the input is split on commas and each token is
// matched case-insensitively against the canonical flag names; an
// unknown token returns an error.
func parseBAMFlagString(s string) (uint16, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return 0, nil
	}
	// Try integer first (decimal/hex/octal via base 0).
	if v, err := strconv.ParseUint(s, 0, 32); err == nil {
		return uint16(v), nil
	}
	var out uint16
	for _, raw := range strings.Split(s, ",") {
		tok := strings.TrimSpace(raw)
		if tok == "" {
			continue
		}
		switch strings.ToUpper(tok) {
		case "PAIRED":
			out |= sam.FlagPaired
		case "PROPER_PAIR":
			out |= sam.FlagProperPair
		case "UNMAP":
			out |= sam.FlagUnmapped
		case "MUNMAP":
			out |= sam.FlagMateUnmapped
		case "REVERSE":
			out |= sam.FlagReverse
		case "MREVERSE":
			out |= sam.FlagMateReverse
		case "READ1":
			out |= sam.FlagRead1
		case "READ2":
			out |= sam.FlagRead2
		case "SECONDARY":
			out |= sam.FlagSecondary
		case "QCFAIL":
			out |= sam.FlagQCFail
		case "DUP":
			out |= sam.FlagDuplicate
		case "SUPPLEMENTARY":
			out |= sam.FlagSupplementary
		default:
			return 0, fmt.Errorf("could not parse flag name %q in %q", tok, s)
		}
	}
	return out, nil
}

// resolveMpileupInputs reads -b/--bam-list (when given) and appends to
// the explicit Inputs slice. Order is inputs-then-list-file, matching
// upstream.
func resolveMpileupInputs(opts MpileupOptions) ([]string, error) {
	out := append([]string{}, opts.Inputs...)
	if opts.BamList == "" {
		return out, nil
	}
	f, err := os.Open(opts.BamList)
	if err != nil {
		return nil, fmt.Errorf("bcftools mpileup: open bam-list: %w", err)
	}
	defer f.Close()
	b, err := io.ReadAll(f)
	if err != nil {
		return nil, fmt.Errorf("bcftools mpileup: read bam-list: %w", err)
	}
	for _, ln := range strings.Split(string(b), "\n") {
		ln = strings.TrimSpace(ln)
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		out = append(out, ln)
	}
	return out, nil
}

// deriveSample picks a sample name for the output. We use the @RG SM
// tag when uniform across the file's @RG lines, falling back to the
// basename of the BAM.
func deriveSample(rd sam.Reader, path string) string {
	hdr := rd.Header()
	var sm string
	for _, rg := range hdr.ReadGroups {
		var rgSM string
		for _, f := range rg.Extra {
			if f.Tag == "SM" {
				rgSM = f.Value
				break
			}
		}
		if rgSM == "" {
			continue
		}
		if sm == "" {
			sm = rgSM
			continue
		}
		if rgSM != sm {
			sm = ""
			break
		}
	}
	if sm != "" {
		return sm
	}
	base := filepath.Base(path)
	if ext := filepath.Ext(base); ext == ".bam" || ext == ".sam" {
		base = strings.TrimSuffix(base, ext)
	}
	return base
}

// mpileupReadBAM pulls every record from rd, applies upstream's
// record-level filters (flag bits, MAPQ), and buckets by RName.
func mpileupReadBAM(rd sam.Reader, opts MpileupOptions) (map[string][]*sam.Record, error) {
	out := map[string][]*sam.Record{}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if !mpileupKeepRecord(rec, opts) {
			continue
		}
		out[rec.RName] = append(out[rec.RName], rec)
	}
	for k := range out {
		sort.SliceStable(out[k], func(i, j int) bool { return out[k][i].Pos < out[k][j].Pos })
	}
	return out, nil
}

// mpileupKeepRecord applies upstream's default read-level filters
// (unmapped, secondary, QCfail, duplicate; orphans unless -A; MAPQ
// floor) and the user-supplied --skip-* BAM-flag masks. The default
// mask (mpileup.c:1392 BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)
// is layered through RflagSkipAnySet by validateMpileupOptions so
// upstream's mplp_func per-read checks (mpileup.c:208-211) drive
// the same predicate.
func mpileupKeepRecord(rec *sam.Record, opts MpileupOptions) bool {
	if rec == nil || rec.Pos <= 0 || rec.RName == "" {
		return false
	}
	if rec.Flag&(sam.FlagUnmapped|sam.FlagSecondary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
		return false
	}
	// --skip-* user filters (mpileup.c:208-211). The `if (mask)` gates
	// match upstream: a zero mask is a no-op, so callers who do not
	// set the flag get default behaviour.
	if m := opts.RflagSkipAnyUnset; m != 0 && (m&rec.Flag) != m {
		return false
	}
	if m := opts.RflagSkipAllSet; m != 0 && (m&rec.Flag) == m {
		return false
	}
	if m := opts.RflagSkipAllUnset; m != 0 && (m&rec.Flag) == 0 {
		return false
	}
	if m := opts.RflagSkipAnySet; m != 0 && (m&rec.Flag) != 0 {
		return false
	}
	if !opts.CountOrphans && rec.Flag&sam.FlagPaired != 0 {
		if rec.Flag&sam.FlagMateUnmapped != 0 {
			return false
		}
		if rec.Flag&sam.FlagProperPair == 0 {
			return false
		}
	}
	if opts.MinMQ > 0 && rec.MapQ < opts.MinMQ {
		return false
	}
	return true
}

// parseMpileupRegions resolves -r/-R/-t/-T into a flat per-chrom list of
// 1-based inclusive windows. Empty result means "no restriction".
func parseMpileupRegions(opts MpileupOptions, chromLen map[string]int) (map[string][][2]int, error) {
	var specs []string
	specs = append(specs, opts.Regions...)
	specs = append(specs, opts.Targets...)
	if opts.RegionsFile != "" {
		extra, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return nil, fmt.Errorf("bcftools mpileup: %w", err)
		}
		specs = append(specs, extra...)
	}
	if opts.TargetsFile != "" {
		extra, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return nil, fmt.Errorf("bcftools mpileup: %w", err)
		}
		specs = append(specs, extra...)
	}
	if len(specs) == 0 {
		return nil, nil
	}
	out := map[string][][2]int{}
	for _, s := range specs {
		chrom, beg, end, err := parseMpileupRegionSpec(s, chromLen)
		if err != nil {
			return nil, fmt.Errorf("bcftools mpileup: %w", err)
		}
		out[chrom] = append(out[chrom], [2]int{beg, end})
	}
	for k, iv := range out {
		sort.Slice(iv, func(i, j int) bool { return iv[i][0] < iv[j][0] })
		merged := iv[:0]
		cur := iv[0]
		for i := 1; i < len(iv); i++ {
			if iv[i][0] <= cur[1]+1 {
				if iv[i][1] > cur[1] {
					cur[1] = iv[i][1]
				}
				continue
			}
			merged = append(merged, cur)
			cur = iv[i]
		}
		merged = append(merged, cur)
		out[k] = merged
	}
	return out, nil
}

// parseMpileupRegionSpec parses a single chr[:beg[-end]] spec into
// 1-based inclusive coordinates.
func parseMpileupRegionSpec(s string, chromLen map[string]int) (chrom string, beg, end int, err error) {
	colon := strings.IndexByte(s, ':')
	if colon < 0 {
		chrom = s
		beg = 1
		end = chromLen[chrom]
		if end == 0 {
			end = 1 << 30
		}
		return
	}
	chrom = s[:colon]
	rest := s[colon+1:]
	dash := strings.IndexByte(rest, '-')
	if dash < 0 {
		beg, err = strconv.Atoi(rest)
		if err != nil {
			return "", 0, 0, fmt.Errorf("bad region %q: %w", s, err)
		}
		end = beg
		return
	}
	beg, err = strconv.Atoi(rest[:dash])
	if err != nil {
		return "", 0, 0, fmt.Errorf("bad region %q: %w", s, err)
	}
	tail := rest[dash+1:]
	if tail == "" {
		end = chromLen[chrom]
		if end == 0 {
			end = 1 << 30
		}
		return
	}
	end, err = strconv.Atoi(tail)
	if err != nil {
		return "", 0, 0, fmt.Errorf("bad region %q: %w", s, err)
	}
	return
}

// regionContains returns true when 1-based pos is inside any of the
// windows associated with chrom. When windows is nil every position
// passes.
func regionContains(windows map[string][][2]int, chrom string, pos1 int) bool {
	if windows == nil {
		return true
	}
	iv, ok := windows[chrom]
	if !ok {
		return false
	}
	for _, r := range iv {
		if pos1 >= r[0] && pos1 <= r[1] {
			return true
		}
	}
	return false
}

// writeMpileupVCF walks every chrom in chromOrder, gathers per-position
// pileup columns from every input, runs the glfgen/combine/2bcf
// pipeline, and writes one record per covered position to out.
func writeMpileupVCF(out io.Writer, opts MpileupOptions, ref *fasta.RandomAccess,
	chromOrder []string, chromLen map[string]int,
	perInputRecs []map[string][]*sam.Record, samples []string,
	regWindows map[string][][2]int) error {

	hdr := buildMpileupHeader(opts, chromOrder, chromLen, samples)
	w, finish, err := openMpileupOutput(out, opts, hdr)
	if err != nil {
		return err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return err
	}

	// The errmod tables are expensive to build, so do it once.
	em := errmod.Init(1.0 - mpileupTheta)

	// Run-level BQBZ / MQSBZ leak (bam2bcf.c:1175-1183 + the lack of a
	// reset in bcf_callaux_clean). Upstream's bcf_call_t is allocated
	// once for the whole run; mwu_bq / mwu_mqs are C floats default-
	// initialised to 0.0, then OVERWRITTEN only by has-alt SNP combines.
	// The indel-branch combine reads whatever the most recent SNP-combine
	// value was, regardless of the chromosome boundary. We mirror that by
	// threading a single biasLeak instance through every chromosome,
	// pre-initialised to "value 0, ok=true" so the very first indel record
	// (before any has-alt SNP combine) sees BQBZ=MQSBZ=0 instead of being
	// silently dropped.
	leak := biasLeak{bq: 0, mqs: 0, bqOK: true, mqsOK: true}

	for _, chrom := range chromOrder {
		if regWindows != nil {
			if _, ok := regWindows[chrom]; !ok {
				continue
			}
		}
		refLen := chromLen[chrom]
		if refLen <= 0 {
			continue
		}
		anyHit := false
		perInputChromRecs := make([][]*sam.Record, len(perInputRecs))
		for i, recs := range perInputRecs {
			perInputChromRecs[i] = recs[chrom]
			if len(perInputChromRecs[i]) > 0 {
				anyHit = true
			}
		}
		if !anyHit {
			continue
		}
		refSlab, err := ref.Fetch(chrom, 0, int64(refLen))
		if err != nil {
			return fmt.Errorf("bcftools mpileup: fetch %s: %w", chrom, err)
		}
		if err := emitChromMpileup(w, em, chrom, refSlab, refLen, perInputChromRecs, opts, regWindows, &leak); err != nil {
			return err
		}
	}
	return w.Flush()
}

// emitChromMpileup walks every covered position on one chromosome and
// writes one record per position that has read coverage. Unlike the
// pre-MAQ port, this emits a record for every covered position (not
// only SNP candidates) with `<*>` as the unseen allele.
func emitChromMpileup(w variantWriter, em *errmod.Errmod, chrom string, refSlab []byte, refLen int,
	perInputChromRecs [][]*sam.Record, opts MpileupOptions,
	regWindows map[string][][2]int, leak *biasLeak) error {

	nIn := len(perInputChromRecs)
	// -d/--max-depth: htslib applies the cap per alignment-start inside
	// bam_plp_push (sam.c:6090), dropping a read when iter->pos ==
	// b->core.pos and the queue already holds maxcnt active reads.
	// applyMpileupDepthCap simulates that on the coordinate-sorted
	// record stream so deeply-piled columns (homopolymers, PCR-dup
	// stacks) drop the same reads upstream does. The cap runs BEFORE
	// the smart-overlap tweak so dropped reads do not leak quality
	// edits into their surviving mates: upstream's overlap_push runs
	// inside bam_plp_push, *after* the cap test, so reads dropped by
	// the cap never reach the overlap-quality merger.
	applyMpileupDepthCap(perInputChromRecs, opts.MaxDepth)

	// MPLP_SMART_OVERLAPS + BAQ. Upstream htslib's bam_plp_push
	// (reference_code/htslib/sam.c:6083-6132) interleaves these per
	// read: when a mapped read enters the pileup queue, overlap_push
	// runs only if its mate is already queued (i.e. when the second
	// mate arrives); BAQ (mplp_realn) then runs at the read's first
	// eligible pileup column. The upshot is:
	//   * mate 1's BAQ runs at its first eligible column; if that
	//     column precedes the second mate's arrival, BAQ sees raw
	//     quals (overlap_push has not run yet) — otherwise BAQ sees
	//     the already-merged quals (overlap_push ran when the second
	//     mate was pushed);
	//   * mate 2's BAQ runs after its own arrival triggers
	//     overlap_push, so it always sees merged quals.
	// The previous batched ordering ran applySmartOverlaps for the
	// entire chromosome before applyMpileupBAQ, which fed every mate
	// 1's BAQ overlap-merged quals — driving the indel-AD.1.out
	// cluster-2 SNP-row I16 BQ drifts at 000000F:446-624. The split
	// below reproduces upstream's per-read ordering: phase 1 BAQs
	// standalones + first-mates whose trigger column precedes the
	// mate's start (raw quals); applySmartOverlaps merges; phase 2
	// BAQs second-mates + any remaining first-mates (merged quals).
	// preMergeQual carries, for each first-mate of an overlapping pair, a
	// snapshot of rec.Qual taken AFTER any pre-merge BAQ but BEFORE
	// applySmartOverlaps mutates the array. accumulateMpileupBases consults
	// it when emitting prevQ/nextQ for columns that upstream's bam_plp_next
	// would have drained BEFORE bam_plp_push pushed the second mate and
	// fired overlap_push (sam.c:5970-5980) — at those columns the
	// delta_baseQ neighbour cap in bcfCallGlfgenCore (bam2bcf.c:428-435)
	// reads raw qual[qpos±1], not the post-merge zeros. Without this
	// snapshot our chromosome-wide applySmartOverlaps zeroes the future
	// neighbour bytes early, the cap over-clips, and the SNP-row I16 BQ
	// sums drift — see the indel-AD.1.out cluster-2 residual (000000F:450
	// and :500). preMergeDrain maps a first-mate record to its drain
	// threshold (max BAM-order intermediate read pos between the two
	// mates) so the snapshot is consulted at column C iff C < threshold.
	var preMergeQual map[*sam.Record][]byte
	var preMergeDrain map[*sam.Record]int
	if opts.NoBAQ {
		// -B disables BAQ entirely. Overlap-merge alone matches the
		// upstream MPLP_NO_BAQ branch (mpileup.c bypasses mplp_realn
		// while keeping bam_plp_push's overlap_push call).
		if !opts.IgnoreOverlaps {
			pairs := classifyMatePairs(perInputChromRecs)
			preMergeQual, preMergeDrain = snapshotFirstMateQuals(perInputChromRecs, pairs)
			applySmartOverlaps(perInputChromRecs)
		}
	} else if opts.IgnoreOverlaps {
		// -x disables overlap-merge. BAQ runs once per read on raw
		// quals — equivalent to the previous batched form. No
		// pre-merge snapshot is needed because overlap_push never
		// fires upstream either.
		applyMpileupBAQ(perInputChromRecs, refSlab, opts, nil, nil)
	} else {
		// Default: per-pair interleaving. Classify mates once, then
		// run BAQ→overlap→BAQ as upstream does. Eligibility uses
		// the trigger column (pos0) to decide whether overlap_push
		// would have run by the time the read's first eligible BAQ
		// column is reached: for first-mates that's pos0 >=
		// mate.Start; for second-mates the overlap always precedes
		// BAQ (their push runs overlap_push at or before their first
		// eligible column). realigned is a shared dedup set so
		// reads BAQ'd in phase 1 are not re-BAQ'd in phase 2.
		pairs := classifyMatePairs(perInputChromRecs)
		// Compute the per-first-mate drain threshold BEFORE phase 1
		// so the phase-1 predicate can mirror upstream's bam_plp_push
		// timing. Upstream drains a column C with the mate already
		// merged when no intermediate read with pos strictly between
		// F.Pos and mate.Pos has been pushed; drainAt captures that
		// timing per read (= max(intermediate.Pos), init F.Pos). The
		// drain threshold depends only on push (= coordinate) order,
		// not on quals, so it is safe to compute pre-BAQ.
		drainAt := computeFirstMateDrainThresholds(perInputChromRecs, pairs)
		realigned := make(map[*sam.Record]bool)
		// Phase 1: standalones + first-mates whose trigger column is
		// strictly before their drain threshold (raw quals). Reads
		// with no drainAt entry (second mates, non-overlapping
		// first-mates, standalones) fall through to the !ok branch
		// above. For first-mates with no intermediate strictly
		// between F.Pos and mate.Pos, drainAt == F.Pos and the gate
		// is always false → phase 1 skips them and phase 2 BAQs the
		// kept mate on merged quals, matching upstream's
		// bam_plp_push (sam.c:6083-6132) firing overlap_push before
		// draining col F.Pos with mplp_realn (mpileup.c:573).
		applyMpileupBAQ(perInputChromRecs, refSlab, opts, func(rec *sam.Record, pos0 int) bool {
			p, ok := pairs[rec]
			if !ok {
				return true
			}
			if p.class != mateClassFirst {
				return false
			}
			d, dOk := drainAt[rec]
			if !dOk {
				d = p.mateStart
			}
			return pos0 < d
		}, realigned)
		// Snapshot first-mate quals AFTER phase 1 BAQ (which mirrors
		// upstream's mplp_realn run at the read's first eligible
		// column) but BEFORE overlap-merge zeroes the future
		// neighbour bytes. This is the same quality array upstream's
		// bcf_call_glfgen reads as qual[qpos±1] at any column C
		// drained before the second mate's push (see
		// snapshotFirstMateQuals for the drain-threshold derivation).
		preMergeQual, preMergeDrain = snapshotFirstMateQuals(perInputChromRecs, pairs)
		applySmartOverlaps(perInputChromRecs)
		// Phase 2: all second-mates + the first-mates phase 1 left
		// untouched (their trigger column lies on or after the drain
		// threshold, so by upstream's bam_plp_push ordering
		// overlap_push has already merged their quals by the time
		// BAQ runs).
		applyMpileupBAQ(perInputChromRecs, refSlab, opts, func(rec *sam.Record, pos0 int) bool {
			p, ok := pairs[rec]
			if !ok {
				return false
			}
			if p.class == mateClassSecond {
				return true
			}
			d, dOk := drainAt[rec]
			if !dOk {
				d = p.mateStart
			}
			return pos0 >= d
		}, realigned)
	}

	// Upstream's pileup engine walks reads' CIGARs and emits a column
	// for every covered reference position, including positions past
	// the FASTA end (those rows carry REF=N). Compute the effective
	// chromosome length as max(refLen, maxReadEnd) so the events array
	// can hold those trailing columns. EndPosition is 1-based inclusive,
	// which equals the 0-based exclusive end matching upstream's
	// bam_endpos.
	effLen := refLen
	for i := 0; i < nIn; i++ {
		for _, rec := range perInputChromRecs[i] {
			if e := int(rec.EndPosition()); e > effLen {
				effLen = e
			}
		}
	}

	// events[input][pos0] is the pileup column for one input at one
	// reference position. The preMergeQual / preMergeDrain side-maps
	// thread the pre-overlap-merge quality snapshot for first-mates into
	// accumulateMpileupBases so the per-column prevQ/nextQ delta_baseQ
	// neighbours match the live qual[qpos±1] upstream's bcf_call_glfgen
	// would read at the same column. Nil maps degrade to "always read
	// rec.Qual" — the correct behaviour with -x/--ignore-overlaps and for
	// non-first-mate reads.
	events := make([][][]pileupBase, nIn)
	for i := 0; i < nIn; i++ {
		events[i] = make([][]pileupBase, effLen)
		for _, rec := range perInputChromRecs[i] {
			accumulateMpileupBases(rec, events[i], preMergeQual, preMergeDrain)
		}
	}

	calls := make([]bcfCallret, nIn)
	// Indel-pass state. bca is the per-call aux struct; leak threads the
	// BQBZ / MQSBZ scalars from the last has_alt SNP combine into the
	// indel combine (upstream's bcf_call_t reuse — only the bca arrays
	// are cleaned between SNP and indel passes).
	bca := newBcfCallauxIndel(opts)
	bca.Chr = chrom
	indelCalls := make([]bcfCallret, nIn)
	piles := make([][]pileupBase, nIn)
	for pos0 := 0; pos0 < effLen; pos0++ {
		pos1 := pos0 + 1
		if !regionContains(regWindows, chrom, pos1) {
			continue
		}
		refB := byte('N')
		if pos0 < len(refSlab) {
			refB = upperByte(refSlab[pos0])
		}
		ref4 := seqNt16Int[baseToNt16(refB)]

		// Per-sample glfgen. Track total coverage so all-empty
		// positions are skipped.
		anyCov := false
		for i := 0; i < nIn; i++ {
			piles[i] = filterMpileupPile(events[i][pos0])
			if len(piles[i]) > 0 {
				anyCov = true
			}
			bcfCallGlfgen(piles[i], ref4, opts, em, &calls[i])
		}
		if !anyCov {
			continue
		}
		call := bcfCallCombine(calls, ref4)
		v := bcfCall2bcf(chrom, pos1, refB, &call, opts.FmtFlag)
		if err := w.Write(v); err != nil {
			return err
		}
		// Update the SNP→indel bias leak: only has_alt SNP combines
		// overwrite BQBZ / MQSBZ on the shared bcf_call_t.
		leak.update(&call)

		// Indel pass (mpileup.c:589-613). Upstream gates this on
		// `total_depth < max_indel_depth`, which our port treats as
		// always-true (MaxIDepth is accepted but unused so far). Skip
		// only when -I/--skip-indels is in force.
		if opts.SkipIndels {
			continue
		}
		// Past the FASTA end there is no reference to anchor an indel
		// call (upstream's indel pass reads ref_fai only within the
		// FASTA), so the N-REF SNP column is emitted alone.
		if pos0 >= refLen {
			continue
		}
		var iret int
		if opts.IndelsCNS {
			// --indels-cns / --indels-2.0: consensus-based caller via
			// the in-tree edlib engine (bam2bcf_indelcns.go, port of
			// reference_code/bcftools/bam2bcf_edlib.c).
			iret = bcfCallGapPrepCNS(piles, pos0, bca, refSlab)
		} else {
			iret = bcfCallGapPrep(piles, pos0, bca, refSlab)
		}
		if iret < 0 {
			continue
		}
		// Per-sample indel-branch glfgen.
		for i := 0; i < nIn; i++ {
			bcfCallGlfgenIndel(piles[i], opts, em, &indelCalls[i])
		}
		icall, ok := bcfCallCombineIndel(indelCalls, calls, bca, *leak)
		if !ok {
			continue
		}
		iv := bcfCall2bcfIndel(chrom, pos0, refSlab, &icall, bca, opts.FmtFlag)
		if err := w.Write(iv); err != nil {
			return err
		}
	}
	return nil
}

// mpileupReadBAQInfo caches the per-read CIGAR facts that mplp_realn's
// realignment heuristic needs, so they are computed once per read
// instead of once per covered column.
type mpileupReadBAQInfo struct {
	rec       *sam.Record
	beg       int  // 0-based reference start
	end       int  // 0-based reference end (exclusive)
	hasIndel  bool // CIGAR contains an I/D/N op (upstream PLP_HAS_INDEL)
	hasClip   bool // CIGAR contains a soft-clip op (PLP_HAS_SOFT_CLIP)
	ncig      int  // number of CIGAR ops
	leadMatch int  // leading consecutive M/=/X reference length (lm)
	tailMatch int  // trailing consecutive M/=/X reference length (rm)
	allMatch  bool // every CIGAR op is M/=/X (nm == ncig)
	realigned bool // BAQ already applied (upstream PLP_IS_REALN)
}

// mpileupBuildBAQInfo derives the mplp_realn heuristic facts for rec. lr
// (long-read) controls whether clip ops are skipped while measuring the
// leading/trailing match runs, mirroring mplp_realn's `lr` branch.
func mpileupBuildBAQInfo(rec *sam.Record) mpileupReadBAQInfo {
	info := mpileupReadBAQInfo{rec: rec, beg: int(rec.Pos) - 1, ncig: len(rec.Cigar)}
	refLen := 0
	for _, op := range rec.Cigar {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			refLen += int(op.Length())
		case sam.CigarDeletion, sam.CigarSkipped:
			refLen += int(op.Length())
			info.hasIndel = true
		case sam.CigarInsertion:
			info.hasIndel = true
		case sam.CigarSoftClip:
			info.hasClip = true
		}
	}
	info.end = info.beg + refLen
	lr := len(rec.Seq) > 500
	// Leading match run.
	nm := 0
	for _, op := range rec.Cigar {
		o := op.Op()
		if lr && (o == sam.CigarHardClip || o == sam.CigarSoftClip) {
			continue
		}
		if o == sam.CigarMatch || o == sam.CigarEqual || o == sam.CigarMismatch {
			info.leadMatch += int(op.Length())
			nm++
		} else {
			break
		}
	}
	info.allMatch = nm == info.ncig
	// Trailing match run.
	for k := len(rec.Cigar) - 1; k >= 0; k-- {
		o := rec.Cigar[k].Op()
		if lr && (o == sam.CigarHardClip || o == sam.CigarSoftClip) {
			continue
		}
		if o == sam.CigarMatch || o == sam.CigarEqual || o == sam.CigarMismatch {
			info.tailMatch += int(rec.Cigar[k].Length())
		} else {
			break
		}
	}
	return info
}

// applyMpileupBAQ ports mpileup.c's mplp_realn: it walks every covered
// reference column and runs baq.SamProbRealn (apply+extend mode) on each
// read the first time a column it covers selects it. A read is realigned
// at most once, matching upstream's PLP_IS_REALN dedup.
//
// By default upstream sets MPLP_REALN_PARTIAL (mpileup.c:1389): the
// per-column has_indel/soft-clip skip heuristic and the per-read
// spanning check both apply. `-D/--full-BAQ` (opts.FullBAQ) clears
// MPLP_REALN_PARTIAL (mpileup.c:1567), so both of those checks are
// bypassed and every read on the chromosome is realigned ("full BAQ").
//
// One faithful-port caveat: upstream's per-column `p->indel` term (an
// indel event adjacent to the column, supplied by the pileup engine) is
// not available without indel detection (slice 4). has_indel here counts
// only reads whose CIGAR carries an I/D/N op (PLP_HAS_INDEL), so for
// indel-bearing inputs the partial heuristic is a slight underestimate;
// for indel-free inputs it is exact.
func applyMpileupBAQ(perInputChromRecs [][]*sam.Record, refSlab []byte, opts MpileupOptions, eligible func(*sam.Record, int) bool, realigned map[*sam.Record]bool) {
	baqFlag := mpileupBAQFlag(opts)
	// max_read_len: upstream default is 500 unless -M overrides it.
	maxReadLen := opts.MaxReadLen
	if maxReadLen <= 0 {
		maxReadLen = 500
	}

	// Build per-read heuristic info and an interval index keyed by
	// covered reference position. The column-heuristic counters
	// (nt/has_indel/has_clip) include ALL reads at the column,
	// matching upstream mplp_realn (mpileup.c:430-441) which sees
	// the whole pileup pile regardless of which read it might end up
	// realigning at this column. eligible (when non-nil) only gates
	// which reads we actually realign — never which reads we count.
	// realigned (when non-nil) is a shared per-record dedup set that
	// persists across multiple BAQ phases — phase 1 marks reads it
	// realigns so phase 2 does not double-BAQ them; upstream's
	// PLP_IS_REALN flag plays the same role.
	var infos []*mpileupReadBAQInfo
	maxPos := 0
	for _, recs := range perInputChromRecs {
		for _, rec := range recs {
			if rec.IsUnmapped() || len(rec.Cigar) == 0 {
				continue
			}
			info := mpileupBuildBAQInfo(rec)
			if realigned != nil && realigned[rec] {
				info.realigned = true
			}
			infos = append(infos, &info)
			if info.end > maxPos {
				maxPos = info.end
			}
		}
	}
	if len(infos) == 0 {
		return
	}
	// column[pos0] lists the reads overlapping that column.
	column := make([][]*mpileupReadBAQInfo, maxPos)
	for _, info := range infos {
		for p := info.beg; p < info.end; p++ {
			if p >= 0 && p < maxPos {
				column[p] = append(column[p], info)
			}
		}
	}

	// partial mirrors upstream's MPLP_REALN_PARTIAL bit: set by default,
	// cleared by -D/--full-BAQ (opts.FullBAQ). When false the per-column
	// skip heuristic and the per-read spanning check are both bypassed.
	partial := !opts.FullBAQ

	for pos0 := 0; pos0 < maxPos; pos0++ {
		col := column[pos0]
		if len(col) == 0 {
			continue
		}
		nt := len(col)
		hasIndel, hasClip := 0, 0
		for _, info := range col {
			if info.hasIndel {
				hasIndel++
			}
			if info.hasClip {
				hasClip++
			}
		}
		// MPLP_REALN_PARTIAL skip heuristic (mpileup.c:445). max_indel
		// and min_indel both collapse to 0 here (no per-column indel
		// term), so max_indel==min_indel is always satisfied. Skipped
		// entirely under -D/--full-BAQ.
		if partial {
			if hasIndel == 0 ||
				(float64(hasClip) < 0.2*float64(nt) &&
					(float64(hasIndel) < 0.1*float64(nt) || hasIndel == 1)) {
				continue
			}
		}
		// realnDist mirrors the REALN_DIST macro.
		realnDist := 40
		if nt < 40 {
			realnDist += 10
		}
		if nt < 20 {
			realnDist += 10
		}
		for _, info := range col {
			if info.realigned {
				continue
			}
			if eligible != nil && !eligible(info.rec, pos0) {
				// Skip but do NOT mark realigned: this read may be
				// processed in a later phase (second-mate BAQ after
				// overlap-merge). The shared realigned map persists
				// the "already BAQ'd" status across phases so a
				// first-mate that was realigned in phase 1 is not
				// re-realigned in phase 2 — matching upstream's
				// per-read PLP_IS_REALN dedup.
				continue
			}
			info.realigned = true
			if realigned != nil {
				realigned[info.rec] = true
			}
			if len(info.rec.Seq) > maxReadLen {
				continue
			}
			// Per-read spanning check (mpileup.c:495). Only when
			// MPLP_REALN_PARTIAL is on, nt > 15 and the read has more
			// than one CIGAR op. Bypassed entirely under -D/--full-BAQ.
			if partial && nt > 15 && info.ncig > 1 && !info.allMatch {
				lm, rm := info.leadMatch, info.tailMatch
				if lm >= realnDist*4 && rm >= realnDist*4 {
					continue
				}
				clipThresh := 0.15
				if nt > 20 {
					clipThresh = 0.20
				}
				if lm >= realnDist && rm >= realnDist &&
					float64(hasClip) < clipThresh*float64(nt) {
					continue
				}
			}
			// Long-read band-width blow-up guard (mpileup.c:540-545):
			// for reads longer than 500bp, skip BAQ when the gap
			// between the read's reference span and its query length
			// would force an expensive wide alignment band. rl is the
			// CIGAR reference length (bam_cigar2rlen) — info.end-info.beg.
			if qseq := len(info.rec.Seq); qseq > 500 {
				rl := info.end - info.beg
				diff := rl - qseq
				if diff < 0 {
					diff = -diff
				}
				if diff*qseq >= 500000 {
					continue
				}
			}
			baq.SamProbRealn(info.rec, refSlab, baqFlag)
		}
	}
}

// accumulateMpileupBases walks rec's CIGAR and appends one pileupBase
// per covered reference position into events[pos0]. The base quality is
// captured raw together with its neighbours so glfgen can apply the
// delta_baseQ cap. CIGAR ops that produce no SNP-candidate base
// (D, N, S, H, P, I) are skipped as far as the SNP pileup is
// concerned, but the column immediately before an I/D op has its
// pileupBase.indel set to match upstream htslib's bam_pileup1_t.indel
// semantics (positive for an insertion of that length, negative for a
// deletion). This is the input the indel calling slice (4c+4d) consumes.
//
// preMergeQual and preMergeDrain (when non-nil) carry the
// pre-overlap-merge quality snapshot for first-mates of overlapping
// pairs together with a per-record drain threshold. For columns C
// strictly less than that threshold, the delta_baseQ neighbour
// qualities (prevQ/nextQ) are read from the snapshot rather than from
// rec.Qual — matching upstream's bcf_call_glfgen at bam2bcf.c:428-435:
// at iter->pos=C, qual[qpos±1] reflects only the overlap-merges that
// fired in bam_plp_push (sam.c:5970-5980) before C was drained, which
// happens when no intermediate read between the two mates lifted
// iter->max_pos past C. snapshotFirstMateQuals derives the threshold as
// max(X.Pos) over intermediate reads X (initial value F.Pos so the
// predicate is unreachable when no intermediate raises max_pos). Reads
// without an entry in preMergeQual (second mates, standalones, or all
// reads when -x is in force) always read rec.Qual.
func accumulateMpileupBases(rec *sam.Record, events [][]pileupBase, preMergeQual map[*sam.Record][]byte, preMergeDrain map[*sam.Record]int) {
	if rec.Pos <= 0 {
		return
	}
	refPos := int(rec.Pos) - 1
	queryPos := 0
	isReverse := rec.Flag&sam.FlagReverse != 0
	qlen := len(rec.Seq)
	// preMerge / drainAt are the pre-overlap-merge neighbour-qual
	// snapshot for this record (nil if rec is not a first-mate of an
	// overlapping pair) and the column threshold below which upstream's
	// pileup engine had not yet fired overlap_push when it drained the
	// column. For columns C < drainAt the neighbour qualities are read
	// from preMerge; from C >= drainAt onwards rec.Qual already carries
	// the post-merge values (or was never merged) and is used directly.
	var preMerge []byte
	drainAt := -1
	if preMergeQual != nil {
		if q, ok := preMergeQual[rec]; ok {
			preMerge = q
			drainAt = preMergeDrain[rec]
		}
	}
	// CIGAR op codes / lengths in BAM order, for get_position's
	// soft-clip-aware read-position and soft-clip-length annotations.
	cigarOps := make([]int, len(rec.Cigar))
	cigarLens := make([]int, len(rec.Cigar))
	hasSC := false
	for k, op := range rec.Cigar {
		cigarOps[k] = int(op.Op())
		cigarLens[k] = int(op.Length())
		if op.Op() == sam.CigarSoftClip {
			hasSC = true
		}
	}
	// Track whether this read ever produces an indel-bearing column so
	// the back-pointer to rec can be set on those columns. The
	// per-column indels themselves are computed by peeking at the next
	// consuming CIGAR op when walking the read.
	for idx, op := range rec.Cigar {
		l := int(op.Length())
		o := op.Op()
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			// nextIndel is the indel value to stamp on the LAST base
			// of this match run (htslib pileup engine semantics:
			// p->indel is read from the CIGAR op that follows the
			// match op covering this column). Upstream reference:
			// htslib sam.c:5463-5491 — same-type indel ops are
			// MERGED (1D2D -> -3, 2I1I -> +3), and a CPAD-leading
			// insertion run (e.g. M P I) accumulates inserts across
			// CPAD ops. We reproduce that here so the per-column
			// p->indel matches bam_plp_next.
			nextIndel := 0
			if idx+1 < len(rec.Cigar) {
				no2 := rec.Cigar[idx+1].Op()
				nl2 := int(rec.Cigar[idx+1].Length())
				switch {
				case no2 == sam.CigarDeletion:
					// Start of a new deletion: merge e.g. 1D2D to
					// 3D by accumulating any further D ops.
					nextIndel = -nl2
					for k := idx + 2; k < len(rec.Cigar); k++ {
						if rec.Cigar[k].Op() == sam.CigarDeletion {
							nextIndel -= int(rec.Cigar[k].Length())
						} else {
							break
						}
					}
				case no2 == sam.CigarInsertion:
					// Insertion run: accumulate consecutive I ops,
					// skipping CPAD ops between them, stop on
					// anything else.
					nextIndel = nl2
					for k := idx + 2; k < len(rec.Cigar); k++ {
						kop := rec.Cigar[k].Op()
						if kop == sam.CigarInsertion {
							nextIndel += int(rec.Cigar[k].Length())
						} else if kop == sam.CigarPadding {
							continue
						} else {
							break
						}
					}
				case no2 == sam.CigarPadding && idx+2 < len(rec.Cigar):
					// Pure CPAD-leading run; sum any I ops until a
					// consuming op terminates the run (D/M/N/=/X).
					l3 := 0
					for k := idx + 2; k < len(rec.Cigar); k++ {
						kop := rec.Cigar[k].Op()
						if kop == sam.CigarInsertion {
							l3 += int(rec.Cigar[k].Length())
						} else if kop == sam.CigarDeletion ||
							kop == sam.CigarMatch ||
							kop == sam.CigarSkipped ||
							kop == sam.CigarEqual ||
							kop == sam.CigarMismatch {
							break
						}
					}
					if l3 > 0 {
						nextIndel = l3
					}
				}
				// Otherwise (N, M/=/X, S, H without P-leading): no
				// indel; nextIndel stays 0, mirroring upstream
				// where p->indel = 0 is the default.
			}
			for k := 0; k < l; k++ {
				p := refPos + k
				q := queryPos + k
				if p < 0 || p >= len(events) {
					continue
				}
				if q >= len(rec.Seq) {
					continue
				}
				base := upperByte(rec.Seq[q])
				b4 := seqNt16Int[baseToNt16(base)]
				var rawQual uint8
				if q < len(rec.Qual) {
					rawQual = rec.Qual[q]
				}
				// neighbourSrc selects the qual array upstream's
				// bcf_call_glfgen would see at iter->pos = p. For a
				// first-mate of an overlapping pair, columns drained
				// before overlap_push fired (p < drainAt) predate the
				// merge and must read the pre-merge snapshot; columns
				// from drainAt onwards see the post-merge rec.Qual
				// (because by the time they were drained, mate2's push
				// had already fired tweak_overlap_quality).
				neighbourSrc := rec.Qual
				if preMerge != nil && p < drainAt {
					neighbourSrc = preMerge
				}
				prevQ := -1
				if q > 0 && q-1 < len(neighbourSrc) {
					prevQ = int(neighbourSrc[q-1])
				}
				nextQ := -1
				if q+1 < qlen && q+1 < len(neighbourSrc) {
					nextQ = int(neighbourSrc[q+1])
				}
				gp := getPosition(cigarOps, cigarLens, q, qlen)
				epos, scLen := biasPositionBins(gp)
				ind := 0
				if k == l-1 {
					ind = nextIndel
				}
				// Always carry the back-pointer to the originating read.
				// bcfCallGapPrep needs it to compute the indel-flavored
				// get_pos result for every read at the column (not just
				// the indel-bearing ones), so the iref_*/ialt_*
				// histograms can be populated.
				recPtr := rec
				events[p] = append(events[p], pileupBase{
					base4:       b4,
					rawQual:     rawQual,
					prevQ:       prevQ,
					nextQ:       nextQ,
					mapq:        rec.MapQ,
					reverse:     isReverse,
					qpos:        q,
					qlen:        qlen,
					qname:       rec.QName,
					epos:        epos,
					scLen:       scLen,
					indel:       ind,
					rec:         recPtr,
					hasSoftClip: hasSC,
				})
			}
			refPos += l
			queryPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			// Spanning deletion / ref-skip: emit one pileupBase per
			// covered reference column so the indel branch's glfgen
			// can iterate these reads (upstream bam2bcf.c:307 lets
			// is_del reads through when is_indel=1; bam2bcf.c:301 has
			// is_refskip skip both branches). The SNP branch ignores
			// these columns via the isDel / isRefskip gate in
			// bcfCallGlfgenCore. We carry the back-pointer to the
			// originating record together with the read's mapq and
			// strand so the indel histograms remain accurate.
			isDel := o == sam.CigarDeletion
			isRefskip := o == sam.CigarSkipped
			// p->qpos in upstream points at the LAST consumed query
			// base before the D/N op (i.e. queryPos-1 here when at
			// least one match preceded). Use that for the
			// bias-position helper so iref/ialt histograms are
			// populated with sensible bins.
			qref := queryPos - 1
			if qref < 0 {
				qref = 0
			}
			var rawQual uint8
			if qref < len(rec.Qual) {
				rawQual = rec.Qual[qref]
			}
			gp := getPosition(cigarOps, cigarLens, qref, qlen)
			epos, scLen := biasPositionBins(gp)
			for k := 0; k < l; k++ {
				p := refPos + k
				if p < 0 || p >= len(events) {
					continue
				}
				events[p] = append(events[p], pileupBase{
					base4:       0,
					rawQual:     rawQual,
					prevQ:       -1,
					nextQ:       -1,
					mapq:        rec.MapQ,
					reverse:     isReverse,
					qpos:        qref,
					qlen:        qlen,
					qname:       rec.QName,
					epos:        epos,
					scLen:       scLen,
					indel:       0,
					rec:         rec,
					isDel:       isDel,
					isRefskip:   isRefskip,
					hasSoftClip: hasSC,
				})
			}
			refPos += l
		case sam.CigarInsertion, sam.CigarSoftClip:
			queryPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// no advance.
		}
	}
}

// wangHash is the Go port of htslib's __ac_Wang_hash (khash.h:528),
// operating on 32-bit unsigned integers like the C khint_t.
func wangHash(key uint32) uint32 {
	key += ^(key << 15)
	key ^= key >> 10
	key += key << 3
	key ^= key >> 6
	key += ^(key << 11)
	key ^= key >> 16
	return key
}

// x31HashString is the Go port of htslib's __ac_X31_hash_string
// (khash.h:454), the read-name hash used to pick which mate keeps its
// quality during overlap merging.
func x31HashString(s string) uint32 {
	if len(s) == 0 {
		return 0
	}
	h := uint32(s[0])
	for i := 1; i < len(s); i++ {
		h = (h << 5) - h + uint32(s[i])
	}
	return h
}

// cigarIter is the Go port of htslib's cigar_iref2iseq_set/next state
// machine (sam.c:5731). It walks a read's CIGAR yielding the sequence
// index of each M/=/X base together with its reference offset relative
// to the read start.
type cigarIter struct {
	ops  []int // BAM op codes
	lens []int // op lengths
	idx  int   // current CIGAR op
	icig int   // position within the current op
	iseq int   // sequence index
	iref int   // reference offset from the read start
}

// cigarSet positions the iterator at reference offset pos (iref). It
// returns true if pos is covered by an M/=/X op, false otherwise.
func (it *cigarIter) cigarSet(pos int) bool {
	if pos < 0 {
		return false
	}
	it.idx, it.icig, it.iseq, it.iref = 0, 0, 0, 0
	for it.idx < len(it.ops) {
		cig := it.ops[it.idx]
		ncig := it.lens[it.idx]
		switch cig {
		case sam.CigarSoftClip:
			it.idx++
			it.iseq += ncig
			it.icig = 0
		case sam.CigarHardClip, sam.CigarPadding:
			it.idx++
			it.icig = 0
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			pos -= ncig
			if pos < 0 {
				it.icig = ncig + pos
				it.iseq += it.icig
				it.iref += it.icig
				return true
			}
			it.idx++
			it.iseq += ncig
			it.icig = 0
			it.iref += ncig
		case sam.CigarInsertion:
			it.idx++
			it.iseq += ncig
			it.icig = 0
		case sam.CigarDeletion, sam.CigarSkipped:
			pos -= ncig
			if pos < 0 {
				pos = 0
			}
			it.idx++
			it.icig = 0
			it.iref += ncig
		default:
			return false
		}
	}
	it.iseq = -1
	return false
}

// cigarNext advances to the next M/=/X base. It returns true while a
// base is available and false when the CIGAR is exhausted.
func (it *cigarIter) cigarNext() bool {
	for it.idx < len(it.ops) {
		cig := it.ops[it.idx]
		ncig := it.lens[it.idx]
		switch cig {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if it.icig >= ncig-1 {
				it.icig = -1
				it.idx++
				continue
			}
			it.iseq++
			it.icig++
			it.iref++
			return true
		case sam.CigarDeletion, sam.CigarSkipped:
			it.idx++
			it.iref += ncig
			it.icig = -1
		case sam.CigarInsertion, sam.CigarSoftClip:
			it.idx++
			it.iseq += ncig
			it.icig = -1
		case sam.CigarHardClip, sam.CigarPadding:
			it.idx++
			it.icig = -1
		default:
			return false
		}
	}
	it.iseq, it.iref = -1, -1
	return false
}

// newCigarIter builds a cigarIter from a record's CIGAR.
func newCigarIter(c sam.Cigar) *cigarIter {
	it := &cigarIter{ops: make([]int, len(c)), lens: make([]int, len(c))}
	for i, op := range c {
		it.ops[i] = int(op.Op())
		it.lens[i] = int(op.Length())
	}
	return it
}

// seqBase returns the uppercase nucleotide at sequence index i of rec,
// or 0 if out of range.
func seqBase(rec *sam.Record, i int) byte {
	if i < 0 || i >= len(rec.Seq) {
		return 0
	}
	return upperByte(rec.Seq[i])
}

// tweakOverlapQuality is the Go port of htslib's tweak_overlap_quality
// (sam.c:5805). For two reads of the same template whose alignments
// overlap, it merges the base qualities in the overlapping reference
// span: matching bases get the summed quality on the kept mate and 0 on
// the other; mismatching bases keep 0.8x the higher quality and zero
// the lower. Which mate keeps quality is chosen by a name hash so the
// choice is deterministic and unbiased.
func tweakOverlapQuality(a, b *sam.Record) {
	aCig := newCigarIter(a.Cigar)
	bCig := newCigarIter(b.Cigar)
	aPos := int(a.Pos) - 1
	bPos := int(b.Pos) - 1
	iref := bPos
	if !aCig.cigarSet(iref - aPos) {
		return
	}
	if !bCig.cigarSet(iref - bPos) {
		return
	}

	// Pick which mate keeps quality, by read-name hash.
	var amul, bmul uint8
	if wangHash(x31HashString(a.QName))&1 != 0 {
		amul, bmul = 1, 0
	} else {
		amul, bmul = 0, 1
	}

	for {
		// Step a and b to the next matching reference position.
		for aCig.iref >= 0 && aCig.iref < iref-aPos {
			if !aCig.cigarNext() {
				return
			}
		}
		if aCig.iref < 0 {
			return
		}
		for bCig.iref >= 0 && bCig.iref < iref-bPos {
			if !bCig.cigarNext() {
				return
			}
		}
		if bCig.iref < 0 {
			return
		}

		if iref < aCig.iref+aPos {
			iref = aCig.iref + aPos
		}
		if iref < bCig.iref+bPos {
			iref = bCig.iref + bPos
		}
		iref++

		// If a or b has a deletion the other catches up. Upstream zeroes
		// (or scales by 0.8 on the kept mate) the caught-up bases, the
		// same rule it uses for mismatches.
		if aCig.iref+aPos != bCig.iref+bPos {
			switch {
			case aCig.iref+aPos < bCig.iref+bPos && bCig.idx > 0 &&
				bCig.ops[bCig.idx-1] == sam.CigarDeletion:
				for aCig.iref+aPos < bCig.iref+bPos {
					if amul != 0 {
						a.Qual[aCig.iseq] = uint8(float64(a.Qual[aCig.iseq]) * 0.8)
					} else {
						a.Qual[aCig.iseq] = 0
					}
					if !aCig.cigarNext() {
						return
					}
				}
				continue
			case aCig.idx > 0 && aCig.ops[aCig.idx-1] == sam.CigarDeletion:
				for bCig.iref+bPos < aCig.iref+aPos {
					if bmul != 0 {
						b.Qual[bCig.iseq] = uint8(float64(b.Qual[bCig.iseq]) * 0.8)
					} else {
						b.Qual[bCig.iseq] = 0
					}
					if !bCig.cigarNext() {
						return
					}
				}
				continue
			default:
				continue
			}
		}

		if aCig.iseq >= len(a.Qual) || bCig.iseq >= len(b.Qual) {
			return
		}

		aq := a.Qual[aCig.iseq]
		bq := b.Qual[bCig.iseq]
		if seqBase(a, aCig.iseq) == seqBase(b, bCig.iseq) {
			// Confident: keep the summed quality, capped at 200.
			qual := int(aq) + int(bq)
			if qual > 200 {
				qual = 200
			}
			a.Qual[aCig.iseq] = uint8(amul) * uint8(qual)
			b.Qual[bCig.iseq] = uint8(bmul) * uint8(qual)
		} else {
			// Mismatch: keep 0.8x the higher quality, zero the lower.
			switch {
			case aq > bq:
				a.Qual[aCig.iseq] = uint8(0.8 * float64(aq))
				b.Qual[bCig.iseq] = 0
			case aq < bq:
				b.Qual[bCig.iseq] = uint8(0.8 * float64(bq))
				a.Qual[aCig.iseq] = 0
			default:
				a.Qual[aCig.iseq] = uint8(float64(amul) * 0.8 * float64(aq))
				b.Qual[bCig.iseq] = uint8(float64(bmul) * 0.8 * float64(bq))
			}
		}
	}
}

// mateClass labels a read's role within an overlap-pair, as identified
// by the same predicate htslib's overlap_push uses (sam.c:5950): the
// first mate to arrive (coordinate-sorted) is mateClassFirst, the
// second is mateClassSecond. Reads that never get paired (no mate in
// the input, mate fails the proper-pair / coord guards, or wild CIGAR)
// are absent from the classifyMatePairs map and treated as standalones.
type mateClass uint8

const (
	mateClassFirst mateClass = iota + 1
	mateClassSecond
)

// matePairInfo records a read's pair role plus the 0-based reference
// start of the read's mate. mateStart lets the BAQ engine answer
// "would overlap_push have run by the time BAQ realigns this read?":
// upstream htslib's overlap_push runs when the second mate is pushed
// (at its alignment start). For a first-mate, if its first eligible
// BAQ column precedes its mate's start, BAQ sees raw quals; otherwise
// BAQ sees the post-merge quals. For a second-mate, BAQ always runs
// after overlap-merge because its own push triggers overlap_push at
// or before its first eligible column.
type matePairInfo struct {
	class     mateClass
	mateStart int
}

// classifyMatePairs walks each input's records in arrival (coordinate)
// order and labels every read that will be paired up by
// applySmartOverlaps. The predicate must stay byte-identical to
// applySmartOverlaps' loop so the two passes agree on which reads are
// "first" and "second" mates. The result is keyed by *sam.Record
// pointer identity, which is stable across BAQ / overlap-merge phases
// because we never copy records.
//
// This pre-classification lets emitChromMpileup interleave BAQ and
// overlap-merge per upstream's bam_plp_push ordering: first-mates
// whose BAQ trigger column precedes their mate's start are BAQ'd on
// raw quals; the remaining first-mates plus all second-mates are
// BAQ'd after overlap-merge.
func classifyMatePairs(perInputChromRecs [][]*sam.Record) map[*sam.Record]matePairInfo {
	out := make(map[*sam.Record]matePairInfo)
	for _, recs := range perInputChromRecs {
		buffered := make(map[string]*sam.Record)
		for _, rec := range recs {
			if rec.IsUnmapped() || rec.IsMateUnmapped() || !rec.IsProperPair() {
				continue
			}
			pos := int(rec.Pos) - 1
			mpos := int(rec.PNext) - 1
			end := int(rec.EndPosition())
			if rec.RNext != "" && rec.RNext != "=" && rec.RNext != rec.RName {
				continue
			}
			isize := int(rec.TLen)
			if isize < 0 {
				isize = -isize
			}
			if isize >= 2*len(rec.Seq) && mpos >= end {
				continue
			}
			if mate, ok := buffered[rec.QName]; ok {
				out[mate] = matePairInfo{class: mateClassFirst, mateStart: pos}
				out[rec] = matePairInfo{class: mateClassSecond, mateStart: int(mate.Pos) - 1}
				delete(buffered, rec.QName)
				continue
			}
			if mpos >= pos || (rec.IsPaired() && rec.PNext == 0) {
				buffered[rec.QName] = rec
			}
		}
	}
	return out
}

// snapshotFirstMateQuals returns a per-record pre-merge copy of rec.Qual
// for each first-mate of an overlapping pair (as labelled by pairs), plus
// a parallel map from each first-mate to the column threshold below
// which its pile entries must consult the snapshot for delta_baseQ
// neighbour quals. The threshold mirrors upstream htslib's pileup-
// iterator timing: in bam_plp_push (sam.c:6083-6132) the overlap_push
// call that merges the pair's quals fires only when the SECOND mate is
// pushed (sam.c:5970-5980), and the engine emits a column C as soon as
// iter->max_pos > C — which is set to the maximum BAM-order position
// pushed so far. So overlap_push fires "in time" for column C only when
// no read pushed strictly between the two mates has Pos > C; otherwise
// that earlier intermediate read already drove iter->max_pos > C and
// drained C with raw quals BEFORE mate2 (and overlap_push) arrived.
//
// drainThreshold for a first-mate F is therefore max(X.Pos) over reads
// X strictly between F and its mate M in BAM order (initial value F.Pos
// so the predicate "C < threshold" stays unreachable when there are no
// intermediates with Pos > F.Pos). Columns C < drainThreshold use the
// snapshot; columns C >= drainThreshold use the post-merge rec.Qual.
//
// Second mates need no snapshot: by the time upstream's iter->pos
// reaches any column they cover, overlap_push has already fired (it
// fired when they were pushed, before iter->max_pos could advance past
// their own start), so accumulateMpileupBases can read rec.Qual
// directly for them.
// computeFirstMateDrainThresholds returns just the per-first-mate drain
// thresholds (see snapshotFirstMateQuals for the derivation) without
// allocating the quality snapshot. It is called before phase-1 BAQ so
// the phase-1 / phase-2 eligibility predicates can gate on the same
// upstream-equivalent timing the snapshot lookup uses; the snapshot
// itself is still produced after phase 1 by snapshotFirstMateQuals
// (the threshold is push-order-only and unaffected by BAQ).
func computeFirstMateDrainThresholds(perInputChromRecs [][]*sam.Record, pairs map[*sam.Record]matePairInfo) map[*sam.Record]int {
	drainThreshold := make(map[*sam.Record]int)
	for _, recs := range perInputChromRecs {
		pending := make(map[*sam.Record]int)
		matePending := make(map[string]*sam.Record)
		for _, rec := range recs {
			pos := int(rec.Pos) - 1
			p, pairOk := pairs[rec]
			if pairOk && p.class == mateClassSecond {
				if f, ok2 := matePending[rec.QName]; ok2 {
					drainThreshold[f] = pending[f]
					delete(pending, f)
					delete(matePending, rec.QName)
				}
			}
			for f, m := range pending {
				if pos > m {
					pending[f] = pos
				}
			}
			if pairOk && p.class == mateClassFirst {
				pending[rec] = pos
				matePending[rec.QName] = rec
			}
		}
		for f, m := range pending {
			drainThreshold[f] = m
		}
	}
	return drainThreshold
}

func snapshotFirstMateQuals(perInputChromRecs [][]*sam.Record, pairs map[*sam.Record]matePairInfo) (map[*sam.Record][]byte, map[*sam.Record]int) {
	snap := make(map[*sam.Record][]byte)
	drainThreshold := make(map[*sam.Record]int)
	for _, recs := range perInputChromRecs {
		// Walk this input's records in coordinate (= push) order so we
		// can compute, for each first-mate F, the running max of any
		// later reads' positions until F's mate is encountered. pending
		// maps first-mate -> current max(X.Pos) for X seen strictly
		// between F and (still-to-arrive) M. The initial max is F.Pos
		// so reads with X.Pos > F.Pos can lift it; reads with
		// X.Pos == F.Pos (same-position bystanders pushed before F's
		// mate) don't change drainage.
		pending := make(map[*sam.Record]int)
		matePending := make(map[string]*sam.Record)
		for _, rec := range recs {
			pos := int(rec.Pos) - 1
			p, pairOk := pairs[rec]
			// Resolve a pending pair FIRST when this rec is its second
			// mate: mate2's own push is what fires overlap_push, so it
			// must NOT be counted as an intermediate that lifts the
			// drain threshold. Once the threshold is frozen we proceed
			// with the lift loop for whichever first-mates are still
			// pending (other pairs whose mate has yet to arrive).
			if pairOk && p.class == mateClassSecond {
				if f, ok2 := matePending[rec.QName]; ok2 {
					drainThreshold[f] = pending[f]
					delete(pending, f)
					delete(matePending, rec.QName)
				}
			}
			// Lift drain thresholds of every still-pending first-mate
			// whose mate has not yet been encountered. This rec's push
			// would have set upstream's iter->max_pos to max(prev, pos),
			// draining cols < pos for those pending first-mates if pos
			// exceeds the previous running max.
			for f, m := range pending {
				if pos > m {
					pending[f] = pos
				}
			}
			if pairOk && p.class == mateClassFirst {
				// Snapshot F's quality array NOW (called just before
				// applySmartOverlaps, after any pre-merge BAQ; this is
				// the post-mplp_realn / pre-overlap_push state upstream
				// keeps in bam1_t->qual until its mate's push fires
				// tweak_overlap_quality at sam.c:5982). The drain
				// threshold starts at F.Pos and grows as intermediates
				// are pushed; it freezes when F's mate appears.
				buf := make([]byte, len(rec.Qual))
				copy(buf, rec.Qual)
				snap[rec] = buf
				pending[rec] = pos
				matePending[rec.QName] = rec
			}
		}
		// First-mates whose mate never showed up (shouldn't happen
		// when classifyMatePairs labelled them, but be defensive)
		// still get a sensible threshold = current running max so the
		// snapshot is never consulted past that point.
		for f, m := range pending {
			drainThreshold[f] = m
		}
	}
	return snap, drainThreshold
}

// applySmartOverlaps ports htslib's overlap_push (sam.c:5950): for each
// input it pairs up the two mates of every proper read pair and calls
// tweakOverlapQuality so a base covered by both mates is not counted
// twice. Reads must be in coordinate order (as the pileup engine sees
// them); the per-input record slices already are.
func applySmartOverlaps(perInputChromRecs [][]*sam.Record) {
	for _, recs := range perInputChromRecs {
		buffered := make(map[string]*sam.Record)
		for _, rec := range recs {
			// Mapped reads in proper pairs only.
			if rec.IsUnmapped() || rec.IsMateUnmapped() || !rec.IsProperPair() {
				continue
			}
			pos := int(rec.Pos) - 1
			mpos := int(rec.PNext) - 1
			end := int(rec.EndPosition()) // 0-based exclusive ref end
			// No overlap possible for mates on a different contig or for
			// wild CIGARs (matching the overlap_push guard).
			if rec.RNext != "" && rec.RNext != "=" && rec.RNext != rec.RName {
				continue
			}
			isize := int(rec.TLen)
			if isize < 0 {
				isize = -isize
			}
			if isize >= 2*len(rec.Seq) && mpos >= end {
				continue
			}
			if mate, ok := buffered[rec.QName]; ok {
				tweakOverlapQuality(mate, rec)
				delete(buffered, rec.QName)
				continue
			}
			// Buffer reads whose mate is still to arrive.
			if mpos >= pos || (rec.IsPaired() && rec.PNext == 0) {
				buffered[rec.QName] = rec
			}
		}
	}
}

// filterMpileupPile finalises a per-position pileup column as the
// per-sample glfgen input. -Q/--min-BQ is applied inside glfgen (after
// the delta_baseQ cap), matching upstream's ordering. Read-pair overlap
// handling is not done here: it is a pre-pileup quality merge
// (applySmartOverlaps), matching upstream where overlap_push runs while
// reads are pushed into the pileup engine. The -d/--max-depth cap is
// applied earlier, in applyMpileupDepthCap, matching htslib's per-
// alignment-start drop semantics in bam_plp_push (sam.c:6090).
func filterMpileupPile(evs []pileupBase) []pileupBase {
	if len(evs) == 0 {
		return nil
	}
	out := make([]pileupBase, len(evs))
	copy(out, evs)
	return out
}

// endHeap is a min-heap of 0-based exclusive reference end positions,
// used by applyMpileupDepthCap to track in-queue read lifetimes when
// porting htslib's per-alignment-start depth cap.
type endHeap []int

func (h endHeap) Len() int            { return len(h) }
func (h endHeap) Less(i, j int) bool  { return h[i] < h[j] }
func (h endHeap) Swap(i, j int)       { h[i], h[j] = h[j], h[i] }
func (h *endHeap) Push(x interface{}) { *h = append(*h, x.(int)) }
func (h *endHeap) Pop() interface{} {
	n := len(*h)
	x := (*h)[n-1]
	*h = (*h)[:n-1]
	return x
}

// applyMpileupDepthCap drops reads using htslib's per-alignment-start
// rule from bam_plp_push (reference_code/htslib/sam.c:6090): when a new
// read is pushed at b->core.pos == iter->pos and mp->cnt > maxcnt, the
// read is silently discarded. Because htslib's mempool keeps one
// sentinel node alive, the predicate is equivalent to "drop when the
// number of accepted, still-active reads already in the queue is at
// least maxcnt".
//
// The cap fires only at the boundary between same-start reads: between
// two pushes the pileup engine has drained columns up to the last
// accepted read's alignment start, so iter->pos at push time equals the
// previous accepted read's start position. The cap therefore differs
// from the previous per-column cap (which dropped reads once a column
// already had MaxDepth bases regardless of where they started) and
// matches upstream's golden output at deep homopolymer columns where
// many reads share an alignment start.
//
// Records must be coordinate-sorted on input (they are, post
// mpileupReadBAM). The slice in perInputChromRecs[i] is rewritten in
// place to retain only the accepted reads.
func applyMpileupDepthCap(perInputChromRecs [][]*sam.Record, maxcnt int) {
	if maxcnt <= 0 {
		return
	}
	for i, recs := range perInputChromRecs {
		if len(recs) == 0 {
			continue
		}
		ends := &endHeap{}
		// iter->pos starts at 0 (calloc-zeroed in bam_plp_init). The
		// cap predicate is iter->pos == b->core.pos; the very first
		// read therefore only triggers the cap when its 0-based start
		// is also 0, which our zero initialisation reproduces.
		iterPos := int32(0)
		// lastAccepted tracks whether any read on this contig has been
		// accepted yet. Before the first accepted read, iter->pos is
		// still its init value, so we honour the upstream check via
		// iterPos == startPos directly.
		out := recs[:0]
		for _, rec := range recs {
			if rec == nil || rec.Pos <= 0 {
				continue
			}
			startPos := rec.Pos - 1 // 0-based
			// Drain reads whose end (0-based exclusive) is <= iterPos.
			// htslib frees these inside bam_plp64_next before counting
			// mp->cnt, so they no longer occupy a queue slot.
			for ends.Len() > 0 && (*ends)[0] <= int(iterPos) {
				heap.Pop(ends)
			}
			// htslib's bam_plp_push cap check: drop the read when
			// iter->pos == b->core.pos and mp->cnt > maxcnt. mp->cnt
			// equals "active reads + 1" (one initial sentinel node),
			// so the >maxcnt test fires once the number of active
			// reads reaches maxcnt.
			if iterPos == startPos && ends.Len() >= maxcnt {
				continue
			}
			out = append(out, rec)
			// htslib only allocates a new tail node (incrementing
			// mp->cnt) when tail->end > iter->pos, i.e. the read is
			// still active. Pure-insertion CIGARs with refLen == 0
			// satisfy end == startPos and therefore do not queue.
			endPos := int(rec.EndPosition())
			if endPos > int(startPos) {
				heap.Push(ends, endPos)
			}
			// After this push the pileup engine drains columns up to
			// the read's alignment start: iter->pos catches up to the
			// new max_pos == startPos before the next bam_plp_push
			// call returns NULL.
			iterPos = startPos
		}
		perInputChromRecs[i] = out
	}
}

// bcfCall2bcf is the Go port of bcf_call2bcf (bam2bcf.c:1200) for the
// SNP path. It turns a combined bcfCall into a vcf.Variant: REF, the
// ALT alleles (including the `<*>` unseen allele), QUAL=0, INFO/DP/I16/
// QS/MQ0F and FORMAT/PL plus FORMAT/AD,ADF,ADR controlled by fmtFlag
// (the resolved `-a/--annotate` bitset).
func bcfCall2bcf(chrom string, pos1 int, refB byte, call *bcfCall, fmtFlag uint32) *vcf.Variant {
	alleles := make([]string, 0, call.nAlleles)
	alleles = append(alleles, string(refB)) // REF
	for i := 1; i < call.nAlleles; i++ {
		if call.unseen == i {
			alleles = append(alleles, "<*>")
		} else {
			alleles = append(alleles, string("ACGTN"[call.alleles[i]]))
		}
	}
	alt := alleles[1:]

	// INFO/DP, I16, QS, MQ0F. The order matches bam2bcf.c:1300-1336.
	info := map[string]string{}
	infoOrder := make([]string, 0, 4)
	info["DP"] = strconv.Itoa(call.oriDepth)
	infoOrder = append(infoOrder, "DP")

	// INFO/SCR is emitted before I16, mirroring bam2bcf.c:1298-1299
	// (after the per-sample AD/ADF/ADR INFO tags — which the current
	// port does not yet emit at the INFO level — and before I16/QS).
	if fmtFlag&B2BInfoSCR != 0 {
		info["SCR"] = strconv.Itoa(call.scrTotal)
		infoOrder = append(infoOrder, "SCR")
	}

	var i16 strings.Builder
	for j, v := range call.anno {
		if j > 0 {
			i16.WriteByte(',')
		}
		i16.WriteString(formatI16Number(v))
	}
	info["I16"] = i16.String()
	infoOrder = append(infoOrder, "I16")

	// INFO/QS carries one value per allele (the coverage-normalised
	// quality sum); the `<*>` allele's qsum is 0. Upstream stores QS as
	// a C float and renders it with %g, so round through float32 for
	// byte-for-byte parity.
	var qs strings.Builder
	for j := 0; j < call.nAlleles; j++ {
		if j > 0 {
			qs.WriteByte(',')
		}
		qs.WriteString(formatFloat32G(call.qsum[j]))
	}
	info["QS"] = qs.String()
	infoOrder = append(infoOrder, "QS")

	// Bias annotations, emitted only for sites with a real ALT allele
	// and only when the value is defined (upstream skips HUGE_VAL). The
	// order matches bam2bcf.c:1306-1336: VDB, SGB, RPBZ, MQBZ, MQSBZ,
	// BQBZ, NMBZ, SCBZ. NMBZ is gated by B2BInfoNMBZ; the others are in
	// the mpileup default set.
	if call.hasAlt {
		addBias := func(tag string, v float64, ok bool) {
			if ok {
				info[tag] = formatFloat32G(v)
				infoOrder = append(infoOrder, tag)
			}
		}
		addBias("VDB", call.vdb, call.vdbOK)
		addBias("SGB", call.segBias, call.sgbOK)
		addBias("RPBZ", call.mwuPos, call.rpbzOK)
		addBias("MQBZ", call.mwuMq, call.mqbzOK)
		addBias("MQSBZ", call.mwuMqs, call.mqsbzOK)
		addBias("BQBZ", call.mwuBq, call.bqbzOK)
		if fmtFlag&B2BInfoNMBZ != 0 {
			addBias("NMBZ", call.mwuNm, call.nmbzOK)
		}
		addBias("SCBZ", call.mwuSc, call.scbzOK)
	}

	mq0f := 0.0
	if call.oriDepth > 0 {
		mq0f = float64(call.mq0) / float64(call.oriDepth)
	}
	info["MQ0F"] = formatFloat32G(mq0f)
	infoOrder = append(infoOrder, "MQ0F")

	// FORMAT — PL first, then optional AD/ADF/ADR controlled by fmtFlag,
	// then SCR. Upstream emits ADF/ADR before AD (bam2bcf.c:1376-1384),
	// but does so for INFO; for FORMAT both are gated by fmtFlag bits
	// and the column order matches the bits' order: AD, ADF, ADR, SCR.
	format := []string{"PL"}
	emitAD := fmtFlag&B2BFmtAD != 0
	emitADF := fmtFlag&B2BFmtADF != 0
	emitADR := fmtFlag&B2BFmtADR != 0
	emitSCR := fmtFlag&B2BFmtSCR != 0
	if emitAD {
		format = append(format, "AD")
	}
	if emitADF {
		format = append(format, "ADF")
	}
	if emitADR {
		format = append(format, "ADR")
	}
	if emitSCR {
		format = append(format, "SCR")
	}
	samplesOut := make([]vcf.Sample, len(call.pl))
	for s := range call.pl {
		var pl strings.Builder
		for k, v := range call.pl[s] {
			if k > 0 {
				pl.WriteByte(',')
			}
			pl.WriteString(strconv.Itoa(v))
		}
		data := map[string]string{"PL": pl.String()}
		if emitAD && s < len(call.adf) {
			data["AD"] = formatPerAlleleSum(call.adf[s], call.adr[s])
		}
		if emitADF && s < len(call.adf) {
			data["ADF"] = formatPerAllele(call.adf[s])
		}
		if emitADR && s < len(call.adr) {
			data["ADR"] = formatPerAllele(call.adr[s])
		}
		if emitSCR && s < len(call.scr) {
			data["SCR"] = strconv.Itoa(call.scr[s])
		}
		samplesOut[s] = vcf.Sample{Data: data}
	}

	return &vcf.Variant{
		Chrom:     chrom,
		Pos:       pos1,
		ID:        ".",
		Ref:       string(refB),
		Alt:       alt,
		Qual:      0,
		Filter:    []string{"."},
		Info:      info,
		InfoOrder: infoOrder,
		Format:    format,
		Samples:   samplesOut,
	}
}

// bcfCall2bcfIndel is the Go port of the `ori_ref < 0` branch of
// bcf_call2bcf (bam2bcf.c:1211-1234, 1257-1283). It turns a combined
// indel bcfCall into a vcf.Variant, building REF/ALT from
// bca.IndelTypes / bca.Inscns / bca.IndelReg, and emitting the
// INDEL flag together with IDV/IMF, plus the standard INFO/DP/I16/QS,
// the bias subset present, and FORMAT/PL.
//
// refSlab is the chromosome reference (0-indexed); pos0 is the 0-based
// reference position of the indel (the position immediately PRIOR to
// the inserted/deleted span — i.e. the anchor base printed as the first
// REF nucleotide).
func bcfCall2bcfIndel(chrom string, pos0 int, refSlab []byte, call *bcfCall,
	bca *bcfCallauxIndel, fmtFlag uint32) *vcf.Variant {

	// Build REF and ALT strings (bam2bcf.c:1211-1234).
	refBase := byte('N')
	if pos0 >= 0 && pos0 < len(refSlab) {
		refBase = upperByte(refSlab[pos0])
	}
	var ref strings.Builder
	ref.WriteByte(refBase)
	for j := 0; j < bca.IndelReg; j++ {
		if pos0+1+j < len(refSlab) {
			ref.WriteByte(upperByte(refSlab[pos0+1+j]))
		}
	}
	alts := make([]string, 0, call.nAlleles-1)
	for i := 1; i < call.nAlleles && i < 4; i++ {
		a := call.alleles[i]
		if a < 0 {
			break
		}
		var b strings.Builder
		b.WriteByte(refBase)
		t := bca.IndelTypes[a]
		switch {
		case t < 0:
			// Deletion: skip the deleted reference bases.
			for j := -t; j < bca.IndelReg; j++ {
				if pos0+1+j < len(refSlab) {
					b.WriteByte(upperByte(refSlab[pos0+1+j]))
				}
			}
		case t > 0:
			// Insertion: emit the per-type consensus payload, then the
			// indelreg ref tail.
			ins := bca.Inscns[a*bca.MaxIns : (a+1)*bca.MaxIns]
			for j := 0; j < t; j++ {
				code := ins[j]
				if int(code) < 5 {
					b.WriteByte("ACGTN"[code])
				} else {
					b.WriteByte('N')
				}
			}
			for j := 0; j < bca.IndelReg; j++ {
				if pos0+1+j < len(refSlab) {
					b.WriteByte(upperByte(refSlab[pos0+1+j]))
				}
			}
		default:
			// t == 0: this should not appear among the ALT slots because
			// allele 0 carries the REF type. Skip defensively.
			continue
		}
		alts = append(alts, b.String())
	}

	info := map[string]string{}
	infoOrder := make([]string, 0, 16)

	// INDEL flag, IDV, IMF (bam2bcf.c:1257-1283).
	info["INDEL"] = ""
	infoOrder = append(infoOrder, "INDEL")
	if fmtFlag&B2BInfoIDV != 0 {
		info["IDV"] = strconv.Itoa(bca.MaxSupport)
		infoOrder = append(infoOrder, "IDV")
	}
	if fmtFlag&B2BInfoIMF != 0 {
		info["IMF"] = formatFloat32G(bca.MaxFrac)
		infoOrder = append(infoOrder, "IMF")
	}

	info["DP"] = strconv.Itoa(call.oriDepth)
	infoOrder = append(infoOrder, "DP")

	// INFO/SCR mirrors the SNP path (bam2bcf.c:1298-1299): emitted
	// before I16 when -a INFO/SCR is selected, using the shared
	// per-column tally produced by bcfCallCombine.
	if fmtFlag&B2BInfoSCR != 0 {
		info["SCR"] = strconv.Itoa(call.scrTotal)
		infoOrder = append(infoOrder, "SCR")
	}

	var i16 strings.Builder
	for j, v := range call.anno {
		if j > 0 {
			i16.WriteByte(',')
		}
		i16.WriteString(formatI16Number(v))
	}
	info["I16"] = i16.String()
	infoOrder = append(infoOrder, "I16")

	var qs strings.Builder
	for j := 0; j < call.nAlleles; j++ {
		if j > 0 {
			qs.WriteByte(',')
		}
		qs.WriteString(formatFloat32G(call.qsum[j]))
	}
	info["QS"] = qs.String()
	infoOrder = append(infoOrder, "QS")

	// Bias annotations: only with a real ALT (always true here), and
	// only when defined. Order matches the SNP path (bam2bcf.c:1306-1336):
	// VDB, SGB, RPBZ, MQBZ, MQSBZ, BQBZ, SCBZ.
	addBias := func(tag string, v float64, ok bool) {
		if ok {
			info[tag] = formatFloat32G(v)
			infoOrder = append(infoOrder, tag)
		}
	}
	addBias("VDB", call.vdb, call.vdbOK)
	addBias("SGB", call.segBias, call.sgbOK)
	addBias("RPBZ", call.mwuPos, call.rpbzOK)
	addBias("MQBZ", call.mwuMq, call.mqbzOK)
	addBias("MQSBZ", call.mwuMqs, call.mqsbzOK)
	addBias("BQBZ", call.mwuBq, call.bqbzOK)
	if fmtFlag&B2BInfoNMBZ != 0 {
		addBias("NMBZ", call.mwuNm, call.nmbzOK)
	}
	addBias("SCBZ", call.mwuSc, call.scbzOK)

	mq0f := 0.0
	if call.oriDepth > 0 {
		mq0f = float64(call.mq0) / float64(call.oriDepth)
	}
	info["MQ0F"] = formatFloat32G(mq0f)
	infoOrder = append(infoOrder, "MQ0F")

	// FORMAT: PL plus optional AD/ADF/ADR/SCR. Column order matches
	// the SNP path: AD, ADF, ADR, SCR (mirrors the bit ordering).
	format := []string{"PL"}
	emitAD := fmtFlag&B2BFmtAD != 0
	emitADF := fmtFlag&B2BFmtADF != 0
	emitADR := fmtFlag&B2BFmtADR != 0
	emitSCR := fmtFlag&B2BFmtSCR != 0
	if emitAD {
		format = append(format, "AD")
	}
	if emitADF {
		format = append(format, "ADF")
	}
	if emitADR {
		format = append(format, "ADR")
	}
	if emitSCR {
		format = append(format, "SCR")
	}
	samplesOut := make([]vcf.Sample, len(call.pl))
	for s := range call.pl {
		var pl strings.Builder
		for k, v := range call.pl[s] {
			if k > 0 {
				pl.WriteByte(',')
			}
			pl.WriteString(strconv.Itoa(v))
		}
		data := map[string]string{"PL": pl.String()}
		if emitAD && s < len(call.adf) {
			data["AD"] = formatPerAlleleSum(call.adf[s], call.adr[s])
		}
		if emitADF && s < len(call.adf) {
			data["ADF"] = formatPerAllele(call.adf[s])
		}
		if emitADR && s < len(call.adr) {
			data["ADR"] = formatPerAllele(call.adr[s])
		}
		if emitSCR && s < len(call.scr) {
			data["SCR"] = strconv.Itoa(call.scr[s])
		}
		samplesOut[s] = vcf.Sample{Data: data}
	}

	return &vcf.Variant{
		Chrom:     chrom,
		Pos:       pos0 + 1,
		ID:        ".",
		Ref:       ref.String(),
		Alt:       alts,
		Qual:      0,
		Filter:    []string{"."},
		Info:      info,
		InfoOrder: infoOrder,
		Format:    format,
		Samples:   samplesOut,
	}
}

// formatFloat32G renders v exactly as upstream renders an INFO float
// field: bcftools stores the value as a C `float` (32-bit) and the VCF
// writer prints it with C's `%g` conversion, which uses six significant
// digits with trailing zeros stripped. Rounding through float32 first
// is what makes INFO/QS, I16 and the bias tags (VDB, SGB, RPBZ, ...)
// byte-for-byte identical to upstream output.
func formatFloat32G(v float64) string {
	f := float64(float32(v))
	if math.IsInf(f, 0) || math.IsNaN(f) {
		return strconv.FormatFloat(f, 'g', 6, 64)
	}
	return strconv.FormatFloat(f, 'g', 6, 64)
}

// formatI16Number renders one I16 slot. Upstream copies the double into
// a C float before writing, so I16 goes through the same float32 + %g
// path as the other INFO floats.
func formatI16Number(v float64) string {
	return formatFloat32G(v)
}

// formatPerAllele renders a comma-separated list of integer counts, one
// per output allele, for FORMAT/ADF or FORMAT/ADR. The slice is already
// reordered to the site's allele order by bcfCallCombine.
func formatPerAllele(vals []int) string {
	if len(vals) == 0 {
		return "0"
	}
	var b strings.Builder
	for i, v := range vals {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(v))
	}
	return b.String()
}

// formatPerAlleleSum renders FORMAT/AD: the elementwise sum of the
// per-allele forward and reverse allelic depths. Both inputs are
// expected to be the same length and ordered to the site's allele list.
func formatPerAlleleSum(adf, adr []int) string {
	n := len(adf)
	if len(adr) > n {
		n = len(adr)
	}
	if n == 0 {
		return "0"
	}
	var b strings.Builder
	for i := 0; i < n; i++ {
		var v int
		if i < len(adf) {
			v += adf[i]
		}
		if i < len(adr) {
			v += adr[i]
		}
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(v))
	}
	return b.String()
}

// buildMpileupHeader builds the VCF header (metadata + sample list) for
// the output. The INFO/FORMAT lines match the SNP-relevant subset of
// upstream's mpileup header.
func buildMpileupHeader(opts MpileupOptions, chroms []string, chromLen map[string]int, samples []string) *vcf.Header {
	meta := []string{"##fileformat=VCFv4.2"}
	if !opts.NoVersion {
		meta = append(meta,
			"##bcftoolsVersion=bio_ai_experiment",
			"##bcftools_mpileupCommand=mpileup",
		)
	}
	meta = append(meta, `##FILTER=<ID=PASS,Description="All filters passed">`)
	for _, c := range chroms {
		meta = append(meta, fmt.Sprintf("##contig=<ID=%s,length=%d>", c, chromLen[c]))
	}
	// INFO/FORMAT declarations, in the same order as upstream mpileup so
	// the header matches byte-for-byte. INDEL/IDV/IMF are declared even
	// though indel calling is not yet ported: upstream always emits them.
	meta = append(meta,
		`##ALT=<ID=*,Description="Represents allele(s) other than observed.">`,
		`##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">`,
		`##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">`,
		`##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">`,
		`##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">`,
		`##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">`,
		`##INFO=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">`,
		`##INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">`,
		`##INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">`,
		`##INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">`,
	)
	// INFO/NMBZ is optional, gated by B2BInfoNMBZ. Upstream emits it
	// between MQSBZ and SCBZ (mpileup.c:816-817 + bam2bcf.c:1326).
	if opts.FmtFlag&B2BInfoNMBZ != 0 {
		meta = append(meta,
			`##INFO=<ID=NMBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better; approximate, experimental, make me localized?)">`)
	}
	meta = append(meta,
		`##INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">`,
		`##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf">`,
		`##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">`,
		`##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">`,
		`##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">`,
		`##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">`,
	)
	// Optional FORMAT tags emitted only when their `-a` bits are set. The
	// header-line text matches upstream verbatim (mpileup.c:843-849).
	if opts.FmtFlag&B2BFmtAD != 0 {
		meta = append(meta,
			`##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">`)
	}
	if opts.FmtFlag&B2BFmtADF != 0 {
		meta = append(meta,
			`##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allelic depths on the forward strand (high-quality bases)">`)
	}
	if opts.FmtFlag&B2BFmtADR != 0 {
		meta = append(meta,
			`##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allelic depths on the reverse strand (high-quality bases)">`)
	}
	// INFO/SCR + FORMAT/SCR, gated independently by their bits. Header
	// text matches upstream verbatim (mpileup.c:858-860).
	if opts.FmtFlag&B2BInfoSCR != 0 {
		meta = append(meta,
			`##INFO=<ID=SCR,Number=1,Type=Integer,Description="Number of soft-clipped reads (at high-quality bases)">`)
	}
	if opts.FmtFlag&B2BFmtSCR != 0 {
		meta = append(meta,
			`##FORMAT=<ID=SCR,Number=1,Type=Integer,Description="Per-sample number of soft-clipped reads (at high-quality bases)">`)
	}
	return &vcf.Header{MetaInfo: meta, Samples: samples}
}

// openMpileupOutput returns a variantWriter for the requested -O format
// plus a cleanup function that flushes/closes any wrapping compressor.
// The caller still owns the underlying writer.
func openMpileupOutput(out io.Writer, opts MpileupOptions, hdr *vcf.Header) (variantWriter, func(), error) {
	switch opts.OutputFormat {
	case OutputVCFGz:
		gw := gzip.NewWriter(out)
		if opts.CompressLevel > 0 {
			if g, err := gzip.NewWriterLevel(out, opts.CompressLevel); err == nil {
				gw = g
			}
		}
		return &vcfVariantWriter{vcf.NewWriter(gw, hdr)}, func() { _ = gw.Close() }, nil
	case OutputBCF:
		bw := bgzip.NewWriter(out)
		w, err := bcf.NewWriterFromVCFHeader(bw, hdr)
		if err != nil {
			_ = bw.Close()
			return nil, func() {}, err
		}
		return &bcfVariantWriter{w}, func() { _ = w.Flush(); _ = bw.Close() }, nil
	case OutputBCFUncompressed:
		w, err := bcf.NewWriterFromVCFHeader(out, hdr)
		if err != nil {
			return nil, func() {}, err
		}
		return &bcfVariantWriter{w}, func() { _ = w.Flush() }, nil
	}
	return &vcfVariantWriter{vcf.NewWriter(out, hdr)}, func() {}, nil
}
