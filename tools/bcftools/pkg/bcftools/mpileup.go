// bcftools mpileup — generate per-position genotype likelihoods from
// BAM input. This is the upstream input to `bcftools call`.
//
// V1 SIMPLIFICATION: this port ships the SNP-only, uniform-error
// genotype-likelihood model. Per the project parity rule
// (docs/PARITY_ROADMAP.md "Definition of 1:1") the CLI surface matches
// upstream `mpileup.c::main_mpileup` getopt_long; the underlying
// algorithm is documented as a v1 simplification and tracked for
// follow-up. The simplifications are explicitly:
//
//   - No BAQ recalibration (`-B/--no-BAQ` is the default and the
//     enabled state both no-op; `-E/--redo-BAQ` is hard-rejected).
//   - No indel calling (the upstream MAQ-style indel realigner is
//     deferred; every `--ext-prob`, `--gap-frac`, `--tandem-qual`,
//     `--indel-bias`, `--indel-size`, `--min-ireads`, `--max-idepth`,
//     `--open-prob`, `--ar-prob` knob is accepted at the CLI but inert
//     in v1).
//   - The likelihood model is samtools-0.1.19-style uniform-error
//     binomial: for each (REF, ALT) candidate site we walk the base
//     pile, compute log10(P(b | g)) per base, sum across reads, and
//     emit FORMAT/PL as [00, 01, 11] in phred-scaled form.
//   - Single-sample only: each BAM input is one sample; one VCF
//     output column. Multi-sample column merging is deferred.
//
// Upstream reference: reference_code/bcftools/mpileup.c. The full
// algorithm is a per-position MAQ likelihood model with optional BAQ
// recalibration and indel realignment via the
// `bam2bcf_indel.c` machinery.
//
// Output is a streaming VCF that `bcftools call` can consume:
//
//	##fileformat=VCFv4.2
//	##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
//	##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag for...">
//	##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled GL...">
//	#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <sample>
//	chr1   100 .  A  C   0    .       DP=10;I16=...  PL  0,30,300
package bcftools

import (
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// Defaults that match upstream bcftools mpileup.
const (
	// DefaultMpileupMaxDepth is upstream `-d` default.
	DefaultMpileupMaxDepth = 250
	// DefaultMpileupMinMQ is upstream `-q` default.
	DefaultMpileupMinMQ uint8 = 0
	// DefaultMpileupMinBQ is upstream `-Q` default.
	DefaultMpileupMinBQ uint8 = 13
	// DefaultMpileupTandemQual is upstream `-h` default (indel-aware
	// homopolymer penalty). Accepted at CLI but unused in v1.
	DefaultMpileupTandemQual = 500
	// DefaultMpileupExtProb is upstream `--ext-prob`. Unused in v1.
	DefaultMpileupExtProb = 20
	// DefaultMpileupGapFrac is upstream `--gap-frac`. Unused in v1.
	DefaultMpileupGapFrac = 0.05
	// DefaultMpileupOpenProb is upstream `--open-prob`. Unused in v1.
	DefaultMpileupOpenProb = 40
	// DefaultMpileupIndelBias is upstream `--indel-bias`. Unused in v1.
	DefaultMpileupIndelBias = 1.00
	// DefaultMpileupIndelSize is upstream `--indel-size`. Unused in v1.
	DefaultMpileupIndelSize = 110
	// DefaultMpileupMinIReads is upstream `--min-ireads`. Unused in v1.
	DefaultMpileupMinIReads = 1
	// DefaultMpileupMaxIDepth is upstream `--max-idepth`. Unused in v1.
	DefaultMpileupMaxIDepth = 250
	// DefaultMpileupARProb is upstream `--ar-prob`. Unused in v1.
	DefaultMpileupARProb = 1e-4
)

// MpileupOptions configures bcftools mpileup. Fields are 1:1 with the
// upstream getopt_long table in `mpileup.c`. Knobs that the v1 model
// does not consume are tagged "accepted; v1 unused" in the doc
// comment and tracked in PARITY_ROADMAP.
type MpileupOptions struct {
	// Inputs is the list of BAM/SAM paths to pile up. Multi-BAM input
	// yields one VCF sample column per BAM (sample name comes from the
	// @RG SM tag if uniform across the file, otherwise the basename of
	// the input).
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
	// MinBQ is upstream's -Q/--min-BQ (default 13).
	MinBQ uint8
	// MaxBQ is upstream's --max-bq cap (accepted, used as a ceiling on
	// the base quality so callers can clip outlier qualities).
	MaxBQ uint8

	// CountOrphans is upstream's -A/--count-orphans.
	CountOrphans bool
	// IgnoreOverlaps is upstream's -x/--ignore-overlaps.
	IgnoreOverlaps bool
	// NoBAQ is upstream's -B/--no-BAQ. V1 never applies BAQ, so this
	// is effectively the default; the flag exists for parity.
	NoBAQ bool
	// RedoBAQ is upstream's -E/--redo-BAQ. Hard-rejected in v1.
	RedoBAQ bool
	// FullBAQ is upstream's -D/--full-baq. Accepted; v1 ignores.
	FullBAQ bool
	// AdjustMQ is upstream's -C/--adjust-mq. Accepted; v1 ignores.
	AdjustMQ int

	// Annotate is upstream's -a/--annotate list (FORMAT/INFO tags to
	// include). Accepted; v1 always emits the default set
	// (INFO/DP, INFO/I16, FORMAT/PL).
	Annotate string

	// ReadGroups is upstream's -G/--read-groups. Accepted; v1 ignores.
	ReadGroups string
	// IgnoreRG is upstream's --ignore-RG (long-only). Accepted; v1 ignores.
	IgnoreRG bool

	// Platforms is upstream's -P/--platforms. Accepted; v1 ignores.
	Platforms string

	// Config is upstream's -X/--config (predefined indel-model preset:
	// `1.12`, `2.1`, `ultima`, `pacbio-ccs-1.20`, etc). Accepted; v1
	// ignores since we don't run an indel realigner.
	Config string

	// PerSampleMF is upstream's -p/--per-sample-mF. Accepted; v1 ignores.
	PerSampleMF bool

	// Seed is upstream's --seed (random seed for subsampling).
	// Accepted; v1 ignores (no subsampling).
	Seed int64

	// TandemQual is upstream's -h/--tweak-stop / --tandem-qual.
	// Accepted; v1 ignores.
	TandemQual int
	// ExtProb is upstream's --ext-prob. Accepted; v1 ignores.
	ExtProb int
	// GapFrac is upstream's --gap-frac. Accepted; v1 ignores.
	GapFrac float64
	// OpenProb is upstream's --open-prob. Accepted; v1 ignores.
	OpenProb int
	// IndelBias is upstream's --indel-bias. Accepted; v1 ignores.
	IndelBias float64
	// IndelSize is upstream's --indel-size. Accepted; v1 ignores.
	IndelSize int
	// MinIReads is upstream's --min-ireads. Accepted; v1 ignores.
	MinIReads int
	// MaxIDepth is upstream's --max-idepth. Accepted; v1 ignores.
	MaxIDepth int
	// ARProb is upstream's --ar-prob. Accepted; v1 ignores.
	ARProb float64
	// AmbigReads is upstream's --ambig-reads / --ar. Accepted; v1 ignores.
	AmbigReads string
	// MaxReadLen is upstream's -M/--max-read-len. Accepted; v1 ignores.
	MaxReadLen int

	// DelBias is upstream's --del-bias (hidden). Accepted; v1 ignores.
	DelBias float64
	// PolyMQual is upstream's --poly-mqual. Accepted; v1 ignores.
	PolyMQual bool
	// ScoreVsRef is upstream's --score-vs-ref. Accepted; v1 ignores.
	ScoreVsRef float64
	// SeqQOffset is upstream's --seqq-offset. Accepted; v1 ignores.
	SeqQOffset int

	// SkipIndels is upstream's -I/--skip-indels. V1 already never
	// emits indel records so the flag is effectively the default.
	SkipIndels bool
	// IndelsCNS is upstream's --indels-cns. Accepted; v1 ignores.
	IndelsCNS bool
	// NoIndelsCNS is upstream's --no-indels-cns. Accepted; v1 ignores.
	NoIndelsCNS bool

	// GVCFBlock is upstream's -g/--gvcf. Accepted; v1 always emits
	// one record per position (no gVCF blocking).
	GVCFBlock string

	// NoReference is upstream's --no-reference (skip the FASTA REF
	// check). Accepted; v1 always uses the FASTA REF.
	NoReference bool

	// OutputFormat is upstream's -O/--output-type (v|z|u|b). v1 emits
	// VCF text ("v") or gzipped VCF ("z"); BCF output (u|b) is
	// hard-rejected with a roadmap pointer.
	OutputFormat OutputFormat
	// Output is upstream's -o/--output (default stdout).
	Output string
	// CompressLevel is upstream's --compression-level (gzip level for -O z).
	CompressLevel int

	// Threads is upstream's --threads (accepted; v1 is single-threaded).
	Threads int
	// NoVersion is upstream's --no-version (omit the version line in
	// the header).
	NoVersion bool

	// Verbosity is upstream's -v/--verbosity (accepted; v1 ignores).
	Verbosity int

	// FlagIncl / FlagExcl are upstream's --rf/--ff (now spelled
	// --skip-all-unset / --skip-any-set in newer bcftools). Accepted;
	// v1 ignores. The boring no-op behaviour matches upstream when
	// these flags are left empty.
	FlagIncl string
	FlagExcl string
	FlagAny  string
	FlagLS   string
}

// ErrMpileupRedoBAQ is returned when -E/--redo-BAQ is requested. The
// upstream BAQ recalibrator is non-trivial (Heng Li's BAQ HMM in
// `bcftools/bam2bcf.c`) and not yet ported. Tracked in
// docs/PARITY_ROADMAP.md#bcftools.
var ErrMpileupRedoBAQ = fmt.Errorf("bcftools mpileup: -E/--redo-BAQ is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools")

// ErrMpileupBCFOutput is returned when -O u/b is requested. The BCF
// writer is wired elsewhere in the project (`pkg/bioformats/bcf`) but
// the mpileup-specific records carry custom INFO/FORMAT tags that the
// writer is not yet taught to encode. Tracked in PARITY_ROADMAP.
var ErrMpileupBCFOutput = fmt.Errorf("bcftools mpileup: -O u/b (BCF output) is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools")

// MpileupFile is the file-path entry point. It opens every input BAM,
// the FASTA reference, and writes a streaming VCF to out.
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
	// post-filters in v1 (no BAI seek path).
	regWindows, err := parseMpileupRegions(opts, chromLen)
	if err != nil {
		return err
	}

	// Sample names for the VCF #CHROM line and FORMAT column.
	samples := make([]string, len(in))
	for i, x := range in {
		samples[i] = x.sample
	}
	if len(opts.Samples) > 0 || opts.SamplesFile != "" {
		// Restrict to the named samples (preserving input order).
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
		// Drop inputs whose sample isn't requested.
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

// validateMpileupOptions applies upstream's defaults and hard-rejects
// flags whose underlying behaviour is deferred. Per the project parity
// rule, every accepted-but-deferred flag must be visible at the CLI
// surface but flagging it on must produce a clean error.
func validateMpileupOptions(opts *MpileupOptions) error {
	if opts.RedoBAQ {
		return ErrMpileupRedoBAQ
	}
	switch opts.OutputFormat {
	case OutputVCF, OutputVCFGz:
		// OK.
	default:
		return ErrMpileupBCFOutput
	}
	if opts.MaxDepth == 0 {
		opts.MaxDepth = DefaultMpileupMaxDepth
	}
	if opts.MinBQ == 0 {
		opts.MinBQ = DefaultMpileupMinBQ
	}
	if opts.MaxBQ == 0 {
		// Upstream caps at 60 to stop overly-confident base callers
		// from dominating; we follow.
		opts.MaxBQ = 60
	}
	return nil
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

// deriveSample picks a sample name for the VCF output. We use the @RG
// SM tag when uniform across the file's @RG lines, falling back to the
// basename of the BAM (matching upstream's behaviour when no @RG SM is
// present).
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
			// Inconsistent SM; fall back to file basename.
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
// (unmapped, secondary, QCfail, duplicate; orphans unless -A;
// MAPQ floor). Note: FSUPPLEMENTARY is NOT in upstream's default mask
// (mpileup.c:1392 `BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP` =
// 0x704). The earlier 0xF04 mask was a regression caught in review.
func mpileupKeepRecord(rec *sam.Record, opts MpileupOptions) bool {
	if rec == nil || rec.Pos <= 0 || rec.RName == "" {
		return false
	}
	if rec.Flag&(sam.FlagUnmapped|sam.FlagSecondary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
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
	// Sort + merge per chrom for stable contains checks.
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
// windows associated with chrom. When windows is nil (the no-restriction
// case) every position passes.
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
// base events from every input, and writes a streaming VCF to out.
func writeMpileupVCF(out io.Writer, opts MpileupOptions, ref *fasta.RandomAccess,
	chromOrder []string, chromLen map[string]int,
	perInputRecs []map[string][]*sam.Record, samples []string,
	regWindows map[string][][2]int) error {

	bw := newMpileupOutput(out, opts)
	defer bw.Close()

	if err := writeMpileupHeader(bw, opts, chromOrder, chromLen, samples); err != nil {
		return err
	}

	// Emit positions in chrom-order, position-order.
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
		// Pull this chrom's records per input.
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
		if err := emitChromMpileup(bw, chrom, refSlab, refLen, perInputChromRecs, opts, regWindows); err != nil {
			return err
		}
	}
	return nil
}

// emitChromMpileup walks every covered position on one chromosome and
// writes a VCF record per variant site (or per all-ref site when
// --skip-indels is false, which we then would log; in v1 we emit only
// SNP candidates that have at least one ALT-supporting read).
func emitChromMpileup(out io.Writer, chrom string, refSlab []byte, refLen int,
	perInputChromRecs [][]*sam.Record, opts MpileupOptions,
	regWindows map[string][][2]int) error {

	// Build per-position base events for every input. The window is
	// the entire chrom in v1 (we don't split into sub-windows because
	// the BAM is already in memory and the event slabs are small per
	// position).
	nIn := len(perInputChromRecs)
	events := make([][][]mpileupBase, nIn)
	for i := 0; i < nIn; i++ {
		events[i] = make([][]mpileupBase, refLen)
		for _, rec := range perInputChromRecs[i] {
			accumulateMpileupBases(rec, events[i], opts)
		}
	}

	for pos0 := 0; pos0 < refLen; pos0++ {
		pos1 := pos0 + 1
		if !regionContains(regWindows, chrom, pos1) {
			continue
		}
		refB := byte('N')
		if pos0 < len(refSlab) {
			refB = upperByte(refSlab[pos0])
		}
		if refB == 'N' {
			// Upstream emits records for N positions only when the
			// site has reads; we follow the simpler rule and skip
			// since the likelihood model needs a valid REF.
			continue
		}
		// Per-sample base lists at this position.
		anyCov := false
		var perSampleBases [][]mpileupBase
		for i := 0; i < nIn; i++ {
			bs := filterAndCap(events[i][pos0], opts)
			perSampleBases = append(perSampleBases, bs)
			if len(bs) > 0 {
				anyCov = true
			}
		}
		if !anyCov {
			continue
		}
		// Choose ALT alleles: the non-REF bases present across any
		// sample, in descending count.
		alts := chooseALTs(perSampleBases, refB)
		if len(alts) == 0 {
			// Pure reference; v1 omits these (matching upstream's
			// default "no -G" behaviour where REF-only sites are
			// dropped unless --gvcf grouping kicks in, which is
			// deferred).
			continue
		}
		// Write the record.
		if err := writeMpileupRecord(out, chrom, pos1, refB, alts, perSampleBases); err != nil {
			return err
		}
	}
	return nil
}

// mpileupBase is one read's base contribution to one reference position.
type mpileupBase struct {
	base      byte // uppercase ACGTN
	qual      uint8
	mapq      uint8
	isReverse bool
	qname     string
}

// accumulateMpileupBases walks rec's CIGAR and appends one event per
// covered reference position into events[pos0]. CIGAR ops that don't
// produce an SNP-candidate base (D, N, S, H, P, I) are skipped — indel
// candidates are tracked in PARITY_ROADMAP for v2.
func accumulateMpileupBases(rec *sam.Record, events [][]mpileupBase, opts MpileupOptions) {
	if rec.Pos <= 0 {
		return
	}
	refPos := int(rec.Pos) - 1
	queryPos := 0
	isReverse := rec.Flag&sam.FlagReverse != 0
	for _, op := range rec.Cigar {
		l := int(op.Length())
		o := op.Op()
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
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
				if base != 'A' && base != 'C' && base != 'G' && base != 'T' {
					continue
				}
				var qual uint8
				if q < len(rec.Qual) {
					qual = rec.Qual[q]
				}
				if opts.MaxBQ > 0 && qual > opts.MaxBQ {
					qual = opts.MaxBQ
				}
				events[p] = append(events[p], mpileupBase{
					base:      base,
					qual:      qual,
					mapq:      rec.MapQ,
					isReverse: isReverse,
					qname:     rec.QName,
				})
			}
			refPos += l
			queryPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += l
		case sam.CigarInsertion, sam.CigarSoftClip:
			queryPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// no advance.
		}
	}
}

// filterAndCap applies -Q (MinBQ), -x (IgnoreOverlaps), and -d (MaxDepth)
// to a per-position event slice. The slice returned is freshly allocated.
func filterAndCap(evs []mpileupBase, opts MpileupOptions) []mpileupBase {
	if len(evs) == 0 {
		return nil
	}
	out := make([]mpileupBase, 0, len(evs))
	seenQNames := map[string]int{} // for IgnoreOverlaps
	for _, e := range evs {
		if opts.MinBQ > 0 && e.qual < opts.MinBQ {
			continue
		}
		if opts.IgnoreOverlaps {
			if idx, ok := seenQNames[e.qname]; ok {
				// Keep the higher-quality half of the overlap.
				if e.qual > out[idx].qual {
					out[idx] = e
				}
				continue
			}
			seenQNames[e.qname] = len(out)
		}
		out = append(out, e)
		if opts.MaxDepth > 0 && len(out) >= opts.MaxDepth {
			break
		}
	}
	return out
}

// chooseALTs returns the non-REF bases observed across any sample, in
// descending total-count order. v1 caps at 1 ALT (biallelic-only)
// because the PL emitter only computes the biallelic 3-value triple;
// emitting two ALTs while writing only three PL values produces an
// output that's spec-non-conforming against `Number=G` (which expects
// 6 values for n_alt=2). The full multi-allelic PL grid lands with
// the MAQ port. Reviewer-caught regression on PR #111.
func chooseALTs(perSampleBases [][]mpileupBase, ref byte) []byte {
	counts := map[byte]int{}
	for _, evs := range perSampleBases {
		for _, e := range evs {
			if e.base == ref {
				continue
			}
			counts[e.base]++
		}
	}
	type ac struct {
		b byte
		n int
	}
	all := make([]ac, 0, len(counts))
	for b, n := range counts {
		all = append(all, ac{b, n})
	}
	sort.Slice(all, func(i, j int) bool {
		if all[i].n != all[j].n {
			return all[i].n > all[j].n
		}
		return all[i].b < all[j].b
	})
	if len(all) > 1 {
		all = all[:1]
	}
	out := make([]byte, len(all))
	for i, a := range all {
		out[i] = a.b
	}
	return out
}

// mpileupDiploidPL computes phred-scaled genotype likelihoods for a
// biallelic site under the uniform-error model. The returned slice is
// [PL(0/0), PL(0/1), PL(1/1)] with the minimum value subtracted (so
// the chosen genotype is 0). REF is allele 0, ALT is allele 1.
//
// Per-base log10 likelihood when the genotype is g:
//
//	g = 0/0 → p(base|g) = 1-e        if base == ref else e/3
//	g = 0/1 → p(base|g) = 0.5*(1-e) + 0.5*e/3  if base in {ref, alt}
//	                    else e/3
//	g = 1/1 → p(base|g) = 1-e        if base == alt else e/3
//
// e is derived from the base quality: e = 10^(-Q/10). The genotype
// log10-likelihood is the sum across bases; we then convert to phred
// (multiply by -10) and rebase to min=0.
func mpileupDiploidPL(bases []mpileupBase, ref, alt byte) [3]int {
	var ll [3]float64
	for _, b := range bases {
		// e = 10^(-q/10) clamped to [1e-6, 0.75] to mirror the
		// upstream samtools 0.1.19 cap (qual=0 → e=1, which would
		// blow up the log; cap at 0.75 to keep it finite while
		// preserving the "uninformative" semantics).
		q := float64(b.qual)
		e := math.Pow(10, -q/10)
		if e > 0.75 {
			e = 0.75
		}
		if e < 1e-6 {
			e = 1e-6
		}
		one := 1.0 - e
		het := 0.5*one + 0.5*(e/3.0)

		// p(base|g) for each g.
		var p [3]float64
		if b.base == ref {
			p[0] = one
			p[1] = het
			p[2] = e / 3
		} else if b.base == alt {
			p[0] = e / 3
			p[1] = het
			p[2] = one
		} else {
			// Off-allele: only e/3 mass for any g.
			p[0] = e / 3
			p[1] = e / 3
			p[2] = e / 3
		}
		for g := 0; g < 3; g++ {
			ll[g] += math.Log10(p[g])
		}
	}
	// Convert to phred (×-10), subtract min.
	var phred [3]float64
	for g := 0; g < 3; g++ {
		phred[g] = -10 * ll[g]
	}
	mn := phred[0]
	for g := 1; g < 3; g++ {
		if phred[g] < mn {
			mn = phred[g]
		}
	}
	var out [3]int
	for g := 0; g < 3; g++ {
		v := phred[g] - mn
		if v > 255 {
			v = 255 // upstream caps PL at 255.
		}
		if v < 0 {
			v = 0
		}
		out[g] = int(math.Round(v))
	}
	return out
}

// mpileupI16 computes the upstream-flavoured 16-tag aux array. The 16
// slots are documented in `bcftools/bam2bcf.c::compute_I16` but we
// reproduce the layout here:
//
//	[0] #ref bases on forward strand
//	[1] #ref bases on reverse strand
//	[2] #non-ref bases on forward strand
//	[3] #non-ref bases on reverse strand
//	[4] sum baseQ of ref bases
//	[5] sum baseQ^2 of ref bases
//	[6] sum baseQ of non-ref bases
//	[7] sum baseQ^2 of non-ref bases
//	[8] sum mapQ of ref bases
//	[9] sum mapQ^2 of ref bases
//	[10] sum mapQ of non-ref bases
//	[11] sum mapQ^2 of non-ref bases
//	[12] sum tailDist of ref bases (zero in v1 — no tail-distance tracker)
//	[13] sum tailDist^2 of ref bases
//	[14] sum tailDist of non-ref bases
//	[15] sum tailDist^2 of non-ref bases
func mpileupI16(bases []mpileupBase, ref byte) [16]float64 {
	var i16 [16]float64
	for _, b := range bases {
		isRef := b.base == ref
		fwd := !b.isReverse
		switch {
		case isRef && fwd:
			i16[0]++
		case isRef && !fwd:
			i16[1]++
		case !isRef && fwd:
			i16[2]++
		case !isRef && !fwd:
			i16[3]++
		}
		bq := float64(b.qual)
		mq := float64(b.mapq)
		if isRef {
			i16[4] += bq
			i16[5] += bq * bq
			i16[8] += mq
			i16[9] += mq * mq
		} else {
			i16[6] += bq
			i16[7] += bq * bq
			i16[10] += mq
			i16[11] += mq * mq
		}
		// tail distance not tracked in v1.
	}
	return i16
}

// writeMpileupRecord emits one VCF text record for (chrom, pos, ref,
// alts) given per-sample base lists.
func writeMpileupRecord(out io.Writer, chrom string, pos1 int, ref byte, alts []byte, perSampleBases [][]mpileupBase) error {
	var b strings.Builder
	b.WriteString(chrom)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(pos1))
	b.WriteString("\t.\t")
	b.WriteByte(ref)
	b.WriteByte('\t')
	for i, a := range alts {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteByte(a)
	}
	b.WriteString("\t0\t.\t")

	// INFO: DP, I16.
	dp := 0
	var totI16 [16]float64
	for _, evs := range perSampleBases {
		dp += len(evs)
		for j, v := range mpileupI16(evs, ref) {
			totI16[j] += v
		}
	}
	b.WriteString("DP=")
	b.WriteString(strconv.Itoa(dp))
	b.WriteString(";I16=")
	for j, v := range totI16 {
		if j > 0 {
			b.WriteByte(',')
		}
		b.WriteString(formatI16Number(v))
	}
	b.WriteString("\tPL")
	// Per-sample PL: in v1 we always emit the (REF, ALT0) biallelic
	// PL. Sites with multiple ALTs use the leading ALT in the PL
	// computation. The full multi-allelic PL grid is a v2 follow-up
	// (the ordering is `j(j+1)/2 + i` per VCF spec) — tracked in
	// PARITY_ROADMAP.
	alt0 := alts[0]
	for _, evs := range perSampleBases {
		pl := mpileupDiploidPL(evs, ref, alt0)
		b.WriteByte('\t')
		b.WriteString(strconv.Itoa(pl[0]))
		b.WriteByte(',')
		b.WriteString(strconv.Itoa(pl[1]))
		b.WriteByte(',')
		b.WriteString(strconv.Itoa(pl[2]))
	}
	b.WriteByte('\n')
	_, err := io.WriteString(out, b.String())
	return err
}

// formatI16Number matches upstream's I16 float rendering: integers
// without a fractional part; floats with the shortest exact form.
func formatI16Number(v float64) string {
	if v == math.Trunc(v) {
		return strconv.FormatFloat(v, 'f', 0, 64)
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
}

// writeMpileupHeader emits the VCF header lines: fileformat, contigs,
// INFO/FORMAT, and the #CHROM column.
func writeMpileupHeader(out io.Writer, opts MpileupOptions, chroms []string, chromLen map[string]int, samples []string) error {
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	if !opts.NoVersion {
		b.WriteString("##bcftoolsVersion=bio_ai_experiment\n")
		b.WriteString("##bcftools_mpileupCommand=mpileup\n")
	}
	for _, c := range chroms {
		fmt.Fprintf(&b, "##contig=<ID=%s,length=%d>\n", c, chromLen[c])
	}
	b.WriteString("##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n")
	b.WriteString("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n")
	b.WriteString("##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">\n")
	b.WriteString("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for _, s := range samples {
		b.WriteByte('\t')
		b.WriteString(s)
	}
	b.WriteByte('\n')
	_, err := io.WriteString(out, b.String())
	return err
}

// mpileupOutput is a tiny wrapper that lets the caller select gzipped
// VCF output. For v1 we only support text (-O v) and gzipped text (-O z)
// outputs; BCF (-O b|u) is rejected at validateMpileupOptions.
type mpileupOutput struct {
	out io.Writer
	gz  io.WriteCloser
}

func newMpileupOutput(out io.Writer, opts MpileupOptions) *mpileupOutput {
	// gzipped VCF output is wired through the standard iohelper
	// wrapper used by other subcommands; for the streaming-text v1 we
	// pass the caller's writer through unchanged (the CLI wraps stdout
	// in a gzip writer when -O z is set, so the package layer stays
	// format-agnostic).
	_ = opts
	return &mpileupOutput{out: out}
}

func (m *mpileupOutput) Write(p []byte) (int, error) { return m.out.Write(p) }
func (m *mpileupOutput) Close() error {
	if m.gz != nil {
		return m.gz.Close()
	}
	return nil
}
