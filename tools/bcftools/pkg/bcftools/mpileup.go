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
	// DefaultMpileupMinIReads is upstream `--min-ireads`. Unused.
	DefaultMpileupMinIReads = 1
	// DefaultMpileupMaxIDepth is upstream `--max-idepth`. Unused.
	DefaultMpileupMaxIDepth = 250
	// DefaultMpileupARProb is upstream `--ar-prob`. Unused.
	DefaultMpileupARProb = 1e-4
	// mpileupTheta is upstream's CALL_DEFTHETA (bam2bcf.c:39); the errmod
	// depth-correlation parameter is 1 - theta.
	mpileupTheta = 0.83
)

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
	// include). Accepted; the default set (INFO/DP, I16, QS, MQ0F,
	// FORMAT/PL) is always emitted.
	Annotate string

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
	// AmbigReads is upstream's --ambig-reads / --ar. Accepted; ignored.
	AmbigReads string
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
	// IndelsCNS is upstream's --indels-cns. Accepted; ignored.
	IndelsCNS bool
	// NoIndelsCNS is upstream's --no-indels-cns. Accepted; ignored.
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

	// FlagIncl / FlagExcl are upstream's --rf/--ff. Accepted; ignored.
	FlagIncl string
	FlagExcl string
	FlagAny  string
	FlagLS   string
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
// floor). FSUPPLEMENTARY is NOT in upstream's default mask
// (mpileup.c:1392 BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP).
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
	em := ErrmodInit(1.0 - mpileupTheta)

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
		if err := emitChromMpileup(w, em, chrom, refSlab, refLen, perInputChromRecs, opts, regWindows); err != nil {
			return err
		}
	}
	return w.Flush()
}

// emitChromMpileup walks every covered position on one chromosome and
// writes one record per position that has read coverage. Unlike the
// pre-MAQ port, this emits a record for every covered position (not
// only SNP candidates) with `<*>` as the unseen allele.
func emitChromMpileup(w variantWriter, em *Errmod, chrom string, refSlab []byte, refLen int,
	perInputChromRecs [][]*sam.Record, opts MpileupOptions,
	regWindows map[string][][2]int) error {

	nIn := len(perInputChromRecs)
	// BAQ realignment. Upstream mpileup.c runs sam_prob_realn on each
	// mapped read against the chromosome reference (BAQ on by default);
	// in apply mode it lowers rec.Qual in place. applyMpileupBAQ ports
	// mplp_realn's column-gated decision and edits rec.Qual before
	// accumulateMpileupBases reads the (now BAQ-adjusted) qualities.
	// -B/--no-BAQ skips it entirely.
	if !opts.NoBAQ {
		applyMpileupBAQ(perInputChromRecs, refSlab, opts)
	}

	// events[input][pos0] is the pileup column for one input at one
	// reference position.
	events := make([][][]pileupBase, nIn)
	for i := 0; i < nIn; i++ {
		events[i] = make([][]pileupBase, refLen)
		for _, rec := range perInputChromRecs[i] {
			accumulateMpileupBases(rec, events[i])
		}
	}

	calls := make([]bcfCallret, nIn)
	for pos0 := 0; pos0 < refLen; pos0++ {
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
			pile := filterMpileupPile(events[i][pos0], opts)
			if len(pile) > 0 {
				anyCov = true
			}
			bcfCallGlfgen(pile, ref4, opts, em, &calls[i])
		}
		if !anyCov {
			continue
		}
		call := bcfCallCombine(calls, ref4)
		v := bcfCall2bcf(chrom, pos1, refB, &call)
		if err := w.Write(v); err != nil {
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
func applyMpileupBAQ(perInputChromRecs [][]*sam.Record, refSlab []byte, opts MpileupOptions) {
	baqFlag := mpileupBAQFlag(opts)
	// max_read_len: upstream default is 500 unless -M overrides it.
	maxReadLen := opts.MaxReadLen
	if maxReadLen <= 0 {
		maxReadLen = 500
	}

	// Build per-read heuristic info and an interval index keyed by
	// covered reference position.
	var infos []*mpileupReadBAQInfo
	maxPos := 0
	for _, recs := range perInputChromRecs {
		for _, rec := range recs {
			if rec.IsUnmapped() || len(rec.Cigar) == 0 {
				continue
			}
			info := mpileupBuildBAQInfo(rec)
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
			info.realigned = true
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
// (D, N, S, H, P, I) are skipped — indel candidates are slice-4 work.
func accumulateMpileupBases(rec *sam.Record, events [][]pileupBase) {
	if rec.Pos <= 0 {
		return
	}
	refPos := int(rec.Pos) - 1
	queryPos := 0
	isReverse := rec.Flag&sam.FlagReverse != 0
	qlen := len(rec.Seq)
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
				b4 := seqNt16Int[baseToNt16(base)]
				var rawQual uint8
				if q < len(rec.Qual) {
					rawQual = rec.Qual[q]
				}
				prevQ := -1
				if q > 0 && q-1 < len(rec.Qual) {
					prevQ = int(rec.Qual[q-1])
				}
				nextQ := -1
				if q+1 < qlen && q+1 < len(rec.Qual) {
					nextQ = int(rec.Qual[q+1])
				}
				events[p] = append(events[p], pileupBase{
					base4:   b4,
					rawQual: rawQual,
					prevQ:   prevQ,
					nextQ:   nextQ,
					mapq:    rec.MapQ,
					reverse: isReverse,
					qpos:    q,
					qlen:    qlen,
					qname:   rec.QName,
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

// filterMpileupPile applies -x (IgnoreOverlaps) and -d (MaxDepth) to a
// per-position pileup column. -Q/--min-BQ is applied inside glfgen
// (after the delta_baseQ cap), matching upstream's ordering.
func filterMpileupPile(evs []pileupBase, opts MpileupOptions) []pileupBase {
	if len(evs) == 0 {
		return nil
	}
	out := make([]pileupBase, 0, len(evs))
	var seenQNames map[string]int
	if opts.IgnoreOverlaps {
		seenQNames = make(map[string]int, len(evs))
	}
	for _, e := range evs {
		if opts.IgnoreOverlaps {
			if idx, ok := seenQNames[e.qname]; ok {
				// Keep the higher-quality half of the overlapping pair.
				if e.rawQual > out[idx].rawQual {
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

// bcfCall2bcf is the Go port of bcf_call2bcf (bam2bcf.c:1200) for the
// SNP path. It turns a combined bcfCall into a vcf.Variant: REF, the
// ALT alleles (including the `<*>` unseen allele), QUAL=0, INFO/DP/I16/
// QS/MQ0F and FORMAT/PL. The bias INFO tags (VDB/SGB/RPBZ/...) are
// slice-4 work and deliberately omitted.
func bcfCall2bcf(chrom string, pos1 int, refB byte, call *bcfCall) *vcf.Variant {
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
	// quality sum); the `<*>` allele's qsum is 0.
	var qs strings.Builder
	for j := 0; j < call.nAlleles; j++ {
		if j > 0 {
			qs.WriteByte(',')
		}
		qs.WriteString(formatQSNumber(call.qsum[j]))
	}
	info["QS"] = qs.String()
	infoOrder = append(infoOrder, "QS")

	mq0f := 0.0
	if call.oriDepth > 0 {
		mq0f = float64(call.mq0) / float64(call.oriDepth)
	}
	info["MQ0F"] = formatQSNumber(mq0f)
	infoOrder = append(infoOrder, "MQ0F")

	// FORMAT/PL — one upper-triangle grid per sample.
	format := []string{"PL"}
	samplesOut := make([]vcf.Sample, len(call.pl))
	for s := range call.pl {
		var pl strings.Builder
		for k, v := range call.pl[s] {
			if k > 0 {
				pl.WriteByte(',')
			}
			pl.WriteString(strconv.Itoa(v))
		}
		samplesOut[s] = vcf.Sample{Data: map[string]string{"PL": pl.String()}}
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

// formatI16Number matches upstream's I16 float rendering: integers
// without a fractional part, otherwise the shortest exact form.
func formatI16Number(v float64) string {
	if v == math.Trunc(v) && !math.IsInf(v, 0) {
		return strconv.FormatFloat(v, 'f', 0, 64)
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
}

// formatQSNumber renders an INFO/QS or MQ0F float: whole numbers print
// as integers, otherwise the shortest exact decimal form is used.
func formatQSNumber(v float64) string {
	if v == math.Trunc(v) && !math.IsInf(v, 0) {
		return strconv.FormatFloat(v, 'f', 0, 64)
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
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
	meta = append(meta,
		`##ALT=<ID=*,Description="Represents allele(s) other than observed.">`,
		`##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">`,
		`##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">`,
		`##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">`,
		`##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">`,
		`##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">`,
	)
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
