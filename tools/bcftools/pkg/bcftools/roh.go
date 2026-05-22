package bcftools

// bcftools roh — runs-of-autozygosity HMM (port).
//
// Upstream reference: reference_code/bcftools/vcfroh.c and HMM.c. The
// port runs a 2-state HMM (HW = Hardy-Weinberg, AZ = autozygous) per
// sample and emits two tables:
//
//	ST  per-site:    Sample, CHROM, POS, State (0/1), Quality
//	RG  per-region:  Sample, CHROM, Start, End, Length, NumMarkers, Quality
//
// AZ-state emission (upstream's pdg[0]*(1-f) + pdg[2]*f) and HW-state
// emission (pdg[0]*(1-f)^2 + 2*pdg[1]*f*(1-f) + pdg[2]*f^2) are
// computed verbatim from the per-genotype likelihood vector pdg.
//
// Transition probabilities scale with the distance between consecutive
// markers. The base matrix is built from the per-bp parameters
// HW->AZ and AZ->HW; the HMM engine (hmm.go) precomputes 10000 matrix
// powers so that a gap of d bp between two sites uses base^d. When a
// genetic map (-G/--genetic-map) or a constant recombination rate
// (-M/--rec-rate) is supplied the off-diagonal transition entries are
// further scaled by the cross-over probability for the interval.
//
// Quality scores are forward-backward phred scores
// (phred_score(1 - posterior)), matching upstream exactly. The
// per-region quality is the mean of the per-site qualities over the
// run.

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// RohOptions controls the behaviour of Roh / RohFile.
type RohOptions struct {
	// AFTag picks an INFO tag for allele frequency. Empty defaults
	// to "AF".
	AFTag string
	// AFDflt is the --AF-dflt value. When nil, missing AF means
	// "skip the site". When non-nil, missing AF (or AF==0) falls back
	// to this value.
	AFDflt *float64
	// AFFile is the upstream --AF-file argument (CHR\tPOS\tREF,ALT\tAF).
	AFFile string

	// GTsOnly is upstream's -G/--GTs-only FLOAT: when >0, hard-GT
	// mode is enabled and the floating-point error probability is
	// pow(10, -GTsOnly/10).
	GTsOnly float64

	// IgnoreHomRef skips hom-ref genotypes (-i/--ignore-homref).
	IgnoreHomRef bool
	// IncludeNoalt allows records with no ALT (default skip).
	IncludeNoalt bool
	// SkipIndels drops indel records (-I/--skip-indels).
	SkipIndels bool

	// HWtoAZ / AZtoHW transition probabilities per bp. Defaults are
	// upstream's literal values: 6.7e-8 and 5e-9.
	HWtoAZ float64
	AZtoHW float64

	// OutputTypes controls which output sections are emitted. Each
	// rune is one of 's' (per-site), 'r' (regions). Default "sr".
	OutputTypes string

	// Samples / Regions / Targets / Include / Exclude — post-filters.
	Samples     []string
	SamplesFile string
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
	IncludeExpr string
	ExcludeExpr string

	// EstimateAF mirrors upstream `-e/--estimate-AF`: an optional
	// "GT," / "PL," tag prefix followed by "-" (all samples) or a
	// comma-separated sample list. The frequency is then estimated
	// from the genotypes (or PLs) of those samples.
	EstimateAF string

	// BufferSize mirrors upstream -b/--buffer-size: "INT" or
	// "INT,INT" (buffer size and overlap). A negative first value is
	// a memory budget in MB. Empty means unlimited.
	BufferSize string
	// GeneticMap is the -m/--genetic-map IMPUTE2-format file (the
	// "{CHROM}" placeholder is honoured).
	GeneticMap string
	// RecRate is the -M/--rec-rate constant recombination rate per bp.
	RecRate float64
	// ViterbiTraining is the -V/--viterbi-training convergence
	// threshold; when >0 the transition probabilities are
	// re-estimated by Baum-Welch.
	ViterbiTraining float64

	RegionsOverlap int
	TargetsOverlap int
}

// RohSite captures a single per-site state assignment.
type RohSite struct {
	Sample string
	Chrom  string
	Pos    int
	State  int     // 0 = HW, 1 = AZ
	Qual   float64 // forward-backward phred score
	AF     float64
}

// RohRegion captures a contiguous AZ run for one sample. Qual is the
// mean of the per-site forward-backward phred scores over the run,
// matching upstream's RG quality column.
type RohRegion struct {
	Sample     string
	Chrom      string
	Start      int
	End        int
	Length     int
	NumMarkers int
	Qual       float64
}

// RohResult is the full pile of (per-site, per-region) outputs.
type RohResult struct {
	Sites   []RohSite
	Regions []RohRegion
}

const (
	stateHW = 0
	stateAZ = 1
)

// DefaultHWtoAZ matches upstream's `HW->AZ = 6.7e-8` per bp.
const DefaultHWtoAZ = 6.7e-8

// DefaultAZtoHW matches upstream's `AZ->HW = 5e-9` per bp.
const DefaultAZtoHW = 5e-9

// RohFile is the file-aware entry point used by the CLI.
func RohFile(path string, out io.Writer, opts RohOptions) (RohResult, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return RohResult{}, fmt.Errorf("bcftools roh: open %s: %w", path, err)
	}
	defer in.Close()
	return Roh(in, out, opts)
}

// rohMarker is one accepted site with its emission probabilities.
type rohMarker struct {
	chrom string
	pos   int // 1-based
	eHW   float64
	eAZ   float64
}

// Roh streams VCF input through the HMM and writes the requested
// output sections to out.
func Roh(in io.Reader, out io.Writer, opts RohOptions) (RohResult, error) {
	if opts.HWtoAZ <= 0 {
		opts.HWtoAZ = DefaultHWtoAZ
	}
	if opts.AZtoHW <= 0 {
		opts.AZtoHW = DefaultAZtoHW
	}
	if opts.OutputTypes == "" {
		opts.OutputTypes = "sr"
	}
	if opts.SamplesFile != "" {
		extra, err := LoadSamplesFile(opts.SamplesFile)
		if err != nil {
			return RohResult{}, fmt.Errorf("bcftools roh: samples-file: %w", err)
		}
		opts.Samples = append(opts.Samples, extra...)
	}
	if opts.RegionsFile != "" {
		extra, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return RohResult{}, fmt.Errorf("bcftools roh: regions-file: %w", err)
		}
		opts.Regions = append(opts.Regions, extra...)
	}
	if opts.TargetsFile != "" {
		extra, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return RohResult{}, fmt.Errorf("bcftools roh: targets-file: %w", err)
		}
		opts.Targets = append(opts.Targets, extra...)
	}
	if opts.AFTag == "" {
		opts.AFTag = "AF"
	}

	// Upstream `roh` default mode scores emissions from FORMAT/PL
	// genotype likelihoods; `-G/--GTs-only` switches to hard-GT mode.
	// PL-based emission is a deliberate, documented deferral in this
	// port (see docs/PARITY_ROADMAP.md): only hard-GT mode is
	// implemented. To keep the deferral honest rather than silently
	// decoding from hard GT, reject a non-`-G` invocation with
	// upstream init_data's exact diagnostic (vcfroh.c:151-158).
	if opts.GTsOnly <= 0 {
		return RohResult{}, fmt.Errorf("Error: The FORMAT/PL tag not found in the header, consider running with -G")
	}

	// Upstream rejects -b with -V (vcfroh.c:1255).
	if opts.ViterbiTraining > 0 && opts.BufferSize != "" {
		return RohResult{}, fmt.Errorf("Error: cannot use -b with -V")
	}

	hdr, vars, err := readAllVariants(in)
	if err != nil {
		return RohResult{}, fmt.Errorf("bcftools roh: %w", err)
	}

	samples := selectSamples(hdr.Samples, opts.Samples)
	if len(samples) == 0 {
		return RohResult{}, fmt.Errorf("bcftools roh: no samples to analyse")
	}
	sampleIdx := indexSamples(hdr)

	afFromFile, err := loadAFFile(opts.AFFile)
	if err != nil {
		return RohResult{}, err
	}

	bufMax, bufOlap, err := parseBufferSize(opts.BufferSize, len(samples))
	if err != nil {
		return RohResult{}, err
	}

	// AF estimation cohort (--estimate-AF).
	afEst, afEstActive, err := parseEstimateAF(opts.EstimateAF, hdr, sampleIdx)
	if err != nil {
		return RohResult{}, err
	}

	hardGTErr := 0.0
	if opts.GTsOnly > 0 {
		hardGTErr = math.Pow(10, -opts.GTsOnly/10.0)
	}

	// Build per-sample marker lists. Markers are split by chromosome.
	type chromBlock struct {
		chrom   string
		markers []rohMarker
	}
	perSample := make([][]chromBlock, len(samples))

	prevChrom := ""
	prevPos := -1
	for _, v := range vars {
		if len(opts.Regions) > 0 && !regionMatches(v, opts.Regions) {
			continue
		}
		if len(opts.Targets) > 0 && !regionMatches(v, opts.Targets) {
			continue
		}
		if opts.SkipIndels && isIndel(v) {
			continue
		}
		// Count real ALT alleles, excluding the symbolic <*>/<NON_REF>.
		nalt := 0
		ial := 0
		for i, a := range v.Alt {
			if a == "" || a == "." {
				continue
			}
			if a == "<*>" || a == "<NON_REF>" {
				continue
			}
			nalt++
			if ial == 0 {
				ial = i + 1
			}
		}
		if nalt == 0 {
			if !opts.IncludeNoalt {
				continue
			}
		} else if nalt > 1 {
			// Multiallelic sites are skipped (upstream considers only
			// biallelic records).
			continue
		}
		// Skip duplicate positions on the same chromosome.
		if v.Chrom == prevChrom && v.Pos == prevPos {
			continue
		}
		prevChrom, prevPos = v.Chrom, v.Pos

		af, ok := resolveAF(v, ial, opts, afFromFile, afEst, afEstActive, sampleIdx, hardGTErr)
		if !ok {
			continue
		}

		for si, sample := range samples {
			idx, ok := sampleIdx[sample]
			if !ok {
				continue
			}
			cls, ok := rohGTClass(v, idx, ial)
			if !ok {
				continue
			}
			pdg := pdgFromClass(cls, hardGTErr)
			sum := pdg[0] + pdg[1] + pdg[2]
			if sum == 0 {
				continue
			}
			for i := range pdg {
				pdg[i] /= sum
			}
			if opts.IgnoreHomRef && pdg[0] > 0.99 {
				continue
			}
			eAZ := pdg[0]*(1-af) + pdg[2]*af
			eHW := pdg[0]*(1-af)*(1-af) + 2*pdg[1]*(1-af)*af + pdg[2]*af*af

			blk := &perSample[si]
			if len(*blk) == 0 || (*blk)[len(*blk)-1].chrom != v.Chrom {
				*blk = append(*blk, chromBlock{chrom: v.Chrom})
			}
			cur := &(*blk)[len(*blk)-1]
			cur.markers = append(cur.markers, rohMarker{
				chrom: v.Chrom, pos: v.Pos, eHW: eHW, eAZ: eAZ,
			})
		}
	}

	// Genetic map / recombination-rate transition modifier.
	gmap, err := loadGenMap(opts.GeneticMap)
	if err != nil {
		return RohResult{}, err
	}

	result := RohResult{}
	for si, sample := range samples {
		blocks := perSample[si]
		for _, blk := range blocks {
			if len(blk.markers) == 0 {
				continue
			}
			runRohSample(&result, sample, blk.chrom, blk.markers, opts, gmap, bufMax, bufOlap)
		}
	}

	if err := writeRoh(out, result, opts.OutputTypes); err != nil {
		return result, err
	}
	return result, nil
}

// baseTprob builds the 2x2 base transition matrix (row-major) from the
// per-bp HW->AZ and AZ->HW parameters.
func baseTprob(t2AZ, t2HW float64) []float64 {
	m := make([]float64, 4)
	matSet(m, 2, stateHW, stateHW, 1-t2AZ)
	matSet(m, 2, stateHW, stateAZ, t2HW)
	matSet(m, 2, stateAZ, stateHW, t2AZ)
	matSet(m, 2, stateAZ, stateAZ, 1-t2HW)
	return m
}

// runRohSample runs the HMM over one chromosome block for one sample
// and appends the per-site and per-region results.
func runRohSample(result *RohResult, sample, chrom string, markers []rohMarker, opts RohOptions, gmap *genMap, bufMax, bufOlap int) {
	n := len(markers)
	eprobs := make([]float64, n*2)
	sites := make([]uint32, n)
	for i, m := range markers {
		eprobs[i*2+stateHW] = m.eHW
		eprobs[i*2+stateAZ] = m.eAZ
		sites[i] = uint32(m.pos - 1) // 0-based, as upstream
	}

	h := newHMM(2, baseTprob(opts.HWtoAZ, opts.AZtoHW), 10000)
	tmod := newTprobModifier(opts, gmap, chrom)
	if tmod != nil {
		h.setTprob = tmod.apply
	}

	if opts.ViterbiTraining > 0 {
		runRohViterbiTraining(result, sample, chrom, markers, eprobs, sites, h, opts)
		return
	}

	// Single-pass Viterbi + forward-backward. With --buffer-size the
	// stream is processed in overlapping windows: the HMM state is
	// snapshotted just before the overlap region and restored for the
	// next window, so the buffered decode matches the unbuffered one.
	acc := &rohRegionAccum{}
	if bufMax <= 0 {
		if tmod != nil {
			tmod.reset()
		}
		h.runViterbi(n, eprobs, sites)
		h.runFwdBwd(n, eprobs, sites)
		appendRohSites(result, acc, sample, chrom, markers, h, opts)
		acc.finish(result)
		return
	}

	pos := 0
	if tmod != nil {
		tmod.reset()
	}
	for pos < n {
		end := n
		if n-pos > bufMax {
			end = pos + bufMax
		}
		segN := end - pos
		emit := segN
		if end < n && bufOlap > 0 && segN > bufOlap {
			emit = segN - bufOlap
			h.requestSnapshot(sites[end-bufOlap-1])
		} else {
			h.snapAtPos = 0
		}
		segE := eprobs[pos*2 : end*2]
		segS := sites[pos:end]
		h.runViterbi(segN, segE, segS)
		h.runFwdBwd(segN, segE, segS)
		appendRohSites(result, acc, sample, chrom, markers[pos:pos+emit], h, opts)

		if end < n && bufOlap > 0 {
			h.restoreSnapshot()
			pos = end - bufOlap
		} else {
			pos = end
		}
	}
	acc.finish(result)
}

// rohRegionAccum tracks the AZ run currently open for one sample so
// that a run can span buffer-flush boundaries, mirroring upstream's
// persistent smpl->rg struct.
type rohRegionAccum struct {
	open    bool
	sample  string
	chrom   string
	start   int
	end     int
	markers int
	qual    float64
}

// finish flushes any still-open region into result.
func (a *rohRegionAccum) finish(result *RohResult) {
	if !a.open {
		return
	}
	result.Regions = append(result.Regions, RohRegion{
		Sample:     a.sample,
		Chrom:      a.chrom,
		Start:      a.start,
		End:        a.end,
		Length:     a.end - a.start + 1,
		NumMarkers: a.markers,
		Qual:       a.qual / float64(a.markers),
	})
	a.open = false
}

// appendRohSites emits per-site states and rolls up AZ regions for one
// decoded segment. The accumulator acc carries an open AZ run across
// segment boundaries so buffered decoding produces the same regions as
// the unbuffered pass.
func appendRohSites(result *RohResult, acc *rohRegionAccum, sample, chrom string, markers []rohMarker, h *hmm, opts RohOptions) {
	wantST := strings.ContainsRune(opts.OutputTypes, 's')
	wantRG := strings.ContainsRune(opts.OutputTypes, 'r')

	for i, m := range markers {
		state := int(h.vpath[i])
		posterior := h.fwd[i*2+state]
		qual := phredScore(1.0 - posterior)
		if wantST {
			result.Sites = append(result.Sites, RohSite{
				Sample: sample, Chrom: chrom, Pos: m.pos,
				State: state, Qual: qual,
			})
		}
		if !wantRG {
			continue
		}
		if state == stateAZ {
			if !acc.open {
				acc.open = true
				acc.sample = sample
				acc.chrom = chrom
				acc.start = m.pos
				acc.markers = 0
				acc.qual = 0
			}
			acc.end = m.pos
			acc.markers++
			acc.qual += qual
		} else if acc.open {
			acc.finish(result)
		}
	}
}

// runRohViterbiTraining performs Baum-Welch re-estimation of the
// transition probabilities, then decodes with the converged matrix.
// It mirrors flush_viterbi's --viterbi-training branch.
func runRohViterbiTraining(result *RohResult, sample, chrom string, markers []rohMarker, eprobs []float64, sites []uint32, h *hmm, opts RohOptions) {
	h.setTprobMatrix(baseTprob(opts.HWtoAZ, opts.AZtoHW), 10000)
	for {
		cur := h.tprobMatrix()
		t2azPrev := matAt(cur, 2, stateAZ, stateHW)
		t2hwPrev := matAt(cur, 2, stateHW, stateAZ)
		next := h.runBaumWelch(len(markers), eprobs, sites)
		nextCopy := make([]float64, 4)
		copy(nextCopy, next)
		h.setTprobMatrix(nextCopy, 10000)
		deltaAZ := math.Abs(matAt(nextCopy, 2, stateAZ, stateHW) - t2azPrev)
		deltaHW := math.Abs(matAt(nextCopy, 2, stateHW, stateAZ) - t2hwPrev)
		if deltaAZ <= opts.ViterbiTraining && deltaHW <= opts.ViterbiTraining {
			break
		}
	}
	h.runViterbi(len(markers), eprobs, sites)
	h.runFwdBwd(len(markers), eprobs, sites)
	acc := &rohRegionAccum{}
	appendRohSites(result, acc, sample, chrom, markers, h, opts)
	acc.finish(result)
}

// AFDfltPtr returns a *float64 pointing at the given value.
func AFDfltPtr(v float64) *float64 { return &v }

// rohGTClass classifies the sample genotype for the site as hom-ref
// (0), het (1) or hom-alt (2). The classification depends only on
// whether the two alleles equal the reference, mirroring upstream's
// fake_PLs branch. The second return value is false for missing or
// non-diploid genotypes.
func rohGTClass(v *vcf.Variant, idx, ial int) (int, bool) {
	if v == nil || idx < 0 || idx >= len(v.Samples) {
		return 0, false
	}
	gt, ok := v.Samples[idx].Data["GT"]
	if !ok || gt == "" || gt == "." {
		return 0, false
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) != 2 {
		return 0, false
	}
	a, errA := strconv.Atoi(parts[0])
	b, errB := strconv.Atoi(parts[1])
	if errA != nil || errB != nil || parts[0] == "." || parts[1] == "." {
		return 0, false
	}
	if a != b {
		return 1, true // het
	}
	if a == 0 {
		return 0, true // hom-ref
	}
	return 2, true // hom-alt
}

// pdgFromClass collapses a genotype class into a 3-vector P(D|G). With
// err>0 the mass is redistributed per upstream vcfroh.c's fake_PLs:
//
//	HOM-REF  PDG = [1 - err - err², err, err²]
//	HET      PDG = [err, 1 - 2*err, err]
//	HOM-ALT  PDG = [err², err, 1 - err - err²]
//
// With err==0 the result is one-hot.
func pdgFromClass(class int, err float64) [3]float64 {
	if err <= 0 {
		switch class {
		case 0:
			return [3]float64{1, 0, 0}
		case 1:
			return [3]float64{0, 1, 0}
		case 2:
			return [3]float64{0, 0, 1}
		}
		return [3]float64{}
	}
	err2 := err * err
	switch class {
	case 0:
		return [3]float64{1 - err - err2, err, err2}
	case 1:
		return [3]float64{err, 1 - 2*err, err}
	case 2:
		return [3]float64{err2, err, 1 - err - err2}
	}
	return [3]float64{}
}

// pdgFromDose is retained for backward compatibility with callers and
// tests that classify by the alt-allele dosage (0/1/2).
func pdgFromDose(dose int, err float64) [3]float64 {
	return pdgFromClass(dose, err)
}

func isIndel(v *vcf.Variant) bool {
	if len(v.Ref) != 1 {
		return true
	}
	for _, a := range v.Alt {
		if a == "<*>" || a == "<NON_REF>" || a == "." || a == "" {
			continue
		}
		if len(a) != 1 {
			return true
		}
	}
	return false
}

// resolveAF determines the alternate-allele frequency for a record,
// following upstream process_line's precedence (vcfroh.c:836-890):
// explicit AF-tag, AF-file, AF-dflt, --estimate-AF, then INFO/AC,AN.
//
// After the precedence chain, upstream applies (vcfroh.c:888-890):
//   - if --AF-dflt>0 and (no AF resolved OR AF==0): use --AF-dflt;
//   - else if no AF resolved: skip the site (counted as nno_af);
//   - else if AF==0: skip the site (counted as nno_af).
//
// An AF==0 site with no --AF-dflt is therefore rejected, not clamped.
func resolveAF(v *vcf.Variant, ial int, opts RohOptions, fileAF map[string]float64, afEst estimateAFCohort, afEstActive bool, sampleIdx map[string]int, hardGTErr float64) (float64, bool) {
	var (
		af  float64
		ret bool
	)
	switch {
	case hasNonAFTag(opts.AFTag):
		af, ret = lookupINFOTag(v, opts.AFTag, ial)
	case opts.AFFile != "":
		af, ret = lookupAFFile(v, fileAF)
	case opts.AFDflt != nil && *opts.AFDflt > 0:
		// --AF-dflt is a primary branch (vcfroh.c:850), taken before
		// estimate-AF and AC/AN. The post-resolution clause below
		// still applies it as a fallback for the AF-tag/AF-file paths.
		af, ret = *opts.AFDflt, true
	case afEstActive:
		af, ret = estimateAF(v, ial, afEst, sampleIdx, hardGTErr)
	default:
		af, ret = lookupINFOTag(v, "AF", ial)
		if !ret {
			af, ret = lookupACAN(v)
		}
	}

	// vcfroh.c:888-890 — apply --AF-dflt as a fallback, otherwise an
	// unresolved or zero AF skips the site.
	if opts.AFDflt != nil && *opts.AFDflt > 0 && (!ret || af == 0.0) {
		af = *opts.AFDflt
	} else if !ret {
		return 0, false
	} else if af == 0.0 {
		return 0, false
	}
	// Clamp only the >=1 side to keep emission maths well-defined; an
	// AF of 0 is already rejected above (upstream never clamps it up).
	if af >= 1 {
		af = 1 - 1e-9
	}
	return af, true
}

// hasNonAFTag reports whether the AF tag was set to a non-default
// value. resolveAF treats the bare default "AF" specially so that
// AC/AN can act as the fallback.
func hasNonAFTag(tag string) bool {
	return tag != "" && tag != "AF"
}

func lookupINFOTag(v *vcf.Variant, tag string, ial int) (float64, bool) {
	if v.Info == nil {
		return 0, false
	}
	val, ok := v.Info[tag]
	if !ok || val == "" || val == "." {
		return 0, false
	}
	parts := strings.Split(val, ",")
	// INFO/AF is Number=A: element ial-1 is the freq of the ial-th
	// allele.
	i := ial - 1
	if i < 0 || i >= len(parts) {
		i = 0
	}
	af, err := strconv.ParseFloat(parts[i], 64)
	if err != nil {
		return 0, false
	}
	return af, true
}

func lookupACAN(v *vcf.Variant) (float64, bool) {
	if v.Info == nil {
		return 0, false
	}
	anStr, ok := v.Info["AN"]
	if !ok {
		return 0, false
	}
	an, err := strconv.Atoi(strings.SplitN(anStr, ",", 2)[0])
	if err != nil || an <= 0 {
		return 0, false
	}
	acStr, ok := v.Info["AC"]
	if !ok {
		return 0, false
	}
	ac, err := strconv.Atoi(strings.SplitN(acStr, ",", 2)[0])
	if err != nil || ac < 0 {
		return 0, false
	}
	return float64(ac) / float64(an), true
}

func lookupAFFile(v *vcf.Variant, fileAF map[string]float64) (float64, bool) {
	if len(fileAF) == 0 {
		return 0, false
	}
	af, ok := fileAF[afFileKey(v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))]
	return af, ok
}

func afFileKey(chrom string, pos int, ref, alt string) string {
	return chrom + "\x00" + strconv.Itoa(pos) + "\x00" + ref + "\x00" + alt
}

func loadAFFile(path string) (map[string]float64, error) {
	if path == "" {
		return nil, nil
	}
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("bcftools roh: --AF-file: %w", err)
	}
	defer in.Close()
	out := map[string]float64{}
	sc := bufio.NewScanner(in)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		// CHR \t POS \t REF,ALT \t AF
		fields := strings.Split(line, "\t")
		if len(fields) < 4 {
			continue
		}
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}
		af, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			continue
		}
		// The third column carries the comma-joined alleles; store it
		// verbatim so it matches strings.Join(v.Alt, ...) with REF
		// prepended.
		alleles := strings.SplitN(fields[2], ",", 2)
		if len(alleles) != 2 {
			continue
		}
		out[afFileKey(fields[0], pos, alleles[0], alleles[1])] = af
	}
	return out, sc.Err()
}

// writeRoh emits the requested output sections. The numeric format
// (%.1f for qualities) matches upstream's vcfroh.c.
func writeRoh(out io.Writer, r RohResult, types string) error {
	w := bufio.NewWriter(out)
	wantST := strings.ContainsRune(types, 's')
	wantRG := strings.ContainsRune(types, 'r')
	if wantRG {
		if _, err := w.WriteString("# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality (average fwd-bwd phred score)\n"); err != nil {
			return err
		}
	}
	if wantST {
		if _, err := w.WriteString("# ST\t[2]Sample\t[3]Chromosome\t[4]Position\t[5]State (0:HW, 1:AZ)\t[6]Quality (fwd-bwd phred score)\n"); err != nil {
			return err
		}
	}
	if wantST {
		for _, s := range r.Sites {
			if _, err := fmt.Fprintf(w, "ST\t%s\t%s\t%d\t%d\t%.1f\n",
				s.Sample, s.Chrom, s.Pos, s.State, s.Qual); err != nil {
				return err
			}
		}
	}
	if wantRG {
		for _, reg := range r.Regions {
			if _, err := fmt.Fprintf(w, "RG\t%s\t%s\t%d\t%d\t%d\t%d\t%.1f\n",
				reg.Sample, reg.Chrom, reg.Start, reg.End, reg.Length, reg.NumMarkers, reg.Qual); err != nil {
				return err
			}
		}
	}
	return w.Flush()
}

// parseBufferSize parses upstream's -b/--buffer-size argument. A
// negative first value is a memory budget in MB; the optional second
// value is the overlap. Returns (bufMax, bufOlap, err); both zero
// means unlimited.
func parseBufferSize(spec string, nsamples int) (int, int, error) {
	if spec == "" {
		return 0, 0, nil
	}
	olap := -1
	mainPart := spec
	if comma := strings.IndexByte(spec, ','); comma >= 0 {
		mainPart = spec[:comma]
		o, err := strconv.Atoi(strings.TrimSpace(spec[comma+1:]))
		if err != nil || o < 0 {
			return 0, 0, fmt.Errorf("bcftools roh: could not parse --buffer-size %q", spec)
		}
		olap = o
	}
	val, err := strconv.ParseFloat(strings.TrimSpace(mainPart), 64)
	if err != nil {
		return 0, 0, fmt.Errorf("bcftools roh: could not parse --buffer-size %q", spec)
	}
	var bufMax int
	if val < 0 {
		if nsamples == 0 {
			nsamples = 1
		}
		bufMax = int(math.Abs(val) * 1e6 / (4 + 8*2) / float64(nsamples))
	} else {
		bufMax = int(val)
	}
	if olap < 0 {
		olap = int(float64(bufMax) * 0.01)
	}
	return bufMax, olap, nil
}
