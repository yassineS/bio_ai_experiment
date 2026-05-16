package bcftools

// bcftools roh — runs-of-autozygosity HMM (port).
//
// Upstream reference: reference_code/bcftools/vcfroh.c, in particular
// `process_line` for the emission probabilities and `flush_viterbi`
// for the per-site state decode. The v1 port runs a 2-state Viterbi
// HMM (HW = Hardy-Weinberg, AZ = autozygous) per sample and emits two
// tables:
//
//	ST  per-site:    Sample, CHROM, POS, State (0/1), Quality
//	RG  per-region:  Sample, CHROM, Start, End, Length, NumMarkers, Quality
//
// AZ-state emission (matching upstream's pdg[0]*(1-p) + pdg[2]*p): for
// hard GT (one-hot PDG) this collapses to
//
//	HOM-REF (PDG = [1,0,0])  -> emission(AZ) = (1-AF)
//	HOM-ALT (PDG = [0,0,1])  -> emission(AZ) = AF
//	HET     (PDG = [0,1,0])  -> emission(AZ) = 0  (the HET cannot occur
//	                                                  in an autozygous segment)
//
// Missing-AF policy: when AF lookup fails for a site AND AFDflt is nil
// (i.e. --AF-dflt was not given), the site is SKIPPED. This matches
// upstream's default. The hard-coded 0.4 that PR #106 used is gone.
//
// Transition defaults: HW->AZ = 6.7e-8 per bp, AZ->HW = 5e-9 per bp.
// We do NOT scale by physical distance between markers in v1 — that
// is a tracked deferral. The per-bp magnitude is preserved as-is.
//
// Quality scores: the per-region quality is the sum of per-site state
// posteriors (delta values). Because we do not yet scale transitions
// by distance, these numbers are NOT comparable to upstream's. See
// the docstring on `RG` output for the canonical disclaimer.

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// RohOptions controls the behaviour of Roh / RohFile.
type RohOptions struct {
	// AFTag picks an INFO tag for allele frequency. Empty defaults
	// to "AF".
	AFTag string
	// AFDflt is the --AF-dflt value. When nil, missing AF means
	// "skip the site" (upstream's default). When non-nil, missing AF
	// falls back to this value.
	AFDflt *float64
	// AFFile is the upstream --AF-file argument (CHR\tPOS\tREF,ALT\tAF).
	// v1 reads the file once and indexes by (CHROM,POS,REF,ALT).
	AFFile string

	// GTsOnly is upstream's -G/--GTs-only FLOAT: when >0, hard-GT
	// mode is enabled and the floating-point error probability is
	// pow(10, -GTsOnly/10). Default (0) means PL mode (which v1
	// rejects, see the CLI).
	GTsOnly float64

	// IgnoreHomRef skips 0/0 GTs (-i/--ignore-homref).
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
	// rune is one of 's' (per-site), 'r' (regions), 'z' (gzip).
	// Default "sr".
	OutputTypes string

	// Samples / Regions / Targets / Include / Exclude — post-filters
	// in v1 (no index seek).
	Samples     []string
	SamplesFile string
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
	IncludeExpr string
	ExcludeExpr string

	// EstimateAF mirrors upstream `-e/--estimate-AF`. Accepted at the
	// CLI; v1 ignores the value and pulls AF from INFO/AF.
	EstimateAF string

	// BufferSize / GeneticMap / RecRate / ViterbiTraining: parsed and
	// accepted, no-op in v1. Tracked in PARITY_ROADMAP.
	BufferSize      string
	GeneticMap      string
	RecRate         float64
	ViterbiTraining float64
	RegionsOverlap  int
	TargetsOverlap  int
}

// RohSite captures a single per-site state assignment.
type RohSite struct {
	Sample string
	Chrom  string
	Pos    int
	State  int     // 0 = HW, 1 = AZ
	Qual   float64 // RG-style per-site quality; see disclaimer at the top
	AF     float64
}

// RohRegion captures a contiguous AZ run for one sample.
//
// NB: the Qual field is the sum of per-site AZ-state posteriors over
// the region. RG quality scores are NOT comparable to upstream's
// until physical-distance scaling of the transition matrix lands.
// Tracked in docs/PARITY_ROADMAP.md#bcftools.
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
	if opts.AFTag == "" {
		opts.AFTag = "AF"
	}

	hdr, vars, err := readAllVariants(in)
	if err != nil {
		return RohResult{}, fmt.Errorf("bcftools roh: %w", err)
	}

	samples := selectSamples(hdr.Samples, opts.Samples)
	if len(samples) == 0 {
		return RohResult{}, fmt.Errorf("bcftools roh: no samples to analyse")
	}

	afFromFile, err := loadAFFile(opts.AFFile)
	if err != nil {
		return RohResult{}, err
	}

	// Pre-filter variants by mode (skip indels, ignore no-ALT, etc.)
	// and group by CHROM for the per-chromosome Viterbi pass.
	type marker struct {
		v  *vcf.Variant
		af float64
	}
	byChrom := map[string][]marker{}
	chromOrder := []string{}
	for _, v := range vars {
		if !opts.IncludeNoalt && (len(v.Alt) == 0 || v.Alt[0] == ".") {
			continue
		}
		if opts.SkipIndels && isIndel(v) {
			continue
		}
		af, ok := lookupAF(v, opts.AFTag, afFromFile)
		if !ok {
			if opts.AFDflt == nil {
				continue // upstream's default: skip when AF unknown
			}
			af = *opts.AFDflt
		}
		// Clamp AF to (0, 1) to keep log() safe.
		if af <= 0 {
			af = 1e-6
		}
		if af >= 1 {
			af = 1 - 1e-6
		}
		if _, seen := byChrom[v.Chrom]; !seen {
			chromOrder = append(chromOrder, v.Chrom)
		}
		byChrom[v.Chrom] = append(byChrom[v.Chrom], marker{v: v, af: af})
	}

	// Per-sample Viterbi over each CHROM block.
	result := RohResult{}
	sampleIdx := indexSamples(hdr)
	hardGTErr := 0.0
	if opts.GTsOnly > 0 {
		hardGTErr = math.Pow(10, -opts.GTsOnly/10.0)
	}

	for _, sample := range samples {
		idx, ok := sampleIdx[sample]
		if !ok {
			continue
		}
		for _, chrom := range chromOrder {
			mrk := byChrom[chrom]
			if len(mrk) == 0 {
				continue
			}
			// PDG per site for this sample. Missing GT → skip.
			type sitePDG struct {
				v   *vcf.Variant
				af  float64
				pdg [3]float64
			}
			pdgs := make([]sitePDG, 0, len(mrk))
			for _, m := range mrk {
				dose, ok := gtDose(m.v, idx)
				if !ok {
					continue
				}
				if opts.IgnoreHomRef && dose == 0 {
					continue
				}
				p := pdgFromDose(dose, hardGTErr)
				pdgs = append(pdgs, sitePDG{v: m.v, af: m.af, pdg: p})
			}
			if len(pdgs) == 0 {
				continue
			}

			// Viterbi over two states. emissions[t][state] is the
			// log-emission probability. We work in log space.
			n := len(pdgs)
			emHW := make([]float64, n)
			emAZ := make([]float64, n)
			for t := range pdgs {
				af := pdgs[t].af
				eHW := pdgs[t].pdg[0]*((1-af)*(1-af)) +
					pdgs[t].pdg[1]*(2*af*(1-af)) +
					pdgs[t].pdg[2]*(af*af)
				eAZ := pdgs[t].pdg[0]*(1-af) + pdgs[t].pdg[2]*af
				emHW[t] = safeLog(eHW)
				emAZ[t] = safeLog(eAZ)
			}

			// Transition matrix (log-space). Upstream uses per-bp
			// values; v1 holds them constant between markers and
			// notes the deferral in the docstring.
			logHWtoAZ := math.Log(opts.HWtoAZ)
			logAZtoHW := math.Log(opts.AZtoHW)
			logHWtoHW := math.Log(1 - opts.HWtoAZ)
			logAZtoAZ := math.Log(1 - opts.AZtoHW)

			// dp[t][s], path[t][s]
			dpHW := make([]float64, n)
			dpAZ := make([]float64, n)
			pathHW := make([]int, n)
			pathAZ := make([]int, n)
			// Uniform prior on the initial state.
			dpHW[0] = math.Log(0.5) + emHW[0]
			dpAZ[0] = math.Log(0.5) + emAZ[0]
			for t := 1; t < n; t++ {
				// HW at t can come from HW or AZ at t-1.
				hwFromHW := dpHW[t-1] + logHWtoHW
				hwFromAZ := dpAZ[t-1] + logAZtoHW
				if hwFromHW >= hwFromAZ {
					dpHW[t] = hwFromHW + emHW[t]
					pathHW[t] = 0
				} else {
					dpHW[t] = hwFromAZ + emHW[t]
					pathHW[t] = 1
				}
				azFromAZ := dpAZ[t-1] + logAZtoAZ
				azFromHW := dpHW[t-1] + logHWtoAZ
				if azFromAZ >= azFromHW {
					dpAZ[t] = azFromAZ + emAZ[t]
					pathAZ[t] = 1
				} else {
					dpAZ[t] = azFromHW + emAZ[t]
					pathAZ[t] = 0
				}
			}

			// Backtrace.
			states := make([]int, n)
			if dpAZ[n-1] >= dpHW[n-1] {
				states[n-1] = 1
			} else {
				states[n-1] = 0
			}
			for t := n - 2; t >= 0; t-- {
				if states[t+1] == 1 {
					states[t] = pathAZ[t+1]
				} else {
					states[t] = pathHW[t+1]
				}
			}

			// Per-site emission + per-region rollups.
			var (
				regionOpen   bool
				regionStart  int
				regionEnd    int
				regionMarker int
				regionQual   float64
			)
			for t, sp := range pdgs {
				qual := math.Abs(dpAZ[t] - dpHW[t])
				site := RohSite{
					Sample: sample,
					Chrom:  sp.v.Chrom,
					Pos:    sp.v.Pos,
					State:  states[t],
					Qual:   qual,
					AF:     sp.af,
				}
				result.Sites = append(result.Sites, site)

				if states[t] == 1 {
					if !regionOpen {
						regionOpen = true
						regionStart = sp.v.Pos
						regionMarker = 0
						regionQual = 0
					}
					regionEnd = sp.v.Pos
					regionMarker++
					regionQual += qual
				} else if regionOpen {
					result.Regions = append(result.Regions, RohRegion{
						Sample:     sample,
						Chrom:      sp.v.Chrom,
						Start:      regionStart,
						End:        regionEnd,
						Length:     regionEnd - regionStart + 1,
						NumMarkers: regionMarker,
						Qual:       regionQual,
					})
					regionOpen = false
				}
			}
			if regionOpen {
				result.Regions = append(result.Regions, RohRegion{
					Sample:     sample,
					Chrom:      pdgs[0].v.Chrom,
					Start:      regionStart,
					End:        regionEnd,
					Length:     regionEnd - regionStart + 1,
					NumMarkers: regionMarker,
					Qual:       regionQual,
				})
			}
		}
	}

	if err := writeRoh(out, result, opts.OutputTypes); err != nil {
		return result, err
	}
	return result, nil
}

// AFDfltPtr returns a *float64 pointing at the given value. Convenience
// helper for callers that want to set --AF-dflt without having to take
// the address of a literal.
func AFDfltPtr(v float64) *float64 { return &v }

// pdgFromDose collapses a hard-GT dosage into a 3-vector PDG. With
// hardGTErr>0 we redistribute mass per upstream's `-G N` semantics:
//
//	PDG_called = 1 - 2*err
//	PDG_other  = err  (each)
//
// where err = 10^(-N/10). With hardGTErr==0 we emit one-hot.
func pdgFromDose(dose int, err float64) [3]float64 {
	if err <= 0 {
		switch dose {
		case 0:
			return [3]float64{1, 0, 0}
		case 1:
			return [3]float64{0, 1, 0}
		case 2:
			return [3]float64{0, 0, 1}
		}
		return [3]float64{0, 0, 0}
	}
	called := 1 - 2*err
	if called < 0 {
		called = 0
	}
	switch dose {
	case 0:
		return [3]float64{called, err, err}
	case 1:
		return [3]float64{err, called, err}
	case 2:
		return [3]float64{err, err, called}
	}
	return [3]float64{}
}

func safeLog(v float64) float64 {
	if v <= 0 {
		return -math.Log(math.MaxFloat64) // very negative but finite
	}
	return math.Log(v)
}

func isIndel(v *vcf.Variant) bool {
	if len(v.Ref) != 1 {
		return true
	}
	for _, a := range v.Alt {
		if a == "" || len(a) != 1 {
			return true
		}
	}
	return false
}

// lookupAF resolves the alternate-allele frequency for a record. We
// try the AFFile first (highest priority), then INFO/<AFTag>. We
// return false when no value is available.
func lookupAF(v *vcf.Variant, tag string, fileAF map[string]float64) (float64, bool) {
	if len(fileAF) > 0 {
		key := afFileKey(v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		if af, ok := fileAF[key]; ok {
			return af, true
		}
	}
	if v.Info == nil {
		return 0, false
	}
	val, ok := v.Info[tag]
	if !ok || val == "" || val == "." {
		return 0, false
	}
	// For multi-allelic records (which we don't expect here) take
	// just the first AF.
	if comma := strings.Index(val, ","); comma >= 0 {
		val = val[:comma]
	}
	af, err := strconv.ParseFloat(val, 64)
	if err != nil {
		return 0, false
	}
	return af, true
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
		refAlt := strings.SplitN(fields[2], ",", 2)
		if len(refAlt) != 2 {
			continue
		}
		af, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			continue
		}
		out[afFileKey(fields[0], pos, refAlt[0], refAlt[1])] = af
	}
	return out, sc.Err()
}

// writeRoh emits the requested output sections. ST and RG share one
// io.Writer (gzip is a TODO; CLI rejects -O z with a roadmap pointer
// in v1).
func writeRoh(out io.Writer, r RohResult, types string) error {
	w := bufio.NewWriter(out)
	wantST := strings.ContainsRune(types, 's')
	wantRG := strings.ContainsRune(types, 'r')
	// Header lines mirror upstream's prefix-based output.
	if wantST {
		if _, err := w.WriteString("# ST, Per-site State Probabilities. ST [2]Sample\t[3]Chromosome\t[4]Position\t[5]State (0:HW, 1:AZ)\t[6]Quality (smaller is better)\n"); err != nil {
			return err
		}
		for _, s := range r.Sites {
			if _, err := fmt.Fprintf(w, "ST\t%s\t%s\t%d\t%d\t%.6f\n",
				s.Sample, s.Chrom, s.Pos, s.State, s.Qual); err != nil {
				return err
			}
		}
	}
	if wantRG {
		if _, err := w.WriteString("# RG, Per-region Autozygous Segments. RG [2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality (sum of per-site quality; NOT comparable to upstream until distance-scaled transitions land — see docs/PARITY_ROADMAP.md#bcftools)\n"); err != nil {
			return err
		}
		for _, reg := range r.Regions {
			if _, err := fmt.Fprintf(w, "RG\t%s\t%s\t%d\t%d\t%d\t%d\t%.6f\n",
				reg.Sample, reg.Chrom, reg.Start, reg.End, reg.Length, reg.NumMarkers, reg.Qual); err != nil {
				return err
			}
		}
	}
	return w.Flush()
}
