package bcftools

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

// CallModel selects the variant-calling algorithm.
type CallModel int

const (
	// CallModelNone means the caller was not requested. Call() rejects it.
	CallModelNone CallModel = iota
	// CallModelConsensus selects the original Li-2011 consensus caller (`-c`).
	CallModelConsensus
	// CallModelMultiallelic selects the multiallelic caller (`-m`). When
	// the mpileup INFO/QS annotation is present, Call runs the faithful
	// port of mcall.c (see callm.go): EM allele-frequency estimation, the
	// per-site QUAL, the max-likelihood GT, and the INFO rewrite
	// (AN/AC/DP4/MQ). The mpileup | call -m pipeline byte-matches upstream.
	// Synthetic PL-only fixtures (no QS) fall through to the heuristic
	// path. Remaining gaps (trios, --constrain, sample groups) are tracked
	// in docs/PARITY_ROADMAP.md (bcftools call section).
	CallModelMultiallelic
)

// PloidySpec selects the per-sample ploidy. Fixed ploidies (1 or 2 across
// every sample/contig) live here directly; the GRCh37/GRCh38 per-contig
// sex-chromosome maps come through CallOptions.PloidyTable (see
// call_ploidy.go).
type PloidySpec int

const (
	// PloidyDiploid treats every sample as diploid (default).
	PloidyDiploid PloidySpec = 2
	// PloidyHaploid treats every sample as haploid (same as upstream `-X`).
	PloidyHaploid PloidySpec = 1
)

// CallOptions controls the behaviour of Call / CallFile.
type CallOptions struct {
	// Model selects the variant caller. Required.
	Model CallModel
	// KeepAlts mirrors upstream `-A`: emit every ALT allele declared in
	// the input, even those with zero supporting reads.
	KeepAlts bool
	// VariantsOnly mirrors upstream `-v`: drop records whose called
	// genotypes are all reference.
	VariantsOnly bool
	// Prior is the per-base mutation rate prior (upstream default 1.1e-3).
	Prior float64
	// PvalThreshold is the cutoff for the variant posterior: a site is
	// emitted as variant when posterior > 1 - PvalThreshold. The upstream
	// default is 0.5, matching `-p 0.5`.
	PvalThreshold float64
	// Ploidy selects diploid (default) or haploid calling. When
	// PloidyTable is set it overrides this for the per-record, per-sample
	// resolution; Ploidy is then used only as the global fallback.
	Ploidy PloidySpec
	// PloidySpec, when non-empty, is the textual ploidy spec ("2", "1",
	// "GRCh37", ...). Setting it to one of the GRCh* aliases populates
	// PloidyTable automatically when ParsePloidySpec is used.
	PloidySpec string
	// PloidyTable, when non-nil, replaces the global Ploidy with a
	// per-region, per-sex map (built from --ploidy GRCh37/GRCh38 or a
	// --ploidy-file argument).
	PloidyTable *PloidyTable
	// SampleSexes maps the input sample index to the sex id registered
	// in PloidyTable. When nil every sample defaults to
	// PloidyTable.DefaultSexID (the last registered sex, F for the
	// GRCh predefs — matching vcfcall.c sample2sex initialisation).
	SampleSexes []int
	// OutputFormat passes through to the writer (see openOutput).
	OutputFormat OutputFormat
	// CompressLevel sets the gzip level for -O z output.
	CompressLevel int
	// Regions / Targets / Samples mirror the same fields on ViewOptions
	// and are applied with the same semantics (Regions are index-aware,
	// Targets are post-filter).
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
	Samples     []string
	SamplesFile string
	// GVCFSpec is the raw textual "--gvcf 0,5,10" value preserved for
	// validation reporting. When non-empty, Call parses it into
	// GVCFRange before streaming starts (see callm_gvcf.go).
	GVCFSpec string
	// GVCFRange is the parsed-and-sorted DP-threshold slice (upstream's
	// _gvcf_t.dp_range). When non-empty, consecutive post-mcall REF-only
	// records that share a per-sample MIN_DP bin are banded into a
	// single INFO/END+MIN_DP record by callGVCFBlocker.
	GVCFRange []int
}

// defaults applies upstream-equivalent defaults for any unset field.
func (o *CallOptions) defaults() {
	if o.Prior == 0 {
		o.Prior = 1.1e-3
	}
	if o.PvalThreshold == 0 {
		o.PvalThreshold = 0.5
	}
	if o.Ploidy == 0 {
		o.Ploidy = PloidyDiploid
	}
}

// ParsePloidySpec turns a "--ploidy" string into the typed value. The
// returned spec string is preserved so Call() can record provenance.
func ParsePloidySpec(s string) (PloidySpec, string, error) {
	switch strings.TrimSpace(s) {
	case "", "2":
		return PloidyDiploid, "2", nil
	case "1":
		return PloidyHaploid, "1", nil
	case "GRCh37", "GRCh38":
		return PloidyDiploid, s, nil
	}
	return 0, "", fmt.Errorf("bcftools call: unknown --ploidy %q (expect 1, 2, GRCh37, GRCh38)", s)
}

// BuildPloidyTableFromSpec returns a PloidyTable for one of the
// recognised --ploidy aliases ("GRCh37", "GRCh38", "1", "2"), or nil
// when the spec resolves to a simple uniform ploidy that doesn't need
// the table machinery.
func BuildPloidyTableFromSpec(spec string) (*PloidyTable, error) {
	body := LookupPredefPloidy(spec)
	if body == "" {
		return nil, nil
	}
	return ParsePloidyTable(body, 2)
}

// resolveSampleSexes returns a sex-id slice of length nsmpl. When
// opts.SampleSexes is set we use it (clamped to the registered range);
// otherwise every sample defaults to PloidyTable.DefaultSexID(), which
// matches vcfcall.c's `sample2sex[i] = args->nsex - 1` initialisation.
func resolveSampleSexes(tbl *PloidyTable, samples []string, opts *CallOptions) []int {
	if tbl == nil {
		return nil
	}
	dflt := tbl.DefaultSexID()
	out := make([]int, len(samples))
	for i := range out {
		out[i] = dflt
	}
	if opts != nil {
		for i, sid := range opts.SampleSexes {
			if i >= len(out) {
				break
			}
			if sid >= 0 && sid < tbl.NSex() {
				out[i] = sid
			}
		}
	}
	return out
}

// perSamplePloidy returns the per-sample ploidy slice for one record.
// When opts.PloidyTable is nil it falls back to opts.Ploidy (the global
// 1- or 2-uniform mode), so existing call sites keep working. nsmpl is
// the number of input samples on the record.
func perSamplePloidy(opts CallOptions, sexes []int, chrom string, pos, nsmpl int) []int {
	if opts.PloidyTable == nil {
		out := make([]int, nsmpl)
		p := int(opts.Ploidy)
		if p == 0 {
			p = 2
		}
		for i := range out {
			out[i] = p
		}
		return out
	}
	if len(sexes) < nsmpl {
		// Pad with the table default (matches vcfcall.c init).
		dflt := opts.PloidyTable.DefaultSexID()
		ext := make([]int, nsmpl)
		copy(ext, sexes)
		for i := len(sexes); i < nsmpl; i++ {
			ext[i] = dflt
		}
		sexes = ext
	}
	return opts.PloidyTable.PerSamplePloidy(chrom, pos, sexes[:nsmpl])
}

// Call streams VCF/BCF input from in, applies the consensus / multiallelic
// variant caller, and writes the called records to out. It is the
// streaming entry point used by `bcftools call` when no region query is
// requested.
func Call(in io.Reader, out io.Writer, opts CallOptions) (int, error) {
	opts.defaults()
	if opts.Model == CallModelNone {
		return 0, fmt.Errorf("bcftools call: a caller must be selected (-c / --consensus-caller or -m / --multiallelic-caller)")
	}
	// Resolve --gvcf spec (library callers may leave the parsed slice
	// nil and pass the textual form only).
	if opts.GVCFSpec != "" && len(opts.GVCFRange) == 0 {
		ranges, err := parseGVCFRanges(opts.GVCFSpec)
		if err != nil {
			return 0, err
		}
		opts.GVCFRange = ranges
	}
	// Upstream vcfcall.c:1190 rejects --variants-only with --gvcf.
	if len(opts.GVCFRange) > 0 && opts.VariantsOnly {
		return 0, fmt.Errorf("bcftools call: The two options cannot be combined: --variants-only and --gvcf")
	}
	// Upstream vcfcall.c:1182 rejects --gvcf with the consensus caller.
	if len(opts.GVCFRange) > 0 && opts.Model == CallModelConsensus {
		return 0, fmt.Errorf("bcftools call: gvcf -g option not functional with -c calling mode yet")
	}
	if opts.PloidyTable == nil && (opts.PloidySpec == "GRCh37" || opts.PloidySpec == "GRCh38") {
		tbl, err := BuildPloidyTableFromSpec(opts.PloidySpec)
		if err != nil {
			return 0, err
		}
		opts.PloidyTable = tbl
	}
	parsedTargets, err := parseRegions(opts.Targets)
	if err != nil {
		return 0, err
	}
	postFilters := parsedTargets
	if len(opts.Regions) > 0 {
		regs, err := parseRegions(opts.Regions)
		if err != nil {
			return 0, err
		}
		postFilters = append(postFilters, regs...)
	}
	return callStreaming(in, out, opts, postFilters)
}

// CallFile is the file-aware entry point for `bcftools call`. Today it
// always streams (no chunk-seek), but it does normalise the input path so
// gzipped / bgzipped / plain files all work. Region queries are evaluated
// as post-filters in v1 (matching the streaming path in `view`).
func CallFile(path string, out io.Writer, opts CallOptions, stderr io.Writer) (int, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	defer in.Close()
	if len(opts.Regions) > 0 && stderr != nil {
		fmt.Fprintln(stderr, "bcftools call: index-backed region queries are deferred; treating -r as a post-filter")
	}
	return Call(in, out, opts)
}

// callStreaming is the inner loop. It consumes a VCF stream, dispatches
// each record through the chosen caller, and writes the results. The
// targets slice is the union of -t and (when no index is available) -r.
func callStreaming(in io.Reader, out io.Writer, opts CallOptions, targets []region) (int, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return 0, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return 0, fmt.Errorf("bcftools call: BCF input is not yet wired through the caller; convert with `bcftools view in.bcf` first (see docs/PARITY_ROADMAP.md bcftools call)")
	}
	r := vcf.NewReader(br)
	hdr, err := r.ReadHeader()
	if err != nil {
		return 0, err
	}
	hdr = filterHeaderSamples(hdr, opts.Samples)
	// For -c, when the input declares INFO/I16 we route through the
	// faithful consensus port (callc.go) and rewrite the header
	// accordingly. Otherwise the heuristic v1 augmentation applies
	// (FORMAT/GT + AC/AN only).
	if opts.Model == CallModelConsensus && headerHasInfo(hdr, "I16") {
		hdr = augmentCallHeaderConsensus(hdr)
	} else {
		hdr = augmentCallHeader(hdr, opts.Model)
	}

	sexes := resolveSampleSexes(opts.PloidyTable, hdr.Samples, &opts)

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	count := 0
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, err
		}
		if len(targets) > 0 && !overlapsAny(v, targets) {
			continue
		}
		if len(opts.Samples) > 0 {
			restrictSamples(v, opts.Samples)
		}
		samplePloidy := perSamplePloidy(opts, sexes, v.Chrom, v.Pos, len(v.Samples))
		called, keep := callVariant(v, opts, samplePloidy)
		if !keep {
			continue
		}
		if err := w.Write(called); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// augmentCallHeader inserts the meta lines that upstream bcftools adds
// when calling: ##INFO/AC, ##INFO/AN, and the standard FORMAT/GT
// declaration (idempotent when already present).
func augmentCallHeader(hdr *vcf.Header, model CallModel) *vcf.Header {
	if hdr == nil {
		return hdr
	}
	out := &vcf.Header{Samples: append([]string(nil), hdr.Samples...)}
	// The multiallelic caller rewrites the INFO block: the auxiliary
	// mpileup tags I16/QS are removed (they are consumed into DP4/MQ/QUAL)
	// and FORMAT/GT plus INFO AC/AN/DP4/MQ are appended in upstream order.
	dropMpileupAux := model == CallModelMultiallelic
	for _, m := range hdr.MetaInfo {
		if dropMpileupAux && (strings.HasPrefix(m, `##INFO=<ID=I16,`) || strings.HasPrefix(m, `##INFO=<ID=QS,`)) {
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
	}
	declarations := []struct {
		marker string
		line   string
	}{
		{`##FORMAT=<ID=GT,`, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`},
		{`##INFO=<ID=AC,`, `##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">`},
		{`##INFO=<ID=AN,`, `##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">`},
		{`##INFO=<ID=DP4,`, `##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">`},
		{`##INFO=<ID=MQ,`, `##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">`},
	}
	for _, d := range declarations {
		// DP4/MQ are only emitted by the multiallelic path.
		if !dropMpileupAux && (strings.HasPrefix(d.marker, `##INFO=<ID=DP4,`) || strings.HasPrefix(d.marker, `##INFO=<ID=MQ,`)) {
			continue
		}
		found := false
		for _, m := range out.MetaInfo {
			if strings.HasPrefix(m, d.marker) {
				found = true
				break
			}
		}
		if !found {
			out.MetaInfo = append(out.MetaInfo, d.line)
		}
	}
	return out
}

// callVariant runs the consensus or multiallelic caller on v and returns
// the called record plus a "keep" flag. v is consumed (its Samples slice
// is rewritten in place with the called GTs); callers should pass a
// per-record value rather than reusing a shared variant.
//
// The decision logic:
//
//   - For each sample we find the most-likely genotype index from PL.
//   - We compute a variant posterior using a Hardy-Weinberg + mutation
//     rate prior (the Li 2011 model — same family upstream's `-c` uses).
//   - The site is "variant" iff posterior > 1 - opts.PvalThreshold OR any
//     called genotype is non-reference. (The latter is a fail-safe that
//     mirrors upstream's behaviour on small-sample inputs where the
//     prior dominates the posterior.)
//   - The site is emitted iff !opts.VariantsOnly OR the site is variant
//     OR opts.KeepAlts is set.
func callVariant(v *vcf.Variant, opts CallOptions, samplePloidy []int) (*vcf.Variant, bool) {
	// The faithful multiallelic caller runs when -m is selected and the
	// mpileup INFO/QS annotation is present. The synthetic PL-only
	// fixtures (no QS) fall through to the heuristic path below.
	if opts.Model == CallModelMultiallelic && hasQS(v) {
		if out, keep, ok := mcallSite(v, opts, samplePloidy); ok {
			return out, keep
		}
	}
	// The faithful consensus caller runs when -c is selected and the
	// mpileup INFO/I16 annotation is present. The synthetic PL-only
	// fixtures (no I16) fall through to the heuristic v1 path below.
	if opts.Model == CallModelConsensus && hasI16(v) {
		if out, keep, ok := ccallSite(v, opts, samplePloidy); ok {
			return out, keep
		}
	}
	nAlts := len(v.Alt)
	if nAlts == 1 && v.Alt[0] == "." {
		nAlts = 0
		v.Alt = nil
	}
	nAlleles := nAlts + 1
	// Effective per-sample ploidy. When --ploidy is the global "1" or
	// "2" this matches opts.Ploidy for every sample; for the GRCh*
	// tables the heuristic path still runs only when all samples on a
	// record agree (mixed-sex PED is plumbed through the mcall path).
	effPloidy := func(i int) PloidySpec {
		if i < len(samplePloidy) {
			if samplePloidy[i] == 1 {
				return PloidyHaploid
			}
			if samplePloidy[i] == 0 {
				return 0
			}
		}
		return PloidyDiploid
	}
	plByGT := make([][]int, len(v.Samples))
	mostLikely := make([]int, len(v.Samples))
	gtPloidy := make([]PloidySpec, len(v.Samples))
	haveGenotypeData := false
	for i, s := range v.Samples {
		p := effPloidy(i)
		gtPloidy[i] = p
		if p == 0 {
			plByGT[i] = nil
			mostLikely[i] = -1
			continue
		}
		pl, ok := decodePL(s.Data["PL"], nAlleles, p)
		if !ok {
			plByGT[i] = nil
			mostLikely[i] = -1
			continue
		}
		haveGenotypeData = true
		plByGT[i] = pl
		mostLikely[i] = argMinIndex(pl)
	}

	ac := make([]int, nAlts)
	an := 0
	for i := range v.Samples {
		if mostLikely[i] < 0 {
			continue
		}
		a1, a2, ok := decomposeGTIndex(mostLikely[i], nAlleles, gtPloidy[i])
		if !ok {
			continue
		}
		if a1 > 0 && a1-1 < nAlts {
			ac[a1-1]++
		}
		an++
		if gtPloidy[i] == PloidyDiploid {
			if a2 > 0 && a2-1 < nAlts {
				ac[a2-1]++
			}
			an++
		}
	}

	logRatio := 0.0
	for _, pl := range plByGT {
		if pl == nil {
			continue
		}
		refPL := pl[0]
		bestNonRef := math.MaxInt32
		for j := 1; j < len(pl); j++ {
			if pl[j] < bestNonRef {
				bestNonRef = pl[j]
			}
		}
		if bestNonRef == math.MaxInt32 {
			continue
		}
		ratio := float64(refPL-bestNonRef) / 10.0
		if ratio > 20 {
			ratio = 20
		}
		if ratio < -20 {
			ratio = -20
		}
		logRatio += ratio
	}
	priorLog10 := math.Log10(opts.Prior)
	posteriorLog10 := logRatio + priorLog10
	var posterior float64
	if posteriorLog10 > 50 {
		posterior = 1
	} else if posteriorLog10 < -50 {
		posterior = math.Pow(10, posteriorLog10)
	} else {
		x := math.Pow(10, posteriorLog10)
		posterior = x / (1 + x)
	}

	anyAltCall := false
	for _, c := range ac {
		if c > 0 {
			anyAltCall = true
			break
		}
	}
	threshold := 1 - opts.PvalThreshold
	isVariant := haveGenotypeData && (posterior > threshold || anyAltCall)

	if opts.VariantsOnly && !isVariant && !opts.KeepAlts {
		return v, false
	}

	out := *v
	out.Samples = make([]vcf.Sample, len(v.Samples))
	for i, s := range v.Samples {
		newSample := vcf.Sample{Name: s.Name, Data: copyStringMap(s.Data)}
		if mostLikely[i] >= 0 {
			newSample.Data["GT"] = encodeGT(mostLikely[i], nAlleles, gtPloidy[i])
		} else if gtPloidy[i] == 0 {
			// Samples with ploidy 0 (e.g. F on chrY) are emitted as
			// missing, matching upstream mcall.c.
			newSample.Data["GT"] = "."
		}
		out.Samples[i] = newSample
	}
	if !hasFormat(out.Format, "GT") {
		out.Format = append([]string{"GT"}, out.Format...)
	}

	if !opts.KeepAlts && nAlts > 0 {
		out.Alt = trimUnsupportedAlts(out.Alt, ac, &out)
	}
	// VCF spec requires ALT="." when no ALT alleles remain.
	if len(out.Alt) == 0 {
		out.Alt = []string{"."}
	}

	if isVariant {
		if posterior >= 1-1e-30 {
			out.Qual = 999
		} else {
			out.Qual = -10 * math.Log10(1-posterior)
			if out.Qual > 999 {
				out.Qual = 999
			}
			if out.Qual < 0 {
				out.Qual = 0
			}
		}
	} else {
		// QUAL is left as "missing" (".") for non-variant sites; upstream
		// emits "." rather than 0 in this case.
		out.Qual = -1
	}

	out.Info = copyStringMap(out.Info)
	out.InfoOrder = append([]string(nil), out.InfoOrder...)
	// Emit AC / AN only when at least one real ALT allele remains.
	realAlts := len(out.Alt)
	if realAlts == 1 && out.Alt[0] == "." {
		realAlts = 0
	}
	if realAlts > 0 {
		acStrs := make([]string, realAlts)
		newAC := computeACFromGT(&out, realAlts)
		for i := 0; i < realAlts; i++ {
			if i < len(newAC) {
				acStrs[i] = strconv.Itoa(newAC[i])
			} else {
				acStrs[i] = "0"
			}
		}
		setInfo(&out, "AC", strings.Join(acStrs, ","))
		setInfo(&out, "AN", strconv.Itoa(totalAN(&out)))
	}
	_ = an

	return &out, true
}

// hasFormat reports whether the FORMAT slice already declares key.
func hasFormat(fmtKeys []string, key string) bool {
	for _, k := range fmtKeys {
		if k == key {
			return true
		}
	}
	return false
}

// copyStringMap returns a shallow copy of m so callers can mutate it
// without disturbing the source.
func copyStringMap(m map[string]string) map[string]string {
	if m == nil {
		return map[string]string{}
	}
	out := make(map[string]string, len(m))
	for k, v := range m {
		out[k] = v
	}
	return out
}

// setInfo writes (or overwrites) an INFO field while preserving the
// InfoOrder slice. New keys are appended; existing keys keep their
// original position.
func setInfo(v *vcf.Variant, key, value string) {
	if v.Info == nil {
		v.Info = make(map[string]string)
	}
	if _, ok := v.Info[key]; !ok {
		v.InfoOrder = append(v.InfoOrder, key)
	}
	v.Info[key] = value
}

// computeACFromGT recounts INFO/AC across the called genotypes. It
// returns a slice of length max(numAlts, observed-max-allele-index).
func computeACFromGT(v *vcf.Variant, numAlts int) []int {
	out := make([]int, numAlts)
	for _, s := range v.Samples {
		gt := s.Data["GT"]
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			if a == "." || a == "" {
				continue
			}
			n, err := strconv.Atoi(a)
			if err != nil || n <= 0 {
				continue
			}
			idx := n - 1
			if idx < len(out) {
				out[idx]++
			}
		}
	}
	return out
}

// totalAN returns the count of called alleles across every sample.
func totalAN(v *vcf.Variant) int {
	total := 0
	for _, s := range v.Samples {
		gt := s.Data["GT"]
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			if a == "." || a == "" {
				continue
			}
			if _, err := strconv.Atoi(a); err == nil {
				total++
			}
		}
	}
	return total
}

// trimUnsupportedAlts drops ALT alleles whose ac[k] is zero. It also
// rewrites the per-sample GT calls to renumber surviving alleles. This
// runs only when -A is not set.
func trimUnsupportedAlts(alts []string, ac []int, v *vcf.Variant) []string {
	if len(alts) == 0 || len(ac) == 0 {
		return alts
	}
	keepIdx := make([]int, 0, len(alts))
	for i := 0; i < len(alts) && i < len(ac); i++ {
		if ac[i] > 0 {
			keepIdx = append(keepIdx, i)
		}
	}
	if len(keepIdx) == len(alts) {
		return alts
	}
	remap := make(map[int]int, len(alts)+1)
	remap[0] = 0
	newAlts := make([]string, 0, len(keepIdx))
	for newI, oldI := range keepIdx {
		remap[oldI+1] = newI + 1
		newAlts = append(newAlts, alts[oldI])
	}
	for i := range v.Samples {
		gt := v.Samples[i].Data["GT"]
		v.Samples[i].Data["GT"] = remapGTByIndex(gt, remap)
	}
	return newAlts
}

// remapGTByIndex renumbers allele indices in a GT string per the remap
// map. Unknown indices become ".". This is the index-aware sibling of
// the per-ALT remapGT used by `bcftools norm`.
func remapGTByIndex(gt string, remap map[int]int) string {
	if gt == "" || gt == "." {
		return gt
	}
	sep := byte('/')
	if strings.Contains(gt, "|") {
		sep = '|'
	}
	parts := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	out := make([]string, 0, len(parts))
	for _, p := range parts {
		if p == "." || p == "" {
			out = append(out, ".")
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			out = append(out, ".")
			continue
		}
		newN, ok := remap[n]
		if !ok {
			out = append(out, ".")
			continue
		}
		out = append(out, strconv.Itoa(newN))
	}
	return strings.Join(out, string(sep))
}

// argMinIndex returns the index of the smallest element in xs.
func argMinIndex(xs []int) int {
	best := 0
	for i := 1; i < len(xs); i++ {
		if xs[i] < xs[best] {
			best = i
		}
	}
	return best
}

// decodePL parses a colon-free PL string ("0,30,255") into a slice of
// Phred-scaled likelihoods, one per possible genotype. The expected
// length follows the bcftools convention:
//
//	diploid:  (n*(n+1))/2 entries for n alleles
//	haploid:  n entries
//
// Missing values ("." or empty) are decoded as 255 to mark them as
// "extremely unlikely" without breaking downstream maths. The bool is
// false when the field is missing entirely.
func decodePL(s string, nAlleles int, ploidy PloidySpec) ([]int, bool) {
	if s == "" || s == "." {
		return nil, false
	}
	parts := strings.Split(s, ",")
	expected := expectedPLLen(nAlleles, ploidy)
	out := make([]int, len(parts))
	for i, p := range parts {
		if p == "." || p == "" {
			out[i] = 255
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, false
		}
		out[i] = n
	}
	for len(out) < expected {
		out = append(out, 255)
	}
	return out, true
}

// expectedPLLen returns the canonical PL vector length for a site with
// nAlleles total alleles (1 REF + (nAlleles-1) ALTs) at the given ploidy.
func expectedPLLen(nAlleles int, ploidy PloidySpec) int {
	switch ploidy {
	case PloidyHaploid:
		return nAlleles
	default:
		return nAlleles * (nAlleles + 1) / 2
	}
}

// decomposeGTIndex turns a PL-vector index back into the (a1, a2)
// genotype it represents. For diploid sites bcftools uses the canonical
// VCF ordering:
//
//	idx 0 -> 0/0   idx 1 -> 0/1   idx 2 -> 1/1
//	idx 3 -> 0/2   idx 4 -> 1/2   idx 5 -> 2/2 ...
//
// For haploid sites the index is the allele number directly.
func decomposeGTIndex(idx, nAlleles int, ploidy PloidySpec) (int, int, bool) {
	if ploidy == PloidyHaploid {
		if idx < 0 || idx >= nAlleles {
			return 0, 0, false
		}
		return idx, 0, true
	}
	k := 0
	for a2 := 0; a2 < nAlleles; a2++ {
		for a1 := 0; a1 <= a2; a1++ {
			if k == idx {
				return a1, a2, true
			}
			k++
		}
	}
	return 0, 0, false
}

// encodeGT renders the (a1, a2) genotype that corresponds to PL index
// idx back into a "a/b" or "a" string.
func encodeGT(idx, nAlleles int, ploidy PloidySpec) string {
	a1, a2, ok := decomposeGTIndex(idx, nAlleles, ploidy)
	if !ok {
		return "."
	}
	if ploidy == PloidyHaploid {
		return strconv.Itoa(a1)
	}
	return fmt.Sprintf("%d/%d", a1, a2)
}
