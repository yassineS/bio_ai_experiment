package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// HapConvertOptions controls the IMPUTE2 HAP/legend/sample conversion
// modes of `bcftools convert` (the --hapsample / --haplegendsample
// family and their *2vcf inverses).
//
// The port mirrors upstream's vcfconvert.c: vcf_to_hapsample,
// vcf_to_haplegendsample, hapsample_to_vcf and haplegendsample_to_vcf.
// Genotype encoding follows convert.c's process_gt_to_hap (and
// process_gt_to_hap2 when Hap2Dip is set).
type HapConvertOptions struct {
	// Prefix is the value passed to --hapsample / --haplegendsample
	// (output) or --hapsample2vcf / --haplegendsample2vcf (input). It
	// may be a single prefix or a comma-separated list of explicit
	// filenames (2 entries for hap/sample, 3 for hap/legend/sample).
	Prefix string

	// Hap2Dip mirrors --haploid2diploid: haploid genotypes are emitted
	// as diploid homozygotes in the .hap output.
	Hap2Dip bool

	// OutputFormat / CompressLevel control the VCF/BCF container emitted
	// by the *2vcf inverse modes.
	OutputFormat  OutputFormat
	CompressLevel int

	// Samples / SamplesFile / ForceSamples restrict the per-sample
	// columns of the input VCF before export (VCF->hap directions only).
	Samples      []string
	SamplesFile  string
	ForceSamples bool

	// Regions / RegionsFile / Targets / TargetsFile apply CHROM[:beg-end]
	// post-filters to the input VCF (VCF->hap directions only).
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string

	// IncludeExpr / ExcludeExpr are the standard -i / -e filter
	// expressions, applied to the input VCF (VCF->hap directions only).
	IncludeExpr string
	ExcludeExpr string
}

// hapOutputNames resolves a --hapsample prefix (or explicit
// comma-separated name list) into the hap and sample file names. A bare
// prefix expands to "<prefix>.hap.gz" and "<prefix>.samples"; a
// two-element list supplies the names directly (with "." or "" meaning
// "skip this file").
func hapOutputNames(prefix string) (hapName, sampleName string, err error) {
	files := strings.Split(prefix, ",")
	switch len(files) {
	case 1:
		return prefix + ".hap.gz", prefix + ".samples", nil
	case 2:
		return resolveName(files[0]), resolveName(files[1]), nil
	default:
		return "", "", fmt.Errorf("error parsing --hapsample filenames: %s", prefix)
	}
}

// hapLegendOutputNames resolves a --haplegendsample prefix (or explicit
// comma-separated name list) into the hap, legend and sample file names.
// A bare prefix expands to "<prefix>.hap.gz", "<prefix>.legend.gz" and
// "<prefix>.samples"; a three-element list supplies the names directly
// (with "." or "" meaning "skip this file").
func hapLegendOutputNames(prefix string) (hapName, legendName, sampleName string, err error) {
	files := strings.Split(prefix, ",")
	switch len(files) {
	case 1:
		return prefix + ".hap.gz", prefix + ".legend.gz", prefix + ".samples", nil
	case 3:
		return resolveName(files[0]), resolveName(files[1]), resolveName(files[2]), nil
	default:
		return "", "", "", fmt.Errorf("error parsing --haplegendsample filenames: %s", prefix)
	}
}

// resolveName maps the upstream "skip" sentinels ("" and ".") to an empty
// string and otherwise returns the name unchanged.
func resolveName(name string) string {
	if name == "" || name == "." {
		return ""
	}
	return name
}

// hapInputNames resolves a --hapsample2vcf argument into the hap and
// sample file names. A bare prefix expands to "<prefix>.hap.gz" and
// "<prefix>.samples"; a comma-separated "hap,sample" supplies them
// directly.
func hapInputNames(arg string) (hapName, sampleName string) {
	if i := strings.IndexByte(arg, ','); i >= 0 {
		return arg[:i], arg[i+1:]
	}
	return arg + ".hap.gz", arg + ".samples"
}

// hapLegendInputNames resolves a --haplegendsample2vcf argument into the
// hap, legend and sample file names. A bare prefix expands to the default
// suffixes; a comma-separated "hap,legend,sample" supplies them directly.
func hapLegendInputNames(arg string) (hapName, legendName, sampleName string, err error) {
	i := strings.IndexByte(arg, ',')
	if i < 0 {
		return arg + ".hap.gz", arg + ".legend.gz", arg + ".samples", nil
	}
	rest := arg[i+1:]
	j := strings.IndexByte(rest, ',')
	if j < 0 {
		return "", "", "", fmt.Errorf("could not parse hap/legend/sample file names: %s", arg)
	}
	// Upstream order on the explicit path is hap,legend,sample.
	return arg[:i], rest[:j], rest[j+1:], nil
}

// isGzName reports whether name ends in ".gz" (case-insensitive), the
// signal upstream uses to decide whether to BGZF-compress an output file.
func isGzName(name string) bool {
	return len(name) >= 3 && strings.EqualFold(name[len(name)-3:], ".gz")
}

// hapWriteCloser bundles an output file with the (optional) BGZF layer so
// callers can flush and close in one step.
type hapWriteCloser struct {
	w      io.Writer
	bgz    *bgzf.Writer
	file   *os.File
	closed bool
}

// Close flushes the BGZF layer (if any) and closes the underlying file.
// It is idempotent so callers may both `defer Close()` and call it
// explicitly to surface the flush error.
func (h *hapWriteCloser) Close() error {
	if h.closed {
		return nil
	}
	h.closed = true
	var err error
	if h.bgz != nil {
		err = h.bgz.Close()
	}
	if cerr := h.file.Close(); err == nil {
		err = cerr
	}
	return err
}

// openHapOutput creates name and wraps it in a BGZF writer when the name
// ends in ".gz", matching upstream's compression choice.
func openHapOutput(name string) (*hapWriteCloser, error) {
	f, err := os.Create(name)
	if err != nil {
		return nil, err
	}
	hc := &hapWriteCloser{file: f}
	if isGzName(name) {
		hc.bgz = bgzf.NewWriter(f)
		hc.w = hc.bgz
	} else {
		hc.w = f
	}
	return hc, nil
}

// gtFields splits a VCF GT string into its allele tokens and reports
// whether the genotype is phased. A genotype is phased iff its separator
// is '|'; haploid genotypes (a single allele) are reported as phased,
// matching how upstream treats a vector-end second allele.
func gtFields(gt string) (alleles []string, phased bool) {
	if gt == "" || gt == "." {
		return []string{"."}, true
	}
	// VCF allows the two alleles to be separated by either '|' or '/'.
	if i := strings.IndexAny(gt, "|/"); i >= 0 {
		phased = gt[i] == '|'
		return []string{gt[:i], gt[i+1:]}, phased
	}
	return []string{gt}, true
}

// gtToHap renders one sample's GT string into the IMPUTE2 hap encoding
// used by upstream's process_gt_to_hap / process_gt_to_hap2. It returns
// the two space-separated haplotype tokens (without a trailing space).
// When hap2Dip is true haploid genotypes are emitted as diploid
// homozygotes.
func gtToHap(gt string, hap2Dip bool) string {
	alleles, phased := gtFields(gt)

	if len(alleles) == 1 {
		a := alleles[0]
		if a == "." {
			// A missing haploid genotype encodes as bcf_gt_missing (0).
			// In hap2dip mode upstream's process_gt_to_hap2 reaches its
			// bcf_gt_is_missing(ptr[0]) branch and emits "? ?". In the
			// plain process_gt_to_hap path the haploid value falls
			// through to kputw(bcf_gt_allele(0)) == -1, emitting "-1 -".
			if hap2Dip {
				return "? ?"
			}
			return "-1 -"
		}
		if hap2Dip {
			// Haploid -> diploid homozygote: "0" -> "0 0".
			return a + " " + a
		}
		// Haploid: second column is the missing-allele marker '-'.
		return a + " -"
	}

	a0, a1 := alleles[0], alleles[1]

	// A missing first allele collapses the whole genotype to "? ?",
	// except the haploid-missing case ("." with a vector-end second
	// allele) which upstream renders "? -" for the non-hap2dip path.
	if a0 == "." {
		if a1 == "" {
			if hap2Dip {
				return "? ?"
			}
			return "? -"
		}
		return "? ?"
	}
	// A missing second allele (e.g. "0/." or "0|.") also yields "? ?".
	if a1 == "." {
		return "? ?"
	}

	if phased {
		return a0 + " " + a1
	}
	// Unphased diploid: each allele is suffixed with '*'.
	return a0 + "* " + a1 + "*"
}

// VCFToHapSample implements `bcftools convert --hapsample`: it reads a VCF
// or BCF from path, applies the sample/region/expression filters in opts,
// and writes the WTCCC-style .hap and .samples files. It returns the
// number of records written. status messages are written to stderr.
func VCFToHapSample(path string, opts HapConvertOptions, stderr io.Writer) (int, error) {
	hdr, variants, err := loadFilteredVCF(path, &opts)
	if err != nil {
		return 0, err
	}

	hapName, sampleName, err := hapOutputNames(opts.Prefix)
	if err != nil {
		return 0, err
	}
	if hapName != "" {
		fmt.Fprintf(stderr, "Hap file: %s\n", hapName)
	}
	if sampleName != "" {
		fmt.Fprintf(stderr, "Sample file: %s\n", sampleName)
	}

	if sampleName != "" {
		if err := writeHapSampleFile(sampleName, hdr.Samples); err != nil {
			return 0, err
		}
	}
	if hapName == "" {
		return 0, nil
	}

	hout, err := openHapOutput(hapName)
	if err != nil {
		return 0, err
	}
	defer hout.Close()

	w := bufio.NewWriter(hout.w)
	nok, skipped, err := streamHapLines(w, hdr, variants, opts.Hap2Dip, func(v *vcf.Variant, gtCols string) error {
		first := firstAlt(v)
		id := fmt.Sprintf("%s:%d_%s_%s", v.Chrom, v.Pos, v.Ref, first)
		// "%CHROM %CHROM:%POS\_%REF\_%FIRST_ALT %POS %REF %FIRST_ALT GT..."
		_, err := fmt.Fprintf(w, "%s %s %d %s %s %s\n", v.Chrom, id, v.Pos, v.Ref, first, gtCols)
		return err
	}, stderr)
	if err != nil {
		return nok, err
	}
	if err := w.Flush(); err != nil {
		return nok, err
	}
	reportHapSkips(stderr, nok, skipped)
	return nok, hout.Close()
}

// VCFToHapLegendSample implements `bcftools convert --haplegendsample`: it
// reads a VCF/BCF from path, applies the filters in opts, and writes the
// IMPUTE2 .hap, .legend and .samples files. It returns the number of
// records written.
func VCFToHapLegendSample(path string, opts HapConvertOptions, stderr io.Writer) (int, error) {
	hdr, variants, err := loadFilteredVCF(path, &opts)
	if err != nil {
		return 0, err
	}

	hapName, legendName, sampleName, err := hapLegendOutputNames(opts.Prefix)
	if err != nil {
		return 0, err
	}
	if hapName != "" {
		fmt.Fprintf(stderr, "Hap file: %s\n", hapName)
	}
	if legendName != "" {
		fmt.Fprintf(stderr, "Legend file: %s\n", legendName)
	}
	if sampleName != "" {
		fmt.Fprintf(stderr, "Sample file: %s\n", sampleName)
	}

	if sampleName != "" {
		if err := writeHapLegendSampleFile(sampleName, hdr.Samples); err != nil {
			return 0, err
		}
	}
	if hapName == "" && legendName == "" {
		return 0, nil
	}

	var hapW, legW *bufio.Writer
	var hout, lout *hapWriteCloser
	if hapName != "" {
		hout, err = openHapOutput(hapName)
		if err != nil {
			return 0, err
		}
		defer hout.Close()
		hapW = bufio.NewWriter(hout.w)
	}
	if legendName != "" {
		lout, err = openHapOutput(legendName)
		if err != nil {
			return 0, err
		}
		defer lout.Close()
		legW = bufio.NewWriter(lout.w)
		if _, err := legW.WriteString("id position a0 a1\n"); err != nil {
			return 0, err
		}
	}

	nok, skipped, err := streamHapLines(hapW, hdr, variants, opts.Hap2Dip, func(v *vcf.Variant, gtCols string) error {
		if hapW != nil {
			if _, err := hapW.WriteString(gtCols); err != nil {
				return err
			}
			if err := hapW.WriteByte('\n'); err != nil {
				return err
			}
		}
		if legW != nil {
			first := firstAlt(v)
			id := fmt.Sprintf("%s:%d_%s_%s", v.Chrom, v.Pos, v.Ref, first)
			if _, err := fmt.Fprintf(legW, "%s %d %s %s\n", id, v.Pos, v.Ref, first); err != nil {
				return err
			}
		}
		return nil
	}, stderr)
	if err != nil {
		return nok, err
	}
	if hapW != nil {
		if err := hapW.Flush(); err != nil {
			return nok, err
		}
	}
	if legW != nil {
		if err := legW.Flush(); err != nil {
			return nok, err
		}
	}
	reportHapSkips(stderr, nok, skipped)
	if hout != nil {
		if err := hout.Close(); err != nil {
			return nok, err
		}
	}
	if lout != nil {
		if err := lout.Close(); err != nil {
			return nok, err
		}
	}
	return nok, nil
}

// hapSkips tracks the three skip reasons upstream reports at the end of a
// VCF->hap conversion.
type hapSkips struct {
	noAlt        int
	nonBiallelic int
	filtered     int
}

// streamHapLines is the shared body of the VCF->hap exporters. For each
// biallelic variant it builds the per-sample hap columns and invokes emit
// with the variant and the rendered column string. Multiallelic and
// no-ALT records are skipped (counted), matching upstream. When hapW is
// nil the hap columns are still rendered (legend-only output still needs
// to walk the records).
func streamHapLines(hapW *bufio.Writer, hdr *vcf.Header, variants []*vcf.Variant, hap2Dip bool, emit func(*vcf.Variant, string) error, stderr io.Writer) (int, hapSkips, error) {
	var skips hapSkips
	nok := 0
	warnedMulti := false
	var sb strings.Builder
	for _, v := range variants {
		// ALT allele is required.
		if len(v.Alt) == 0 || (len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "")) {
			skips.noAlt++
			continue
		}
		// Biallelic required.
		if len(v.Alt) > 1 {
			if !warnedMulti {
				fmt.Fprintln(stderr, "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'bcftools norm -m-'.")
				warnedMulti = true
			}
			skips.nonBiallelic++
			continue
		}

		sb.Reset()
		gtIdx := formatIndex(v.Format, "GT")
		if gtIdx < 0 {
			return nok, skips, fmt.Errorf("FORMAT/GT tag not present at %s:%d", v.Chrom, v.Pos)
		}
		for i, s := range v.Samples {
			if i > 0 {
				sb.WriteByte(' ')
			}
			sb.WriteString(gtToHap(s.Data["GT"], hap2Dip))
		}
		if err := emit(v, sb.String()); err != nil {
			return nok, skips, err
		}
		nok++
	}
	return nok, skips, nil
}

// reportHapSkips prints upstream's per-conversion summary line to stderr.
func reportHapSkips(stderr io.Writer, nok int, s hapSkips) {
	fmt.Fprintf(stderr, "%d records written, %d skipped: %d/%d/%d no-ALT/non-biallelic/filtered\n",
		nok, s.noAlt+s.nonBiallelic+s.filtered, s.noAlt, s.nonBiallelic, s.filtered)
}

// firstAlt returns the first ALT allele of v (the only ALT for the
// biallelic records the hap exporters accept).
func firstAlt(v *vcf.Variant) string {
	if len(v.Alt) == 0 {
		return "."
	}
	return v.Alt[0]
}

// formatIndex returns the position of tag in the FORMAT list, or -1.
func formatIndex(format []string, tag string) int {
	for i, f := range format {
		if f == tag {
			return i
		}
	}
	return -1
}

// writeHapSampleFile writes the SHAPEIT-style .sample file used by
// --hapsample (header "ID_1 ID_2 missing", then "0 0 0", then one row per
// sample).
func writeHapSampleFile(name string, samples []string) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	var w io.Writer = f
	var bgz *bgzf.Writer
	if isGzName(name) {
		bgz = bgzf.NewWriter(f)
		w = bgz
	}
	bw := bufio.NewWriter(w)
	if _, err := bw.WriteString("ID_1 ID_2 missing\n0 0 0\n"); err != nil {
		return err
	}
	for _, s := range samples {
		if _, err := fmt.Fprintf(bw, "%s %s 0\n", s, s); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	if bgz != nil {
		if err := bgz.Close(); err != nil {
			return err
		}
	}
	return f.Close()
}

// writeHapLegendSampleFile writes the IMPUTE2-style .sample file used by
// --haplegendsample (header "sample population group sex", then one row
// per sample with sex defaulting to '2').
func writeHapLegendSampleFile(name string, samples []string) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	var w io.Writer = f
	var bgz *bgzf.Writer
	if isGzName(name) {
		bgz = bgzf.NewWriter(f)
		w = bgz
	}
	bw := bufio.NewWriter(w)
	if _, err := bw.WriteString("sample population group sex\n"); err != nil {
		return err
	}
	for _, s := range samples {
		if _, err := fmt.Fprintf(bw, "%s %s %s 2\n", s, s, s); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	if bgz != nil {
		if err := bgz.Close(); err != nil {
			return err
		}
	}
	return f.Close()
}

// loadFilteredVCF reads path, applies opts' sample restriction, region
// post-filter and include/exclude expressions, and returns the resulting
// header and variants. It is shared by the two VCF->hap exporters.
func loadFilteredVCF(path string, opts *HapConvertOptions) (*vcf.Header, []*vcf.Variant, error) {
	if opts.SamplesFile != "" {
		names, err := LoadSamplesFile(opts.SamplesFile)
		if err != nil {
			return nil, nil, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Samples = append(opts.Samples, names...)
	}
	if opts.RegionsFile != "" {
		regs, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return nil, nil, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Regions = append(opts.Regions, regs...)
	}
	if opts.TargetsFile != "" {
		regs, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return nil, nil, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Targets = append(opts.Targets, regs...)
	}

	r, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, nil, fmt.Errorf("bcftools convert: open %s: %w", path, err)
	}
	defer r.Close()
	hdr, variants, err := readAllVariants(r)
	if err != nil {
		return nil, nil, err
	}

	if len(opts.Samples) > 0 {
		missing := missingSamples(hdr.Samples, opts.Samples)
		if len(missing) > 0 && !opts.ForceSamples {
			return nil, nil, fmt.Errorf("bcftools convert: requested samples missing from input: %s (use --force-samples to ignore)", strings.Join(missing, ", "))
		}
		hdr = filterHeaderSamples(hdr, opts.Samples)
	}

	postFilters := append([]string{}, opts.Targets...)
	postFilters = append(postFilters, opts.Regions...)
	parsedTargets, err := parseRegions(postFilters)
	if err != nil {
		return nil, nil, fmt.Errorf("bcftools convert: %w", err)
	}

	include, exclude, err := compileExpressions(ViewOptions{
		IncludeExpr: opts.IncludeExpr,
		ExcludeExpr: opts.ExcludeExpr,
	})
	if err != nil {
		return nil, nil, fmt.Errorf("bcftools convert: %w", err)
	}

	kept := variants[:0]
	for _, v := range variants {
		if len(parsedTargets) > 0 && !overlapsAny(v, parsedTargets) {
			continue
		}
		if include != nil && !include.Eval(v) {
			continue
		}
		if exclude != nil && exclude.Eval(v) {
			continue
		}
		if len(opts.Samples) > 0 {
			restrictSamples(v, opts.Samples)
		}
		kept = append(kept, v)
	}
	return hdr, kept, nil
}

// HapSampleToVCF implements `bcftools convert --hapsample2vcf`: it reads a
// SHAPEIT-style .hap (+.samples) pair and emits a VCF/BCF with filled GT.
// It returns the number of records processed.
func HapSampleToVCF(arg string, out io.Writer, opts HapConvertOptions, stderr io.Writer) (int, error) {
	hapName, sampleName := hapInputNames(arg)

	sampleNames, err := readHapSampleNames(sampleName, 2)
	if err != nil {
		return 0, err
	}

	hr, err := iohelper.OpenReader(hapName)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: open %s: %w", hapName, err)
	}
	defer hr.Close()
	sc := bufio.NewScanner(hr)
	sc.Buffer(make([]byte, 0, 64<<10), 64<<20)

	if !sc.Scan() {
		return 0, fmt.Errorf("empty file: %s", hapName)
	}
	first := sc.Text()
	chrom, err := chromFromHapLine(first)
	if err != nil {
		return 0, fmt.Errorf("%w in %s", err, hapName)
	}

	hdr := newHapVCFHeader(chrom, sampleNames)
	w, finish, err := openOutput(out, ViewOptions{OutputFormat: opts.OutputFormat, CompressLevel: opts.CompressLevel}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	total := 0
	line := first
	for {
		v, perr := parseHapSampleLine(line, chrom, sampleNames)
		if perr != nil {
			return total, perr
		}
		if err := w.Write(v); err != nil {
			return total, err
		}
		total++
		if !sc.Scan() {
			break
		}
		line = sc.Text()
		if strings.TrimSpace(line) == "" {
			break
		}
	}
	if err := sc.Err(); err != nil {
		return total, err
	}
	fmt.Fprintf(stderr, "Number of processed rows: \t%d\n", total)
	return total, w.Flush()
}

// HapLegendSampleToVCF implements `bcftools convert --haplegendsample2vcf`:
// it reads an IMPUTE2 .hap/.legend/.samples triple and emits a VCF/BCF
// with filled GT. It returns the number of records processed.
func HapLegendSampleToVCF(arg string, out io.Writer, opts HapConvertOptions, stderr io.Writer) (int, error) {
	hapName, legendName, sampleName, err := hapLegendInputNames(arg)
	if err != nil {
		return 0, err
	}

	sampleNames, err := readHapSampleNames(sampleName, 1)
	if err != nil {
		return 0, err
	}

	lr, err := iohelper.OpenReader(legendName)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: open %s: %w", legendName, err)
	}
	defer lr.Close()
	legSc := bufio.NewScanner(lr)
	legSc.Buffer(make([]byte, 0, 64<<10), 64<<20)

	// Eat the legend header line.
	if !legSc.Scan() {
		return 0, fmt.Errorf("empty file: %s", legendName)
	}
	if !legSc.Scan() {
		return 0, fmt.Errorf("empty file: %s", legendName)
	}
	firstLeg := legSc.Text()
	chrom, err := chromFromLegendLine(firstLeg)
	if err != nil {
		return 0, fmt.Errorf("%w in %s", err, legendName)
	}

	hr, err := iohelper.OpenReader(hapName)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: open %s: %w", hapName, err)
	}
	defer hr.Close()
	hapSc := bufio.NewScanner(hr)
	hapSc.Buffer(make([]byte, 0, 64<<10), 64<<20)

	hdr := newHapVCFHeader(chrom, sampleNames)
	w, finish, err := openOutput(out, ViewOptions{OutputFormat: opts.OutputFormat, CompressLevel: opts.CompressLevel}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	total := 0
	legLine := firstLeg
	for {
		v, revAls, perr := parseLegendLine(legLine, chrom)
		if perr != nil {
			return total, perr
		}
		if !hapSc.Scan() {
			return total, fmt.Errorf("different number of records in %s and %s", legendName, hapName)
		}
		if err := fillHapGenotypes(v, strings.Fields(hapSc.Text()), sampleNames, revAls); err != nil {
			return total, err
		}
		if err := w.Write(v); err != nil {
			return total, err
		}
		total++

		if !legSc.Scan() {
			if hapSc.Scan() {
				return total, fmt.Errorf("different number of records in %s and %s", legendName, hapName)
			}
			break
		}
		legLine = legSc.Text()
		if strings.TrimSpace(legLine) == "" {
			break
		}
	}
	if err := legSc.Err(); err != nil {
		return total, err
	}
	fmt.Fprintf(stderr, "Number of processed rows: \t%d\n", total)
	return total, w.Flush()
}

// readHapSampleNames reads a .sample/.samples file and returns the sample
// IDs. skipRows is the number of leading header rows to drop (2 for the
// SHAPEIT format used by --hapsample, 1 for the IMPUTE2 format used by
// --haplegendsample). The sample ID is the first whitespace-delimited
// column of each remaining row.
func readHapSampleNames(name string, skipRows int) ([]string, error) {
	r, err := iohelper.OpenReader(name)
	if err != nil {
		return nil, fmt.Errorf("bcftools convert: open %s: %w", name, err)
	}
	defer r.Close()
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var rows []string
	for sc.Scan() {
		rows = append(rows, sc.Text())
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if len(rows) < skipRows {
		return nil, fmt.Errorf("could not read %s", name)
	}
	var names []string
	for _, row := range rows[skipRows:] {
		fields := strings.Fields(row)
		if len(fields) == 0 {
			continue
		}
		names = append(names, fields[0])
	}
	return names, nil
}

// newHapVCFHeader builds the minimal VCF header upstream emits for the
// hap*2vcf modes: an END INFO tag, a GT FORMAT tag, a single contig with
// the maximum CSI length, and the sample columns.
func newHapVCFHeader(chrom string, samples []string) *vcf.Header {
	return &vcf.Header{
		MetaInfo: []string{
			"##fileformat=VCFv4.2",
			`##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">`,
			`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
			fmt.Sprintf("##contig=<ID=%s,length=%d>", chrom, 0x7fffffff),
		},
		Samples: samples,
	}
}

// chromFromHapLine extracts the chromosome name from the first column of a
// SHAPEIT .hap line. The first column is the bare CHROM and the second is
// CHROM:POS_REF_ALT; upstream reads CHROM from the leading column.
func chromFromHapLine(line string) (string, error) {
	fields := strings.Fields(line)
	if len(fields) < 2 {
		return "", fmt.Errorf("could not determine CHROM")
	}
	// The CHROM[:POS_REF_ALT] string is in the second column.
	sb := fields[1]
	c := strings.IndexByte(sb, ':')
	if c < 0 {
		return "", fmt.Errorf("could not determine CHROM in the second column")
	}
	return sb[:c], nil
}

// chromFromLegendLine extracts the chromosome name from a legend data line
// whose first column is CHROM:POS_REF_ALT.
func chromFromLegendLine(line string) (string, error) {
	c := strings.IndexByte(line, ':')
	if c < 0 {
		return "", fmt.Errorf("expected CHROM:POS_REF_ALT in first column")
	}
	return line[:c], nil
}

// parseChromPosRefAlt splits a CHROM:POS_REF_ALT[_END] token into its
// components. END, when present, is returned in end (0 means absent).
func parseChromPosRefAlt(tok string) (chrom string, pos int, ref, alt string, end int, err error) {
	colon := strings.IndexByte(tok, ':')
	if colon < 0 {
		return "", 0, "", "", 0, fmt.Errorf("could not parse the CHROM:POS_REF_ALT[_END] string: %s", tok)
	}
	chrom = tok[:colon]
	rest := tok[colon+1:]
	// POS is the run of digits up to the first '_'.
	us := strings.IndexByte(rest, '_')
	if us < 0 {
		return "", 0, "", "", 0, fmt.Errorf("could not parse the CHROM:POS_REF_ALT[_END] string: %s", tok)
	}
	pos, err = strconv.Atoi(rest[:us])
	if err != nil {
		return "", 0, "", "", 0, fmt.Errorf("could not parse the CHROM:POS_REF_ALT[_END] string: %s", tok)
	}
	rest = rest[us+1:]
	us = strings.IndexByte(rest, '_')
	if us < 0 {
		return "", 0, "", "", 0, fmt.Errorf("could not parse the CHROM:POS_REF_ALT[_END] string: %s", tok)
	}
	ref = rest[:us]
	rest = rest[us+1:]
	// ALT runs to the next '_' (END) or end of token.
	if us = strings.IndexByte(rest, '_'); us >= 0 {
		alt = rest[:us]
		end, err = strconv.Atoi(rest[us+1:])
		if err != nil || end < 1 {
			return "", 0, "", "", 0, fmt.Errorf("could not parse the CHROM:POS_REF_ALT[_END] string: %s", tok)
		}
	} else {
		alt = rest
	}
	return chrom, pos, ref, alt, end, nil
}

// parseLegendLine parses a legend data line ("id position a0 a1") into a
// Variant, cross-checking the id's POS/REF/ALT against the dedicated
// columns the way upstream's tsv_setter_verify_pos / verify_ref_alt do.
// It also reports whether the a0/a1 columns are reversed relative to the
// id-derived REF/ALT (rev_als), so the haplotype 0/1 labels can be
// swapped to match.
func parseLegendLine(line, chrom string) (*vcf.Variant, bool, error) {
	fields := strings.Fields(line)
	if len(fields) < 4 {
		return nil, false, fmt.Errorf("could not parse legend line: %s", line)
	}
	c, pos, ref, alt, end, err := parseChromPosRefAlt(fields[0])
	if err != nil {
		return nil, false, err
	}
	if c != chrom {
		// Upstream emits all records under the first chromosome it saw;
		// a different CHROM is a genuine inconsistency.
		return nil, false, fmt.Errorf("CHROM mismatch in legend line: %s", line)
	}
	vpos, err := strconv.Atoi(fields[1])
	if err != nil {
		return nil, false, fmt.Errorf("could not parse POS: %s", fields[1])
	}
	if vpos != pos {
		return nil, false, fmt.Errorf("POS mismatch: %s", line)
	}
	revAls, err := refAltOrientation(ref, alt, fields[2], fields[3])
	if err != nil {
		return nil, false, err
	}
	return newHapVariant(chrom, pos, ref, alt, end), revAls, nil
}

// refAltOrientation compares the (a0,a1) columns of a hap/legend file
// against the REF/ALT taken from the CHROM:POS_REF_ALT id, mirroring
// upstream's tsv_setter_verify_ref_alt. It returns false when they agree,
// true when they are swapped (rev_als), and an error when neither matches.
func refAltOrientation(ref, alt, a0, a1 string) (bool, error) {
	if a0 == ref && a1 == alt {
		return false, nil
	}
	if a0 == alt && a1 == ref {
		return true, nil
	}
	return false, fmt.Errorf("REF/ALT mismatch: [%s][%s]", a0, a1)
}

// parseHapSampleLine parses a SHAPEIT .hap line (CHROM CHROM:POS_REF_ALT
// POS REF ALT hap...) into a Variant with GT filled from the haplotype
// columns.
func parseHapSampleLine(line, chrom string, samples []string) (*vcf.Variant, error) {
	fields := strings.Fields(line)
	if len(fields) < 5 {
		return nil, fmt.Errorf("could not parse: %s", line)
	}
	c, pos, ref, alt, end, err := parseChromPosRefAlt(fields[1])
	if err != nil {
		return nil, err
	}
	if c != chrom {
		return nil, fmt.Errorf("CHROM mismatch: %s", line)
	}
	vpos, err := strconv.Atoi(fields[2])
	if err != nil {
		return nil, fmt.Errorf("could not parse POS: %s", fields[2])
	}
	if vpos != pos {
		return nil, fmt.Errorf("POS mismatch: %s", line)
	}
	revAls, err := refAltOrientation(ref, alt, fields[3], fields[4])
	if err != nil {
		return nil, err
	}
	v := newHapVariant(chrom, pos, ref, alt, end)
	if err := fillHapGenotypes(v, fields[5:], samples, revAls); err != nil {
		return nil, err
	}
	return v, nil
}

// newHapVariant constructs a biallelic Variant with the END INFO tag set
// when end>0, the shape produced by the hap*2vcf importers.
func newHapVariant(chrom string, pos int, ref, alt string, end int) *vcf.Variant {
	v := &vcf.Variant{
		Chrom:   chrom,
		Pos:     pos,
		ID:      ".",
		Ref:     ref,
		Alt:     []string{alt},
		Qual:    -1,
		Filter:  []string{"."},
		Info:    map[string]string{},
		Format:  []string{"GT"},
		Samples: nil,
	}
	if end > 0 {
		v.Info["END"] = strconv.Itoa(end)
		v.InfoOrder = []string{"END"}
	}
	return v
}

// fillHapGenotypes converts 2*nsamples haplotype tokens into per-sample GT
// strings. revAls swaps the 0/1 allele labels (set when the .hap file's
// REF/ALT are reversed relative to the legend). The token grammar mirrors
// upstream's tsv_setter_haps: '0'/'1' are phased alleles, '?' is a phased
// missing allele, '-' is a vector-end (haploid) marker, and a trailing '*'
// on both alleles downgrades the genotype to unphased.
func fillHapGenotypes(v *vcf.Variant, tokens, samples []string, revAls bool) error {
	if len(tokens) != 2*len(samples) {
		return fmt.Errorf("wrong number of hap fields: got %d, want %d", len(tokens), 2*len(samples))
	}
	a0, a1 := "0", "1"
	if revAls {
		a0, a1 = "1", "0"
	}
	v.Samples = make([]vcf.Sample, len(samples))
	for i, name := range samples {
		t0, t1 := tokens[2*i], tokens[2*i+1]
		al0, un0, err := hapToken(t0, a0, a1)
		if err != nil {
			return err
		}
		al1, un1, err := hapToken(t1, a0, a1)
		if err != nil {
			return err
		}
		// Upstream requires the unphased marker on both alleles or
		// neither.
		if un0 != un1 {
			return fmt.Errorf("missing unphased marker '*'")
		}
		var gt string
		switch {
		case al1 == "": // haploid (second allele was '-')
			gt = al0
		case un0:
			gt = al0 + "/" + al1
		default:
			gt = al0 + "|" + al1
		}
		v.Samples[i] = vcf.Sample{Name: name, Data: map[string]string{"GT": gt}}
	}
	return nil
}

// hapToken decodes a single haplotype token into its allele label and
// whether it carried the trailing '*' unphased marker. An empty allele
// label signals the '-' vector-end (haploid) marker.
func hapToken(tok, a0, a1 string) (allele string, unphased bool, err error) {
	if strings.HasSuffix(tok, "*") {
		unphased = true
		tok = tok[:len(tok)-1]
	}
	switch tok {
	case "0":
		return a0, unphased, nil
	case "1":
		return a1, unphased, nil
	case "?":
		return ".", unphased, nil
	case "-":
		return "", unphased, nil
	default:
		return "", false, fmt.Errorf("could not parse hap token: %q", tok)
	}
}
