// Native port of the upstream `contrast` plugin (plugins/contrast.c). It runs a
// basic per-site association test between two groups of samples (controls -0 and
// cases -1) and flags novel alleles and genotypes seen only in the case group,
// adding INFO annotations PASSOC, FASSOC, NASSOC, NOVELAL and NOVELGT. The VCF
// output is emitted unchanged except for the added INFO fields; a one-line
// summary is written to stderr at the end.
//
// A single -i/--include or -e/--exclude filter expression is supported as a
// site-level pre-filter (upstream's contrast.c calls filter_test with a NULL
// sample mask and drops the record entirely before annotating or writing it).
// The -o/-O options select the output file and container, and -W/--write-index
// indexes that file (a CSI by default, a TBI for -W=tbi on VCF.gz) exactly as
// upstream does; -W to stdout is rejected with upstream's "failed to initialise
// index for -" error.
//
// The rare-allele enrichment mode (-f/--max-allele-freq) is supported: the
// per-site VCF output and PASSOC/FASSOC/NASSOC/NOVEL* annotations are unchanged,
// and in addition the region-wide pooled minor-allele counts (folded over the
// records whose minor allele is at or below the -f threshold, exactly as
// contrast.c does) feed a second stderr summary line "max_AC/PASSOC/FASSOC/
// NASSOC:" with the Fisher's exact probability and control/case non-REF
// fractions. An integer -f argument is a raw allele-count threshold; a float in
// [0,1] is an allele-frequency threshold scaled by the total sample count.
//
// The --regions-overlap / --targets-overlap region-matching modes remain
// unsupported (the native region/target filter does not replicate htslib's
// overlap semantics).
package bcftools

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("contrast", func() NativePlugin { return &contrastPlugin{} }) }

// contrast annotation flags, matching PRINT_* in contrast.c.
const (
	contrastPASSOC  = 1 << 0
	contrastFASSOC  = 1 << 1
	contrastNASSOC  = 1 << 2
	contrastNOVELAL = 1 << 3
	contrastNOVELGT = 1 << 4
)

// contrastPlugin implements the `contrast` plugin (per-record, but with running
// totals for the stderr summary, so it runs serially).
type contrastPlugin struct {
	annots      int
	controlIdx  []int
	caseIdx     []int
	forceSample bool

	ntotal, nskipped, ntested, ncaseAl, ncaseGt int

	// Rare-allele enrichment mode (-f/--max-allele-freq). maxAC is the resolved
	// minor-allele-count threshold (an integer argument is taken verbatim; a
	// float in [0,1] is multiplied by the total sample count, floored, min 1).
	// maxACSet records that -f was given. enrichNals accumulates, region-wide,
	// the control-ref/control-alt/case-ref/case-alt allele counts of the records
	// whose minor allele is at or below the threshold, mirroring contrast.c's
	// args->nals folding. The accumulated counts feed the extra stderr summary
	// line printed at end of run.
	maxAC      int
	maxACSet   bool
	enrichNals [4]int

	filter *pluginFilter // compiled -i/-e site-level pre-filter, nil if none

	outputFile string        // -o/--output FILE; "" or "-" means stdout
	format     OutputFormat  // -O/--output-type; defaults to VCF
	clevel     int           // -O z/b level; -1 for the package default
	writeIndex writeIndexFmt // -W/--write-index[=FMT]; writeIndexOff when unset
	outputSet  bool          // whether any of -o/-O/-W was given

	stderr io.Writer
}

// OutputControl reports the -o/-O/-W selection contrast parsed from its argv so
// the framework's stage-3 writer targets the chosen file and container and
// indexes it. It implements pluginOutputControl.
func (p *contrastPlugin) OutputControl() (string, OutputFormat, int, writeIndexFmt, bool) {
	return p.outputFile, p.format, p.clevel, p.writeIndex, p.outputSet
}

// Name returns the plugin name.
func (p *contrastPlugin) Name() string { return "contrast" }

// RegionTargetCaps opts contrast into the shared -r/-R/-t/-T region/target
// filter, applied to the records before the case/control contrast accounting.
func (p *contrastPlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching contrast.c about().
func (p *contrastPlugin) About() string {
	return "Find novel alleles and genotypes in two groups of samples.\n"
}

// Parallel reports false: the running summary totals are updated serially.
func (p *contrastPlugin) Parallel() bool { return false }

// RunStyle reports that contrast is a run()-style plugin (options precede the
// file, no `--` separator).
func (p *contrastPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of contrast's value-taking flags consumes
// the following token, used by the host to split the input-file positional.
func (p *contrastPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-a", "--annots", "-0", "--control-samples", "--bg-samples",
		"-1", "--case-samples", "--novel-samples", "-i", "--include", "-e", "--exclude",
		"-f", "--max-allele-freq", "-o", "--output", "-O", "--output-type",
		"-r", "--regions", "-R", "--regions-file", "-t", "--targets", "-T", "--targets-file",
		"--regions-overlap", "--targets-overlap", "-v", "--verbosity":
		return true
	}
	return false
}

// SetStderr wires the host stderr writer the summary line is printed to.
func (p *contrastPlugin) SetStderr(w io.Writer) { p.stderr = w }

// Init parses the options, builds the control/case sample-index lists, appends
// the requested INFO header lines, and rejects the unsupported modes.
func (p *contrastPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	annotsStr := "PASSOC,FASSOC"
	var controlStr, caseStr string
	var filterExpr string
	var filterExclude, haveFilter bool
	var maxACStr string
	p.format = OutputVCF
	p.clevel = -1
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("contrast: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		if sel, handled, werr := parseWriteIndexArg(a); handled {
			if werr != nil {
				return nil, fmt.Errorf("contrast: %w", werr)
			}
			p.writeIndex = sel
			p.outputSet = true
			continue
		}
		switch a {
		case "-a", "--annots":
			v, err := next()
			if err != nil {
				return nil, err
			}
			annotsStr = v
		case "-0", "--control-samples", "--bg-samples":
			v, err := next()
			if err != nil {
				return nil, err
			}
			controlStr = v
		case "-1", "--case-samples", "--novel-samples":
			v, err := next()
			if err != nil {
				return nil, err
			}
			caseStr = v
		case "--force-samples":
			p.forceSample = true
		case "-i", "--include", "-e", "--exclude":
			if haveFilter {
				return nil, fmt.Errorf("contrast: only one -i/--include or -e/--exclude expression can be given, and they cannot be combined")
			}
			v, err := next()
			if err != nil {
				return nil, err
			}
			filterExpr = v
			filterExclude = a == "-e" || a == "--exclude"
			haveFilter = true
		case "-f", "--max-allele-freq":
			v, err := next()
			if err != nil {
				return nil, err
			}
			maxACStr = v
		case "--regions-overlap", "--targets-overlap":
			return nil, fmt.Errorf("contrast: %s is not supported by the native plugin", a)
		case "-o", "--output":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.outputFile = v
			p.outputSet = true
		case "-O", "--output-type":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if err := p.parseOutputType(v); err != nil {
				return nil, err
			}
			p.outputSet = true
		case "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		default:
			// Attached -O<x> / -o<path> getopt forms (e.g. `-Oz`, `-oout.vcf`).
			if strings.HasPrefix(a, "-O") && len(a) > 2 {
				if err := p.parseOutputType(a[2:]); err != nil {
					return nil, err
				}
				p.outputSet = true
				continue
			}
			if strings.HasPrefix(a, "-o") && len(a) > 2 {
				p.outputFile = a[2:]
				p.outputSet = true
				continue
			}
			return nil, fmt.Errorf("contrast: unsupported option %q", a)
		}
	}
	if controlStr == "" {
		return nil, fmt.Errorf("contrast: missing the -0, --control-samples option")
	}
	if caseStr == "" {
		return nil, fmt.Errorf("contrast: missing the -1, --case-samples option")
	}

	for _, tok := range strings.Split(annotsStr, ",") {
		switch strings.ToUpper(strings.TrimSpace(tok)) {
		case "PASSOC":
			p.annots |= contrastPASSOC
		case "FASSOC":
			p.annots |= contrastFASSOC
		case "NASSOC":
			p.annots |= contrastNASSOC
		case "NOVELAL":
			p.annots |= contrastNOVELAL
		case "NOVELGT":
			p.annots |= contrastNOVELGT
		default:
			return nil, fmt.Errorf("contrast: the annotation is not recognised: %s", tok)
		}
	}

	var err error
	p.controlIdx, err = resolveSampleGroup(hdr, controlStr, p.forceSample)
	if err != nil {
		return nil, err
	}
	p.caseIdx, err = resolveSampleGroup(hdr, caseStr, p.forceSample)
	if err != nil {
		return nil, err
	}

	if haveFilter {
		f, ferr := newPluginFilterWithHeader(filterExpr, filterExclude, hdr)
		if ferr != nil {
			return nil, fmt.Errorf("contrast: %w", ferr)
		}
		p.filter = f
	}

	if maxACStr != "" {
		ac, err := parseContrastMaxAC(maxACStr, len(hdr.Samples))
		if err != nil {
			return nil, err
		}
		p.maxAC = ac
		p.maxACSet = true
	}

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	if p.annots&contrastPASSOC != 0 {
		out.MetaInfo = append(out.MetaInfo, `##INFO=<ID=PASSOC,Number=1,Type=Float,Description="Fisher's exact test probability of genotypic association (REF vs non-REF allele)">`)
	}
	if p.annots&contrastFASSOC != 0 {
		out.MetaInfo = append(out.MetaInfo, `##INFO=<ID=FASSOC,Number=2,Type=Float,Description="Proportion of non-REF allele in controls and cases">`)
	}
	if p.annots&contrastNASSOC != 0 {
		out.MetaInfo = append(out.MetaInfo, `##INFO=<ID=NASSOC,Number=4,Type=Integer,Description="Number of control-ref, control-alt, case-ref and case-alt alleles">`)
	}
	if p.annots&contrastNOVELAL != 0 {
		out.MetaInfo = append(out.MetaInfo, `##INFO=<ID=NOVELAL,Number=.,Type=String,Description="List of samples with novel alleles. Note that samples listed here are not listed in NOVELGT again.">`)
	}
	if p.annots&contrastNOVELGT != 0 {
		out.MetaInfo = append(out.MetaInfo, `##INFO=<ID=NOVELGT,Number=.,Type=String,Description="List of samples with novel genotypes">`)
	}
	return out, nil
}

// parseContrastMaxAC resolves the -f/--max-allele-freq argument into a
// minor-allele-count threshold, mirroring contrast.c init_data(): a clean
// integer is taken verbatim; otherwise a float in [0,1] is multiplied by the
// total sample count and floored, with a floor of 1 when the product rounds to
// zero. A value that is neither a clean integer nor a float in [0,1] is an
// error.
func parseContrastMaxAC(s string, nsamples int) (int, error) {
	// strtol-style: accept a clean base-10 integer (the whole string).
	if n, err := strconv.Atoi(s); err == nil {
		return n, nil
	}
	val, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, fmt.Errorf("contrast: could not parse the argument: -f, --max-allele-freq %s", s)
	}
	if val < 0 || val > 1 {
		return 0, fmt.Errorf("contrast: expected integer or float from the range [0,1]: -f, --max-allele-freq %s", s)
	}
	ac := int(val * float64(nsamples))
	if ac == 0 {
		ac = 1
	}
	return ac, nil
}

// parseOutputType parses upstream's -O u|b|v|z[0-9] spelling into the contrast
// plugin's container/level selection.
func (p *contrastPlugin) parseOutputType(s string) error {
	if s == "" {
		return fmt.Errorf("contrast: empty -O argument")
	}
	switch s[0] {
	case 'b':
		p.format = OutputBCF
	case 'u':
		p.format = OutputBCFUncompressed
	case 'z':
		p.format = OutputVCFGz
	case 'v':
		p.format = OutputVCF
	default:
		n, err := strconv.Atoi(s)
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("contrast: the output type %q not recognised", s)
		}
		p.clevel = n
		return nil
	}
	if len(s) > 1 {
		n, err := strconv.Atoi(s[1:])
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("contrast: could not parse compression level %q", s[1:])
		}
		p.clevel = n
	}
	return nil
}

// Process annotates one record with the requested INFO fields, mirroring
// contrast.c process_record(). The record is always emitted (the upstream
// run-loop writes every record); skipped records carry no new INFO.
func (p *contrastPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	// Upstream applies the -i/-e site-level filter before process_record and
	// before writing, so a filtered-out record is dropped from the output and
	// excluded from every summary counter (including ntotal).
	if !p.filter.testSite(v) {
		return nil, nil
	}
	p.ntotal++

	// Control group: gather the union of alleles, the per-genotype bitmask set
	// (for NOVELGT), and ref/alt allele counts.
	controlAls := uint32(0)
	var nals [4]int // ctrl-ref, ctrl-alt, case-ref, case-alt
	controlGTs := map[uint32]struct{}{}
	for _, si := range p.controlIdx {
		gt, ok := sampleGT(v, si)
		if !ok {
			continue
		}
		var bits uint32
		for _, a := range gt.alleles {
			if a == missingAllele {
				continue
			}
			if a > 31 {
				p.nskipped++
				return []*vcf.Variant{v}, nil
			}
			controlAls |= 1 << uint(a)
			bits |= 1 << uint(a)
			if a != 0 {
				nals[1]++
			} else {
				nals[0]++
			}
		}
		if p.annots&contrastNOVELGT != 0 {
			controlGTs[bits] = struct{}{}
		}
	}
	if controlAls == 0 && len(p.controlIdx) > 0 {
		p.nskipped++
		return []*vcf.Variant{v}, nil
	}

	var caseAlSmpl, caseGtSmpl []string
	hasGT := false
	for _, si := range p.caseIdx {
		gt, ok := sampleGT(v, si)
		if !ok {
			continue
		}
		caseAl := false
		var bits uint32
		for _, a := range gt.alleles {
			if a == missingAllele {
				continue
			}
			if a > 31 {
				p.nskipped++
				return []*vcf.Variant{v}, nil
			}
			if controlAls&(1<<uint(a)) == 0 {
				caseAl = true
			}
			bits |= 1 << uint(a)
			if a != 0 {
				nals[3]++
			} else {
				nals[2]++
			}
		}
		if bits == 0 {
			continue
		}
		hasGT = true
		name := v.Samples[si].Name
		if caseAl && p.annots&contrastNOVELAL != 0 {
			caseAlSmpl = append(caseAlSmpl, name)
		} else if p.annots&contrastNOVELGT != 0 {
			if _, seen := controlGTs[bits]; !seen {
				caseGtSmpl = append(caseGtSmpl, name)
			}
		}
	}
	if !hasGT && len(p.caseIdx) > 0 {
		p.nskipped++
		return []*vcf.Variant{v}, nil
	}

	// Rare-allele enrichment (-f) folding, mirroring contrast.c. The minor allele
	// (the rarer of REF / non-REF across the trio of control+case counts) is
	// pooled into the region-wide enrichNals only when its count is at or below
	// the threshold. When the non-REF allele is the minor one the counts are
	// added verbatim; when REF is the minor one the ref/alt columns are swapped
	// so enrichNals always tracks "minor=alt".
	if p.maxACSet {
		if nals[0]+nals[2] > nals[1]+nals[3] {
			if nals[1]+nals[3] <= p.maxAC {
				for i := 0; i < 4; i++ {
					p.enrichNals[i] += nals[i]
				}
			}
		} else {
			if nals[0]+nals[2] <= p.maxAC {
				p.enrichNals[0] += nals[1]
				p.enrichNals[1] += nals[0]
				p.enrichNals[2] += nals[3]
				p.enrichNals[3] += nals[2]
			}
		}
	}

	if p.annots&contrastPASSOC != 0 && len(p.controlIdx) > 0 && len(p.caseIdx) > 0 {
		_, _, fisher := mpileupFisherExact(int64(nals[0]), int64(nals[1]), int64(nals[2]), int64(nals[3]))
		setInfo(v, "PASSOC", formatVCFFloat(fisher))
	}
	if p.annots&contrastFASSOC != 0 && len(p.controlIdx) > 0 && len(p.caseIdx) > 0 {
		var c0, c1 string
		if nals[0]+nals[1] != 0 {
			c0 = formatVCFFloat(float64(float32(nals[1]) / float32(nals[0]+nals[1])))
		} else {
			c0 = "."
		}
		if nals[2]+nals[3] != 0 {
			c1 = formatVCFFloat(float64(float32(nals[3]) / float32(nals[2]+nals[3])))
		} else {
			c1 = "."
		}
		setInfo(v, "FASSOC", c0+","+c1)
	}
	if p.annots&contrastNASSOC != 0 {
		setInfo(v, "NASSOC", fmt.Sprintf("%d,%d,%d,%d", nals[0], nals[1], nals[2], nals[3]))
	}
	if len(caseAlSmpl) > 0 {
		setInfo(v, "NOVELAL", strings.Join(caseAlSmpl, ","))
		p.ncaseAl++
	}
	if len(caseGtSmpl) > 0 {
		setInfo(v, "NOVELGT", strings.Join(caseGtSmpl, ","))
		p.ncaseGt++
	}
	p.ntested++
	return []*vcf.Variant{v}, nil
}

// Destroy writes the one-line summary to stderr, mirroring contrast.c run().
// With -f/--max-allele-freq it also writes the rare-allele enrichment line: the
// region-wide Fisher's exact probability over the pooled minor-allele counts and
// the control/case non-REF fractions.
func (p *contrastPlugin) Destroy() error {
	if p.stderr != nil {
		fmt.Fprintf(p.stderr, "Total/processed/skipped/case_allele/case_gt:\t%d\t%d\t%d\t%d\t%d\n",
			p.ntotal, p.ntested, p.nskipped, p.ncaseAl, p.ncaseGt)
		if p.maxACSet {
			n := &p.enrichNals
			_, _, fisher := mpileupFisherExact(int64(n[0]), int64(n[1]), int64(n[2]), int64(n[3]))
			val1, val2 := 0.0, 0.0
			if n[0]+n[1] != 0 {
				val1 = float64(float32(n[1]) / float32(n[0]+n[1]))
			}
			if n[2]+n[3] != 0 {
				val2 = float64(float32(n[3]) / float32(n[2]+n[3]))
			}
			fmt.Fprintf(p.stderr, "max_AC/PASSOC/FASSOC/NASSOC:\t%d\t%e\t%f,%f\t%d,%d,%d,%d\n",
				p.maxAC, fisher, val1, val2, n[0], n[1], n[2], n[3])
		}
	}
	return nil
}

// resolveSampleGroup resolves a comma-separated sample list or @file into the
// sorted list of VCF sample indices, mirroring read_sample_list_or_file() in
// contrast.c. A leading "@" forces file interpretation; otherwise the string is
// first tried as a sample list and, on a miss, as a file. Upstream sorts the
// indices with a comparator that is effectively order-preserving for distinct
// values, so the user-supplied order is kept (matching observed output).
func resolveSampleGroup(hdr *vcf.Header, str string, force bool) ([]int, error) {
	names, err := readSampleListOrFile(hdr, str)
	if err != nil {
		return nil, err
	}
	idx := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		idx[s] = i
	}
	var out []int
	var skipped int
	for _, n := range names {
		i, ok := idx[n]
		if !ok {
			if force {
				skipped++
				continue
			}
			return nil, fmt.Errorf("contrast: the sample \"%s\" is not present in the VCF", n)
		}
		out = append(out, i)
	}
	if len(out) == 0 && !force {
		return nil, fmt.Errorf("contrast: none of the samples are present in the VCF: %s", str)
	}
	_ = skipped
	return out, nil
}

// readSampleListOrFile returns the sample names for a -0/-1 argument. A "@"
// prefix (or a plain string that is not a known VCF sample but names a readable
// file) is read as a one-name-per-line file; otherwise the string is split on
// commas. This mirrors hts_readlist's list-or-file behaviour used by upstream.
func readSampleListOrFile(hdr *vcf.Header, str string) ([]string, error) {
	known := make(map[string]bool, len(hdr.Samples))
	for _, s := range hdr.Samples {
		known[s] = true
	}
	readFile := func(path string) ([]string, error) {
		data, err := os.ReadFile(path)
		if err != nil {
			return nil, err
		}
		var names []string
		for _, line := range strings.Split(string(data), "\n") {
			line = strings.TrimRight(line, "\r")
			line = strings.TrimSpace(line)
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			// Take the first whitespace-delimited token, matching hts_readlist.
			names = append(names, strings.Fields(line)[0])
		}
		return names, nil
	}
	if strings.HasPrefix(str, "@") {
		return readFile(str[1:])
	}
	// Comma-separated list: keep it as a list unless none of the entries are
	// known samples and the string names a readable file (matching upstream's
	// list-then-file fallback).
	list := strings.Split(str, ",")
	anyKnown := false
	for _, n := range list {
		if known[n] {
			anyKnown = true
			break
		}
	}
	if !anyKnown {
		if names, err := readFile(str); err == nil && len(names) > 0 {
			return names, nil
		}
	}
	return list, nil
}
