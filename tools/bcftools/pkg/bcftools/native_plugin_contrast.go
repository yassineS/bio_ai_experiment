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
// The rare-allele enrichment mode (-f/--max-allele-freq) and the index/region
// jump options require htslib machinery the native pipeline does not provide
// and remain unsupported.
package bcftools

import (
	"fmt"
	"io"
	"os"
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

	filter *pluginFilter // compiled -i/-e site-level pre-filter, nil if none

	stderr io.Writer
}

// Name returns the plugin name.
func (p *contrastPlugin) Name() string { return "contrast" }

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
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("contrast: %s requires an argument", a)
			}
			i++
			return args[i], nil
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
			return nil, fmt.Errorf("contrast: the rare-allele enrichment mode (-f) is not supported by the native plugin")
		case "-r", "--regions", "-R", "--regions-file", "-t", "--targets", "-T", "--targets-file":
			return nil, fmt.Errorf("contrast: region/target selection (%s) is not supported by the native plugin", a)
		case "--regions-overlap", "--targets-overlap":
			return nil, fmt.Errorf("contrast: %s is not supported by the native plugin", a)
		case "-o", "--output":
			return nil, fmt.Errorf("contrast: writing to a file (-o) is not supported by the native plugin; use stdout")
		case "-O", "--output-type":
			if _, err := next(); err != nil {
				return nil, err
			}
			return nil, fmt.Errorf("contrast: -O is handled by the host -O; not a plugin option here")
		case "-W", "--write-index":
			return nil, fmt.Errorf("contrast: --write-index is not supported by the native plugin")
		case "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		default:
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
		f, ferr := newPluginFilter(filterExpr, filterExclude)
		if ferr != nil {
			return nil, fmt.Errorf("contrast: %w", ferr)
		}
		p.filter = f
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
func (p *contrastPlugin) Destroy() error {
	if p.stderr != nil {
		fmt.Fprintf(p.stderr, "Total/processed/skipped/case_allele/case_gt:\t%d\t%d\t%d\t%d\t%d\n",
			p.ntotal, p.ntested, p.nskipped, p.ncaseAl, p.ncaseGt)
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
