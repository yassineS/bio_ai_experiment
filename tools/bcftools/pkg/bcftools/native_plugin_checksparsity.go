// Native port of the upstream `check-sparsity` plugin
// (plugins/check-sparsity.c). It reports samples that lack a sufficient number
// of genotyped markers within a chromosome (the default) and prints, for each
// chromosome, the samples that did not reach the -n threshold. The plugin
// suppresses the VCF/BCF output and writes its "<region>\t<sample>" report to
// stdout. It needs to see the whole stream grouped by chromosome, so it is a
// serial bufferedPlugin.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("check-sparsity", func() NativePlugin { return &checkSparsityPlugin{} })
}

// checkSparsityPlugin implements the `check-sparsity` plugin. In its default
// mode it groups the stream by chromosome and reports the samples that did not
// reach -n genotyped markers per chromosome. With -r/-R it instead processes
// each region independently — overlapping records only — and reports the still
// sparse samples once, labelled by the verbatim region token, exactly as
// upstream's indexed test_region() does (no intra-region chromosome grouping,
// because the index query is per-region). check-sparsity therefore self-applies
// the shared region selection (it is a regionTargetSink) rather than letting the
// framework pre-filter and collapse the per-region labels.
type checkSparsityPlugin struct {
	hdr      *vcf.Header
	minSites int
	out      io.Writer
	// regionSpecs holds the -r/-R region tokens (verbatim label + parsed window)
	// in order. Empty means the default per-chromosome mode.
	regionSpecs []regionSpec
}

// SuppressVCF reports true: `+check-sparsity` emits no VCF/BCF output, only its
// textual report on stdout (upstream's run() prints and returns).
func (p *checkSparsityPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *checkSparsityPlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *checkSparsityPlugin) Name() string { return "check-sparsity" }

// About returns the one-line description, matching check-sparsity.c about().
func (p *checkSparsityPlugin) About() string {
	return "Print samples without genotypes in a region or chromosome"
}

// RunStyle reports that check-sparsity is a run()-style plugin: upstream
// exports a `run` symbol, so its options precede the input file with no `--`
// separator (e.g. `+check-sparsity FILE -n 2`).
func (p *checkSparsityPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of check-sparsity's own flags consumes the
// following CLI token, used by the host to split the input-file positional from
// the plugin options.
func (p *checkSparsityPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-n", "--n-markers", "-r", "--regions", "-R", "--regions-file", "-v", "--verbosity":
		return true
	}
	return false
}

// RegionTargetCaps reports that check-sparsity exposes NEITHER family to the
// shared region/target filter. Its -r/-R are region/chromosome SELECTORS that
// also group and label the report by region, and — unlike every other plugin —
// upstream's check-sparsity reads its -R file with hts_readlist + tbx_itr_querys
// (verbatim region-list strings, colon syntax), NOT the synced reader's TSV
// format. check-sparsity therefore parses -r/-R itself in Init (it has no -t/-T
// option), so the shared filter must leave both families alone.
func (p *checkSparsityPlugin) RegionTargetCaps() regionTargetCaps {
	return regionTargetCaps{regions: false, targets: false}
}

// Init parses -n/--n-markers and rejects the index-based region modes (-r/-R),
// which require tabix/BCF index jumping not available in the native pipeline.
func (p *checkSparsityPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	p.minSites = 1
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("check-sparsity: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-n", "--n-markers":
			v, err := next()
			if err != nil {
				return nil, err
			}
			var n int
			if _, err := fmt.Sscanf(v, "%d", &n); err != nil {
				return nil, fmt.Errorf("check-sparsity: could not parse: -n %s", v)
			}
			p.minSites = n
		case "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		case "-r", "--regions":
			v, err := next()
			if err != nil {
				return nil, err
			}
			specs, perr := parseRegionSpecs(SplitCommaList(v))
			if perr != nil {
				return nil, fmt.Errorf("check-sparsity: %w", perr)
			}
			p.regionSpecs = append(p.regionSpecs, specs...)
		case "-R", "--regions-file":
			v, err := next()
			if err != nil {
				return nil, err
			}
			// Upstream reads -R via hts_readlist + tbx_itr_querys: each line is a
			// verbatim region-list string (colon syntax), NOT the synced reader's
			// TSV. A tab-separated / BED line therefore fails tbx_itr_querys and
			// upstream silently produces no output for it. We reproduce the
			// region-list parse and document the BED quirk in docs/UPSTREAM_BUGS.md.
			lines, perr := loadCheckSparsityRegionFile(v)
			if perr != nil {
				return nil, fmt.Errorf("check-sparsity: %w", perr)
			}
			specs, perr := parseRegionSpecs(lines)
			if perr != nil {
				return nil, fmt.Errorf("check-sparsity: %w", perr)
			}
			p.regionSpecs = append(p.regionSpecs, specs...)
		default:
			return nil, fmt.Errorf("check-sparsity: unsupported option %q", a)
		}
	}
	if !hasFormatHeader(hdr.MetaInfo, "GT") {
		return nil, fmt.Errorf("check-sparsity: GT field is not present")
	}
	return hdr, nil
}

// ProcessAll reports the sparse samples. With -r/-R it processes each region
// independently (records overlapping it) and reports once per region labelled by
// the verbatim region token, mirroring upstream's indexed test_region(). Without
// regions it groups the whole stream by chromosome, mirroring
// test_region(reg==NULL).
func (p *checkSparsityPlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	if len(p.hdr.Samples) == 0 {
		return nil, nil
	}
	if len(p.regionSpecs) > 0 {
		// Per-region: each region token is its own group, labelled by the token,
		// with no intra-region chromosome boundary (the index query never sees a
		// chromosome change). Regions are processed in the order they were given.
		for _, spec := range p.regionSpecs {
			region := []region{spec.region}
			sub := variants[:0:0]
			for _, v := range variants {
				if overlapsAny(v, region) {
					sub = append(sub, v)
				}
			}
			// Upstream's indexed test_region() returns without reporting when the
			// region query yields no records (the itr is empty), so an empty
			// region prints nothing rather than "all samples sparse".
			if len(sub) == 0 {
				continue
			}
			p.scanGroup(sub, spec.label, false)
		}
		return nil, nil
	}

	// Default per-chromosome grouping over the whole stream. This mirrors
	// upstream's single non-indexed test_region(reg==NULL): report at each
	// chromosome boundary, reset per chromosome, and break the ENTIRE scan once
	// every sample has reached min_sites (so later chromosomes are not reported).
	nsmpl := len(p.hdr.Samples)
	reset := func() ([]int, []int) {
		smpl := make([]int, nsmpl)
		for i := range smpl {
			smpl[i] = i
		}
		return smpl, make([]int, nsmpl)
	}
	smpl, nsites := reset()
	report := func(reg string) {
		if p.out != nil {
			for _, s := range smpl {
				fmt.Fprintf(p.out, "%s\t%s\n", reg, p.hdr.Samples[s])
			}
		}
		smpl, nsites = reset()
	}
	curChrom := ""
	haveChrom := false
	nread := false
	for _, v := range variants {
		if haveChrom && v.Chrom != curChrom {
			report(curChrom)
			nread = false
		}
		curChrom = v.Chrom
		haveChrom = true
		if !formatHasTag(v, "GT") {
			continue
		}
		for i := 0; i < len(smpl); i++ {
			gt, ok := sampleGT(v, smpl[i])
			if !ok {
				continue
			}
			if len(gt.alleles) == 0 || gt.alleles[0] == missingAllele {
				continue
			}
			nsites[i]++
			if nsites[i] < p.minSites {
				continue
			}
			smpl = append(smpl[:i], smpl[i+1:]...)
			nsites = append(nsites[:i], nsites[i+1:]...)
			i--
		}
		nread = true
		if len(smpl) == 0 {
			break
		}
	}
	if nread {
		report(curChrom)
	}
	return nil, nil
}

// scanGroup runs the sparsity scan over one region's records and reports the
// still-sparse samples labelled by the region token, mirroring the final
// report() of upstream's indexed test_region() (rid stays -1, so the label is
// the verbatim region string and there is no intra-region chromosome grouping).
func (p *checkSparsityPlugin) scanGroup(variants []*vcf.Variant, label string, _ bool) {
	nsmpl := len(p.hdr.Samples)
	smpl := make([]int, nsmpl)
	for i := range smpl {
		smpl[i] = i
	}
	nsites := make([]int, nsmpl)
	for _, v := range variants {
		if !formatHasTag(v, "GT") {
			continue
		}
		for i := 0; i < len(smpl); i++ {
			gt, ok := sampleGT(v, smpl[i])
			if !ok {
				continue
			}
			if len(gt.alleles) == 0 || gt.alleles[0] == missingAllele {
				continue
			}
			nsites[i]++
			if nsites[i] < p.minSites {
				continue
			}
			smpl = append(smpl[:i], smpl[i+1:]...)
			nsites = append(nsites[:i], nsites[i+1:]...)
			i--
		}
		if len(smpl) == 0 {
			break
		}
	}
	if p.out != nil {
		for _, s := range smpl {
			fmt.Fprintf(p.out, "%s\t%s\n", label, p.hdr.Samples[s])
		}
	}
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *checkSparsityPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return nil, nil
}

// Destroy releases resources (none held).
func (p *checkSparsityPlugin) Destroy() error { return nil }

// loadCheckSparsityRegionFile reads a -R file the way upstream check-sparsity
// does: hts_readlist returns each non-comment line verbatim, and that string is
// handed to tbx_itr_querys. Only the colon region-list syntax ("chr",
// "chr:beg-end") parses there; a tab-separated / BED line fails to parse and
// upstream silently produces no output for it. We reproduce that by keeping only
// the single-column lines (a verbatim region token) and dropping multi-column
// (TSV/BED) lines — see docs/UPSTREAM_BUGS.md for the reproducer.
func loadCheckSparsityRegionFile(path string) ([]string, error) {
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var out []string
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		// A line with whitespace/tabs is a TSV/BED line that upstream's
		// tbx_itr_querys cannot parse: silently drop it (upstream prints nothing).
		if strings.IndexFunc(line, func(r rune) bool { return r == ' ' || r == '\t' }) >= 0 {
			continue
		}
		out = append(out, line)
	}
	return out, sc.Err()
}
