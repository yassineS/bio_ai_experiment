// Native port of the upstream `check-sparsity` plugin
// (plugins/check-sparsity.c). It reports samples that lack a sufficient number
// of genotyped markers within a chromosome (the default) and prints, for each
// chromosome, the samples that did not reach the -n threshold. The plugin
// suppresses the VCF/BCF output and writes its "<region>\t<sample>" report to
// stdout. It needs to see the whole stream grouped by chromosome, so it is a
// serial bufferedPlugin.
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("check-sparsity", func() NativePlugin { return &checkSparsityPlugin{} })
}

// checkSparsityPlugin implements the `check-sparsity` plugin in its default
// per-chromosome mode. The index-based -r/-R region modes are reported as
// unsupported because they depend on htslib's tabix/BCF index region jumping,
// which the native pipeline does not expose to plugins.
type checkSparsityPlugin struct {
	hdr      *vcf.Header
	minSites int
	out      io.Writer
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
		case "-r", "--regions", "-R", "--regions-file":
			return nil, fmt.Errorf("check-sparsity: the index-based region modes (-r/-R) are not supported by the native plugin; run upstream for indexed region jumping")
		default:
			return nil, fmt.Errorf("check-sparsity: unsupported option %q", a)
		}
	}
	if !hasFormatHeader(hdr.MetaInfo, "GT") {
		return nil, fmt.Errorf("check-sparsity: GT field is not present")
	}
	return hdr, nil
}

// ProcessAll streams the records grouped by chromosome and prints, at each
// chromosome boundary (and at EOF), the samples that did not accumulate
// min_sites genotyped markers, mirroring test_region(reg==NULL) and report().
func (p *checkSparsityPlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	nsmpl := len(p.hdr.Samples)
	if nsmpl == 0 {
		return nil, nil
	}

	// smpl holds the indices still being tracked (not yet reaching min_sites),
	// nsites their genotyped-marker counts; both are reset per chromosome.
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

	// rid==-1 mirrors upstream's "no record seen yet"; nread tracks whether at
	// least one GT-bearing record was processed since the last report, which
	// gates the final report().
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

		// Skip records with no GT data (upstream `continue` on missing fmt_gt),
		// which leaves nread unchanged for those sites.
		if !formatHasTag(v, "GT") {
			continue
		}

		// Update the tracked samples: a sample with a present (non-missing)
		// genotype at this site advances its counter; once it reaches min_sites
		// it is dropped from tracking. Samples whose GT is missing are skipped.
		for i := 0; i < len(smpl); i++ {
			gt, ok := sampleGT(v, smpl[i])
			if !ok {
				continue
			}
			// ial==0 means the first allele is missing/end => treated as missing.
			if len(gt.alleles) == 0 || gt.alleles[0] == missingAllele {
				continue
			}
			nsites[i]++
			if nsites[i] < p.minSites {
				continue
			}
			// Remove sample i from the tracking lists.
			smpl = append(smpl[:i], smpl[i+1:]...)
			nsites = append(nsites[:i], nsites[i+1:]...)
			i--
		}
		nread = true
		if len(smpl) == 0 {
			// Upstream breaks the entire read loop once every sample reached
			// min_sites; no further chromosomes are processed or reported.
			break
		}
	}
	if nread {
		report(curChrom)
	}
	return nil, nil
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *checkSparsityPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return nil, nil
}

// Destroy releases resources (none held).
func (p *checkSparsityPlugin) Destroy() error { return nil }
