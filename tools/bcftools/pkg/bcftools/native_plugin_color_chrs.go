// Native port of the upstream `color-chrs` plugin (plugins/color-chrs.c). It
// colors shared chromosomal segments in a trio (mother, father, child) or a
// pair of unrelated samples by running an HMM Viterbi decode over the phased
// genotypes, writing the segmentation and per-sample switch counts to a
// separate `<prefix>.dat` file (NOT stdout).
//
// The model is a fixed-state HMM whose transition matrix is built once from a
// crossover probability pij=2e-8 (plus a genotype-error term pgt_err=1e-9), then
// pre-multiplied to up to 10000 positions and matrix-power-jumped for larger
// gaps — exactly htslib/HMM.c's hmm_init / hmm_set_tprob / _set_tprob /
// hmm_run_viterbi. Every operation on this path is IEEE-754 +,-,*,/ (the
// emission probabilities are products of the constant pgt_err / 1-pgt_err and,
// for the unrelated mode, the fixed af=0.5 terms); there is NO transcendental
// in the decode, so the Go port is byte-reproducible against the C HMM. The
// `.dat` output therefore matches upstream exactly.
//
// color-chrs is a generic init/process plugin (its options follow `--`). It is
// registered as a multiOutputPlugin so the host hands it the whole invocation
// and it owns its single `<prefix>.dat` writer.
package bcftools

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("color-chrs", func() NativePlugin { return &colorChrsPlugin{} })
}

// color-chrs mode constants, matching C_TRIO / C_UNRL in color-chrs.c.
const (
	ccTrio = 1
	ccUnrl = 2
)

// Trio state indices, matching the TRIO_* defines in color-chrs.c.
const (
	trioAC = 0
	trioAD = 1
	trioBC = 2
	trioBD = 3
	trioCA = 4
	trioDA = 5
	trioCB = 6
	trioDB = 7
)

// Unrelated state indices, matching the UNRL_* defines in color-chrs.c.
const (
	unrlXXXX = 0
	unrl0x0x = 1
	unrl0xx0 = 2
	unrlX00x = 3
	unrlX0x0 = 4
	unrl0101 = 5
	unrl0110 = 6
)

// Switch flags, matching SW_MOTHER / SW_FATHER in color-chrs.c.
const (
	swMother = 1
	swFather = 2
)

// colorChrsPlugin implements the `color-chrs` plugin end to end.
type colorChrsPlugin struct{}

// Name returns the plugin name.
func (p *colorChrsPlugin) Name() string { return "color-chrs" }

// About returns the one-line description, matching color-chrs.c about().
func (p *colorChrsPlugin) About() string {
	return "Color shared chromosomal segments, requires phased GTs.\n"
}

// RunMulti executes color-chrs: parse options, read the input, run the per-
// chromosome Viterbi decode, and write the `<prefix>.dat` file, matching run()
// in color-chrs.c.
func (p *colorChrsPlugin) RunMulti(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	var prefix, trioSamples, unrelatedSamples string
	args := opts.Args
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("color-chrs: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-p", "--prefix":
			v, err := next()
			if err != nil {
				return err
			}
			prefix = v
		case "-t", "--trio":
			v, err := next()
			if err != nil {
				return err
			}
			trioSamples = v
		case "-u", "--unrelated":
			v, err := next()
			if err != nil {
				return err
			}
			unrelatedSamples = v
		default:
			return fmt.Errorf("color-chrs: unsupported option %q", a)
		}
	}
	if trioSamples != "" && unrelatedSamples != "" {
		return fmt.Errorf("color-chrs: expected only one of the -t/-u options")
	}
	if trioSamples == "" && unrelatedSamples == "" {
		return fmt.Errorf("color-chrs: expected one of the -t/-u options")
	}
	if prefix == "" {
		return fmt.Errorf("color-chrs: expected the -p option")
	}

	hdr, variants, err := readPluginInput(opts, stderr)
	if err != nil {
		return fmt.Errorf("color-chrs: %w", err)
	}
	idx := sampleIndex(hdr)

	st := &colorChrsState{hdr: hdr, pij: 2e-8, pgtErr: 1e-9, prevRID: ""}
	if trioSamples != "" {
		names := strings.Split(trioSamples, ",")
		if len(names) != 3 {
			return fmt.Errorf("color-chrs: expected three sample names with -t")
		}
		for i, name := range names {
			j, ok := idx[name]
			if !ok {
				return fmt.Errorf("color-chrs: %d-th sample not found: %s", i+1, name)
			}
			switch i {
			case 0:
				st.imother = j
			case 1:
				st.ifather = j
			case 2:
				st.ichild = j
			}
		}
		st.mode = ccTrio
		st.initHMMTrio()
	} else {
		names := strings.Split(unrelatedSamples, ",")
		if len(names) != 2 {
			return fmt.Errorf("color-chrs: expected two sample names with -u")
		}
		for i, name := range names {
			j, ok := idx[name]
			if !ok {
				return fmt.Errorf("color-chrs: %d-th sample not found: %s", i+1, name)
			}
			if i == 0 {
				st.isample = j
			} else {
				st.jsample = j
			}
		}
		st.mode = ccUnrl
		st.initHMMUnrelated()
	}

	// Open the .dat output. The header lines are written lazily on the first
	// flush (matching the !args->fp guard in flush_viterbi).
	f, err := os.Create(prefix + ".dat")
	if err != nil {
		return fmt.Errorf("color-chrs: %s.dat: %w", prefix, err)
	}
	defer f.Close()
	st.fp = f

	// process(): on a new chromosome, flush the previous one; accumulate the
	// observed probabilities for the current site.
	for _, v := range variants {
		if st.prevRID == "" {
			st.prevRID = v.Chrom
		}
		if st.prevRID != v.Chrom {
			st.flushViterbi()
		}
		st.prevRID = v.Chrom
		if st.mode == ccTrio {
			st.setObservedProbTrio(v)
		} else {
			st.setObservedProbUnrelated(v)
		}
	}
	// destroy(): flush the final chromosome.
	st.flushViterbi()
	return nil
}

// Init satisfies NativePlugin; the real work runs in RunMulti.
func (p *colorChrsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return hdr, nil
}

// Process satisfies NativePlugin; never reached (RunMulti owns the run).
func (p *colorChrsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy satisfies NativePlugin.
func (p *colorChrsPlugin) Destroy() error { return nil }
