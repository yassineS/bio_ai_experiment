// Native registration for the upstream `vrfs` plugin (plugins/vrfs.c). vrfs
// (variant read frequency score) assesses site noisiness by pileup of a large
// set of BAM/CRAM alignments against a FASTA reference: for every site it walks
// the reads via htslib's mpileup2 engine (BAQ/realignment, MIN_MQ/MAX_BQ
// tuning, legacy pileup buffers), bins the per-sample variant-allele
// frequencies into a regidx-indexed histogram, and derives a per-site score
// from the across-sample variance. Reproducing this requires the mpileup2
// pileup machinery, BAM/CRAM read access, and the regidx region cursor — none
// of which the native (stdlib-only, VCF-record) plugin framework exposes.
//
// Rather than silently diverge, the plugin is registered but Init returns a
// clean "unsupported" error explaining the reason.
package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("vrfs", func() NativePlugin { return &vrfsPlugin{} }) }

// vrfsPlugin is the unsupported native entry for `vrfs`.
type vrfsPlugin struct{}

// Name returns the plugin name.
func (p *vrfsPlugin) Name() string { return "vrfs" }

// About returns the one-line description, matching vrfs.c about().
func (p *vrfsPlugin) About() string {
	return "Localised assessment of sequencing artefacts, estimate site noisiness (variant read frequency score)"
}

// RunStyle reports that vrfs is a run()-style plugin: upstream exports a `run`
// symbol, so its options precede any input with no `--` separator.
func (p *vrfsPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of vrfs's own flags consumes the following
// CLI token as its value. Most vrfs options take a value; -i/--use-index is a
// boolean toggle and the verbosity flag is optional-argument.
func (p *vrfsPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-a", "--alns", "-f", "--fasta-ref", "-s", "--sites",
		"-d", "--min-depth", "-n", "--nbins", "-r", "--recalc",
		"-b", "--batch", "-m", "--merge-batches", "-M", "--merge-files",
		"-o", "--output", "-O", "--output-type":
		return true
	}
	return false
}

// Init reports that the native port is unsupported and why.
func (p *vrfsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return nil, fmt.Errorf("vrfs: the native plugin is not supported; it requires htslib's mpileup2 pileup engine (BAQ/realignment, legacy pileup buffers), BAM/CRAM read access, and the regidx region cursor to build per-site variant-read-frequency histograms from many alignments, none of which are available in the native plugin framework")
}

// Process is never reached (Init fails first).
func (p *vrfsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *vrfsPlugin) Destroy() error { return nil }
