// Native port of the upstream `variantkey-hex` plugin
// (plugins/variantkey-hex.c). It generates three unsorted VariantKey lookup
// table files into a directory (the first plugin positional, default "./"):
//
//   - vkrs.unsorted.hex : VariantKey(16 hex) \t rsID(8 hex)
//   - rsvk.unsorted.hex : rsID(8 hex) \t VariantKey(16 hex)
//   - nrvk.unsorted.tsv : VariantKey(16 hex) \t REF \t ALT  (non-reversible keys only)
//
// and prints the variant / non-reversible counts to stdout. The VariantKey is
// computed from CHROM, 0-based POS, REF and the first ALT (as in
// add-variantkey). This is a multiOutputPlugin: it owns the three file writers.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("variantkey-hex", func() NativePlugin { return &variantKeyHexPlugin{} })
}

// variantKeyHexPlugin implements variantkey-hex.
type variantKeyHexPlugin struct {
	dir string
}

// Name returns the plugin name.
func (p *variantKeyHexPlugin) Name() string { return "variantkey-hex" }

// About returns the one-line description, matching variantkey-hex.c about().
func (p *variantKeyHexPlugin) About() string {
	return "Generate VariantKey index files\n"
}

// Init records the output directory (the first plugin positional, default
// "./"), matching the C init() which reads argv[1].
func (p *variantKeyHexPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.dir = "./"
	for _, a := range args {
		if a == "" {
			continue
		}
		if a[0] == '-' {
			return nil, fmt.Errorf("variantkey-hex: unsupported option %q", a)
		}
		p.dir = a
		break
	}
	return hdr, nil
}

// Process is unused: variantkey-hex is a multiOutputPlugin.
func (p *variantKeyHexPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *variantKeyHexPlugin) Destroy() error { return nil }

// RunMulti reads the input, writes the three index files, and prints the counts
// to out, mirroring process()/destroy().
func (p *variantKeyHexPlugin) RunMulti(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	if _, err := p.Init(opts.Args, nil); err != nil {
		return err
	}
	hdr, variants, err := readPluginInput(opts, stderr)
	if err != nil {
		return err
	}
	_ = hdr

	vkrsPath := filepath.Join(p.dir, "vkrs.unsorted.hex")
	rsvkPath := filepath.Join(p.dir, "rsvk.unsorted.hex")
	nrvkPath := filepath.Join(p.dir, "nrvk.unsorted.tsv")

	fVkrs, err := os.Create(vkrsPath)
	if err != nil {
		return fmt.Errorf("variantkey-hex: %s: %w", vkrsPath, err)
	}
	defer fVkrs.Close()
	fRsvk, err := os.Create(rsvkPath)
	if err != nil {
		return fmt.Errorf("variantkey-hex: %s: %w", rsvkPath, err)
	}
	defer fRsvk.Close()
	fNrvk, err := os.Create(nrvkPath)
	if err != nil {
		return fmt.Errorf("variantkey-hex: %s: %w", nrvkPath, err)
	}
	defer fNrvk.Close()

	wVkrs := bufio.NewWriter(fVkrs)
	wRsvk := bufio.NewWriter(fRsvk)
	wNrvk := bufio.NewWriter(fNrvk)

	var numvar, nrv uint64
	for _, v := range variants {
		alt := ""
		if len(v.Alt) > 0 {
			alt = v.Alt[0]
		}
		vk := variantKey(v.Chrom, uint32(v.Pos-1), v.Ref, alt)
		rs := parseRSID(v.ID)
		vkx := variantKeyHex(vk)
		rsx := rsHex(rs)
		fmt.Fprintf(wVkrs, "%s\t%s\n", vkx, rsx)
		fmt.Fprintf(wRsvk, "%s\t%s\n", rsx, vkx)
		if vk&1 != 0 {
			fmt.Fprintf(wNrvk, "%s\t%s\t%s\n", vkx, v.Ref, alt)
			nrv++
		}
		numvar++
	}
	if err := wVkrs.Flush(); err != nil {
		return err
	}
	if err := wRsvk.Flush(); err != nil {
		return err
	}
	if err := wNrvk.Flush(); err != nil {
		return err
	}

	// destroy() prints the two count lines to stdout.
	fmt.Fprintf(out, "VariantKeys: %d\n", numvar)
	fmt.Fprintf(out, "Non-reversible VariantKeys: %d\n", nrv)
	return nil
}

// parseRSID reproduces upstream's rsID extraction: drop the first two characters
// of the ID (the "rs" prefix) and parse the leading decimal digits as a uint32.
func parseRSID(id string) uint32 {
	var rs uint32
	if len(id) > 2 {
		s := id[2:]
		for i := 0; i < len(s) && s[i] >= '0' && s[i] <= '9'; i++ {
			rs = rs*10 + uint32(s[i]-'0')
		}
	}
	return rs
}
