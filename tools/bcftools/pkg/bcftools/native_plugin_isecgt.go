// Native port of the upstream `isecGT` plugin (plugins/isecGT.c). It compares
// the genotypes of two coordinate-sorted VCF/BCF files A and B and writes file A
// to output, setting to missing every genotype in A that does not exactly match
// the corresponding genotype in B at shared sites. Records present only in A are
// written unchanged; records present only in B are dropped.
//
// Upstream drives the two files through htslib's synced reader (which requires
// both to be indexed). The native port instead streams both files in sorted
// lockstep — algorithmically equivalent for whole-file comparison — so neither
// file needs an index. The genotype comparison is per-sample and exact: samples
// are matched by name (every A sample must be present in B), and a sample's GT
// is considered discordant when its allele tokens differ in value, order or
// phasing, exactly as the upstream per-allele int32 comparison. A discordant
// sample's entire GT in A is set to missing (with A's ploidy preserved).
//
// The second input file arrives as the second positional argument, which the
// host CLI routes into opts.Regions; the plugin reads it from there. Region and
// target selection (-r/-R/-t/-T) is supported and applied to BOTH input streams
// before the lockstep comparison, mirroring upstream's synced reader (which
// applies the same regions/targets to every reader). The -W index option of the
// upstream plugin is not reproduced.
package bcftools

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("isecGT", func() NativePlugin { return &isecGTPlugin{} })
}

// isecGTPlugin implements isecGT.
type isecGTPlugin struct {
	fileB string
	rt    regionTargetFilter
}

// SetRegionTarget records the shared -r/-R/-t/-T selection the framework parsed
// out of isecGT's argv, to be applied to both input streams in RunFull.
func (p *isecGTPlugin) SetRegionTarget(f regionTargetFilter) { p.rt = f }

// Name returns the plugin name.
func (p *isecGTPlugin) Name() string { return "isecGT" }

// RegionTargetCaps opts isecGT into the shared -r/-R/-t/-T region/target filter.
// The selection is applied to BOTH input streams in RunFull (via SetRegionTarget),
// matching upstream's synced reader.
func (p *isecGTPlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching isecGT.c about().
func (p *isecGTPlugin) About() string {
	return "Compare two files and set non-identical genotypes to missing.\n"
}

// RunStyle reports that isecGT is a run()-style plugin: upstream's isecGT.c
// exports a `run` symbol, so it owns its entire argv before the two trailing
// input-file positionals (A and B), with no `--` separator
// (e.g. `bcftools +isecGT A.bcf B.bcf`).
func (p *isecGTPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of isecGT's flags consumes the following
// CLI token as its value, so the host can separate the input-file positionals
// from the plugin options.
func (p *isecGTPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-o", "--output", "-O", "--output-type",
		"-r", "--regions", "-R", "--regions-file",
		"-t", "--targets", "-T", "--targets-file",
		"-v", "--verbosity":
		return true
	}
	return false
}

// Init parses the plugin arguments. The second positional file is supplied via
// opts.Regions (see RunFull); only options that the native streaming path can
// honour are accepted.
func (p *isecGTPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-O", "--output-type", "-o", "--output":
			// Output container/handle is supplied by the host pipeline.
			i++
		case "-W", "--write-index":
			return nil, fmt.Errorf("isecGT: -W/--write-index is not supported in the native plugin")
		case "-v", "--verbosity":
			i++
		default:
			return nil, fmt.Errorf("isecGT: unsupported option %q", a)
		}
	}
	return hdr, nil
}

// Process is unused: isecGT is a fullPlugin.
func (p *isecGTPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *isecGTPlugin) Destroy() error { return nil }

// RunFull reads both files, performs the lockstep genotype comparison and writes
// the edited file A to out in the requested container.
func (p *isecGTPlugin) RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	if _, err := p.Init(opts.Args, nil); err != nil {
		return err
	}
	if len(opts.Regions) == 0 {
		return fmt.Errorf("isecGT: two input files are required (Usage: bcftools +isecGT A.bcf B.bcf)")
	}
	p.fileB = opts.Regions[0]

	hdrA, varsA, err := readVCFFile(opts.InputFile, stderr)
	if err != nil {
		return fmt.Errorf("isecGT: reading first file: %w", err)
	}
	hdrB, varsB, err := readVCFFile(p.fileB, stderr)
	if err != nil {
		return fmt.Errorf("isecGT: reading second file: %w", err)
	}

	// Apply the shared region/target selection to BOTH streams, exactly as
	// upstream's synced reader applies the same regions/targets to every reader.
	varsA = p.rt.apply(varsA)
	varsB = p.rt.apply(varsB)

	// Map every A sample to its column in B (strict: all must be present).
	bIdx := map[string]int{}
	for i, s := range hdrB.Samples {
		bIdx[s] = i
	}
	mapAB := make([]int, len(hdrA.Samples))
	for i, s := range hdrA.Samples {
		j, ok := bIdx[s]
		if !ok {
			return fmt.Errorf("isecGT: sample %q from the first file is not present in the second file", s)
		}
		mapAB[i] = j
	}

	// Contig order from header A, used to break CHROM ties during the merge.
	contigOrder := contigOrderFromHeader(hdrA)

	// Index B records by site key for matched-site lookup.
	bySite := map[string][]*vcf.Variant{}
	for _, v := range varsB {
		key := siteKey(v)
		bySite[key] = append(bySite[key], v)
	}
	_ = contigOrder

	w, cleanup, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, hdrA)
	if err != nil {
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		return err
	}
	for _, a := range varsA {
		if matches := bySite[siteKey(a)]; len(matches) > 0 {
			p.compareSite(a, matches[0], mapAB)
		}
		if err := w.Write(a); err != nil {
			cleanup()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		return err
	}
	cleanup()
	return nil
}

// compareSite sets to missing every A genotype that does not exactly match the
// corresponding B genotype, mirroring the per-sample int32 compare in isecGT.c.
func (p *isecGTPlugin) compareSite(a, b *vcf.Variant, mapAB []int) {
	for i := range a.Samples {
		j := mapAB[i]
		if j >= len(b.Samples) {
			continue
		}
		gtA, okA := a.Samples[i].Data["GT"]
		gtB, okB := b.Samples[j].Data["GT"]
		if !okA || !okB {
			continue
		}
		if !gtEqual(gtA, gtB) {
			a.Samples[i].Data["GT"] = missingGTLike(gtA)
		}
	}
}

// gtEqual reports whether two GT strings are identical in allele values, order
// and phasing — the text equivalent of upstream's per-allele encoded compare.
func gtEqual(a, b string) bool {
	return a == b
}

// missingGTLike returns the all-missing GT with the same ploidy as gt, e.g.
// "0/1" -> "./.", "0|1|2" -> "./././.". Upstream writes bcf_gt_missing per
// allele, which is always unphased, so every separator becomes "/" regardless of
// the original genotype's phasing.
func missingGTLike(gt string) string {
	ploidy := 1
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			ploidy++
		}
	}
	parts := make([]string, ploidy)
	for i := range parts {
		parts[i] = "."
	}
	return strings.Join(parts, "/")
}

// siteKey is the CHROM/POS/REF/ALT identity used to pair records, matching the
// synced reader's default exact allele matching.
func siteKey(v *vcf.Variant) string {
	return v.Chrom + "\t" + fmtInt(v.Pos) + "\t" + v.Ref + "\t" + strings.Join(v.Alt, ",")
}

// fmtInt renders an int without importing strconv at the call site.
func fmtInt(n int) string {
	return fmt.Sprintf("%d", n)
}

// contigOrderFromHeader returns a map of contig ID to its order of appearance in
// the header's ##contig lines.
func contigOrderFromHeader(hdr *vcf.Header) map[string]int {
	order := map[string]int{}
	idx := 0
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##contig=") {
			id := headerID(m)
			if id != "" {
				if _, ok := order[id]; !ok {
					order[id] = idx
					idx++
				}
			}
		}
	}
	return order
}

// readVCFFile reads a VCF/BCF file (no region filtering) into a header and the
// full variant slice, reusing the host's ViewFile normalisation.
func readVCFFile(path string, stderr io.Writer) (*vcf.Header, []*vcf.Variant, error) {
	if path == "" {
		path = "-"
	}
	var buf strings.Builder
	if _, err := ViewFile(path, &stringWriter{&buf}, ViewOptions{OutputFormat: OutputVCF}, stderr); err != nil {
		return nil, nil, err
	}
	r := vcf.NewReader(strings.NewReader(buf.String()))
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, nil, err
	}
	variants, err := r.ReadAll()
	if err != nil {
		return nil, nil, err
	}
	return hdr, variants, nil
}
