package bcftools

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ConvertOptions controls the behaviour of Convert / ConvertFile.
//
// The v1 port covers the common-case "round-trip VCF/BCF through a
// different container format with optional sample / region filtering"
// pipeline plus the GEN/HAP/HAPLEGEND/TSV/gVCF modes wired in
// subcmds_convert.go. Upstream's `vcfconvert.c` advertises PLINK
// (`--plink`/`--tped`) options too, but those lines are commented out
// upstream with no implementation, so there is nothing to port there.
type ConvertOptions struct {
	// OutputFormat is the requested container (-O v|z|b|u).
	OutputFormat OutputFormat
	// CompressLevel forwards the -l gzip level when OutputFormat is OutputVCFGz.
	CompressLevel int
	// Threads is upstream's -@/--threads. When greater than 1 it enables
	// parallel BGZF compression of -O z and -O b output via bgzf.MultiWriter;
	// the framed result decodes byte-identically regardless of thread count.
	Threads int

	// Samples / SamplesFile restrict the per-sample columns to the named
	// set (in the requested order). Missing names are skipped silently
	// unless ForceSamples is false (the default), in which case a missing
	// requested sample is reported via an error.
	Samples      []string
	SamplesFile  string
	ForceSamples bool

	// Regions / RegionsFile apply CHROM[:beg-end] filters at the record
	// level. v1 always treats them as a post-filter (no index seek).
	Regions     []string
	RegionsFile string
	// Targets / TargetsFile are identical to Regions in v1; upstream
	// distinguishes them by whether they imply an indexed seek.
	Targets     []string
	TargetsFile string

	// IncludeExpr / ExcludeExpr are the standard bcftools filter
	// expressions (-i / -e). They re-use the Filter type from view.go.
	IncludeExpr string
	ExcludeExpr string
}

// Convert streams a VCF/BCF source through opts and writes the requested
// format to out. The (path, reader) pair is split from ConvertFile so the
// streaming entry point can be used in tests without a sibling file.
func Convert(in io.Reader, out io.Writer, opts ConvertOptions) (int, error) {
	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return 0, err
	}
	return writeConverted(hdr, variants, out, opts)
}

// ConvertFile opens the named input through iohelper (transparent gzip
// + BCF auto-detect) and emits the converted stream to out.
func ConvertFile(path string, out io.Writer, opts ConvertOptions) (int, error) {
	// Apply samples-file before reading so a bad path fails fast.
	if opts.SamplesFile != "" {
		names, err := LoadSamplesFile(opts.SamplesFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Samples = append(opts.Samples, names...)
	}
	if opts.RegionsFile != "" {
		regs, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Regions = append(opts.Regions, regs...)
	}
	if opts.TargetsFile != "" {
		regs, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Targets = append(opts.Targets, regs...)
	}

	r, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: open %s: %w", path, err)
	}
	defer r.Close()
	return Convert(r, out, opts)
}

// writeConverted is shared between Convert (in-memory pipeline) and any
// future streaming entry point. It applies the sample / region /
// expression filters and emits the result via the unified openOutput
// helper from view.go.
func writeConverted(hdr *vcf.Header, variants []*vcf.Variant, out io.Writer, opts ConvertOptions) (int, error) {
	// 1) Sample restriction. We validate the requested samples against
	//    the input header before narrowing it.
	if len(opts.Samples) > 0 {
		missing := missingSamples(hdr.Samples, opts.Samples)
		if len(missing) > 0 && !opts.ForceSamples {
			return 0, fmt.Errorf("bcftools convert: requested samples missing from input: %s (use --force-samples to ignore)", strings.Join(missing, ", "))
		}
	}
	if len(opts.Samples) > 0 {
		hdr = filterHeaderSamples(hdr, opts.Samples)
	}

	// 2) Build the region post-filter list (Targets + Regions both act
	//    as post-filters in v1).
	postFilters := append([]string{}, opts.Targets...)
	postFilters = append(postFilters, opts.Regions...)
	parsedTargets, err := parseRegions(postFilters)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}

	// 3) Filter expressions.
	include, exclude, err := compileExpressions(ViewOptions{
		IncludeExpr: opts.IncludeExpr,
		ExcludeExpr: opts.ExcludeExpr,
	})
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}

	// 4) Open writer in the requested format.
	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, hdr)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	count := 0
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
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// missingSamples returns the entries in requested that do not appear in
// have. Used for the --force-samples gate.
func missingSamples(have []string, requested []string) []string {
	if len(requested) == 0 {
		return nil
	}
	seen := make(map[string]bool, len(have))
	for _, s := range have {
		seen[s] = true
	}
	var miss []string
	for _, w := range requested {
		if !seen[w] {
			miss = append(miss, w)
		}
	}
	return miss
}
