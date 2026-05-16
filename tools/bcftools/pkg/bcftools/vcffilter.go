package bcftools

import (
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// FilterMode controls how the -s/--soft-filter name is folded into
// each record's FILTER column. It mirrors the upstream `-m +|x` flag.
//
//	FilterModeReplace  ("" default) — overwrite the FILTER column with the
//	                                  soft-filter name when the record fails.
//	FilterModeAdd      ("+")        — append the soft-filter name (preserving
//	                                  any pre-existing names).
//	FilterModeReset    ("x")        — clear FILTER (set to PASS) at sites that
//	                                  pass the include/exclude expression.
//
// Both `+` and `x` may be combined (`+x`), in which case the failing-record
// behaviour is "add" and the passing-record behaviour is "reset".
type FilterMode int

// FilterMode values; see the type doc for the upstream `-m +|x` mapping.
const (
	FilterModeReplace FilterMode = 0
	FilterModeAdd     FilterMode = 1 << iota
	FilterModeReset
)

// ParseFilterMode parses the upstream -m/--mode value. The empty string
// yields FilterModeReplace (the upstream default). Any combination of '+'
// and 'x' is accepted; unknown characters return an error.
func ParseFilterMode(s string) (FilterMode, error) {
	mode := FilterModeReplace
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case '+':
			mode |= FilterModeAdd
		case 'x':
			mode |= FilterModeReset
		default:
			return 0, fmt.Errorf("bcftools filter: bad --mode %q (allowed: +, x, +x)", s)
		}
	}
	return mode, nil
}

// SetGTsMode controls how -S/--set-GTs rewrites the GT field of samples
// belonging to a failing record. Upstream accepts the literal strings "."
// (missing) or "0" (homozygous reference).
type SetGTsMode int

// SetGTsMode values; see the type doc for the upstream `-S .|0` mapping.
const (
	SetGTsOff     SetGTsMode = 0
	SetGTsMissing SetGTsMode = 1
	SetGTsRef     SetGTsMode = 2
)

// ParseSetGTsMode parses the upstream -S/--set-GTs value. The empty string
// yields SetGTsOff (the upstream default).
func ParseSetGTsMode(s string) (SetGTsMode, error) {
	switch s {
	case "":
		return SetGTsOff, nil
	case ".":
		return SetGTsMissing, nil
	case "0":
		return SetGTsRef, nil
	}
	return 0, fmt.Errorf("bcftools filter: bad --set-GTs %q (allowed: . or 0)", s)
}

// VCFFilterOptions controls the behaviour of VCFFilter / VCFFilterFile.
// The v1 port covers the common-case "soft-filter records by include /
// exclude expression and optionally rewrite GT on failing samples" path.
//
// SnpGap and IndelGap mirror the upstream `-g/--SnpGap` and `-G/--IndelGap`
// flags: failing variants are tagged with the soft-filter name when a
// neighbouring indel (for SnpGap) or another indel cluster member (for
// IndelGap) is within INT base pairs on the same contig.
type VCFFilterOptions struct {
	// OutputFormat / CompressLevel forward to the shared openOutput helper.
	OutputFormat  OutputFormat
	CompressLevel int

	// IncludeExpr / ExcludeExpr drive the soft-filter decision. A record
	// passes the test when (include matches) OR (exclude does not match).
	// Both empty means every record passes; in that case the soft-filter
	// is never set.
	IncludeExpr string
	ExcludeExpr string

	// SoftFilter is the upstream -s/--soft-filter NAME. The literal "+"
	// asks the library to invent a unique name ("FilterN") per call;
	// matching upstream's behaviour.
	SoftFilter string

	// Mode is the combination of `+` (add to existing FILTER) and `x`
	// (reset FILTER on passing records). Default is replace-on-fail.
	Mode FilterMode

	// SetGTs rewrites the GT entry of every sample in a failing record.
	SetGTs SetGTsMode

	// SnpGap and IndelGap correspond to -g/--SnpGap and -G/--IndelGap.
	SnpGap   int
	IndelGap int

	// Regions / RegionsFile and Targets / TargetsFile run as post-filters in v1.
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string

	// NoVersion suppresses the ##bcftools_filterCommand= header line we
	// would otherwise inject.
	NoVersion bool
}

// VCFFilter streams a VCF/BCF source through opts and writes the
// soft-filtered output to out. It is the in-memory entry point for tests
// and CLI tools that already have an io.Reader.
func VCFFilter(in io.Reader, out io.Writer, opts VCFFilterOptions) (int, error) {
	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return 0, err
	}
	return writeFiltered(hdr, variants, out, opts)
}

// VCFFilterFile opens the named input through iohelper and emits the
// soft-filtered stream to out.
func VCFFilterFile(path string, out io.Writer, opts VCFFilterOptions) (int, error) {
	if opts.RegionsFile != "" {
		regs, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools filter: %w", err)
		}
		opts.Regions = append(opts.Regions, regs...)
	}
	if opts.TargetsFile != "" {
		regs, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools filter: %w", err)
		}
		opts.Targets = append(opts.Targets, regs...)
	}
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools filter: open %s: %w", path, err)
	}
	defer r.Close()
	return VCFFilter(r, out, opts)
}

// writeFiltered applies the soft-filter pipeline. The algorithm:
//
//  1. Compile -i / -e expressions.
//  2. Pre-pass the records to find indel positions per chromosome (so we
//     can answer "is there an indel within snpGap bp?" in step 4).
//  3. Inject the ##FILTER=<ID=name,...> header line if it's not already
//     present.
//  4. For each record, evaluate the expression. If it FAILS, set the
//     FILTER column according to Mode; optionally rewrite GTs. If it
//     PASSES with mode 'x', force FILTER to PASS.
func writeFiltered(hdr *vcf.Header, variants []*vcf.Variant, out io.Writer, opts VCFFilterOptions) (int, error) {
	include, exclude, err := compileExpressions(ViewOptions{
		IncludeExpr: opts.IncludeExpr,
		ExcludeExpr: opts.ExcludeExpr,
	})
	if err != nil {
		return 0, fmt.Errorf("bcftools filter: %w", err)
	}

	postFilters := append([]string{}, opts.Targets...)
	postFilters = append(postFilters, opts.Regions...)
	parsedTargets, err := parseRegions(postFilters)
	if err != nil {
		return 0, fmt.Errorf("bcftools filter: %w", err)
	}

	// Pre-compute indel positions per chromosome for SnpGap/IndelGap.
	indelPos := indelIndex(variants)

	// Resolve the soft-filter name. The literal "+" means "auto-pick a
	// unique name"; we follow upstream's "FilterN" convention but emit
	// just one name per invocation (every failing record uses it).
	filterName := opts.SoftFilter
	if filterName == "+" {
		filterName = uniqueFilterName(hdr)
	}
	if filterName != "" {
		ensureSoftFilterHeader(hdr, filterName)
	}

	if !opts.NoVersion {
		hdr.MetaInfo = append(hdr.MetaInfo, "##bcftools_filterCommand=filter")
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return 0, fmt.Errorf("bcftools filter: %w", err)
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

		fails := evaluateFilter(v, include, exclude)
		if !fails && (opts.SnpGap > 0 || opts.IndelGap > 0) {
			if violatesGap(v, indelPos, opts.SnpGap, opts.IndelGap) {
				fails = true
			}
		}

		applyFilterDecision(v, fails, filterName, opts)
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// evaluateFilter returns true when the record FAILS the expression
// (matching upstream's "soft-filter when this is true" semantics):
//   - With -e: fails when the exclude expression evaluates true.
//   - With -i: fails when the include expression evaluates false.
//   - With neither: passes (returns false).
func evaluateFilter(v *vcf.Variant, include, exclude *Filter) bool {
	if exclude != nil && exclude.Eval(v) {
		return true
	}
	if include != nil && !include.Eval(v) {
		return true
	}
	return false
}

// applyFilterDecision mutates v.Filter (and optionally v.Samples) per
// the upstream -m/--mode and -S/--set-GTs rules.
func applyFilterDecision(v *vcf.Variant, fails bool, filterName string, opts VCFFilterOptions) {
	switch {
	case fails && filterName != "":
		if opts.Mode&FilterModeAdd != 0 {
			// Append, preserving existing FILTER entries other than PASS / "."
			out := make([]string, 0, len(v.Filter)+1)
			for _, f := range v.Filter {
				if f == "" || f == "." || f == "PASS" {
					continue
				}
				if f == filterName {
					// Already present; nothing to add later.
					return
				}
				out = append(out, f)
			}
			out = append(out, filterName)
			v.Filter = out
		} else {
			// Replace.
			v.Filter = []string{filterName}
		}
		if opts.SetGTs != SetGTsOff {
			rewriteFailingGTs(v, opts.SetGTs)
		}
	case !fails && opts.Mode&FilterModeReset != 0:
		// Reset on pass: force PASS.
		v.Filter = []string{"PASS"}
	case !fails && isEmptyFilter(v.Filter):
		// Upstream vcffilter.c:715 also forces FILTER=PASS on any
		// passing record whose FILTER is empty/`.`, regardless of
		// mode. Preserve this default-PASS-on-pass behaviour.
		v.Filter = []string{"PASS"}
	}
}

// isEmptyFilter returns true for the missing-FILTER encodings: empty
// slice, single empty string, or single ".".
func isEmptyFilter(f []string) bool {
	if len(f) == 0 {
		return true
	}
	if len(f) == 1 && (f[0] == "" || f[0] == ".") {
		return true
	}
	return false
}

// rewriteFailingGTs replaces every per-sample GT entry per the upstream
// -S/--set-GTs rule. Sample data without a GT entry is left untouched.
func rewriteFailingGTs(v *vcf.Variant, mode SetGTsMode) {
	if len(v.Samples) == 0 {
		return
	}
	for si := range v.Samples {
		gt, ok := v.Samples[si].Data["GT"]
		if !ok {
			continue
		}
		// Preserve the phase separator from the original GT.
		sep := byte('/')
		if strings.ContainsRune(gt, '|') {
			sep = '|'
		}
		// Count allele slots (= ploidy) by splitting on either separator.
		slots := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
		if len(slots) == 0 {
			continue
		}
		var rep string
		switch mode {
		case SetGTsMissing:
			rep = "."
		case SetGTsRef:
			rep = "0"
		default:
			// SetGTsOff or an unrecognised mode: don't rewrite.
			continue
		}
		newSlots := make([]string, len(slots))
		for i := range newSlots {
			newSlots[i] = rep
		}
		v.Samples[si].Data["GT"] = strings.Join(newSlots, string(sep))
	}
}

// indelIndex returns, per chromosome, the sorted set of POSitions occupied
// by indel records (REF or any ALT length differs from 1). We use a
// sorted slice so the gap check can binary-search.
func indelIndex(variants []*vcf.Variant) map[string][]int {
	out := make(map[string][]int)
	for _, v := range variants {
		if !isIndel(v) {
			continue
		}
		out[v.Chrom] = append(out[v.Chrom], v.Pos)
	}
	for k := range out {
		sort.Ints(out[k])
	}
	return out
}

// violatesGap returns true if v sits within snpGap bp of an indel (when
// v itself is a SNP) or within indelGap bp of another indel (when v is
// an indel). Position comparison is on the leftmost POS only — this
// matches upstream's first-pass behaviour for non-overlap indel modes.
func violatesGap(v *vcf.Variant, indels map[string][]int, snpGap, indelGap int) bool {
	positions := indels[v.Chrom]
	if len(positions) == 0 {
		return false
	}
	gap := 0
	if isIndel(v) {
		gap = indelGap
	} else if snpGap > 0 {
		gap = snpGap
	}
	if gap <= 0 {
		return false
	}
	// Binary-search for the nearest indel.
	lo, hi := 0, len(positions)
	for lo < hi {
		mid := (lo + hi) / 2
		if positions[mid] < v.Pos {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	candidates := []int{}
	if lo < len(positions) {
		candidates = append(candidates, positions[lo])
	}
	if lo > 0 {
		candidates = append(candidates, positions[lo-1])
	}
	for _, p := range candidates {
		if p == v.Pos {
			// Same record; skip.
			continue
		}
		d := p - v.Pos
		if d < 0 {
			d = -d
		}
		if d <= gap {
			return true
		}
	}
	return false
}

// ensureSoftFilterHeader injects a ##FILTER=<ID=name,Description=...>
// header line if one with the same ID isn't present yet.
func ensureSoftFilterHeader(hdr *vcf.Header, name string) {
	if hdr == nil || name == "" {
		return
	}
	for _, m := range hdr.MetaInfo {
		k, id := structuredID(m)
		if k == "FILTER" && id == name {
			return
		}
	}
	line := fmt.Sprintf(`##FILTER=<ID=%s,Description="Set if not true: see bcftools filter expression">`, name)
	hdr.MetaInfo = append(hdr.MetaInfo, line)
}

// uniqueFilterName mirrors upstream's "+" sentinel: when -s + is used,
// pick the smallest "FilterN" that doesn't already appear in the header.
func uniqueFilterName(hdr *vcf.Header) string {
	used := make(map[string]bool)
	if hdr != nil {
		for _, m := range hdr.MetaInfo {
			k, id := structuredID(m)
			if k == "FILTER" {
				used[id] = true
			}
		}
	}
	for i := 1; i < 1<<30; i++ {
		cand := fmt.Sprintf("Filter%d", i)
		if !used[cand] {
			return cand
		}
	}
	return "Filter"
}
