// Shared region/target (-r/-R/-t/-T) selection for the native plugin
// framework. Upstream bcftools applies these options uniformly across every
// in-tree plugin via htslib's bcf_sr_set_regions / bcf_sr_set_targets; this
// port re-creates the same selection in pure Go so a single filter can be
// applied by the native-plugin host before any record reaches a plugin's
// Process, instead of every plugin re-implementing (or rejecting) the flags.
//
// Semantics replicated byte-for-byte against upstream bcftools 1.23.1:
//
//   - -r/--regions and -R/--regions-file are REGION mode: a record is kept if
//     its span [POS, POS+len(REF)-1] OVERLAPS any region. (Upstream uses the
//     index to jump; for a fully decoded stream an overlap test is equivalent.)
//
//   - -t/--targets and -T/--targets-file are TARGET mode: a streaming
//     POSITIONAL filter keyed on the record START position — a record is kept
//     if its POS falls within [beg, end] of any target. A leading '^' on the
//     target list / file path NEGATES the selection (exclude matches).
//
//   - The -t vs -r difference is exactly upstream's: -r is span-overlap based,
//     -t is start-position based. An indel at POS=100 spanning 100..104 is
//     INCLUDED by `-r chr:102-102` (overlap) but EXCLUDED by `-t chr:102-102`
//     (its start, 100, is not in 102..102), matching upstream verbatim.
//
//   - File-format autodetection for -R/-T mirrors htslib's synced-reader
//     bcf_sr_regions_init (the path every in-tree plugin uses): a path ending in
//     .bed or .bed.gz is parsed as BED (0-based half-open, converted to 1-based
//     inclusive); any other file is TAB-separated with two columns meaning
//     chr<TAB>pos (a single 1-based position) and three or more meaning
//     chr<TAB>beg<TAB>end (1-based inclusive). Unlike `bcftools view`'s regidx,
//     the synced reader does NOT accept a bare contig or a chr:beg-end region
//     string inside a file — such a single-column line is a parse error, which
//     this loader reproduces. Inline -r/-t region STRINGS (chr, chr:pos,
//     chr:beg-end, chr:beg-) are still 1-based inclusive and colon-spelled.
package bcftools

import (
	"bufio"
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// regionSpec is one -r region token (or one -R file line) together with the
// verbatim label upstream uses when it reports per-region results. For an inline
// -r token the label is the token itself ("chr1", "chr1:200-405"); for an -R
// file line it is the synthesised "chr:beg-end" form. check-sparsity is the one
// plugin whose output is grouped and labelled by region, so it consults
// RegionSpecs directly rather than the flattened region slice.
type regionSpec struct {
	label  string
	region region
}

// Overlap modes for --regions-overlap / --targets-overlap. They mirror
// htslib's BCF_SR_REGIONS_OVERLAP / BCF_SR_TARGETS_OVERLAP integer codes and
// bcftools' parse_overlap_option spellings ("pos"/0, "record"/1, "variant"/2).
// The mode selects the record interval [beg,end] used for the region/target
// overlap test:
//
//   - overlapPos     (0): beg = end = POS                       (POS-in-region)
//   - overlapRecord  (1): beg = POS, end = POS + rlen - 1       (VCF line span)
//   - overlapVariant (2): beg = POS + off, end = POS + rlen - 1 (true variant
//     span, where off is the longest leading run of bases common to REF and
//     every ALT, capped at rlen — so a left-anchored indel's shared first base
//     is excluded from the span). See htslib _set_variant_boundaries.
//
// Upstream's synced reader defaults are regions=1 (record) and targets=0 (pos),
// set in bcf_sr_init: BCF_SR_REGIONS_OVERLAP=1, BCF_SR_TARGETS_OVERLAP=0.
const (
	overlapPos     = 0
	overlapRecord  = 1
	overlapVariant = 2
)

// regionTargetFilter captures the parsed -r/-R/-t/-T selection. The zero value
// is a no-op filter (active() reports false), so plugins that received no
// region/target option pay nothing. Note the zero value also encodes the
// upstream overlap defaults (regionsOverlap == overlapRecord == 1 only after
// parseRegionTargetArgs initialises it; the targets default overlapPos == 0
// coincides with the zero value). parseRegionTargetArgs always sets both modes,
// so any filter produced there carries the correct defaults.
type regionTargetFilter struct {
	// regions are -r/-R windows, matched by span overlap.
	regions []region
	// regionSpecs preserves each region's verbatim label, in order, for the few
	// plugins (check-sparsity) that group/label output by region.
	regionSpecs []regionSpec
	// targets are -t/-T windows, matched by record start position.
	targets []region
	// targetsNegated inverts the target match (the leading '^').
	targetsNegated bool
	// hasRegions / hasTargets record whether the corresponding option was
	// supplied at all (an option supplied but matching nothing must still drop
	// every record, distinct from "option absent").
	hasRegions bool
	hasTargets bool
	// regionsOverlap / targetsOverlap are the --regions-overlap /
	// --targets-overlap modes (overlapPos/Record/Variant). parseRegionTargetArgs
	// initialises them to the upstream defaults (record for regions, pos for
	// targets) and overrides them when the option is supplied.
	regionsOverlap int
	targetsOverlap int
}

// active reports whether any region/target selection was requested.
func (f *regionTargetFilter) active() bool {
	return f != nil && (f.hasRegions || f.hasTargets)
}

// keep reports whether the variant passes the region AND target selection.
// When both -r/-R and -t/-T are given they are ANDed, matching upstream
// (regions restrict the index scan, targets then post-filter the survivors).
// The region and target match each use the variant interval implied by their
// overlap mode (regionsOverlap / targetsOverlap), tested against the windows.
func (f *regionTargetFilter) keep(v *vcf.Variant) bool {
	if f.hasRegions {
		beg, end := overlapBoundaries(v, f.regionsOverlap)
		if !intervalOverlapsAny(v.Chrom, beg, end, f.regions) {
			return false
		}
	}
	if f.hasTargets {
		beg, end := overlapBoundaries(v, f.targetsOverlap)
		match := intervalOverlapsAny(v.Chrom, beg, end, f.targets)
		if match == f.targetsNegated {
			return false
		}
	}
	return true
}

// overlapBoundaries returns the 1-based inclusive interval [beg,end] of the
// variant under the given overlap mode, replicating htslib's per-mode interval
// selection (synced_bcf_reader.c). rlen is len(REF): for concrete (non-symbolic)
// alleles this equals htslib's rec->rlen exactly.
//
//   - overlapPos:     beg = end = POS.
//   - overlapRecord:  beg = POS, end = POS + rlen - 1.
//   - overlapVariant: beg = POS + off, end = POS + rlen - 1, where off is the
//     longest run of leading bases common to REF and EVERY ALT allele, capped at
//     rlen (htslib _set_variant_boundaries). With no ALT alleles off is 0.
func overlapBoundaries(v *vcf.Variant, mode int) (beg, end int) {
	rlen := len(v.Ref)
	if rlen < 1 {
		rlen = 1
	}
	switch mode {
	case overlapPos:
		return v.Pos, v.Pos
	case overlapVariant:
		off := variantLeadingOffset(v.Ref, v.Alt, rlen)
		return v.Pos + off, v.Pos + rlen - 1
	default: // overlapRecord
		return v.Pos, v.Pos + rlen - 1
	}
}

// variantLeadingOffset returns the length of the longest prefix shared by ref
// and every alt allele, capped at rlen. It mirrors the loop in htslib's
// _set_variant_boundaries: off starts at rlen and is reduced to the shortest
// common-prefix length across all ALT alleles (stopping early at 0). A symbolic
// or empty ALT (one starting with '<' or '*', or REF mismatch at position 0)
// yields off 0, leaving the full record span.
func variantLeadingOffset(ref string, alts []string, rlen int) int {
	off := rlen
	for _, alt := range alts {
		j := 0
		for j < len(ref) && j < len(alt) && ref[j] == alt[j] {
			j++
		}
		if off > j {
			off = j
		}
		if off == 0 {
			break
		}
	}
	if off > rlen {
		off = rlen
	}
	return off
}

// intervalOverlapsAny reports whether the 1-based inclusive interval [beg,end]
// on chrom overlaps any of the windows, the same predicate htslib uses
// (region.start <= end && region.end >= beg).
func intervalOverlapsAny(chrom string, beg, end int, regions []region) bool {
	for _, r := range regions {
		if r.chrom != chrom {
			continue
		}
		if r.beg <= end && r.end >= beg {
			return true
		}
	}
	return false
}

// apply returns the subset of variants that pass the filter, preserving order.
// A nil/inactive filter returns the input unchanged.
func (f *regionTargetFilter) apply(variants []*vcf.Variant) []*vcf.Variant {
	if !f.active() {
		return variants
	}
	out := variants[:0:0]
	for _, v := range variants {
		if f.keep(v) {
			out = append(out, v)
		}
	}
	return out
}

// startInAny reports whether the variant's START position (1-based POS) falls
// within any of the listed 1-based inclusive windows. This is the -t/-T
// positional semantics, distinct from overlapsAny (span overlap) used by -r.
func startInAny(v *vcf.Variant, regions []region) bool {
	for _, r := range regions {
		if r.chrom != v.Chrom {
			continue
		}
		if v.Pos >= r.beg && v.Pos <= r.end {
			return true
		}
	}
	return false
}

// regionTargetCaps reports which of the region (-r/-R) and target (-t/-T)
// option families a plugin actually exposes to the shared filter. The default is
// the empty {false,false}: the shared filter consumes NOTHING and leaves every
// flag for the plugin's own parser. This is deliberately conservative because
// several plugins repurpose these letters (tag2tag's -r is --replace, -t is
// --tags; guess-ploidy's -t is --tag), so a plugin opts IN to the shared filter
// by declaring the families it genuinely uses for region/target selection.
type regionTargetCaps struct {
	regions bool // plugin exposes -r/-R (region overlap)
	targets bool // plugin exposes -t/-T (target start-position)
	// overlapOption reports whether the plugin exposes --regions-overlap /
	// --targets-overlap. Only the plugins that upstream genuinely wires to
	// bcf_sr_set_opt(BCF_SR_*_OVERLAP) accept these flags (contrast, split,
	// scatter), plus mendelian2 (which advertises them — upstream has a getopt
	// bug that drops them, fixed-on-port) and trio-dnm3 (which upstream
	// consumes-and-ignores — honoured here). For every other plugin the option is
	// passed through to the plugin's own parser, which rejects it just as upstream
	// does for an unknown option. The overlap MODES still default to the upstream
	// synced-reader values (record for regions, pos for targets) regardless, so a
	// plugin that opts out still applies the correct default predicate.
	overlapOption bool
}

// allRegionTargetCaps is the common opt-in: both -r/-R and -t/-T are
// region/target selection (the htslib synced-reader convention). It does NOT
// claim the --regions-overlap/--targets-overlap flags, since most plugins do not
// accept them upstream; a plugin that does opts in with overlapRegionTargetCaps.
var allRegionTargetCaps = regionTargetCaps{regions: true, targets: true}

// overlapRegionTargetCaps opts a plugin into -r/-R/-t/-T AND the
// --regions-overlap/--targets-overlap mode flags. Used by the plugins that
// genuinely honour the overlap option (contrast, split, scatter, mendelian2,
// trio-dnm3).
var overlapRegionTargetCaps = regionTargetCaps{regions: true, targets: true, overlapOption: true}

// regionsOnlyCaps opts a plugin into -r/-R region selection while leaving -t/-T
// to the plugin's own parser (guess-ploidy, whose -t is --tag).
var regionsOnlyCaps = regionTargetCaps{regions: true, targets: false}

// regionTargetCapabler is implemented by native plugins that opt into the shared
// -r/-R/-t/-T filter, declaring which families they expose. A plugin that does
// not implement it gets the empty caps: the shared filter consumes nothing and
// the plugin parses (or rejects) the flags itself, exactly as before.
type regionTargetCapabler interface {
	// RegionTargetCaps reports which region/target option families this plugin
	// exposes, so the shared filter consumes only those.
	RegionTargetCaps() regionTargetCaps
}

// parseRegionTargetArgs extracts the -r/-R/-t/-T options (and their values)
// from args, returning the remaining args (with those options removed) and the
// assembled filter. Unknown options and positionals are passed through
// untouched so each plugin's own parser still sees them. Both the joined-value
// (-rREGION) and separate-value (-r REGION) spellings are accepted, matching
// getopt. caps restricts which option families are consumed; a flag for a
// family the plugin does not expose (e.g. -t on guess-ploidy, where it means
// --tag) is passed through untouched.
func parseRegionTargetArgs(args []string, caps regionTargetCaps) (remaining []string, filter regionTargetFilter, err error) {
	// Seed the upstream synced-reader overlap defaults (regions=record,
	// targets=pos). They are honoured whether or not the explicit option is
	// supplied.
	filter.regionsOverlap = overlapRecord
	filter.targetsOverlap = overlapPos
	i := 0
	value := func(flag, joined string) (string, error) {
		if joined != "" {
			return joined, nil
		}
		if i+1 >= len(args) {
			return "", fmt.Errorf("option %q requires an argument", flag)
		}
		i++
		return args[i], nil
	}
	for i = 0; i < len(args); i++ {
		a := args[i]
		var flag, joined string
		switch {
		case a == "-r" || a == "--regions":
			flag = "-r"
		case a == "-R" || a == "--regions-file":
			flag = "-R"
		case a == "-t" || a == "--targets":
			flag = "-t"
		case a == "-T" || a == "--targets-file":
			flag = "-T"
		case strings.HasPrefix(a, "-r") && len(a) > 2 && a[1] != '-':
			flag, joined = "-r", a[2:]
		case strings.HasPrefix(a, "-R") && len(a) > 2 && a[1] != '-':
			flag, joined = "-R", a[2:]
		case strings.HasPrefix(a, "-t") && len(a) > 2 && a[1] != '-':
			flag, joined = "-t", a[2:]
		case strings.HasPrefix(a, "-T") && len(a) > 2 && a[1] != '-':
			flag, joined = "-T", a[2:]
		case strings.HasPrefix(a, "--regions="):
			flag, joined = "-r", a[len("--regions="):]
		case strings.HasPrefix(a, "--regions-file="):
			flag, joined = "-R", a[len("--regions-file="):]
		case strings.HasPrefix(a, "--targets="):
			flag, joined = "-t", a[len("--targets="):]
		case strings.HasPrefix(a, "--targets-file="):
			flag, joined = "-T", a[len("--targets-file="):]
		case a == "--regions-overlap":
			flag = "--regions-overlap"
		case a == "--targets-overlap":
			flag = "--targets-overlap"
		case strings.HasPrefix(a, "--regions-overlap="):
			flag, joined = "--regions-overlap", a[len("--regions-overlap="):]
		case strings.HasPrefix(a, "--targets-overlap="):
			flag, joined = "--targets-overlap", a[len("--targets-overlap="):]
		default:
			remaining = append(remaining, a)
			continue
		}

		// Respect the plugin's declared capabilities: a flag for a family the
		// plugin does not expose (e.g. -t on guess-ploidy = --tag) is left for
		// the plugin's own parser. The overlap options ride along with the family
		// they configure: --regions-overlap with the region family, and
		// --targets-overlap with the target family.
		if (flag == "-r" || flag == "-R") && !caps.regions {
			remaining = append(remaining, a)
			continue
		}
		if (flag == "-t" || flag == "-T") && !caps.targets {
			remaining = append(remaining, a)
			continue
		}
		// The overlap-mode flags are only consumed by plugins that genuinely
		// expose them (caps.overlapOption); otherwise they pass through to the
		// plugin's parser, which rejects them exactly as upstream rejects an
		// unknown option.
		if flag == "--regions-overlap" && (!caps.regions || !caps.overlapOption) {
			remaining = append(remaining, a)
			continue
		}
		if flag == "--targets-overlap" && (!caps.targets || !caps.overlapOption) {
			remaining = append(remaining, a)
			continue
		}

		val, verr := value(flag, joined)
		if verr != nil {
			return nil, regionTargetFilter{}, verr
		}
		switch flag {
		case "-r":
			tokens := SplitCommaList(val)
			specs, rerr := parseRegionSpecs(tokens)
			if rerr != nil {
				return nil, regionTargetFilter{}, rerr
			}
			filter.addRegionSpecs(specs)
		case "-R":
			tokens, rerr := loadRegionTargetFile(val)
			if rerr != nil {
				return nil, regionTargetFilter{}, rerr
			}
			specs, rerr := parseRegionSpecs(tokens)
			if rerr != nil {
				return nil, regionTargetFilter{}, rerr
			}
			filter.addRegionSpecs(specs)
		case "-t":
			neg, body := splitNegation(val)
			regs, rerr := parseRegions(SplitCommaList(body))
			if rerr != nil {
				return nil, regionTargetFilter{}, rerr
			}
			filter.targets = append(filter.targets, regs...)
			filter.hasTargets = true
			if neg {
				filter.targetsNegated = true
			}
		case "-T":
			neg, path := splitNegation(val)
			specs, rerr := loadRegionTargetFile(path)
			if rerr != nil {
				return nil, regionTargetFilter{}, rerr
			}
			regs, rerr := parseRegions(specs)
			if rerr != nil {
				return nil, regionTargetFilter{}, rerr
			}
			filter.targets = append(filter.targets, regs...)
			filter.hasTargets = true
			if neg {
				filter.targetsNegated = true
			}
		case "--regions-overlap":
			mode, perr := parseOverlapOption(val)
			if perr != nil {
				return nil, regionTargetFilter{}, perr
			}
			filter.regionsOverlap = mode
		case "--targets-overlap":
			mode, perr := parseOverlapOption(val)
			if perr != nil {
				return nil, regionTargetFilter{}, perr
			}
			filter.targetsOverlap = mode
		}
	}
	return remaining, filter, nil
}

// parseOverlapOption maps a --regions-overlap/--targets-overlap MODE string to
// its integer code, replicating bcftools' parse_overlap_option (version.c):
// "pos"/0, "record"/1, "variant"/2, case-insensitive on the words. An
// unrecognised value is an error.
func parseOverlapOption(arg string) (int, error) {
	switch {
	case strings.EqualFold(arg, "pos") || arg == "0":
		return overlapPos, nil
	case strings.EqualFold(arg, "record") || arg == "1":
		return overlapRecord, nil
	case strings.EqualFold(arg, "variant") || arg == "2":
		return overlapVariant, nil
	default:
		return -1, fmt.Errorf("could not parse the overlap option: %q", arg)
	}
}

// parseRegionSpecs turns region tokens ("chr1", "chr1:200-405", ...) into
// labelled specs, keeping the verbatim token as the label so region-grouped
// plugins (check-sparsity) can reproduce upstream's per-region report keys.
func parseRegionSpecs(tokens []string) ([]regionSpec, error) {
	specs := make([]regionSpec, 0, len(tokens))
	for _, tok := range tokens {
		regs, err := parseRegions([]string{tok})
		if err != nil {
			return nil, err
		}
		for _, r := range regs {
			specs = append(specs, regionSpec{label: tok, region: r})
		}
	}
	return specs, nil
}

// addRegionSpecs records the labelled region specs and their flattened regions,
// marking the region family active.
func (f *regionTargetFilter) addRegionSpecs(specs []regionSpec) {
	for _, s := range specs {
		f.regionSpecs = append(f.regionSpecs, s)
		f.regions = append(f.regions, s.region)
	}
	f.hasRegions = true
}

// splitNegation strips a leading '^' (upstream's target-negation marker) from a
// targets list or targets-file path, returning whether it was present.
func splitNegation(s string) (negated bool, rest string) {
	if strings.HasPrefix(s, "^") {
		return true, s[1:]
	}
	return false, s
}

// loadRegionTargetFile reads a -R/-T file into a slice of 1-based-inclusive
// region-list specs ("chr:beg-end"), autodetecting the format the way htslib's
// regidx does: a .bed/.bed.gz path is BED (0-based half-open); any other
// tab/whitespace-separated file is 1-based inclusive, with two columns meaning
// a single position and three+ meaning beg..end. Comment (#) and blank lines
// are skipped. The file may be gzip-compressed (handled transparently).
func loadRegionTargetFile(path string) ([]string, error) {
	isBED := strings.HasSuffix(path, ".bed") || strings.HasSuffix(path, ".bed.gz")
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
		fields := strings.Fields(line)
		switch {
		case len(fields) == 1:
			// The synced reader (the plugin region/target path) requires at
			// least the chr<TAB>pos columns; a single-column line (a bare contig
			// or a chr:beg-end region string) is a parse error, exactly as
			// upstream's bcf_sr_regions_init reports it.
			return nil, fmt.Errorf("could not parse line %q of %q, expected columns chr,pos[,end]", line, path)
		case len(fields) == 2:
			// chr<TAB>pos: a single 1-based position. htslib treats a
			// two-column tab file as chr,pos regardless of a .bed suffix.
			pos, perr := strconv.Atoi(fields[1])
			if perr != nil {
				return nil, fmt.Errorf("bad position in %q: %q", path, line)
			}
			out = append(out, fmt.Sprintf("%s:%d-%d", fields[0], pos, pos))
		default:
			beg, perr := strconv.Atoi(fields[1])
			if perr != nil {
				return nil, fmt.Errorf("bad start in %q: %q", path, line)
			}
			end, perr := strconv.Atoi(fields[2])
			if perr != nil {
				return nil, fmt.Errorf("bad end in %q: %q", path, line)
			}
			if isBED {
				// BED is 0-based half-open [beg,end); convert to 1-based
				// inclusive [beg+1, end].
				beg++
			}
			out = append(out, fmt.Sprintf("%s:%d-%d", fields[0], beg, end))
		}
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// LoadRegionTargetFile is the exported entry point used by the host CLI to
// resolve a -R/-T file path to 1-based region-list specs with the same
// autodetection the native plugin filter uses. It is a thin wrapper over
// loadRegionTargetFile.
func LoadRegionTargetFile(path string) ([]string, error) {
	return loadRegionTargetFile(path)
}
