package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// OutputFormat describes how `bcftools view` should emit records.
type OutputFormat int

const (
	// OutputVCF is uncompressed VCF text (the default).
	OutputVCF OutputFormat = iota
	// OutputVCFGz is gzip-compressed VCF.
	OutputVCFGz
	// OutputBCF marks compressed BCF output. The csq path writes it via
	// openCSQOutput; the view runner still defers it (see viewBCFStream).
	OutputBCF
	// OutputBCFUncompressed marks uncompressed BCF output. As with
	// OutputBCF, csq writes it but view defers it.
	OutputBCFUncompressed
)

// ParseOutputFormat turns the short bcftools -O codes into the typed value.
func ParseOutputFormat(s string) (OutputFormat, error) {
	switch strings.ToLower(s) {
	case "", "v":
		return OutputVCF, nil
	case "z":
		return OutputVCFGz, nil
	case "b":
		return OutputBCF, nil
	case "u":
		return OutputBCFUncompressed, nil
	}
	return 0, fmt.Errorf("bcftools view: unknown output format %q (expect v, z, u, b)", s)
}

// ViewOptions controls the behaviour of View / ViewFile.
type ViewOptions struct {
	OutputFormat   OutputFormat
	HeaderOnly     bool
	NoHeader       bool
	DropGenotypes  bool
	MinAlleleCount int     // -c
	MaxAlleleCount int     // -C (-1 disables)
	MinAlleleFreq  float64 // -q
	MaxAlleleFreq  float64 // -Q (-1 disables)
	IncludeExpr    string  // -i
	ExcludeExpr    string  // -e
	ApplyFilters   []string
	Regions        []string // -r
	RegionsFile    string   // -R
	Targets        []string // -t (post-filter)
	TargetsFile    string   // -T
	Samples        []string // -s
	SamplesFile    string   // -S
	// Private (-x/--private) keeps only sites whose non-reference alleles are
	// carried exclusively by the subset samples: the non-reference allele
	// count within the subset is greater than zero AND equals the
	// non-reference allele count across all samples (none outside the subset).
	// It is meaningful only together with a sample subset (-s/-S).
	Private bool
	// ExcludePrivate (-X/--exclude-private) is the inverse of Private: it drops
	// sites whose non-reference alleles are carried exclusively by the subset
	// samples. Like Private, it requires a sample subset to take effect.
	ExcludePrivate bool
	CompressLevel  int // -l
	// Threads is the value of -@/--threads. When greater than 1 and the
	// selected output format is BGZF-framed (-O z VCF.gz or -O b BCF), the
	// output is compressed by a pool of that many worker goroutines via
	// bgzf.MultiWriter. The framed result decodes byte-identically regardless
	// of the thread count. A value of 0 or 1 uses the single-threaded writer.
	Threads int
	// IncludeTypes (-v/--types) keeps only records that include at least one
	// of the named variant types; ExcludeTypes (-V/--exclude-types) drops
	// records that include any of the named types. Type names are the
	// upstream lowercase set: snps, indels, mnps, ref, bnd, other.
	IncludeTypes []string
	ExcludeTypes []string
	// NoUpdate (-I/--no-update) suppresses the recomputation of INFO/AC and
	// INFO/AN after a sample subset (-s/-S). By default, like upstream
	// bcftools view, AC and AN are recomputed from the kept genotypes (and
	// added to the header/record when absent).
	NoUpdate bool
	// SuppressPASSFilter prevents openOutput from re-injecting the implicit
	// ##FILTER=<ID=PASS> header line. It is set by `annotate -x FILTER`, which
	// (like upstream remove_hdr_lines(BCF_HL_FLT)) removes every FILTER header
	// line — PASS included — and they must not reappear at write time.
	SuppressPASSFilter bool
	// CalcAC mirrors upstream vcfview.c's calc_ac flag. It is set by the CLI
	// when any of -c/-C/-q/-Q (allele count/frequency selectors) or -x/-X
	// (private) is requested. When CalcAC is true and NoUpdate is false,
	// INFO/AC and INFO/AN are (re)computed and appended to every output
	// record exactly as upstream does — even when no records are dropped and
	// no sample subset is applied. A sample subset (-s/-S) without -I also
	// triggers recomputation; that is handled directly by the Samples path.
	CalcAC bool
}

// recalcAC reports whether INFO/AC and INFO/AN should be (re)computed and
// written to output records. It mirrors upstream vcfview.c, where
// calc_ac (set by -c/-C/-q/-Q/-x/-X) OR a sample subset (-s/-S) enables the
// update, and -I/--no-update suppresses it.
func (o ViewOptions) recalcAC() bool {
	if o.NoUpdate {
		return false
	}
	return o.CalcAC || len(o.Samples) > 0
}

// variantTypeMask returns the OR of the per-allele variant-type bits for v.
// A record with no ALT (or only the missing/ref allele) has mask 0, which is
// the "ref" type.
func variantTypeMask(v *vcf.Variant) int {
	mask := 0
	for _, alt := range v.Alt {
		mask |= variantTypeBit(v.Ref, alt)
	}
	return mask
}

// typeNameBits maps an upstream --types name to its variant-type bit. The
// "ref" type is the special mask-0 case and is reported via the bool.
func typeNameBits(name string) (bit int, isRef bool) {
	switch strings.ToLower(name) {
	case "snps", "snp":
		return vtSNP, false
	case "indels", "indel":
		return vtINDEL, false
	case "mnps", "mnp":
		return vtMNP, false
	case "bnd":
		return vtBND, false
	case "other":
		return vtOTHER, false
	case "overlap":
		return vtOVERLAP, false
	case "ref":
		return 0, true
	}
	return 0, false
}

// matchesTypeSet reports whether v's variant types intersect the named set.
func matchesTypeSet(v *vcf.Variant, names []string) bool {
	mask := variantTypeMask(v)
	for _, n := range names {
		bit, isRef := typeNameBits(n)
		if isRef {
			// "ref" matches a record whose only allele classifies as the
			// reference (alt equals ref) — not a missing/absent ALT ("."),
			// which upstream does not treat as a ref-type record.
			if mask == 0 && len(v.Alt) > 0 && v.Alt[0] != "." && v.Alt[0] != "" {
				return true
			}
			continue
		}
		if bit != 0 && mask&bit != 0 {
			return true
		}
	}
	return false
}

// passesTypeFilter applies the -v/--types (include) and -V/--exclude-types
// (exclude) selectors, mirroring upstream bcftools view: a record is kept by
// -v when it includes at least one requested type, and dropped by -V when it
// includes any excluded type.
func (o ViewOptions) passesTypeFilter(v *vcf.Variant) bool {
	if len(o.IncludeTypes) > 0 && !matchesTypeSet(v, o.IncludeTypes) {
		return false
	}
	if len(o.ExcludeTypes) > 0 && matchesTypeSet(v, o.ExcludeTypes) {
		return false
	}
	return true
}

// recomputeACAN recomputes INFO/AC (per-ALT, Number=A) and INFO/AN (total
// called alleles) from v's current sample genotypes and writes them back into
// v.Info, mirroring upstream `bcftools view -s` (which updates AC/AN after a
// sample subset). The tags are appended to the INFO order when absent. AF and
// other INFO tags are left untouched, matching upstream.
//
// It always computes from the per-sample GT data, which is the upstream
// behaviour for the sample-subset path (bcf_calc_ac with BCF_UN_FMT only).
func recomputeACAN(v *vcf.Variant) {
	ac, an, ok := acanFromGT(v)
	if !ok {
		return
	}
	writeACAN(v, ac, an)
}

// updateACAN mirrors upstream vcfview.c's non-subset calc_ac path. Upstream
// calls bcf_calc_ac(hdr,line,ac,BCF_UN_INFO|BCF_UN_FMT), which prefers the
// pre-existing INFO/AC and INFO/AN when both are present and the AC arity
// matches the number of ALT alleles, and otherwise computes them from the GT
// FORMAT field. The resulting values are then written back to INFO/AC and
// INFO/AN whenever update_info is set. updateACAN reproduces that exactly:
// existing INFO/AC and INFO/AN are preserved when consistent, otherwise they
// are recomputed from GT.
func updateACAN(v *vcf.Variant) {
	if ac, an, ok := acanFromINFO(v); ok {
		writeACAN(v, ac, an)
		return
	}
	if ac, an, ok := acanFromGT(v); ok {
		writeACAN(v, ac, an)
	}
}

// acanFromINFO reads per-ALT allele counts and the total allele number from
// the pre-existing INFO/AC and INFO/AN tags. It returns ok=false (mirroring
// upstream bcf_calc_ac) when either tag is absent, AN cannot be parsed, or the
// number of AC entries does not equal the number of ALT alleles.
func acanFromINFO(v *vcf.Variant) (ac []int, an int, ok bool) {
	acStr := v.Info["AC"]
	anStr := v.Info["AN"]
	if acStr == "" || anStr == "" {
		return nil, 0, false
	}
	an, err := strconv.Atoi(anStr)
	if err != nil || an < 0 {
		return nil, 0, false
	}
	parts := strings.Split(acStr, ",")
	if len(parts) != len(v.Alt) {
		return nil, 0, false
	}
	ac = make([]int, len(parts))
	for i, p := range parts {
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, 0, false
		}
		ac[i] = n
	}
	return ac, an, true
}

// acanFromGT computes per-ALT allele counts and the total called-allele number
// (AN) from v's per-sample GT data. It returns ok=false when the record has no
// ALT alleles, leaving the caller to skip the update as upstream does.
func acanFromGT(v *vcf.Variant) (ac []int, an int, ok bool) {
	if len(v.Alt) == 0 {
		return nil, 0, false
	}
	ac = make([]int, len(v.Alt))
	for _, s := range v.Samples {
		gt, has := s.Data["GT"]
		if !has {
			continue
		}
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			if a == "." || a == "" {
				continue
			}
			n, err := strconv.Atoi(a)
			if err != nil {
				continue
			}
			an++
			if n >= 1 && n <= len(ac) {
				ac[n-1]++
			}
		}
	}
	return ac, an, true
}

// applyACAN performs the per-record AC/AN handling for a kept variant,
// restricting samples first when a subset is requested and then (re)computing
// INFO/AC and INFO/AN as upstream vcfview.c does: from the kept genotypes after
// a subset, or from INFO-preferred-then-GT otherwise. It is a no-op when
// recomputation is disabled (NoUpdate, or neither calc_ac nor a subset).
func applyACAN(v *vcf.Variant, opts ViewOptions) {
	if len(opts.Samples) > 0 {
		restrictSamples(v, opts.Samples)
		if !opts.NoUpdate {
			recomputeACAN(v)
		}
		return
	}
	if opts.recalcAC() {
		updateACAN(v)
	}
}

// writeACAN stores per-ALT allele counts in INFO/AC and the total allele
// number in INFO/AN, appending the tags to the INFO order when absent.
func writeACAN(v *vcf.Variant, ac []int, an int) {
	acParts := make([]string, len(ac))
	for i, c := range ac {
		acParts[i] = strconv.Itoa(c)
	}
	setInfo(v, "AC", strings.Join(acParts, ","))
	setInfo(v, "AN", strconv.Itoa(an))
}

// ensureACANHeader appends the ##INFO header lines for AC and AN when they are
// absent, using the exact definitions upstream bcftools emits (via
// bcf_hdr_append, which adds new lines at the end of the meta block), so a
// sample-subset that introduces AC/AN produces a valid, parity-matching header.
func ensureACANHeader(hdr *vcf.Header) {
	if hdr == nil {
		return
	}
	have := map[string]bool{}
	for _, m := range hdr.MetaInfo {
		if k, id := structuredID(m); k == "INFO" {
			have[id] = true
		}
	}
	if !have["AC"] {
		hdr.MetaInfo = append(hdr.MetaInfo, `##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">`)
	}
	if !have["AN"] {
		hdr.MetaInfo = append(hdr.MetaInfo, `##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">`)
	}
}

// applyAlleleFilters returns true if the variant passes the AC/AF filters.
func (o ViewOptions) applyAlleleFilters(v *vcf.Variant) bool {
	if o.MinAlleleCount == 0 && o.MaxAlleleCount <= 0 && o.MinAlleleFreq == 0 && o.MaxAlleleFreq <= 0 {
		return true
	}
	ac, an := computeAC(v)
	if o.MinAlleleCount > 0 && ac < o.MinAlleleCount {
		return false
	}
	if o.MaxAlleleCount > 0 && ac > o.MaxAlleleCount {
		return false
	}
	if an > 0 {
		af := float64(ac) / float64(an)
		if o.MinAlleleFreq > 0 && af < o.MinAlleleFreq {
			return false
		}
		if o.MaxAlleleFreq > 0 && af > o.MaxAlleleFreq {
			return false
		}
	}
	return true
}

// computeAC returns the total non-reference allele count and total called
// allele count across all samples. If no GT information is available it
// falls back to the INFO AC/AN tags.
func computeAC(v *vcf.Variant) (ac, an int) {
	for _, s := range v.Samples {
		gt, ok := s.Data["GT"]
		if !ok {
			continue
		}
		// Split on either separator.
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			if a == "." || a == "" {
				continue
			}
			an++
			if n, err := strconv.Atoi(a); err == nil && n > 0 {
				ac++
			}
		}
	}
	if an == 0 {
		// Fall back to INFO/AC and INFO/AN if no per-sample GTs.
		if v.Info["AC"] != "" {
			parts := strings.Split(v.Info["AC"], ",")
			for _, p := range parts {
				if n, err := strconv.Atoi(p); err == nil {
					ac += n
				}
			}
		}
		if v.Info["AN"] != "" {
			if n, err := strconv.Atoi(v.Info["AN"]); err == nil {
				an = n
			}
		}
	}
	return ac, an
}

// nonRefACOver returns the total non-reference allele count summed across the
// samples in v whose names appear in want. When want is nil every sample is
// counted. The second return value reports whether any of the counted samples
// carried a GT field at all (mirroring upstream's requirement that a GT FORMAT
// field be present before the private filter is evaluated).
func nonRefACOver(v *vcf.Variant, want map[string]bool) (ac int, hasGT bool) {
	for _, s := range v.Samples {
		if want != nil && !want[s.Name] {
			continue
		}
		gt, ok := s.Data["GT"]
		if !ok {
			continue
		}
		hasGT = true
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			if a == "." || a == "" {
				continue
			}
			if n, err := strconv.Atoi(a); err == nil && n > 0 {
				ac++
			}
		}
	}
	return ac, hasGT
}

// passesPrivateFilter implements the -x/--private and -X/--exclude-private
// selectors. A site is "private" to the subset when the non-reference allele
// count within the subset is greater than zero AND equals the non-reference
// allele count across all samples (i.e. no non-reference allele is carried
// outside the subset). This mirrors upstream bcftools vcfview.c, which applies
// the test after sample subsetting and only when a GT FORMAT field is present.
//
// With no sample subset, or no GT data, the filter is a no-op (the variant is
// kept) — matching upstream, where the private check is gated on the sample
// subset being non-empty and the allele counts being recomputed from GTs.
func passesPrivateFilter(v *vcf.Variant, opts ViewOptions) bool {
	if !opts.Private && !opts.ExcludePrivate {
		return true
	}
	if len(opts.Samples) == 0 {
		return true
	}
	want := make(map[string]bool, len(opts.Samples))
	for _, name := range opts.Samples {
		want[name] = true
	}
	acFull, fullHasGT := nonRefACOver(v, nil)
	acSub, subHasGT := nonRefACOver(v, want)
	if !fullHasGT || !subHasGT {
		return true
	}
	isPrivate := acSub > 0 && acFull == acSub
	if opts.Private {
		return isPrivate
	}
	// ExcludePrivate: drop private sites, keep everything else.
	return !isPrivate
}

// applyFilterColumnFilters returns true if v passes the -f filter list.
func (o ViewOptions) applyFilterColumnFilters(v *vcf.Variant) bool {
	if len(o.ApplyFilters) == 0 {
		return true
	}
	for _, want := range o.ApplyFilters {
		for _, f := range v.Filter {
			if f == want {
				return true
			}
		}
	}
	return false
}

// dropGenotypes strips all per-sample data from v.
func dropGenotypes(v *vcf.Variant) {
	v.Format = nil
	v.Samples = nil
}

// restrictSamples keeps only the named samples in v, in the order they were
// requested. Samples present in v but not in the request list are dropped;
// requested names not present in v are silently skipped.
func restrictSamples(v *vcf.Variant, wanted []string) {
	if len(wanted) == 0 || len(v.Samples) == 0 {
		return
	}
	bySample := make(map[string]vcf.Sample, len(v.Samples))
	for _, s := range v.Samples {
		bySample[s.Name] = s
	}
	kept := make([]vcf.Sample, 0, len(wanted))
	for _, name := range wanted {
		if s, ok := bySample[name]; ok {
			kept = append(kept, s)
		}
	}
	v.Samples = kept
}

// View streams VCF or BCF from in, applies opts, and writes to out. It is
// the streaming entry point used by the `view` subcommand when no region
// query is requested (regions need a seekable file, see ViewFile).
func View(in io.Reader, out io.Writer, opts ViewOptions) (int, error) {
	return viewStreaming(in, out, opts, false, nil)
}

// ViewFile is the file-aware entry point for `bcftools view`. It supports
// region queries via a sibling .tbi index on bgzipped VCF inputs (and a
// sibling .csi index on BCF inputs) and falls back to a streaming scan
// otherwise.
func ViewFile(path string, out io.Writer, opts ViewOptions, stderr io.Writer) (int, error) {
	if len(opts.Regions) > 0 && HasCSI(path) {
		bcfHead, _ := looksLikeBCF(path)
		if bcfHead {
			return viewBCFRegions(path, out, opts, stderr)
		}
	}
	if len(opts.Regions) > 0 && hasTabixIndex(path) {
		return viewRegions(path, out, opts, stderr)
	}
	if len(opts.Regions) > 0 && stderr != nil {
		fmt.Fprintln(stderr, "bcftools view: no .tbi/.csi index found; treating -r as a post-filter (slower)")
	}
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	defer in.Close()
	// Promote -r to -t when no index is available.
	postFilters := opts.Targets
	if len(opts.Regions) > 0 {
		postFilters = append([]string{}, opts.Targets...)
		postFilters = append(postFilters, opts.Regions...)
	}
	parsedTargets, err := parseRegions(postFilters)
	if err != nil {
		return 0, err
	}
	return viewStreaming(in, out, opts, true, parsedTargets)
}

// hasTabixIndex returns true if path has a sibling .tbi file. For remote URLs
// the probe is an existence check through hfile (a HEAD/ranged request that
// 404s when the index is absent), matching htslib's remote index discovery.
func hasTabixIndex(path string) bool {
	return siblingExists(path + ".tbi")
}

// region captures one parsed `chr:start-end` window in 1-based inclusive
// coordinates (matching VCF text).
type region struct {
	chrom string
	beg   int
	end   int
}

// parseRegions turns "chr:start-end" specs into a slice of regions. If a
// spec lacks a start/end it defaults to the whole contig.
func parseRegions(specs []string) ([]region, error) {
	var out []region
	for _, s := range specs {
		colon := strings.IndexByte(s, ':')
		if colon < 0 {
			out = append(out, region{chrom: s, beg: 1, end: 1 << 30})
			continue
		}
		chrom := s[:colon]
		rest := s[colon+1:]
		dash := strings.IndexByte(rest, '-')
		var beg, end int
		var err error
		if dash < 0 {
			beg, err = strconv.Atoi(rest)
			if err != nil {
				return nil, fmt.Errorf("bcftools view: bad region %q", s)
			}
			end = beg
		} else {
			beg, err = strconv.Atoi(rest[:dash])
			if err != nil {
				return nil, fmt.Errorf("bcftools view: bad region %q", s)
			}
			if rest[dash+1:] == "" {
				end = 1 << 30
			} else {
				end, err = strconv.Atoi(rest[dash+1:])
				if err != nil {
					return nil, fmt.Errorf("bcftools view: bad region %q", s)
				}
			}
		}
		out = append(out, region{chrom: chrom, beg: beg, end: end})
	}
	return out, nil
}

// overlapsAny returns true if v's coordinates fall in any of the listed
// regions (1-based inclusive). End is taken as POS + len(REF) - 1.
func overlapsAny(v *vcf.Variant, regions []region) bool {
	if len(regions) == 0 {
		return true
	}
	vEnd := v.Pos + len(v.Ref) - 1
	if vEnd < v.Pos {
		vEnd = v.Pos
	}
	for _, r := range regions {
		if r.chrom != v.Chrom {
			continue
		}
		if vEnd >= r.beg && v.Pos <= r.end {
			return true
		}
	}
	return false
}

// viewStreaming handles VCF (gzipped or plain) and BCF inputs as a sequential
// stream. If applyTargets is true the variant must fall within targets to be
// kept; this is how the -t / -T flags plus the -r fallback path work.
func viewStreaming(in io.Reader, out io.Writer, opts ViewOptions, applyTargets bool, targets []region) (int, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return 0, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return viewBCFStream(br, out, opts, applyTargets, targets)
	}
	return viewVCFStream(br, out, opts, applyTargets, targets)
}

// viewVCFStream is the VCF (text or gzip-decoded) path. It uses the existing
// vcf.Reader + vcf.Writer types.
func viewVCFStream(in io.Reader, out io.Writer, opts ViewOptions, applyTargets bool, targets []region) (int, error) {
	r := vcf.NewReader(in)
	hdr, err := r.ReadHeader()
	if err != nil {
		return 0, err
	}

	includeF, excludeF, err := compileExpressions(opts, hdr)
	if err != nil {
		return 0, err
	}
	hdr = filterHeaderSamples(hdr, opts.Samples)
	if opts.recalcAC() {
		ensureACANHeader(hdr)
	}
	if opts.DropGenotypes {
		hdr = stripFormatLines(hdr)
		hdr.Samples = nil
	}

	w, finish, err := openOutput(out, opts, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if !opts.NoHeader {
		if err := w.WriteHeader(); err != nil {
			return 0, err
		}
	}
	if opts.HeaderOnly {
		return 0, w.Flush()
	}

	count := 0
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, err
		}
		if !keepVariant(v, opts, includeF, excludeF, applyTargets, targets) {
			continue
		}
		applyACAN(v, opts)
		if opts.DropGenotypes {
			dropGenotypes(v)
		}
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// viewBCFStream consumes a BCF stream and emits VCF text. Writing BCF is
// deferred — when callers ask for OutputBCF/OutputBCFUncompressed the
// outer runner has already rejected the request.
func viewBCFStream(in io.Reader, out io.Writer, opts ViewOptions, applyTargets bool, targets []region) (int, error) {
	br, err := bcf.NewReader(in)
	if err != nil {
		return 0, err
	}
	inputHdr := br.Header().VCF
	hdr := filterHeaderSamples(inputHdr, opts.Samples)
	if opts.recalcAC() {
		ensureACANHeader(hdr)
	}
	if opts.DropGenotypes {
		hdr = stripFormatLines(hdr)
		hdr.Samples = nil
	}
	includeF, excludeF, err := compileExpressions(opts, inputHdr)
	if err != nil {
		return 0, err
	}

	w, finish, err := openOutput(out, opts, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if !opts.NoHeader {
		if err := w.WriteHeader(); err != nil {
			return 0, err
		}
	}
	if opts.HeaderOnly {
		return 0, w.Flush()
	}

	count := 0
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, err
		}
		v := rec.ToVariant(br.Header())
		if !keepVariant(v, opts, includeF, excludeF, applyTargets, targets) {
			continue
		}
		applyACAN(v, opts)
		if opts.DropGenotypes {
			dropGenotypes(v)
		}
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// viewBCFRegions executes CSI-backed region queries on a BGZF-wrapped BCF.
func viewBCFRegions(path string, out io.Writer, opts ViewOptions, _ io.Writer) (int, error) {
	regions, err := parseRegions(opts.Regions)
	if err != nil {
		return 0, err
	}
	hdr, recs, err := ReadBCFRegions(path, regions)
	if err != nil {
		return 0, err
	}
	inputHdr := hdr.VCF
	vhdr := filterHeaderSamples(inputHdr, opts.Samples)
	if opts.recalcAC() {
		ensureACANHeader(vhdr)
	}
	if opts.DropGenotypes {
		vhdr = stripFormatLines(vhdr)
		vhdr.Samples = nil
	}

	includeF, excludeF, err := compileExpressions(opts, inputHdr)
	if err != nil {
		return 0, err
	}

	w, finish, err := openOutput(out, opts, vhdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if !opts.NoHeader {
		if err := w.WriteHeader(); err != nil {
			return 0, err
		}
	}
	if opts.HeaderOnly {
		return 0, w.Flush()
	}

	count := 0
	for _, rec := range recs {
		v := rec.ToVariant(hdr)
		if !keepVariant(v, opts, includeF, excludeF, true, regions) {
			continue
		}
		applyACAN(v, opts)
		if opts.DropGenotypes {
			dropGenotypes(v)
		}
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// viewRegions executes index-backed region queries on a bgzipped VCF.
func viewRegions(path string, out io.Writer, opts ViewOptions, stderr io.Writer) (int, error) {
	idx, err := readTabixIndex(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools view: load .tbi: %w", err)
	}
	regions, err := parseRegions(opts.Regions)
	if err != nil {
		return 0, err
	}

	// Open the bgzipped data file once for the whole batch of region queries.
	// openSeekable returns a *os.File for local paths and a ranged remote
	// handle (http(s)/s3/gs) wrapped in an io.SectionReader for URLs, so the
	// same index-driven seek path serves both.
	data, err := openSeekable(path)
	if err != nil {
		return 0, err
	}
	defer data.Close()

	// Read the header through iohelper for the metadata. We use a separate
	// stream for that because the tabix queries deliver record bytes only.
	hdrIn, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	hdrReader := vcf.NewReader(hdrIn)
	hdr, err := hdrReader.ReadHeader()
	if err != nil {
		hdrIn.Close()
		return 0, err
	}
	hdrIn.Close()
	inputHdr := hdr
	hdr = filterHeaderSamples(hdr, opts.Samples)
	if opts.recalcAC() {
		ensureACANHeader(hdr)
	}
	if opts.DropGenotypes {
		hdr = stripFormatLines(hdr)
		hdr.Samples = nil
	}

	includeF, excludeF, err := compileExpressions(opts, inputHdr)
	if err != nil {
		return 0, err
	}

	w, finish, err := openOutput(out, opts, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if !opts.NoHeader {
		if err := w.WriteHeader(); err != nil {
			return 0, err
		}
	}
	if opts.HeaderOnly {
		return 0, w.Flush()
	}

	// For each region: query bytes, parse each line through vcf.NewReader,
	// then push through the keep/filter pipeline.
	count := 0
	for _, reg := range regions {
		// Tabix queries use 0-based half-open; VCF specs are 1-based inclusive.
		// Translate accordingly.
		beg := reg.beg - 1
		if beg < 0 {
			beg = 0
		}
		end := reg.end
		lines, qerr := idx.QueryBytesReader(data, reg.chrom, beg, end)
		if qerr != nil {
			return count, qerr
		}
		for _, line := range lines {
			v, perr := parseVCFLine(line, hdr)
			if perr != nil {
				if stderr != nil {
					fmt.Fprintf(stderr, "bcftools view: skipping bad record: %v\n", perr)
				}
				continue
			}
			if !keepVariant(v, opts, includeF, excludeF, true, []region{reg}) {
				continue
			}
			applyACAN(v, opts)
			if opts.DropGenotypes {
				dropGenotypes(v)
			}
			if err := w.Write(v); err != nil {
				return count, err
			}
			count++
		}
	}
	return count, w.Flush()
}

// parseVCFLine wraps the existing vcf.Reader so we can decode one isolated
// record line (no header parsing needed). It mirrors the layout used by
// tabix queries (raw line bytes).
func parseVCFLine(line []byte, hdr *vcf.Header) (*vcf.Variant, error) {
	// Splice the line into a minimal stream with the header attached so the
	// vcf.Reader can re-use its existing logic. Building a header on every
	// line is cheap because the vcf.Reader doesn't allocate per record.
	var buf strings.Builder
	for _, m := range hdr.MetaInfo {
		buf.WriteString(m)
		buf.WriteByte('\n')
	}
	buf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	if len(hdr.Samples) > 0 {
		buf.WriteString("\tFORMAT\t")
		buf.WriteString(strings.Join(hdr.Samples, "\t"))
	}
	buf.WriteByte('\n')
	buf.Write(line)
	buf.WriteByte('\n')
	r := vcf.NewReader(strings.NewReader(buf.String()))
	if _, err := r.ReadHeader(); err != nil {
		return nil, err
	}
	return r.Read()
}

// keepVariant runs every filter in opts. Returns true when the variant
// should be emitted.
func keepVariant(v *vcf.Variant, opts ViewOptions, includeF, excludeF *Filter, applyTargets bool, targets []region) bool {
	if applyTargets && len(targets) > 0 && !overlapsAny(v, targets) {
		return false
	}
	if !opts.applyFilterColumnFilters(v) {
		return false
	}
	if !opts.applyAlleleFilters(v) {
		return false
	}
	if !passesPrivateFilter(v, opts) {
		return false
	}
	if !opts.passesTypeFilter(v) {
		return false
	}
	if includeF != nil && !includeF.Eval(v) {
		return false
	}
	if excludeF != nil && excludeF.Eval(v) {
		return false
	}
	return true
}

// compileExpressions parses the -i / -e flags into Filter trees, resolving
// bare identifiers against hdr so a FORMAT-only tag (e.g. `GQ>30`) is treated
// as FORMAT and a tag declared as both INFO and FORMAT is rejected as
// ambiguous, matching upstream filter.c. hdr should be the input header (before
// any sample restriction), as upstream resolves against the input header.
func compileExpressions(opts ViewOptions, hdr *vcf.Header) (include, exclude *Filter, err error) {
	if opts.IncludeExpr != "" {
		include, err = CompileFilterWithHeader(opts.IncludeExpr, hdr)
		if err != nil {
			return nil, nil, err
		}
	}
	if opts.ExcludeExpr != "" {
		exclude, err = CompileFilterWithHeader(opts.ExcludeExpr, hdr)
		if err != nil {
			return nil, nil, err
		}
	}
	return include, exclude, nil
}

// filterHeaderSamples returns a copy of hdr with Samples narrowed to the
// requested set (in the requested order). If samples is empty the header
// is returned unchanged.
func filterHeaderSamples(hdr *vcf.Header, samples []string) *vcf.Header {
	if len(samples) == 0 {
		return hdr
	}
	out := &vcf.Header{MetaInfo: append([]string{}, hdr.MetaInfo...)}
	seen := make(map[string]bool, len(hdr.Samples))
	for _, s := range hdr.Samples {
		seen[s] = true
	}
	for _, want := range samples {
		if seen[want] {
			out.Samples = append(out.Samples, want)
		}
	}
	return out
}

// stripFormatLines removes every ##FORMAT=... meta line from hdr. Upstream
// bcftools drops these lines when -G/--drop-genotypes is given because the
// resulting record has no FORMAT column to describe.
func stripFormatLines(hdr *vcf.Header) *vcf.Header {
	if hdr == nil {
		return hdr
	}
	out := &vcf.Header{Samples: hdr.Samples}
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##FORMAT=") {
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
	}
	return out
}

// variantWriter is the small interface View uses to send variants downstream.
// vcfWriter wraps *vcf.Writer; bcfWriter wraps *bcf.Writer; both speak the
// same WriteHeader / Write / Flush dance.
type variantWriter interface {
	WriteHeader() error
	Write(*vcf.Variant) error
	Flush() error
}

// bgzfFlusher is the minimal interface a vcf/bcfVariantWriter needs from the
// underlying BGZF writer for the header-block flush. It is satisfied by
// *bgzf.Writer (serial block compression) and by *bgzf.MultiWriter (parallel
// block compression, returned when -@/--threads > 1). Salvaged from PR #219 so
// the view/call/mpileup -Ob/-Oz output paths can close the BGZF block after the
// header, matching upstream htslib's vcf_hdr_write / bcf_hdr_write semantics.
type bgzfFlusher interface {
	Flush() error
	Close() error
	Write(p []byte) (int, error)
}

// vcfVariantWriter is the trivial pass-through adapter for the VCF / VCF.gz
// output paths. When bgzf is non-nil (the -Oz path), WriteHeader closes the
// BGZF block after the header so the header occupies its own block, matching
// upstream bcftools / htslib's vcf_hdr_write which calls bgzf_flush after the
// header.
type vcfVariantWriter struct {
	w    *vcf.Writer
	bgzf bgzfFlusher
}

func (a *vcfVariantWriter) WriteHeader() error {
	if err := a.w.WriteHeader(); err != nil {
		return err
	}
	if a.bgzf != nil {
		if err := a.w.Flush(); err != nil {
			return err
		}
		return a.bgzf.Flush()
	}
	return nil
}
func (a *vcfVariantWriter) Write(v *vcf.Variant) error { return a.w.Write(v) }
func (a *vcfVariantWriter) Flush() error               { return a.w.Flush() }

// bcfVariantWriter wraps a bcf.Writer so View can treat both output formats
// the same way. As with vcfVariantWriter, a non-nil bgzf closes the BGZF block
// after the header for the -Ob / -Ou paths.
type bcfVariantWriter struct {
	w    *bcf.Writer
	bgzf bgzfFlusher
}

func (a *bcfVariantWriter) WriteHeader() error {
	if err := a.w.WriteHeader(); err != nil {
		return err
	}
	if a.bgzf != nil {
		if err := a.w.Flush(); err != nil {
			return err
		}
		return a.bgzf.Flush()
	}
	return nil
}
func (a *bcfVariantWriter) Write(v *vcf.Variant) error { return a.w.Write(v) }
func (a *bcfVariantWriter) Flush() error               { return a.w.Flush() }

// ensurePASSFilter inserts ##FILTER=<ID=PASS,...> immediately after the
// ##fileformat line (or at the top of the header if no fileformat line is
// present) when the header lacks an explicit PASS definition. This mirrors
// upstream htslib's bcf_hdr_parse_line behaviour and is required for
// byte-for-byte parity with bcftools text output.
func ensurePASSFilter(hdr *vcf.Header) {
	if hdr == nil {
		return
	}
	for _, m := range hdr.MetaInfo {
		k, id := structuredID(m)
		if k == "FILTER" && id == "PASS" {
			return
		}
	}
	passLine := `##FILTER=<ID=PASS,Description="All filters passed">`
	insertAt := 0
	for i, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##fileformat=") {
			insertAt = i + 1
			break
		}
	}
	out := make([]string, 0, len(hdr.MetaInfo)+1)
	out = append(out, hdr.MetaInfo[:insertAt]...)
	out = append(out, passLine)
	out = append(out, hdr.MetaInfo[insertAt:]...)
	hdr.MetaInfo = out
}

// openOutput returns a variantWriter plus a cleanup function. The cleanup
// closes any wrapping compressor that needs an explicit Close (BGZF for -O z
// and -O b). For the BGZF paths the variantWriter is given the bgzf flusher so
// WriteHeader closes the header into its own BGZF block, matching upstream
// htslib's vcf_hdr_write / bcf_hdr_write (which call bgzf_flush after the
// header) and keeping tabix/.csi offsets clean. We deliberately do not close
// `out` itself — the caller still owns it.
func openOutput(out io.Writer, opts ViewOptions, hdr *vcf.Header) (variantWriter, func(), error) {
	if !opts.SuppressPASSFilter {
		ensurePASSFilter(hdr)
	}
	switch opts.OutputFormat {
	case OutputVCFGz:
		bw, err := newBGZFOutput(out, opts.CompressLevel, opts.Threads)
		if err != nil {
			return nil, func() {}, err
		}
		return &vcfVariantWriter{w: vcf.NewWriter(bw, hdr), bgzf: bw}, func() { _ = bw.Close() }, nil
	case OutputBCF:
		bw, err := newBGZFOutput(out, opts.CompressLevel, opts.Threads)
		if err != nil {
			return nil, func() {}, err
		}
		w, err := bcf.NewWriterFromVCFHeader(bw, hdr)
		if err != nil {
			_ = bw.Close()
			return nil, func() {}, err
		}
		return &bcfVariantWriter{w: w, bgzf: bw}, func() { _ = w.Flush(); _ = bw.Close() }, nil
	case OutputBCFUncompressed:
		w, err := bcf.NewWriterFromVCFHeader(out, hdr)
		if err != nil {
			return nil, func() {}, err
		}
		return &bcfVariantWriter{w: w}, func() { _ = w.Flush() }, nil
	}
	return &vcfVariantWriter{w: vcf.NewWriter(out, hdr)}, func() {}, nil
}

// LoadSamplesFilePairs reads a samples file like LoadSamplesFile but
// also returns the optional second whitespace-separated column as a
// rename target. Each result entry is {name, alias}; alias is "" when
// the line carries only one field. Mirrors bam_smpl_add_samples
// (bam_sample.c) which threads the second column as a per-sample
// rename map. Lines starting with `#` and blank lines are ignored.
// Salvaged from PR #219 for the mpileup -S/--samples-file path.
func LoadSamplesFilePairs(path string) ([][2]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var pairs [][2]string
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		var name, alias string
		if i := strings.IndexAny(line, "\t "); i >= 0 {
			name = line[:i]
			alias = strings.TrimSpace(line[i:])
		} else {
			name = line
		}
		if name == "" {
			continue
		}
		pairs = append(pairs, [2]string{name, alias})
	}
	return pairs, sc.Err()
}

// newBGZFOutput returns a bgzfFlusher (a BGZF io.WriteCloser that also exposes
// Flush) for compressed VCF.gz (-O z) and BCF (-O b) output. It mirrors
// upstream bcftools, which writes BGZF (a gzip-compatible block format) for
// both. A level < 0 selects the package default compression level.
//
// When threads > 1 the returned writer is a bgzf.MultiWriter that compresses
// blocks across that many worker goroutines; otherwise it is the
// single-threaded bgzf.Writer. Because every BGZF block is an independent gzip
// member, both paths produce a stream that decodes to byte-identical plaintext
// regardless of the thread count. The caller must Close the returned writer to
// flush the final block and emit the BGZF EOF marker.
//
// The concrete result (*bgzf.Writer or *bgzf.MultiWriter) satisfies
// bgzfFlusher, so callers can wire it into a vcf/bcfVariantWriter to flush the
// header block (matching upstream htslib's vcf_hdr_write / bcf_hdr_write).
func newBGZFOutput(out io.Writer, level, threads int) (bgzfFlusher, error) {
	if level < 0 {
		level = bgzip.DefaultCompression
	}
	if threads > 1 {
		return bgzip.NewMultiWriter(out, level, threads)
	}
	return bgzip.NewWriterLevel(out, level)
}

// LoadSamplesFile reads one sample name per line. Blank lines and comment
// lines beginning with '#' are skipped.
func LoadSamplesFile(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var names []string
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		// `samples-file` may include extra tab-separated fields (e.g. for
		// `bcftools view -S`); the sample ID is always the first column.
		if i := strings.IndexAny(line, "\t "); i >= 0 {
			line = line[:i]
		}
		names = append(names, line)
	}
	return names, sc.Err()
}

// LoadRegionsFile reads BED-style regions (CHROM \t BEG \t END) one per line,
// converting them to 1-based inclusive ranges that match VCF spec text.
func LoadRegionsFile(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var out []string
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) == 1 {
			out = append(out, fields[0])
			continue
		}
		if len(fields) < 3 {
			return nil, fmt.Errorf("bcftools view: bad regions-file line %q", line)
		}
		beg, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("bcftools view: bad start in regions-file: %w", err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("bcftools view: bad end in regions-file: %w", err)
		}
		// BED is 0-based half-open; VCF/region spec is 1-based inclusive.
		out = append(out, fmt.Sprintf("%s:%d-%d", fields[0], beg+1, end))
	}
	return out, sc.Err()
}

// SplitCommaList splits "a,b,c" honoring trailing/leading whitespace.
func SplitCommaList(s string) []string {
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	for i := range parts {
		parts[i] = strings.TrimSpace(parts[i])
	}
	return parts
}
