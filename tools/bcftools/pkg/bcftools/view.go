package bcftools

import (
	"bufio"
	"compress/flate"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// OutputFormat describes how `bcftools view` should emit records.
type OutputFormat int

const (
	// OutputVCF is uncompressed VCF text (the default).
	OutputVCF OutputFormat = iota
	// OutputVCFGz is gzip-compressed VCF.
	OutputVCFGz
	// OutputBCF marks BGZF-compressed BCF output (-Ob): BCF bytes streamed
	// through a bgzf.Writer.
	OutputBCF
	// OutputBCFUncompressed marks "uncompressed" BCF output (-Ou). This
	// mirrors htslib's "wbu" open mode: the BCF byte stream is still wrapped
	// in BGZF framing, but each deflate block is stored uncompressed
	// (compression level 0). It is therefore BGZF-readable like -Ob, just
	// without the deflate cost — not a raw, frameless BCF.
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
	CompressLevel  int      // -l
	// IncludeTypes / ExcludeTypes are the comma lists driving
	// `-v/--types` and `-V/--exclude-types`. Accepted token set:
	// snps, indels, mnps, ref, bnd, other. Both empty means no filter;
	// supplying both is an error (matches vcfview.c:170).
	IncludeTypes []string
	ExcludeTypes []string
	// NoUpdateINFO is the `-I/--no-update` flag (vcfview.c:567+669).
	// When true the per-record INFO/AC and INFO/AN are NOT recomputed
	// after a sample subset; the original values pass through. Default
	// (false) matches upstream's behaviour of always recomputing when
	// samples are restricted.
	NoUpdateINFO bool
	// SkipPASSInjection disables the automatic ##FILTER=<ID=PASS,...>
	// header insertion. Internal use only: `bcftools reheader` is a
	// pure header rewrite — upstream does not synthesise a PASS line
	// when one isn't already present, so neither should we.
	SkipPASSInjection bool
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

// recomputeACAN walks the per-sample GT field and rewrites v.Info["AC"]
// and v.Info["AN"] from scratch. AC is the per-ALT non-reference count
// (so a multi-allelic site emits a comma-separated list of length
// len(v.Alt)); AN is the total called-allele count. Mirrors the path
// upstream vcfview.c takes when `args->update_info` is on and samples
// have been subset (vcfview.c:355-364, 452-454).
//
// If no sample carries a GT field the function leaves Info untouched
// (mirroring vcfview.c:354 `if ( ... && !bcf_get_fmt(... "GT")) update_ac = 0`).
func recomputeACAN(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}
	hasGT := false
	for _, s := range v.Samples {
		if _, ok := s.Data["GT"]; ok {
			hasGT = true
			break
		}
	}
	if !hasGT {
		return
	}
	ac := make([]int, len(v.Alt))
	an := 0
	for _, s := range v.Samples {
		gt, ok := s.Data["GT"]
		if !ok {
			continue
		}
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			if a == "." || a == "" {
				continue
			}
			an++
			n, err := strconv.Atoi(a)
			if err != nil || n <= 0 {
				continue
			}
			if n-1 < len(ac) {
				ac[n-1]++
			}
		}
	}
	if v.Info == nil {
		v.Info = make(map[string]string)
	}
	// AC may not have been declared in the source INFO map at all. We
	// still update it; the header line is expected to exist (upstream
	// also writes the tag unconditionally and relies on the header
	// declaration). InfoOrder is touched only if needed so the source
	// key order is otherwise preserved.
	acParts := make([]string, len(ac))
	for i, n := range ac {
		acParts[i] = strconv.Itoa(n)
	}
	acStr := strings.Join(acParts, ",")
	if _, present := v.Info["AC"]; !present {
		v.InfoOrder = append(v.InfoOrder, "AC")
	}
	v.Info["AC"] = acStr
	if _, present := v.Info["AN"]; !present {
		v.InfoOrder = append(v.InfoOrder, "AN")
	}
	v.Info["AN"] = strconv.Itoa(an)
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

// hasTabixIndex returns true if path has a sibling .tbi file.
func hasTabixIndex(path string) bool {
	_, err := os.Stat(path + ".tbi")
	return err == nil
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
		if len(opts.Samples) > 0 {
			restrictSamples(v, opts.Samples)
			if !opts.NoUpdateINFO {
				recomputeACAN(v)
			}
		}
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

// viewBCFStream consumes a BCF stream and re-emits records in the
// requested output format (VCF, gzipped VCF, raw BCF, or BGZF BCF) via
// openOutput.
func viewBCFStream(in io.Reader, out io.Writer, opts ViewOptions, applyTargets bool, targets []region) (int, error) {
	br, err := bcf.NewReader(in)
	if err != nil {
		return 0, err
	}
	origHdr := br.Header().VCF
	hdr := filterHeaderSamples(origHdr, opts.Samples)
	if opts.DropGenotypes {
		hdr = stripFormatLines(hdr)
		hdr.Samples = nil
	}
	includeF, excludeF, err := compileExpressions(opts, origHdr)
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
		if len(opts.Samples) > 0 {
			restrictSamples(v, opts.Samples)
			if !opts.NoUpdateINFO {
				recomputeACAN(v)
			}
		}
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
	origHdr := hdr.VCF
	vhdr := filterHeaderSamples(origHdr, opts.Samples)
	if opts.DropGenotypes {
		vhdr = stripFormatLines(vhdr)
		vhdr.Samples = nil
	}

	includeF, excludeF, err := compileExpressions(opts, origHdr)
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
		if len(opts.Samples) > 0 {
			restrictSamples(v, opts.Samples)
			if !opts.NoUpdateINFO {
				recomputeACAN(v)
			}
		}
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
	idx, err := tabix.ReadFile(path + ".tbi")
	if err != nil {
		return 0, fmt.Errorf("bcftools view: load .tbi: %w", err)
	}
	regions, err := parseRegions(opts.Regions)
	if err != nil {
		return 0, err
	}

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
	origHdr := hdr
	hdr = filterHeaderSamples(hdr, opts.Samples)
	if opts.DropGenotypes {
		hdr = stripFormatLines(hdr)
		hdr.Samples = nil
	}

	includeF, excludeF, err := compileExpressions(opts, origHdr)
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
		lines, qerr := idx.QueryBytes(path, reg.chrom, beg, end)
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
			if len(opts.Samples) > 0 {
				restrictSamples(v, opts.Samples)
				if !opts.NoUpdateINFO {
					recomputeACAN(v)
				}
			}
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
	if !opts.applyVariantTypeFilters(v) {
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

// applyVariantTypeFilters returns true if v passes the -v/-V type
// selectors. A variant is kept when at least one of its per-ALT types
// is in the include set, and none of its per-ALT types are in the
// exclude set. Mirrors vcfview.c:325-329 — but operates on the
// already-decoded ALT strings rather than the bcf bitmask.
func (o ViewOptions) applyVariantTypeFilters(v *vcf.Variant) bool {
	if len(o.IncludeTypes) == 0 && len(o.ExcludeTypes) == 0 {
		return true
	}
	types := variantTypesPerALT(v)
	if len(o.IncludeTypes) > 0 {
		ok := false
		for _, t := range o.IncludeTypes {
			for _, vt := range types {
				if vt == t {
					ok = true
					break
				}
			}
			if ok {
				break
			}
		}
		if !ok {
			return false
		}
	}
	if len(o.ExcludeTypes) > 0 {
		for _, t := range o.ExcludeTypes {
			for _, vt := range types {
				if vt == t {
					return false
				}
			}
		}
	}
	return true
}

// variantTypesPerALT classifies each ALT allele of v into one of the
// upstream bucket names: "snps", "indels", "mnps", "bnd", "ref",
// "other". A site with N ALTs returns N labels. Mirrors htslib's
// bcf_set_variant_type (vcf.c:5380-5439) — REF-length, ALT-length and
// breakend bracket syntax drive the bucket choice.
func variantTypesPerALT(v *vcf.Variant) []string {
	if len(v.Alt) == 0 {
		return []string{"ref"}
	}
	out := make([]string, 0, len(v.Alt))
	refLen := len(v.Ref)
	for _, a := range v.Alt {
		switch {
		case a == "" || a == "." || a == "*":
			out = append(out, "other")
		case strings.ContainsAny(a, "[]"):
			out = append(out, "bnd")
		case strings.HasPrefix(a, "<") && strings.HasSuffix(a, ">"):
			out = append(out, "other")
		case len(a) == 1 && refLen == 1:
			out = append(out, "snps")
		case len(a) == refLen:
			out = append(out, "mnps")
		default:
			out = append(out, "indels")
		}
	}
	return out
}

// compileExpressions parses the -i / -e flags into Filter trees, resolving
// bare tags against hdr (which should be the unmodified input header so INFO
// vs FORMAT resolution matches htslib). A nil hdr falls back to the
// header-less observable resolution rule.
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

// vcfVariantWriter is the trivial pass-through adapter for the VCF / VCF.gz
// output paths. When bgzf is non-nil (the -Oz path), WriteHeader closes the
// BGZF block after the header so the header occupies its own block, matching
// upstream bcftools / htslib's vcf_hdr_write which calls bgzf_flush after the
// header.
type vcfVariantWriter struct {
	w    *vcf.Writer
	bgzf *bgzip.Writer
}

func (a *vcfVariantWriter) WriteHeader() error {
	if err := a.w.WriteHeader(); err != nil {
		return err
	}
	if a.bgzf != nil {
		// Push the header bytes out of the vcf.Writer's bufio buffer so they
		// reach the BGZF layer, then close the BGZF block.
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
	bgzf *bgzip.Writer
}

func (a *bcfVariantWriter) WriteHeader() error {
	if err := a.w.WriteHeader(); err != nil {
		return err
	}
	if a.bgzf != nil {
		// Push the header bytes out of the bcf.Writer's bufio buffer so they
		// reach the BGZF layer, then close the BGZF block.
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
// closes any wrapping compressor that needs an explicit Close (gzip for -O z,
// bgzip for -O b). We deliberately do not close `out` itself — the caller
// still owns it.
func openOutput(out io.Writer, opts ViewOptions, hdr *vcf.Header) (variantWriter, func(), error) {
	if !opts.SkipPASSInjection {
		ensurePASSFilter(hdr)
	}
	switch opts.OutputFormat {
	case OutputVCFGz:
		// Upstream -Oz emits BGZF (block-gzip with the BC subfield), not
		// plain gzip, so the output is tabix-indexable and byte-matches
		// `bgzip`. Use the BGZF writer rather than compress/gzip.
		var gw *bgzip.Writer
		if opts.CompressLevel < 0 {
			gw = bgzip.NewWriter(out)
		} else {
			var err error
			gw, err = bgzip.NewWriterLevel(out, opts.CompressLevel)
			if err != nil {
				return nil, func() {}, err
			}
		}
		return &vcfVariantWriter{w: vcf.NewWriter(gw, hdr), bgzf: gw}, func() { _ = gw.Close() }, nil
	case OutputBCF:
		bw := bgzip.NewWriter(out)
		w, err := bcf.NewWriterFromVCFHeader(bw, hdr)
		if err != nil {
			_ = bw.Close()
			return nil, func() {}, err
		}
		return &bcfVariantWriter{w: w, bgzf: bw}, func() { _ = w.Flush(); _ = bw.Close() }, nil
	case OutputBCFUncompressed:
		// htslib's "wbu" mode wraps BCF in BGZF but stores each block
		// uncompressed (deflate level 0). Match that framing so the output
		// is structurally identical to genuine `bcftools view -Ou`.
		bw, err := bgzip.NewWriterLevel(out, flate.NoCompression)
		if err != nil {
			return nil, func() {}, err
		}
		w, err := bcf.NewWriterFromVCFHeader(bw, hdr)
		if err != nil {
			_ = bw.Close()
			return nil, func() {}, err
		}
		return &bcfVariantWriter{w: w, bgzf: bw}, func() { _ = w.Flush(); _ = bw.Close() }, nil
	}
	return &vcfVariantWriter{w: vcf.NewWriter(out, hdr)}, func() {}, nil
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

// LoadSamplesFilePairs reads a samples file like LoadSamplesFile but
// also returns the optional second whitespace-separated column as a
// rename target. Each result entry is {name, alias}; alias is "" when
// the line carries only one field. Mirrors bam_smpl_add_samples
// (bam_sample.c) which threads the second column as a per-sample
// rename map. Lines starting with `#` and blank lines are ignored.
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
