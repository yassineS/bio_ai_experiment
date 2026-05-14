package bcftools

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
	"github.com/yassineS/bio_ai_experiment/tools/tabix/pkg/tabix"
)

// OutputFormat describes how `bcftools view` should emit records.
type OutputFormat int

const (
	// OutputVCF is uncompressed VCF text (the default).
	OutputVCF OutputFormat = iota
	// OutputVCFGz is gzip-compressed VCF.
	OutputVCFGz
	// OutputBCF marks compressed BCF output. NOT IMPLEMENTED in this slice;
	// the runner returns an explanatory error so callers do not silently get
	// the wrong format.
	OutputBCF
	// OutputBCFUncompressed marks uncompressed BCF output. Also deferred.
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

// View streams VCF or BCF from in, applies opts, and writes to out. It is
// the streaming entry point used by the `view` subcommand when no region
// query is requested (regions need a seekable file, see ViewFile).
func View(in io.Reader, out io.Writer, opts ViewOptions) (int, error) {
	return viewStreaming(in, out, opts, false, nil)
}

// ViewFile is the file-aware entry point for `bcftools view`. It supports
// region queries via a sibling .tbi index on bgzipped VCF inputs and falls
// back to a streaming scan otherwise.
func ViewFile(path string, out io.Writer, opts ViewOptions, stderr io.Writer) (int, error) {
	if len(opts.Regions) > 0 && hasTabixIndex(path) {
		return viewRegions(path, out, opts, stderr)
	}
	if len(opts.Regions) > 0 && stderr != nil {
		fmt.Fprintln(stderr, "bcftools view: no .tbi index found; treating -r as a post-filter (slower)")
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

	includeF, excludeF, err := compileExpressions(opts)
	if err != nil {
		return 0, err
	}
	hdr = filterHeaderSamples(hdr, opts.Samples)
	if opts.DropGenotypes {
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

// viewBCFStream consumes a BCF stream and emits VCF text. Writing BCF is
// deferred — when callers ask for OutputBCF/OutputBCFUncompressed the
// outer runner has already rejected the request.
func viewBCFStream(in io.Reader, out io.Writer, opts ViewOptions, applyTargets bool, targets []region) (int, error) {
	br, err := bcf.NewReader(in)
	if err != nil {
		return 0, err
	}
	hdr := br.Header().VCF
	hdr = filterHeaderSamples(hdr, opts.Samples)
	if opts.DropGenotypes {
		hdr.Samples = nil
	}
	includeF, excludeF, err := compileExpressions(opts)
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
	hdr = filterHeaderSamples(hdr, opts.Samples)
	if opts.DropGenotypes {
		hdr.Samples = nil
	}

	includeF, excludeF, err := compileExpressions(opts)
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
	if includeF != nil && !includeF.Eval(v) {
		return false
	}
	if excludeF != nil && excludeF.Eval(v) {
		return false
	}
	return true
}

// compileExpressions parses the -i / -e flags into Filter trees.
func compileExpressions(opts ViewOptions) (include, exclude *Filter, err error) {
	if opts.IncludeExpr != "" {
		include, err = CompileFilter(opts.IncludeExpr)
		if err != nil {
			return nil, nil, err
		}
	}
	if opts.ExcludeExpr != "" {
		exclude, err = CompileFilter(opts.ExcludeExpr)
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

// openOutput returns a vcf.Writer plus a cleanup function. The cleanup
// closes the gzip writer if one was created. We deliberately do not close
// out itself — the caller still owns it.
func openOutput(out io.Writer, opts ViewOptions, hdr *vcf.Header) (*vcf.Writer, func(), error) {
	switch opts.OutputFormat {
	case OutputVCFGz:
		gw, err := gzipWriter(out, opts.CompressLevel)
		if err != nil {
			return nil, func() {}, err
		}
		return vcf.NewWriter(gw, hdr), func() { _ = gw.Close() }, nil
	case OutputBCF, OutputBCFUncompressed:
		return nil, func() {}, fmt.Errorf("bcftools view: -O b/u (BCF output) is not yet implemented; use -O v or -O z")
	}
	return vcf.NewWriter(out, hdr), func() {}, nil
}

// gzipWriter returns a gzip writer at the requested level (or default if
// level < 0).
func gzipWriter(out io.Writer, level int) (*gzip.Writer, error) {
	if level < 0 {
		return gzip.NewWriter(out), nil
	}
	return gzip.NewWriterLevel(out, level)
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
