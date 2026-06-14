package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// readTabixLines runs a tabix-backed region query and returns the raw record
// bytes. It is shared by the query and concat region paths. The .tbi index and
// the data file are opened through the hfile-aware helpers, so path may be a
// local file or a remote http(s)://, s3:// or gs:// URL.
func readTabixLines(path string, regs []region) ([][]byte, error) {
	idx, err := readTabixIndex(path)
	if err != nil {
		return nil, fmt.Errorf("load .tbi: %w", err)
	}
	src, err := openSeekable(path)
	if err != nil {
		return nil, err
	}
	defer src.Close()
	var out [][]byte
	for _, reg := range regs {
		beg := reg.beg - 1
		if beg < 0 {
			beg = 0
		}
		end := reg.end
		lines, qerr := idx.QueryBytesReader(src, reg.chrom, beg, end)
		if qerr != nil {
			return out, qerr
		}
		out = append(out, lines...)
	}
	return out, nil
}

// QueryOptions controls the behaviour of Query / QueryFile.
type QueryOptions struct {
	// Format is the bcftools-style format string. Tokens supported:
	//   %CHROM %POS %POS0 %END %END0 %REF %ALT %FIRST_ALT %QUAL %ID
	//   %FILTER %TYPE %IS_TS %TGT %GT
	//   %INFO        (bare)                      — entire INFO column
	//   %INFO/<TAG>                              — INFO field by name
	//   [%FIELD]                                 — sample-level field, tab-repeated
	//   \n \t                                    — newline / tab
	// Anything else is emitted verbatim.
	Format string
	// PrintHeader emits a "# " prefixed header row derived from the format
	// string (the `-H` / `--print-header` flag).
	PrintHeader bool
	// ListSamples short-circuits formatting and emits one sample name per
	// line (the `-l` / `--list-samples` flag).
	ListSamples bool
	// Samples narrows the per-sample output to these names, in order.
	Samples []string
	// SamplesFile is a sibling file that supplies sample names (one per
	// line). It is loaded and concatenated to Samples by the runner.
	SamplesFile string
	// Regions are index-backed region filters (`-r`).
	Regions []string
	// RegionsFile is the BED-style sibling for Regions (`-R`).
	RegionsFile string
	// Targets are post-filter region specs (`-t`), always applied after
	// reading regardless of indexing.
	Targets []string
	// TargetsFile is the BED-style sibling for Targets (`-T`).
	TargetsFile string
	// IncludeExpr / ExcludeExpr reuse the View filter evaluator.
	IncludeExpr string
	ExcludeExpr string
	// ApplyFilters is the FILTER name list to keep (`-F` / `--apply-filters`).
	// Upstream bcftools query maps `-F` here because `-f` is the format
	// string. An empty list disables the filter.
	ApplyFilters []string
}

// Query streams VCF or BCF from in, applies opts, and writes the formatted
// output to out. Records that do not pass the filters are silently skipped.
// The number of records written is returned.
func Query(in io.Reader, out io.Writer, opts QueryOptions) (int, error) {
	return queryStreaming(in, out, opts, false, nil)
}

// QueryFile is the file-aware entry point: it loads regions from a sibling
// CSI/TBI when available and otherwise falls back to a streaming scan with
// regions promoted to post-filters.
func QueryFile(path string, out io.Writer, opts QueryOptions, stderr io.Writer) (int, error) {
	if opts.ListSamples {
		return listSamplesFromFile(path, out)
	}
	if len(opts.Regions) > 0 && HasCSI(path) {
		bcfHead, _ := looksLikeBCF(path)
		if bcfHead {
			return queryBCFRegions(path, out, opts)
		}
	}
	if len(opts.Regions) > 0 && hasTabixIndex(path) {
		return queryVCFRegions(path, out, opts, stderr)
	}
	if len(opts.Regions) > 0 && stderr != nil {
		fmt.Fprintln(stderr, "bcftools query: no .tbi/.csi index found; treating -r as a post-filter (slower)")
	}
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	defer in.Close()
	postFilters := opts.Targets
	if len(opts.Regions) > 0 {
		postFilters = append([]string{}, opts.Targets...)
		postFilters = append(postFilters, opts.Regions...)
	}
	targets, err := parseRegions(postFilters)
	if err != nil {
		return 0, err
	}
	return queryStreaming(in, out, opts, len(targets) > 0, targets)
}

// listSamplesFromFile prints one sample name per line and is the
// implementation behind `-l`.
func listSamplesFromFile(path string, out io.Writer) (int, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	defer in.Close()
	samples, err := readHeaderSamples(in)
	if err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(out)
	for _, s := range samples {
		if _, err := fmt.Fprintln(bw, s); err != nil {
			return 0, err
		}
	}
	return len(samples), bw.Flush()
}

// readHeaderSamples reads enough of in to extract the sample list — either
// from a BCF text header or a VCF #CHROM line.
func readHeaderSamples(in io.Reader) ([]string, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return nil, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		r, err := bcf.NewReader(br)
		if err != nil {
			return nil, err
		}
		return append([]string{}, r.Header().VCF.Samples...), nil
	}
	r := vcf.NewReader(br)
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, err
	}
	return append([]string{}, hdr.Samples...), nil
}

// queryStreaming dispatches by magic to the VCF or BCF reader and drives the
// format-string formatter.
func queryStreaming(in io.Reader, out io.Writer, opts QueryOptions, applyTargets bool, targets []region) (int, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return 0, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return queryBCFStream(br, out, opts, applyTargets, targets)
	}
	return queryVCFStream(br, out, opts, applyTargets, targets)
}

// queryVCFStream drives the VCF (text or gzip-decoded) path.
func queryVCFStream(in io.Reader, out io.Writer, opts QueryOptions, applyTargets bool, targets []region) (int, error) {
	r := vcf.NewReader(in)
	hdr, err := r.ReadHeader()
	if err != nil {
		return 0, err
	}
	if opts.ListSamples {
		return writeSampleList(out, hdr.Samples)
	}
	tokens, err := ParseFormatString(opts.Format)
	if err != nil {
		return 0, err
	}
	includeF, excludeF, err := compileQueryExpressions(opts)
	if err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()
	if opts.PrintHeader {
		if err := writeHeaderRow(bw, tokens, opts.Samples, hdr.Samples); err != nil {
			return 0, err
		}
	}
	sampleFilter := buildSampleFilter(opts.Samples, hdr.Samples)
	count := 0
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, err
		}
		if !keepQueryVariant(v, opts, includeF, excludeF, applyTargets, targets) {
			continue
		}
		if err := emitRecord(bw, tokens, v, sampleFilter, hdr.Samples); err != nil {
			return count, err
		}
		count++
	}
	return count, nil
}

// queryBCFStream is the BCF mirror of queryVCFStream.
func queryBCFStream(in io.Reader, out io.Writer, opts QueryOptions, applyTargets bool, targets []region) (int, error) {
	r, err := bcf.NewReader(in)
	if err != nil {
		return 0, err
	}
	hdr := r.Header()
	if opts.ListSamples {
		return writeSampleList(out, hdr.VCF.Samples)
	}
	tokens, err := ParseFormatString(opts.Format)
	if err != nil {
		return 0, err
	}
	includeF, excludeF, err := compileQueryExpressions(opts)
	if err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()
	if opts.PrintHeader {
		if err := writeHeaderRow(bw, tokens, opts.Samples, hdr.VCF.Samples); err != nil {
			return 0, err
		}
	}
	sampleFilter := buildSampleFilter(opts.Samples, hdr.VCF.Samples)
	count := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, err
		}
		v := rec.ToVariant(hdr)
		if !keepQueryVariant(v, opts, includeF, excludeF, applyTargets, targets) {
			continue
		}
		if err := emitRecord(bw, tokens, v, sampleFilter, hdr.VCF.Samples); err != nil {
			return count, err
		}
		count++
	}
	return count, nil
}

// queryBCFRegions handles CSI-backed region queries on BCF inputs.
func queryBCFRegions(path string, out io.Writer, opts QueryOptions) (int, error) {
	regs, err := parseRegions(opts.Regions)
	if err != nil {
		return 0, err
	}
	hdr, recs, err := ReadBCFRegions(path, regs)
	if err != nil {
		return 0, err
	}
	tokens, err := ParseFormatString(opts.Format)
	if err != nil {
		return 0, err
	}
	includeF, excludeF, err := compileQueryExpressions(opts)
	if err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()
	if opts.PrintHeader {
		if err := writeHeaderRow(bw, tokens, opts.Samples, hdr.VCF.Samples); err != nil {
			return 0, err
		}
	}
	sampleFilter := buildSampleFilter(opts.Samples, hdr.VCF.Samples)
	count := 0
	for _, rec := range recs {
		v := rec.ToVariant(hdr)
		if !keepQueryVariant(v, opts, includeF, excludeF, true, regs) {
			continue
		}
		if err := emitRecord(bw, tokens, v, sampleFilter, hdr.VCF.Samples); err != nil {
			return count, err
		}
		count++
	}
	return count, nil
}

// queryVCFRegions handles TBI-backed region queries on bgzipped VCF inputs.
// It reuses the same parseVCFLine helper as `view`.
func queryVCFRegions(path string, out io.Writer, opts QueryOptions, stderr io.Writer) (int, error) {
	// Reuse the existing viewRegions plumbing by capturing its output through
	// an intermediary buffer is overkill; mirror the index logic here.
	hdrIn, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	r := vcf.NewReader(hdrIn)
	hdr, err := r.ReadHeader()
	if err != nil {
		hdrIn.Close()
		return 0, err
	}
	hdrIn.Close()
	tokens, err := ParseFormatString(opts.Format)
	if err != nil {
		return 0, err
	}
	includeF, excludeF, err := compileQueryExpressions(opts)
	if err != nil {
		return 0, err
	}
	regs, err := parseRegions(opts.Regions)
	if err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()
	if opts.PrintHeader {
		if err := writeHeaderRow(bw, tokens, opts.Samples, hdr.Samples); err != nil {
			return 0, err
		}
	}
	sampleFilter := buildSampleFilter(opts.Samples, hdr.Samples)
	// We delegate to the same path used by view for the tabix indexed read.
	count, err := queryTBI(path, hdr, regs, opts, includeF, excludeF, sampleFilter, tokens, bw, stderr)
	return count, err
}

// queryTBI is a thin wrapper around the shared tabix indexing logic.
func queryTBI(path string, hdr *vcf.Header, regs []region, opts QueryOptions, includeF, excludeF *Filter, sampleFilter []int, tokens []FormatToken, w io.Writer, stderr io.Writer) (int, error) {
	lines, err := readTabixLines(path, regs)
	if err != nil {
		return 0, err
	}
	count := 0
	for _, line := range lines {
		v, perr := parseVCFLine(line, hdr)
		if perr != nil {
			if stderr != nil {
				fmt.Fprintf(stderr, "bcftools query: skipping bad record: %v\n", perr)
			}
			continue
		}
		if !keepQueryVariant(v, opts, includeF, excludeF, true, regs) {
			continue
		}
		if err := emitRecord(w, tokens, v, sampleFilter, hdr.Samples); err != nil {
			return count, err
		}
		count++
	}
	return count, nil
}

// writeSampleList emits one sample name per line.
func writeSampleList(out io.Writer, names []string) (int, error) {
	bw := bufio.NewWriter(out)
	for _, s := range names {
		if _, err := fmt.Fprintln(bw, s); err != nil {
			return 0, err
		}
	}
	return len(names), bw.Flush()
}

// compileQueryExpressions parses -i/-e using the shared expression compiler.
func compileQueryExpressions(opts QueryOptions) (include, exclude *Filter, err error) {
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

// keepQueryVariant evaluates the filter pipeline for query. The set of
// filters is intentionally a subset of view's: targets, ApplyFilters, -i/-e.
func keepQueryVariant(v *vcf.Variant, opts QueryOptions, includeF, excludeF *Filter, applyTargets bool, targets []region) bool {
	if applyTargets && len(targets) > 0 && !overlapsAny(v, targets) {
		return false
	}
	if len(opts.ApplyFilters) > 0 {
		ok := false
		for _, want := range opts.ApplyFilters {
			for _, f := range v.Filter {
				if f == want {
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
	if includeF != nil && !includeF.Eval(v) {
		return false
	}
	if excludeF != nil && excludeF.Eval(v) {
		return false
	}
	return true
}

// buildSampleFilter returns a list of indexes into v.Samples to emit, in the
// requested order. Names not present in the header are silently skipped. An
// empty selection means "all samples in header order" and is signalled with a
// nil slice.
func buildSampleFilter(wanted, headerSamples []string) []int {
	if len(wanted) == 0 {
		return nil
	}
	pos := make(map[string]int, len(headerSamples))
	for i, s := range headerSamples {
		pos[s] = i
	}
	out := make([]int, 0, len(wanted))
	for _, w := range wanted {
		if idx, ok := pos[w]; ok {
			out = append(out, idx)
		}
	}
	return out
}

// FormatToken is one piece of a parsed query format string.
type FormatToken struct {
	// Kind is one of: literal, placeholder, sample.
	Kind FormatTokenKind
	// Text is the literal text (Kind=literal) or placeholder name (Kind=placeholder).
	Text string
	// Inner is the parsed body of a [%...] sample-repeated group. The inner
	// tokens may themselves contain literals and placeholders but not nested
	// sample groups.
	Inner []FormatToken
}

// FormatTokenKind tags FormatToken variants.
type FormatTokenKind int

const (
	// TokenLiteral is verbatim text.
	TokenLiteral FormatTokenKind = iota
	// TokenPlaceholder is a `%NAME` placeholder (e.g. CHROM, INFO/DP).
	TokenPlaceholder
	// TokenSample is a `[%...]` sample-repeated group.
	TokenSample
)

// ParseFormatString tokenises a bcftools-query format string. It returns a
// flat list of tokens; sample-repeated groups appear as a single TokenSample
// with their inner tokens populated.
func ParseFormatString(s string) ([]FormatToken, error) {
	if s == "" {
		return nil, fmt.Errorf("bcftools query: empty format string")
	}
	return tokenize(s, false)
}

// tokenize walks src; when inSample is true the `]` closing bracket terminates
// the current group. It returns an error if any structural problem is seen
// (unbalanced brackets, trailing backslash, ...).
func tokenize(src string, inSample bool) ([]FormatToken, error) {
	var out []FormatToken
	var lit strings.Builder
	flushLit := func() {
		if lit.Len() > 0 {
			out = append(out, FormatToken{Kind: TokenLiteral, Text: lit.String()})
			lit.Reset()
		}
	}
	i := 0
	for i < len(src) {
		c := src[i]
		switch c {
		case '\\':
			if i+1 >= len(src) {
				return nil, fmt.Errorf("bcftools query: trailing backslash in format string")
			}
			switch src[i+1] {
			case 'n':
				lit.WriteByte('\n')
			case 't':
				lit.WriteByte('\t')
			case '\\':
				lit.WriteByte('\\')
			case '"':
				lit.WriteByte('"')
			case '%':
				lit.WriteByte('%')
			case '[':
				lit.WriteByte('[')
			case ']':
				lit.WriteByte(']')
			default:
				// Unknown escape — keep the second character verbatim, matching
				// upstream's permissive behaviour.
				lit.WriteByte(src[i+1])
			}
			i += 2
		case '%':
			flushLit()
			name, consumed, err := readPlaceholder(src[i+1:])
			if err != nil {
				return nil, err
			}
			out = append(out, FormatToken{Kind: TokenPlaceholder, Text: name})
			i += 1 + consumed
		case '[':
			if inSample {
				return nil, fmt.Errorf("bcftools query: nested '[' in format string")
			}
			flushLit()
			end := strings.IndexByte(src[i+1:], ']')
			if end < 0 {
				return nil, fmt.Errorf("bcftools query: missing ']' in format string")
			}
			body := src[i+1 : i+1+end]
			if strings.IndexByte(body, '[') >= 0 {
				return nil, fmt.Errorf("bcftools query: nested '[' in format string")
			}
			// The inner body is parsed in outer (non-inSample) mode because
			// it is a separate substring; the caller has already stripped the
			// matching ']'.
			inner, err := tokenize(body, false)
			if err != nil {
				return nil, err
			}
			out = append(out, FormatToken{Kind: TokenSample, Inner: inner})
			i = i + 1 + end + 1
		case ']':
			if !inSample {
				return nil, fmt.Errorf("bcftools query: stray ']' in format string")
			}
			flushLit()
			return out, nil
		default:
			lit.WriteByte(c)
			i++
		}
	}
	if inSample {
		return nil, fmt.Errorf("bcftools query: unterminated '['")
	}
	flushLit()
	return out, nil
}

// readPlaceholder consumes one placeholder name starting at src. It accepts
// the optional INFO/<TAG> form. Returns the name, number of bytes consumed,
// and an error if the placeholder is malformed.
func readPlaceholder(src string) (string, int, error) {
	if len(src) == 0 {
		return "", 0, fmt.Errorf("bcftools query: trailing '%%' in format string")
	}
	// Allow %INFO/<TAG> — uppercase letters plus `/` for the INFO group, then
	// letters/digits/underscore.
	i := 0
	for i < len(src) {
		c := src[i]
		if c >= 'A' && c <= 'Z' {
			i++
			continue
		}
		if c >= 'a' && c <= 'z' {
			i++
			continue
		}
		if c >= '0' && c <= '9' {
			i++
			continue
		}
		if c == '_' || c == '/' {
			i++
			continue
		}
		break
	}
	if i == 0 {
		return "", 0, fmt.Errorf("bcftools query: empty placeholder")
	}
	return src[:i], i, nil
}

// writeHeaderRow emits the `-H` derived header line. The format string is
// scanned and each placeholder contributes its name (with the leading `#` only
// on the first column). Sample-repeated groups expand to one column per sample
// per placeholder.
func writeHeaderRow(w io.Writer, tokens []FormatToken, requested, headerSamples []string) error {
	samples := headerSamples
	if len(requested) > 0 {
		samples = filterSamplesByName(requested, headerSamples)
	}
	col := 1
	var sb strings.Builder
	for _, t := range tokens {
		switch t.Kind {
		case TokenLiteral:
			sb.WriteString(t.Text)
		case TokenPlaceholder:
			writeHeaderColumn(&sb, t.Text, &col)
		case TokenSample:
			for sIdx, name := range samples {
				_ = sIdx
				for _, inner := range t.Inner {
					switch inner.Kind {
					case TokenLiteral:
						sb.WriteString(inner.Text)
					case TokenPlaceholder:
						writeHeaderColumn(&sb, inner.Text+":"+name, &col)
					}
				}
			}
		}
	}
	s := sb.String()
	if !strings.HasSuffix(s, "\n") {
		s += "\n"
	}
	_, err := w.Write([]byte(s))
	return err
}

// writeHeaderColumn appends one column header to sb, prefixing the first
// column with the conventional "# " marker.
func writeHeaderColumn(sb *strings.Builder, name string, col *int) {
	if *col == 1 {
		sb.WriteString("# [1]" + name)
	} else {
		sb.WriteString("[" + strconv.Itoa(*col) + "]" + name)
	}
	*col++
}

// filterSamplesByName narrows headerSamples to the wanted set, preserving the
// requested order. Missing names are silently dropped.
func filterSamplesByName(wanted, headerSamples []string) []string {
	present := make(map[string]bool, len(headerSamples))
	for _, s := range headerSamples {
		present[s] = true
	}
	out := make([]string, 0, len(wanted))
	for _, w := range wanted {
		if present[w] {
			out = append(out, w)
		}
	}
	return out
}

// emitRecord formats one variant per the token list and writes it to w.
// headerSamples supplies the per-sample names used by the %SAMPLE token
// inside `[ ... ]` groups; an index into it is taken from sampleFilter (or
// the natural order when sampleFilter is nil).
func emitRecord(w io.Writer, tokens []FormatToken, v *vcf.Variant, sampleFilter []int, headerSamples []string) error {
	var sb strings.Builder
	for _, t := range tokens {
		switch t.Kind {
		case TokenLiteral:
			sb.WriteString(t.Text)
		case TokenPlaceholder:
			sb.WriteString(formatPlaceholder(t.Text, v, -1, headerSamples))
		case TokenSample:
			indexes := sampleFilter
			if indexes == nil {
				indexes = make([]int, len(v.Samples))
				for i := range v.Samples {
					indexes[i] = i
				}
			}
			// Upstream bcftools repeats the inner pattern verbatim for each
			// sample with no auto-inserted separator — any inter-sample
			// delimiter must be expressed in the inner pattern itself
			// (e.g. `[\t%GT]` puts a tab in front of each sample's GT).
			for _, idx := range indexes {
				for _, inner := range t.Inner {
					switch inner.Kind {
					case TokenLiteral:
						sb.WriteString(inner.Text)
					case TokenPlaceholder:
						sb.WriteString(formatPlaceholder(inner.Text, v, idx, headerSamples))
					}
				}
			}
		}
	}
	_, err := w.Write([]byte(sb.String()))
	return err
}

// formatPlaceholder resolves a single placeholder against v. sampleIdx is the
// sample index when evaluating inside `[ ... ]`, or -1 when in the outer
// scope (in which case sample-only placeholders fall back to a missing
// value). headerSamples carries the VCF sample names so the %SAMPLE token
// can resolve sampleIdx to a name.
func formatPlaceholder(name string, v *vcf.Variant, sampleIdx int, headerSamples []string) string {
	switch name {
	case "CHROM":
		return v.Chrom
	case "POS":
		return strconv.Itoa(v.Pos)
	case "POS0":
		// 0-based start, mirroring convert.c process_pos0 (line->pos).
		return strconv.Itoa(v.Pos - 1)
	case "END":
		// 1-based inclusive end = pos0 + rlen, where rlen follows INFO/END
		// when present (convert.c process_end, htslib line->rlen).
		return strconv.Itoa(queryVariantEnd(v))
	case "END0":
		// 0-based end (convert.c process_end0 = pos + rlen - 1).
		return strconv.Itoa(queryVariantEnd(v) - 1)
	case "REF":
		return v.Ref
	case "ALT":
		if len(v.Alt) == 0 {
			return "."
		}
		return strings.Join(v.Alt, ",")
	case "FIRST_ALT":
		// First ALT allele, or "." when the site has no ALT (convert.c
		// process_first_alt). Not publicly advertised but supported upstream.
		if len(v.Alt) == 0 {
			return "."
		}
		return v.Alt[0]
	case "QUAL":
		if v.Qual < 0 {
			return "."
		}
		// Match upstream's "%g"-ish default; trim trailing zeros for prettiness.
		return strconv.FormatFloat(v.Qual, 'g', -1, 64)
	case "ID":
		if v.ID == "" {
			return "."
		}
		return v.ID
	case "FILTER":
		if len(v.Filter) == 0 {
			return "."
		}
		return strings.Join(v.Filter, ";")
	case "TYPE":
		return variantType(v)
	case "IS_TS":
		// 1 for a transition (A<->G / C<->T) SNP/MNP, else 0 (convert.c
		// process_is_ts: |acgt2int(ref[0]) - acgt2int(alt0[0])| == 2).
		return variantIsTransition(v)
	case "INFO":
		// Bare %INFO (no /TAG) emits the entire INFO column, mirroring
		// upstream convert.c's process_info with a NULL key.
		return formatWholeInfo(v)
	case "SAMPLE":
		// %SAMPLE resolves to the per-sample name inside a [ ... ] group;
		// outside one it has no sample and prints ".".
		if sampleIdx < 0 || sampleIdx >= len(headerSamples) {
			return "."
		}
		return headerSamples[sampleIdx]
	case "GT":
		if sampleIdx < 0 {
			return "."
		}
		return sampleField(v, sampleIdx, "GT")
	case "TGT":
		if sampleIdx < 0 {
			return "."
		}
		return translatedGenotype(v, sampleIdx)
	case "TBCSQ":
		if sampleIdx < 0 {
			return "."
		}
		return expandTBCSQ(v, sampleIdx, "BCSQ")
	}
	if strings.HasPrefix(name, "INFO/") {
		key := name[len("INFO/"):]
		if val, ok := v.Info[key]; ok {
			if val == "" {
				return "1"
			}
			return val
		}
		return "."
	}
	if strings.HasPrefix(name, "FMT/") || strings.HasPrefix(name, "FORMAT/") {
		key := name[strings.IndexByte(name, '/')+1:]
		if sampleIdx < 0 {
			return "."
		}
		return sampleField(v, sampleIdx, key)
	}
	// Sample-context FORMAT field referenced bare (e.g. [%DP]).
	if sampleIdx >= 0 {
		return sampleField(v, sampleIdx, name)
	}
	return "."
}

// formatWholeInfo renders the entire INFO column, mirroring upstream
// convert.c's process_info when the format tag is a bare %INFO (no /TAG).
// Key order follows the record's original INFO order.
func formatWholeInfo(v *vcf.Variant) string {
	return v.InfoString()
}

// expandTBCSQ ports upstream convert.c's process_tbcsq: it decodes the
// per-sample FORMAT/BCSQ bitmask into the comma-separated list of the
// referenced INFO/BCSQ entries. The output is `hap1\thap2` when both
// haplotypes carry consequences (or "." per side when empty), matching
// `bcftools query -f'[%TBCSQ\n]'`.
//
// Bit layout (mirrors upstream's csq.c emission and convert.c expansion):
//   - each int32 value carries 30 effective bits (top 2 reserved);
//   - bit (2*k + ihap) within a value at array offset j corresponds to
//     INFO/BCSQ entry index (j*30 + 2*k + ihap) / 2 = j*15 + k.
func expandTBCSQ(v *vcf.Variant, sampleIdx int, tag string) string {
	if sampleIdx < 0 || sampleIdx >= len(v.Samples) {
		return "."
	}
	info, ok := v.Info[tag]
	if !ok || info == "" {
		return ".\t."
	}
	csqs := strings.Split(info, ",")
	raw := sampleField(v, sampleIdx, tag)
	if raw == "" || raw == "." {
		return ".\t."
	}
	parts := strings.Split(raw, ",")
	vals := make([]uint32, len(parts))
	for i, p := range parts {
		// Treat negative/missing as zero.
		if p == "" || p == "." {
			continue
		}
		n, err := strconv.ParseInt(p, 10, 64)
		if err != nil || n < 0 {
			continue
		}
		vals[i] = uint32(n)
	}
	var hap1, hap2 strings.Builder
	appendCSQ := func(sb *strings.Builder, idx int) {
		if idx < 0 || idx >= len(csqs) {
			return
		}
		if sb.Len() > 0 {
			sb.WriteByte(',')
		}
		sb.WriteString(csqs[idx])
	}
	const nbits = 30 // 30 effective bits per int32
	for j, val := range vals {
		if val == 0 {
			continue
		}
		// hap1 lives in even bits (0,2,4,...), hap2 in odd (1,3,5,...).
		for k := 0; k < nbits; k += 2 {
			idx := (j*nbits + k) / 2
			if val&(1<<uint(k)) != 0 {
				appendCSQ(&hap1, idx)
			}
			if val&(1<<uint(k+1)) != 0 {
				appendCSQ(&hap2, idx)
			}
		}
	}
	h1 := hap1.String()
	if h1 == "" {
		h1 = "."
	}
	h2 := hap2.String()
	if h2 == "" {
		h2 = "."
	}
	return h1 + "\t" + h2
}

// sampleField returns the FORMAT key value for a sample, or "." when missing.
func sampleField(v *vcf.Variant, idx int, key string) string {
	if idx < 0 || idx >= len(v.Samples) {
		return "."
	}
	val, ok := v.Samples[idx].Data[key]
	if !ok || val == "" {
		return "."
	}
	return val
}

// variantType returns SNP / MNP / INDEL / OTHER, following the bcftools
// definition: SNP when REF and every ALT are single bases; MNP when REF and
// every ALT have the same (greater than one) length; INDEL when any ALT
// length differs from REF; OTHER otherwise (no ALT, structural, ...).
func variantType(v *vcf.Variant) string {
	if len(v.Alt) == 0 {
		return "OTHER"
	}
	refLen := len(v.Ref)
	hasIndel := false
	allSNP := refLen == 1
	allSameLen := true
	for _, a := range v.Alt {
		if a == "" || a == "." || a == "*" {
			return "OTHER"
		}
		if strings.ContainsAny(a, "<>[]") {
			return "OTHER"
		}
		if len(a) != refLen {
			hasIndel = true
			allSameLen = false
		}
		if len(a) != 1 {
			allSNP = false
		}
	}
	if hasIndel {
		return "INDEL"
	}
	if allSNP {
		return "SNP"
	}
	if allSameLen && refLen > 1 {
		return "MNP"
	}
	return "OTHER"
}

// queryVariantEnd returns the 1-based inclusive end position for the %END
// query token, mirroring htslib's line->rlen accounting: the reference span
// is len(REF) unless an integer INFO/END widens it (symbolic/long records),
// in which case END is that value. Equivalent to convert.c's process_end
// (line->pos + line->rlen). It differs from the annotate-side variantEnd in
// honouring INFO/END.
func queryVariantEnd(v *vcf.Variant) int {
	rlen := len(v.Ref)
	if es, ok := v.Info["END"]; ok {
		if e, err := strconv.Atoi(strings.TrimSpace(es)); err == nil {
			rlen = e - (v.Pos - 1)
		}
	}
	return (v.Pos - 1) + rlen
}

// acgt2int maps a nucleotide byte to htslib's bcf_acgt2int encoding
// (A=0, C=1, G=2, T=3); any other byte returns -1.
func acgt2int(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	}
	return -1
}

// variantIsTransition returns "1" when the site is a SNP/MNP substitution
// whose first REF/ALT bases form a transition (A<->G or C<->T), else "0",
// mirroring convert.c process_is_ts.
func variantIsTransition(v *vcf.Variant) string {
	t := variantType(v)
	if t != "SNP" && t != "MNP" {
		return "0"
	}
	if len(v.Ref) == 0 || len(v.Alt) == 0 || len(v.Alt[0]) == 0 {
		return "0"
	}
	r := acgt2int(v.Ref[0])
	a := acgt2int(v.Alt[0][0])
	if r < 0 || a < 0 {
		return "0"
	}
	if r-a == 2 || a-r == 2 {
		return "1"
	}
	return "0"
}

// translatedGenotype maps the indices in a GT field to their allele strings.
// e.g. REF=A ALT=T,G and GT=1/2 → "T/G". Missing alleles stay as ".".
func translatedGenotype(v *vcf.Variant, idx int) string {
	if idx < 0 || idx >= len(v.Samples) {
		return "."
	}
	gt, ok := v.Samples[idx].Data["GT"]
	if !ok || gt == "" || gt == "." {
		return "."
	}
	alleles := append([]string{v.Ref}, v.Alt...)
	var sb strings.Builder
	tok := ""
	flush := func() {
		if tok == "" || tok == "." {
			sb.WriteString(".")
			return
		}
		n, err := strconv.Atoi(tok)
		if err != nil || n < 0 || n >= len(alleles) {
			sb.WriteString(tok)
			return
		}
		sb.WriteString(alleles[n])
	}
	for i := 0; i < len(gt); i++ {
		c := gt[i]
		switch c {
		case '/', '|':
			flush()
			sb.WriteByte(c)
			tok = ""
		default:
			tok += string(c)
		}
	}
	flush()
	return sb.String()
}
