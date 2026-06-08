// Package bcftools — `bcftools annotate` subcommand.
//
// `bcftools annotate` enriches a VCF/BCF stream with extra INFO/ID/FILTER
// fields drawn from a sidecar annotation source. The annotation source can
// be one of:
//
//   - A TAB-delimited table (`.tab`, `.tab.gz`, plain) whose first columns
//     identify the record (CHROM[, POS[, REF[, ALT]]]) and remaining columns
//     hold the annotation values. The `-c/--columns` flag maps each
//     annotation column to the destination field, e.g.
//     `CHROM,POS,REF,ALT,INFO/AC,INFO/AN`.
//   - A VCF/BCF (`.vcf`, `.vcf.gz`, `.bcf`) — a "transfer" annotation: each
//     INFO/FILTER tag named in `-c` is copied from a matching record in the
//     annotation file to the input record.
//
// Other supported flags:
//
//   - `-h/--header-lines FILE` injects extra `##...` lines into the output
//     header (in front of the existing header, after `##fileformat=`).
//   - `-x/--remove FIELD,...` drops the named fields from each record. Each
//     entry is `INFO/TAG`, `FILTER`, `ID`, or `INFO` (drop all INFO).
//   - `--rename-chrs FILE` rewrites CHROM via a two-column tab map.
//   - `-r/--regions LIST` is a post-filter (no index seek in v1).
//
// What is NOT covered in this v1 port (tracked in PARITY_ROADMAP):
//
//   - `-c CHROM,FROM,TO,...` BED-style ranges (we always require a single
//     POS column rather than FROM/TO; for full BED-range overlap rewrite,
//     use `bedtools annotate`).
//   - `-i/--include` and `-e/--exclude` expressions (the same engine works
//     here, but the v1 wiring focuses on the column mapping).
//   - `--mark-sites` and `--set-id` extras.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// AnnotateOptions controls the `bcftools annotate` subcommand.
type AnnotateOptions struct {
	// Annotations is the source path. Empty disables the column-transfer
	// step (useful when the goal is just `-x` field removal or `--rename-chrs`).
	Annotations string
	// Columns is the comma list described in the package docstring.
	Columns string
	// ColumnsFile is `-C/--columns-file FILE`: read column names (one per
	// row) from a file, equivalent to passing them via `-c`. An optional
	// second whitespace-separated token on each row selects a merge type
	// for `-l/--merge-logic` (currently rejected with rejection-parity).
	ColumnsFile string
	// HeaderLines is `-h/--header-lines FILE`: a file whose `##...` lines
	// are injected into the output header.
	HeaderLines string
	// HeaderLine is `-H/--header-line STR`: one or more literal `##...`
	// header lines appended to the output header. Repeatable.
	HeaderLine []string
	// Remove is the comma-list `-x` argument.
	Remove string
	// Regions is a post-filter on the input records.
	Regions []string
	// RegionsFile is the BED-like sidecar for `-R/--regions-file`.
	RegionsFile string
	// RegionsOverlap is upstream `--regions-overlap 0|1|2`. v1 accepts
	// the flag for parity but applies a post-filter based on simple
	// record-overlap (matching the default 1).
	RegionsOverlap int
	// RenameChromMap is the two-column tab file driving `--rename-chrs`.
	RenameChromMap string
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
	// SetID is the `-I/--set-id [+]FORMAT` value: a bcftools-query-style
	// format string used to populate the ID column. A leading '+' means
	// "only set when ID is missing (.)"; without the '+' the value is
	// unconditionally replaced. An empty string leaves the ID column
	// untouched. Mirrors vcfannotate.c:3250-3253.
	SetID string
	// IncludeExpr / ExcludeExpr are upstream `-i/--include` and
	// `-e/--exclude` filter expressions. Records are dropped (or, with
	// KeepSites, passed through unmodified) when ExcludeExpr is true or
	// IncludeExpr is false.
	IncludeExpr string
	ExcludeExpr string
	// KeepSites mirrors upstream `-k/--keep-sites`: instead of discarding
	// sites that fail -i/-e, leave them unchanged in the output.
	KeepSites bool
	// MarkSites is upstream `-m/--mark-sites [+-]TAG`: tag sites that
	// match (+) or do not match (-) the -a source with INFO/TAG.
	MarkSites string
	// Samples is upstream `-s/--samples [^]LIST`: comma-separated list
	// of samples to annotate (or exclude when prefixed with "^").
	Samples []string
	// SamplesFile is upstream `-S/--samples-file [^]FILE`: read sample
	// names from a file, "^"-prefix inverts.
	SamplesFile string
	// SamplesExclude is set when -s/-S used the `^` prefix; the named
	// samples are the exclusion set rather than the inclusion set.
	SamplesExclude bool
	// Force mirrors upstream `--force`: continue past malformed records
	// instead of failing.
	Force bool
	// NoVersion suppresses the `##bcftools_annotateVersion` /
	// `##bcftools_annotateCommand` provenance lines.
	NoVersion bool
	// PGCommand is the verbatim command line stamped into the
	// provenance header line when NoVersion is false.
	PGCommand string
	// WriteIndex is upstream `-W/--write-index[=FMT]`. Non-empty
	// values are "csi" or "tbi"; empty disables indexing.
	WriteIndex string
	// Verbosity mirrors upstream `--verbosity INT`.
	Verbosity int
}

// AnnotateFile is the file-aware entry point. It opens path through
// iohelper, reads the records, applies the requested annotations / removals,
// and writes the result to out.
func AnnotateFile(path string, out io.Writer, opts AnnotateOptions) (int, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools annotate: open %s: %w", path, err)
	}
	defer in.Close()
	return Annotate(in, out, opts)
}

// Annotate is the streaming entry point. It reads the input fully (matching
// upstream's behaviour where the annotation source dictates record-by-record
// random access), applies the transformations, and writes the result.
func Annotate(in io.Reader, out io.Writer, opts AnnotateOptions) (int, error) {
	hdr, recs, err := readAllVariants(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools annotate: %w", err)
	}

	// Apply --rename-chrs first so downstream key matching uses new names.
	var chromMap map[string]string
	if opts.RenameChromMap != "" {
		chromMap, err = loadChromRenameMap(opts.RenameChromMap)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: --rename-chrs %s: %w", opts.RenameChromMap, err)
		}
		for _, v := range recs {
			if n, ok := chromMap[v.Chrom]; ok {
				v.Chrom = n
			}
		}
		// Header contig lines: rewrite IDs in place.
		for i, m := range hdr.MetaInfo {
			if !strings.HasPrefix(m, "##contig=") {
				continue
			}
			k, id := structuredID(m)
			if k != "contig" {
				continue
			}
			if n, ok := chromMap[id]; ok {
				hdr.MetaInfo[i] = strings.Replace(m, "ID="+id, "ID="+n, 1)
			}
		}
	}

	// -h/--header-lines: append extra meta lines to the end of the
	// header (matches upstream vcfannotate.c::init_header_lines).
	if opts.HeaderLines != "" {
		extra, err := readHeaderLines(opts.HeaderLines)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: -h %s: %w", opts.HeaderLines, err)
		}
		hdr = appendMetaLines(hdr, extra)
	}
	// -H/--header-line STR (repeatable): each entry is a literal `##...`
	// line appended to the output header.
	if len(opts.HeaderLine) > 0 {
		extra := make([]string, 0, len(opts.HeaderLine))
		for _, line := range opts.HeaderLine {
			line = strings.TrimRight(line, "\r\n")
			if line == "" {
				continue
			}
			if !strings.HasPrefix(line, "##") {
				return 0, fmt.Errorf("bcftools annotate: -H %q: must begin with ##", line)
			}
			extra = append(extra, line)
		}
		hdr = appendMetaLines(hdr, extra)
	}

	// -s/-S samples restriction. Upstream's `annotate -s` restricts
	// the set of samples whose per-sample data is overwritten when -a
	// transfers FORMAT/* columns; it does NOT drop sample columns from
	// the output. We currently parse and accept the list (so it round-
	// trips alongside provenance) but do not project the per-sample
	// payload because the FORMAT-transfer path is not yet implemented
	// (tracked in PARITY_ROADMAP). The sample list will be wired into
	// the FORMAT-transfer path when that lands.
	_ = opts.Samples

	// -C/--columns-file FILE: read columns from a file (one per row).
	// Lines may contain an optional second whitespace-separated token,
	// which selects the merge-logic type — currently rejected.
	if opts.ColumnsFile != "" {
		extra, err := readColumnsFile(opts.ColumnsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: -C %s: %w", opts.ColumnsFile, err)
		}
		if opts.Columns != "" {
			opts.Columns = opts.Columns + "," + extra
		} else {
			opts.Columns = extra
		}
	}

	// -a + -c: column mapping.
	cols, err := parseAnnColumns(opts.Columns)
	if err != nil {
		return 0, fmt.Errorf("bcftools annotate: -c: %w", err)
	}
	matched := map[*vcf.Variant]bool{}
	if opts.Annotations != "" && len(cols) > 0 {
		switch {
		case strings.HasSuffix(opts.Annotations, ".vcf") ||
			strings.HasSuffix(opts.Annotations, ".vcf.gz") ||
			strings.HasSuffix(opts.Annotations, ".bcf"):
			if err := applyVCFAnnotations(opts.Annotations, recs, cols, matched); err != nil {
				return 0, fmt.Errorf("bcftools annotate: %w", err)
			}
		default:
			if err := applyTableAnnotations(opts.Annotations, recs, cols, hdr, matched); err != nil {
				return 0, fmt.Errorf("bcftools annotate: %w", err)
			}
		}
	}
	// -m/--mark-sites [+-]TAG: tag sites that match (+) or do not match
	// (-) the -a source with INFO/TAG=1.
	if opts.MarkSites != "" {
		if err := applyMarkSites(recs, hdr, opts.MarkSites, matched, opts.Annotations != ""); err != nil {
			return 0, fmt.Errorf("bcftools annotate: --mark-sites: %w", err)
		}
	}

	// -x: field removal.
	if opts.Remove != "" {
		applyRemovals(recs, hdr, opts.Remove)
	}

	// -I/--set-id: rewrite the ID column from a query-style format string.
	if opts.SetID != "" {
		if err := applySetID(recs, opts.SetID); err != nil {
			return 0, fmt.Errorf("bcftools annotate: --set-id: %w", err)
		}
	}

	// -i/--include and -e/--exclude expressions: filter records (or pass
	// through when -k/--keep-sites is set).
	var (
		incF, excF *Filter
		failed     map[*vcf.Variant]bool
	)
	if opts.IncludeExpr != "" {
		f, err := CompileFilterWithHeader(opts.IncludeExpr, hdr)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: --include: %w", err)
		}
		incF = f
	}
	if opts.ExcludeExpr != "" {
		f, err := CompileFilterWithHeader(opts.ExcludeExpr, hdr)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: --exclude: %w", err)
		}
		excF = f
	}
	if incF != nil || excF != nil {
		failed = map[*vcf.Variant]bool{}
		for _, v := range recs {
			ok := true
			if incF != nil && !incF.Eval(v) {
				ok = false
			}
			if ok && excF != nil && excF.Eval(v) {
				ok = false
			}
			if !ok {
				failed[v] = true
			}
		}
	}

	// --no-version: stamp provenance unless suppressed.
	if !opts.NoVersion {
		stampAnnotateProvenance(hdr, opts.PGCommand)
	}

	// Region post-filter. Merge -r and -R into one set.
	regionsSpec := opts.Regions
	if opts.RegionsFile != "" {
		regs, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: -R %s: %w", opts.RegionsFile, err)
		}
		regionsSpec = append(regionsSpec, regs...)
	}
	regions, err := parseRegions(regionsSpec)
	if err != nil {
		return 0, err
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}
	count := 0
	for _, v := range recs {
		if len(regions) > 0 && !overlapsAny(v, regions) {
			continue
		}
		if failed[v] && !opts.KeepSites {
			continue
		}
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// readColumnsFile reads `-C/--columns-file` entries: one column name per
// non-blank, non-comment line. A second whitespace-separated token (the
// merge-logic type) is rejected with rejection-parity to flag the
// unimplemented behaviour up front.
func readColumnsFile(path string) (string, error) {
	f, err := os.Open(path)
	if err != nil {
		return "", err
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
		fields := strings.Fields(line)
		if len(fields) > 1 {
			return "", fmt.Errorf("--merge-logic is to be implemented, please open an issue on github")
		}
		names = append(names, fields[0])
	}
	if err := sc.Err(); err != nil {
		return "", err
	}
	return strings.Join(names, ","), nil
}

// stampAnnotateProvenance appends the upstream-style
// `##bcftools_annotateVersion` / `##bcftools_annotateCommand` lines to
// hdr. Mirrors appendNormProvenance.
func stampAnnotateProvenance(hdr *vcf.Header, cmdLine string) {
	if hdr == nil {
		return
	}
	hdr.MetaInfo = append(hdr.MetaInfo,
		`##bcftools_annotateVersion=bio_ai_experiment`,
	)
	if cmdLine == "" {
		cmdLine = "annotate"
	}
	hdr.MetaInfo = append(hdr.MetaInfo,
		`##bcftools_annotateCommand=`+cmdLine,
	)
}

// annColumn is one entry in the parsed -c list.
type annColumn struct {
	// Kind: "CHROM", "POS", "REF", "ALT", "ID", "FILTER", "INFO" or "-"
	// (the "skip" sentinel).
	Kind string
	// Tag is the INFO tag when Kind == "INFO".
	Tag string
}

// parseAnnColumns parses a `-c CHROM,POS,REF,ALT,INFO/AC,INFO/AN,-` style
// argument. The "-" entry means "skip this column".
func parseAnnColumns(spec string) ([]annColumn, error) {
	if spec == "" {
		return nil, nil
	}
	var out []annColumn
	for _, p := range strings.Split(spec, ",") {
		p = strings.TrimSpace(p)
		switch p {
		case "":
			return nil, fmt.Errorf("empty entry in -c %q", spec)
		case "-":
			out = append(out, annColumn{Kind: "-"})
		case "CHROM", "POS", "REF", "ALT", "ID", "FILTER":
			out = append(out, annColumn{Kind: p})
		default:
			if strings.HasPrefix(p, "INFO/") {
				out = append(out, annColumn{Kind: "INFO", Tag: p[len("INFO/"):]})
				continue
			}
			// Bare tag — assume INFO.
			out = append(out, annColumn{Kind: "INFO", Tag: p})
		}
	}
	return out, nil
}

// applyTableAnnotations reads a TAB-delimited annotation file and copies the
// chosen columns onto each matching record. Matching uses the (CHROM, POS)
// key by default; (CHROM, POS, REF) and (CHROM, POS, REF, ALT) tighten the
// match when those columns are present in the column spec.
func applyTableAnnotations(path string, recs []*vcf.Variant, cols []annColumn, hdr *vcf.Header, matched map[*vcf.Variant]bool) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return err
	}
	defer in.Close()

	// Build a lookup from key -> map[tag]value.
	// Determine which columns of the table form the key.
	keyCols := []int{}
	for i, c := range cols {
		switch c.Kind {
		case "CHROM", "POS", "REF", "ALT":
			keyCols = append(keyCols, i)
		}
	}
	if len(keyCols) < 2 {
		return fmt.Errorf("annotation table needs at least CHROM and POS columns")
	}

	// Walk the table, accumulating INFO assignments.
	type assignment struct {
		tag, value string
	}
	type tableRow struct {
		assigns []assignment
		id      string
		filter  string
	}
	keyed := make(map[string]tableRow)
	br := bufio.NewReader(in)
	for {
		line, err := br.ReadString('\n')
		if line == "" && err != nil {
			break
		}
		line = strings.TrimRight(line, "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			if err != nil {
				break
			}
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < len(cols) {
			if err != nil {
				break
			}
			continue
		}
		key := buildTableKey(fields, cols)
		row := tableRow{}
		for i, c := range cols {
			switch c.Kind {
			case "ID":
				row.id = fields[i]
			case "FILTER":
				row.filter = fields[i]
			case "INFO":
				if fields[i] != "." && fields[i] != "" {
					row.assigns = append(row.assigns, assignment{tag: c.Tag, value: fields[i]})
				}
			}
		}
		keyed[key] = row
		if err != nil {
			break
		}
	}

	// Apply.
	for _, v := range recs {
		key := buildVariantKey(v, cols)
		row, ok := keyed[key]
		if !ok {
			continue
		}
		if matched != nil {
			matched[v] = true
		}
		if row.id != "" && row.id != "." {
			v.ID = row.id
		}
		if row.filter != "" && row.filter != "." {
			v.Filter = strings.Split(row.filter, ";")
		}
		for _, a := range row.assigns {
			if v.Info == nil {
				v.Info = map[string]string{}
			}
			if _, exists := v.Info[a.tag]; !exists {
				v.InfoOrder = append(v.InfoOrder, a.tag)
			}
			v.Info[a.tag] = a.value
		}
	}
	// Inject ##INFO header lines for tags we wrote that are not yet declared.
	for _, c := range cols {
		if c.Kind != "INFO" {
			continue
		}
		if !hasInfoLine(hdr, c.Tag) {
			hdr.MetaInfo = append([]string{hdr.MetaInfo[0],
				fmt.Sprintf(`##INFO=<ID=%s,Number=.,Type=String,Description="Added by bcftools annotate">`, c.Tag)},
				hdr.MetaInfo[1:]...)
		}
	}
	return nil
}

// applyVCFAnnotations transfers the named INFO/FILTER fields from the
// matching records of a VCF/BCF annotation file to the input records.
func applyVCFAnnotations(path string, recs []*vcf.Variant, cols []annColumn, matched map[*vcf.Variant]bool) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return err
	}
	defer in.Close()
	_, annRecs, err := readAllVariants(in)
	if err != nil {
		return err
	}
	// Index by (CHROM, POS, REF, ALT).
	idx := map[string]*vcf.Variant{}
	for _, v := range annRecs {
		idx[strictKey(v)] = v
	}
	wantInfo := map[string]bool{}
	wantFilter := false
	wantID := false
	for _, c := range cols {
		switch c.Kind {
		case "INFO":
			if c.Tag != "" {
				wantInfo[c.Tag] = true
			}
		case "FILTER":
			wantFilter = true
		case "ID":
			wantID = true
		}
	}
	for _, v := range recs {
		src, ok := idx[strictKey(v)]
		if !ok {
			continue
		}
		if matched != nil {
			matched[v] = true
		}
		for tag := range wantInfo {
			val, has := src.Info[tag]
			if !has {
				continue
			}
			if v.Info == nil {
				v.Info = map[string]string{}
			}
			if _, exists := v.Info[tag]; !exists {
				v.InfoOrder = append(v.InfoOrder, tag)
			}
			v.Info[tag] = val
		}
		if wantFilter && len(src.Filter) > 0 {
			v.Filter = append([]string{}, src.Filter...)
		}
		if wantID && src.ID != "" && src.ID != "." {
			v.ID = src.ID
		}
	}
	return nil
}

// applyMarkSites implements `-m/--mark-sites [+-]TAG`: when the leading
// sign is '+' (or absent) the TAG flag is set on records that matched
// the -a source; when the sign is '-' it is set on records that did NOT
// match. Mirrors vcfannotate.c:mark_sites.
func applyMarkSites(recs []*vcf.Variant, hdr *vcf.Header, spec string, matched map[*vcf.Variant]bool, haveAnnSource bool) error {
	mark := true
	tag := spec
	switch spec[0] {
	case '+':
		tag = spec[1:]
		mark = true
	case '-':
		tag = spec[1:]
		mark = false
	}
	if tag == "" {
		return fmt.Errorf("--mark-sites: missing TAG name in %q", spec)
	}
	if !haveAnnSource {
		// Without -a we have no notion of "matched"; treat all sites as
		// unmatched (mark="-" => all sites tagged; mark="+" => no
		// sites tagged), matching upstream's interpretation.
	}
	if !hasInfoLine(hdr, tag) {
		hdr.MetaInfo = append([]string{hdr.MetaInfo[0],
			fmt.Sprintf(`##INFO=<ID=%s,Number=0,Type=Flag,Description="Site %s -a source (added by bcftools annotate --mark-sites)">`,
				tag,
				map[bool]string{true: "matched", false: "did not match"}[mark])},
			hdr.MetaInfo[1:]...)
	}
	for _, v := range recs {
		isMatch := matched[v]
		if (mark && isMatch) || (!mark && !isMatch) {
			if v.Info == nil {
				v.Info = map[string]string{}
			}
			if _, exists := v.Info[tag]; !exists {
				v.InfoOrder = append(v.InfoOrder, tag)
			}
			v.Info[tag] = "" // Flag (Number=0): no value, just presence.
		}
	}
	return nil
}

// applyRemovals processes -x/--remove. Each entry is one of `INFO/TAG`,
// `INFO` (drop all INFO), `FILTER`, `FILTER/NAME`, or `ID`.
func applyRemovals(recs []*vcf.Variant, hdr *vcf.Header, spec string) {
	for _, raw := range strings.Split(spec, ",") {
		ent := strings.TrimSpace(raw)
		switch {
		case ent == "INFO":
			for _, v := range recs {
				v.Info = map[string]string{}
				v.InfoOrder = nil
			}
		case strings.HasPrefix(ent, "INFO/"):
			tag := ent[len("INFO/"):]
			for _, v := range recs {
				if v.Info != nil {
					delete(v.Info, tag)
				}
				if len(v.InfoOrder) > 0 {
					out := v.InfoOrder[:0]
					for _, k := range v.InfoOrder {
						if k != tag {
							out = append(out, k)
						}
					}
					v.InfoOrder = out
				}
			}
			// Strip the matching ##INFO header line.
			out := hdr.MetaInfo[:0]
			for _, m := range hdr.MetaInfo {
				k, id := structuredID(m)
				if k == "INFO" && id == tag {
					continue
				}
				out = append(out, m)
			}
			hdr.MetaInfo = out
		case ent == "FILTER":
			for _, v := range recs {
				v.Filter = []string{"."}
			}
		case strings.HasPrefix(ent, "FILTER/"):
			name := ent[len("FILTER/"):]
			for _, v := range recs {
				out := v.Filter[:0]
				for _, f := range v.Filter {
					if f == name {
						continue
					}
					out = append(out, f)
				}
				if len(out) == 0 {
					out = []string{"."}
				}
				v.Filter = out
			}
		case ent == "ID":
			for _, v := range recs {
				v.ID = "."
			}
		case ent == "FORMAT" || ent == "FMT":
			// Drop everything except GT (matches upstream: removing all
			// FORMAT keeps GT as the genotype column is part of the
			// record contract).
			for _, v := range recs {
				v.Format = filterFormatKeep(v.Format, "GT")
				for i := range v.Samples {
					if v.Samples[i].Data != nil {
						for k := range v.Samples[i].Data {
							if k != "GT" {
								delete(v.Samples[i].Data, k)
							}
						}
					}
				}
			}
			hdr.MetaInfo = filterFormatHeaderKeep(hdr.MetaInfo, "GT")
		case strings.HasPrefix(ent, "FORMAT/") || strings.HasPrefix(ent, "FMT/"):
			tag := ent[strings.IndexByte(ent, '/')+1:]
			for _, v := range recs {
				v.Format = filterFormatDrop(v.Format, tag)
				for i := range v.Samples {
					if v.Samples[i].Data != nil {
						delete(v.Samples[i].Data, tag)
					}
				}
			}
			// Strip the matching ##FORMAT header line.
			out := hdr.MetaInfo[:0]
			for _, m := range hdr.MetaInfo {
				k, id := structuredID(m)
				if k == "FORMAT" && id == tag {
					continue
				}
				out = append(out, m)
			}
			hdr.MetaInfo = out
		}
	}
}

// filterFormatDrop returns a copy of fmts with the named tag removed
// (preserving order). The empty input case returns nil.
func filterFormatDrop(fmts []string, tag string) []string {
	out := fmts[:0]
	for _, k := range fmts {
		if k == tag {
			continue
		}
		out = append(out, k)
	}
	return out
}

// filterFormatKeep returns a copy of fmts with everything except keep
// removed.
func filterFormatKeep(fmts []string, keep string) []string {
	out := fmts[:0]
	for _, k := range fmts {
		if k == keep {
			out = append(out, k)
		}
	}
	return out
}

// filterFormatHeaderKeep drops every ##FORMAT=<ID=...> meta line whose
// ID is not in keep. Other meta lines are passed through unchanged.
func filterFormatHeaderKeep(meta []string, keep string) []string {
	out := meta[:0]
	for _, m := range meta {
		k, id := structuredID(m)
		if k == "FORMAT" && id != keep {
			continue
		}
		out = append(out, m)
	}
	return out
}

// readHeaderLines reads a file of `##...` meta lines, ignoring blank lines
// and `#`-comments that don't start with `##`.
func readHeaderLines(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var out []string
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := sc.Text()
		if !strings.HasPrefix(line, "##") {
			continue
		}
		out = append(out, line)
	}
	return out, sc.Err()
}

// loadChromRenameMap reads a two-column tab file `OLD\tNEW` (one per line)
// and returns a string->string map.
func loadChromRenameMap(path string) (map[string]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	out := map[string]string{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.SplitN(line, "\t", 2)
		if len(fields) != 2 {
			return nil, fmt.Errorf("bad rename-chrs line %q (need OLD<TAB>NEW)", line)
		}
		out[strings.TrimSpace(fields[0])] = strings.TrimSpace(fields[1])
	}
	return out, sc.Err()
}

// appendMetaLines appends extra ##... lines to the end of hdr,
// de-duplicating by exact-string equality. Matches upstream's behaviour
// of placing -h/-H lines after the existing meta block.
func appendMetaLines(hdr *vcf.Header, extra []string) *vcf.Header {
	seen := map[string]bool{}
	for _, m := range hdr.MetaInfo {
		seen[m] = true
	}
	for _, m := range extra {
		if seen[m] {
			continue
		}
		hdr.MetaInfo = append(hdr.MetaInfo, m)
		seen[m] = true
	}
	return hdr
}

// injectMetaLines inserts extra ##... lines into hdr after the
// ##fileformat line (or at the top if none is present), de-duplicating by
// exact-string equality.
func injectMetaLines(hdr *vcf.Header, extra []string) *vcf.Header {
	out := &vcf.Header{Samples: hdr.Samples}
	insertAt := 0
	for i, m := range hdr.MetaInfo {
		out.MetaInfo = append(out.MetaInfo, m)
		if strings.HasPrefix(m, "##fileformat=") {
			insertAt = i + 1
		}
	}
	seen := map[string]bool{}
	for _, m := range out.MetaInfo {
		seen[m] = true
	}
	added := []string{}
	for _, m := range extra {
		if seen[m] {
			continue
		}
		added = append(added, m)
		seen[m] = true
	}
	out.MetaInfo = append(out.MetaInfo[:insertAt], append(added, out.MetaInfo[insertAt:]...)...)
	return out
}

// hasInfoLine reports whether the header already declares an ##INFO line
// for the named tag.
func hasInfoLine(hdr *vcf.Header, tag string) bool {
	for _, m := range hdr.MetaInfo {
		k, id := structuredID(m)
		if k == "INFO" && id == tag {
			return true
		}
	}
	return false
}

// strictKey is the canonical key for VCF-source annotation matching.
func strictKey(v *vcf.Variant) string {
	return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
}

// buildTableKey assembles the lookup key for a TSV row. The set of key
// columns matches the variant-side buildVariantKey.
func buildTableKey(fields []string, cols []annColumn) string {
	parts := make([]string, 0, 4)
	for i, c := range cols {
		switch c.Kind {
		case "CHROM":
			parts = append(parts, fields[i])
		case "POS":
			parts = append(parts, fields[i])
		case "REF":
			parts = append(parts, fields[i])
		case "ALT":
			parts = append(parts, fields[i])
		}
	}
	return strings.Join(parts, "\t")
}

// buildVariantKey assembles the same key from a Variant.
func buildVariantKey(v *vcf.Variant, cols []annColumn) string {
	parts := make([]string, 0, 4)
	for _, c := range cols {
		switch c.Kind {
		case "CHROM":
			parts = append(parts, v.Chrom)
		case "POS":
			parts = append(parts, fmt.Sprintf("%d", v.Pos))
		case "REF":
			parts = append(parts, v.Ref)
		case "ALT":
			parts = append(parts, strings.Join(v.Alt, ","))
		}
	}
	return strings.Join(parts, "\t")
}
