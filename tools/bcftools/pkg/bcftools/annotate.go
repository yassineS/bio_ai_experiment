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
	"strconv"
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
	// HeaderLines is a path whose `##...` lines are injected into the
	// output header.
	HeaderLines string
	// Remove is the comma-list `-x` argument.
	Remove string
	// Regions is a post-filter on the input records.
	Regions []string
	// RegionsFile is the BED-like sidecar.
	RegionsFile string
	// RenameChromMap is the two-column tab file driving `--rename-chrs`.
	RenameChromMap string
	// SetID is the `--set-id [+]<FORMAT>` macro string. Empty disables it.
	SetID string
	// MergeLogic is the raw `--merge-logic` argument (TAG:LOGIC[,...]).
	MergeLogic string
	// MinOverlap is the raw `--min-overlap ANN:VCF` argument.
	MinOverlap string
	// PairLogic is the raw `--pair-logic` argument (default "some").
	PairLogic string
	// SingleOverlaps enables `--single-overlaps` (apply only the first
	// overlapping annotation row).
	SingleOverlaps bool
	// RenameAnnots is the `--rename-annots` map file path.
	RenameAnnots string
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
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

	// -h/--header-lines: prepend extra meta lines (after fileformat).
	if opts.HeaderLines != "" {
		extra, err := readHeaderLines(opts.HeaderLines)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: -h %s: %w", opts.HeaderLines, err)
		}
		hdr = injectMetaLines(hdr, extra)
	}

	// --rename-annots: rename existing INFO/FORMAT/FILTER tags. Upstream
	// applies this before the source annotations are transferred.
	if opts.RenameAnnots != "" {
		maps, err := loadRenameAnnots(opts.RenameAnnots)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: --rename-annots %s: %w", opts.RenameAnnots, err)
		}
		applyRenameAnnots(recs, hdr, maps)
	}

	// Parse the advanced overlap/merge/pair options.
	mergeLogic, err := ParseMergeLogicSpec(opts.MergeLogic)
	if err != nil {
		return 0, fmt.Errorf("bcftools annotate: %w", err)
	}
	if opts.SingleOverlaps && opts.MergeLogic != "" {
		return 0, fmt.Errorf("bcftools annotate: the options --merge-logic and --single-overlaps cannot be combined")
	}
	minOverlap, err := ParseMinOverlap(opts.MinOverlap)
	if err != nil {
		return 0, fmt.Errorf("bcftools annotate: %w", err)
	}
	if opts.MinOverlap != "" && opts.SingleOverlaps {
		return 0, fmt.Errorf("bcftools annotate: the options --single-overlaps and --min-overlap cannot be combined")
	}
	pairLogic := PairSome
	if opts.PairLogic != "" {
		pairLogic, err = ParsePairLogic(opts.PairLogic)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: %w", err)
		}
	}

	// -a + -c: column mapping.
	cols, err := parseAnnColumns(opts.Columns)
	if err != nil {
		return 0, fmt.Errorf("bcftools annotate: -c: %w", err)
	}
	isVCFSource := strings.HasSuffix(opts.Annotations, ".vcf") ||
		strings.HasSuffix(opts.Annotations, ".vcf.gz") ||
		strings.HasSuffix(opts.Annotations, ".bcf")
	if opts.MinOverlap != "" && isVCFSource {
		return 0, fmt.Errorf("bcftools annotate: the --min-overlap option cannot be used when annotating from a VCF")
	}
	if opts.MergeLogic != "" && isVCFSource {
		return 0, fmt.Errorf("bcftools annotate: the --merge-logic is intended for use with BED or TAB-delimited files only.")
	}
	if opts.Annotations != "" && len(cols) > 0 {
		switch {
		case isVCFSource:
			if err := applyVCFAnnotations(opts.Annotations, recs, cols, pairLogic, hdr); err != nil {
				return 0, fmt.Errorf("bcftools annotate: %w", err)
			}
		case hasRangeColumns(cols):
			if err := applyTableRangeAnnotations(opts.Annotations, recs, cols, hdr, mergeLogic, minOverlap, opts.SingleOverlaps); err != nil {
				return 0, fmt.Errorf("bcftools annotate: %w", err)
			}
		default:
			if err := applyTableAnnotations(opts.Annotations, recs, cols, hdr); err != nil {
				return 0, fmt.Errorf("bcftools annotate: %w", err)
			}
		}
	}

	// -x: field removal.
	if opts.Remove != "" {
		applyRemovals(recs, hdr, opts.Remove)
	}

	// --set-id: macro-expand the ID column. Upstream applies this last, after
	// every other transformation.
	if opts.SetID != "" {
		prog, err := ParseSetID(opts.SetID)
		if err != nil {
			return 0, fmt.Errorf("bcftools annotate: --set-id: %w", err)
		}
		for _, v := range recs {
			val, ok := prog.expand(v)
			if !ok {
				continue
			}
			if prog.onlyIfEmpty && v.ID != "" && v.ID != "." {
				continue
			}
			v.ID = val
		}
	}

	// Region post-filter.
	regions, err := parseRegions(opts.Regions)
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
		if err := w.Write(v); err != nil {
			return count, err
		}
		count++
	}
	return count, w.Flush()
}

// annColumn is one entry in the parsed -c list.
type annColumn struct {
	// Kind: "CHROM", "POS", "BEG", "END", "REF", "ALT", "ID", "FILTER",
	// "INFO" or "-" (the "skip" sentinel). BEG/END (aliases FROM/TO) mark
	// the interval columns of a range-style annotation table.
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
		case "BEG", "FROM":
			out = append(out, annColumn{Kind: "BEG"})
		case "END", "TO":
			out = append(out, annColumn{Kind: "END"})
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
func applyTableAnnotations(path string, recs []*vcf.Variant, cols []annColumn, hdr *vcf.Header) error {
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

// hasRangeColumns reports whether the -c spec includes BEG/END (FROM/TO)
// interval columns, selecting the range-overlap annotation path instead of
// the point-key path.
func hasRangeColumns(cols []annColumn) bool {
	beg, end := false, false
	for _, c := range cols {
		switch c.Kind {
		case "BEG":
			beg = true
		case "END":
			end = true
		}
	}
	return beg && end
}

// annRangeRow is one parsed interval row of a range-style annotation table.
type annRangeRow struct {
	chrom    string
	beg, end int
	id       string
	filter   string
	assigns  []struct{ tag, value string }
}

// applyTableRangeAnnotations reads a TAB-delimited annotation table whose
// columns include CHROM,BEG,END (range form) and copies the chosen columns
// onto every overlapping VCF record. Multi-overlap merging is governed by
// mergeLogic (default first), the per-row overlap thresholds by minOverlap,
// and singleOverlaps restricts to the first overlapping row. It mirrors the
// observable behaviour of annotate_from_regidx.
func applyTableRangeAnnotations(path string, recs []*vcf.Variant, cols []annColumn, hdr *vcf.Header, mergeLogic map[string]MergeLogic, minOverlap MinOverlap, singleOverlaps bool) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return err
	}
	defer in.Close()

	chromIdx, begIdx, endIdx := -1, -1, -1
	for i, c := range cols {
		switch c.Kind {
		case "CHROM":
			chromIdx = i
		case "BEG":
			begIdx = i
		case "END":
			endIdx = i
		}
	}
	if chromIdx < 0 || begIdx < 0 || endIdx < 0 {
		return fmt.Errorf("annotation table needs CHROM, BEG and END columns")
	}

	// Read all interval rows, grouped by chromosome.
	byChrom := map[string][]annRangeRow{}
	br := bufio.NewReader(in)
	for {
		line, rerr := br.ReadString('\n')
		if line == "" && rerr != nil {
			break
		}
		line = strings.TrimRight(line, "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			if rerr != nil {
				break
			}
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < len(cols) {
			if rerr != nil {
				break
			}
			continue
		}
		beg, err1 := strconv.Atoi(strings.TrimSpace(fields[begIdx]))
		end, err2 := strconv.Atoi(strings.TrimSpace(fields[endIdx]))
		if err1 != nil || err2 != nil {
			if rerr != nil {
				break
			}
			continue
		}
		row := annRangeRow{chrom: fields[chromIdx], beg: beg, end: end}
		for i, c := range cols {
			switch c.Kind {
			case "ID":
				row.id = fields[i]
			case "FILTER":
				row.filter = fields[i]
			case "INFO":
				row.assigns = append(row.assigns, struct{ tag, value string }{c.Tag, fields[i]})
			}
		}
		byChrom[row.chrom] = append(byChrom[row.chrom], row)
		if rerr != nil {
			break
		}
	}

	for _, v := range recs {
		vbeg := v.Pos
		vend := variantEnd(v)
		// Collect overlapping rows in table order.
		var overlaps []annRangeRow
		for _, row := range byChrom[v.Chrom] {
			if vbeg > row.end || vend < row.beg {
				continue
			}
			if !passesMinOverlap(minOverlap, row.beg, row.end, vbeg, vend) {
				continue
			}
			overlaps = append(overlaps, row)
			if singleOverlaps {
				break
			}
		}
		if len(overlaps) == 0 {
			continue
		}

		// ID and FILTER follow first-overlap semantics.
		first := overlaps[0]
		if first.id != "" && first.id != "." {
			v.ID = first.id
		}
		if first.filter != "" && first.filter != "." {
			v.Filter = strings.Split(first.filter, ";")
		}

		// Merge INFO assignments per the chosen logic.
		accs := map[string]*mergeAccumulator{}
		var tagOrder []string
		for _, row := range overlaps {
			for _, a := range row.assigns {
				acc, ok := accs[a.tag]
				if !ok {
					logic := mergeLogic[a.tag]
					acc = newMergeAccumulator(logic, infoTagIsInteger(hdr, a.tag))
					accs[a.tag] = acc
					tagOrder = append(tagOrder, a.tag)
				}
				acc.add(a.value)
			}
		}
		for _, tag := range tagOrder {
			val, ok := accs[tag].reduce()
			if !ok {
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
// matching records of a VCF/BCF annotation file to the input records. The
// pairLogic mode controls how an annotation record is paired with an input
// record (see --pair-logic).
func applyVCFAnnotations(path string, recs []*vcf.Variant, cols []annColumn, pairLogic PairLogic, hdr *vcf.Header) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return err
	}
	defer in.Close()
	annHdr, annRecs, err := readAllVariants(in)
	if err != nil {
		return err
	}
	// Group annotation records by (CHROM, POS) so pair-logic can scan the
	// candidates at a site.
	byPos := map[string][]*vcf.Variant{}
	for _, v := range annRecs {
		k := v.Chrom + "\t" + strconv.Itoa(v.Pos)
		byPos[k] = append(byPos[k], v)
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
		var src *vcf.Variant
		for _, cand := range byPos[v.Chrom+"\t"+strconv.Itoa(v.Pos)] {
			if pairLogicMatches(v, cand, pairLogic) {
				src = cand
				break
			}
		}
		if src == nil {
			continue
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
	// Carry the source's ##INFO header definitions for each transferred tag
	// into the output header (matching upstream, which copies the matching
	// header records from the annotation file). We iterate cols (not the
	// wantInfo map) so the appended lines follow the -c column order.
	for _, c := range cols {
		if c.Kind != "INFO" || c.Tag == "" {
			continue
		}
		if hasInfoLine(hdr, c.Tag) {
			continue
		}
		def := findMetaLineByID(annHdr, "INFO", c.Tag)
		if def == "" {
			def = fmt.Sprintf(`##INFO=<ID=%s,Number=.,Type=String,Description="Added by bcftools annotate">`, c.Tag)
		}
		hdr.MetaInfo = append(hdr.MetaInfo, def)
	}
	return nil
}

// findMetaLineByID returns the verbatim ##<kind>=<...ID=id...> meta line from
// hdr, or "" when no such line exists.
func findMetaLineByID(hdr *vcf.Header, kind, id string) string {
	for _, m := range hdr.MetaInfo {
		k, mid := structuredID(m)
		if k == kind && mid == id {
			return m
		}
	}
	return ""
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
		}
	}
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

// injectMetaLines appends extra ##... lines to the end of hdr's meta block
// (just before the #CHROM line), de-duplicating by exact-string equality.
// This mirrors upstream's bcf_hdr_append, which adds new header records after
// the existing ones.
func injectMetaLines(hdr *vcf.Header, extra []string) *vcf.Header {
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	seen := map[string]bool{}
	for _, m := range out.MetaInfo {
		seen[m] = true
	}
	for _, m := range extra {
		if seen[m] {
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
		seen[m] = true
	}
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

// infoTagIsInteger reports whether the ##INFO line for tag declares
// Type=Integer. It is used so the numeric merge reductions truncate integer
// tags exactly as upstream does.
func infoTagIsInteger(hdr *vcf.Header, tag string) bool {
	for _, m := range hdr.MetaInfo {
		k, id := structuredID(m)
		if k != "INFO" || id != tag {
			continue
		}
		return strings.Contains(m, "Type=Integer")
	}
	return false
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
