// Native port of the upstream `split-vep` plugin (plugins/split-vep.c). It
// queries the structured INFO/CSQ (VEP) or INFO/BCSQ (bcftools/csq) annotation,
// splitting it into its pipe-delimited subfields. Two output modes are
// supported, exactly as upstream:
//
//   - annotate mode (default): -c/--columns extracts the named/indexed subfields
//     into new INFO tags on each record (re-emitting VCF/BCF).
//   - text mode: -f/--format prints a `bcftools query`-style line per record (or
//     per transcript with -d), dropping records without a passing consequence.
//
// Transcript and consequence selection (-s/--select) supports the all/worst
// transcript modes and the any / term[+|-] severity filter using upstream's
// default severity scale (overridable subset). The annot-prefix (-p), all-fields
// expansion (-A), duplicate (-d), drop/keep-sites (-x/-X) and allow-undef (-u)
// options are honoured.
//
// Parts of upstream that require the bcftools filter and convert engines, the
// gene-list machinery, or canonical-transcript expression selection are not
// reproduced byte-for-byte and are reported as a clean unsupported Init error
// rather than emitting silently divergent output: -i/-e expressions, -g/--gene-list,
// the EXPRESSION / primary / pick / mane transcript selectors, --columns-types FILE,
// and -S FILE severity overrides.
//
// Unlike the stats/contrast/split plugins (whose -i/-e are now wired to the
// shared filter engine as a plain VCF site/sample pre-filter), split-vep's -i/-e
// cannot be a simple pre-filter: upstream registers the expanded per-transcript
// CSQ subfields as synthetic INFO tags in the OUTPUT header (filter_init runs on
// args->hdr_out) so an expression like -i 'gnomAD_AF<0.1' or -e 'IMPACT="LOW"'
// resolves names that exist only after split-vep's per-transcript expansion, not
// in the input VCF columns. Re-enabling it would require teaching the filter
// engine about those derived columns, which is out of scope here, so it stays a
// clean unsupported error. The query format engine supports the common
// %CHROM/%POS/%ID/%REF/%ALT/%QUAL/%FILTER/%INFO-tag and %CSQ-subfield tokens with
// literal text, \t and \n; other convert directives are rejected.
package bcftools

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("split-vep", func() NativePlugin { return &splitVepPlugin{} })
}

// split-vep column value types, mirroring the BCF_HT_* values used upstream.
const (
	svTypeStr = iota
	svTypeInt
	svTypeReal
)

// split-vep transcript-selection modes (the supported subset of SELECT_TR_*).
const (
	svTrAll = iota
	svTrWorst
)

// svSelectAny is the sentinel min/max severity meaning "any consequence",
// matching SELECT_CSQ_ANY (-1) upstream.
const svSelectAny = -1

// svAnnot is one requested output column: the subfield name, its index in the
// CSQ subfield list (-1 means the whole transcript string), the output INFO tag
// and the value type.
type svAnnot struct {
	field string
	idx   int
	tag   string
	typ   int
}

// splitVepPlugin implements split-vep. It is a fullPlugin: it owns input reading
// and either text or VCF/BCF output, because the text (-f) mode does not fit the
// per-record re-emit pipeline.
type splitVepPlugin struct {
	annotation  string // -a INFO tag to parse ("" => auto CSQ/BCSQ)
	annotPrefix string // -p prefix prepended to all subfield names
	columnStr   string // -c raw column spec
	formatStr   string // -f format string
	allFields   string // -A delimiter (expands %CSQ)
	selectStr   string // -s raw select spec
	duplicate   bool   // -d
	listFields  bool   // -l
	printHeader int    // -H count (1 => header, 2 => omit indices)
	allowUndef  bool   // -u
	dropSites   int    // -1 unset, 0 keep, 1 drop

	// Resolved state (filled by resolve()).
	vepTag    string
	fields    []string       // subfield names in header order
	field2idx map[string]int // first index wins
	csqIdx    int
	annots    []svAnnot
	selTr     int
	minSev    int
	maxSev    int

	// Severity scale: scale entries are substrings, csq2severity maps a token to
	// its tier. Both grow lazily as new consequence tokens are seen.
	scale       []string
	csq2sev     map[string]int
	formatItems []svFmtItem
}

// Name returns the plugin name.
func (p *splitVepPlugin) Name() string { return "split-vep" }

// About returns the one-line description, matching split-vep.c about() verbatim.
func (p *splitVepPlugin) About() string {
	return "Query structured annotations such as INFO/CSQ produced by VEP of INFO/BCSQ produced by bctools/csq.\n"
}

// RunStyle reports that split-vep is a run()-style plugin (options precede the
// input file with no `--` separator).
func (p *splitVepPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether a split-vep flag consumes the following token,
// so the host can split the input-file positional out of the plugin options.
func (p *splitVepPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-a", "--annotation", "-A", "--all-fields", "-c", "--columns",
		"--columns-types", "-f", "--format", "-g", "--gene-list",
		"--gene-list-fields", "-p", "--annot-prefix", "-s", "--select",
		"-S", "--severity", "-i", "--include", "-e", "--exclude",
		"-o", "--output", "-O", "--output-type", "-r", "--regions",
		"-R", "--regions-file", "-t", "--targets", "-T", "--targets-file",
		"--regions-overlap", "--targets-overlap", "-v", "--verbosity":
		return true
	}
	return false
}

// Init is not used: split-vep is a fullPlugin, so runNativePlugin delegates the
// whole invocation to RunFull (which parses the plugin arguments itself). Init
// is retained to satisfy the NativePlugin interface.
func (p *splitVepPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return hdr, nil
}

// parseArgs walks the plugin argv, rejecting the modes that cannot be made
// byte-exact without the upstream filter/convert/gene engines.
func (p *splitVepPlugin) parseArgs(args []string) error {
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("split-vep: option %q requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-a", "--annotation":
			v, err := next()
			if err != nil {
				return err
			}
			p.annotation = v
		case "-p", "--annot-prefix":
			v, err := next()
			if err != nil {
				return err
			}
			p.annotPrefix = v
		case "-c", "--columns":
			v, err := next()
			if err != nil {
				return err
			}
			p.columnStr = v
		case "-f", "--format":
			v, err := next()
			if err != nil {
				return err
			}
			p.formatStr = v
		case "-A", "--all-fields":
			v, err := next()
			if err != nil {
				return err
			}
			p.allFields = v
		case "-s", "--select":
			v, err := next()
			if err != nil {
				return err
			}
			p.selectStr = v
		case "-d", "--duplicate":
			p.duplicate = true
		case "-l", "--list":
			p.listFields = true
		case "-H", "--print-header":
			p.printHeader++
		case "-u", "--allow-undef-tags":
			p.allowUndef = true
		case "-x", "--drop-sites":
			p.dropSites = 1
		case "-X", "--keep-sites":
			p.dropSites = 0
		case "-i", "--include", "-e", "--exclude":
			return fmt.Errorf("split-vep: the -i/-e filter expressions filter over split-vep's expanded per-transcript CSQ subfields (registered on the output header by filter_init), not just the input VCF columns, and so are not supported in the native plugin; run upstream bcftools for that")
		case "-g", "--gene-list", "--gene-list-fields":
			return fmt.Errorf("split-vep: the -g/--gene-list gene-restriction machinery is not supported in the native plugin; run upstream bcftools for that")
		case "-S", "--severity":
			return fmt.Errorf("split-vep: -S/--severity scale overrides are not supported in the native plugin; the built-in default scale is used")
		case "--columns-types":
			return fmt.Errorf("split-vep: --columns-types overrides are not supported in the native plugin")
		default:
			return fmt.Errorf("split-vep: unsupported option %q", a)
		}
	}
	if p.printHeader > 0 && p.formatStr == "" {
		return fmt.Errorf("split-vep: -H/--print-header requires -f/--format")
	}
	if p.allFields != "" && p.formatStr == "" {
		return fmt.Errorf("split-vep: -A/--all-fields requires -f/--format")
	}
	return nil
}

// Process is unused: split-vep is a fullPlugin.
func (p *splitVepPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *splitVepPlugin) Destroy() error { return nil }

// RunFull reads the input, resolves the CSQ header, then dispatches to the list,
// text, or annotate path. It owns its own output writer.
func (p *splitVepPlugin) RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	p.dropSites = -1
	p.csqIdx = -1
	if err := p.parseArgs(opts.Args); err != nil {
		return err
	}
	hdr, variants, err := readPluginInput(opts, stderr)
	if err != nil {
		return err
	}
	if err := p.resolveHeader(hdr); err != nil {
		return err
	}

	if p.listFields {
		var b strings.Builder
		for i, f := range p.fields {
			fmt.Fprintf(&b, "%d\t%s\n", i, f)
		}
		_, err := io.WriteString(out, b.String())
		return err
	}

	if err := p.resolveSelect(); err != nil {
		return err
	}

	if p.formatStr != "" {
		return p.runFormat(variants, out)
	}
	return p.runAnnotate(opts, hdr, variants, out)
}

// resolveHeader discovers the CSQ/BCSQ tag, parses its "Format: a|b|c"
// Description into subfield names, and builds the name->index map.
func (p *splitVepPlugin) resolveHeader(hdr *vcf.Header) error {
	tag := p.annotation
	if tag == "" {
		hasCSQ := infoHeaderExists(hdr, "CSQ")
		hasBCSQ := infoHeaderExists(hdr, "BCSQ")
		switch {
		case hasCSQ:
			tag = "CSQ"
		case hasBCSQ:
			tag = "BCSQ"
		default:
			return fmt.Errorf("split-vep: expected INFO/CSQ or INFO/BCSQ annotation")
		}
	}
	desc := infoHeaderDescription(hdr, tag)
	if desc == "" {
		return fmt.Errorf("split-vep: could not find INFO/%s in the header", tag)
	}
	idx := strings.Index(desc, "Format: ")
	if idx < 0 {
		return fmt.Errorf("split-vep: could not parse the %s Format from the header", tag)
	}
	format := desc[idx+len("Format: "):]
	format = strings.TrimSuffix(format, `"`)
	p.vepTag = tag
	p.field2idx = map[string]int{}
	for _, raw := range splitVepFields(format) {
		name := p.sanitizeField(raw)
		if _, ok := p.field2idx[name]; !ok {
			p.field2idx[name] = len(p.fields)
		}
		p.fields = append(p.fields, name)
	}
	if ci, ok := p.field2idx["Consequence"]; ok {
		p.csqIdx = ci
	} else if p.annotPrefix != "" {
		if ci, ok := p.field2idx[p.annotPrefix+"Consequence"]; ok {
			p.csqIdx = ci
		}
	}
	return nil
}

// splitVepFields tokenizes the "a|b(1-based)|c" format string: fields are
// separated by '|', and a '(' ends the current field name (the bracketed text up
// to the next '|' is discarded), matching upstream's tokenizer.
func splitVepFields(format string) []string {
	var fields []string
	i := 0
	for i < len(format) {
		start := i
		for i < len(format) && format[i] != '|' && format[i] != '"' && format[i] != '(' {
			i++
		}
		fields = append(fields, format[start:i])
		for i < len(format) && format[i] != '|' {
			i++
		}
		if i < len(format) {
			i++ // skip '|'
		}
	}
	return fields
}

// sanitizeField applies the annot-prefix then upstream's sanitize_field_name:
// "1000G" is exempt; a leading '.' or digit gets a '_' prefix; any character
// outside [A-Za-z0-9_.] becomes '_'.
func (p *splitVepPlugin) sanitizeField(name string) string {
	if p.annotPrefix != "" {
		name = p.annotPrefix + name
	}
	if name == "1000G" {
		return name
	}
	var b strings.Builder
	if len(name) > 0 && (name[0] == '.' || (name[0] >= '0' && name[0] <= '9')) {
		b.WriteByte('_')
	}
	for i := 0; i < len(name); i++ {
		c := name[i]
		if (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' {
			b.WriteByte(c)
		} else {
			b.WriteByte('_')
		}
	}
	return b.String()
}

// resolveSelect parses -s into the transcript mode and the severity range, then
// resolves the column spec and the format string (whichever applies). It also
// finalises drop_sites.
func (p *splitVepPlugin) resolveSelect() error {
	p.initSeverityScale()
	p.selTr = svTrAll
	p.minSev, p.maxSev = svSelectAny, svSelectAny

	sel := p.selectStr
	if sel == "" {
		sel = "all:any"
	}
	parts := strings.Split(sel, ":")
	selTr, selCsq := "all", "any"
	if len(parts) > 0 && parts[0] != "" {
		selTr = parts[0]
	}
	if len(parts) > 1 && parts[1] != "" {
		selCsq = parts[1]
	}
	if len(parts) > 2 && parts[2] != "" && parts[2] != "all" {
		return fmt.Errorf("split-vep: the PRN selection %q is not supported in the native plugin (only :all)", parts[2])
	}
	switch selTr {
	case "all":
		p.selTr = svTrAll
	case "worst":
		p.selTr = svTrWorst
	default:
		return fmt.Errorf("split-vep: the transcript selector %q (primary/pick/mane/EXPRESSION) is not supported in the native plugin", selTr)
	}
	if selCsq != "any" {
		modifier := byte('=')
		term := selCsq
		if n := len(term); n > 0 && (term[n-1] == '+' || term[n-1] == '-') {
			modifier = term[n-1]
			term = term[:n-1]
		}
		sev, ok := p.csq2sev[strings.ToLower(term)]
		if !ok {
			return fmt.Errorf("split-vep: unknown consequence %q (see bcftools +split-vep -S -)", term)
		}
		switch modifier {
		case '=':
			p.minSev, p.maxSev = sev, sev
		case '+':
			p.minSev, p.maxSev = sev, int(^uint(0)>>1)
		case '-':
			p.minSev, p.maxSev = 0, sev
		}
	}

	if p.formatStr != "" {
		if err := p.resolveFormat(); err != nil {
			return err
		}
	} else {
		if err := p.resolveColumns(); err != nil {
			return err
		}
		if len(p.annots) == 0 && p.minSev == svSelectAny {
			return fmt.Errorf("split-vep: nothing selected to do; use -c, -f or -s")
		}
	}
	if p.dropSites == -1 {
		if p.formatStr != "" {
			p.dropSites = 1
		} else {
			p.dropSites = 0
		}
	}
	return nil
}

// resolveColumns parses -c into the ordered annot list, resolving names/indexes,
// ranges and the optional :TYPE suffix, and assigns default types.
func (p *splitVepPlugin) resolveColumns() error {
	if p.columnStr == "" {
		return nil
	}
	seen := map[string]bool{}
	add := func(idx, typ int) error {
		if idx < 0 || idx >= len(p.fields) {
			return fmt.Errorf("split-vep: column index %d out of range", idx)
		}
		field := p.fields[idx]
		if seen[field] {
			return nil
		}
		seen[field] = true
		t := typ
		if t < 0 {
			t = p.defaultColumnType(field)
		}
		tag := field
		p.annots = append(p.annots, svAnnot{field: field, idx: idx, tag: tag, typ: t})
		return nil
	}
	for _, item := range strings.Split(p.columnStr, ",") {
		item = strings.TrimSpace(item)
		if item == "" {
			continue
		}
		if item == "-" {
			for i := range p.fields {
				if err := add(i, -1); err != nil {
					return err
				}
			}
			continue
		}
		typ := -1
		name := item
		if colon := strings.LastIndex(item, ":"); colon >= 0 {
			suffix := item[colon+1:]
			if t, ok := parseSvType(suffix); ok {
				typ = t
				name = item[:colon]
			}
		}
		if idx, ok := p.field2idx[name]; ok {
			if err := add(idx, typ); err != nil {
				return err
			}
			continue
		}
		if sname := p.sanitizeField(strings.TrimPrefix(name, p.annotPrefix)); sname != name {
			if idx, ok := p.field2idx[sname]; ok {
				if err := add(idx, typ); err != nil {
					return err
				}
				continue
			}
		}
		// numeric index or range "a-b"
		if lo, hi, ok := parseIndexRange(name); ok {
			for i := lo; i <= hi; i++ {
				if err := add(i, typ); err != nil {
					return err
				}
			}
			continue
		}
		return fmt.Errorf("split-vep: no such column %q", name)
	}
	return nil
}

// parseSvType maps a :TYPE suffix to a column type, matching upstream's accepted
// spellings.
func parseSvType(s string) (int, bool) {
	switch strings.ToLower(s) {
	case "string", "str":
		return svTypeStr, true
	case "integer", "int":
		return svTypeInt, true
	case "float", "real":
		return svTypeReal, true
	}
	return 0, false
}

// parseIndexRange parses "N" or "N-M" into an inclusive index range.
func parseIndexRange(s string) (int, int, bool) {
	if dash := strings.IndexByte(s, '-'); dash > 0 {
		lo, err1 := strconv.Atoi(s[:dash])
		hi, err2 := strconv.Atoi(s[dash+1:])
		if err1 == nil && err2 == nil {
			return lo, hi, true
		}
		return 0, 0, false
	}
	if n, err := strconv.Atoi(s); err == nil {
		return n, n, true
	}
	return 0, 0, false
}

// defaultColumnType returns the default value type for a subfield name, porting
// the default_column_types presets (first matching pattern wins, else String).
func (p *splitVepPlugin) defaultColumnType(field string) int {
	// Strip the annot-prefix when matching, as upstream matches the raw VEP
	// field name against the presets.
	name := strings.TrimPrefix(field, p.annotPrefix)
	switch name {
	case "DISTANCE", "STRAND", "TSL", "GENE_PHENO", "HGVS_OFFSET",
		"MOTIF_POS", "existing_InFrame_oORFs", "existing_OutOfFrame_oORFs",
		"existing_uORFs":
		return svTypeInt
	case "AF", "MAX_AF", "MOTIF_SCORE_CHANGE":
		return svTypeReal
	}
	if strings.HasSuffix(name, "_POPS") {
		return svTypeStr
	}
	if strings.HasSuffix(name, "_AF") {
		return svTypeReal
	}
	if strings.HasPrefix(name, "MAX_AF_") {
		return svTypeReal
	}
	if strings.HasPrefix(name, "SpliceAI_pred_DP_") {
		return svTypeInt
	}
	if strings.HasPrefix(name, "SpliceAI_pred_DS_") {
		return svTypeReal
	}
	return svTypeStr
}

// readPluginInput reads opts.InputFile (with regions) into a header and a slice
// of variants, reusing the host's ViewFile normalisation. It is shared by the
// fullPlugin-style native plugins.
func readPluginInput(opts PluginOptions, stderr io.Writer) (*vcf.Header, []*vcf.Variant, error) {
	regions := append([]string{}, opts.Regions...)
	if opts.RegionsFile != "" {
		regs, rerr := LoadRegionsFile(opts.RegionsFile)
		if rerr != nil {
			return nil, nil, rerr
		}
		regions = append(regions, regs...)
	}
	input := opts.InputFile
	if input == "" {
		input = "-"
	}
	var buf strings.Builder
	if _, err := ViewFile(input, &stringWriter{&buf}, ViewOptions{OutputFormat: OutputVCF, Regions: regions}, stderr); err != nil {
		return nil, nil, fmt.Errorf("reading plugin input: %w", err)
	}
	r := vcf.NewReader(strings.NewReader(buf.String()))
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, nil, err
	}
	variants, err := r.ReadAll()
	if err != nil {
		return nil, nil, fmt.Errorf("malformed VCF input: %w", err)
	}
	return hdr, variants, nil
}

// stringWriter adapts a strings.Builder to io.Writer for ViewFile.
type stringWriter struct{ b *strings.Builder }

func (w *stringWriter) Write(p []byte) (int, error) { return w.b.Write(p) }
