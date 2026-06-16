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
// Transcript and consequence selection (-s/--select) supports the full upstream
// surface: the all/worst transcript modes, the EXPRESSION selectors
// (primary => CANONICAL=YES, pick => PICK=1, mane => MANE_SELECT!="", or an
// arbitrary <FIELD><OP><VALUE> with the =, !=, ~, !~ operators), the
// any / term[+|-] severity filter using the default severity scale (overridable
// via -S FILE), and the PRN qualifier (:all / :worst, where :worst rewrites the
// printed Consequence to its single worst term). Gene restriction
// (-g/--gene-list with the optional leading "+" prioritise mode and
// --gene-list-fields) and the file-based column-type overrides
// (--columns-types -|FILE) are also reproduced. The annot-prefix (-p),
// all-fields expansion (-A), duplicate (-d), drop/keep-sites (-x/-X) and
// allow-undef (-u) options are honoured.
//
// The -i/-e filter expressions are supported: they evaluate against the
// expanded per-transcript CSQ subfields. Upstream registers those subfields as
// synthetic INFO tags on the OUTPUT header and runs filter_init(args->hdr_out)
// after parse_column_str, so an expression like -i 'gnomAD_AF<0.1' or
// -e 'IMPACT="LOW"' resolves names that exist only after split-vep's
// per-transcript expansion. The native port mirrors this (see
// native_plugin_splitvep_filter.go): it auto-registers any CSQ subfield the
// expression references as an extra column, compiles the filter against the
// augmented header via the shared filter engine, and applies it where upstream's
// filter_and_output does — per collapsed record, or per transcript with -d. The
// expression sees the same aggregated INFO tags upstream produces, including the
// array-OR ("any element matches") semantics for the non -d case.
//
// The query format engine supports the common
// %CHROM/%POS/%ID/%REF/%ALT/%QUAL/%FILTER/%INFO-tag and %CSQ-subfield tokens with
// literal text, \t and \n; other convert directives are rejected.
package bcftools

import (
	"fmt"
	"io"
	"regexp"
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

// split-vep transcript-selection modes, mirroring SELECT_TR_* upstream.
const (
	svTrAll   = iota // SELECT_TR_ALL: list every transcript
	svTrWorst        // SELECT_TR_WORST: only the worst-consequence transcript
	svTrExpr         // SELECT_TR_EXPR: transcripts matching a <FIELD><OP><VALUE> rule
)

// split-vep PRN (print-consequence) modes, mirroring PRN_CSQ_* upstream.
const (
	svPrnAll   = iota // PRN_CSQ_ALL: print all consequence terms per transcript
	svPrnWorst        // PRN_CSQ_WORST: print only the worst term per transcript
)

// split-vep --select EXPRESSION operators, mirroring TR_OP_* upstream.
const (
	svTrOpEq = iota // =  string equality
	svTrOpNe        // != string inequality
	svTrOpRe        // ~  regex match
	svTrOpNr        // !~ regex non-match
)

// split-vep -g/--gene-list modes, mirroring GENES_* upstream.
const (
	svGenesRestrict   = iota // restrict to transcripts whose gene is listed
	svGenesPrioritize        // keep all, but move listed-gene transcripts first
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

// svTrExprSpec is a parsed --select EXPRESSION (<FIELD><OP><VALUE>), mirroring
// select_tr_t upstream: the CSQ subfield index to consult, the operator, the
// literal comparison value (for =/!=) and the compiled regex (for ~/!~).
type svTrExprSpec struct {
	field string
	idx   int
	op    int
	value string
	regex *regexp.Regexp
}

// svCol2Type is one --columns-types override rule: a regex matched against the
// (unprefixed) VEP subfield name and the value type it maps to. Mirrors
// col2type_t upstream.
type svCol2Type struct {
	regex *regexp.Regexp
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
	severity    string // -S severity scale (- or FILE; "" => default)
	columnTypes string // --columns-types (- or FILE; "" => default presets)
	genesFname  string // -g/--gene-list spec (may have a leading '+')
	geneFields  string // --gene-list-fields LIST (default SYMBOL,Gene,gene)
	duplicate   bool   // -d
	listFields  bool   // -l
	printHeader int    // -H count (1 => header, 2 => omit indices)
	allowUndef  bool   // -u
	dropSites   int    // -1 unset, 0 keep, 1 drop

	// Resolved state (filled by resolve()).
	inHdr     *vcf.Header // the input header (for filter undef-tag resolution)
	vepTag    string
	fields    []string       // subfield names in header order
	field2idx map[string]int // first index wins
	csqIdx    int
	annots    []svAnnot
	selTr     int
	prnCsq    int
	trExpr    svTrExprSpec // valid when selTr == svTrExpr
	minSev    int
	maxSev    int

	// Gene restriction state (filled by initGeneList when -g is given).
	genesMode    int             // svGenesRestrict or svGenesPrioritize
	genes        map[string]bool // the hashed --gene-list gene names
	geneFieldIdx []int           // resolved CSQ subfield indices to match against

	// Column-type overrides (--columns-types FILE). When non-nil these regex
	// rules replace the built-in default_column_types presets. column2typeErr
	// records a deferred parse/read failure surfaced by the first untyped-column
	// lookup, matching upstream where init_column2type errors fatally the first
	// time get_column_type is reached.
	column2type    []svCol2Type
	column2typeErr error

	// Severity scale: scale entries are substrings, csq2severity maps a token to
	// its tier. Both grow lazily as new consequence tokens are seen.
	scale       []string
	csq2sev     map[string]int
	formatItems []svFmtItem

	// annotVals is the per-record value accumulator (one slice of transcript
	// values per annot), reused across records to mirror upstream's annot_t.str.
	annotVals [][]string

	// Filter (-i/-e). filterStr is the raw expression text, filterExclude is
	// true for -e/--exclude, filterSet records that a -i/-e was given (so an
	// empty expression is distinguished from "no filter"). The compiled filter
	// is built against the split-vep-augmented output header so it can reference
	// the per-transcript CSQ subfield columns, matching upstream
	// filter_init(args->hdr_out).
	filterStr     string
	filterExclude bool
	filterSet     bool
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
		case "-i", "--include":
			v, err := next()
			if err != nil {
				return err
			}
			if p.filterSet {
				return fmt.Errorf("Error: only one -i or -e expression can be given, and they cannot be combined")
			}
			p.filterStr = v
			p.filterExclude = false
			p.filterSet = true
		case "-e", "--exclude":
			v, err := next()
			if err != nil {
				return err
			}
			if p.filterSet {
				return fmt.Errorf("Error: only one -i or -e expression can be given, and they cannot be combined")
			}
			p.filterStr = v
			p.filterExclude = true
			p.filterSet = true
		case "-g", "--gene-list":
			v, err := next()
			if err != nil {
				return err
			}
			p.genesFname = v
		case "--gene-list-fields":
			v, err := next()
			if err != nil {
				return err
			}
			p.geneFields = v
		case "-S", "--severity":
			v, err := next()
			if err != nil {
				return err
			}
			p.severity = v
		case "--columns-types":
			v, err := next()
			if err != nil {
				return err
			}
			p.columnTypes = v
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

// svScaleDumpError carries the default severity-scale or column-types dump that
// upstream prints to stderr when invoked with -S - / -S ? / --columns-types -.
// Upstream emits the text via error() and exits non-zero; the host renders the
// Message to stderr and returns a non-zero exit, matching that behaviour.
type svScaleDumpError struct{ msg string }

func (e *svScaleDumpError) Error() string { return e.msg }

// dumpRequest returns a non-nil error carrying the default severity scale or the
// default column-types table when -S -|? or --columns-types - was requested,
// mirroring upstream's pre-init_data checks in run(). The caller renders it to
// stderr and exits non-zero.
func (p *splitVepPlugin) dumpRequest() error {
	if p.severity == "-" || p.severity == "?" {
		return &svScaleDumpError{msg: svDefaultSeverityText}
	}
	if p.columnTypes == "-" {
		return &svScaleDumpError{msg: svDefaultColumnTypesText}
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
	// -S -|? and --columns-types - print the default tables to stderr and exit
	// non-zero, before any input is read, exactly as upstream's run() does.
	if dump := p.dumpRequest(); dump != nil {
		io.WriteString(stderr, dump.Error())
		return dump
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
	p.inHdr = hdr
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
	if err := p.initSeverityScale(); err != nil {
		return err
	}
	p.selTr = svTrAll
	p.prnCsq = svPrnAll
	p.minSev, p.maxSev = svSelectAny, svSelectAny

	sel := p.selectStr
	if sel == "" {
		sel = "all:any"
	}
	parts := strings.Split(sel, ":")
	selTr, selCsq, prnCsq := "all", "any", "all"
	if len(parts) > 0 && parts[0] != "" {
		selTr = parts[0]
	}
	if len(parts) > 1 && parts[1] != "" {
		selCsq = parts[1]
	}
	if len(parts) > 2 && parts[2] != "" {
		prnCsq = parts[2]
	}
	switch strings.ToLower(selTr) {
	case "all":
		p.selTr = svTrAll
	case "worst":
		p.selTr = svTrWorst
	case "primary":
		if err := p.initSelectTrExpr("CANONICAL=YES"); err != nil {
			return err
		}
	case "pick":
		if err := p.initSelectTrExpr("PICK=1"); err != nil {
			return err
		}
	case "mane":
		if err := p.initSelectTrExpr(`MANE_SELECT!=""`); err != nil {
			return err
		}
	default:
		if err := p.initSelectTrExpr(selTr); err != nil {
			return err
		}
	}
	switch strings.ToLower(prnCsq) {
	case "all":
		p.prnCsq = svPrnAll
	case "worst":
		p.prnCsq = svPrnWorst
	default:
		return fmt.Errorf("Error: could not parse \"%s\" in the expression \"%s\"", prnCsq, sel)
	}
	if selCsq != "any" {
		modifier := byte('=')
		term := selCsq
		if n := len(term); n > 0 && (term[n-1] == '+' || term[n-1] == '-') {
			modifier = term[n-1]
			term = term[:n-1]
		}
		// Upstream looks the term up case-sensitively against the lowercased
		// severity-scale keys, so an upper-case term like ":MISSENSE" is rejected.
		sev, ok := p.csq2sev[term]
		if !ok {
			return fmt.Errorf("Error: the consequence \"%s\" is not recognised. Run \"bcftools +split-vep -S ?\" to see the default list.", term)
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
		// Auto-add CSQ subfields referenced by the -i/-e expression as columns so
		// they are registered on the output header and populated as INFO tags
		// before the filter is evaluated, mirroring parse_filter_str.
		if err := p.addFilterColumns(); err != nil {
			return err
		}
	} else {
		if err := p.resolveColumns(); err != nil {
			return err
		}
		if err := p.addFilterColumns(); err != nil {
			return err
		}
	}
	// Surface a deferred --columns-types FILE read/parse failure, which upstream
	// reports fatally the first time a column needs its type resolved.
	if p.column2typeErr != nil {
		return p.column2typeErr
	}

	// init_gene_list runs after the column/filter parsing upstream.
	if p.genesFname != "" {
		if err := p.initGeneList(); err != nil {
			return err
		}
	}

	// The "why not use bcftools view" guard, ported from run() lines 1684-1703.
	// When none of -c/-f selected a column and no severity range is active, the
	// invocation does nothing unless an EXPRESSION transcript selector is given,
	// in which case upstream defaults to drop_sites=1 (keep only sites hitting a
	// matching transcript). -X is then a no-op error.
	if p.formatStr == "" && len(p.annots) == 0 {
		if p.minSev == svSelectAny && p.maxSev == svSelectAny {
			if p.selTr != svTrExpr {
				return fmt.Errorf("Error: none of the -c,-f,-s options was given, why not use \"bcftools view\" instead?")
			}
			if p.dropSites == -1 {
				p.dropSites = 1
			} else if p.dropSites == 0 {
				return fmt.Errorf("Error: the option -X has no effect without -c,-f, why not use \"bcftools view\" instead?")
			}
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

// defaultColumnType returns the value type for a subfield name. It mirrors
// upstream's get_column_type: the column-type regex table (built lazily from the
// --columns-types FILE, or from the built-in default_column_types presets when no
// override is given) is matched in order against the raw (unprefixed) VEP field
// name; the first match wins, and an unmatched name is String. Each pattern is
// anchored with ^...$ exactly as upstream does.
func (p *splitVepPlugin) defaultColumnType(field string) int {
	// Strip the annot-prefix when matching, as upstream matches the raw VEP
	// field name against the presets.
	name := strings.TrimPrefix(field, p.annotPrefix)
	if p.column2type == nil && p.column2typeErr == nil {
		// Build the regex table on first use, exactly when upstream's
		// get_column_type triggers init_column2type. A read/parse failure is
		// recorded and surfaced by the caller (resolveColumns/resolveFormat).
		if err := p.initColumn2Type(); err != nil {
			p.column2typeErr = err
			return svTypeStr
		}
	}
	for _, ct := range p.column2type {
		if ct.regex.MatchString(name) {
			return ct.typ
		}
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
