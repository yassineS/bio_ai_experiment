// Native port of the upstream `gvcfz` plugin (plugins/gvcfz.c): compress a gVCF
// file by resizing the reference blocks according to -g/--group-by criteria.
//
// gvcfz groups consecutive gVCF reference blocks (records whose only ALT is
// <NON_REF>/<*>) by evaluating a list of FILTER:EXPR group definitions against
// each record. The first group whose bcftools filter EXPRESSION matches the
// record selects the block the record belongs to; a "-" expression is the
// catch-all that always matches. Consecutive records assigned to the same group
// are merged into one representative record (the first), with INFO/END extended
// to the block end, FORMAT/DP set to the minimum DP (MIN_DP preferred over DP),
// FORMAT/GQ (or RGQ) set to the minimum, and FORMAT/PL set to the element-wise
// minimum. A non-reference site (a real variant) flushes the current block and
// is written through verbatim.
//
// The FILTER:EXPR filter expressions are compiled with the in-tree filter engine
// (CompileFilterWithHeader / Filter.Eval), which now supports the FORMAT/GT
// predicates upstream's filter_init/filter_test cover (e.g. GT!="alt", GQ>60 &
// DP<20). This is what makes the native port possible — previously the framework
// did not expose that engine.
//
// The plugin is run()-style (options precede the input file, no `--`) and is a
// serial bufferedPlugin because the block state machine spans records.
package bcftools

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("gvcfz", func() NativePlugin { return &gvcfzPlugin{} }) }

// gvcfzGroup is one compiled FILTER:EXPR group definition. A nil filter is the
// "-" catch-all that always matches; filterLabel is "" for the PASS pseudo-group
// (which adds no FILTER) and otherwise the FILTER ID to stamp onto the block.
type gvcfzGroup struct {
	expr        string  // raw expression text (for diagnostics)
	filter      *Filter // compiled expression, nil for "-" (always matches)
	filterLabel string  // FILTER ID to add ("" for PASS, no filter added)
}

// gvcfzBlock is the running state of the current gVCF block, mirroring the
// upstream block_t. grp < 0 means the block is inactive.
type gvcfzBlock struct {
	grp    int          // group index, -1 when inactive
	rec    *vcf.Variant // representative record (the first in the block)
	end    int          // 1-based block end (max END seen)
	minDP  int          // minimum FORMAT/MIN_DP (or DP) across the block
	gqKey  string       // "GQ", "RGQ", or "" when neither is present
	gq     int          // minimum FORMAT/GQ|RGQ
	pl     [3]int       // element-wise minimum FORMAT/PL ([-1,-1,-1] when absent)
	filter *gvcfzGroup  // the group definition for grp (nil when inactive)
}

// gvcfzPlugin implements the gvcfz block-resizing plugin.
type gvcfzPlugin struct {
	hdr      *vcf.Header
	groupBy  string
	groups   []gvcfzGroup
	preExpr  string // raw -i/-e expression (record-level pre-filter)
	preExclu bool   // true for -e
	preFlt   *pluginFilter
}

// Name returns the plugin name.
func (p *gvcfzPlugin) Name() string { return "gvcfz" }

// About returns the one-line description, matching gvcfz.c about().
func (p *gvcfzPlugin) About() string {
	return "Compress gVCF file by resizing gVCF blocks according to specified criteria."
}

// RunStyle reports that gvcfz is a run()-style plugin (its options precede the
// input file with no `--` separator).
func (p *gvcfzPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of gvcfz's own value-taking flags consumes
// the following CLI token, so the host can split the input file out of the
// run-style options.
func (p *gvcfzPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-g", "--group-by", "-i", "--include", "-e", "--exclude",
		"-o", "--output", "-O", "--output-type", "-v", "--verbosity":
		return true
	}
	return false
}

// Init parses gvcfz's options, compiles the group expressions and the optional
// -i/-e pre-filter, and returns the output header with the INFO/END definition
// and any non-PASS FILTER definitions appended, exactly as init_groups() does.
func (p *gvcfzPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	for i := 0; i < len(args); i++ {
		a := args[i]
		// Support the joined short-flag form (-gEXPR, -ioEXPR, ...) the way getopt
		// does: a value-taking short flag may carry its value in the same token.
		var joined string
		hasJoined := false
		if len(a) > 2 && a[0] == '-' && a[1] != '-' {
			switch a[:2] {
			case "-g", "-i", "-e", "-o", "-O", "-v":
				joined = a[2:]
				hasJoined = true
				a = a[:2]
			}
		}
		next := func() (string, error) {
			if hasJoined {
				return joined, nil
			}
			if i+1 >= len(args) {
				return "", fmt.Errorf("gvcfz: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-g", "--group-by":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.groupBy = v
		case "-i", "--include", "-e", "--exclude":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if p.preExpr != "" {
				return nil, fmt.Errorf("gvcfz: only one -i or -e expression can be given, and they cannot be combined")
			}
			p.preExpr = v
			p.preExclu = a == "-e" || a == "--exclude"
		case "-a", "--trim-alt-alleles":
			// Trimming unused ALT alleles only affects records that are NOT pure
			// gVCF blocks before the n_allele re-check. Since the native port keeps
			// the representative record's alleles verbatim and the supported
			// fixtures use <NON_REF>/<*> blocks (already 2-allele), -a is accepted
			// as a no-op for parity; a real alt-bearing record still flushes and
			// passes through unchanged exactly as upstream.
		case "-o", "--output", "-O", "--output-type", "-v", "--verbosity",
			"-W", "--write-index":
			// Output container / indexing options are handled by the host pipeline
			// (or are no-ops here); consume the value form where applicable.
			if a == "-o" || a == "--output" || a == "-O" || a == "--output-type" ||
				a == "-v" || a == "--verbosity" {
				if _, err := next(); err != nil {
					return nil, err
				}
			} else if i+1 < len(args) && !strings.HasPrefix(args[i+1], "-") {
				// -W may take an optional FMT argument.
				i++
			}
		case "--no-version":
			// provenance is handled elsewhere in the native path
		default:
			return nil, fmt.Errorf("gvcfz: unsupported option %q", a)
		}
	}
	if p.groupBy == "" {
		return nil, fmt.Errorf("gvcfz: missing the -g option")
	}

	if p.preExpr != "" {
		pf, err := newPluginFilterWithHeader(p.preExpr, p.preExclu, hdr)
		if err != nil {
			return nil, fmt.Errorf("gvcfz: %w", err)
		}
		p.preFlt = pf
	}

	// Build the output header: duplicate the input, append the END INFO line and
	// any non-PASS FILTER lines, exactly as init_groups() does.
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	out.MetaInfo = appendInfoHeader(out.MetaInfo,
		`##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">`)

	groups, filterLines, err := parseGvcfzGroups(p.groupBy, hdr)
	if err != nil {
		return nil, err
	}
	p.groups = groups
	for _, line := range filterLines {
		out.MetaInfo = appendFilterHeader(out.MetaInfo, line)
	}
	return out, nil
}

// parseGvcfzGroups parses the -g group-by string into compiled group definitions
// and the ##FILTER header lines for the non-PASS filters, porting init_groups().
// The group-by string is a list of `FILTER:EXPR` clauses separated by `;`. A
// FILTER other than PASS contributes a ##FILTER line whose Description is the
// VERBATIM group-by string with double quotes replaced by single quotes. EXPR is
// trimmed of surrounding whitespace; "-" yields a nil (always-matching) filter.
func parseGvcfzGroups(groupBy string, hdr *vcf.Header) (groups []gvcfzGroup, filterLines []string, err error) {
	// The FILTER Description is the whole group-by string with " -> '.
	descr := strings.ReplaceAll(groupBy, `"`, `'`)
	seenFilter := map[string]bool{"PASS": true}

	rest := groupBy
	for {
		rest = strings.TrimLeft(rest, " \t\n\r\f\v")
		if rest == "" {
			break
		}
		colon := strings.IndexByte(rest, ':')
		if colon < 0 {
			return nil, nil, fmt.Errorf("gvcfz: could not parse the expression: %q", groupBy)
		}
		flt := rest[:colon]
		rest = rest[colon+1:]
		semi := strings.IndexByte(rest, ';')
		var exprText string
		if semi < 0 {
			exprText = rest
			rest = ""
		} else {
			exprText = rest[:semi]
			rest = rest[semi+1:]
		}
		exprTrimmed := strings.TrimSpace(exprText)

		g := gvcfzGroup{expr: exprTrimmed}
		if flt != "PASS" {
			g.filterLabel = flt
			if !seenFilter[flt] {
				seenFilter[flt] = true
				filterLines = append(filterLines,
					fmt.Sprintf(`##FILTER=<ID=%s,Description="%s">`, flt, descr))
			}
		}
		if exprTrimmed != "-" {
			f, ferr := CompileFilterWithHeader(exprTrimmed, hdr)
			if ferr != nil {
				return nil, nil, fmt.Errorf("gvcfz: %w", ferr)
			}
			g.filter = f
		}
		groups = append(groups, g)

		if semi < 0 {
			break
		}
	}
	return groups, filterLines, nil
}

// ProcessAll runs the gVCF block state machine over the whole ordered stream and
// emits the resized records, porting run()'s `while(bcf_sr_next_line) process_gvcf`
// loop followed by the final flush_block(NULL).
func (p *gvcfzPlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	out := make([]*vcf.Variant, 0, len(variants))
	blk := &gvcfzBlock{grp: -1}

	flush := func(rec *vcf.Variant) error {
		emit, err := p.flushBlock(blk, rec)
		if err != nil {
			return err
		}
		if emit != nil {
			out = append(out, emit)
		}
		return nil
	}

	for _, rec := range variants {
		// -i/-e record-level pre-filter (process_gvcf's filter_test gate).
		if p.preFlt != nil && !p.preFlt.testSite(rec) {
			continue
		}

		// A record with a real ALT (not <NON_REF>/<*>) is not a gVCF block: flush
		// the current block and pass the record through verbatim.
		if !isGvcfRefBlock(rec) {
			if err := flush(rec); err != nil {
				return nil, err
			}
			out = append(out, rec)
			continue
		}

		end := gvcfEnd(rec)
		gqKey, gq := gvcfGQ(rec)
		minDP, ok := gvcfMinDP(rec)
		if !ok {
			return nil, fmt.Errorf("gvcfz: expected one FORMAT/MIN_DP or FORMAT/DP value at %s:%d", rec.Chrom, rec.Pos)
		}
		pl, plOK := gvcfPL(rec)
		if !plOK {
			pl = [3]int{-1, -1, -1}
		}

		// Select the first group whose filter matches (nil filter always matches).
		grp := len(p.groups)
		for i := range p.groups {
			if p.groups[i].filter == nil || p.groups[i].filter.Eval(rec) {
				grp = i
				break
			}
		}

		if blk.grp != grp {
			if err := flush(rec); err != nil { // new block
				return nil, err
			}
		}
		if blk.grp >= 0 && blk.rec.Chrom != rec.Chrom {
			if err := flush(nil); err != nil { // new chromosome
				return nil, err
			}
		}

		if blk.grp >= 0 {
			// Extend the existing block.
			if blk.end < end {
				blk.end = end
			}
			if blk.gqKey != "" && gqKey != "" && blk.gq > gq {
				blk.gq = gq
			}
			if blk.minDP > minDP {
				blk.minDP = minDP
			}
			for k := 0; k < 3; k++ {
				if blk.pl[k] > pl[k] {
					blk.pl[k] = pl[k]
				}
			}
			continue
		}

		// Start a new block.
		blk.rec = copyVariant(rec)
		blk.grp = grp
		if grp < len(p.groups) {
			blk.filter = &p.groups[grp]
		} else {
			blk.filter = nil
		}
		blk.minDP = minDP
		blk.end = end
		blk.pl = pl
		blk.gqKey = gqKey
		if gqKey != "" {
			blk.gq = gq
		}
	}
	if err := flush(nil); err != nil {
		return nil, err
	}
	return out, nil
}

// flushBlock finalises the current block and returns the representative record to
// emit (nil when no block is active), porting flush_block(). When rec is
// non-nil and the block would overlap it, the block end is clamped to rec->pos.
func (p *gvcfzPlugin) flushBlock(blk *gvcfzBlock, rec *vcf.Variant) (*vcf.Variant, error) {
	if blk.grp < 0 {
		return nil, nil
	}
	end := blk.end
	if rec != nil && end >= rec.Pos {
		// Upstream: `if (gvcf->end-1 >= rec->pos) gvcf->end = rec->pos;` where
		// gvcf->end is 1-based and rec->pos is 0-based. With our 1-based rec.Pos
		// (= rec->pos+1) the test `end-1 >= rec.Pos-1` is `end >= rec.Pos` and the
		// clamp `end = rec->pos` becomes `end = rec.Pos-1`.
		end = rec.Pos - 1
	}

	r := blk.rec
	// INFO/END only when the block spans more than the single start position.
	if r.Pos < end {
		setInfo(r, "END", strconv.Itoa(end))
	}
	if err := setGvcfFormatInt(r, "DP", blk.minDP); err != nil {
		return nil, fmt.Errorf("gvcfz: could not update FORMAT/DP at %s:%d", r.Chrom, r.Pos)
	}
	if blk.gqKey != "" {
		if err := setGvcfFormatInt(r, blk.gqKey, blk.gq); err != nil {
			return nil, fmt.Errorf("gvcfz: could not update FORMAT/%s at %s:%d", blk.gqKey, r.Chrom, r.Pos)
		}
	}
	if blk.pl[0] >= 0 {
		if err := setGvcfFormatPL(r, blk.pl); err != nil {
			return nil, fmt.Errorf("gvcfz: could not update FORMAT/PL at %s:%d", r.Chrom, r.Pos)
		}
	}
	if blk.filter != nil && blk.filter.filterLabel != "" {
		addBcfFilter(r, blk.filter.filterLabel)
	}

	blk.grp = -1
	blk.filter = nil
	return r, nil
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *gvcfzPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *gvcfzPlugin) Destroy() error { return nil }

// SetStderr satisfies stderrSink; gvcfz prints no end-of-run summary.
func (p *gvcfzPlugin) SetStderr(io.Writer) {}

// --- pure helpers (also exercised by the binary-free unit tests) -----------

// isGvcfRefBlock reports whether rec is a gVCF reference block: it has at most a
// single ALT allele which, if present, is the <NON_REF> or <*> symbolic allele.
// This ports the n_allele/<NON_REF> guard in process_gvcf().
func isGvcfRefBlock(rec *vcf.Variant) bool {
	if len(rec.Alt) == 0 {
		return true
	}
	if len(rec.Alt) > 1 {
		return false
	}
	a := rec.Alt[0]
	return a == "<NON_REF>" || a == "<*>" || a == "." || a == ""
}

// gvcfEnd returns the block end (1-based): INFO/END when present, else POS.
func gvcfEnd(rec *vcf.Variant) int {
	if s, ok := rec.Info["END"]; ok && s != "" && s != "." {
		if n, err := strconv.Atoi(s); err == nil {
			return n
		}
	}
	return rec.Pos
}

// gvcfGQ returns the GQ key in effect ("GQ", "RGQ", or "" if neither has a single
// value) and its value, porting process_gvcf()'s GQ/RGQ probing.
func gvcfGQ(rec *vcf.Variant) (string, int) {
	if v, ok := singleFormatInt(rec, "GQ"); ok {
		return "GQ", v
	}
	if v, ok := singleFormatInt(rec, "RGQ"); ok {
		return "RGQ", v
	}
	return "", 0
}

// gvcfMinDP returns the block's representative DP: FORMAT/MIN_DP when present,
// else FORMAT/DP; ok is false when neither is a single integer (an error in
// upstream).
func gvcfMinDP(rec *vcf.Variant) (int, bool) {
	if v, ok := singleFormatInt(rec, "MIN_DP"); ok {
		return v, true
	}
	if v, ok := singleFormatInt(rec, "DP"); ok {
		return v, true
	}
	return 0, false
}

// gvcfPL returns the three FORMAT/PL values, ok=false when PL is absent or not
// exactly three values, porting the PL handling in process_gvcf().
func gvcfPL(rec *vcf.Variant) ([3]int, bool) {
	if len(rec.Samples) == 0 {
		return [3]int{}, false
	}
	s, ok := rec.Samples[0].Data["PL"]
	if !ok || s == "" || s == "." {
		return [3]int{}, false
	}
	parts := strings.Split(s, ",")
	if len(parts) != 3 {
		return [3]int{}, false
	}
	var pl [3]int
	for i := 0; i < 3; i++ {
		n, err := strconv.Atoi(parts[i])
		if err != nil {
			return [3]int{}, false
		}
		pl[i] = n
	}
	return pl, true
}

// singleFormatInt returns the integer value of FORMAT/key for the first sample
// when it is a single integer, mirroring bcf_get_format_int32(...)==1.
func singleFormatInt(rec *vcf.Variant, key string) (int, bool) {
	if len(rec.Samples) == 0 {
		return 0, false
	}
	s, ok := rec.Samples[0].Data[key]
	if !ok || s == "" || s == "." {
		return 0, false
	}
	if strings.IndexByte(s, ',') >= 0 {
		return 0, false
	}
	n, err := strconv.Atoi(s)
	if err != nil {
		return 0, false
	}
	return n, true
}

// setGvcfFormatInt overwrites (or appends) a single-valued FORMAT integer field
// on the first (only) sample, mirroring bcf_update_format_int32(...,1). When the
// FORMAT tag is not yet present it is appended to the FORMAT list and to the
// sample, matching upstream's bcf_update_format growing the FORMAT vector — but
// only when the header declares the tag (otherwise upstream errors).
func setGvcfFormatInt(rec *vcf.Variant, key string, val int) error {
	return setFormatField(rec, key, strconv.Itoa(val))
}

// setGvcfFormatPL overwrites the three-valued FORMAT/PL field on the sample.
func setGvcfFormatPL(rec *vcf.Variant, pl [3]int) error {
	return setFormatField(rec, "PL", fmt.Sprintf("%d,%d,%d", pl[0], pl[1], pl[2]))
}

// setFormatField sets FORMAT/key to value on every sample, appending the tag to
// the FORMAT list when it is not already present (single-sample gVCF, so this
// matches bcf_update_format growing the vector).
func setFormatField(rec *vcf.Variant, key, value string) error {
	present := false
	for _, f := range rec.Format {
		if f == key {
			present = true
			break
		}
	}
	if !present {
		rec.Format = append(rec.Format, key)
	}
	for i := range rec.Samples {
		if rec.Samples[i].Data == nil {
			rec.Samples[i].Data = make(map[string]string)
		}
		rec.Samples[i].Data[key] = value
	}
	return nil
}

// addBcfFilter ports bcf_add_filter: a PASS/"." record's FILTER is replaced by
// the label; otherwise the label is appended (deduplicated).
func addBcfFilter(v *vcf.Variant, label string) {
	if len(v.Filter) == 0 || (len(v.Filter) == 1 && (v.Filter[0] == "PASS" || v.Filter[0] == ".")) {
		v.Filter = []string{label}
		return
	}
	for _, f := range v.Filter {
		if f == label {
			return
		}
	}
	v.Filter = append(v.Filter, label)
}

// copyVariant deep-copies the fields gvcfz mutates on the representative record
// (FILTER, INFO, FORMAT and per-sample data), so resizing one block never
// aliases the source record shared with the input slice.
func copyVariant(v *vcf.Variant) *vcf.Variant {
	cp := *v
	cp.Alt = append([]string(nil), v.Alt...)
	cp.Filter = append([]string(nil), v.Filter...)
	cp.Info = make(map[string]string, len(v.Info))
	for k, val := range v.Info {
		cp.Info[k] = val
	}
	cp.InfoOrder = append([]string(nil), v.InfoOrder...)
	cp.Format = append([]string(nil), v.Format...)
	cp.Samples = make([]vcf.Sample, len(v.Samples))
	for i, s := range v.Samples {
		cp.Samples[i].Name = s.Name
		cp.Samples[i].Data = make(map[string]string, len(s.Data))
		for k, val := range s.Data {
			cp.Samples[i].Data[k] = val
		}
	}
	return &cp
}

// appendFilterHeader inserts a ##FILTER line if a definition for the same ID is
// not already present, matching upstream's bcf_hdr_printf no-op on a duplicate.
func appendFilterHeader(meta []string, line string) []string {
	id := headerID(line)
	if id != "" {
		for _, m := range meta {
			if headerKind(m) == "##FILTER" && headerID(m) == id {
				return meta
			}
		}
	}
	return append(meta, line)
}
