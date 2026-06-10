// Package bcftools — advanced `bcftools annotate` features.
//
// This file implements the option-tail of `bcftools annotate` that the v1
// port deferred (tracked in docs/PARITY_ROADMAP.md):
//
//   - `--set-id [+]<FORMAT>` — set the ID column from a `bcftools query`-like
//     macro string. The supported macros mirror the non-FORMAT subset of
//     upstream's `convert.c`: %CHROM, %POS, %POS0, %END, %END0, %ID, %REF,
//     %ALT, %FIRST_ALT, %QUAL, %FILTER, %TYPE and %INFO/TAG (also written
//     %TAG). A leading `+` only sets the ID when it is currently missing.
//   - `--merge-logic <tag:logic>` — how to combine values when several
//     annotation rows overlap one VCF record (first|append|append-missing|
//     unique|sum|avg|min|max).
//   - `--min-overlap <ann:vcf>` — minimum reciprocal overlap (as a fraction
//     of the annotation interval, the VCF record, or both) required before an
//     overlapping annotation row is applied.
//   - `--pair-logic <mode>` — allele-pairing mode (exact|some|all|...) used
//     when matching records against a VCF annotation source.
//   - `--single-overlaps` — apply only the first overlapping annotation row
//     (the default with merge-logic unset is also "first", but this flag
//     additionally disables the multi-interval bookkeeping).
//   - `--rename-annots <file>` — rename INFO/FORMAT/FILTER tags per a map
//     file (`TYPE/old<whitespace>new`).
//
// The upstream implementation in vcfannotate.c is built on htslib's streaming
// synced-reader and regidx interval index. This port reads all records up
// front (see annotate.go) and matches them against in-memory annotation
// tables, so the code here mirrors the *observable* upstream semantics rather
// than its internal data structures.
package bcftools

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// MergeLogic enumerates the multi-overlap merge strategies of
// `--merge-logic`. The zero value MergeFirst matches upstream's default
// (keep the first overlapping annotation row, discard the rest).
type MergeLogic int

// Merge-logic strategies, mirroring the MM_* constants of vcfannotate.c.
const (
	// MergeFirst keeps the value from the first overlapping annotation row.
	MergeFirst MergeLogic = iota
	// MergeAppend concatenates every overlapping value (Number=. tags only).
	MergeAppend
	// MergeAppendMissing is MergeAppend but also carries over missing values.
	MergeAppendMissing
	// MergeUnique appends only values not already seen.
	MergeUnique
	// MergeSum sums the numeric values of overlapping rows.
	MergeSum
	// MergeAvg averages the numeric values of overlapping rows.
	MergeAvg
	// MergeMin keeps the minimum numeric value across overlapping rows.
	MergeMin
	// MergeMax keeps the maximum numeric value across overlapping rows.
	MergeMax
)

// PairLogic enumerates the `--pair-logic` allele-pairing modes used when the
// annotation source is a VCF/BCF.
type PairLogic int

// Pair-logic modes, mirroring the BCF_SR_PAIR_* flags upstream sets.
const (
	// PairSome (the default) matches when the REFs are compatible and at
	// least one ALT allele is shared.
	PairSome PairLogic = iota
	// PairExact requires REF and the full ALT set to be identical.
	PairExact
	// PairAll matches on REF alone, ignoring ALT.
	PairAll
	// PairSNPs matches SNP records sharing a REF.
	PairSNPs
	// PairIndels matches indel records sharing a REF.
	PairIndels
	// PairBoth matches SNP or indel records sharing a REF.
	PairBoth
	// PairID additionally requires the ID columns to match.
	PairID
)

// ParseMergeLogic parses one `--merge-logic TAG:LOGIC` argument into the tag
// name and the MergeLogic value. The recognised logic keywords are
// first, append, append-missing, unique, sum, avg, min and max (matching
// vcfannotate.c). Multiple comma-separated tag:logic pairs may be supplied;
// callers split on commas before invoking this for each pair.
func ParseMergeLogic(spec string) (tag string, logic MergeLogic, err error) {
	i := strings.LastIndexByte(spec, ':')
	if i < 0 {
		return "", MergeFirst, fmt.Errorf("could not parse --merge-logic %q, expected TAG:LOGIC", spec)
	}
	tag = strings.TrimSpace(spec[:i])
	name := strings.TrimSpace(spec[i+1:])
	if tag == "" {
		return "", MergeFirst, fmt.Errorf("could not parse --merge-logic %q, empty tag", spec)
	}
	switch strings.ToLower(name) {
	case "first":
		logic = MergeFirst
	case "append":
		logic = MergeAppend
	case "append-missing":
		logic = MergeAppendMissing
	case "unique":
		logic = MergeUnique
	case "sum":
		logic = MergeSum
	case "avg":
		logic = MergeAvg
	case "min":
		logic = MergeMin
	case "max":
		logic = MergeMax
	default:
		return "", MergeFirst, fmt.Errorf("could not parse --merge-logic %q, the logic %q is not recognised", spec, name)
	}
	return tag, logic, nil
}

// ParseMergeLogicSpec parses a full `--merge-logic` argument that may carry
// several comma-separated TAG:LOGIC pairs into a tag->logic map.
func ParseMergeLogicSpec(spec string) (map[string]MergeLogic, error) {
	if spec == "" {
		return nil, nil
	}
	out := map[string]MergeLogic{}
	for _, part := range strings.Split(spec, ",") {
		part = strings.TrimSpace(part)
		if part == "" {
			continue
		}
		tag, logic, err := ParseMergeLogic(part)
		if err != nil {
			return nil, err
		}
		out[tag] = logic
	}
	return out, nil
}

// ParsePairLogic parses a `--pair-logic` argument. The recognised modes are
// snps, indels, both, all/any, some, none/exact and id, matching the switch
// in vcfannotate.c.
func ParsePairLogic(s string) (PairLogic, error) {
	switch s {
	case "snps":
		return PairSNPs, nil
	case "indels":
		return PairIndels, nil
	case "both":
		return PairBoth, nil
	case "all", "any":
		return PairAll, nil
	case "some":
		return PairSome, nil
	case "none", "exact":
		return PairExact, nil
	case "id":
		return PairID, nil
	default:
		return PairSome, fmt.Errorf("the --pair-logic string %q is not recognised", s)
	}
}

// MinOverlap holds the parsed `--min-overlap ANN:VCF` fractions. A zero value
// for either side disables that side's check (matching upstream where 0
// means "no requirement").
type MinOverlap struct {
	// Ann is the minimum overlap as a fraction of the annotation interval.
	Ann float64
	// Vcf is the minimum overlap as a fraction of the VCF record.
	Vcf float64
}

// ParseMinOverlap parses a `--min-overlap` argument. Accepted shapes are
// `ANN`, `ANN:VCF` and `:VCF`, each value in the range [0,1], mirroring the
// strtod-based parser in vcfannotate.c.
func ParseMinOverlap(s string) (MinOverlap, error) {
	var mo MinOverlap
	if s == "" {
		return mo, nil
	}
	annStr, vcfStr := s, ""
	hasVcf := false
	if i := strings.IndexByte(s, ':'); i >= 0 {
		annStr = s[:i]
		vcfStr = s[i+1:]
		hasVcf = true
	}
	if annStr != "" {
		v, err := strconv.ParseFloat(annStr, 64)
		if err != nil || v < 0 || v > 1 {
			return mo, fmt.Errorf("could not parse \"--min-overlap %s\", expected value(s) between 0-1", s)
		}
		mo.Ann = v
	}
	if hasVcf && vcfStr != "" {
		v, err := strconv.ParseFloat(vcfStr, 64)
		if err != nil || v < 0 || v > 1 {
			return mo, fmt.Errorf("could not parse \"--min-overlap %s\", expected value(s) between 0-1", s)
		}
		mo.Vcf = v
	}
	return mo, nil
}

// passesMinOverlap reports whether an annotation interval [annBeg,annEnd]
// overlapping a VCF record [vcfBeg,vcfEnd] (all 1-based, inclusive) meets the
// configured minimum-overlap thresholds. It mirrors the isec/len arithmetic
// of annotate_from_regidx.
func passesMinOverlap(mo MinOverlap, annBeg, annEnd, vcfBeg, vcfEnd int) bool {
	lenAnn := annEnd - annBeg + 1
	lenVcf := vcfEnd - vcfBeg + 1
	hi := annEnd
	if vcfEnd < hi {
		hi = vcfEnd
	}
	lo := annBeg
	if vcfBeg > lo {
		lo = vcfBeg
	}
	isec := hi - lo + 1
	if isec <= 0 {
		return false
	}
	if mo.Ann > 0 && lenAnn > 0 && mo.Ann > float64(isec)/float64(lenAnn) {
		return false
	}
	if mo.Vcf > 0 && lenVcf > 0 && mo.Vcf > float64(isec)/float64(lenVcf) {
		return false
	}
	return true
}

// mergeAccumulator collects the candidate values for a single INFO tag across
// every overlapping annotation row, then reduces them per the chosen logic.
type mergeAccumulator struct {
	logic  MergeLogic
	values []string
	seen   map[string]bool
	// isInt is set when the tag is declared Type=Integer, so the numeric
	// reductions (notably avg) truncate to an integer like upstream.
	isInt bool
}

func newMergeAccumulator(logic MergeLogic, isInt bool) *mergeAccumulator {
	return &mergeAccumulator{logic: logic, seen: map[string]bool{}, isInt: isInt}
}

// add records one overlapping value. Empty/"." values are dropped unless the
// append-missing logic is in force.
func (m *mergeAccumulator) add(value string) {
	if value == "" || value == "." {
		if m.logic != MergeAppendMissing {
			return
		}
	}
	if m.logic == MergeUnique {
		if m.seen[value] {
			return
		}
		m.seen[value] = true
	}
	m.values = append(m.values, value)
	if m.logic == MergeFirst && len(m.values) > 1 {
		m.values = m.values[:1]
	}
}

// reduce produces the final INFO value string, or ("", false) when there is
// nothing to write.
func (m *mergeAccumulator) reduce() (string, bool) {
	if len(m.values) == 0 {
		return "", false
	}
	switch m.logic {
	case MergeFirst:
		return m.values[0], true
	case MergeAppend, MergeAppendMissing, MergeUnique:
		return strings.Join(m.values, ","), true
	case MergeSum, MergeAvg, MergeMin, MergeMax:
		return reduceNumeric(m.logic, m.values, m.isInt)
	default:
		return m.values[0], true
	}
}

// reduceNumeric applies sum/avg/min/max to the numeric components of each
// value. Multi-valued (comma-separated) annotation values are reduced
// element-wise, matching upstream which operates per array slot. When isInt
// is set (the tag is Type=Integer) the avg is truncated toward zero like
// upstream's integer setter.
func reduceNumeric(logic MergeLogic, values []string, isInt bool) (string, bool) {
	var acc []float64
	var cnt []int
	allInt := true
	for _, v := range values {
		for i, tok := range strings.Split(v, ",") {
			f, err := strconv.ParseFloat(tok, 64)
			if err != nil {
				continue
			}
			if f != float64(int64(f)) {
				allInt = false
			}
			for len(acc) <= i {
				acc = append(acc, 0)
				cnt = append(cnt, 0)
			}
			switch logic {
			case MergeSum, MergeAvg:
				acc[i] += f
			case MergeMin:
				if cnt[i] == 0 || f < acc[i] {
					acc[i] = f
				}
			case MergeMax:
				if cnt[i] == 0 || f > acc[i] {
					acc[i] = f
				}
			}
			cnt[i]++
		}
	}
	if len(acc) == 0 {
		return "", false
	}
	out := make([]string, len(acc))
	for i := range acc {
		val := acc[i]
		if logic == MergeAvg && cnt[i] > 0 {
			val = val / float64(cnt[i])
		}
		out[i] = formatMergeNumber(val, allInt, isInt)
	}
	return strings.Join(out, ","), true
}

// formatMergeNumber renders a merged numeric value. For an Integer-typed tag
// the value is truncated toward zero. Otherwise integral values print without
// a decimal point and fractional values use the shortest %g representation.
func formatMergeNumber(v float64, allInt, isInt bool) string {
	if isInt {
		return strconv.FormatInt(int64(v), 10)
	}
	if allInt && v == float64(int64(v)) {
		return strconv.FormatInt(int64(v), 10)
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
}

// renameTag is one parsed entry of a `--rename-annots` map.
type renameTag struct {
	// Type is one of "INFO", "FORMAT" or "FILTER".
	Type string
	// Old and New are the source and destination tag names (without the
	// TYPE/ prefix).
	Old string
	New string
}

// loadRenameAnnots reads a `--rename-annots` map file. Each non-empty line is
// `TYPE/old<whitespace>new`, where TYPE is one of INFO, FORMAT/FMT or FILTER
// (case-insensitive). It mirrors rename_annots_init/rename_annots_init1.
func loadRenameAnnots(path string) ([]renameTag, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var out []renameTag
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("could not parse: %s", line)
		}
		rt, err := parseRenameEntry(fields[0], fields[1])
		if err != nil {
			return nil, err
		}
		out = append(out, rt)
	}
	return out, sc.Err()
}

// parseRenameEntry splits the `TYPE/old` source and the `[TYPE/]new`
// destination into a renameTag, validating that a typed destination matches
// the source type.
func parseRenameEntry(src, dst string) (renameTag, error) {
	typ, old, ok := splitTypePrefix(src)
	if !ok {
		return renameTag{}, fmt.Errorf("cannot rename %q: missing INFO/FORMAT/FILTER prefix", src)
	}
	if dtyp, dnew, dok := splitTypePrefix(dst); dok {
		if dtyp != typ {
			return renameTag{}, fmt.Errorf("cannot transfer %s to %s", old, dtyp)
		}
		dst = dnew
	}
	return renameTag{Type: typ, Old: old, New: dst}, nil
}

// splitTypePrefix recognises a leading INFO/, FORMAT/, FMT/ or FILTER/ prefix
// (case-insensitive) and returns the canonical type plus the remainder.
func splitTypePrefix(s string) (typ, rest string, ok bool) {
	lower := strings.ToLower(s)
	switch {
	case strings.HasPrefix(lower, "info/"):
		return "INFO", s[5:], true
	case strings.HasPrefix(lower, "format/"):
		return "FORMAT", s[7:], true
	case strings.HasPrefix(lower, "fmt/"):
		return "FORMAT", s[4:], true
	case strings.HasPrefix(lower, "filter/"):
		return "FILTER", s[7:], true
	default:
		return "", "", false
	}
}

// applyRenameAnnots renames INFO/FORMAT/FILTER tags in the header and in every
// record, mirroring rename_annots. The header line for the old tag has its
// ID= rewritten to the new name (the rest of the line is preserved).
func applyRenameAnnots(recs []*vcf.Variant, hdr *vcf.Header, maps []renameTag) {
	for _, m := range maps {
		if m.Old == m.New {
			continue
		}
		// Rewrite the matching ##TYPE=<ID=old,...> header line.
		for i, line := range hdr.MetaInfo {
			k, id := structuredID(line)
			if k == m.Type && id == m.Old {
				hdr.MetaInfo[i] = replaceMetaID(line, m.Old, m.New)
			}
		}
		switch m.Type {
		case "INFO":
			for _, v := range recs {
				renameInfoKey(v, m.Old, m.New)
			}
		case "FORMAT":
			for _, v := range recs {
				renameFormatKey(v, m.Old, m.New)
			}
		case "FILTER":
			for _, v := range recs {
				renameFilterName(v, m.Old, m.New)
			}
		}
	}
}

// replaceMetaID rewrites the `ID=old` token of a structured meta line to
// `ID=new`, leaving every other field untouched.
func replaceMetaID(line, oldID, newID string) string {
	return strings.Replace(line, "ID="+oldID, "ID="+newID, 1)
}

// renameInfoKey renames an INFO tag on one record, preserving its position in
// the INFO column.
func renameInfoKey(v *vcf.Variant, oldTag, newTag string) {
	if v.Info == nil {
		return
	}
	val, ok := v.Info[oldTag]
	if !ok {
		return
	}
	delete(v.Info, oldTag)
	v.Info[newTag] = val
	for i, k := range v.InfoOrder {
		if k == oldTag {
			v.InfoOrder[i] = newTag
		}
	}
}

// renameFormatKey renames a FORMAT tag on one record (both the FORMAT list and
// every sample's data map).
func renameFormatKey(v *vcf.Variant, oldTag, newTag string) {
	for i, k := range v.Format {
		if k == oldTag {
			v.Format[i] = newTag
		}
	}
	for si := range v.Samples {
		if v.Samples[si].Data == nil {
			continue
		}
		if val, ok := v.Samples[si].Data[oldTag]; ok {
			delete(v.Samples[si].Data, oldTag)
			v.Samples[si].Data[newTag] = val
		}
	}
}

// renameFilterName renames a FILTER entry on one record.
func renameFilterName(v *vcf.Variant, oldName, newName string) {
	for i, f := range v.Filter {
		if f == oldName {
			v.Filter[i] = newName
		}
	}
}

// setIDProgram is a parsed `--set-id` macro string.
type setIDProgram struct {
	// onlyIfEmpty is set when the format string began with `+`.
	onlyIfEmpty bool
	tokens      []setIDToken
}

// setIDToken is one element of a setIDProgram: either literal text or a macro.
type setIDToken struct {
	// literal holds verbatim text when macro is empty.
	literal string
	// macro is the macro name (e.g. "CHROM", "INFO/AC"); empty for literals.
	macro string
}

// ParseSetID parses a `--set-id [+]<FORMAT>` argument into a setIDProgram.
// The macro grammar mirrors the non-FORMAT subset of convert.c: a `%` starts
// a macro, `%INFO/TAG` (or a bare `%TAG`) reads an INFO field, and `\t`,
// `\n` plus `\<char>` escapes are honoured in literal runs.
func ParseSetID(format string) (*setIDProgram, error) {
	p := &setIDProgram{}
	if strings.HasPrefix(format, "+") {
		p.onlyIfEmpty = true
		format = format[1:]
	}
	var lit strings.Builder
	flush := func() {
		if lit.Len() > 0 {
			p.tokens = append(p.tokens, setIDToken{literal: lit.String()})
			lit.Reset()
		}
	}
	i := 0
	for i < len(format) {
		c := format[i]
		switch c {
		case '%':
			flush()
			macro, next, err := parseSetIDMacro(format, i+1)
			if err != nil {
				return nil, err
			}
			p.tokens = append(p.tokens, setIDToken{macro: macro})
			i = next
		case '\\':
			i++
			if i >= len(format) {
				lit.WriteByte('\\')
				break
			}
			switch format[i] {
			case 'n':
				lit.WriteByte('\n')
			case 't':
				lit.WriteByte('\t')
			default:
				lit.WriteByte(format[i])
			}
			i++
		default:
			lit.WriteByte(c)
			i++
		}
	}
	flush()
	return p, nil
}

// parseSetIDMacro reads one macro name starting at offset i (just past the
// `%`). It returns the macro name, the offset just past it, and an error if
// the name is empty. `%INFO/TAG` is normalised to the macro name "INFO/TAG";
// a bare `%TAG/SUB` (vcf-column form) keeps its slash.
func parseSetIDMacro(format string, i int) (string, int, error) {
	// `%/TAG` is the explicit vcf-column form in convert.c; accept it by
	// skipping the leading slash and treating the remainder as a column.
	vcfColumn := false
	if i < len(format) && format[i] == '/' {
		vcfColumn = true
		i++
	}
	start := i
	for i < len(format) && (isWordByte(format[i])) {
		i++
	}
	name := format[start:i]
	if name == "" {
		return "", i, fmt.Errorf("could not parse format string near %%")
	}
	// %INFO/TAG — consume the slash and the tag name.
	if !vcfColumn && name == "INFO" && i < len(format) && format[i] == '/' {
		i++
		ts := i
		for i < len(format) && isWordByte(format[i]) {
			i++
		}
		tag := format[ts:i]
		if tag == "" {
			return "", i, fmt.Errorf("could not parse format string near %%INFO/")
		}
		return "INFO/" + tag, i, nil
	}
	return name, i, nil
}

// isWordByte reports whether b is part of a macro/tag name (alphanumeric,
// underscore or dot), matching convert.c's parse_tag loop.
func isWordByte(b byte) bool {
	return b == '_' || b == '.' ||
		(b >= '0' && b <= '9') ||
		(b >= 'a' && b <= 'z') ||
		(b >= 'A' && b <= 'Z')
}

// expand renders the program for one record. It returns ("", false) when the
// program emits nothing (an empty result), in which case the ID is left
// unchanged (matching upstream's `if ( args->tmpks.l )` guard).
func (p *setIDProgram) expand(v *vcf.Variant) (string, bool) {
	var b strings.Builder
	for _, tok := range p.tokens {
		if tok.macro == "" {
			b.WriteString(tok.literal)
			continue
		}
		b.WriteString(expandSetIDMacro(tok.macro, v))
	}
	s := b.String()
	if s == "" {
		return "", false
	}
	return s, true
}

// expandSetIDMacro renders a single macro for one record.
func expandSetIDMacro(macro string, v *vcf.Variant) string {
	if strings.HasPrefix(macro, "INFO/") {
		return infoMacroValue(v, macro[len("INFO/"):])
	}
	switch macro {
	case "CHROM":
		return v.Chrom
	case "POS":
		return strconv.Itoa(v.Pos)
	case "POS0":
		return strconv.Itoa(v.Pos - 1)
	case "END":
		return strconv.Itoa(variantEnd(v))
	case "END0":
		return strconv.Itoa(variantEnd(v) - 1)
	case "ID":
		return v.ID
	case "REF":
		return v.Ref
	case "ALT":
		return altMacroValue(v)
	case "FIRST_ALT":
		return firstAltValue(v)
	case "QUAL":
		return qualMacroValue(v)
	case "FILTER":
		return filterMacroValue(v)
	case "TYPE":
		return typeMacroValue(v)
	default:
		// A bare %TAG is interpreted as an INFO tag, matching the
		// fall-through in convert.c.
		return infoMacroValue(v, macro)
	}
}

// infoMacroValue returns the value of an INFO tag, "." if absent, or "1" for
// a flag (a present tag with no value), matching process_info.
func infoMacroValue(v *vcf.Variant, tag string) string {
	if v.Info == nil {
		return "."
	}
	val, ok := v.Info[tag]
	if !ok {
		return "."
	}
	if val == "" {
		return "1"
	}
	return val
}

// altMacroValue returns the comma-joined ALT alleles, or "." for a
// no-ALT record, matching process_alt.
func altMacroValue(v *vcf.Variant) string {
	if len(v.Alt) == 0 || (len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "")) {
		return "."
	}
	return strings.Join(v.Alt, ",")
}

// firstAltValue returns the first ALT allele, or "." when there is none,
// matching process_first_alt.
func firstAltValue(v *vcf.Variant) string {
	if len(v.Alt) == 0 || v.Alt[0] == "" || v.Alt[0] == "." {
		return "."
	}
	return v.Alt[0]
}

// qualMacroValue renders QUAL, "." when missing, matching process_qual.
func qualMacroValue(v *vcf.Variant) string {
	if v.Qual < 0 {
		return "."
	}
	if v.Qual == float64(int64(v.Qual)) {
		return strconv.FormatInt(int64(v.Qual), 10)
	}
	return strconv.FormatFloat(v.Qual, 'g', -1, 64)
}

// filterMacroValue renders the FILTER column joined by ';', "." when empty,
// matching process_filter.
func filterMacroValue(v *vcf.Variant) string {
	if len(v.Filter) == 0 {
		return "."
	}
	parts := make([]string, 0, len(v.Filter))
	for _, f := range v.Filter {
		if f == "" {
			continue
		}
		parts = append(parts, f)
	}
	if len(parts) == 0 {
		return "."
	}
	return strings.Join(parts, ";")
}

// variantEnd returns the 1-based inclusive end coordinate (POS + rlen - 1),
// where rlen is the REF length, matching line->pos+line->rlen with the END
// processor's pos+rlen reading (which uses the 0-based POS internally).
func variantEnd(v *vcf.Variant) int {
	rlen := len(v.Ref)
	if rlen < 1 {
		rlen = 1
	}
	return v.Pos + rlen - 1
}

// typeMacroValue renders %TYPE as the comma-separated set of allele types,
// in upstream's fixed order (SNP,MNP,INDEL,OTHER,BND,OVERLAP), or REF for a
// no-variant record. It mirrors process_type over bcf_get_variant_types.
func typeMacroValue(v *vcf.Variant) string {
	mask := 0
	for _, alt := range v.Alt {
		mask |= variantTypeBit(v.Ref, alt)
	}
	if mask == 0 {
		return "REF"
	}
	var parts []string
	if mask&vtSNP != 0 {
		parts = append(parts, "SNP")
	}
	if mask&vtMNP != 0 {
		parts = append(parts, "MNP")
	}
	if mask&vtINDEL != 0 {
		parts = append(parts, "INDEL")
	}
	if mask&vtOTHER != 0 {
		parts = append(parts, "OTHER")
	}
	if mask&vtBND != 0 {
		parts = append(parts, "BND")
	}
	if mask&vtOVERLAP != 0 {
		parts = append(parts, "OVERLAP")
	}
	if len(parts) == 0 {
		return "REF"
	}
	return strings.Join(parts, ",")
}

// Variant-type bits, matching the VCF_* flags used by process_type (after
// masking to ORIG_VAR_TYPES).
const (
	vtSNP = 1 << iota
	vtMNP
	vtINDEL
	vtOTHER
	vtBND
	vtOVERLAP
)

// variantTypeBit classifies one REF/ALT pair, porting htslib's
// bcf_set_variant_type. Only the original (masked) types are returned so the
// result feeds process_type directly.
func variantTypeBit(ref, alt string) int {
	if alt == "*" {
		return vtOVERLAP
	}
	if ref == "" || alt == "" || alt == "." {
		return 0 // REF
	}
	// Single-base REF and ALT: the common SNP/REF case. Upstream compares the
	// bases case-sensitively here (`*ref==*alt`), so "A" vs "a" is a SNP.
	if len(ref) == 1 && len(alt) == 1 {
		if alt == "." || ref[0] == alt[0] || alt == "X" {
			return 0
		}
		return vtSNP
	}
	if alt[0] == '<' {
		return vtOTHER
	}
	if alt[0] == ']' || alt[0] == '[' {
		return vtBND
	}

	r, a := 0, 0
	for r < len(ref) && a < len(alt) && equalFold1(ref[r], alt[a]) {
		r++
		a++
	}
	switch {
	case a < len(alt) && r == len(ref):
		if alt[len(alt)-1] == ']' || alt[len(alt)-1] == '[' {
			return vtBND
		}
		return vtINDEL
	case r < len(ref) && a == len(alt):
		return vtINDEL
	case r == len(ref) && a == len(alt):
		return 0
	}

	re, ae := len(ref)-1, len(alt)-1
	if alt[ae] == ']' || alt[ae] == '[' {
		return vtBND
	}
	for re > r && ae > a && equalFold1(ref[re], alt[ae]) {
		re--
		ae--
	}
	if ae == a {
		if re == r {
			return vtSNP
		}
		if equalFold1(ref[re], alt[ae]) {
			return vtINDEL
		}
		return vtOTHER
	}
	if re == r {
		if equalFold1(ref[re], alt[ae]) {
			return vtINDEL
		}
		return vtOTHER
	}
	if re-r == ae-a {
		return vtMNP
	}
	return vtOTHER
}

// equalFold1 compares two bytes case-insensitively (ASCII only), matching
// htslib's toupper_c usage in bcf_set_variant_type.
func equalFold1(x, y byte) bool {
	if x >= 'a' && x <= 'z' {
		x -= 'a' - 'A'
	}
	if y >= 'a' && y <= 'z' {
		y -= 'a' - 'A'
	}
	return x == y
}

// SR-type bits, mirroring the SR_REF/SR_SNP/SR_INDEL/SR_OTHER classes that
// bcf_sr_sort.c collapses the full VCF_* mask into for pairing.
const (
	srREF = 1 << iota
	srSNP
	srINDEL
	srOTHER
)

// recordSRType reduces a record's per-allele variant types to the SR_*
// pairing class used by htslib's synced reader. VCF_MNP folds into SR_SNP and
// a no-ALT record is SR_REF, matching bcf_sr_sort's classification.
func recordSRType(v *vcf.Variant) int {
	mask := 0
	for _, alt := range v.Alt {
		mask |= variantTypeBit(v.Ref, alt)
	}
	if mask == 0 {
		return srREF
	}
	t := 0
	if mask&(vtSNP|vtMNP) != 0 {
		t |= srSNP
	}
	if mask&vtINDEL != 0 {
		t |= srINDEL
	}
	if mask&vtOTHER != 0 {
		t |= srOTHER
	}
	if t == 0 {
		// BND/OVERLAP-only records have no SR class; treat as OTHER so they
		// only pair under exact/shared-allele matches.
		t = srOTHER
	}
	return t
}

// alleleTokens builds the REF>ALT token set htslib uses for allele matching
// (see bcf_sr_sort.c). Tokens are upper-cased so comparison is
// case-insensitive, matching strncasecmp. A no-ALT record yields "REF>.".
func alleleTokens(v *vcf.Variant) []string {
	if len(v.Alt) == 0 {
		return []string{strings.ToUpper(v.Ref) + ">."}
	}
	out := make([]string, 0, len(v.Alt))
	for _, alt := range v.Alt {
		out = append(out, strings.ToUpper(v.Ref)+">"+strings.ToUpper(alt))
	}
	return out
}

// tokenSetsEqual reports whether two token slices contain the same multiset of
// tokens (order-insensitive), the multi_is_exact criterion.
func tokenSetsEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for _, ta := range a {
		found := false
		for _, tb := range b {
			if ta == tb {
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	return true
}

// tokensShared reports whether the two token slices share at least one token,
// the multi_is_subset criterion.
func tokensShared(a, b []string) bool {
	for _, ta := range a {
		for _, tb := range b {
			if ta == tb {
				return true
			}
		}
	}
	return false
}

// srScore returns the pairing score for two SR types under the given mode,
// porting the SR_SCORE table built by bcf_sr_init_scores. A score of 0 means
// the two records cannot be paired on type alone.
func srScore(logic PairLogic, ti, tj int) int {
	// PAIR_ANY enables every type combination.
	if logic == PairAll {
		return 1
	}
	best := 0
	consider := func(a, b, score int) {
		if score > best && ti&a != 0 && tj&b != 0 {
			best = score
		}
	}
	switch logic {
	case PairSNPs:
		consider(srSNP, srREF, 2)
		consider(srREF, srSNP, 2)
	case PairIndels:
		consider(srINDEL, srREF, 2)
		consider(srREF, srINDEL, 2)
	case PairBoth:
		consider(srSNP, srSNP, 3)
		consider(srINDEL, srINDEL, 3)
		consider(srSNP, srREF, 2)
		consider(srREF, srSNP, 2)
		consider(srINDEL, srREF, 2)
		consider(srREF, srINDEL, 2)
	}
	return best
}

// pairLogicMatches reports whether a VCF annotation record src may be paired
// with the input record v under the given pair-logic mode. It mirrors the
// single-varset case of pairing_score in htslib's bcf_sr_sort.c: an exact or
// shared-allele match always pairs, otherwise the SR type-score for the mode
// must be non-zero.
func pairLogicMatches(v, src *vcf.Variant, logic PairLogic) bool {
	if v.Chrom != src.Chrom || v.Pos != src.Pos {
		return false
	}
	ti := recordSRType(v)
	tj := recordSRType(src)
	at := alleleTokens(v)
	bt := alleleTokens(src)

	if logic == PairExact {
		if ti != tj {
			return false
		}
		return tokenSetsEqual(at, bt)
	}

	if logic == PairID && !idsMatch(v, src) {
		return false
	}

	// Exact allele-string match or a shared allele always pairs.
	if ti == tj && tokenSetsEqual(at, bt) {
		return true
	}
	if ti&tj != 0 && tokensShared(at, bt) {
		return true
	}

	// PairID falls back to the "some" type semantics once IDs agree, which
	// (like the default) has no extra type score and so only pairs on a
	// shared allele — handled above.
	if logic == PairID || logic == PairSome {
		return false
	}

	return srScore(logic, ti, tj) > 0
}

// idsMatch reports whether the input ID contains the annotation ID as a
// semicolon-delimited member, mirroring upstream's ID overlap requirement.
func idsMatch(v, src *vcf.Variant) bool {
	if src.ID == "" || src.ID == "." {
		return true
	}
	for _, id := range strings.Split(v.ID, ";") {
		if id == src.ID {
			return true
		}
	}
	return false
}
