package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// filterSet is a comma-separated list of FILTER tag names parsed from a CLI
// flag (e.g. --remove-filtered "q10,s50").
type filterSet map[string]struct{}

// parseFilterList splits a comma-separated string into a set of names. Empty
// strings produce an empty set (meaning "no filter to match").
func parseFilterList(s string) filterSet {
	if s == "" {
		return nil
	}
	out := make(filterSet)
	for _, f := range strings.Split(s, ",") {
		f = strings.TrimSpace(f)
		if f != "" {
			out[f] = struct{}{}
		}
	}
	return out
}

// passRemoveFilteredNames returns true when the variant SHOULD be kept under
// --remove-filtered <names>. A site fails (returns false) if any of its
// FILTER entries appears in the set.
func passRemoveFilteredNames(v *vcf.Variant, fset filterSet) bool {
	if len(fset) == 0 {
		return true
	}
	for _, name := range v.Filter {
		if _, hit := fset[name]; hit {
			return false
		}
	}
	return true
}

// passKeepFilteredNames returns true when the variant SHOULD be kept under
// --keep-filtered <names>. A site passes only if at least one of its FILTER
// entries appears in the set. Sites with no FILTER (or FILTER==PASS) are
// rejected, matching upstream behaviour: --keep-filtered selects sites that
// failed at least one of the named filters.
func passKeepFilteredNames(v *vcf.Variant, fset filterSet) bool {
	if len(fset) == 0 {
		return true
	}
	for _, name := range v.Filter {
		if _, hit := fset[name]; hit {
			return true
		}
	}
	return false
}

// infoTagSet is a set of INFO tag names parsed from a CLI flag.
type infoTagSet = filterSet

func parseInfoTagList(s string) infoTagSet {
	return parseFilterList(s)
}

// filterRecodeInfo applies the `--recode-INFO TAG` recode-column selector to
// a recoded variant's INFO map. keepInfo (when non-empty) restricts INFO to
// the listed tags. Mirrors upstream `parameters.cpp:319` (recode_INFO_to_keep).
//
// Returns a fresh map; the caller is responsible for assigning it back.
func filterRecodeInfo(info map[string]string, keepInfo infoTagSet) map[string]string {
	if len(keepInfo) == 0 {
		out := make(map[string]string, len(info))
		for k, v := range info {
			out[k] = v
		}
		return out
	}
	out := make(map[string]string, len(keepInfo))
	for k, v := range info {
		if _, hit := keepInfo[k]; !hit {
			continue
		}
		out[k] = v
	}
	return out
}

// getInfoRunner streams variants and writes <prefix>.INFO with chosen INFO
// tags. The columns are CHROM POS REF ALT followed by one column per
// requested tag, in the order the user supplied them. Missing tags emit ".".
//
// vcftools refers to the input flag as `--get-INFO TAG`, repeatable. Our CLI
// accepts a comma-separated list; both styles map to the same `tags` slice.
type getInfoRunner struct {
	w    *bufio.Writer
	f    io.WriteCloser
	tags []string
}

func newGetInfoRunner(prefix string, tags []string) (*getInfoRunner, error) {
	f, err := iohelper.OpenWriter(prefix + ".INFO")
	if err != nil {
		return nil, fmt.Errorf("opening %s.INFO: %w", prefix, err)
	}
	w := bufio.NewWriter(f)
	// Header.
	if _, err := w.WriteString("CHROM\tPOS\tREF\tALT"); err != nil {
		f.Close()
		return nil, err
	}
	for _, tag := range tags {
		if _, err := w.WriteString("\t" + tag); err != nil {
			f.Close()
			return nil, err
		}
	}
	if _, err := w.WriteString("\n"); err != nil {
		f.Close()
		return nil, err
	}
	return &getInfoRunner{w: w, f: f, tags: append([]string(nil), tags...)}, nil
}

func (g *getInfoRunner) addVariant(v *vcf.Variant) error {
	if g == nil {
		return nil
	}
	altStr := "."
	if len(v.Alt) > 0 {
		altStr = strings.Join(v.Alt, ",")
	}
	if _, err := fmt.Fprintf(g.w, "%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, altStr); err != nil {
		return err
	}
	for _, tag := range g.tags {
		val, ok := v.Info[tag]
		if !ok || val == "" {
			val = "."
		}
		if _, err := g.w.WriteString("\t" + val); err != nil {
			return err
		}
	}
	if _, err := g.w.WriteString("\n"); err != nil {
		return err
	}
	return nil
}

func (g *getInfoRunner) close() error {
	if g == nil {
		return nil
	}
	if err := g.w.Flush(); err != nil {
		g.f.Close()
		return err
	}
	return g.f.Close()
}

// infoHeaderMeta captures the bits of an `##INFO=<ID=X,Number=N,Type=T,...>`
// declaration that the SITE-filter path needs. We only inspect Type and
// ignore Number/Description, mirroring upstream `entry_filters.cpp:1053`,
// which switches on `INFO_map[...].Type != Flag` alone.
type infoHeaderMeta struct {
	ID   string
	Type string
}

// lookupInfoMeta scans a VCF header for an `##INFO=<ID=tag,...>` declaration
// and returns the parsed type field. Returns (meta, false) if no matching
// declaration exists. Used by the `--keep-INFO` site filter to mirror
// upstream's "must be Type=Flag" check at entry_filters.cpp:1053.
//
// The parser is intentionally tolerant — VCF spec allows arbitrary
// field ordering inside the `<...>` envelope, and quoted Description
// fields may embed commas. We pull out the body once and then split on
// commas only outside double-quoted strings.
func lookupInfoMeta(h *vcf.Header, tag string) (infoHeaderMeta, bool) {
	if h == nil {
		return infoHeaderMeta{}, false
	}
	prefix := "##INFO=<"
	for _, line := range h.MetaInfo {
		if !strings.HasPrefix(line, prefix) {
			continue
		}
		// Strip the leading `##INFO=<` and trailing `>` (best-effort).
		body := strings.TrimPrefix(line, prefix)
		if strings.HasSuffix(body, ">") {
			body = body[:len(body)-1]
		}
		fields := splitInfoHeaderFields(body)
		meta := infoHeaderMeta{}
		for _, f := range fields {
			kv := strings.SplitN(f, "=", 2)
			if len(kv) != 2 {
				continue
			}
			key := strings.TrimSpace(kv[0])
			val := strings.TrimSpace(kv[1])
			switch key {
			case "ID":
				meta.ID = val
			case "Type":
				meta.Type = val
			}
		}
		if meta.ID == tag {
			return meta, true
		}
	}
	return infoHeaderMeta{}, false
}

// splitInfoHeaderFields splits the body of an `##INFO=<...>` declaration on
// commas that are NOT inside a double-quoted region. VCF Description fields
// are allowed to embed commas, so a naive `strings.Split(...,",")` would
// chop them apart.
func splitInfoHeaderFields(body string) []string {
	var out []string
	var cur strings.Builder
	inQuote := false
	for i := 0; i < len(body); i++ {
		c := body[i]
		switch {
		case c == '"':
			inQuote = !inQuote
			cur.WriteByte(c)
		case c == ',' && !inQuote:
			out = append(out, cur.String())
			cur.Reset()
		default:
			cur.WriteByte(c)
		}
	}
	if cur.Len() > 0 {
		out = append(out, cur.String())
	}
	return out
}

// validateFlagTypeINFO mirrors upstream's per-tag Type=Flag header check
// in entry_filters.cpp:1053 (keep) and :1072 (remove). It returns a
// descriptive error if any tag in `tags` is undeclared or non-Flag-typed
// in the header. `flagName` is the caller's CLI flag label
// (`"--keep-INFO"` or `"--remove-INFO"`) for inclusion in the error
// message — upstream prefixes its `LOG.error` with the same string. The
// check is header-invariant, so we run it once at Run start rather than
// per-site as upstream does.
func validateFlagTypeINFO(flagName string, tags infoTagSet, h *vcf.Header) error {
	if len(tags) == 0 {
		return nil
	}
	for tag := range tags {
		meta, ok := lookupInfoMeta(h, tag)
		if !ok {
			return fmt.Errorf("%s: INFO tag %q is not declared in the VCF header", flagName, tag)
		}
		if !strings.EqualFold(meta.Type, "Flag") {
			return fmt.Errorf("%s: using INFO flag filtering on non flag type %s will not work correctly", flagName, tag)
		}
	}
	return nil
}

// passKeepINFOSite implements upstream's `--keep-INFO TAG` site filter
// (entry_filters.cpp:1033-1063). A site passes when at least one of the
// named INFO Flag tags is present in its INFO field. When `flags` is empty
// the filter is inactive (caller should not invoke).
//
// "Present" matches upstream's `get_INFO_value(tag) == "1"` test: the tag
// is in the INFO column, regardless of whether it appears as a bare flag
// (`MYFLAG`) or with a `=1` suffix.
func passKeepINFOSite(v *vcf.Variant, flags infoTagSet) bool {
	if len(flags) == 0 {
		return true
	}
	for tag := range flags {
		if infoFlagPresent(v, tag) {
			return true
		}
	}
	return false
}

// passRemoveINFOSite implements upstream's `--remove-INFO TAG` site filter
// (entry_filters.cpp:1068-1086). A site is DROPPED (returns false) when any
// of the named INFO Flag tags is present in its INFO field — i.e. the
// polarity-inverted complement of passKeepINFOSite. When `flags` is empty
// the filter is inactive (caller should not invoke).
//
// Composes with passKeepINFOSite per upstream's
// filter_sites_by_INFO ordering: keep narrows first, remove vetoes the
// survivors.
func passRemoveINFOSite(v *vcf.Variant, flags infoTagSet) bool {
	if len(flags) == 0 {
		return true
	}
	for tag := range flags {
		if infoFlagPresent(v, tag) {
			return false
		}
	}
	return true
}

// infoFlagPresent mirrors upstream's `get_INFO_value(tag) == "1"` presence
// test for Flag-type INFO tags. The parser stores bare-flag form
// (`MYFLAG;OTHER=...`) with an empty string value, so both `""` and `"1"`
// count as present; anything else (e.g. `MYFLAG=0`) is treated as absent.
func infoFlagPresent(v *vcf.Variant, tag string) bool {
	val, ok := v.Info[tag]
	if !ok {
		return false
	}
	return val == "" || val == "1"
}
