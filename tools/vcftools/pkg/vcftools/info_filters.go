package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
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

// filterRecodeInfo applies --keep-INFO / --remove-INFO to a recoded variant's
// INFO map. keepInfo (when non-empty) restricts INFO to the listed tags;
// removeInfo (when non-empty) strips the listed tags. If both are set,
// upstream's behaviour is "intersect both restrictions", which we mirror.
//
// Returns a fresh map; the caller is responsible for assigning it back.
func filterRecodeInfo(info map[string]string, keepInfo, removeInfo infoTagSet) map[string]string {
	out := make(map[string]string, len(info))
	for k, v := range info {
		if len(keepInfo) > 0 {
			if _, hit := keepInfo[k]; !hit {
				continue
			}
		}
		if len(removeInfo) > 0 {
			if _, hit := removeInfo[k]; hit {
				continue
			}
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
