package matrix

import (
	"sort"
	"strings"
)

// Flag describes one flag of a subcommand and the values worth exercising. A
// boolean flag (Values empty) contributes a single "+flag" variant; a
// value-bearing flag contributes one variant per value.
type Flag struct {
	// Name is the flag token as passed on the command line, e.g. "-q" or "-f".
	Name string

	// Values are the value choices to try. Empty means the flag is boolean and
	// is emitted with no argument.
	Values []string

	// Bool reports a no-argument flag even if Values is nil; kept explicit so a
	// value-bearing flag with a deliberately empty value set is distinguishable.
	Bool bool
}

// Combo is a curated multi-flag interaction: a named set of flags (each with a
// chosen value) that exercise a meaningful combination. Combos are how callers
// add specific interactions on top of the auto-generated single-flag cases.
type Combo struct {
	// Name labels the interaction in reports (e.g. "filter+count").
	Name string
	// Flags is the ordered flag tokens+values for this combo, e.g.
	// {"-f", "0x2", "-c"}.
	Flags []string
}

// ExpandSpec describes how to expand one subcommand into a curated set of
// entries.
//
// Combinatorics policy (IMPORTANT — read before extending):
//
// We deliberately do NOT generate the 2^N power set of flags. For a subcommand
// with N flags that is intractable (samtools view alone has dozens), produces
// mostly meaningless or mutually-exclusive combinations, and would make the
// report unreadable. Instead Expand produces:
//
//  1. A baseline entry (no extra flags) — the "does it run at all" case.
//  2. One entry per single flag value (every flag exercised in isolation).
//  3. Exactly the multi-flag interactions listed in Combos — hand-curated by
//     the matrix author to cover flags that genuinely interact (e.g. a region
//     filter combined with a count, or a format flag combined with a header
//     flag).
//
// This keeps the matrix size linear-plus-curated (N + |Combos|) rather than
// exponential, while still making it trivial to add a specific interaction:
// append a Combo. If a future tool truly needs a small full cross-product over
// a handful of flags, build it explicitly as Combos (a helper, CrossProduct,
// is provided for that bounded case).
type ExpandSpec struct {
	Tool           string
	Subcommand     string
	UpstreamTool   string
	UsesSubcommand bool
	Input          InputKind
	Compare        CompareMode
	Heavy          bool

	// BaseArgs are prepended to every generated entry (e.g. the input
	// placeholder and any always-on positional like a region file).
	BaseArgs []string

	// Flags are exercised one at a time (single-flag cases).
	Flags []Flag

	// Combos are the curated multi-flag interactions.
	Combos []Combo
}

// Expand turns an ExpandSpec into the curated set of entries described above.
func (s ExpandSpec) Expand() []Entry {
	mk := func(name string, extra []string) Entry {
		args := make([]string, 0, len(s.BaseArgs)+len(extra))
		args = append(args, s.BaseArgs...)
		args = append(args, extra...)
		return Entry{
			Tool:           s.Tool,
			Subcommand:     s.Subcommand,
			UpstreamTool:   s.UpstreamTool,
			UsesSubcommand: s.UsesSubcommand,
			Name:           name,
			Args:           args,
			Input:          s.Input,
			Compare:        s.Compare,
			Heavy:          s.Heavy,
		}
	}

	var out []Entry
	// 1. Baseline.
	out = append(out, mk("base", nil))

	// 2. Single-flag cases.
	for _, f := range s.Flags {
		if f.Bool || len(f.Values) == 0 {
			out = append(out, mk("flag"+sanitize(f.Name), []string{f.Name}))
			continue
		}
		for _, v := range f.Values {
			out = append(out, mk("flag"+sanitize(f.Name)+"_"+sanitize(v), []string{f.Name, v}))
		}
	}

	// 3. Curated multi-flag interactions.
	for _, c := range s.Combos {
		out = append(out, mk("combo_"+sanitize(c.Name), c.Flags))
	}
	return out
}

// CrossProduct builds the full cross-product of the given value-bearing flags
// as a slice of Combos. It exists for the BOUNDED case where a matrix author
// has deliberately chosen a small handful of flags whose every combination is
// meaningful — it is NOT used by Expand automatically, precisely to keep the
// no-power-set policy explicit. Callers pass the result as ExpandSpec.Combos.
func CrossProduct(flags ...Flag) []Combo {
	combos := []Combo{{Name: "x", Flags: nil}}
	for _, f := range flags {
		var next []Combo
		vals := f.Values
		if len(vals) == 0 {
			vals = []string{""} // boolean: present/absent handled by the empty combo
		}
		for _, c := range combos {
			for _, v := range vals {
				nf := append([]string{}, c.Flags...)
				if v == "" {
					nf = append(nf, f.Name)
				} else {
					nf = append(nf, f.Name, v)
				}
				next = append(next, Combo{
					Name:  strings.TrimPrefix(c.Name+"_"+sanitize(f.Name)+sanitize(v), "x_"),
					Flags: nf,
				})
			}
		}
		combos = next
	}
	// Stable order for reproducible reports.
	sort.Slice(combos, func(i, j int) bool { return combos[i].Name < combos[j].Name })
	return combos
}

// sanitize makes a flag/value safe and compact for use in an entry name.
func sanitize(s string) string {
	s = strings.TrimLeft(s, "-")
	r := strings.NewReplacer("/", "_", " ", "_", ".", "_", "{", "", "}", "", "0x", "x")
	return r.Replace(s)
}
