package bcftools

import (
	"fmt"
	"strconv"
	"strings"
)

// PloidyTable encodes the per-region, per-sex ploidy map used by
// `bcftools call --ploidy`. It mirrors the data structure built by
// reference_code/bcftools/ploidy.c on top of regidx. We keep it as a
// flat list of intervals because the GRCh37/GRCh38 predefs contain only
// a handful of entries and `call` does at most one lookup per record.
//
// Each entry covers [from, to] (1-based, inclusive) on a chromosome
// (possibly "*") for a single sex (possibly "*") and assigns a ploidy.
// The "*" sex/chrom rows act as catch-all defaults.
type PloidyTable struct {
	// sexes is the ordered list of registered sex names. The last entry
	// is the "default" sex assigned to samples without an explicit PED
	// sex (upstream: see vcfcall.c around lines 700-718).
	sexes []string
	// entries lists every per-region rule in declaration order.
	entries []ploidyEntry
	// dflt is the per-sex fallback ploidy used when no explicit entry
	// matches a (chrom, pos, sex) query AND no "*" rule covers that
	// sex. Indexed by sex id.
	dflt []int
	// builtinDflt is the global default ploidy supplied to ploidy_init
	// (upstream "2" for --ploidy {GRCh37,GRCh38}).
	builtinDflt int
}

// ploidyEntry is one row of the ploidy table.
type ploidyEntry struct {
	chrom  string // "*" means any chromosome
	from   int    // 1-based, inclusive; -1 means -infinity
	to     int    // 1-based, inclusive; -1 means +infinity
	sexID  int    // -1 means any sex
	ploidy int
}

// PredefPloidyAliases lists the alias names the CLI accepts and the
// strings ploidy_init_string parses. Mirrors vcfcall.c
// `ploidy_predefs`.
var PredefPloidyAliases = []struct {
	Alias string
	About string
	Body  string
}{
	{
		Alias: "GRCh37",
		About: "Human Genome reference assembly GRCh37 / hg19",
		Body: "X 1 60000 M 1\n" +
			"X 2699521 154931043 M 1\n" +
			"Y 1 59373566 M 1\n" +
			"Y 1 59373566 F 0\n" +
			"MT 1 16569 M 1\n" +
			"MT 1 16569 F 1\n" +
			"chrX 1 60000 M 1\n" +
			"chrX 2699521 154931043 M 1\n" +
			"chrY 1 59373566 M 1\n" +
			"chrY 1 59373566 F 0\n" +
			"chrM 1 16569 M 1\n" +
			"chrM 1 16569 F 1\n" +
			"*  * *     M 2\n" +
			"*  * *     F 2\n",
	},
	{
		Alias: "GRCh38",
		About: "Human Genome reference assembly GRCh38 / hg38",
		Body: "X 1 9999 M 1\n" +
			"X 2781480 155701381 M 1\n" +
			"Y 1 57227415 M 1\n" +
			"Y 1 57227415 F 0\n" +
			"MT 1 16569 M 1\n" +
			"MT 1 16569 F 1\n" +
			"chrX 1 9999 M 1\n" +
			"chrX 2781480 155701381 M 1\n" +
			"chrY 1 57227415 M 1\n" +
			"chrY 1 57227415 F 0\n" +
			"chrM 1 16569 M 1\n" +
			"chrM 1 16569 F 1\n" +
			"*  * *     M 2\n" +
			"*  * *     F 2\n",
	},
}

// LookupPredefPloidy returns the body string for an alias (case-insensitive)
// or "" if not found.
func LookupPredefPloidy(alias string) string {
	for _, p := range PredefPloidyAliases {
		if strings.EqualFold(p.Alias, alias) {
			return p.Body
		}
	}
	return ""
}

// ParsePloidyTable parses the body of a predefined ploidy spec (one rule
// per line: "CHROM FROM TO SEX PLOIDY"). It registers sex names in the
// order they appear, so the upstream behaviour (default sex = last
// registered, which for GRCh37/38 is F) is preserved.
func ParsePloidyTable(body string, dflt int) (*PloidyTable, error) {
	tbl := &PloidyTable{builtinDflt: dflt}
	sexIdx := func(name string) int {
		for i, s := range tbl.sexes {
			if s == name {
				return i
			}
		}
		tbl.sexes = append(tbl.sexes, name)
		tbl.dflt = append(tbl.dflt, dflt)
		return len(tbl.sexes) - 1
	}
	for ln, raw := range strings.Split(body, "\n") {
		line := strings.TrimSpace(raw)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) != 5 {
			return nil, fmt.Errorf("ploidy: line %d: expected 5 fields, got %d", ln+1, len(fields))
		}
		chrom := fields[0]
		fromS, toS, sexS, ploidyS := fields[1], fields[2], fields[3], fields[4]
		ploidy, err := strconv.Atoi(ploidyS)
		if err != nil {
			return nil, fmt.Errorf("ploidy: line %d: invalid ploidy %q", ln+1, ploidyS)
		}
		entry := ploidyEntry{chrom: chrom, ploidy: ploidy, sexID: -1}
		if sexS != "*" {
			entry.sexID = sexIdx(sexS)
		}
		if fromS == "*" {
			entry.from = -1
		} else {
			n, err := strconv.Atoi(fromS)
			if err != nil {
				return nil, fmt.Errorf("ploidy: line %d: invalid from %q", ln+1, fromS)
			}
			entry.from = n
		}
		if toS == "*" {
			entry.to = -1
		} else {
			n, err := strconv.Atoi(toS)
			if err != nil {
				return nil, fmt.Errorf("ploidy: line %d: invalid to %q", ln+1, toS)
			}
			entry.to = n
		}
		// A "* * *" row updates the per-sex default rather than a region
		// entry — this matches ploidy.c (it stores the "*" rows separately
		// in dflt_ploidy[]).
		if entry.chrom == "*" && entry.from == -1 && entry.to == -1 {
			if entry.sexID == -1 {
				for i := range tbl.dflt {
					tbl.dflt[i] = ploidy
				}
				tbl.builtinDflt = ploidy
			} else {
				tbl.dflt[entry.sexID] = ploidy
			}
			continue
		}
		tbl.entries = append(tbl.entries, entry)
	}
	if len(tbl.sexes) == 0 {
		// No sex was ever registered (e.g. body uses only "* * *" rows
		// with "*" sex). Register a single "*" sex so callers have an
		// entry to map samples to.
		tbl.sexes = []string{"*"}
		tbl.dflt = []int{tbl.builtinDflt}
	}
	return tbl, nil
}

// NSex returns the number of registered sexes.
func (t *PloidyTable) NSex() int {
	if t == nil {
		return 0
	}
	return len(t.sexes)
}

// SexID returns the numeric id for the named sex, or -1 if unknown.
func (t *PloidyTable) SexID(name string) int {
	if t == nil {
		return -1
	}
	for i, s := range t.sexes {
		if s == name {
			return i
		}
	}
	return -1
}

// DefaultSexID returns the id of the "default" sex (the last one
// registered). For the GRCh37/38 predefs this is "F", matching
// vcfcall.c's `sample2sex[i] = args->nsex - 1`.
func (t *PloidyTable) DefaultSexID() int {
	if t == nil || len(t.sexes) == 0 {
		return -1
	}
	return len(t.sexes) - 1
}

// SexName returns the name for a registered sex id, or "" if unknown.
func (t *PloidyTable) SexName(id int) string {
	if t == nil || id < 0 || id >= len(t.sexes) {
		return ""
	}
	return t.sexes[id]
}

// Query returns the ploidy for (chrom, pos, sexID). pos is 1-based.
// When no explicit entry matches, the per-sex default is returned.
func (t *PloidyTable) Query(chrom string, pos, sexID int) int {
	if t == nil {
		return 2
	}
	// Walk entries in declaration order, last match wins (matches
	// ploidy.c: it stores intervals in regidx which is order-preserving
	// and the query uses regitr_loop's last hit).
	best := -1
	for _, e := range t.entries {
		if e.chrom != "*" && e.chrom != chrom {
			continue
		}
		if e.from != -1 && pos < e.from {
			continue
		}
		if e.to != -1 && pos > e.to {
			continue
		}
		if e.sexID != -1 && e.sexID != sexID {
			continue
		}
		best = e.ploidy
	}
	if best >= 0 {
		return best
	}
	if sexID >= 0 && sexID < len(t.dflt) {
		return t.dflt[sexID]
	}
	return t.builtinDflt
}

// MaxPloidy returns the largest ploidy that appears in the table or
// per-sex defaults. Mirrors ploidy_max().
func (t *PloidyTable) MaxPloidy() int {
	if t == nil {
		return 2
	}
	m := t.builtinDflt
	for _, e := range t.entries {
		if e.ploidy > m {
			m = e.ploidy
		}
	}
	for _, d := range t.dflt {
		if d > m {
			m = d
		}
	}
	return m
}

// PerSamplePloidy returns the per-sample ploidy slice for one record.
// sexIDs[i] is the registered sex id of sample i; pass DefaultSexID for
// each unmapped sample (matches upstream's default).
func (t *PloidyTable) PerSamplePloidy(chrom string, pos int, sexIDs []int) []int {
	out := make([]int, len(sexIDs))
	if t == nil {
		for i := range out {
			out[i] = 2
		}
		return out
	}
	for i, sid := range sexIDs {
		out[i] = t.Query(chrom, pos, sid)
	}
	return out
}
