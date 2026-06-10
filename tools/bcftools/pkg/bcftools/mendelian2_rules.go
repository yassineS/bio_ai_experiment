// Package bcftools — see doc.go. This file implements the per-contig
// ploidy / inheritance rule machinery for `bcftools +mendelian2`,
// mirroring the `--rules ASSEMBLY` / `--rules-file FILE` support in
// upstream `reference_code/bcftools/plugins/mendelian2.c` (the
// `rules_predef_t` table, `parse_rules`, and `init_rules`).
//
// # Why this exists
//
// The v1 mendelian2 port only understood a single chrX heuristic (male
// children inherit their X from the mother). Upstream is much richer:
// it loads a table of `SEX_ID CHROM:BEG-END INHERITED_FROM` lines —
// either one of the built-in GRCh37 / GRCh38 assemblies or a custom
// `--rules-file` — and, for every site, derives a per-sex (1X / 2X)
// ploidy and inheritance mode. The pseudo-autosomal regions (PAR) on
// X stay diploid; the male-specific X stretch is haploid maternal; Y
// is haploid paternal for males and absent for females; MT is haploid
// maternal for everyone. This file reproduces that table-driven model
// exactly so the `--rules` flag is no longer a deferred stub.
//
// The rules grammar (one line per region):
//
//		SEX_ID  CHROM:BEG-END  INHERITED_FROM
//
//	  - SEX_ID is a free-form token; the human assemblies use "1X"
//	    (male child) and "2X" (female child). The set of distinct
//	    SEX_IDs seen across all lines defines the sex dictionary, in
//	    first-seen order (matching upstream's str2sex_id khash).
//	  - CHROM:BEG-END are 1-based inclusive coordinates.
//	  - INHERITED_FROM is a combination of M (allele inherited from the
//	    mother), F (from the father), "MF"/"FM" (both — diploid), or
//	    "." (neither — ploidy 0, e.g. Y in a female). The number of
//	    M/F letters is the ploidy (max 2).
//
// Any region not covered by a rule defaults to diploid, biparental
// (MF) inheritance for every sex — exactly upstream's behaviour for
// the autosomes.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Inheritance bits. These mirror the `1<<iMOM` / `1<<iDAD` masks in
// mendelian2.c (iMOM and iDAD are the mother/father slots of the trio
// index array). They are bit positions in MendelianRule.Inherits.
const (
	inheritMother = 1 << 0 // M — allele inherited from the mother
	inheritFather = 1 << 1 // F — allele inherited from the father
)

// MendelianRule is the resolved per-region, per-sex inheritance rule:
// how many alleles the child carries here (Ploidy, 0/1/2) and which
// parents may contribute them (Inherits, a combination of inheritMother
// and inheritFather). It mirrors upstream's `rule_t`.
type MendelianRule struct {
	// SexID is the dictionary index of the SEX_ID token this rule
	// applies to (e.g. 0 for the first SEX_ID seen, typically "1X").
	SexID int
	// Ploidy is the number of alleles the child carries in this
	// region for this sex: 0 (absent), 1 (haploid) or 2 (diploid).
	Ploidy int
	// Inherits is the bitmask of permitted parental contributions
	// (inheritMother | inheritFather). Zero means "no applicable
	// rule" (the child carries nothing here).
	Inherits int
}

// mendelianRuleRegion is one parsed `SEX_ID CHROM:BEG-END FROM` line.
// Beg/End are 0-based inclusive (converted from the 1-based inclusive
// coordinates in the file), matching upstream's regidx storage.
type mendelianRuleRegion struct {
	chrom string
	beg   int // 0-based inclusive
	end   int // 0-based inclusive
	rule  MendelianRule
}

// MendelianRules is the table-driven inheritance model loaded from a
// predefined assembly or a custom rules file. It maps every site to a
// per-sex MendelianRule, defaulting to diploid biparental inheritance
// outside any listed region. It is the Go analogue of upstream's
// `regidx_t *rules` plus the `str2sex_id` dictionary.
type MendelianRules struct {
	regions []mendelianRuleRegion
	// sexID maps a SEX_ID token (e.g. "1X") to its dictionary index,
	// assigned in first-seen order. PED/PFM sex values are translated
	// to these tokens ("1X" for male, "2X" for female) and looked up
	// here.
	sexID map[string]int
	// nSex is the number of distinct SEX_ID tokens (len(sexID)).
	nSex int
}

// NumSexes reports how many distinct SEX_ID tokens the rules define
// (two — "1X" and "2X" — for the human assemblies).
func (r *MendelianRules) NumSexes() int { return r.nSex }

// SexIDFor maps an integer sex code (1=male, 2=female, anything else
// treated as female, matching upstream's PED default) to the rule
// dictionary's SEX_ID index. It returns an error when the required
// SEX_ID token ("1X" or "2X") is not present in the rules — exactly
// upstream's `error("Missing the sex ...")`.
func (r *MendelianRules) SexIDFor(sex int) (int, error) {
	token := "2X"
	if sex == 1 {
		token = "1X"
	}
	id, ok := r.sexID[token]
	if !ok {
		return 0, fmt.Errorf("missing the sex %q, it's not in the rules", token)
	}
	return id, nil
}

// rulesFor returns, for the given site (chrom and 0-based inclusive
// [beg,end] span) and sex dictionary index, the applicable
// MendelianRule. The default — when no listed region overlaps — is
// diploid biparental inheritance, mirroring the `rule[i]` initialisation
// at the top of upstream's collect_stats.
func (r *MendelianRules) rulesFor(chrom string, beg, end, sexID int) MendelianRule {
	rule := MendelianRule{SexID: sexID, Ploidy: 2, Inherits: inheritMother | inheritFather}
	for i := range r.regions {
		reg := &r.regions[i]
		if reg.rule.SexID != sexID {
			continue
		}
		if reg.chrom != chrom {
			continue
		}
		if reg.end < beg || reg.beg > end {
			continue
		}
		// Last overlapping rule wins, matching upstream's
		// `while (regitr_overlap(itr)) rule[sex_id] = *rule;`.
		rule = reg.rule
		rule.SexID = sexID
	}
	return rule
}

// mendelian2PredefinedRules holds the built-in assembly tables. The
// rule text is copied verbatim from the `rules_predefs[]` array in
// mendelian2.c so reviewers can diff the coordinates line-for-line.
var mendelian2PredefinedRules = []struct {
	alias string
	about string
	rules string
}{
	{
		alias: "GRCh37",
		about: "Human Genome reference assembly GRCh37 / hg19, both chr naming conventions",
		rules: `   # Unlisted regions inherit from both parents (MF)
   1X  X:1-60000               M
   1X  X:2699521-154931043     M
   1X  Y:1-59373566            F
   2X  Y:1-59373566            .
   1X  MT:1-16569              M
   2X  MT:1-16569              M

   1X  chrX:1-60000            M
   1X  chrX:2699521-154931043  M
   1X  chrY:1-59373566         F
   2X  chrY:1-59373566         .
   1X  chrM:1-16569            M
   2X  chrM:1-16569            M
`,
	},
	{
		alias: "GRCh38",
		about: "Human Genome reference assembly GRCh38 / hg38, both chr naming conventions",
		rules: `   # Unlisted regions inherit from both parents (MF)
   1X  X:1-9999                M
   1X  X:2781480-155701381     M
   1X  Y:1-57227415            F
   2X  Y:1-57227415            .
   1X  MT:1-16569              M
   2X  MT:1-16569              M

   1X  chrX:1-9999             M
   1X  chrX:2781480-155701381  M
   1X  chrY:1-57227415         F
   2X  chrY:1-57227415         .
   1X  chrMT:1-16569           M
   2X  chrMT:1-16569           M
`,
	},
}

// MendelianRulesList renders the human-readable catalogue printed by
// upstream's `--rules list` (and, when detailed is true, the full
// per-region tables printed by `--rules list?`). It is returned as a
// string rather than written to stderr so the caller owns the I/O.
func MendelianRulesList(detailed bool) string {
	var b strings.Builder
	b.WriteString("\nPRE-DEFINED INHERITANCE RULES\n\n")
	b.WriteString(" * Columns are: SEX_ID CHROM:BEG-END INHERITED_FROM\n")
	b.WriteString(" * Coordinates are 1-based inclusive.\n\n")
	for _, p := range mendelian2PredefinedRules {
		fmt.Fprintf(&b, "%s\n   .. %s\n\n", p.alias, p.about)
		if detailed {
			fmt.Fprintf(&b, "%s\n", p.rules)
		}
	}
	b.WriteString("Run as --rules <alias> (e.g. --rules GRCh37).\n")
	b.WriteString("To see the detailed ploidy definition, append a question mark (e.g. --rules GRCh37?).\n\n")
	return b.String()
}

// ErrMendelianRulesList is returned by LoadMendelianRulesByName when
// the caller asked for a rules listing rather than a usable table. The
// CLI prints Listing and exits, mirroring upstream init_rules' two
// `exit(-1)` paths:
//
//   - `--rules list` / `list?` and any unknown alias print the whole
//     catalogue (Listing == MendelianRulesList output).
//   - `--rules GRCh38?` (a known alias with a trailing "?") prints just
//     that assembly's detailed per-region table.
type ErrMendelianRulesList struct {
	// Detailed is true for the "list?" form (full per-region tables in
	// the catalogue).
	Detailed bool
	// Listing is the exact text the CLI should print before exiting.
	Listing string
}

func (e *ErrMendelianRulesList) Error() string { return "mendelian2: rules listing requested" }

// LoadMendelianRulesByName resolves a predefined assembly alias
// (case-insensitive, e.g. "GRCh38") into a MendelianRules table. A
// trailing "?" on a concrete alias requests that assembly's detailed
// table; the aliases "list"/"list?" (and any unknown alias) request the
// full catalogue. For any listing request it returns an
// *ErrMendelianRulesList so the caller can print and exit, mirroring
// upstream init_rules.
func LoadMendelianRulesByName(alias string) (*MendelianRules, error) {
	if alias == "" {
		alias = "GRCh37"
	}
	detailed := false
	if strings.HasSuffix(alias, "?") {
		detailed = true
		alias = alias[:len(alias)-1]
	}
	if strings.EqualFold(alias, "list") {
		return nil, &ErrMendelianRulesList{Detailed: detailed, Listing: MendelianRulesList(detailed)}
	}
	for _, p := range mendelian2PredefinedRules {
		if strings.EqualFold(alias, p.alias) {
			if detailed {
				// `--rules GRCh38?`: upstream prints only this
				// assembly's table and exits.
				return nil, &ErrMendelianRulesList{Detailed: true, Listing: p.rules}
			}
			return parseMendelianRules(strings.NewReader(p.rules))
		}
	}
	// Unknown alias: upstream prints the catalogue and exits.
	return nil, &ErrMendelianRulesList{Detailed: detailed, Listing: MendelianRulesList(detailed)}
}

// LoadMendelianRulesFile parses a custom ploidy / inheritance rules
// file (the `--rules-file` format) into a MendelianRules table. The
// grammar matches LoadMendelianRulesByName's built-in tables.
func LoadMendelianRulesFile(path string) (*MendelianRules, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("could not read rules file: %w", err)
	}
	defer f.Close()
	return parseMendelianRules(f)
}

// parseMendelianRules parses the rules grammar from r. It assigns
// SEX_ID dictionary indices in first-seen order (matching upstream's
// str2sex_id), converts the 1-based inclusive coordinates to 0-based
// inclusive, and decodes the M/F/MF/. inheritance token.
func parseMendelianRules(r io.Reader) (*MendelianRules, error) {
	out := &MendelianRules{sexID: map[string]int{}}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		reg, err := parseMendelianRuleLine(line, out)
		if err != nil {
			return nil, err
		}
		out.regions = append(out.regions, reg)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	out.nSex = len(out.sexID)
	// Deterministic region order for reproducible overlap resolution;
	// the "last overlapping rule wins" precedence is preserved because
	// the predefined tables never list two overlapping regions for the
	// same (sex, chrom).
	sort.SliceStable(out.regions, func(i, j int) bool {
		a, b := out.regions[i], out.regions[j]
		if a.chrom != b.chrom {
			return a.chrom < b.chrom
		}
		if a.rule.SexID != b.rule.SexID {
			return a.rule.SexID < b.rule.SexID
		}
		return a.beg < b.beg
	})
	return out, nil
}

// parseMendelianRuleLine parses one `SEX_ID CHROM:BEG-END FROM` line,
// registering the SEX_ID token in out.sexID. It mirrors upstream's
// parse_rules field-by-field.
func parseMendelianRuleLine(line string, out *MendelianRules) (mendelianRuleRegion, error) {
	fields := strings.Fields(line)
	if len(fields) < 3 {
		return mendelianRuleRegion{}, fmt.Errorf("could not parse the rules line: %q", line)
	}
	sexTok := fields[0]
	sexID, ok := out.sexID[sexTok]
	if !ok {
		sexID = len(out.sexID)
		out.sexID[sexTok] = sexID
	}

	chrom, beg, end, err := parseRuleRegion(fields[1])
	if err != nil {
		return mendelianRuleRegion{}, fmt.Errorf("could not parse the region %q: %w", line, err)
	}

	inherits, ploidy, err := parseInheritanceToken(fields[2])
	if err != nil {
		return mendelianRuleRegion{}, fmt.Errorf("could not parse the region %q: %w", line, err)
	}

	return mendelianRuleRegion{
		chrom: chrom,
		beg:   beg,
		end:   end,
		rule:  MendelianRule{SexID: sexID, Ploidy: ploidy, Inherits: inherits},
	}, nil
}

// parseRuleRegion splits "CHROM:BEG-END" into the chromosome name and a
// 0-based inclusive [beg,end] span (the 1-based file coordinates minus
// one), matching upstream's `*beg = strtol(...)-1` / `*end = ...-1`.
func parseRuleRegion(s string) (string, int, int, error) {
	colon := strings.LastIndex(s, ":")
	if colon < 0 {
		return "", 0, 0, fmt.Errorf("missing ':' in region %q", s)
	}
	chrom := s[:colon]
	span := s[colon+1:]
	dash := strings.Index(span, "-")
	if dash < 0 {
		return "", 0, 0, fmt.Errorf("missing '-' in region %q", s)
	}
	beg1, err := strconv.Atoi(span[:dash])
	if err != nil {
		return "", 0, 0, fmt.Errorf("bad BEG in region %q", s)
	}
	end1, err := strconv.Atoi(span[dash+1:])
	if err != nil {
		return "", 0, 0, fmt.Errorf("bad END in region %q", s)
	}
	if chrom == "" {
		return "", 0, 0, fmt.Errorf("empty chromosome in region %q", s)
	}
	return chrom, beg1 - 1, end1 - 1, nil
}

// parseInheritanceToken decodes the INHERITED_FROM column (M / F / MF /
// FM / "."). The number of M/F letters is the ploidy; "." resets both
// to zero. It mirrors the per-character loop in upstream parse_rules.
func parseInheritanceToken(s string) (inherits, ploidy int, err error) {
	for _, c := range s {
		switch c {
		case 'M':
			inherits |= inheritMother
			ploidy++
		case 'F':
			inherits |= inheritFather
			ploidy++
		case '.':
			inherits = 0
			ploidy = 0
		default:
			return 0, 0, fmt.Errorf("unexpected inheritance character %q", string(c))
		}
		if ploidy > 2 {
			return 0, 0, fmt.Errorf("ploidy > 2 is not supported: %q", s)
		}
	}
	return inherits, ploidy, nil
}

// defaultMendelianRules returns the GRCh37 table, matching upstream's
// `init_rules(args, NULL)` default when neither --rules nor
// --rules-file is given.
func defaultMendelianRules() *MendelianRules {
	r, _ := parseMendelianRules(strings.NewReader(mendelian2PredefinedRules[0].rules))
	return r
}
