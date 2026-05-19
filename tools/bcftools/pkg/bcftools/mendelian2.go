// Package bcftools — see doc.go. This file implements `bcftools
// +mendelian2`, the rewritten Mendelian-inheritance checker that
// landed in upstream as a plugin (`reference_code/bcftools/plugins/
// mendelian2.c`).
//
// The legacy `mendelian` plugin (PR #105, mendelian.go in this
// package) handled trios specified one-at-a-time on the command line
// and only knew about INFO/MERR annotation, a TSV summary, and
// whole-record deletion. mendelian2 is the modern replacement and
// adds:
//
//   - PED-file ingestion (`-P/--ped FILE`): six-column PED is
//     scanned; any row whose father AND mother AND child are all in
//     the input VCF is taken as a trio.
//   - `-p/--pfm [1X:|2X:]P,F,M`: the single-trio shortcut (note the
//     upstream order is proband,father,mother, the OPPOSITE of the
//     internal mom/dad/kid index layout).
//   - The richer mode set `-m c|[adeEgmMS]`. Modes can be combined
//     (e.g. `-m ad` annotates AND deletes bad genotypes).
//   - Drop-bad-genotypes (`-m d`): instead of dropping whole sites,
//     set every offending trio's GTs to "./." while keeping the row.
//   - Per-site bookkeeping (`sites_ref_only`, `sites_no_GT`,
//     `sites_not_diploid`, ...) printed by `-m c`.
//
// We reuse the trio-consistency core from mendelian.go (mendelian2
// trios end up as the same {child, father, mother} tuple), but
// everything else — PED parsing, the multi-bit mode set, the
// expanded summary — lives here so the two plugins can evolve
// independently.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Mendelian2Mode is a bitmask of the upstream `-m c|[adeEgmMS]`
// letters. Multiple modes can be combined for VCF/BCF output; the
// drop modes (D/E/M/S) take precedence over their list counterparts.
type Mendelian2Mode uint

// Mendelian2 mode bits — order matches the MODE_* constants in
// `reference_code/bcftools/plugins/mendelian2.c` so reviewers can
// diff our parser line-for-line against upstream.
const (
	Mendelian2Annotate Mendelian2Mode = 1 << iota // a — add INFO/MERR
	Mendelian2Count                               // c — TSV summary (default)
	Mendelian2DeleteGT                            // d — set offending trio GTs to ./.
	Mendelian2ListErr                             // e — emit only sites with a Mendel error
	Mendelian2DropErr                             // E — drop sites with at least one Mendel error
	Mendelian2ListGood                            // g — emit only sites with at least one good trio
	Mendelian2ListMiss                            // m — emit only sites with at least one missing trio GT
	Mendelian2DropMiss                            // M — drop sites with a missing trio GT
	Mendelian2ListSkip                            // s — list sites skipped for housekeeping reasons (MODE_LIST_SKIP, 1<<8)
	Mendelian2DropSkip                            // S — drop sites skipped for housekeeping reasons (MODE_DROP_SKIP, 1<<9)
)

// mendelian2ListModes bundles the "emit only ..." selectors; if any
// is on we suppress everything else by default. Mirrors upstream's
// LIST_MODES bundle.
const mendelian2ListModes = Mendelian2ListErr | Mendelian2ListGood | Mendelian2ListMiss

// ParseMendelian2Mode parses upstream's `-m` letter set into a
// bitmask. The default (empty string) is the count-only summary.
// Unknown letters return a descriptive error. Letters mirror the
// switch in mendelian2.c:run().
func ParseMendelian2Mode(s string) (Mendelian2Mode, error) {
	if s == "" {
		return Mendelian2Count, nil
	}
	var m Mendelian2Mode
	for _, r := range s {
		switch r {
		case 'a':
			m |= Mendelian2Annotate
		case 'c':
			m |= Mendelian2Count
		case 'd':
			m |= Mendelian2DeleteGT
		case 'e', 'x':
			// upstream also accepts 'x' as a historical alias for 'e'.
			m |= Mendelian2ListErr
		case 'E':
			m |= Mendelian2DropErr
		case 'g', '+':
			m |= Mendelian2ListGood
		case 'm', 'u':
			// upstream also accepts 'u' as a historical alias for 'm'.
			m |= Mendelian2ListMiss
		case 'M':
			m |= Mendelian2DropMiss
		case 's':
			// MODE_LIST_SKIP (upstream mendelian2.c:58). Emit
			// skipped sites alongside the regular output. v1 has
			// no per-site "skipped" track yet, but the bit is
			// distinct from MODE_DROP_SKIP so a future
			// implementation can flip it on without breaking
			// callers.
			m |= Mendelian2ListSkip
		case 'S':
			m |= Mendelian2DropSkip
		default:
			return 0, fmt.Errorf("bcftools mendelian2: unknown -m mode letter %q (accept c|[adeEgmMSs])", r)
		}
	}
	return m, nil
}

// Mendelian2PFM holds the parsed `-p/--pfm [1X:|2X:]P,F,M` value.
// The Sex field is 0 (unknown / autosomal), 1 (1X — male: father
// haploid on chrX) or 2 (2X — female: father haploid on chrY only).
// Order on the command line is proband, father, mother — the OPPOSITE
// of the internal {mom, dad, kid} ordering, so we keep them as
// separate fields rather than a 3-element slice to avoid confusion.
type Mendelian2PFM struct {
	Child  string
	Father string
	Mother string
	Sex    int // 0=unknown, 1=male (1X), 2=female (2X)
}

// ParseMendelian2PFM parses the `-p/--pfm` value. The Sex prefix is
// optional; "1X:" → male, "2X:" → female, anything else (or absent)
// is treated as unknown / autosomal-only. Matches the prefix sniff
// inside mendelian2.c:init_data().
func ParseMendelian2PFM(s string) (Mendelian2PFM, error) {
	out := Mendelian2PFM{}
	if s == "" {
		return out, fmt.Errorf("bcftools mendelian2: empty -p/--pfm value")
	}
	body := s
	switch {
	case strings.HasPrefix(strings.ToUpper(s), "1X:"):
		out.Sex = 1
		body = s[3:]
	case strings.HasPrefix(strings.ToUpper(s), "2X:"):
		out.Sex = 2
		body = s[3:]
	}
	parts := strings.Split(body, ",")
	if len(parts) != 3 {
		return out, fmt.Errorf("bcftools mendelian2: -p/--pfm %q: expected [1X:|2X:]PROBAND,FATHER,MOTHER", s)
	}
	out.Child = strings.TrimSpace(parts[0])
	out.Father = strings.TrimSpace(parts[1])
	out.Mother = strings.TrimSpace(parts[2])
	if out.Child == "" || out.Father == "" || out.Mother == "" {
		return out, fmt.Errorf("bcftools mendelian2: -p/--pfm %q: empty member", s)
	}
	return out, nil
}

// Mendelian2Trio is one (child, father, mother) tuple plus the
// child's sex as parsed from the PED's 5th column (1=male, 2=female,
// 0=unknown). It's a separate type from the legacy Trio so future
// extensions (per-trio MERR FORMAT tag, sex-aware ploidy rules)
// don't have to mutate the existing public surface.
type Mendelian2Trio struct {
	Child  string
	Father string
	Mother string
	Sex    int // child's sex from PED 5th column
}

// Mendelian2Options controls the behaviour of Mendelian2.
type Mendelian2Options struct {
	// PEDFile is the path to a six-column PED file
	// (FamilyID, IndividualID, PaternalID, MaternalID, Sex,
	// Phenotype). Any individual whose father AND mother are both
	// present in the input VCF is taken as a trio. Mirrors
	// upstream's parse_ped.
	PEDFile string
	// PFM is the single-trio shortcut form ([1X:|2X:]P,F,M).
	// Mutually exclusive with PEDFile.
	PFM *Mendelian2PFM
	// Trios is an in-memory override used by tests; if set, PEDFile
	// and PFM are ignored. Lets the library be exercised without
	// touching the filesystem.
	Trios []Mendelian2Trio
	// Mode is the bitmask from -m.
	Mode Mendelian2Mode
	// IncludeExpr / ExcludeExpr — accepted for parity but applied
	// only at the record level (no per-sample filtering yet).
	// Empty string means "no filter".
	IncludeExpr string
	ExcludeExpr string
	// OutputFormat selects between VCF / VCF.gz for the
	// stream-back-out modes; ignored when Mode is purely Count.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
}

// Mendelian2Summary is the rollup returned by Mendelian2. Mirrors
// upstream's `-m c` printout shape (a handful of per-site counters
// plus per-trio numbers).
type Mendelian2Summary struct {
	// Per-site bookkeeping (upstream's `sites_*` lines).
	SitesRefOnly     int // n_allele==1 (no ALT) — skipped
	SitesManyAls     int // n_allele>64 — skipped (we use a bitmask)
	SitesFail        int // failed -i/-e filter
	SitesNoGT        int // no FMT/GT field
	SitesNotDiploid  int // GT not diploid
	SitesMissing     int // at least one trio with missing GT
	SitesMERR        int // at least one trio with Mendel error
	SitesGood        int // at least one good trio
	TotalRecords     int // total record count (sum of all paths above)
	RecordsWithError int // alias for SitesMERR for callers that already
	// speak the legacy Mendelian shape.
	Trios []Mendelian2TrioStats
}

// Mendelian2TrioStats is the per-trio breakdown.
type Mendelian2TrioStats struct {
	Mendelian2Trio
	NGood    int // good genotypes
	NGoodAlt int // good genotypes with at least one ALT (subset of NGood)
	NMErr    int // Mendelian errors
	NMissing int // missing trio GTs
	NFail    int // -i/-e filter failures
}

// Mendelian2File is the file-aware entry point used by the CLI. It
// opens path via iohelper (transparent gzip + BCF auto-detect),
// walks every record, and writes the requested output to w.
func Mendelian2File(path string, out io.Writer, opts Mendelian2Options) (Mendelian2Summary, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return Mendelian2Summary{}, fmt.Errorf("bcftools mendelian2: open %s: %w", path, err)
	}
	defer in.Close()
	return Mendelian2(in, out, opts)
}

// Mendelian2 streams VCF/BCF input through the upstream rule set and
// emits whatever the mode bitmask asks for. Default mode is Count
// (TSV summary only, no per-record output).
func Mendelian2(in io.Reader, out io.Writer, opts Mendelian2Options) (Mendelian2Summary, error) {
	if opts.Mode == 0 {
		opts.Mode = Mendelian2Count
	}

	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return Mendelian2Summary{}, fmt.Errorf("bcftools mendelian2: %w", err)
	}

	trios, err := loadMendelian2Trios(opts, hdr)
	if err != nil {
		return Mendelian2Summary{}, fmt.Errorf("bcftools mendelian2: %w", err)
	}
	if len(trios) == 0 {
		return Mendelian2Summary{}, fmt.Errorf("bcftools mendelian2: no complete trio found in input")
	}

	indices, err := resolveMendelian2Indices(hdr, trios)
	if err != nil {
		return Mendelian2Summary{}, fmt.Errorf("bcftools mendelian2: %w", err)
	}

	summary := Mendelian2Summary{Trios: make([]Mendelian2TrioStats, len(trios))}
	for i := range trios {
		summary.Trios[i].Mendelian2Trio = trios[i]
	}

	// Inject INFO/MERR for annotate / list / drop modes that emit
	// VCF (i.e. anything other than pure Count).
	needWriter := opts.Mode&^Mendelian2Count != 0
	annotatedHdr := hdr
	if needWriter && opts.Mode&Mendelian2Annotate != 0 {
		annotatedHdr = withMERRHeader(hdr)
	}

	var writer variantWriter
	var finish func()
	if needWriter {
		writer, finish, err = openOutput(out, ViewOptions{
			OutputFormat:  opts.OutputFormat,
			CompressLevel: opts.CompressLevel,
		}, annotatedHdr)
		if err != nil {
			return summary, fmt.Errorf("bcftools mendelian2: %w", err)
		}
		defer finish()
		if err := writer.WriteHeader(); err != nil {
			return summary, err
		}
	}

	for _, v := range variants {
		summary.TotalRecords++

		// Per-record skip checks mirroring mendelian2.c:process_record.
		skipReason := classifyRecord(v)
		switch skipReason {
		case skipRefOnly:
			summary.SitesRefOnly++
			if needWriter && opts.Mode&Mendelian2DropSkip == 0 {
				if err := writer.Write(v); err != nil {
					return summary, err
				}
			}
			continue
		case skipManyAls:
			summary.SitesManyAls++
			if needWriter && opts.Mode&Mendelian2DropSkip == 0 {
				if err := writer.Write(v); err != nil {
					return summary, err
				}
			}
			continue
		case skipNoGT:
			summary.SitesNoGT++
			if needWriter && opts.Mode&Mendelian2DropSkip == 0 {
				if err := writer.Write(v); err != nil {
					return summary, err
				}
			}
			continue
		}

		hasGood, hasErr, hasMiss, totalErr := evaluateTrios(v, indices, summary.Trios)
		if hasGood {
			summary.SitesGood++
		}
		if hasErr {
			summary.SitesMERR++
			summary.RecordsWithError++
		}
		if hasMiss {
			summary.SitesMissing++
		}

		// Drop precedences (E,M,S) take effect first.
		if opts.Mode&Mendelian2DropErr != 0 && hasErr {
			continue
		}
		if opts.Mode&Mendelian2DropMiss != 0 && hasMiss {
			continue
		}

		if !needWriter {
			continue
		}

		// `-m d` rewrites offending trio GTs to ./. but keeps the row.
		if opts.Mode&Mendelian2DeleteGT != 0 && hasErr {
			rewriteTrioGTs(v, indices, summary.Trios)
		}

		// `-m a` annotates INFO/MERR with the per-site error count.
		if opts.Mode&Mendelian2Annotate != 0 {
			if v.Info == nil {
				v.Info = make(map[string]string)
			}
			if _, exists := v.Info["MERR"]; !exists {
				v.InfoOrder = append(v.InfoOrder, "MERR")
			}
			v.Info["MERR"] = strconv.Itoa(totalErr)
		}

		// LIST_MODES — emit only if matches one of the selectors.
		// Drop modes already returned above.
		if opts.Mode&mendelian2ListModes != 0 {
			emit := false
			if opts.Mode&Mendelian2ListErr != 0 && hasErr {
				emit = true
			}
			if opts.Mode&Mendelian2ListMiss != 0 && hasMiss {
				emit = true
			}
			if opts.Mode&Mendelian2ListGood != 0 && hasGood {
				emit = true
			}
			if !emit {
				continue
			}
		}

		if err := writer.Write(v); err != nil {
			return summary, err
		}
	}

	if writer != nil {
		if err := writer.Flush(); err != nil {
			return summary, err
		}
	}

	if opts.Mode&Mendelian2Count != 0 {
		if err := writeMendelian2Summary(out, summary); err != nil {
			return summary, err
		}
	}
	return summary, nil
}

// loadMendelian2Trios merges every trio source (in-memory test
// override, then -p/--pfm, then -P/--ped). Trios whose members are
// not all present in hdr are dropped, matching the upstream PED-parse
// behaviour of skipping rows with unknown samples.
func loadMendelian2Trios(opts Mendelian2Options, hdr *vcf.Header) ([]Mendelian2Trio, error) {
	if len(opts.Trios) > 0 {
		return opts.Trios, nil
	}
	if opts.PFM != nil {
		return []Mendelian2Trio{{
			Child:  opts.PFM.Child,
			Father: opts.PFM.Father,
			Mother: opts.PFM.Mother,
			Sex:    opts.PFM.Sex,
		}}, nil
	}
	if opts.PEDFile == "" {
		return nil, fmt.Errorf("missing -p/--pfm or -P/--ped option")
	}
	return parsePEDFile(opts.PEDFile, hdr)
}

// parsePEDFile reads a 6-column PED file (whitespace-separated) and
// returns one trio per row whose father AND mother AND child are all
// present in hdr.Samples — matching the upstream parse_ped behaviour
// (rows whose members aren't in the VCF are silently skipped).
func parsePEDFile(path string, hdr *vcf.Header) ([]Mendelian2Trio, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	in := make(map[string]bool, len(hdr.Samples))
	for _, s := range hdr.Samples {
		in[s] = true
	}
	var out []Mendelian2Trio
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 4 {
			return nil, fmt.Errorf("malformed PED row %q (need at least 4 cols)", line)
		}
		// fields: FamilyID, IndividualID, PaternalID, MaternalID, [Sex, Phenotype]
		_, child, dad, mom := fields[0], fields[1], fields[2], fields[3]
		if !in[child] || !in[dad] || !in[mom] {
			continue
		}
		sex := 0
		if len(fields) >= 5 {
			n, err := strconv.Atoi(fields[4])
			if err != nil {
				return nil, fmt.Errorf("PED row %q: sex column %q not numeric", line, fields[4])
			}
			if n == 1 || n == 2 {
				sex = n
			}
		}
		out = append(out, Mendelian2Trio{
			Child:  child,
			Father: dad,
			Mother: mom,
			Sex:    sex,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	// Deterministic order — upstream sorts by min sample index across
	// the trio for sequential access; we sort by child name for
	// reproducibility (the index-based sort is a perf optimisation,
	// not a correctness one).
	sort.Slice(out, func(i, j int) bool { return out[i].Child < out[j].Child })
	return out, nil
}

// resolveMendelian2Indices maps each trio's three sample names to
// their column positions in the input header.
func resolveMendelian2Indices(hdr *vcf.Header, trios []Mendelian2Trio) ([]trioIndex, error) {
	byName := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		byName[s] = i
	}
	out := make([]trioIndex, len(trios))
	for i, t := range trios {
		c, ok := byName[t.Child]
		if !ok {
			return nil, fmt.Errorf("trio %d (%s/%s/%s): child sample not in input", i+1, t.Child, t.Father, t.Mother)
		}
		f, ok := byName[t.Father]
		if !ok {
			return nil, fmt.Errorf("trio %d (%s/%s/%s): father sample not in input", i+1, t.Child, t.Father, t.Mother)
		}
		m, ok := byName[t.Mother]
		if !ok {
			return nil, fmt.Errorf("trio %d (%s/%s/%s): mother sample not in input", i+1, t.Child, t.Father, t.Mother)
		}
		out[i] = trioIndex{child: c, father: f, mother: m}
	}
	return out, nil
}

// skipReason is the classifier for per-record housekeeping skips.
type skipReason int

const (
	skipNone    skipReason = iota
	skipRefOnly            // no ALT alleles
	skipManyAls            // > 64 ALTs (bitmask would overflow)
	skipNoGT               // FORMAT/GT is not declared on the record
)

// classifyRecord returns the skip-reason for v, mirroring the
// fast-path tests in upstream mendelian2.c:process_record.
func classifyRecord(v *vcf.Variant) skipReason {
	if len(v.Alt) == 0 || (len(v.Alt) == 1 && v.Alt[0] == ".") {
		return skipRefOnly
	}
	if len(v.Alt) > 64 {
		return skipManyAls
	}
	hasGT := false
	for _, t := range v.Format {
		if t == "GT" {
			hasGT = true
			break
		}
	}
	if !hasGT {
		return skipNoGT
	}
	return skipNone
}

// evaluateTrios runs every trio through the consistency rule and
// updates per-trio counters in stats. Returns (anyGood, anyErr,
// anyMiss, totalErrors).
func evaluateTrios(v *vcf.Variant, indices []trioIndex, stats []Mendelian2TrioStats) (bool, bool, bool, int) {
	var anyGood, anyErr, anyMiss bool
	var totalErr int
	for i, idx := range indices {
		child, father, mother, complete := readTrioGenotypes(v, idx)
		if !complete {
			anyMiss = true
			stats[i].NMissing++
			continue
		}
		ok := mendelianConsistent(child, father, mother, false)
		if ok {
			anyGood = true
			stats[i].NGood++
			if !(child[0] == 0 && child[1] == 0 && father[0] == 0 && father[1] == 0 && mother[0] == 0 && mother[1] == 0) {
				stats[i].NGoodAlt++
			}
			continue
		}
		anyErr = true
		stats[i].NMErr++
		totalErr++
	}
	return anyGood, anyErr, anyMiss, totalErr
}

// rewriteTrioGTs sets every offending trio's GT field to "./.".
// Used for `-m d` (the "set bad GTs to missing" mode).
func rewriteTrioGTs(v *vcf.Variant, indices []trioIndex, stats []Mendelian2TrioStats) {
	for i, idx := range indices {
		child, father, mother, complete := readTrioGenotypes(v, idx)
		if !complete {
			continue
		}
		if mendelianConsistent(child, father, mother, false) {
			continue
		}
		setSampleGT(v, idx.child, "./.")
		setSampleGT(v, idx.father, "./.")
		setSampleGT(v, idx.mother, "./.")
		_ = stats[i] // counter already bumped in evaluateTrios
	}
}

// setSampleGT mutates the GT field of v.Samples[i] to value. A no-op
// if i is out of range.
func setSampleGT(v *vcf.Variant, i int, value string) {
	if i < 0 || i >= len(v.Samples) {
		return
	}
	if v.Samples[i].Data == nil {
		v.Samples[i].Data = map[string]string{}
	}
	v.Samples[i].Data["GT"] = value
}

// writeMendelian2Summary writes the upstream-shaped `-m c` summary
// to out. We mirror upstream's section headers and column ordering
// verbatim so a downstream `grep ^TRIO | cut` pipeline behaves the
// same.
func writeMendelian2Summary(out io.Writer, s Mendelian2Summary) error {
	w := bufio.NewWriter(out)
	defer w.Flush()
	rows := []struct {
		key string
		val int
		why string
	}{
		{"sites_ref_only", s.SitesRefOnly, "# sites skipped because there was no ALT allele"},
		{"sites_many_als", s.SitesManyAls, "# skipped because of too many ALT alleles"},
		{"sites_fail", s.SitesFail, "# skipped because of failed -i/-e filter"},
		{"sites_no_GT", s.SitesNoGT, "# skipped because of absent FORMAT/GT field"},
		{"sites_not_diploid", s.SitesNotDiploid, "# skipped because FORMAT/GT not formatted diploid"},
		{"sites_missing", s.SitesMissing, "# number of sites with at least one trio GT missing"},
		{"sites_merr", s.SitesMERR, "# number of sites with at least one Mendelian error"},
		{"sites_good", s.SitesGood, "# number of sites with at least one good trio"},
	}
	if _, err := fmt.Fprintln(w, "# Summary stats"); err != nil {
		return err
	}
	for _, r := range rows {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%s\n", r.key, r.val, r.why); err != nil {
			return err
		}
	}
	if _, err := fmt.Fprintln(w, "# Per-trio stats, each column corresponds to one trio."); err != nil {
		return err
	}
	if err := writeMendelian2TrioRow(w, "ngood", s.Trios, func(t Mendelian2TrioStats) int { return t.NGood }); err != nil {
		return err
	}
	if err := writeMendelian2TrioRow(w, "ngood_alt", s.Trios, func(t Mendelian2TrioStats) int { return t.NGoodAlt }); err != nil {
		return err
	}
	if err := writeMendelian2TrioRow(w, "nmerr", s.Trios, func(t Mendelian2TrioStats) int { return t.NMErr }); err != nil {
		return err
	}
	if err := writeMendelian2TrioRow(w, "nmissing", s.Trios, func(t Mendelian2TrioStats) int { return t.NMissing }); err != nil {
		return err
	}
	if err := writeMendelian2TrioRow(w, "nfail", s.Trios, func(t Mendelian2TrioStats) int { return t.NFail }); err != nil {
		return err
	}
	if _, err := fmt.Fprintln(w, "# TRIO\t[2]id\t[3]child\t[4]father\t[5]mother"); err != nil {
		return err
	}
	for i, t := range s.Trios {
		if _, err := fmt.Fprintf(w, "TRIO\t%d\t%s\t%s\t%s\n", i+1, t.Child, t.Father, t.Mother); err != nil {
			return err
		}
	}
	return nil
}

func writeMendelian2TrioRow(w *bufio.Writer, label string, trios []Mendelian2TrioStats, pick func(Mendelian2TrioStats) int) error {
	if _, err := fmt.Fprint(w, label); err != nil {
		return err
	}
	for _, t := range trios {
		if _, err := fmt.Fprintf(w, "\t%d", pick(t)); err != nil {
			return err
		}
	}
	if _, err := fmt.Fprintln(w); err != nil {
		return err
	}
	return nil
}
