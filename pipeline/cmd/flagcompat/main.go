// Command flagcompat measures, per ported tool, the flag-compatibility
// percentage: of an upstream tool's documented CLI flags, how many our Go port
// accepts under the same name/spelling.
//
// This is a SURFACE-acceptance metric. It checks whether our port registers a
// flag with the same spelling as upstream — not whether the flag behaves
// identically. Semantic equivalence is covered separately by the byte-exact
// upstream-parity tests (the *Upstream* tests run in CI). See
// docs/manuscript/results/flag_compat.md for the rendered results and caveats.
//
// Method:
//
//   - The UPSTREAM flag set per tool (or tool/subcommand) is a curated list in
//     upstream_flags.txt, embedded at build time. Each entry was derived from a
//     primary source under reference_code/ (getopt optstrings, longopts arrays,
//     GetOptions lists, or documented usage text). Provenance is documented in
//     the header of that file.
//   - OUR flag set per tool is extracted automatically from the port's Go
//     source via go/ast (see extract.go), scanning every non-test .go file
//     under tools/<tool>/ for pkg/cliflag and stdlib flag registrations. This
//     needs no execution and is robust to formatting.
//   - compat% = |upstream ∩ ours| / |upstream|, matching on the bare flag token
//     (a name with its leading dashes removed). For subcommanded tools we pool
//     our port's flags across the whole tool family and test each subcommand's
//     upstream flags against that pool, i.e. "does the port accept this
//     spelling anywhere in the tool".
//
// Usage:
//
//	go run ./pipeline/cmd/flagcompat            # print the markdown report
//	go run ./pipeline/cmd/flagcompat -tools DIR # override the tools/ directory
package main

import (
	_ "embed"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

//go:embed upstream_flags.txt
var upstreamData string

// entry is one upstream flag record: a tool (or tool/subcommand) and the set of
// upstream flag tokens it documents.
type entry struct {
	tool  string   // family name, e.g. "samtools"
	sub   string   // subcommand, e.g. "view"; empty for single-command tools
	flags []string // distinct upstream flag tokens (bare, no dashes)
}

func (e entry) label() string {
	if e.sub == "" {
		return e.tool
	}
	return e.tool + " " + e.sub
}

// result holds the computed compatibility for one entry.
type result struct {
	entry
	matched int      // |upstream ∩ ours|
	total   int      // |upstream|
	missing []string // upstream flags our port does not register
}

func (r result) pct() float64 {
	if r.total == 0 {
		return 100.0
	}
	return 100.0 * float64(r.matched) / float64(r.total)
}

func main() {
	toolsDir := flag.String("tools", "", "path to the tools/ directory (default: auto-detected from the repo root)")
	flag.Parse()

	dir, err := resolveToolsDir(*toolsDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "flagcompat:", err)
		os.Exit(1)
	}

	entries, err := parseUpstream(upstreamData)
	if err != nil {
		fmt.Fprintln(os.Stderr, "flagcompat:", err)
		os.Exit(1)
	}

	// Extract our flags once per tool family.
	ourFlags := map[string]map[string]struct{}{}
	for _, e := range entries {
		if _, done := ourFlags[e.tool]; done {
			continue
		}
		fs, err := extractOurFlags(filepath.Join(dir, e.tool))
		if err != nil {
			fmt.Fprintf(os.Stderr, "flagcompat: extracting %s: %v\n", e.tool, err)
			os.Exit(1)
		}
		ourFlags[e.tool] = fs
	}

	results := make([]result, 0, len(entries))
	for _, e := range entries {
		ours := ourFlags[e.tool]
		r := result{entry: e, total: len(e.flags)}
		for _, f := range e.flags {
			if _, ok := ours[f]; ok {
				r.matched++
			} else {
				r.missing = append(r.missing, f)
			}
		}
		results = append(results, r)
	}

	fmt.Print(render(results))
}

// resolveToolsDir returns the tools/ directory, either from the override or by
// walking up from the working directory to the module root (where go.mod is).
func resolveToolsDir(override string) (string, error) {
	if override != "" {
		return override, nil
	}
	wd, err := os.Getwd()
	if err != nil {
		return "", err
	}
	for dir := wd; ; {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return filepath.Join(dir, "tools"), nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", fmt.Errorf("could not locate go.mod above %s", wd)
		}
		dir = parent
	}
}

// parseUpstream parses the embedded upstream_flags.txt. Lines beginning with #
// are comments; a trailing backslash continues a record onto the next line. The
// first whitespace token is the tool[/subcommand] key; the rest are flag
// tokens.
func parseUpstream(data string) ([]entry, error) {
	var entries []entry
	var logical []string
	// Stitch backslash-continued lines into logical records first.
	var buf strings.Builder
	for _, raw := range strings.Split(data, "\n") {
		line := raw
		if i := strings.IndexByte(line, '#'); i >= 0 {
			line = line[:i]
		}
		line = strings.TrimRight(line, " \t\r")
		if strings.HasSuffix(line, "\\") {
			buf.WriteString(strings.TrimSuffix(line, "\\"))
			buf.WriteByte(' ')
			continue
		}
		buf.WriteString(line)
		logical = append(logical, buf.String())
		buf.Reset()
	}

	for _, rec := range logical {
		fields := strings.Fields(rec)
		if len(fields) == 0 {
			continue
		}
		key := fields[0]
		tool, sub, _ := strings.Cut(key, "/")
		seen := map[string]struct{}{}
		var flags []string
		for _, f := range fields[1:] {
			if _, dup := seen[f]; dup {
				continue
			}
			seen[f] = struct{}{}
			flags = append(flags, f)
		}
		entries = append(entries, entry{tool: tool, sub: sub, flags: flags})
	}
	return entries, nil
}

// render produces the markdown report: a per-entry table, a per-family rollup,
// and a flag-count-weighted aggregate.
func render(results []result) string {
	var b strings.Builder

	// Per-entry table.
	b.WriteString("| Tool / subcommand | Upstream flags | Accepted | Compat % | Missing (upstream flags not accepted) |\n")
	b.WriteString("|---|---:|---:|---:|---|\n")
	for _, r := range results {
		miss := "—"
		if len(r.missing) > 0 {
			sort.Strings(r.missing)
			quoted := make([]string, len(r.missing))
			for i, m := range r.missing {
				quoted[i] = "`" + m + "`"
			}
			miss = strings.Join(quoted, " ")
		}
		fmt.Fprintf(&b, "| %s | %d | %d | %.1f%% | %s |\n",
			r.label(), r.total, r.matched, r.pct(), miss)
	}

	// Per-family rollup (sum of matched / sum of total across the family's
	// entries).
	type fam struct {
		matched, total int
	}
	families := map[string]*fam{}
	var order []string
	for _, r := range results {
		f, ok := families[r.tool]
		if !ok {
			f = &fam{}
			families[r.tool] = f
			order = append(order, r.tool)
		}
		f.matched += r.matched
		f.total += r.total
	}
	b.WriteString("\n## Per-tool-family rollup\n\n")
	b.WriteString("| Tool family | Upstream flag slots | Accepted | Compat % |\n")
	b.WriteString("|---|---:|---:|---:|\n")
	var grandMatched, grandTotal int
	for _, name := range order {
		f := families[name]
		pct := 100.0
		if f.total > 0 {
			pct = 100.0 * float64(f.matched) / float64(f.total)
		}
		fmt.Fprintf(&b, "| %s | %d | %d | %.1f%% |\n", name, f.total, f.matched, pct)
		grandMatched += f.matched
		grandTotal += f.total
	}

	// Aggregate, weighted by upstream flag count.
	aggPct := 100.0
	if grandTotal > 0 {
		aggPct = 100.0 * float64(grandMatched) / float64(grandTotal)
	}
	b.WriteString("\n## Aggregate (weighted by upstream flag count)\n\n")
	fmt.Fprintf(&b, "- Upstream flag slots measured: **%d**\n", grandTotal)
	fmt.Fprintf(&b, "- Accepted under the same spelling: **%d**\n", grandMatched)
	fmt.Fprintf(&b, "- Weighted flag-compatibility: **%.1f%%**\n", aggPct)

	return b.String()
}
