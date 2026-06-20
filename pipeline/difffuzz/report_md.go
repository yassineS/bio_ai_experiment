package difffuzz

import (
	"fmt"
	"os"
	"sort"
	"strings"
)

// WriteMarkdown renders the report as a human-readable Markdown document.
func (r *Report) WriteMarkdown(path string) error {
	var b strings.Builder
	fmt.Fprintf(&b, "# Differential fuzzing report\n\n")
	fmt.Fprintf(&b, "Seed: `%d`  Iterations/target: `%d`\n\n", r.Seed, r.Iterations)

	// Summary table.
	fmt.Fprintf(&b, "## Summary\n\n")
	fmt.Fprintf(&b, "| Target | Tool | Status | Inputs | Divergences | Coverage |\n")
	fmt.Fprintf(&b, "|---|---|---|---:|---|---|\n")
	for _, t := range r.Targets {
		status := "ran"
		if t.Skipped {
			status = "SKIP"
		}
		cov := "-"
		if t.Coverage.Enabled {
			cov = fmt.Sprintf("%.1f%%", t.Coverage.Percent)
		}
		fmt.Fprintf(&b, "| %s | %s | %s | %d | %s | %s |\n",
			t.Name, t.Tool, status, t.Inputs, classSummary(t.ByClass), cov)
	}
	b.WriteString("\n")

	// Per-target detail.
	for _, t := range r.Targets {
		fmt.Fprintf(&b, "## %s (`%s`)\n\n", t.Name, t.Tool)
		if t.Skipped {
			fmt.Fprintf(&b, "SKIPPED: %s\n\n", t.SkipReason)
			continue
		}
		fmt.Fprintf(&b, "- inputs: %d\n", t.Inputs)
		fmt.Fprintf(&b, "- divergences by class: %s\n", classSummary(t.ByClass))
		if t.Coverage.Enabled {
			fmt.Fprintf(&b, "- our-code coverage over the run: %.1f%% of statements\n", t.Coverage.Percent)
		} else if t.Coverage.Detail != "" {
			fmt.Fprintf(&b, "- coverage: not captured (%s)\n", firstLine(t.Coverage.Detail))
		}
		b.WriteString("\n")

		if len(t.Reproducers) == 0 {
			fmt.Fprintf(&b, "No divergences found.\n\n")
			continue
		}
		fmt.Fprintf(&b, "### Minimized reproducers\n\n")
		for i, rp := range t.Reproducers {
			fmt.Fprintf(&b, "**%d. %s** (origin=%s, %d→%d bytes)\n\n", i+1, rp.Class, rp.Origin, rp.OrigLen, rp.InputLen)
			fmt.Fprintf(&b, "%s\n\n", indent(rp.Detail))
			if rp.Printable != "" {
				fmt.Fprintf(&b, "Input:\n\n```\n%s\n```\n\n", rp.Printable)
			} else {
				fmt.Fprintf(&b, "Input (base64): `%s`\n\n", rp.Input)
			}
		}
	}

	return os.WriteFile(path, []byte(b.String()), 0o644)
}

// classSummary renders the per-class divergence counts in a stable order.
func classSummary(by map[DivergenceClass]int) string {
	if len(by) == 0 {
		return "none"
	}
	keys := make([]string, 0, len(by))
	for k := range by {
		keys = append(keys, string(k))
	}
	sort.Strings(keys)
	parts := make([]string, 0, len(keys))
	for _, k := range keys {
		parts = append(parts, fmt.Sprintf("%s=%d", k, by[DivergenceClass(k)]))
	}
	return strings.Join(parts, ", ")
}

func indent(s string) string {
	lines := strings.Split(strings.TrimRight(s, "\n"), "\n")
	for i := range lines {
		lines[i] = "> " + lines[i]
	}
	return strings.Join(lines, "\n")
}

func firstLine(s string) string {
	if i := strings.IndexByte(s, '\n'); i >= 0 {
		return s[:i]
	}
	return s
}
