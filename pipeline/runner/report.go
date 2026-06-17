package runner

import (
	"encoding/json"
	"fmt"
	"os"
	"sort"
	"strings"
	"time"
)

// Report aggregates the results of a pipeline run.
type Report struct {
	Scale     string    `json:"scale"`
	Seed      int64     `json:"seed"`
	Generated time.Time `json:"generated"`
	Results   []Result  `json:"results"`
	Summary   Summary   `json:"summary"`
}

// Summary is the run-level tally.
type Summary struct {
	Total   int `json:"total"`
	Pass    int `json:"pass"`
	Similar int `json:"similar"`
	Diverge int `json:"diverge"`
	Skip    int `json:"skip"`
	Error   int `json:"error"`
}

// Build assembles a Report from results, computing the summary.
func Build(scale string, seed int64, results []Result) Report {
	r := Report{Scale: scale, Seed: seed, Generated: time.Now().UTC(), Results: results}
	for _, res := range results {
		r.Summary.Total++
		switch res.Status {
		case StatusPass:
			r.Summary.Pass++
		case StatusSimilar:
			r.Summary.Similar++
		case StatusDiverge:
			r.Summary.Diverge++
		case StatusSkip:
			r.Summary.Skip++
		case StatusError:
			r.Summary.Error++
		}
	}
	return r
}

// Failed reports whether the run had any real failures (DIVERGE or ERROR).
// Similarity matches and skips do not count as failures.
func (r Report) Failed() bool {
	return r.Summary.Diverge > 0 || r.Summary.Error > 0
}

// WriteJSON writes the machine-readable report.
func (r Report) WriteJSON(path string) error {
	b, err := json.MarshalIndent(r, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}

// WriteMarkdown writes the human-readable report grouped by tool, with PASS/
// DIVERGE/SIMILAR/SKIP counts and, for heavy entries, the ours-vs-upstream
// timing ratio.
func (r Report) WriteMarkdown(path string) error {
	var b strings.Builder
	fmt.Fprintf(&b, "# Parity pipeline report\n\n")
	fmt.Fprintf(&b, "- Scale: `%s`\n- Seed: `%d`\n- Generated: %s\n\n",
		r.Scale, r.Seed, r.Generated.Format(time.RFC3339))
	fmt.Fprintf(&b, "## Summary\n\n")
	fmt.Fprintf(&b, "| total | PASS | SIMILAR | DIVERGE | SKIP | ERROR |\n")
	fmt.Fprintf(&b, "|------:|-----:|--------:|--------:|-----:|------:|\n")
	fmt.Fprintf(&b, "| %d | %d | %d | %d | %d | %d |\n\n",
		r.Summary.Total, r.Summary.Pass, r.Summary.Similar, r.Summary.Diverge, r.Summary.Skip, r.Summary.Error)

	// Group by tool, stable order.
	byTool := map[string][]Result{}
	var tools []string
	for _, res := range r.Results {
		if _, ok := byTool[res.Tool]; !ok {
			tools = append(tools, res.Tool)
		}
		byTool[res.Tool] = append(byTool[res.Tool], res)
	}
	sort.Strings(tools)

	for _, tool := range tools {
		rs := byTool[tool]
		var p, s, d, sk, e int
		for _, res := range rs {
			switch res.Status {
			case StatusPass:
				p++
			case StatusSimilar:
				s++
			case StatusDiverge:
				d++
			case StatusSkip:
				sk++
			case StatusError:
				e++
			}
		}
		fmt.Fprintf(&b, "## %s\n\n", tool)
		fmt.Fprintf(&b, "PASS %d · SIMILAR %d · DIVERGE %d · SKIP %d · ERROR %d\n\n", p, s, d, sk, e)
		fmt.Fprintf(&b, "| entry | input | compare | status | ours(ms) | upstream(ms) | ratio | detail |\n")
		fmt.Fprintf(&b, "|-------|-------|---------|--------|---------:|-------------:|------:|--------|\n")
		for _, res := range rs {
			ratio := ""
			if res.Heavy && res.TimingRatio > 0 {
				ratio = fmt.Sprintf("%.2fx", res.TimingRatio)
			}
			detail := strings.ReplaceAll(res.Detail, "\n", " ")
			if len(detail) > 120 {
				detail = detail[:120] + "..."
			}
			fmt.Fprintf(&b, "| %s | %s | %s | %s | %d | %d | %s | %s |\n",
				res.Name, res.Input, res.Compare, res.Status, res.OursMillis, res.UpstreamMs, ratio, detail)
		}
		b.WriteString("\n")
	}
	return os.WriteFile(path, []byte(b.String()), 0o644)
}
