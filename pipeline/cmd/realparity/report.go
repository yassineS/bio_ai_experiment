package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"
)

// inputInfo records one resolved input file: its path, on-disk size, and the
// contig count read from its header (BAM @SQ count / VCF ##contig count). It
// anchors the manuscript claim that this experiment ran on multi-contig
// (whole-genome) data.
type inputInfo struct {
	Role    string `json:"role"`              // "bam" | "vcf" | "ref"
	Path    string `json:"path"`              // absolute path
	SizeB   int64  `json:"size_bytes"`        // on-disk size
	Contigs int    `json:"contigs,omitempty"` // @SQ / ##contig count (0 if N/A)
}

// sideMeas is one binary's measured resource usage on a cell.
type sideMeas struct {
	WallMS  float64 `json:"wall_ms"`
	CPUMS   float64 `json:"cpu_ms"`
	MaxRSSK int64   `json:"max_rss_kb"`
}

func measToSide(m Measurement) sideMeas {
	return sideMeas{
		WallMS:  float64(m.Wall) / float64(time.Millisecond),
		CPUMS:   float64(m.CPUTotal()) / float64(time.Millisecond),
		MaxRSSK: m.MaxRSSKB,
	}
}

// cellResult is the outcome of one battery cell.
type cellResult struct {
	Name   string `json:"name"`
	Tool   string `json:"tool"`
	Status string `json:"status"` // "PASS" | "DIVERGE" | "SKIP" | "ERROR"
	Detail string `json:"detail,omitempty"`
	Multi  bool   `json:"cross_contig,omitempty"`

	Ours     *sideMeas `json:"ours,omitempty"`
	Upstream *sideMeas `json:"upstream,omitempty"`

	// Ratios are ours/upstream (>1 means ours is slower / heavier). Nil when a
	// side did not run.
	WallRatio *float64 `json:"wall_ratio,omitempty"`
	CPURatio  *float64 `json:"cpu_ratio,omitempty"`
	RSSRatio  *float64 `json:"rss_ratio,omitempty"`

	// DiffSnippet is a short provenance-stripped first-diff excerpt on DIVERGE.
	DiffSnippet string `json:"diff_snippet,omitempty"`
}

// report is the whole-experiment record written as report.{json,md}.
type report struct {
	Generated time.Time    `json:"generated"`
	Region    string       `json:"region,omitempty"`
	Reps      int          `json:"reps"`
	OurBin    string       `json:"our_bin_dir"`
	UpBin     string       `json:"upstream_bin_dir"`
	Inputs    []inputInfo  `json:"inputs"`
	Cells     []cellResult `json:"cells"`

	Pass    int `json:"pass"`
	Diverge int `json:"diverge"`
	Skip    int `json:"skip"`
	Errored int `json:"error"`
}

// verdict reports whether the parity gate held (zero DIVERGE among the cells
// that actually ran on both sides).
func (r *report) verdict() string {
	if r.Diverge > 0 {
		return "FAIL"
	}
	if r.Pass == 0 {
		return "NO-DATA"
	}
	return "PASS"
}

func ratio(num, den float64) *float64 {
	if den <= 0 {
		return nil
	}
	v := num / den
	return &v
}

// finalize tallies the per-status counts and derives the ours/upstream ratios.
func (r *report) finalize() {
	r.Pass, r.Diverge, r.Skip, r.Errored = 0, 0, 0, 0
	for i := range r.Cells {
		c := &r.Cells[i]
		switch c.Status {
		case "PASS":
			r.Pass++
		case "DIVERGE":
			r.Diverge++
		case "SKIP":
			r.Skip++
		default:
			r.Errored++
		}
		if c.Ours != nil && c.Upstream != nil {
			c.WallRatio = ratio(c.Ours.WallMS, c.Upstream.WallMS)
			c.CPURatio = ratio(c.Ours.CPUMS, c.Upstream.CPUMS)
			c.RSSRatio = ratio(float64(c.Ours.MaxRSSK), float64(c.Upstream.MaxRSSK))
		}
	}
}

// writeReports writes report.json and report.md into dir, returning their paths.
func writeReports(r *report, dir string) (jsonPath, mdPath string, err error) {
	if dir == "" {
		dir = "."
	}
	if err = os.MkdirAll(dir, 0o755); err != nil {
		return "", "", err
	}
	jsonPath = filepath.Join(dir, "report.json")
	mdPath = filepath.Join(dir, "report.md")
	b, err := json.MarshalIndent(r, "", "  ")
	if err != nil {
		return "", "", err
	}
	b = append(b, '\n')
	if err := os.WriteFile(jsonPath, b, 0o644); err != nil {
		return "", "", err
	}
	if err := os.WriteFile(mdPath, []byte(r.markdown()), 0o644); err != nil {
		return "", "", err
	}
	return jsonPath, mdPath, nil
}

// markdown renders the human-readable report.
func (r *report) markdown() string {
	var b strings.Builder
	w := func(format string, args ...any) { fmt.Fprintf(&b, format+"\n", args...) }

	w("# Real-world differential parity + performance (whole-genome, multi-contig)")
	w("")
	w("Each cell runs the same command against OUR port and the UPSTREAM binary on")
	w("real multi-contig input, comparing provenance-stripped output (the project's")
	w("exact parity definition, `runner.CompareByteExact`) and timing each side.")
	w("")
	w("- Generated: %s", r.Generated.Format("2006-01-02T15:04:05Z"))
	w("- Reps per cell: %d (wall/CPU = min, RSS = max)", r.Reps)
	w("- Region: %s", orWhole(r.Region))
	w("- Our binaries: `%s`", orNone(r.OurBin))
	w("- Upstream binaries: `%s`", orNone(r.UpBin))
	w("")

	w("## Inputs")
	w("")
	w("| Role | Path | Size | Contigs |")
	w("|------|------|------|---------|")
	if len(r.Inputs) == 0 {
		w("| _(none provided)_ | | | |")
	}
	for _, in := range r.Inputs {
		cc := "—"
		if in.Contigs > 0 {
			cc = fmt.Sprintf("%d", in.Contigs)
		}
		w("| %s | `%s` | %s | %s |", in.Role, in.Path, humanSize(in.SizeB), cc)
	}
	w("")

	w("## Verdict: **%s**", r.verdict())
	w("")
	w("Parity gate = zero DIVERGE among cells that ran on both sides. "+
		"PASS=%d DIVERGE=%d SKIP=%d ERROR=%d.", r.Pass, r.Diverge, r.Skip, r.Errored)
	w("")

	w("## Cells")
	w("")
	w("| Cell | Parity | our wall/CPU/RSS | up wall/CPU/RSS | wall× | CPU× | RSS× |")
	w("|------|--------|------------------|-----------------|-------|------|------|")
	for _, c := range r.Cells {
		w("| %s | %s | %s | %s | %s | %s | %s |",
			cellLabel(c), c.Status,
			sideCell(c.Ours), sideCell(c.Upstream),
			ratioCell(c.WallRatio), ratioCell(c.CPURatio), ratioCell(c.RSSRatio))
	}
	w("")

	// Divergence detail with embedded provenance-stripped diff snippets.
	var diverged []cellResult
	for _, c := range r.Cells {
		if c.Status == "DIVERGE" || c.Status == "ERROR" {
			diverged = append(diverged, c)
		}
	}
	if len(diverged) > 0 {
		w("## Divergences / errors")
		w("")
		for _, c := range diverged {
			w("### %s — %s", c.Name, c.Status)
			w("")
			if c.Detail != "" {
				w("%s", c.Detail)
				w("")
			}
			if c.DiffSnippet != "" {
				w("```")
				w("%s", c.DiffSnippet)
				w("```")
				w("")
			}
		}
	}
	return b.String()
}

// cellLabel marks cross-contig cells so the report makes the multi-ref coverage
// explicit.
func cellLabel(c cellResult) string {
	if c.Multi {
		return c.Name + " †"
	}
	return c.Name
}

func sideCell(s *sideMeas) string {
	if s == nil {
		return "—"
	}
	return fmt.Sprintf("%.1f / %.1f ms / %s", s.WallMS, s.CPUMS, humanSizeKB(s.MaxRSSK))
}

func ratioCell(r *float64) string {
	if r == nil {
		return "—"
	}
	return fmt.Sprintf("%.2f", *r)
}

func orNone(s string) string {
	if s == "" {
		return "(none)"
	}
	return s
}

func orWhole(s string) string {
	if s == "" {
		return "(whole-genome, all contigs)"
	}
	return s
}

func humanSize(n int64) string {
	const unit = 1024
	if n < unit {
		return fmt.Sprintf("%d B", n)
	}
	div, exp := int64(unit), 0
	for v := n / unit; v >= unit; v /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.1f %ciB", float64(n)/float64(div), "KMGTPE"[exp])
}

func humanSizeKB(kb int64) string { return humanSize(kb * 1024) }
