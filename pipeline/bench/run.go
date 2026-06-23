package bench

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pipeline/stats"
)

// ratioCIAlpha is the two-sided significance level for the bootstrap CI on the
// wall-time ours/upstream ratio: 0.05 gives a 95% interval, matching the
// binomial CIs used elsewhere in the manuscript.
const ratioCIAlpha = 0.05

// RepSamples holds the RAW per-repetition measurements for one side of a cell,
// in execution order. Every element is one process run; the slices are parallel
// (index i is the same run across all three axes). Times are milliseconds,
// memory is mebibytes. These are what the manuscript's distribution-level
// statistics (median, IQR, bootstrap ratio CI) are computed from — the reduced
// min/max fields on CellResult are a lossy summary kept for backward
// compatibility.
type RepSamples struct {
	WallMs []float64 `json:"wall_ms"`
	CPUMs  []float64 `json:"cpu_ms"`
	RSSMB  []float64 `json:"rss_mb"`
}

// CellResult is the measured cost of one (cell, scale) pair for both sides plus
// the ours/upstream ratios. Times are milliseconds, memory is mebibytes. A
// ratio < 1 means OUR binary used less of that resource than upstream.
//
// Three layers of detail are recorded:
//   - the reduced min/max point fields (Our*/Up*) and their min-based ratios —
//     the original, backward-compatible summary;
//   - the RAW per-rep samples (OurSamples/UpSamples) — the full distribution;
//   - the robust distribution statistics derived from those samples
//     (median, IQR, and a bootstrap CI on the median ours/upstream ratio).
type CellResult struct {
	Cell  string `json:"cell"`
	Group string `json:"group"`
	Scale string `json:"scale"`
	Reps  int    `json:"reps"`

	OurWallMs float64 `json:"our_wall_ms"`
	UpWallMs  float64 `json:"up_wall_ms"`
	OurCPUMs  float64 `json:"our_cpu_ms"`
	UpCPUMs   float64 `json:"up_cpu_ms"`
	OurRSSMB  float64 `json:"our_rss_mb"`
	UpRSSMB   float64 `json:"up_rss_mb"`

	WallRatio float64 `json:"wall_ratio"` // ours/upstream (min-based, legacy)
	CPURatio  float64 `json:"cpu_ratio"`
	RSSRatio  float64 `json:"rss_ratio"`

	// Raw per-rep samples for both sides (every rep, every axis).
	OurSamples RepSamples `json:"our_samples"`
	UpSamples  RepSamples `json:"up_samples"`

	// Robust distribution statistics over the raw samples. Median/IQR are the
	// type-7 quantiles (pipeline/stats.MedianIQR); the *RatioMed fields are
	// median(ours)/median(upstream); the wall-ratio bootstrap CI bounds the
	// headline speed claim. See computeStats.
	OurWallMed float64 `json:"our_wall_med"`
	OurWallIQR float64 `json:"our_wall_iqr"`
	UpWallMed  float64 `json:"up_wall_med"`
	UpWallIQR  float64 `json:"up_wall_iqr"`
	OurCPUMed  float64 `json:"our_cpu_med"`
	OurCPUIQR  float64 `json:"our_cpu_iqr"`
	UpCPUMed   float64 `json:"up_cpu_med"`
	UpCPUIQR   float64 `json:"up_cpu_iqr"`
	OurRSSMed  float64 `json:"our_rss_med"`
	OurRSSIQR  float64 `json:"our_rss_iqr"`
	UpRSSMed   float64 `json:"up_rss_med"`
	UpRSSIQR   float64 `json:"up_rss_iqr"`

	WallRatioMed  float64 `json:"wall_ratio_med"` // median(ours)/median(up)
	CPURatioMed   float64 `json:"cpu_ratio_med"`
	RSSRatioMed   float64 `json:"rss_ratio_med"`
	WallRatioCILo float64 `json:"wall_ratio_ci_lo"` // 95% bootstrap CI on wall ratio
	WallRatioCIHi float64 `json:"wall_ratio_ci_hi"`

	Err string `json:"err,omitempty"`
}

// RunConfig controls a benchmark sweep.
type RunConfig struct {
	Scales   []fixtures.Scale // tiers to sweep (small/medium/large) — the scalability axis
	Reps     int              // repetitions per side; wall/CPU take the min, RSS the max
	Seed     int64            // fixture seed
	CellGlob string           // substring filter on cell name ("" = all)
	GroupSel string           // file-type group filter ("" = all)
	Log      func(string, ...any)
}

// Run executes the bench matrix across the configured scale tiers and returns
// one CellResult per (cell, scale). Fixtures for each tier are generated (and
// cached) on demand, exactly like the parity runner, so the two share inputs.
func Run(cfg RunConfig) ([]CellResult, error) {
	if cfg.Reps < 1 {
		cfg.Reps = 3
	}
	if cfg.Log == nil {
		cfg.Log = func(string, ...any) {}
	}
	root, err := upstream.RepoRoot()
	if err != nil {
		return nil, err
	}
	cache := filepath.Join(root, "pipeline", ".fixtures", "bin")

	cells := BenchMatrix()
	var out []CellResult
	for _, scale := range cfg.Scales {
		cfg.Log("generating %s fixtures…", scale)
		man, err := fixtures.Generate(fixtures.Options{Scale: scale, Seed: cfg.Seed})
		if err != nil {
			return nil, fmt.Errorf("fixtures %s: %w", scale, err)
		}
		for _, c := range cells {
			if cfg.CellGlob != "" && !strings.Contains(c.Name, cfg.CellGlob) {
				continue
			}
			if cfg.GroupSel != "" && !strings.EqualFold(c.Group, cfg.GroupSel) {
				continue
			}
			res := runCell(c, man, scale, cfg.Reps, cache)
			status := fmt.Sprintf("wall x%.2f [%.2f,%.2f] (min x%.2f) cpu x%.2f rss x%.2f",
				res.WallRatioMed, res.WallRatioCILo, res.WallRatioCIHi,
				res.WallRatio, res.CPURatioMed, res.RSSRatioMed)
			if res.Err != "" {
				status = "ERROR: " + res.Err
			}
			cfg.Log("[%-6s] %-22s %s", scale, c.Name, status)
			out = append(out, res)
		}
	}
	return out, nil
}

// runCell resolves both binaries, builds the invocation, and measures each side.
func runCell(c BenchCell, man *fixtures.Manifest, scale fixtures.Scale, reps int, cache string) CellResult {
	r := CellResult{Cell: c.Name, Group: c.Group, Scale: string(scale), Reps: reps}
	tmp, err := os.MkdirTemp("", "bench-"+c.Name+"-")
	if err != nil {
		r.Err = err.Error()
		return r
	}
	defer os.RemoveAll(tmp)

	p := c.Build(man, tmp)
	ourBin, err := upstream.OurBinary(p.ourTool, cache)
	if err != nil {
		r.Err = "our binary: " + err.Error()
		return r
	}
	upBin, err := upstream.Binary(p.upKey)
	if err != nil {
		r.Err = "upstream binary: " + err.Error()
		return r
	}

	om, ourSamples, err := repeatMeasuredSamples(reps, ourBin, p.ourArgs, p.stdin, p.ourStdout)
	if err != nil {
		r.Err = "ours: " + err.Error()
		return r
	}
	um, upSamples, err := repeatMeasuredSamples(reps, upBin, p.upArgs, p.stdin, p.upStdout)
	if err != nil {
		r.Err = "upstream: " + err.Error()
		return r
	}

	// Legacy min/max point summary (kept for backward compatibility).
	r.OurWallMs = ms(om.Wall)
	r.UpWallMs = ms(um.Wall)
	r.OurCPUMs = ms(om.CPUTotal())
	r.UpCPUMs = ms(um.CPUTotal())
	r.OurRSSMB = mb(om.MaxRSSKB)
	r.UpRSSMB = mb(um.MaxRSSKB)
	r.WallRatio = ratio(r.OurWallMs, r.UpWallMs)
	r.CPURatio = ratio(r.OurCPUMs, r.UpCPUMs)
	r.RSSRatio = ratio(r.OurRSSMB, r.UpRSSMB)

	// Raw per-rep samples + robust distribution statistics.
	r.OurSamples = toRepSamples(ourSamples)
	r.UpSamples = toRepSamples(upSamples)
	r.computeStats()
	return r
}

// toRepSamples converts the raw per-rep Measurements into parallel float slices
// in report units (ms, ms, MiB).
func toRepSamples(ms_ []Measurement) RepSamples {
	rs := RepSamples{
		WallMs: make([]float64, len(ms_)),
		CPUMs:  make([]float64, len(ms_)),
		RSSMB:  make([]float64, len(ms_)),
	}
	for i, m := range ms_ {
		rs.WallMs[i] = ms(m.Wall)
		rs.CPUMs[i] = ms(m.CPUTotal())
		rs.RSSMB[i] = mb(m.MaxRSSKB)
	}
	return rs
}

// computeStats fills the median/IQR fields for both sides and all three axes,
// the median-based ratios, and the bootstrap CI on the wall-time ratio, from
// the already-populated raw samples. Median and IQR use the type-7 quantile
// (pipeline/stats.MedianIQR); the wall-ratio CI is the percentile bootstrap of
// pipeline/stats.RatioCI (fixed-seed, hence deterministic).
func (r *CellResult) computeStats() {
	r.OurWallMed, r.OurWallIQR = medIQR(r.OurSamples.WallMs)
	r.UpWallMed, r.UpWallIQR = medIQR(r.UpSamples.WallMs)
	r.OurCPUMed, r.OurCPUIQR = medIQR(r.OurSamples.CPUMs)
	r.UpCPUMed, r.UpCPUIQR = medIQR(r.UpSamples.CPUMs)
	r.OurRSSMed, r.OurRSSIQR = medIQR(r.OurSamples.RSSMB)
	r.UpRSSMed, r.UpRSSIQR = medIQR(r.UpSamples.RSSMB)

	r.CPURatioMed = ratio(r.OurCPUMed, r.UpCPUMed)
	r.RSSRatioMed = ratio(r.OurRSSMed, r.UpRSSMed)
	point, lo, hi := stats.RatioCI(r.OurSamples.WallMs, r.UpSamples.WallMs, ratioCIAlpha)
	r.WallRatioMed = point
	r.WallRatioCILo = lo
	r.WallRatioCIHi = hi
}

// medIQR returns the median and the inter-quartile range (Q3-Q1) of xs.
func medIQR(xs []float64) (med, iqr float64) {
	m, q1, q3 := stats.MedianIQR(xs)
	return m, q3 - q1
}

func ms(d time.Duration) float64 { return float64(d.Microseconds()) / 1000.0 }
func mb(kib int64) float64       { return float64(kib) / 1024.0 }
func ratio(our, up float64) float64 {
	if up == 0 {
		return 0
	}
	return our / up
}

// WriteJSON serialises the results to path, including the RAW per-rep samples
// (our_samples / up_samples) and the derived median/IQR and bootstrap ratio-CI
// fields alongside the legacy min/max point summary.
func WriteJSON(path string, results []CellResult) error {
	b, err := json.MarshalIndent(results, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}

// WriteMarkdown renders a human-readable report: a per-(group,scale) table of
// ours-vs-upstream wall / CPU / RSS with ratios, then a scalability section that
// shows how each cell's wall-time grows across the swept tiers.
func WriteMarkdown(path string, results []CellResult, scales []fixtures.Scale) error {
	var b strings.Builder
	b.WriteString("# Performance & scalability report\n\n")
	b.WriteString("Generated by `pipeline/bench` (`cmd/parity-bench`). Each cell times OUR binary\n")
	b.WriteString("against the vendored UPSTREAM binary over N repetitions. Two summaries are\n")
	b.WriteString("reported: the legacy point estimate (wall/CPU = minimum over reps, RSS = peak)\n")
	b.WriteString("and — the headline for the manuscript — the robust distribution\n")
	b.WriteString("(median ± IQR per side, the median ours/upstream ratio, and a 95% percentile\n")
	b.WriteString("bootstrap CI on the wall-time ratio). `ratio = ours / upstream` (< 1.0 = we use\n")
	b.WriteString("less). The raw per-rep samples for every cell are in the companion `bench.json`.\n")
	b.WriteString("Slow cells (ratio > 1) are reported plainly, not hidden.\n\n")

	// Per-scale resource tables, grouped by file format.
	for _, scale := range scales {
		rows := filterScale(results, string(scale))
		if len(rows) == 0 {
			continue
		}
		fmt.Fprintf(&b, "## Scale: %s\n\n", scale)
		groups := orderedGroups(rows)
		for _, g := range groups {
			fmt.Fprintf(&b, "### %s\n\n", g)
			// Robust (median/IQR/CI) table — the headline.
			b.WriteString("**Robust (median ± IQR over reps; 95% bootstrap CI on wall ratio):**\n\n")
			b.WriteString("| cell | wall ms ours (med±IQR) | wall ms up (med±IQR) | wall× (med) | wall× 95% CI | CPU× (med) | RSS× (med) |\n")
			b.WriteString("|---|---|---|---|---|---|---|\n")
			for _, r := range rows {
				if r.Group != g {
					continue
				}
				if r.Err != "" {
					fmt.Fprintf(&b, "| %s | ERROR: %s | | | | | |\n", r.Cell, r.Err)
					continue
				}
				slow := ""
				if r.WallRatioMed > 1.0 {
					slow = " ⚠" // flag regressions plainly
				}
				fmt.Fprintf(&b, "| %s | %.1f ± %.1f | %.1f ± %.1f | %.2f%s | [%.2f, %.2f] | %.2f | %.2f |\n",
					r.Cell,
					r.OurWallMed, r.OurWallIQR, r.UpWallMed, r.UpWallIQR,
					r.WallRatioMed, slow, r.WallRatioCILo, r.WallRatioCIHi,
					r.CPURatioMed, r.RSSRatioMed)
			}
			b.WriteString("\n")
			// Legacy point-estimate table (min wall/CPU, peak RSS) — kept for
			// backward compatibility with earlier reports.
			b.WriteString("**Point estimate (min wall/CPU, peak RSS) — legacy:**\n\n")
			b.WriteString("| cell | wall ms (ours/up) | wall× | CPU ms (ours/up) | CPU× | RSS MB (ours/up) | RSS× |\n")
			b.WriteString("|---|---|---|---|---|---|---|\n")
			for _, r := range rows {
				if r.Group != g {
					continue
				}
				if r.Err != "" {
					fmt.Fprintf(&b, "| %s | ERROR | | %s | | | |\n", r.Cell, r.Err)
					continue
				}
				fmt.Fprintf(&b, "| %s | %.1f / %.1f | %.2f | %.1f / %.1f | %.2f | %.1f / %.1f | %.2f |\n",
					r.Cell, r.OurWallMs, r.UpWallMs, r.WallRatio,
					r.OurCPUMs, r.UpCPUMs, r.CPURatio,
					r.OurRSSMB, r.UpRSSMB, r.RSSRatio)
			}
			b.WriteString("\n")
		}
	}

	// Scalability: wall-time per cell across tiers (only when >1 tier was swept).
	if len(scales) > 1 {
		b.WriteString("## Scalability — wall-clock ms across tiers (ours | upstream)\n\n")
		cells := orderedCells(results)
		b.WriteString("| cell |")
		for _, s := range scales {
			fmt.Fprintf(&b, " %s |", s)
		}
		b.WriteString("\n|---|")
		for range scales {
			b.WriteString("---|")
		}
		b.WriteString("\n")
		for _, cell := range cells {
			fmt.Fprintf(&b, "| %s |", cell)
			for _, s := range scales {
				r := find(results, cell, string(s))
				if r == nil || r.Err != "" {
					b.WriteString(" — |")
					continue
				}
				fmt.Fprintf(&b, " %.0f \\| %.0f |", r.OurWallMs, r.UpWallMs)
			}
			b.WriteString("\n")
		}
		b.WriteString("\n")
	}

	return os.WriteFile(path, []byte(b.String()), 0o644)
}

func filterScale(rs []CellResult, scale string) []CellResult {
	var out []CellResult
	for _, r := range rs {
		if r.Scale == scale {
			out = append(out, r)
		}
	}
	return out
}

func orderedGroups(rs []CellResult) []string {
	seen := map[string]bool{}
	var out []string
	for _, r := range rs {
		if !seen[r.Group] {
			seen[r.Group] = true
			out = append(out, r.Group)
		}
	}
	sort.Strings(out)
	return out
}

func orderedCells(rs []CellResult) []string {
	seen := map[string]bool{}
	var out []string
	for _, r := range rs {
		if !seen[r.Cell] {
			seen[r.Cell] = true
			out = append(out, r.Cell)
		}
	}
	return out
}

func find(rs []CellResult, cell, scale string) *CellResult {
	for i := range rs {
		if rs[i].Cell == cell && rs[i].Scale == scale {
			return &rs[i]
		}
	}
	return nil
}
