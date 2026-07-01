package realbench

import (
	"encoding/json"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"time"
)

// SideMeasure is one binary's measured resource usage on a cell, in the
// report's JSON units (seconds for wall/CPU, KiB for RSS).
type SideMeasure struct {
	WallS float64 `json:"wall_s"`
	CPUS  float64 `json:"cpu_s"`
	RSSKB int64   `json:"rss_kb"`
}

func measToSide(m Measurement) *SideMeasure {
	return &SideMeasure{
		WallS: m.Wall.Seconds(),
		CPUS:  m.CPUTotal().Seconds(),
		RSSKB: m.MaxRSSKB,
	}
}

// CellRecord is the per-cell result row. Parity is one of PASS/DIFF/SKIP/ERROR.
// wall_x/cpu_x/rss_x are ours/upstream (>1 means ours is slower / heavier); they
// are nil when a side did not run.
type CellRecord struct {
	Tool       string `json:"tool"`
	Name       string `json:"name"`
	Subcommand string `json:"subcommand"`
	Tier       string `json:"tier"`
	Parity     string `json:"parity"` // PASS | DIFF | SKIP | ERROR

	Ours *SideMeasure `json:"ours,omitempty"`
	Up   *SideMeasure `json:"up,omitempty"`

	WallX *float64 `json:"wall_x,omitempty"`
	CPUX  *float64 `json:"cpu_x,omitempty"`
	RSSX  *float64 `json:"rss_x,omitempty"`

	Note string `json:"note,omitempty"`
}

// Machine records the host the benchmark ran on, so results from different
// tiers/runs are interpretable when aggregated.
type Machine struct {
	OS        string `json:"os"`
	Arch      string `json:"arch"`
	NumCPU    int    `json:"num_cpu"`
	GoVersion string `json:"go_version"`
	Hostname  string `json:"hostname,omitempty"`
}

// machineInfo captures the current host.
func machineInfo() Machine {
	h, _ := os.Hostname()
	return Machine{
		OS:        runtime.GOOS,
		Arch:      runtime.GOARCH,
		NumCPU:    runtime.NumCPU(),
		GoVersion: runtime.Version(),
		Hostname:  h,
	}
}

// Report is the whole-experiment record written as realbench.<tier>.{json,md}.
// It is a single object so per-tier reports aggregate trivially by concatenating
// their Cells slices.
type Report struct {
	Tier      string       `json:"tier"`
	Generated time.Time    `json:"generated"`
	Reps      int          `json:"reps"`
	Machine   Machine      `json:"machine"`
	Cells     []CellRecord `json:"cells"`

	Pass    int `json:"pass"`
	Diff    int `json:"diff"`
	Skip    int `json:"skip"`
	Errored int `json:"error"`
}

// ratio returns num/den (>0 den) or nil.
func ratio(num, den float64) *float64 {
	if den <= 0 {
		return nil
	}
	v := num / den
	return &v
}

// finalize tallies per-parity counts and derives the ours/upstream ratios.
func (r *Report) finalize() {
	r.Pass, r.Diff, r.Skip, r.Errored = 0, 0, 0, 0
	for i := range r.Cells {
		c := &r.Cells[i]
		switch c.Parity {
		case "PASS":
			r.Pass++
		case "DIFF":
			r.Diff++
		case "SKIP":
			r.Skip++
		default:
			r.Errored++
		}
		if c.Ours != nil && c.Up != nil {
			c.WallX = ratio(c.Ours.WallS, c.Up.WallS)
			c.CPUX = ratio(c.Ours.CPUS, c.Up.CPUS)
			c.RSSX = ratio(float64(c.Ours.RSSKB), float64(c.Up.RSSKB))
		}
	}
}

// WriteReports writes realbench.<tier>.json and realbench.<tier>.md into dir,
// returning their paths.
func WriteReports(r *Report, dir string) (jsonPath, mdPath string, err error) {
	if dir == "" {
		dir = "."
	}
	// The working directory (".") always exists, so skip creating it. This also
	// avoids a quirk on some FUSE-backed filesystems (e.g. Fusion on AWS Batch),
	// where MkdirAll(".") spuriously fails with "file exists" — which would throw
	// away an entire completed run at the final report-writing step. Only create
	// a real subdirectory, and tolerate it already existing.
	if dir != "." {
		if err = os.MkdirAll(dir, 0o755); err != nil && !errors.Is(err, fs.ErrExist) {
			return "", "", err
		}
	}
	jsonPath = filepath.Join(dir, "realbench."+r.Tier+".json")
	mdPath = filepath.Join(dir, "realbench."+r.Tier+".md")
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

// markdown renders the human-readable table, sorted by rss_x descending (cells
// without an rss_x ratio sink to the bottom).
func (r *Report) markdown() string {
	var b strings.Builder
	w := func(format string, args ...any) { fmt.Fprintf(&b, format+"\n", args...) }

	w("# realbench — real-data parity + performance (tier: %s)", r.Tier)
	w("")
	w("Each cell runs the same operation against OUR port and the UPSTREAM binary")
	w("on real %s-tier data, comparing provenance-stripped output (the project's", r.Tier)
	w("byte-exact parity definition) and timing each side (wall/CPU = min, RSS = max).")
	w("")
	w("- Generated: %s", r.Generated.Format("2006-01-02T15:04:05Z"))
	w("- Reps per cell: %d", r.Reps)
	w("- Machine: %s/%s, %d CPU, %s", r.Machine.OS, r.Machine.Arch, r.Machine.NumCPU, r.Machine.GoVersion)
	w("")
	w("PASS=%d DIFF=%d SKIP=%d ERROR=%d (of %d cells)", r.Pass, r.Diff, r.Skip, r.Errored, len(r.Cells))
	w("")

	// Sort a copy by rss_x desc; nil ratios sort last.
	rows := make([]CellRecord, len(r.Cells))
	copy(rows, r.Cells)
	sort.SliceStable(rows, func(i, j int) bool {
		return rssKey(rows[i]) > rssKey(rows[j])
	})

	w("| Cell | Tool | Sub | Parity | our wall/cpu/rss | up wall/cpu/rss | wall× | cpu× | rss× | note |")
	w("|------|------|-----|--------|------------------|-----------------|-------|------|------|------|")
	for _, c := range rows {
		w("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
			c.Name, c.Tool, c.Subcommand, c.Parity,
			sideCell(c.Ours), sideCell(c.Up),
			ratioCell(c.WallX), ratioCell(c.CPUX), ratioCell(c.RSSX),
			truncNote(c.Note))
	}
	w("")
	return b.String()
}

// rssKey is the sort key for the MD table: the rss_x ratio, or -1 for cells
// without one (so they sink to the bottom).
func rssKey(c CellRecord) float64 {
	if c.RSSX != nil {
		return *c.RSSX
	}
	return -1
}

func sideCell(s *SideMeasure) string {
	if s == nil {
		return "—"
	}
	return fmt.Sprintf("%.3f / %.3f s / %s", s.WallS, s.CPUS, humanSizeKB(s.RSSKB))
}

func ratioCell(r *float64) string {
	if r == nil {
		return "—"
	}
	return fmt.Sprintf("%.2f", *r)
}

func truncNote(s string) string {
	s = strings.ReplaceAll(s, "\n", " ")
	s = strings.ReplaceAll(s, "|", "/")
	if len(s) > 120 {
		return s[:120] + "…"
	}
	return s
}

func humanSizeKB(kb int64) string {
	const unit = 1024
	n := kb
	if n < unit {
		return fmt.Sprintf("%d KiB", n)
	}
	div, exp := int64(unit), 0
	for v := n / unit; v >= unit; v /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.1f %ciB", float64(n)/float64(div), "MGTPE"[exp])
}
