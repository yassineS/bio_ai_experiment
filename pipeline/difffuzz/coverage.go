package difffuzz

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// CoverageResult reports the Go statement/branch coverage our tool binary
// exercised over the fuzz run, when coverage capture was enabled.
type CoverageResult struct {
	Enabled bool    `json:"enabled"`
	Tool    string  `json:"tool"`
	Percent float64 `json:"percent"` // statements covered, percent
	Detail  string  `json:"detail,omitempty"`
}

// coverageBuilder builds our tool with -cover and accumulates coverage data
// into a GOCOVERDIR across the whole fuzz run, then renders a percentage.
//
// Mechanism (Go 1.20+ binary coverage, available in the project's toolchain):
//   - build the tool binary with `go build -cover -o <bin> <pkg>`;
//   - run that binary with GOCOVERDIR pointed at a per-run directory, so each
//     invocation appends coverage counters there;
//   - after the run, `go tool covdata percent -i=<dir>` prints per-package
//     coverage, from which we extract the total.
//
// If anything in this path is unavailable (older toolchain, build failure) the
// builder degrades gracefully: it returns a non-instrumented binary and the
// final CoverageResult reports Enabled=false with a Detail explaining why, so
// the rest of the fuzz run is unaffected. This is the "leave a hook, document
// how to get it" requirement satisfied with a working default.
type coverageBuilder struct {
	enabled  bool
	tool     string
	bin      string // instrumented binary (falls back to plain bin on failure)
	covDir   string // GOCOVERDIR
	buildErr string
}

// newCoverageBuilder builds an instrumented binary for tool under workDir. When
// enable is false (or the build fails) it returns a builder that yields the
// ordinary (uninstrumented) binary and reports Enabled=false.
func newCoverageBuilder(enable bool, tool, plainBin, workDir string) *coverageBuilder {
	cb := &coverageBuilder{tool: tool, bin: plainBin}
	if !enable {
		return cb
	}
	covDir := filepath.Join(workDir, "covdata-"+tool)
	if err := os.MkdirAll(covDir, 0o755); err != nil {
		cb.buildErr = "mkdir covdir: " + err.Error()
		return cb
	}
	root, err := upstream.RepoRoot()
	if err != nil {
		cb.buildErr = err.Error()
		return cb
	}
	out := filepath.Join(workDir, tool+"-cover")
	pkg := fmt.Sprintf("github.com/yassineS/bio_ai_experiment/tools/%s/cmd/%s", tool, tool)
	cmd := exec.Command("go", "build", "-cover", "-o", out, pkg)
	cmd.Dir = root
	if b, err := cmd.CombinedOutput(); err != nil {
		cb.buildErr = fmt.Sprintf("build -cover failed (need Go >=1.20): %v\n%s", err, strings.TrimSpace(string(b)))
		return cb
	}
	cb.enabled = true
	cb.bin = out
	cb.covDir = covDir
	return cb
}

// env returns the environment to run the instrumented binary with (GOCOVERDIR
// set). Returns nil when coverage is disabled (caller then leaves env default).
func (cb *coverageBuilder) env() []string {
	if !cb.enabled {
		return nil
	}
	return append(os.Environ(), "GOCOVERDIR="+cb.covDir)
}

// result renders the accumulated coverage into a CoverageResult.
func (cb *coverageBuilder) result() CoverageResult {
	if !cb.enabled {
		detail := "coverage capture disabled"
		if cb.buildErr != "" {
			detail = "coverage capture unavailable: " + cb.buildErr +
				"\nTo capture manually: build the tool with `go build -cover`, run it with " +
				"GOCOVERDIR set to a dir, then `go tool covdata percent -i=<dir>`."
		}
		return CoverageResult{Enabled: false, Tool: cb.tool, Detail: detail}
	}
	cmd := exec.Command("go", "tool", "covdata", "percent", "-i="+cb.covDir)
	out, err := cmd.CombinedOutput()
	if err != nil {
		return CoverageResult{Enabled: false, Tool: cb.tool,
			Detail: fmt.Sprintf("covdata percent failed: %v\n%s", err, strings.TrimSpace(string(out)))}
	}
	pct, detail := parseCovdataPercent(out)
	return CoverageResult{Enabled: true, Tool: cb.tool, Percent: pct, Detail: detail}
}

// parseCovdataPercent extracts the overall coverage percent from `go tool
// covdata percent` output. Each line looks like:
//
//	github.com/.../tools/seqtk/pkg/seqtk	coverage: 42.1% of statements
//
// We average the per-package percentages weighted equally (covdata does not
// emit a single total), which is a representative figure for a single tool's
// packages; the raw per-package lines are kept in Detail for transparency.
func parseCovdataPercent(out []byte) (float64, string) {
	lines := strings.Split(strings.TrimSpace(string(out)), "\n")
	var sum float64
	var n int
	for _, ln := range lines {
		idx := strings.Index(ln, "coverage:")
		if idx < 0 {
			continue
		}
		rest := ln[idx+len("coverage:"):]
		pctStr := strings.TrimSpace(rest)
		if sp := strings.IndexByte(pctStr, '%'); sp >= 0 {
			pctStr = pctStr[:sp]
		}
		pctStr = strings.TrimSpace(pctStr)
		if v, err := strconv.ParseFloat(pctStr, 64); err == nil {
			sum += v
			n++
		}
	}
	detail := strings.TrimSpace(string(out))
	if n == 0 {
		return 0, detail
	}
	return sum / float64(n), detail
}
