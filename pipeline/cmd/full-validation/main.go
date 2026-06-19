// Command full-validation is the end-to-end validation flow: for each scale
// tier it (1) runs the entire parity matrix — every command × its flag-combo
// cells — comparing our output to upstream byte-for-byte (or within a numeric
// tolerance), (2) runs the format round-trip checks (lossless encode/decode for
// BAM/CRAM/BGZF/VCF-BCF/FASTQ), and (3) runs the performance/scalability bench
// (wall, CPU, peak RSS per cell). It writes a consolidated report per scale and
// exits non-zero if anything diverges, so it can gate a release.
//
// It is designed to be run unattended on a fat node for the large tier (which
// needs ~30 GB+ scratch). Typical invocations:
//
//	# everything, large tier (manuscript / release gate):
//	go run ./pipeline/cmd/full-validation -scales=large -reps=5
//
//	# quick local sanity across the small tiers:
//	go run ./pipeline/cmd/full-validation -scales=smoke,small -reps=2
//
//	# one tool, one tier while iterating:
//	go run ./pipeline/cmd/full-validation -scales=medium -tools=bcftools
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/bench"
	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
	"github.com/yassineS/bio_ai_experiment/pipeline/roundtrip"
	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

func main() {
	var (
		scalesFlag = flag.String("scales", "large", "comma-separated tiers: smoke|small|medium|large")
		toolsFlag  = flag.String("tools", "", "comma-separated tool filter for the parity matrix (default: all)")
		repsFlag   = flag.Int("reps", 5, "bench repetitions per side")
		seedFlag   = flag.Int64("seed", fixtures.DefaultSeed, "fixture RNG seed")
		updateFix  = flag.Bool("update-fixtures", false, "force fixture regeneration")
		skipBench  = flag.Bool("skip-bench", false, "skip the performance/RSS bench stage")
		outFlag    = flag.String("out", "", "report root (default: pipeline/.fixtures/<scale>/full-validation)")
	)
	flag.Parse()

	root, err := upstream.RepoRoot()
	if err != nil {
		fatal(err)
	}
	cacheDir := filepath.Join(root, "pipeline", ".fixtures", "bin")

	var scales []fixtures.Scale
	for _, s := range strings.Split(*scalesFlag, ",") {
		s = strings.TrimSpace(s)
		if s == "" {
			continue
		}
		sc, err := fixtures.ParseScale(s)
		if err != nil {
			fatal(err)
		}
		scales = append(scales, sc)
	}
	if len(scales) == 0 {
		fatal(fmt.Errorf("no scales given"))
	}

	allowed := map[string]bool{}
	for _, t := range strings.Split(*toolsFlag, ",") {
		if t = strings.TrimSpace(t); t != "" {
			allowed[t] = true
		}
	}

	overallOK := true
	for _, scale := range scales {
		ok := runScale(root, cacheDir, scale, allowed, *seedFlag, *repsFlag, *updateFix, *skipBench, *outFlag)
		overallOK = overallOK && ok
	}

	if !overallOK {
		fmt.Fprintln(os.Stderr, "\nFULL VALIDATION: FAIL — see the per-scale reports above.")
		os.Exit(1)
	}
	fmt.Fprintln(os.Stderr, "\nFULL VALIDATION: PASS")
}

func runScale(root, cacheDir string, scale fixtures.Scale, allowed map[string]bool, seed int64, reps int, force, skipBench bool, outRoot string) bool {
	logf := func(f string, a ...any) { fmt.Fprintf(os.Stderr, f+"\n", a...) }
	logf("\n======== scale=%s ========", scale)

	// --- Fixtures ------------------------------------------------------
	t0 := time.Now()
	man, err := fixtures.Generate(fixtures.Options{Scale: scale, Seed: seed, Force: force, Logf: logf})
	if err != nil {
		fatal(fmt.Errorf("fixtures %s: %w", scale, err))
	}
	logf("fixtures ready in %s", time.Since(t0).Round(time.Millisecond))

	// --- 1. Parity matrix ---------------------------------------------
	entries := matrix.Default().FilterTools(allowed)
	logf("[parity] running %d matrix entries...", len(entries))
	t1 := time.Now()
	cfg := runner.Config{Manifest: man, CacheDir: cacheDir, Logf: func(string, ...any) {}}
	var results []runner.Result
	var diverge, errors, similar int
	for _, e := range entries {
		res := runner.RunEntry(cfg, e)
		results = append(results, res)
		switch res.Status {
		case runner.StatusDiverge:
			diverge++
			logf("  DIVERGE %s/%s: %s", res.Tool, res.Name, res.Detail)
		case runner.StatusError:
			errors++
			logf("  ERROR   %s/%s: %s", res.Tool, res.Name, res.Detail)
		case runner.StatusSimilar:
			similar++
		}
	}
	report := runner.Build(string(scale), seed, results)
	parityDur := time.Since(t1)

	// --- 2. Round-trip -------------------------------------------------
	logf("[round-trip] format encode/decode checks...")
	t2 := time.Now()
	rt := roundtrip.Run(man, cacheDir)
	rtFail := 0
	for _, r := range rt {
		if r.Status == roundtrip.Fail {
			rtFail++
			logf("  FAIL %s (%s): %s", r.Name, r.Format, r.Detail)
		}
	}
	rtDur := time.Since(t2)

	// --- 3. Performance / RSS bench -----------------------------------
	var cells []bench.CellResult
	var benchDur time.Duration
	if !skipBench {
		logf("[bench] wall/CPU/RSS sweep (reps=%d)...", reps)
		t3 := time.Now()
		cells, err = bench.Run(bench.RunConfig{Scales: []fixtures.Scale{scale}, Reps: reps, Seed: seed, Log: func(string, ...any) {}})
		if err != nil {
			logf("  bench error: %v", err)
		}
		benchDur = time.Since(t3)
	}

	// --- Consolidated report ------------------------------------------
	outDir := outRoot
	if outDir == "" {
		outDir = filepath.Join(root, "pipeline", ".fixtures", string(scale), "full-validation")
	}
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		fatal(err)
	}
	_ = writeReports(outDir, scale, report, rt, cells)

	// --- Verdict + resource summary -----------------------------------
	var mem runtime.MemStats
	runtime.ReadMemStats(&mem)
	pass := diverge == 0 && errors == 0 && rtFail == 0
	logf("\n  scale=%s verdict: %s", scale, verdict(pass))
	logf("    parity      : %d entries — %d pass, %d similar, %d diverge, %d error  [%s]",
		report.Summary.Total, report.Summary.Pass, report.Summary.Similar, report.Summary.Diverge, report.Summary.Error, parityDur.Round(time.Millisecond))
	logf("    round-trip  : %d checks — %d fail  [%s]", len(rt), rtFail, rtDur.Round(time.Millisecond))
	if !skipBench {
		logf("    bench       : %d cells captured (wall/CPU/RSS)  [%s]", len(cells), benchDur.Round(time.Millisecond))
	}
	logf("    orchestrator peak heap: %d MB", mem.Sys/(1<<20))
	logf("    reports     : %s/{report.json,report.md,roundtrip.md,bench.md}", outDir)
	return pass
}

func writeReports(outDir string, scale fixtures.Scale, report runner.Report, rt []roundtrip.Result, cells []bench.CellResult) error {
	// Parity (reuse the runner's own JSON/MD writers).
	_ = report.WriteJSON(filepath.Join(outDir, "report.json"))
	_ = report.WriteMarkdown(filepath.Join(outDir, "report.md"))
	// Round-trip.
	var b strings.Builder
	fmt.Fprintf(&b, "# Round-trip validation — scale %s\n\n| check | format | status | detail |\n|---|---|---|---|\n", scale)
	for _, r := range rt {
		fmt.Fprintf(&b, "| %s | %s | %s | %s |\n", r.Name, r.Format, r.Status, r.Detail)
	}
	_ = os.WriteFile(filepath.Join(outDir, "roundtrip.md"), []byte(b.String()), 0o644)
	// Bench.
	if len(cells) > 0 {
		_ = bench.WriteJSON(filepath.Join(outDir, "bench.json"), cells)
		_ = bench.WriteMarkdown(filepath.Join(outDir, "bench.md"), cells, []fixtures.Scale{scale})
	}
	return nil
}

func verdict(pass bool) string {
	if pass {
		return "PASS"
	}
	return "FAIL"
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, "full-validation:", err)
	os.Exit(2)
}
