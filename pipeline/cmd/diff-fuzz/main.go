// Command diff-fuzz is the CLI for the differential-fuzzing harness
// (pipeline/difffuzz). It generates fuzzed inputs for a set of tool/subcommand
// targets, runs BOTH our binary and the vendored upstream binary on each,
// diffs stdout/stderr/exit-code (after the parity harness's provenance
// normalization), classifies and minimizes divergences, optionally captures Go
// coverage of our binaries, and writes difffuzz.{json,md}.
//
// Usage:
//
//	# quick local smoke (seconds, 3 targets) against the vendored binaries
//	go run ./pipeline/cmd/diff-fuzz -quick
//
//	# a fuller run with coverage capture
//	go run ./pipeline/cmd/diff-fuzz -iters=2000 -coverage -targets=bcftools-view,bedtools-merge
//
// Exit status is 1 if any STRICT parity divergence is found (stdout/stderr/
// exit-code/one-crashed), 0 otherwise (both-crashed and SKIPs do not fail).
// See docs/DIFFERENTIAL_FUZZING.md.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/difffuzz"
	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

func main() {
	var (
		quick     = flag.Bool("quick", false, "quick mode: few iterations on a small target subset (runs in seconds)")
		iters     = flag.Int("iters", 0, "iterations per target (default: 200, or 40 in -quick)")
		seed      = flag.Int64("seed", 1, "RNG seed (reproducible)")
		targets   = flag.String("targets", "", "comma-separated target-name filter (default: all default targets)")
		coverage  = flag.Bool("coverage", false, "capture Go statement coverage of our binaries over the run")
		timeout   = flag.Duration("timeout", 10*time.Second, "per-invocation timeout")
		scaleFlag = flag.String("scale", "small", "fixture scale for seed inputs: smoke|small|medium")
		outFlag   = flag.String("out", "", "output dir for difffuzz.{json,md} (default: pipeline/.fixtures/<scale>/difffuzz)")
		minSteps  = flag.Int("minimize-steps", 500, "max delta-debugging steps per reproducer")
		verbose   = flag.Bool("v", false, "verbose progress logging")
	)
	flag.Parse()

	iterations := *iters
	if iterations == 0 {
		if *quick {
			iterations = 40
		} else {
			iterations = 200
		}
	}

	root, err := upstream.RepoRoot()
	if err != nil {
		fatal(err)
	}

	scale, err := fixtures.ParseScale(*scaleFlag)
	if err != nil {
		fatal(err)
	}

	logf := func(string, ...any) {}
	if *verbose {
		logf = func(f string, a ...any) { fmt.Fprintf(os.Stderr, f+"\n", a...) }
	}

	// Generate / reuse seed fixtures (the mutation strategy perturbs these).
	man, err := fixtures.Generate(fixtures.Options{
		Scale: scale,
		Seed:  fixtures.DefaultSeed,
		Logf:  func(f string, a ...any) { fmt.Fprintf(os.Stderr, f+"\n", a...) },
	})
	if err != nil {
		fatal(fmt.Errorf("generating fixtures (needed for seed inputs): %w", err))
	}

	// Select targets.
	var tgts []difffuzz.Target
	if *quick && *targets == "" {
		tgts = difffuzz.QuickTargets()
	} else {
		tgts = difffuzz.DefaultTargets()
	}
	if *targets != "" {
		allow := map[string]bool{}
		for _, t := range strings.Split(*targets, ",") {
			if t = strings.TrimSpace(t); t != "" {
				allow[t] = true
			}
		}
		var filtered []difffuzz.Target
		for _, t := range tgts {
			if allow[t.Name] {
				filtered = append(filtered, t)
			}
		}
		if len(filtered) == 0 {
			fatal(fmt.Errorf("no targets match -targets=%q (available: %s)", *targets, targetNames()))
		}
		tgts = filtered
	}

	cacheDir := filepath.Join(root, "pipeline", ".fixtures", "bin")
	workDir := filepath.Join(root, "pipeline", ".fixtures", string(scale), "difffuzz-work")

	cfg := difffuzz.Config{
		Targets:       tgts,
		Iterations:    iterations,
		Seed:          *seed,
		Timeout:       *timeout,
		CacheDir:      cacheDir,
		WorkDir:       workDir,
		Manifest:      man,
		Coverage:      *coverage,
		MinimizeSteps: *minSteps,
		Logf:          logf,
	}

	fmt.Fprintf(os.Stderr, "diff-fuzz: %d target(s), %d iters each, seed=%d, coverage=%v\n",
		len(tgts), iterations, *seed, *coverage)

	rep, err := difffuzz.Run(cfg)
	if err != nil {
		fatal(err)
	}
	defer os.RemoveAll(workDir)

	outDir := *outFlag
	if outDir == "" {
		outDir = filepath.Join(root, "pipeline", ".fixtures", string(scale), "difffuzz")
	}
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		fatal(err)
	}
	jsonPath := filepath.Join(outDir, "difffuzz.json")
	mdPath := filepath.Join(outDir, "difffuzz.md")
	if err := rep.WriteJSON(jsonPath); err != nil {
		fatal(err)
	}
	if err := rep.WriteMarkdown(mdPath); err != nil {
		fatal(err)
	}

	// Console summary.
	fmt.Fprintln(os.Stderr)
	for _, t := range rep.Targets {
		if t.Skipped {
			fmt.Fprintf(os.Stderr, "  [SKIP] %-22s %s\n", t.Name, firstLine(t.SkipReason))
			continue
		}
		div := 0
		for c, n := range t.ByClass {
			if c != difffuzz.ClassNone {
				div += n
			}
		}
		cov := ""
		if t.Coverage.Enabled {
			cov = fmt.Sprintf("  cov=%.1f%%", t.Coverage.Percent)
		}
		fmt.Fprintf(os.Stderr, "  [%-4d] %-22s divergences=%d%s\n", t.Inputs, t.Name, div, cov)
	}
	fmt.Fprintf(os.Stderr, "\nreports: %s\n         %s\n", jsonPath, mdPath)

	if rep.HasParityDivergence() {
		fmt.Fprintln(os.Stderr, "\nFAIL: strict parity divergence(s) found (see report).")
		os.Exit(1)
	}
}

func targetNames() string {
	var names []string
	for _, t := range difffuzz.DefaultTargets() {
		names = append(names, t.Name)
	}
	return strings.Join(names, ", ")
}

func firstLine(s string) string {
	if i := strings.IndexByte(s, '\n'); i >= 0 {
		return s[:i]
	}
	return s
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, "diff-fuzz:", err)
	os.Exit(2)
}
