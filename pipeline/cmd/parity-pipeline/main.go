// Command parity-pipeline runs the combinatorics parity matrix: it generates
// (or reuses) real-sized fixtures at a chosen scale, runs every matrix entry
// with OUR tool binary and the vendored UPSTREAM binary, compares their output,
// and writes a machine-readable report.json plus a human report.md.
//
// Usage:
//
//	go run ./pipeline/cmd/parity-pipeline -scale=smoke
//	go run ./pipeline/cmd/parity-pipeline -scale=small -tools=samtools,bcftools
//	go run ./pipeline/cmd/parity-pipeline -update-fixtures -scale=medium
//
// Exit status is non-zero if any entry DIVERGES or ERRORS (SIMILAR matches and
// SKIPs do not fail the run). See pipeline/README.md.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

func main() {
	var (
		scaleFlag   = flag.String("scale", "", "fixture scale: smoke|small|medium|large (default: $PIPELINE_SCALE or small)")
		toolsFlag   = flag.String("tools", "", "comma-separated tool filter (default: all)")
		benchFlag   = flag.Bool("bench", false, "informational: point to the bench harness (see -h)")
		updateFix   = flag.Bool("update-fixtures", false, "force regeneration of fixtures even if cached")
		outFlag     = flag.String("out", "", "output directory for report.json/report.md (default: pipeline/.fixtures/<scale>/report)")
		seedFlag    = flag.Int64("seed", fixtures.DefaultSeed, "RNG seed for fixture generation")
		verboseFlag = flag.Bool("v", false, "verbose logging")
	)
	flag.Parse()

	if *benchFlag {
		fmt.Fprintln(os.Stderr, "Performance benchmarks live in pipeline/bench; run them with:\n"+
			"  PIPELINE_SCALE=medium go test -bench=. ./pipeline/bench")
	}

	scale, err := fixtures.ParseScale(*scaleFlag)
	if err != nil {
		fatal(err)
	}

	logf := func(format string, args ...any) {}
	if *verboseFlag {
		logf = func(format string, args ...any) { fmt.Fprintf(os.Stderr, format+"\n", args...) }
	}

	root, err := upstream.RepoRoot()
	if err != nil {
		fatal(err)
	}

	// Generate / reuse fixtures.
	man, err := fixtures.Generate(fixtures.Options{
		Scale: scale,
		Seed:  *seedFlag,
		Force: *updateFix,
		Logf:  func(f string, a ...any) { fmt.Fprintf(os.Stderr, f+"\n", a...) },
	})
	if err != nil {
		fatal(err)
	}

	// Filter the matrix by tool if requested.
	allowed := map[string]bool{}
	for _, t := range strings.Split(*toolsFlag, ",") {
		t = strings.TrimSpace(t)
		if t != "" {
			allowed[t] = true
		}
	}
	entries := matrix.Default().FilterTools(allowed)
	if len(entries) == 0 {
		fatal(fmt.Errorf("no matrix entries match -tools=%q", *toolsFlag))
	}

	cacheDir := filepath.Join(root, "pipeline", ".fixtures", "bin")
	cfg := runner.Config{Manifest: man, CacheDir: cacheDir, Logf: logf}

	fmt.Fprintf(os.Stderr, "running %d matrix entries at scale=%s...\n", len(entries), scale)
	var results []runner.Result
	for _, e := range entries {
		res := runner.RunEntry(cfg, e)
		results = append(results, res)
		fmt.Fprintf(os.Stderr, "  [%-7s] %s/%s\n", res.Status, res.Tool, res.Name)
	}

	report := runner.Build(string(scale), *seedFlag, results)

	outDir := *outFlag
	if outDir == "" {
		outDir = filepath.Join(root, "pipeline", ".fixtures", string(scale), "report")
	}
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		fatal(err)
	}
	jsonPath := filepath.Join(outDir, "report.json")
	mdPath := filepath.Join(outDir, "report.md")
	if err := report.WriteJSON(jsonPath); err != nil {
		fatal(err)
	}
	if err := report.WriteMarkdown(mdPath); err != nil {
		fatal(err)
	}

	s := report.Summary
	fmt.Fprintf(os.Stderr, "\nsummary: total=%d PASS=%d SIMILAR=%d DIVERGE=%d SKIP=%d ERROR=%d\n",
		s.Total, s.Pass, s.Similar, s.Diverge, s.Skip, s.Error)
	fmt.Fprintf(os.Stderr, "reports: %s\n         %s\n", jsonPath, mdPath)

	if report.Failed() {
		os.Exit(1)
	}
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, "parity-pipeline:", err)
	os.Exit(2)
}
