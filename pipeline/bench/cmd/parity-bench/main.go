// Command parity-bench runs the performance & scalability matrix: it times OUR
// tool binaries against the vendored UPSTREAM binaries over the small / medium /
// large fixture tiers, recording wall-clock, CPU time (user+sys) and peak RSS
// for each, and emits a JSON + Markdown report. Sweeping the tiers turns the
// point measurements into scalability curves versus input size (read count,
// coverage, variant count, interval count grow per tier — see
// pipeline/fixtures/scale.go).
//
// Examples:
//
//	# the headline sweep for the manuscript (heavy; run on a beefy box / HPC node)
//	go run ./pipeline/bench/cmd/parity-bench -scales medium,large -reps 5
//
//	# a quick local sanity pass
//	go run ./pipeline/bench/cmd/parity-bench -scales small -reps 2
//
//	# one file-format group, one tier
//	go run ./pipeline/bench/cmd/parity-bench -scales medium -group BED
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/bench"
	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
)

func main() {
	scalesFlag := flag.String("scales", "small", "comma-separated fixture tiers to sweep: smoke|small|medium|large")
	reps := flag.Int("reps", 3, "repetitions per side (wall/CPU take the min, RSS the max)")
	seed := flag.Int64("seed", 1, "fixture generation seed")
	cellGlob := flag.String("cell", "", "substring filter on cell name (default: all)")
	group := flag.String("group", "", "file-format group filter: BAM|CRAM|VCF|BED|FASTQ (default: all)")
	out := flag.String("out", "", "report output dir (default: pipeline/.fixtures/<lastScale>/bench)")
	flag.Parse()

	var scales []fixtures.Scale
	for _, s := range strings.Split(*scalesFlag, ",") {
		s = strings.TrimSpace(s)
		if s == "" {
			continue
		}
		sc, err := fixtures.ParseScale(s)
		if err != nil {
			log.Fatalf("bad -scales: %v", err)
		}
		scales = append(scales, sc)
	}
	if len(scales) == 0 {
		log.Fatal("no scales given")
	}

	results, err := bench.Run(bench.RunConfig{
		Scales:   scales,
		Reps:     *reps,
		Seed:     *seed,
		CellGlob: *cellGlob,
		GroupSel: *group,
		Log:      func(f string, a ...any) { fmt.Printf(f+"\n", a...) },
	})
	if err != nil {
		log.Fatalf("bench: %v", err)
	}

	outDir := *out
	if outDir == "" {
		outDir = filepath.Join("pipeline", ".fixtures", string(scales[len(scales)-1]), "bench")
	}
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		log.Fatalf("mkdir %s: %v", outDir, err)
	}
	jsonPath := filepath.Join(outDir, "bench.json")
	mdPath := filepath.Join(outDir, "bench.md")
	if err := bench.WriteJSON(jsonPath, results); err != nil {
		log.Fatalf("write json: %v", err)
	}
	if err := bench.WriteMarkdown(mdPath, results, scales); err != nil {
		log.Fatalf("write md: %v", err)
	}

	// One-line summary so a CI/cluster log shows the headline at a glance.
	var ok, errs int
	for _, r := range results {
		if r.Err != "" {
			errs++
		} else {
			ok++
		}
	}
	fmt.Printf("\nbench done: %d ok, %d errored across %d tier(s) → %s\n", ok, errs, len(scales), outDir)
}
