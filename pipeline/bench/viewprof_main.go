//go:build ignore

// viewprof_main.go is a throwaway CPU-profiling harness for `samtools view`
// region/BED queries. Run with:
//
//	go run pipeline/bench/viewprof_main.go -bam <bam> -bed <bed> -cpuprofile out.pprof
//
// It drives samtools.ViewFile the same way the CLI does and writes SAM to
// /dev/null so the profile reflects the read/filter/serialise hot path.
package main

import (
	"flag"
	"io"
	"log"
	"os"
	"runtime/pprof"

	"github.com/yassineS/bio_ai_experiment/tools/samtools/pkg/samtools"
)

func main() {
	bam := flag.String("bam", "", "input BAM")
	bed := flag.String("bed", "", "BED file for -L")
	region := flag.String("region", "", "region spec chr:start-end")
	cpu := flag.String("cpuprofile", "", "cpu profile output")
	reps := flag.Int("reps", 1, "repetitions")
	flag.Parse()

	if *cpu != "" {
		f, err := os.Create(*cpu)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	for i := 0; i < *reps; i++ {
		opts := samtools.ViewOptions{}
		if *bed != "" {
			opts.BedPath = *bed
		}
		if *region != "" {
			opts.Regions = []string{*region}
		}
		n, err := samtools.ViewFile(*bam, io.Discard, opts, os.Stderr)
		if err != nil {
			log.Fatal(err)
		}
		_ = n
	}
}
