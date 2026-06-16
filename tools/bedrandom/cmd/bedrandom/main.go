// bedrandom generates random fixed-length intervals across a genome and writes
// them as BED6 records (mirrors `bedtools random`).
package main

import (
	"flag"
	"fmt"
	"os"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedrandom/pkg/bedrandom"
)

const version = "1.0.0"

const usage = `bedrandom - Generate random intervals across a genome

Usage:
  bedrandom [options] -g <genome>

Description:
  Generates -n random intervals, each of length -l, placed uniformly at random
  across the supplied genome. Each interval is emitted as a BED6 record:

      chrom  start  end  index  length  strand

  where index counts 1..n and the score column holds the interval length. The
  strand is chosen at random (+ or -). With -seed the output is fully
  reproducible and byte-for-byte identical to upstream 'bedtools random'.

Options:
  -l, --length NUM         Length of each interval (default: 100)
  -n, --number NUM         Number of intervals to generate (default: 1000000)
  -seed NUM                Seed for the random number generator. If unset, a
                           seed is derived from the current time and pid (the
                           result is then non-reproducible).
  -g, --genome FILE        Genome (chrom-sizes) file: 'chrom<TAB>size' per
                           line. Required.
  -o, --output FILE        Output file ('-' for stdout, default: stdout)
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # 100 intervals of length 1000, reproducible
  bedrandom -g hg38.sizes -l 1000 -n 100 -seed 42 > random.bed

  # default: 1,000,000 intervals of length 100
  bedrandom -g hg38.sizes > random.bed
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var genomeFile, outputFile string
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Genome (chrom-sizes) file (required)")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output file")

	var length, number int
	cliflag.IntVar(fs, &length, "l", "length", bedrandom.DefaultLength, "Length of each interval")
	cliflag.IntVar(fs, &number, "n", "number", bedrandom.DefaultNum, "Number of intervals")

	// "seed" is a single-dash long flag in bedtools; Go's flag package treats
	// -seed and --seed identically. Use a sentinel to detect whether it was set.
	const seedUnset = int(^uint(0) >> 1) // max int sentinel
	seed := seedUnset
	fs.IntVar(&seed, "seed", seedUnset, "Seed for the RNG")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedrandom version %s\n", version)
		os.Exit(0)
	}

	if genomeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -g/--genome is required")
		os.Exit(2)
	}

	opts := bedrandom.Options{Length: length, Num: number}
	if seed != seedUnset {
		opts.Seed = seed
		opts.HaveSeed = true
	} else {
		// Match upstream: seed = time(0) + getpid().
		opts.Seed = int(time.Now().Unix()) + os.Getpid()
		opts.HaveSeed = false
	}

	gr, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening genome file: %v\n", err)
		os.Exit(1)
	}
	genome, err := bedrandom.ParseGenome(gr)
	gr.Close()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
		os.Exit(1)
	}

	out, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedrandom.Generate(genome, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
