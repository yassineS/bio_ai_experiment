// bedshuffle randomly relocates BED intervals across a genome, mirroring
// `bedtools shuffle`.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedshuffle/pkg/bedshuffle"
)

const usage = `bedshuffle - Randomly relocate BED intervals across a genome

Usage:
  bedshuffle [options] -i <input.bed> -g <genome.txt>

Description:
  For each input interval, sample a new (chrom, start) coordinate uniformly
  at random from the genome (weighted by chromosome length), preserving the
  interval length and any extra columns. Several flags constrain the draw:

    -incl   restrict placements to regions in an include BED
    -excl   forbid placements that overlap regions in an exclude BED
    -chrom  keep each shuffled interval on its original chromosome

Options:
  -i,   --input FILE       Input BED (required; '-' for stdin)
  -g,   --genome FILE      Chrom-size table (tab-separated chrom<TAB>size)
        --output FILE      Output file (default: stdout)
  -incl, --include FILE    Include-region BED
  -excl, --exclude FILE    Exclude-region BED
  -chrom,--chromOnly       Keep each interval on its original chromosome
  -seed, --seed N          Deterministic seed (default: 0)
  -maxTries N              Placement retries per interval (default: 1000)
  -h,    --help            Show this help
  -v,    --version         Show version

Examples:
  # Basic shuffle on hg19 chrom sizes:
  bedshuffle -i input.bed -g hg19.genome > shuffled.bed

  # Restrict placements to regions in incl.bed:
  bedshuffle -i input.bed -g hg19.genome -incl incl.bed > shuffled.bed

  # Avoid centromeres (exclude file):
  bedshuffle -i input.bed -g hg19.genome -excl centromeres.bed > shuffled.bed

  # Keep intervals on their original chromosome:
  bedshuffle -i input.bed -g hg19.genome -chrom > shuffled.bed
`

const version = "bedshuffle 0.1.0"

func main() {
	fs := flag.CommandLine

	var inputFile, genomeFile, outputFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED")
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Chrom-size file")
	cliflag.StringVar(fs, &outputFile, "", "output", "", "Output (default stdout)")

	var inclFile, exclFile string
	cliflag.StringVar(fs, &inclFile, "", "incl", "", "Include BED")
	cliflag.StringVar(fs, &exclFile, "", "excl", "", "Exclude BED")
	// Also accept the long-form names from bedtools shuffle.
	cliflag.StringVar(fs, &inclFile, "", "include", "", "Include BED (alias of -incl)")
	cliflag.StringVar(fs, &exclFile, "", "exclude", "", "Exclude BED (alias of -excl)")

	var keepChrom bool
	cliflag.BoolVar(fs, &keepChrom, "", "chrom", false, "Keep on original chrom")
	cliflag.BoolVar(fs, &keepChrom, "", "chromOnly", false, "Keep on original chrom (alias of -chrom)")

	var seed int
	cliflag.IntVar(fs, &seed, "", "seed", 0, "RNG seed")

	var maxTries int
	cliflag.IntVar(fs, &maxTries, "", "maxTries", 1000, "Placement retries per interval")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version")

	flag.Parse()
	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Println(version)
		os.Exit(0)
	}

	if inputFile == "" || genomeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -i and -g are required")
		fmt.Fprintln(os.Stderr, "Use -h for help.")
		os.Exit(1)
	}

	gf, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening genome: %v\n", err)
		os.Exit(1)
	}
	defer gf.Close()
	genome, err := bedshuffle.ParseGenome(gf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading genome: %v\n", err)
		os.Exit(1)
	}

	opts := bedshuffle.Options{
		Genome:     genome,
		Chrom:      keepChrom,
		Seed:       int64(seed),
		MaxRetries: maxTries,
	}

	if inclFile != "" {
		f, err := iohelper.OpenReader(inclFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening include: %v\n", err)
			os.Exit(1)
		}
		recs, err := bedshuffle.ParseBED(f)
		f.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading include: %v\n", err)
			os.Exit(1)
		}
		opts.Include = recs
	}
	if exclFile != "" {
		f, err := iohelper.OpenReader(exclFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening exclude: %v\n", err)
			os.Exit(1)
		}
		recs, err := bedshuffle.ParseBED(f)
		f.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading exclude: %v\n", err)
			os.Exit(1)
		}
		opts.Exclude = recs
	}

	in, err := iohelper.OpenReader(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()
	out, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedshuffle.Shuffle(in, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}
}
