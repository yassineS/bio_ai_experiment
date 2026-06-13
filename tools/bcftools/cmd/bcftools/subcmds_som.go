// CLI runner for `bcftools som` (Self-Organizing Map variant classifier).
// It follows the same parse-validate-dispatch shape as the other runners.
// Upstream's getopt list (vcfsom.c:677) is reproduced here, with the
// train→classify pipeline fixed so the `.som` map is actually usable (the
// upstream som_write_map fwrite-return bug is corrected in the library).
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const somUsage = `bcftools som - Self-Organizing Map (Kohonen map) variant classifier.

Usage:
  bcftools som --train    [options] <in.vcf[.gz]|in.bcf>
  bcftools som --classify [options] <in.vcf[.gz]|in.bcf>

Modes:
  -t, --train                    Train a map from the input VCF's INFO
                                 annotations and write <prefix>.som.
  -c, --classify                 Load <prefix>.som and print one SOM score
                                 per site (1 = most training-like).

Options:
  -p, --prefix STRING            Prefix for the map file (<prefix>.som). Required.
  -T, --training-annots LIST     Comma list of INFO tags forming the per-site
                                 vector (QUAL is also accepted). Default:
                                 QUAL,MQ,MQ0F,BQB,MQB,RPB,SGB.
  -s, --size INT                 Map edge length (nbin) [20]. The 2-D map has
                                 size² nodes.
  -d, --som-dimension INT        Map dimensionality [2]. Upstream requires >=2.
  -n, --ntrain-sites INT         Effective number of training iterations used
                                 in the learning-rate decay [number of sites].
  -l, --learning-rate FLOAT      Learning rate [1.0].
  -b, --bmu-threshold FLOAT      Count threshold selecting classification
                                 nodes [0.9].
  -r, --random-seed INT          RNG seed for weight init [1].
  -o, --output PATH              Write classify scores here (default stdout).
  -?, --help                     Show this help.
      --version                  Show version.

Notes:
  - Unlike upstream (which reads a pre-extracted annots.tab.gz), this port
    reads a VCF/BCF directly and extracts the annotation vector from INFO.
  - The on-disk map format is our own clean, versioned binary format
    (magic "SOMGO1"); upstream's is unusable due to a write bug. See
    docs/UPSTREAM_BUGS.md#bcftools-som-write-map.
  - Upstream's experimental -f/--nfold cross-validation, -m/--merge,
    and -e/--exclude-bad knobs are accepted-as-no-op surface only (v1
    trains a single map).
`

func runSom(args []string) int {
	fs := flag.NewFlagSet("bcftools som", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		train          bool
		classify       bool
		prefix         string
		trainingAnnots string
		size           int
		somDimension   int
		ntrainSites    int
		learningRate   float64
		bmuThreshold   float64
		randomSeed     int
		outputPath     string
		// Accepted upstream surface that v1 treats as a single-map no-op.
		nfold      int
		mergeAlg   string
		excludeBad bool
		showHelp   bool
		showVer    bool
	)

	cliflag.BoolVar(fs, &train, "t", "train", false, "Train a map")
	cliflag.BoolVar(fs, &classify, "c", "classify", false, "Classify with a map")
	cliflag.StringVar(fs, &prefix, "p", "prefix", "", "Map file prefix")
	cliflag.StringVar(fs, &trainingAnnots, "T", "training-annots", "", "INFO annotation list")
	cliflag.IntVar(fs, &size, "s", "size", bcftools.SomDefaultSize, "Map edge length")
	cliflag.IntVar(fs, &somDimension, "d", "som-dimension", bcftools.SomDefaultDimension, "Map dimensionality")
	cliflag.IntVar(fs, &ntrainSites, "n", "ntrain-sites", 0, "Effective number of training iterations")
	cliflag.Float64Var(fs, &learningRate, "l", "learning-rate", bcftools.SomDefaultLearn, "Learning rate")
	cliflag.Float64Var(fs, &bmuThreshold, "b", "bmu-threshold", bcftools.SomDefaultBmuThreshold, "Best-matching-unit count threshold")
	cliflag.IntVar(fs, &randomSeed, "r", "random-seed", bcftools.SomDefaultRandomSeed, "RNG seed")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path for classify scores")
	// Accepted-but-no-op upstream surface (single-map v1).
	cliflag.IntVar(fs, &nfold, "f", "nfold", 1, "n-fold cross-validation (accepted; v1 trains a single map)")
	cliflag.StringVar(fs, &mergeAlg, "m", "merge", "avg", "-f merge algorithm (accepted; v1 trains a single map)")
	cliflag.BoolVar(fs, &excludeBad, "e", "exclude-bad", false, "Exclude bad sites from training (accepted; v1 trains on all sites)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Surface-only flags are intentionally not threaded into the library
	// (v1 trains a single map on all sites).
	_ = nfold
	_ = mergeAlg
	_ = excludeBad

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, somUsage)
		return 2
	}
	if showHelp {
		fmt.Print(somUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	if train == classify {
		// Both unset or both set: upstream requires exactly one mode.
		fmt.Fprintln(os.Stderr, "bcftools som: exactly one of --train or --classify is required")
		fmt.Fprint(os.Stderr, somUsage)
		return 2
	}
	if prefix == "" {
		fmt.Fprintln(os.Stderr, "bcftools som: -p/--prefix is required")
		return 2
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools som: missing input VCF/BCF")
		fmt.Fprint(os.Stderr, somUsage)
		return 2
	}
	input := rest[0]

	opts := bcftools.SomOptions{
		Prefix:       prefix,
		Size:         size,
		NDim:         somDimension,
		NTrain:       ntrainSites,
		Learn:        learningRate,
		BmuThreshold: bmuThreshold,
		RandomSeed:   int64(randomSeed),
	}
	if trainingAnnots != "" {
		opts.TrainingAnnots = bcftools.SplitCommaList(trainingAnnots)
	}

	if train {
		opts.Action = bcftools.SomActionTrain
		res, err := bcftools.SomTrain(input, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools som: %v\n", err)
			return 1
		}
		fmt.Fprintf(os.Stderr, "bcftools som: trained %d-node map on %d sites (%d annots) -> %s.som\n",
			res.MapSize, res.NSites, res.KDim, prefix)
		return 0
	}

	// classify
	opts.Action = bcftools.SomActionClassify
	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools som: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.SomClassify(input, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools som: %v\n", err)
		return 1
	}
	return 0
}
