// Command realparity runs OUR bioinformatics ports against the UPSTREAM
// binaries on real, multi-contig (whole-genome) input files — a GIAB-class
// BAM/CRAM/VCF plus an indexed reference FASTA — over a battery of
// representative samtools/bcftools commands. For each cell it reports BOTH:
//
//   - parity: byte-exact-after-provenance-stripping equality of the two outputs
//     (the project's exact parity definition, runner.CompareByteExact /
//     runner.StripProvenance — the SAME notion used everywhere else in the repo),
//     and
//   - performance: wall-clock, CPU (user+sys), and peak RSS for each side, plus
//     the ours/upstream ratio.
//
// This is the manuscript's "real-world differential parity + performance on
// whole-genome data" experiment (claims C2/C3), replacing the synthetic-only
// large tier with real multi-contig inputs. It is purely our-vs-upstream
// differential testing: no truth set, hap.py, or vcfeval is needed.
//
// Every input is optional. A cell whose required input is absent SKIPs (with a
// clear reason) rather than failing, so the command runs partially with whatever
// data is present and runs unchanged — exiting 0 with an all-SKIP report — on a
// machine with no real data (the common CI case). The heavy run happens on a
// machine that has GIAB-class data.
//
// Usage:
//
//	realparity -bam HG002.bam -vcf HG002.vcf.gz -ref GRCh38.fa -out ./reports
//	realparity -bam x.bam -region chr20 -reps 5 -v
//
// Cross-contig behaviour (BAI/CSI multi-ref bins, RNEXT '=' vs mate-on-other-
// contig, coordinate sort across contigs, per-contig idxstats counts) is exactly
// what the multi-contig input exercises; `view` of the whole file and `idxstats`
// are in the battery on purpose. The parity gate is zero DIVERGE.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
)

func main() {
	os.Exit(run(os.Args[1:]))
}

func run(argv []string) int {
	fs := flag.NewFlagSet("realparity", flag.ContinueOnError)
	var (
		ourBin  = fs.String("our-bin", "", "dir of our built binaries (samtools, bcftools); if empty, build them into a temp dir from tools/*/cmd")
		upBin   = fs.String("upstream-bin", "", "dir with upstream samtools/bcftools; default resolves the vendored reference_code/ build locations")
		ref     = fs.String("ref", "", "indexed reference FASTA (.fa with .fai)")
		bam     = fs.String("bam", "", "multi-contig BAM (indexed)")
		vcf     = fs.String("vcf", "", "bgzipped + indexed VCF (.vcf.gz with .tbi/.csi)")
		region  = fs.String("region", "", "optional region (e.g. chr20 or chr20:1-1000000); empty = whole-genome, all contigs")
		reps    = fs.Int("reps", 3, "measurement repetitions (wall/CPU = min, RSS = max)")
		outDir  = fs.String("out", "", "report output directory (default: current directory)")
		verbose = fs.Bool("v", false, "verbose progress logging to stderr")
		showVer = fs.Bool("version", false, "print version and exit")
	)
	fs.BoolVar(showVer, "V", false, "print version and exit")
	if err := fs.Parse(argv); err != nil {
		return 2
	}
	if *showVer {
		fmt.Println("realparity (bio_ai_experiment) v0.1.0")
		return 0
	}
	if *reps < 1 {
		*reps = 1
	}

	tmpDir, err := os.MkdirTemp("", "realparity-")
	if err != nil {
		fmt.Fprintln(os.Stderr, "realparity: creating temp dir:", err)
		return 1
	}
	defer os.RemoveAll(tmpDir)
	ourCache := filepath.Join(tmpDir, "ourbins")

	bins, notes, err := resolveBins(*ourBin, *upBin, ourCache)
	if err != nil {
		fmt.Fprintln(os.Stderr, "realparity: resolving binaries:", err)
		return 1
	}
	if *verbose {
		for _, n := range notes {
			fmt.Fprintln(os.Stderr, "realparity: note:", n)
		}
		fmt.Fprintf(os.Stderr, "realparity: ours samtools=%q bcftools=%q\n", bins.oursSamtools, bins.oursBcftools)
		fmt.Fprintf(os.Stderr, "realparity: upstream samtools=%q bcftools=%q\n", bins.upSamtools, bins.upBcftools)
	}

	cfg := config{
		bins:    bins,
		in:      inputs{ref: abs(*ref), bam: abs(*bam), vcf: abs(*vcf)},
		region:  *region,
		reps:    *reps,
		outDir:  *outDir,
		tmpDir:  tmpDir,
		verbose: *verbose,
	}

	rep, err := runBattery(cfg)
	if err != nil {
		fmt.Fprintln(os.Stderr, "realparity:", err)
		return 1
	}

	jsonPath, mdPath, err := writeReports(rep, *outDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "realparity: writing reports:", err)
		return 1
	}

	fmt.Printf("realparity verdict: %s (PASS=%d DIVERGE=%d SKIP=%d ERROR=%d)\n",
		rep.verdict(), rep.Pass, rep.Diverge, rep.Skip, rep.Errored)
	fmt.Printf("wrote %s\n", jsonPath)
	fmt.Printf("wrote %s\n", mdPath)

	// Exit non-zero only on a genuine parity divergence, so this can gate CI on a
	// data-bearing machine. A run with no data (all SKIP) exits 0.
	if rep.Diverge > 0 {
		return 1
	}
	return 0
}

// abs resolves p to an absolute path, leaving "" untouched and falling back to
// the original on error.
func abs(p string) string {
	if p == "" {
		return ""
	}
	if a, err := filepath.Abs(p); err == nil {
		return a
	}
	return p
}
