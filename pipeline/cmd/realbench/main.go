// Command realbench runs OUR bioinformatics ports against the UPSTREAM binaries
// on REAL data (a chr20 / exome / wgs tier of GIAB-class BAM/CRAM/VCF + indexed
// reference, FASTQ pairs, intervals BED and a genes GFF) over a matrix covering
// the FULL ported surface — every subcommand of samtools, bcftools, the ~41
// bed* tools, seqtk, fastp, sickle, skewer, prinseq, vcftools, mosdepth, bgzip,
// tabix and htsfile, plus the principal and memory/format-relevant flags per
// subcommand.
//
// For each cell it reports BOTH:
//
//   - parity: byte-exact-after-provenance-stripping equality of the two outputs
//     (the project's exact parity definition, runner.StripProvenance —
//     re-decoding BAM/CRAM/BGZF through a view/gunzip step so framing
//     differences do not count), as PASS / DIFF / SKIP / ERROR, and
//   - performance: wall-clock, CPU (user+sys), and peak RSS for each side, plus
//     the ours/upstream ratio (wall_x / cpu_x / rss_x).
//
// It supersedes the synthetic pipeline/bench harness: the data is real and the
// surface is the whole tool set rather than a curated micro-benchmark.
//
// A cell whose required input or required binary is absent SKIPs (recording a
// reason) rather than crashing, so the command runs partially with whatever data
// and binaries are present and exits cleanly on a machine with no real data.
//
// Usage:
//
//	realbench -tier chr20 -ref GRCh38.fa -bam HG002.bam -cram HG002.cram \
//	  -vcf HG002.vcf.gz -fastq1 R1.fq.gz -fastq2 R2.fq.gz -bed intervals.bed \
//	  -gff genes.gff.gz -our-bin ./ourbin -upstream-bin ./upbin -reps 3 -out ./reports
//
// It writes <out>/realbench.<tier>.json and <out>/realbench.<tier>.md.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pipeline/realbench"
)

func main() {
	os.Exit(run(os.Args[1:]))
}

func run(argv []string) int {
	fs := flag.NewFlagSet("realbench", flag.ContinueOnError)
	var (
		tier    = fs.String("tier", "chr20", "data tier label: chr20 | exome | wgs")
		ref     = fs.String("ref", "", "indexed reference FASTA (.fa with .fai)")
		bam     = fs.String("bam", "", "indexed BAM (.bam with .bai)")
		cram    = fs.String("cram", "", "indexed CRAM (.cram with .crai)")
		vcf     = fs.String("vcf", "", "bgzipped + indexed VCF (.vcf.gz with .tbi)")
		fastq1  = fs.String("fastq1", "", "read-1 FASTQ (R1.fq.gz)")
		fastq2  = fs.String("fastq2", "", "read-2 FASTQ (R2.fq.gz)")
		bed     = fs.String("bed", "", "intervals BED")
		gff     = fs.String("gff", "", "genes GFF (.gff.gz)")
		ourBin  = fs.String("our-bin", "", "dir of our built tool binaries; if empty, build them into a temp dir from tools/*/cmd")
		upBin   = fs.String("upstream-bin", "", "dir of upstream tool binaries; default resolves the vendored reference_code/ build locations")
		reps    = fs.Int("reps", 3, "measurement repetitions (wall/CPU = min, RSS = max)")
		outDir  = fs.String("out", "", "report output directory (default: current directory)")
		tmpRoot = fs.String("tmp", "", "scratch directory for intermediate files (default: system temp). Point this at the task work dir on cloud/Fusion so scratch is fast, visible, and off the small instance root volume.")
		reportO = fs.Bool("report-only", false, "always exit 0 even when parity divergences (DIFF) are found; use for benchmark collection (the default exits 1 on DIFF to gate CI).")
		verbose = fs.Bool("v", false, "verbose progress logging to stderr")
		showVer = fs.Bool("version", false, "print version and exit")
	)
	fs.BoolVar(showVer, "V", false, "print version and exit")
	if err := fs.Parse(argv); err != nil {
		return 2
	}
	if *showVer {
		fmt.Println("realbench (bio_ai_experiment) v0.1.0")
		return 0
	}
	if *reps < 1 {
		*reps = 1
	}

	// Scratch root: default to the system temp dir, but honour -tmp so cloud runs
	// can keep intermediates in the (large, fast, S3-backed) task work dir rather
	// than the instance's small root /tmp. If -tmp names a dir that does not yet
	// exist, create it first.
	if *tmpRoot != "" {
		if err := os.MkdirAll(*tmpRoot, 0o755); err != nil && !os.IsExist(err) {
			fmt.Fprintln(os.Stderr, "realbench: creating -tmp dir:", err)
			return 1
		}
	}
	tmpDir, err := os.MkdirTemp(*tmpRoot, "realbench-")
	if err != nil {
		fmt.Fprintln(os.Stderr, "realbench: creating temp dir:", err)
		return 1
	}
	defer os.RemoveAll(tmpDir)
	ourCache := filepath.Join(tmpDir, "ourbins")

	var notes []string
	resolver := realbench.NewBinResolver(*ourBin, *upBin, ourCache, &notes)

	cfg := realbench.Config{
		Tier:    *tier,
		Inputs:  realbench.ResolveInputs(*ref, *bam, *cram, *vcf, *fastq1, *fastq2, *bed, *gff),
		Reps:    *reps,
		OutDir:  *outDir,
		TmpDir:  tmpDir,
		Verbose: *verbose,
	}
	cfg = realbench.WithResolver(cfg, resolver)

	if *verbose {
		fmt.Fprintf(os.Stderr, "realbench: tier=%s reps=%d out=%q\n", *tier, *reps, *outDir)
	}

	rep, err := realbench.Run(cfg)
	if err != nil {
		fmt.Fprintln(os.Stderr, "realbench:", err)
		return 1
	}

	if *verbose {
		for _, n := range notes {
			fmt.Fprintln(os.Stderr, "realbench: note:", n)
		}
	}

	jsonPath, mdPath, err := realbench.WriteReports(rep, *outDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "realbench: writing reports:", err)
		return 1
	}

	fmt.Printf("realbench tier=%s: PASS=%d DIFF=%d SKIP=%d ERROR=%d (of %d cells)\n",
		rep.Tier, rep.Pass, rep.Diff, rep.Skip, rep.Errored, len(rep.Cells))
	fmt.Printf("wrote %s\n", jsonPath)
	fmt.Printf("wrote %s\n", mdPath)

	// Exit non-zero on a genuine parity divergence, so this can gate CI on a
	// data-bearing machine. A run with no data (all SKIP) exits 0. With
	// -report-only the run always exits 0 (the reports are written either way),
	// so a benchmark pipeline can collect DIFF/perf results without the harness
	// tripping the orchestrator's error handling.
	if rep.Diff > 0 && !*reportO {
		return 1
	}
	return 0
}
