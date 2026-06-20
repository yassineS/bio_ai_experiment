// Command giab-concordance is the turnkey GIAB variant-calling concordance
// harness (manuscript experiment C2/P1). Given a config describing a GIAB
// sample (reference FASTA, aligned reads BAM, GIAB truth VCF + high-confidence
// BED, optional GA4GH stratification BEDs, and the two bcftools binaries to
// drive), it:
//
//  1. Produces a call set with OUR bcftools and with the UPSTREAM bcftools
//     (bcftools mpileup -f REF BAM | bcftools call -mv -O z).
//  2. Compares the two call sets record-by-record within the high-confidence
//     BED, including a ULP-flip detector that confirms QUAL/PL last-place
//     differences do not flip a genotype or PASS/FAIL verdict.
//  3. If hap.py or vcfeval is available, scores both call sets against the GIAB
//     truth set, stratified by the provided GA4GH BEDs, and reports precision/
//     recall/F1 for SNVs and indels per stratum, ours vs upstream.
//  4. Writes giab_concordance.{json,md}.
//
// Every external prerequisite is checked; a missing one yields a clear SKIP with
// a pointer to docs/GIAB_CONCORDANCE.md and the command exits 0 ("nothing to
// do"). It therefore runs unchanged on a machine with no GIAB data — the common
// CI case. The heavy run happens on an external machine that has the GIAB data.
//
// Usage:
//
//	giab-concordance -config run.json
//	giab-concordance -config run.json -out ./reports -v
//	giab-concordance -print-config   # emit a template config to stdout
//
// See docs/GIAB_CONCORDANCE.md for data-acquisition URLs and the full recipe.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pipeline/giab"
)

func main() {
	var (
		configFlag = flag.String("config", "", "path to the JSON run config (see -print-config and docs/GIAB_CONCORDANCE.md)")
		outFlag    = flag.String("out", "", "output directory for giab_concordance.{json,md} (default: config out_dir, else current directory)")
		printCfg   = flag.Bool("print-config", false, "print an annotated template config to stdout and exit")
		verbose    = flag.Bool("v", false, "verbose logging to stderr")
		versionF   = flag.Bool("version", false, "print version and exit")
	)
	flag.BoolVar(versionF, "V", false, "print version and exit")
	flag.Usage = usage
	flag.Parse()

	if *versionF {
		fmt.Println("giab-concordance (bio_ai_experiment) v0.1.0")
		return
	}
	if *printCfg {
		fmt.Print(templateConfig)
		return
	}

	if *configFlag == "" {
		fmt.Fprintln(os.Stderr, "giab-concordance: no -config provided.")
		fmt.Fprintln(os.Stderr)
		fmt.Fprintln(os.Stderr, "This harness needs a JSON config naming the GIAB data and tool binaries.")
		fmt.Fprintln(os.Stderr, "Get a template with:")
		fmt.Fprintln(os.Stderr, "    giab-concordance -print-config > run.json")
		fmt.Fprintln(os.Stderr)
		fmt.Fprintln(os.Stderr, "Then fill in the paths (reference, reads BAM, GIAB truth VCF + high-conf BED,")
		fmt.Fprintln(os.Stderr, "stratification BEDs, our/upstream bcftools) and run:")
		fmt.Fprintln(os.Stderr, "    giab-concordance -config run.json -out ./reports")
		fmt.Fprintln(os.Stderr)
		fmt.Fprintln(os.Stderr, "Data-acquisition URLs and the full recipe are in docs/GIAB_CONCORDANCE.md.")
		// Nothing to do, but not an error condition: exit 0.
		return
	}

	cfg, err := giab.LoadConfig(*configFlag)
	if err != nil {
		fmt.Fprintln(os.Stderr, "giab-concordance:", err)
		os.Exit(2)
	}
	if *outFlag != "" {
		cfg.OutDir = *outFlag
	}

	var logf giab.Logf
	if *verbose {
		logf = func(format string, args ...any) { fmt.Fprintf(os.Stderr, format+"\n", args...) }
	}

	rep, err := giab.Run(cfg, logf)
	if err != nil {
		fmt.Fprintln(os.Stderr, "giab-concordance:", err)
		os.Exit(2)
	}

	jsonPath, mdPath, err := giab.WriteReports(rep, cfg.OutDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "giab-concordance: writing reports:", err)
		os.Exit(2)
	}

	fmt.Fprintf(os.Stderr, "stages: call=%s concordance=%s biological=%s\n",
		rep.CallStatus, rep.ConcordanceStatus, rep.BiologicalStatus)
	if rep.Concordance != nil {
		fmt.Fprintln(os.Stderr, rep.Concordance.Headline())
	}
	fmt.Fprintf(os.Stderr, "reports: %s\n         %s\n", jsonPath, mdPath)

	if rep.Failed() {
		os.Exit(1)
	}
}

func usage() {
	fmt.Fprintln(os.Stderr, "giab-concordance — GIAB variant-calling concordance harness (experiment C2/P1)")
	fmt.Fprintln(os.Stderr)
	fmt.Fprintln(os.Stderr, "Usage:")
	fmt.Fprintln(os.Stderr, "  giab-concordance -config run.json [-out DIR] [-v]")
	fmt.Fprintln(os.Stderr, "  giab-concordance -print-config > run.json")
	fmt.Fprintln(os.Stderr)
	fmt.Fprintln(os.Stderr, "It drives OUR and UPSTREAM bcftools over a GIAB sample, compares the two call")
	fmt.Fprintln(os.Stderr, "sets record-by-record (with a ULP-flip detector) within the high-confidence")
	fmt.Fprintln(os.Stderr, "BED, and — when hap.py/vcfeval is present — scores both against the GIAB truth")
	fmt.Fprintln(os.Stderr, "set per GA4GH stratification. Missing data/tools SKIP cleanly (exit 0).")
	fmt.Fprintln(os.Stderr)
	fmt.Fprintln(os.Stderr, "Flags:")
	flag.PrintDefaults()
	fmt.Fprintln(os.Stderr)
	fmt.Fprintln(os.Stderr, "Exit status: 0 = all PASS/SKIP, 1 = a genotype/PASS-FAIL flip or benchmarking")
	fmt.Fprintln(os.Stderr, "error, 2 = a usage/IO error. See docs/GIAB_CONCORDANCE.md.")
}

// templateConfig is the annotated config emitted by -print-config. JSON does not
// permit comments, so the annotations live in docs/GIAB_CONCORDANCE.md; this is
// a ready-to-edit skeleton with the canonical field set.
var templateConfig = mustTemplate()

func mustTemplate() string {
	c := giab.Config{
		Sample:           "HG002",
		Build:            "GRCh38",
		Reference:        "/data/ref/GRCh38.fa",
		ReadsBAM:         "/data/HG002/HG002.GRCh38.bam",
		TruthVCF:         "/data/giab/HG002_GRCh38_v4.2.1_benchmark.vcf.gz",
		HighConfBED:      "/data/giab/HG002_GRCh38_v4.2.1_benchmark_noinconsistent.bed",
		OurBcftools:      "",
		UpstreamBcftools: "",
		HappyBin:         "",
		VcfevalBin:       "",
		SDFTemplate:      "/data/ref/GRCh38.sdf",
		OutDir:           "./reports",
		Stratifications: []giab.Stratification{
			{Name: "CMRG", BED: "/data/strat/GRCh38_CMRG_benchmark.bed"},
			{Name: "alldifficultregions", BED: "/data/strat/GRCh38_alldifficultregions.bed.gz"},
		},
	}
	b, _ := json.MarshalIndent(c, "", "  ")
	return string(b) + "\n"
}
