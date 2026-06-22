package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

// binset holds the resolved binary paths for one tool on each side.
type binset struct {
	oursSamtools, upSamtools string
	oursBcftools, upBcftools string
}

func (b binset) ours(tool string) string {
	if tool == "bcftools" {
		return b.oursBcftools
	}
	return b.oursSamtools
}

func (b binset) upstream(tool string) string {
	if tool == "bcftools" {
		return b.upBcftools
	}
	return b.upSamtools
}

// inputs holds the resolved input file paths (empty when not provided).
type inputs struct {
	ref, bam, vcf string
}

// config is the fully-resolved run configuration.
type config struct {
	bins    binset
	in      inputs
	region  string
	reps    int
	outDir  string
	tmpDir  string
	verbose bool
}

// runBattery executes every cell in the battery, returning the assembled
// report. It never returns an error for an individual cell failure (those
// become ERROR/DIVERGE/SKIP rows); it only returns an error for a setup
// problem that prevents producing a report at all.
func runBattery(cfg config) (*report, error) {
	rep := &report{
		Generated: time.Now().UTC(),
		Region:    cfg.region,
		Reps:      cfg.reps,
		OurBin:    binDir(cfg.bins.oursSamtools, cfg.bins.oursBcftools),
		UpBin:     binDir(cfg.bins.upSamtools, cfg.bins.upBcftools),
	}
	rep.Inputs = describeInputs(cfg)

	for _, spec := range battery() {
		rep.Cells = append(rep.Cells, runCell(cfg, spec))
	}
	rep.finalize()
	return rep, nil
}

// runCell runs one cell: SKIP if its required input/binaries are absent,
// otherwise time and compare both sides. The returned cellResult is fully
// populated (status, measurements, ratios via finalize, diff snippet).
func runCell(cfg config, spec cellSpec) cellResult {
	res := cellResult{Name: spec.Name, Tool: spec.Tool, Multi: spec.Multi}

	ourBin := cfg.bins.ours(spec.Tool)
	upBin := cfg.bins.upstream(spec.Tool)
	if ourBin == "" || upBin == "" {
		res.Status = "SKIP"
		res.Detail = fmt.Sprintf("missing %s binary (ours=%q upstream=%q)", spec.Tool, ourBin, upBin)
		return res
	}
	if reason, ok := inputAvailable(cfg.in, spec.Need); !ok {
		res.Status = "SKIP"
		res.Detail = reason
		return res
	}

	ourOut, ourMeas, ourErr := runSide(cfg, spec, ourBin)
	upOut, upMeas, upErr := runSide(cfg, spec, upBin)

	if ourErr != nil || upErr != nil {
		// A quickcheck-style postNone cell encodes pass/fail as exit status, so a
		// non-nil error there is meaningful data, compared below. Any other
		// command failing is a hard cell error.
		if spec.Post != postNone {
			res.Status = "ERROR"
			res.Detail = fmt.Sprintf("execution failed: ours_err=%v upstream_err=%v", ourErr, upErr)
			if ourMeas != nil {
				res.Ours = ptrSide(measToSide(*ourMeas))
			}
			if upMeas != nil {
				res.Upstream = ptrSide(measToSide(*upMeas))
			}
			return res
		}
	}

	if ourMeas != nil {
		res.Ours = ptrSide(measToSide(*ourMeas))
	}
	if upMeas != nil {
		res.Upstream = ptrSide(measToSide(*upMeas))
	}

	// For postNone (quickcheck) parity is "same exit verdict".
	if spec.Post == postNone {
		ourOK := ourErr == nil
		upOK := upErr == nil
		if ourOK == upOK {
			res.Status = "PASS"
			res.Detail = fmt.Sprintf("both %s", verdictWord(ourOK))
		} else {
			res.Status = "DIVERGE"
			res.Detail = fmt.Sprintf("exit verdict differs: ours_ok=%v upstream_ok=%v", ourOK, upOK)
		}
		return res
	}

	cmp := runner.CompareByteExact(ourOut, upOut)
	if cmp.Equal {
		res.Status = "PASS"
		return res
	}
	res.Status = "DIVERGE"
	res.Detail = cmp.Detail
	res.DiffSnippet = diffSnippet(ourOut, upOut)
	return res
}

// runSide builds the argv for one binary, runs it reps times, and returns the
// comparable text stream (stdout, or the re-decoded SAM for file-producing
// cells), the reduced measurement, and any execution error.
func runSide(cfg config, spec cellSpec, bin string) ([]byte, *Measurement, error) {
	args, outPath, cleanup, err := buildArgs(cfg, spec)
	if cleanup != nil {
		defer cleanup()
	}
	if err != nil {
		return nil, nil, err
	}

	var env []string
	if spec.Need == needBAMRef && spec.WriteOut == ".cram" {
		env = cramEnv(os.Environ())
	}

	r, err := repeatRun(cfg.reps, bin, args, "", env)
	if err != nil {
		return r.Stdout, &r.Meas, err
	}

	switch spec.Post {
	case postStdout, postNone:
		return r.Stdout, &r.Meas, nil
	case postViewSAM:
		// Re-decode the written BAM/CRAM through the SAME binary's `view -h` so
		// the two outputs are compared by decoded records, not BGZF/CRAM framing
		// (the repo-documented caveat). The decode itself is not timed.
		decoded, derr := decodeAlignment(cfg, bin, outPath, spec)
		if derr != nil {
			return nil, &r.Meas, fmt.Errorf("re-decoding %s output: %w", spec.Name, derr)
		}
		return decoded, &r.Meas, nil
	default:
		return r.Stdout, &r.Meas, nil
	}
}

// decodeAlignment runs `bin view -h <file>` (with -T ref for CRAM) and returns
// stdout, so two alignment files are compared by their decoded SAM records.
func decodeAlignment(cfg config, bin, file string, spec cellSpec) ([]byte, error) {
	args := []string{"view", "-h"}
	if spec.WriteOut == ".cram" {
		args = append(args, "-T", cfg.in.ref)
	}
	args = append(args, file)
	var env []string
	if spec.WriteOut == ".cram" {
		env = cramEnv(os.Environ())
	}
	r, err := runOnce(bin, args, "", env)
	return r.Stdout, err
}

// buildArgs substitutes placeholders into the cell's argv, allocating a temp
// output path for file-producing cells and appending -region where the command
// accepts it. It returns the argv, the {out} path (empty if none), and a
// cleanup func for any temp file.
func buildArgs(cfg config, spec cellSpec) (args []string, outPath string, cleanup func(), err error) {
	repl := strings.NewReplacer(
		"{bam}", cfg.in.bam,
		"{vcf}", cfg.in.vcf,
		"{ref}", cfg.in.ref,
	)
	if spec.WriteOut != "" {
		f, ferr := os.CreateTemp(cfg.tmpDir, "realparity-*"+spec.WriteOut)
		if ferr != nil {
			return nil, "", nil, ferr
		}
		outPath = f.Name()
		f.Close()
		cleanup = func() { os.Remove(outPath) }
	}
	for _, a := range spec.Args {
		a = repl.Replace(a)
		a = strings.ReplaceAll(a, "{out}", outPath)
		args = append(args, a)
	}
	// Append -region where the command accepts a trailing region and one was
	// requested: samtools view / depth and bcftools view / query / stats take a
	// trailing region on an indexed input. norm/sort/flagstat/idxstats/stats(sam)
	// do not, so only the explicitly region-friendly cells get it.
	if cfg.region != "" && acceptsRegion(spec) {
		args = append(args, cfg.region)
	}
	return args, outPath, cleanup, nil
}

// acceptsRegion reports whether appending a trailing region to this cell's argv
// is accepted by both binaries.
func acceptsRegion(spec cellSpec) bool {
	switch spec.Name {
	case "samtools_view_sam", "samtools_view_sam_header", "samtools_depth_a",
		"bcftools_view", "bcftools_view_body", "bcftools_stats", "bcftools_query":
		return true
	}
	return false
}

// inputAvailable reports whether the required inputs for a cell are present.
func inputAvailable(in inputs, need inputKind) (string, bool) {
	switch need {
	case needBAM:
		if in.bam == "" {
			return "no -bam provided", false
		}
	case needBAMRef:
		if in.bam == "" {
			return "no -bam provided", false
		}
		if in.ref == "" {
			return "no -ref provided (CRAM path)", false
		}
	case needVCF:
		if in.vcf == "" {
			return "no -vcf provided", false
		}
	case needVCFRef:
		if in.vcf == "" {
			return "no -vcf provided", false
		}
		if in.ref == "" {
			return "no -ref provided (norm -f)", false
		}
	}
	return "", true
}

func ptrSide(s sideMeas) *sideMeas { return &s }

func verdictWord(ok bool) string {
	if ok {
		return "OK"
	}
	return "FAIL"
}

// diffSnippet returns a short provenance-stripped first-diff excerpt for the
// report. It reuses runner.StripProvenance so the embedded snippet shows the
// SAME normalized bytes the parity comparison saw.
func diffSnippet(ours, upstream []byte) string {
	a := strings.Split(string(runner.StripProvenance(ours)), "\n")
	b := strings.Split(string(runner.StripProvenance(upstream)), "\n")
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			lo := i - 1
			if lo < 0 {
				lo = 0
			}
			var sb strings.Builder
			for j := lo; j <= i && j < n; j++ {
				fmt.Fprintf(&sb, "  ours[%d]: %s\n", j+1, trunc(a[j]))
				fmt.Fprintf(&sb, "  upst[%d]: %s\n", j+1, trunc(b[j]))
			}
			return strings.TrimRight(sb.String(), "\n")
		}
	}
	if len(a) != len(b) {
		return fmt.Sprintf("line count differs: ours=%d upstream=%d", len(a), len(b))
	}
	return "streams differ (no line-level diff located)"
}

func trunc(s string) string {
	if len(s) > 160 {
		return s[:160] + "…"
	}
	return s
}

// binDir returns the common parent directory of the two resolved binaries, for
// the report header. Falls back to the samtools binary's dir.
func binDir(sam, bcf string) string {
	if sam != "" {
		return filepath.Dir(sam)
	}
	if bcf != "" {
		return filepath.Dir(bcf)
	}
	return ""
}
