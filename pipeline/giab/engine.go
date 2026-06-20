package giab

import (
	"bytes"
	"fmt"
	"os/exec"
	"path/filepath"
	"strings"
)

// resolvedEngine pairs an engine identity with its resolved binary path.
type resolvedEngine struct {
	engine Engine
	bin    string
}

// resolveEngine determines which benchmarking engine to use. Explicit config
// wins; otherwise it probes PATH for hap.py, then rtg (vcfeval). It returns
// (engine, "") with bin empty if none is available.
func (r *runner) resolveEngine() resolvedEngine {
	if r.cfg.HappyBin != "" {
		return resolvedEngine{EngineHappy, r.cfg.HappyBin}
	}
	if r.cfg.VcfevalBin != "" {
		return resolvedEngine{EngineVcfeval, r.cfg.VcfevalBin}
	}
	if p, err := exec.LookPath("hap.py"); err == nil {
		return resolvedEngine{EngineHappy, p}
	}
	if p, err := exec.LookPath("rtg"); err == nil {
		return resolvedEngine{EngineVcfeval, p}
	}
	if p, err := exec.LookPath("vcfeval"); err == nil {
		return resolvedEngine{EngineVcfeval, p}
	}
	return resolvedEngine{}
}

// biological runs the GA4GH benchmarking stage: it scores OURS and UPSTREAM
// against the GIAB truth set within the high-conf BED, once over the whole
// region (stratum "*") and once per provided stratification BED. It SKIPs when
// the truth set, the high-conf BED, or an engine is unavailable.
func (r *runner) biological(rep *Report, oursVCF, upVCF string) {
	var miss []string
	if missingFile(r.cfg.TruthVCF) {
		miss = append(miss, "GIAB truth VCF ("+orNone(r.cfg.TruthVCF)+")")
	}
	if missingFile(r.cfg.HighConfBED) {
		miss = append(miss, "high-confidence BED ("+orNone(r.cfg.HighConfBED)+")")
	}
	eng := r.resolveEngine()
	if eng.bin == "" {
		miss = append(miss, "a benchmarking engine (hap.py or vcfeval) on PATH")
	}
	if len(miss) > 0 {
		rep.BiologicalStatus = StatusSkip
		rep.BiologicalDetail = "biological concordance skipped, missing: " + strings.Join(miss, ", ") + "; " + DocPointer
		r.logf.printf("SKIP biological stage: %s", rep.BiologicalDetail)
		return
	}

	// Whole-region run (stratum "*") plus one run per stratification BED.
	strata := []Stratification{{Name: "*", BED: ""}}
	strata = append(strata, r.cfg.Stratifications...)

	var results []EngineResult
	anyFail := false
	for _, s := range strata {
		res := r.scoreStratum(eng, oursVCF, upVCF, s)
		if res.Status == StatusFail {
			anyFail = true
		}
		results = append(results, res)
	}
	rep.Biological = results
	switch {
	case anyFail:
		rep.BiologicalStatus = StatusFail
		rep.BiologicalDetail = "one or more strata failed to score"
	default:
		rep.BiologicalStatus = StatusPass
		rep.BiologicalDetail = fmt.Sprintf("scored %d stratum/strata with %s", len(results), eng.engine)
	}
}

// scoreStratum scores both call sets for one stratum and packages the metrics.
func (r *runner) scoreStratum(eng resolvedEngine, oursVCF, upVCF string, s Stratification) EngineResult {
	res := EngineResult{Engine: eng.engine, Stratum: s.Name}
	oursM, err := r.scoreOne(eng, oursVCF, s, "ours")
	if err != nil {
		res.Status = StatusFail
		res.Detail = "scoring ours: " + err.Error()
		return res
	}
	upM, err := r.scoreOne(eng, upVCF, s, "up")
	if err != nil {
		res.Status = StatusFail
		res.Detail = "scoring upstream: " + err.Error()
		return res
	}
	res.Ours = stampStratum(oursM, s.Name)
	res.Up = stampStratum(upM, s.Name)
	res.Status = StatusPass
	return res
}

// scoreOne invokes the engine on one query VCF and parses its summary.
func (r *runner) scoreOne(eng resolvedEngine, queryVCF string, s Stratification, tag string) ([]BenchMetrics, error) {
	prefix := filepath.Join(r.tmp, fmt.Sprintf("bench_%s_%s", tag, sanitize(s.Name)))
	switch eng.engine {
	case EngineHappy:
		return r.runHappy(eng.bin, queryVCF, s, prefix)
	case EngineVcfeval:
		return r.runVcfeval(eng.bin, queryVCF, s, prefix)
	default:
		return nil, fmt.Errorf("unknown engine %q", eng.engine)
	}
}

// runHappy shells out to hap.py and parses <prefix>.summary.csv.
//
//	hap.py TRUTH QUERY -r REF -f HIGHCONF [-R STRATBED] -o PREFIX
//
// Stratification by a single BED is expressed with -f (confident regions) when
// the stratum is the whole high-conf, and additionally restricted with -R/-T
// when a stratification BED is supplied. hap.py's richer --stratification TSV
// mode is documented in docs/GIAB_CONCORDANCE.md; here we drive one BED at a
// time so each stratum is an independent, parseable run.
func (r *runner) runHappy(bin, queryVCF string, s Stratification, prefix string) ([]BenchMetrics, error) {
	args := []string{
		r.cfg.TruthVCF, queryVCF,
		"-r", r.cfg.Reference,
		"-f", r.cfg.HighConfBED,
		"-o", prefix,
	}
	if s.BED != "" {
		// Restrict the comparison to the stratification region.
		args = append(args, "-T", s.BED)
	}
	if _, err := runEngineCmd(bin, args...); err != nil {
		return nil, err
	}
	return ParseHappySummaryFile(prefix + ".summary.csv")
}

// runVcfeval shells out to RTG vcfeval and parses <outdir>/summary.txt.
//
//	rtg vcfeval -b TRUTH -c QUERY -t SDF -e HIGHCONF [--bed-regions STRAT] -o OUTDIR
//
// vcfeval requires an SDF reference template (cfg.SDFTemplate); without it the
// run errors and the stage reports the actionable reason.
func (r *runner) runVcfeval(bin, queryVCF string, s Stratification, prefix string) ([]BenchMetrics, error) {
	if r.cfg.SDFTemplate == "" {
		return nil, fmt.Errorf("vcfeval needs an SDF reference template (set sdf_template); %s", DocPointer)
	}
	outDir := prefix + ".vcfeval"
	args := []string{}
	// `rtg vcfeval ...` vs a standalone `vcfeval` binary: if the resolved bin is
	// rtg, the subcommand is prepended.
	if filepath.Base(bin) == "rtg" {
		args = append(args, "vcfeval")
	}
	args = append(args,
		"-b", r.cfg.TruthVCF,
		"-c", queryVCF,
		"-t", r.cfg.SDFTemplate,
		"-e", r.cfg.HighConfBED,
		"-o", outDir,
	)
	if s.BED != "" {
		args = append(args, "--bed-regions", s.BED)
	}
	if _, err := runEngineCmd(bin, args...); err != nil {
		return nil, err
	}
	return ParseVcfevalSummaryFile(filepath.Join(outDir, "summary.txt"))
}

// runEngineCmd runs an engine command, capturing stderr for diagnostics.
func runEngineCmd(bin string, args ...string) ([]byte, error) {
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout, cmd.Stderr = &out, &errb
	if err := cmd.Run(); err != nil {
		return out.Bytes(), fmt.Errorf("%s %v: %w: %s", filepath.Base(bin), args, err, strings.TrimSpace(errb.String()))
	}
	return out.Bytes(), nil
}

// stampStratum overwrites the Stratum label on each metric (engines report "*"
// by default).
func stampStratum(ms []BenchMetrics, name string) []BenchMetrics {
	for i := range ms {
		ms[i].Stratum = name
	}
	return ms
}

// sanitize makes a stratum name safe for use in a filename.
func sanitize(s string) string {
	if s == "*" {
		return "all"
	}
	repl := func(r rune) rune {
		switch {
		case r >= 'a' && r <= 'z', r >= 'A' && r <= 'Z', r >= '0' && r <= '9', r == '-', r == '_':
			return r
		default:
			return '_'
		}
	}
	return strings.Map(repl, s)
}
