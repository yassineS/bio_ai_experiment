package giab

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// Status is the verdict of a single stage of the harness.
type Status string

const (
	// StatusPass means the stage ran and met its expectation.
	StatusPass Status = "PASS"
	// StatusFail means the stage ran and a substantive discrepancy was found
	// (e.g. a genotype/PASS-FAIL flip, or a benchmarking engine error).
	StatusFail Status = "FAIL"
	// StatusSkip means a prerequisite (data file or external tool) was absent.
	StatusSkip Status = "SKIP"
)

// EngineResult holds the per-call-set biological concordance metrics for one
// engine run (hap.py or vcfeval), keyed implicitly by the call set it scored.
type EngineResult struct {
	Engine  Engine         `json:"engine"`
	Stratum string         `json:"stratum"`
	Ours    []BenchMetrics `json:"ours"`
	Up      []BenchMetrics `json:"upstream"`
	Status  Status         `json:"status"`
	Detail  string         `json:"detail,omitempty"`
}

// Report is the full output of a concordance run, serialised to
// giab_concordance.{json,md}.
type Report struct {
	Sample    string    `json:"sample"`
	Build     string    `json:"build"`
	Generated time.Time `json:"generated"`

	// CallStage records whether the two call sets were produced.
	CallStatus Status `json:"call_status"`
	CallDetail string `json:"call_detail,omitempty"`

	// Concordance is the ours-vs-upstream record comparison (nil if skipped).
	Concordance       *Concordance `json:"concordance,omitempty"`
	ConcordanceStatus Status       `json:"concordance_status"`
	ConcordanceDetail string       `json:"concordance_detail,omitempty"`

	// Biological holds the per-stratum hap.py/vcfeval results (empty if skipped).
	Biological       []EngineResult `json:"biological,omitempty"`
	BiologicalStatus Status         `json:"biological_status"`
	BiologicalDetail string         `json:"biological_detail,omitempty"`

	// OursVCF / UpVCF are the produced call-set paths (when produced).
	OursVCF string `json:"ours_vcf,omitempty"`
	UpVCF   string `json:"upstream_vcf,omitempty"`
}

// Failed reports whether the run found a substantive failure (a genotype or
// PASS/FAIL flip in the ours-vs-upstream comparison, or a benchmarking error).
// SKIPs never count as failures.
func (r *Report) Failed() bool {
	if r.ConcordanceStatus == StatusFail {
		return true
	}
	if r.BiologicalStatus == StatusFail {
		return true
	}
	return false
}

// Logf is an optional structured logger; nil is a no-op.
type Logf func(format string, args ...any)

func (l Logf) printf(format string, args ...any) {
	if l != nil {
		l(format, args...)
	}
}

// runner carries resolved binaries and a scratch directory through a run.
type runner struct {
	cfg  *Config
	ours string // our bcftools
	up   string // upstream bcftools
	tmp  string
	logf Logf
}

// Run executes the harness for cfg. It never panics on missing data: each stage
// is independently guarded and contributes a SKIP with a reason when its
// prerequisites are absent. The returned Report is always non-nil.
func Run(cfg *Config, logf Logf) (*Report, error) {
	rep := &Report{
		Sample:    cfg.Sample,
		Build:     cfg.Build,
		Generated: time.Now().UTC(),
	}

	tmp, err := os.MkdirTemp("", "giab-concordance")
	if err != nil {
		return rep, err
	}
	defer os.RemoveAll(tmp)

	r := &runner{cfg: cfg, tmp: tmp, logf: logf}

	// Resolve the two bcftools binaries.
	r.ours, r.up = r.resolveBcftools(rep)

	// Stage 1: produce both call sets (or SKIP).
	oursVCF, upVCF, ok := r.produceCallSets(rep)
	if !ok {
		// Call sets unavailable: the ours-vs-upstream stage cannot run, but the
		// biological stage might still compare pre-existing call sets if a
		// truth set + engine are present. We treat absence as a SKIP cascade.
		rep.ConcordanceStatus = StatusSkip
		if rep.ConcordanceDetail == "" {
			rep.ConcordanceDetail = "call sets not produced; " + DocPointer
		}
		rep.BiologicalStatus = StatusSkip
		rep.BiologicalDetail = "no call sets to benchmark; " + DocPointer
		return rep, nil
	}
	rep.OursVCF, rep.UpVCF = oursVCF, upVCF

	// Stage 2: ours-vs-upstream record concordance (the core claim).
	r.compareCallSets(rep, oursVCF, upVCF)

	// Stage 3: biological concordance vs GIAB truth via hap.py / vcfeval.
	r.biological(rep, oursVCF, upVCF)

	return rep, nil
}

// resolveBcftools resolves OUR and UPSTREAM bcftools, recording a call-stage
// skip reason on the report when one cannot be found. It returns ("","") if
// neither resolves; callers still attempt the half they have.
func (r *runner) resolveBcftools(rep *Report) (ours, up string) {
	ours = r.cfg.OurBcftools
	if ours == "" {
		cache := filepath.Join(r.tmp, "bin")
		if p, err := upstream.OurBinary("bcftools", cache); err == nil {
			ours = p
		} else {
			rep.CallDetail = "our bcftools unavailable: " + err.Error()
		}
	} else if missingFile(ours) {
		rep.CallDetail = "configured our_bcftools does not exist: " + ours
		ours = ""
	}

	up = r.cfg.UpstreamBcftools
	if up == "" {
		if p, err := upstream.Binary("bcftools"); err == nil {
			up = p
		} else {
			if rep.CallDetail != "" {
				rep.CallDetail += "; "
			}
			rep.CallDetail += "upstream bcftools unavailable: " + err.Error()
		}
	} else if missingFile(up) {
		if rep.CallDetail != "" {
			rep.CallDetail += "; "
		}
		rep.CallDetail += "configured upstream_bcftools does not exist: " + up
		up = ""
	}
	return ours, up
}

// produceCallSets runs `bcftools mpileup -f REF BAM | bcftools call -mv -Oz`
// with each binary, returning the two bgzipped VCF paths. It SKIPs (returns
// ok=false) when the reference, reads BAM, or either binary is missing.
func (r *runner) produceCallSets(rep *Report) (oursVCF, upVCF string, ok bool) {
	var miss []string
	if missingFile(r.cfg.Reference) {
		miss = append(miss, "reference FASTA ("+orNone(r.cfg.Reference)+")")
	}
	if missingFile(r.cfg.ReadsBAM) {
		miss = append(miss, "reads BAM ("+orNone(r.cfg.ReadsBAM)+")")
	}
	if r.ours == "" {
		miss = append(miss, "our bcftools")
	}
	if r.up == "" {
		miss = append(miss, "upstream bcftools")
	}
	if len(miss) > 0 {
		rep.CallStatus = StatusSkip
		rep.CallDetail = "cannot produce call sets, missing: " + strings.Join(miss, ", ") + "; " + DocPointer
		r.logf.printf("SKIP call stage: %s", rep.CallDetail)
		return "", "", false
	}

	oursVCF = filepath.Join(r.tmp, "ours.vcf.gz")
	upVCF = filepath.Join(r.tmp, "up.vcf.gz")
	if err := r.callOne(r.ours, oursVCF); err != nil {
		rep.CallStatus = StatusFail
		rep.CallDetail = "our call set failed: " + err.Error()
		return "", "", false
	}
	if err := r.callOne(r.up, upVCF); err != nil {
		rep.CallStatus = StatusFail
		rep.CallDetail = "upstream call set failed: " + err.Error()
		return "", "", false
	}
	rep.CallStatus = StatusPass
	rep.CallDetail = "produced ours.vcf.gz and up.vcf.gz via bcftools mpileup | call -mv -Oz"
	return oursVCF, upVCF, true
}

// callOne runs the two-stage pipeline `bcftools mpileup -f REF BAM` piped into
// `bcftools call -mv -O z -o out` with the given bcftools binary, then indexes
// the result.
func (r *runner) callOne(bcftools, outVCF string) error {
	mpileup := exec.Command(bcftools, "mpileup", "-f", r.cfg.Reference, r.cfg.ReadsBAM)
	call := exec.Command(bcftools, "call", "-mv", "-O", "z", "-o", outVCF)

	pipe, err := mpileup.StdoutPipe()
	if err != nil {
		return err
	}
	call.Stdin = pipe

	var mErr, cErr bytes.Buffer
	mpileup.Stderr = &mErr
	call.Stderr = &cErr

	if err := call.Start(); err != nil {
		return fmt.Errorf("starting call: %w", err)
	}
	if err := mpileup.Run(); err != nil {
		return fmt.Errorf("mpileup: %w: %s", err, mErr.String())
	}
	if err := call.Wait(); err != nil {
		return fmt.Errorf("call: %w: %s", err, cErr.String())
	}
	// Index (tbi); a failure here is non-fatal for the record comparison (we
	// re-read the VCF directly) but vcfeval/hap.py need it, so surface it.
	idx := exec.Command(bcftools, "index", "-t", outVCF)
	var iErr bytes.Buffer
	idx.Stderr = &iErr
	if err := idx.Run(); err != nil {
		return fmt.Errorf("index %s: %w: %s", filepath.Base(outVCF), err, iErr.String())
	}
	return nil
}

// compareCallSets decompresses both call sets, restricts to the high-conf BED
// when present, runs the comparator + ULP-flip detector, and records the
// verdict. The stage FAILs only if a genotype or PASS/FAIL flip is found.
func (r *runner) compareCallSets(rep *Report, oursVCF, upVCF string) {
	oursRecs, err := r.readVCF(oursVCF)
	if err != nil {
		rep.ConcordanceStatus = StatusFail
		rep.ConcordanceDetail = "reading ours.vcf.gz: " + err.Error()
		return
	}
	upRecs, err := r.readVCF(upVCF)
	if err != nil {
		rep.ConcordanceStatus = StatusFail
		rep.ConcordanceDetail = "reading up.vcf.gz: " + err.Error()
		return
	}

	var region *RegionSet
	if !missingFile(r.cfg.HighConfBED) {
		region, err = ParseBEDFile(r.cfg.HighConfBED)
		if err != nil {
			rep.ConcordanceStatus = StatusFail
			rep.ConcordanceDetail = "parsing high-conf BED: " + err.Error()
			return
		}
	} else {
		r.logf.printf("high-confidence BED absent (%s); comparing over all sites", orNone(r.cfg.HighConfBED))
	}

	con := CompareCallSets(oursRecs, upRecs, region, r.cfg.effectiveQualULP())
	// Cap embedded diffs.
	if cap := r.cfg.effectiveMaxDiffs(); cap > 0 && len(con.Diffs) > cap {
		con.Diffs = con.Diffs[:cap]
	}
	rep.Concordance = &con
	if con.GenotypeOrFilterFlips > 0 {
		rep.ConcordanceStatus = StatusFail
		rep.ConcordanceDetail = fmt.Sprintf("%d differing site(s) flip a genotype or PASS/FAIL verdict", con.GenotypeOrFilterFlips)
	} else {
		rep.ConcordanceStatus = StatusPass
		rep.ConcordanceDetail = con.Headline()
	}
}

// readVCF decompresses a (possibly .gz) VCF via our bcftools `view` and parses
// the records. Using bcftools view avoids an in-repo gzip dependency here and
// keeps the decode identical to what the engines see.
func (r *runner) readVCF(path string) ([]VCFRecord, error) {
	if !strings.HasSuffix(path, ".gz") {
		return ParseVCFFile(path)
	}
	bin := r.ours
	if bin == "" {
		bin = r.up
	}
	if bin == "" {
		return nil, fmt.Errorf("no bcftools available to decompress %s", filepath.Base(path))
	}
	cmd := exec.Command(bin, "view", path)
	var out, errb bytes.Buffer
	cmd.Stdout, cmd.Stderr = &out, &errb
	if err := cmd.Run(); err != nil {
		return nil, fmt.Errorf("bcftools view %s: %w: %s", filepath.Base(path), err, errb.String())
	}
	return ParseVCF(&out)
}

// orNone renders an empty path as "<none>" for log/skip messages.
func orNone(p string) string {
	if p == "" {
		return "<none>"
	}
	return p
}
