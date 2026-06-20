package giab

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestLoadConfig(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "run.json")
	body := `{
  "sample": "HG002",
  "build": "GRCh38",
  "reference": "/data/ref.fa",
  "reads_bam": "/data/reads.bam",
  "truth_vcf": "/data/truth.vcf.gz",
  "high_conf_bed": "/data/hc.bed",
  "qual_ulp": 0.25,
  "stratifications": [
    {"name": "CMRG", "bed": "/data/cmrg.bed"}
  ]
}`
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg, err := LoadConfig(path)
	if err != nil {
		t.Fatalf("LoadConfig: %v", err)
	}
	if cfg.Sample != "HG002" || cfg.Build != "GRCh38" {
		t.Fatalf("bad sample/build: %+v", cfg)
	}
	if cfg.effectiveQualULP() != 0.25 {
		t.Fatalf("qual_ulp override not applied: %v", cfg.effectiveQualULP())
	}
	if len(cfg.Stratifications) != 1 || cfg.Stratifications[0].Name != "CMRG" {
		t.Fatalf("bad stratifications: %+v", cfg.Stratifications)
	}
}

func TestLoadConfig_BadJSON(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "bad.json")
	if err := os.WriteFile(path, []byte("{not json"), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := LoadConfig(path); err == nil {
		t.Fatal("expected error on malformed JSON")
	}
}

func TestEffectiveDefaults(t *testing.T) {
	c := &Config{}
	if c.effectiveQualULP() != DefaultQualULP {
		t.Fatalf("default ULP: %v", c.effectiveQualULP())
	}
	if c.effectiveMaxDiffs() != 200 {
		t.Fatalf("default max diffs: %d", c.effectiveMaxDiffs())
	}
	c.MaxDiffs = -1
	if c.effectiveMaxDiffs() != 0 {
		t.Fatalf("unlimited sentinel: %d", c.effectiveMaxDiffs())
	}
}

// TestRun_SkipsCleanlyWithNoData is the central CI guarantee: with an empty
// config (no GIAB data, no binaries configured), Run must not error and must
// SKIP every data-dependent stage with a doc pointer.
func TestRun_SkipsCleanlyWithNoData(t *testing.T) {
	cfg := &Config{
		Sample:           "HG002",
		Reference:        "/nonexistent/ref.fa",
		ReadsBAM:         "/nonexistent/reads.bam",
		TruthVCF:         "/nonexistent/truth.vcf.gz",
		HighConfBED:      "/nonexistent/hc.bed",
		OurBcftools:      "/nonexistent/our-bcftools",
		UpstreamBcftools: "/nonexistent/up-bcftools",
	}
	rep, err := Run(cfg, nil)
	if err != nil {
		t.Fatalf("Run should not error on missing data: %v", err)
	}
	if rep.CallStatus != StatusSkip {
		t.Fatalf("call stage should SKIP, got %s (%s)", rep.CallStatus, rep.CallDetail)
	}
	if rep.ConcordanceStatus != StatusSkip {
		t.Fatalf("concordance stage should SKIP, got %s", rep.ConcordanceStatus)
	}
	if rep.BiologicalStatus != StatusSkip {
		t.Fatalf("biological stage should SKIP, got %s", rep.BiologicalStatus)
	}
	if rep.Failed() {
		t.Fatal("a skip-only run must not be a failure")
	}
	if !strings.Contains(rep.CallDetail, DocPointer) {
		t.Fatalf("skip detail should point at the doc: %q", rep.CallDetail)
	}
}

func TestWriteReports(t *testing.T) {
	dir := t.TempDir()
	con := Concordance{Common: 5, Identical: 4, Differ: 1, QualULPOnly: 1}
	rep := &Report{
		Sample:            "HG002",
		Build:             "GRCh38",
		CallStatus:        StatusPass,
		Concordance:       &con,
		ConcordanceStatus: StatusPass,
		ConcordanceDetail: con.Headline(),
		BiologicalStatus:  StatusSkip,
		BiologicalDetail:  "no engine; " + DocPointer,
		Biological: []EngineResult{{
			Engine:  EngineHappy,
			Stratum: "*",
			Status:  StatusPass,
			Ours:    []BenchMetrics{{VarType: "SNP", Stratum: "*", Recall: 0.99, Precision: 0.999, F1: 0.994}},
			Up:      []BenchMetrics{{VarType: "SNP", Stratum: "*", Recall: 0.99, Precision: 0.999, F1: 0.994}},
		}},
	}
	jsonPath, mdPath, err := WriteReports(rep, dir)
	if err != nil {
		t.Fatalf("WriteReports: %v", err)
	}
	jb, err := os.ReadFile(jsonPath)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(jb), "\"concordance\"") {
		t.Fatalf("json missing concordance: %s", jb)
	}
	mb, err := os.ReadFile(mdPath)
	if err != nil {
		t.Fatal(err)
	}
	md := string(mb)
	if !strings.Contains(md, "ULP-flip result") {
		t.Fatalf("markdown missing ULP-flip section:\n%s", md)
	}
	if !strings.Contains(md, "Biological concordance") {
		t.Fatalf("markdown missing biological section")
	}
	if !strings.Contains(md, "Recall (ours / up)") {
		t.Fatalf("markdown missing side-by-side P/R/F1 header")
	}
}

func TestReportFailed(t *testing.T) {
	r := &Report{ConcordanceStatus: StatusFail}
	if !r.Failed() {
		t.Fatal("concordance FAIL should mark the run failed")
	}
	r2 := &Report{BiologicalStatus: StatusFail}
	if !r2.Failed() {
		t.Fatal("biological FAIL should mark the run failed")
	}
	r3 := &Report{ConcordanceStatus: StatusSkip, BiologicalStatus: StatusSkip}
	if r3.Failed() {
		t.Fatal("all-skip run should not be failed")
	}
}
