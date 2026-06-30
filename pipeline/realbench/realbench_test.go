package realbench

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestSubstituteArgs checks that every placeholder is replaced with the resolved
// input path, that {fai} derives from the reference, and that {out}/{outdir} are
// substituted from the per-side temp targets.
func TestSubstituteArgs(t *testing.T) {
	in := Inputs{
		Ref:    "/data/ref.fa",
		BAM:    "/data/x.bam",
		CRAM:   "/data/x.cram",
		VCF:    "/data/x.vcf.gz",
		Fastq1: "/data/R1.fq.gz",
		Fastq2: "/data/R2.fq.gz",
		BED:    "/data/i.bed",
		GFF:    "/data/g.gff.gz",
	}
	args := []string{
		"view", "-C", "-T", phRef, "-o", phOut, phBAM, phCRAM, phVCF,
		phFastq1, phFastq2, phBED, phGFF, phFai, phOutdir,
	}
	got := substituteArgs(args, in, "/tmp/out.cram", "/tmp/wd")
	want := []string{
		"view", "-C", "-T", "/data/ref.fa", "-o", "/tmp/out.cram",
		"/data/x.bam", "/data/x.cram", "/data/x.vcf.gz",
		"/data/R1.fq.gz", "/data/R2.fq.gz", "/data/i.bed", "/data/g.gff.gz",
		"/data/ref.fa.fai", "/tmp/wd",
	}
	if len(got) != len(want) {
		t.Fatalf("arg count: got %d want %d (%v)", len(got), len(want), got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("arg[%d]: got %q want %q", i, got[i], want[i])
		}
	}
}

// TestSubstituteArgs_EmptyFai checks that {fai} stays empty when no reference is
// provided (so a cell needing it would SKIP rather than substitute a bare
// ".fai").
func TestSubstituteArgs_EmptyFai(t *testing.T) {
	got := substituteArgs([]string{"-g", phFai}, Inputs{}, "", "")
	if got[1] != "" {
		t.Errorf("{fai} with no ref: got %q want empty", got[1])
	}
}

// TestMissingInput verifies the SKIP-on-missing-input decision: a cell whose
// required inputs are all present has no missing input; absent ones are named.
func TestMissingInput(t *testing.T) {
	full := Inputs{Ref: "r", BAM: "b", CRAM: "c", VCF: "v", Fastq1: "1", Fastq2: "2", BED: "e", GFF: "g"}
	if m := full.missing(NeedBAM | NeedRef | NeedVCF); m != "" {
		t.Errorf("all present: got missing=%q want empty", m)
	}
	cases := []struct {
		name string
		in   Inputs
		need InputKind
		want string
	}{
		{"no bam", Inputs{Ref: "r"}, NeedBAM, "-bam"},
		{"no ref", Inputs{BAM: "b"}, NeedBAM | NeedRef, "-ref"},
		{"no vcf", Inputs{}, NeedVCF, "-vcf"},
		{"no bed", Inputs{}, NeedBED, "-bed"},
		{"no gff", Inputs{}, NeedGFF, "-gff"},
		{"no fastq2", Inputs{Fastq1: "1"}, NeedFastq1 | NeedFastq2, "-fastq2"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if m := c.in.missing(c.need); m != c.want {
				t.Errorf("missing(%s): got %q want %q", c.name, m, c.want)
			}
		})
	}
}

// TestRunCell_SkipMissingInput drives a real cell through runCell with NO inputs
// and a resolver that needs no binaries: the cell must SKIP (never crash) and
// record the missing-input reason, with no measurements.
func TestRunCell_SkipMissingInput(t *testing.T) {
	var notes []string
	r := NewBinResolver("/nonexistent/our", "/nonexistent/up", t.TempDir(), &notes)
	cfg := WithResolver(Config{Tier: "chr20", Reps: 1, TmpDir: t.TempDir()}, r)

	spec := stdout("samtools", "samtools_view_sam", "view", NeedBAM, "view", phBAM)
	res := runCell(cfg, spec)

	if res.Parity != "SKIP" {
		t.Fatalf("parity: got %q want SKIP", res.Parity)
	}
	if !strings.Contains(res.Note, "-bam") {
		t.Errorf("note: got %q, expected it to mention -bam", res.Note)
	}
	if res.Ours != nil || res.Up != nil {
		t.Errorf("a skipped cell should carry no measurements: ours=%v up=%v", res.Ours, res.Up)
	}
}

// TestRunCell_SkipMissingBinary checks that, with the required input present but
// our binary unresolvable (a dir with no binaries and no build), the cell SKIPs
// with a clear reason rather than crashing.
func TestRunCell_SkipMissingBinary(t *testing.T) {
	bam := filepath.Join(t.TempDir(), "x.bam")
	if err := os.WriteFile(bam, []byte("not a real bam"), 0o644); err != nil {
		t.Fatal(err)
	}
	var notes []string
	// ourDir is a real (empty) dir: ourBinary looks there and finds nothing, so
	// it never tries to build — the cell SKIPs.
	r := NewBinResolver(t.TempDir(), t.TempDir(), t.TempDir(), &notes)
	cfg := WithResolver(Config{Tier: "chr20", Reps: 1, TmpDir: t.TempDir(), Inputs: Inputs{BAM: bam}}, r)

	res := runCell(cfg, stdout("samtools", "samtools_view_sam", "view", NeedBAM, "view", phBAM))
	if res.Parity != "SKIP" {
		t.Fatalf("parity: got %q want SKIP", res.Parity)
	}
	if !strings.Contains(res.Note, "samtools") {
		t.Errorf("note: got %q, expected it to mention the samtools binary", res.Note)
	}
}

// TestReportJSONShape builds a small report by hand, runs finalize, writes it,
// and verifies the JSON shape: a single object {tier, machine, cells:[...]} with
// the per-cell fields and the derived ours/upstream ratios.
func TestReportJSONShape(t *testing.T) {
	rep := &Report{
		Tier: "exome",
		Reps: 3,
		Cells: []CellRecord{
			{
				Tool: "samtools", Name: "samtools_view_sam", Subcommand: "view", Tier: "exome",
				Parity: "PASS",
				Ours:   &SideMeasure{WallS: 2.0, CPUS: 1.5, RSSKB: 200},
				Up:     &SideMeasure{WallS: 1.0, CPUS: 1.0, RSSKB: 100},
			},
			{Tool: "bedrandom", Name: "bedrandom", Subcommand: "bedtools random", Tier: "exome", Parity: "SKIP", Note: "ours-only perf cell"},
		},
	}
	rep.finalize()

	if rep.Pass != 1 || rep.Skip != 1 {
		t.Errorf("counts: pass=%d skip=%d (want 1/1)", rep.Pass, rep.Skip)
	}
	c := rep.Cells[0]
	if c.WallX == nil || *c.WallX != 2.0 {
		t.Errorf("wall_x: got %v want 2.0", c.WallX)
	}
	if c.RSSX == nil || *c.RSSX != 2.0 {
		t.Errorf("rss_x: got %v want 2.0", c.RSSX)
	}
	if c.CPUX == nil || *c.CPUX != 1.5 {
		t.Errorf("cpu_x: got %v want 1.5", c.CPUX)
	}

	dir := t.TempDir()
	jsonPath, mdPath, err := WriteReports(rep, dir)
	if err != nil {
		t.Fatal(err)
	}
	if filepath.Base(jsonPath) != "realbench.exome.json" {
		t.Errorf("json path: got %q", jsonPath)
	}
	if filepath.Base(mdPath) != "realbench.exome.md" {
		t.Errorf("md path: got %q", mdPath)
	}

	raw, err := os.ReadFile(jsonPath)
	if err != nil {
		t.Fatal(err)
	}
	// Decode into a generic map to assert the aggregate object shape.
	var obj map[string]json.RawMessage
	if err := json.Unmarshal(raw, &obj); err != nil {
		t.Fatalf("report JSON does not parse: %v", err)
	}
	for _, key := range []string{"tier", "machine", "cells"} {
		if _, ok := obj[key]; !ok {
			t.Errorf("report JSON missing top-level key %q", key)
		}
	}
	// Round-trip into the typed Report and re-check a per-cell field shape.
	var back Report
	if err := json.Unmarshal(raw, &back); err != nil {
		t.Fatalf("typed re-parse: %v", err)
	}
	if back.Tier != "exome" || len(back.Cells) != 2 {
		t.Fatalf("round-trip: tier=%q cells=%d", back.Tier, len(back.Cells))
	}
	if back.Cells[0].WallX == nil {
		t.Errorf("round-trip lost wall_x on cell 0")
	}

	// The MD table must be sorted by rss_x desc: the PASS cell (rss_x=2.0) comes
	// before the SKIP cell (no ratio).
	mdRaw, err := os.ReadFile(mdPath)
	if err != nil {
		t.Fatal(err)
	}
	md := string(mdRaw)
	iPass := strings.Index(md, "samtools_view_sam")
	iSkip := strings.Index(md, "bedrandom")
	if iPass < 0 || iSkip < 0 || iPass > iSkip {
		t.Errorf("MD not sorted by rss_x desc: pass@%d skip@%d", iPass, iSkip)
	}
}

// TestMatrixCoversToolSurface asserts the matrix actually exercises every ported
// tool family and a healthy number of cells, so a future refactor that silently
// drops a tool is caught.
func TestMatrixCoversToolSurface(t *testing.T) {
	cells := Matrix("chr20")
	if len(cells) < 100 {
		t.Errorf("matrix has only %d cells; expected the full surface (>=100)", len(cells))
	}
	tools := map[string]bool{}
	names := map[string]bool{}
	for _, c := range cells {
		tools[c.Tool] = true
		if names[c.Name] {
			t.Errorf("duplicate cell name %q", c.Name)
		}
		names[c.Name] = true
	}
	for _, want := range []string{
		"samtools", "bcftools", "seqtk", "fastp", "sickle", "skewer",
		"prinseq", "vcftools", "mosdepth", "bgzip", "tabix", "htsfile",
		"bedintersect", "bedmerge", "bedgenomecov", "bedmakewindows",
	} {
		if !tools[want] {
			t.Errorf("matrix does not cover tool %q", want)
		}
	}
	// Every bed* tool with an upstream pair should appear.
	for tool := range bedSubcommand {
		if !tools[tool] {
			t.Errorf("matrix missing bed tool %q", tool)
		}
	}
}

// TestBedUpstreamMapping checks that each bed cell resolves to the right upstream
// `bedtools <sub>` stub via the resolver (using a fake upstream dir so no real
// bedtools is needed).
func TestBedUpstreamMapping(t *testing.T) {
	upDir := t.TempDir()
	// Provide a fake bedtools so locateUpstream finds it.
	fakeBedtools := filepath.Join(upDir, "bedtools")
	if err := os.WriteFile(fakeBedtools, []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatal(err)
	}
	var notes []string
	r := NewBinResolver(t.TempDir(), upDir, t.TempDir(), &notes)

	rb := r.upstreamBinary("bedintersect")
	if rb.Path != fakeBedtools {
		t.Errorf("bedintersect upstream path: got %q want %q", rb.Path, fakeBedtools)
	}
	if len(rb.UpStub) != 1 || rb.UpStub[0] != "intersect" {
		t.Errorf("bedintersect upstream stub: got %v want [intersect]", rb.UpStub)
	}

	rb2 := r.upstreamBinary("bedmakewindows")
	if len(rb2.UpStub) != 1 || rb2.UpStub[0] != "makewindows" {
		t.Errorf("bedmakewindows stub: got %v want [makewindows]", rb2.UpStub)
	}
}
