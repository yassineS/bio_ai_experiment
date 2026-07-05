package realbench

import (
	"compress/gzip"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
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

// TestWriteReportsCurrentDir guards the "-out ." case the Nextflow harness uses:
// WriteReports must succeed when the output directory is "." or "". Previously it
// called MkdirAll("."), which on some FUSE-backed filesystems (Fusion on AWS
// Batch) spuriously returns "file exists" and threw away a completed run.
func TestWriteReportsCurrentDir(t *testing.T) {
	rep := &Report{Tier: "chr20", Reps: 1}
	rep.finalize()

	// Run inside a temp dir so the "." outputs land somewhere disposable.
	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Chdir(t.TempDir()); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { _ = os.Chdir(wd) })

	for _, dir := range []string{".", ""} {
		jsonPath, mdPath, err := WriteReports(rep, dir)
		if err != nil {
			t.Fatalf("WriteReports(%q) errored: %v", dir, err)
		}
		if _, err := os.Stat(jsonPath); err != nil {
			t.Errorf("WriteReports(%q): json missing: %v", dir, err)
		}
		if _, err := os.Stat(mdPath); err != nil {
			t.Errorf("WriteReports(%q): md missing: %v", dir, err)
		}
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

// TestBedCellArgWiring guards the seven previously-broken bed* cells: each must
// now have valid arg wiring (no BED3-invalid column reference, a real
// paired/BEDPE/BED4 input where the subcommand demands it, and a work dir for
// the -p prefix cell). These were HARNESS bugs (BED3 placeholder used where the
// subcommand needs more than BED3), not tool bugs; this test locks in the fix.
func TestBedCellArgWiring(t *testing.T) {
	byName := map[string]CellSpec{}
	for _, c := range Matrix("chr20") {
		byName[c.Name] = c
	}
	contains := func(args []string, want string) bool {
		for _, a := range args {
			if a == want {
				return true
			}
		}
		return false
	}
	joined := func(args []string) string { return strings.Join(args, " ") }

	// 1 & 2: bedmap / bedexpand must NOT reference col 4 against a BED3.
	for _, name := range []string{"bedmap", "bedexpand"} {
		c, ok := byName[name]
		if !ok {
			t.Fatalf("%s cell missing", name)
		}
		for i, a := range c.OurArgs {
			if a == "-c" && i+1 < len(c.OurArgs) && c.OurArgs[i+1] == "4" {
				t.Errorf("%s still uses -c 4 against BED3: %v", name, c.OurArgs)
			}
		}
		if !contains(c.OurArgs, "-c") || !contains(c.OurArgs, "3") {
			t.Errorf("%s expected -c 3 wiring, got %v", name, c.OurArgs)
		}
	}

	// 3: bedoverlap must run on the windowed input with four columns.
	if c := byName["bedoverlap"]; true {
		if !contains(c.OurArgs, phWindow) {
			t.Errorf("bedoverlap must use the derived window input, got %v", c.OurArgs)
		}
		if c.Need&NeedWindow == 0 {
			t.Errorf("bedoverlap must require NeedWindow")
		}
		if !contains(c.OurArgs, "2,3,6,7") {
			t.Errorf("bedoverlap must select four position columns, got %v", c.OurArgs)
		}
	}

	// 4 & 5: pairtopair / pairtobed must reference the BEDPE input.
	for _, name := range []string{"bedpairtopair", "bedpairtobed"} {
		c := byName[name]
		if !contains(c.OurArgs, phBEDPE) {
			t.Errorf("%s must reference the BEDPE input, got %v", name, c.OurArgs)
		}
		if c.Need&NeedBEDPE == 0 {
			t.Errorf("%s must require NeedBEDPE", name)
		}
		if c.Post != PostOursOnly {
			t.Errorf("%s should be ours-only, got post=%v", name, c.Post)
		}
	}
	// pairtobed's -a is the BEDPE and -b is the BED3.
	if c := byName["bedpairtobed"]; !contains(c.OurArgs, phBED) {
		t.Errorf("bedpairtobed -b must be the BED, got %v", c.OurArgs)
	}

	// 6: bedsplit must run in a work dir with a prefix inside it.
	if c := byName["bedsplit"]; true {
		if !c.WorkDirOut {
			t.Errorf("bedsplit must set WorkDirOut so the -p prefix dir exists")
		}
		if !strings.Contains(joined(c.OurArgs), phOutdir) {
			t.Errorf("bedsplit -p must live under {outdir}, got %v", c.OurArgs)
		}
	}

	// 7: bedtobam must feed a named (BED4+) input, not the BED3.
	if c := byName["bedtobam"]; true {
		if !contains(c.OurArgs, phBED4) {
			t.Errorf("bedtobam must feed the BED4 (named) input, got %v", c.OurArgs)
		}
		if c.Need&NeedBED4 == 0 {
			t.Errorf("bedtobam must require NeedBED4")
		}
		if contains(c.OurArgs, phBED) {
			t.Errorf("bedtobam must not feed the bare BED3, got %v", c.OurArgs)
		}
	}

	// 8: bedunionbedg must feed the derived 4-col BedGraph, not BED3 (upstream
	// SIGABRTs on BED3), and require NeedBedGraph.
	if c := byName["bedunionbedg"]; true {
		if !contains(c.OurArgs, phBedGraph) {
			t.Errorf("bedunionbedg must feed the derived bedgraph, got %v", c.OurArgs)
		}
		if c.Need&NeedBedGraph == 0 {
			t.Errorf("bedunionbedg must require NeedBedGraph")
		}
		if contains(c.OurArgs, phBED) {
			t.Errorf("bedunionbedg must not feed the bare BED3, got %v", c.OurArgs)
		}
	}

	// 9: bedtag must NOT pass both -labels and -names (upstream rejects them as
	// mutually exclusive); -names must be gone.
	if c := byName["bedtag"]; true {
		if contains(c.OurArgs, "-names") {
			t.Errorf("bedtag must not pass -names alongside -labels, got %v", c.OurArgs)
		}
		if !contains(c.OurArgs, "-labels") {
			t.Errorf("bedtag must still pass -labels, got %v", c.OurArgs)
		}
	}
}

// TestSamtoolsBcftoolsCellArgWiring locks in the harness fixes for the samtools
// fixmate/markdup and bcftools reheader cells: each must consume the derived
// prerequisite input (name-collated BAM, fixmate'd BAM, sample-rename map) and
// require the matching Need bit, so the stricter upstream oracle gets an input
// it accepts.
func TestSamtoolsBcftoolsCellArgWiring(t *testing.T) {
	byName := map[string]CellSpec{}
	for _, c := range Matrix("chr20") {
		byName[c.Name] = c
	}
	contains := func(args []string, want string) bool {
		for _, a := range args {
			if a == want {
				return true
			}
		}
		return false
	}

	// fixmate must feed the name-collated BAM and require NeedNameSortBAM, not
	// the raw {bam}.
	if c := byName["samtools_fixmate"]; true {
		if !contains(c.OurArgs, phNameBAM) {
			t.Errorf("samtools_fixmate must feed the name-collated BAM, got %v", c.OurArgs)
		}
		if c.Need&NeedNameSortBAM == 0 {
			t.Errorf("samtools_fixmate must require NeedNameSortBAM")
		}
		if contains(c.OurArgs, phBAM) {
			t.Errorf("samtools_fixmate must not feed the raw {bam}, got %v", c.OurArgs)
		}
	}

	// markdup must feed the fixmate'd (markdup-ready) BAM and require
	// NeedFixmateBAM, not the raw {bam}.
	if c := byName["samtools_markdup"]; true {
		if !contains(c.OurArgs, phFixmateBAM) {
			t.Errorf("samtools_markdup must feed the fixmate'd BAM, got %v", c.OurArgs)
		}
		if c.Need&NeedFixmateBAM == 0 {
			t.Errorf("samtools_markdup must require NeedFixmateBAM")
		}
		if contains(c.OurArgs, phBAM) {
			t.Errorf("samtools_markdup must not feed the raw {bam}, got %v", c.OurArgs)
		}
	}

	// reheader must supply a modification directive (-s <rename file>) and
	// require NeedSampleRename.
	if c := byName["bcftools_reheader"]; true {
		if !contains(c.OurArgs, "-s") || !contains(c.OurArgs, phSampleRename) {
			t.Errorf("bcftools_reheader must pass -s <sample-rename>, got %v", c.OurArgs)
		}
		if c.Need&NeedSampleRename == 0 {
			t.Errorf("bcftools_reheader must require NeedSampleRename")
		}
	}
}

// TestDeriveInputs checks that the synthetic bed* inputs are derived
// deterministically from a BED3: a BED4 with a name column, an 8-field windowed
// file, and a 10-field BEDPE. It also verifies an empty BED derives nothing.
func TestDeriveInputs(t *testing.T) {
	dir := t.TempDir()
	bed := filepath.Join(dir, "in.bed")
	// Four BED3 records (two pairs) plus a comment and a track line to skip.
	content := "track name=x\n#comment\nchr20\t100\t200\nchr20\t300\t400\nchr20\t500\t600\nchr20\t700\t800\n"
	if err := os.WriteFile(bed, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(dir, "work")
	if err := os.MkdirAll(out, 0o755); err != nil {
		t.Fatal(err)
	}
	got, err := deriveInputs(Inputs{BED: bed}, out, "")
	if err != nil {
		t.Fatalf("deriveInputs: %v", err)
	}
	if got.BED4 == "" || got.BEDPE == "" || got.Window == "" {
		t.Fatalf("derived paths empty: %+v", got)
	}
	if got.BedGraph == "" {
		t.Fatalf("derived bedgraph path empty: %+v", got)
	}

	bg := readAll(t, got.BedGraph)
	// The 4-col bedgraph has chrom/start/end plus an integer value column.
	for _, line := range nonEmptyLines(bg) {
		f := strings.Split(line, "\t")
		if len(f) != 4 {
			t.Errorf("bedgraph line has %d fields, want 4: %q", len(f), line)
		}
		if _, err := strconv.Atoi(f[3]); err != nil {
			t.Errorf("bedgraph value column %q is not an integer: %q", f[3], line)
		}
	}
	// First record 100..200 ⇒ value = length 100.
	if !strings.Contains(bg, "chr20\t100\t200\t100\n") {
		t.Errorf("bedgraph missing expected length-valued row, got:\n%s", bg)
	}

	bed4 := readAll(t, got.BED4)
	// BED4 has one line per record, each with four tab fields and a region_ name.
	for _, line := range nonEmptyLines(bed4) {
		if n := len(strings.Split(line, "\t")); n != 4 {
			t.Errorf("BED4 line has %d fields, want 4: %q", n, line)
		}
	}
	if !strings.Contains(bed4, "region_1") {
		t.Errorf("BED4 missing synthetic name column: %q", bed4)
	}

	win := readAll(t, got.Window)
	for _, line := range nonEmptyLines(win) {
		if n := len(strings.Split(line, "\t")); n != 8 {
			t.Errorf("window line has %d fields, want 8: %q", n, line)
		}
	}

	bedpe := readAll(t, got.BEDPE)
	for _, line := range nonEmptyLines(bedpe) {
		if n := len(strings.Split(line, "\t")); n != 10 {
			t.Errorf("BEDPE line has %d fields, want 10: %q", n, line)
		}
	}

	// Determinism: re-derive into a fresh dir and compare bytes.
	out2 := filepath.Join(dir, "work2")
	if err := os.MkdirAll(out2, 0o755); err != nil {
		t.Fatal(err)
	}
	got2, err := deriveInputs(Inputs{BED: bed}, out2, "")
	if err != nil {
		t.Fatal(err)
	}
	if readAll(t, got.BEDPE) != readAll(t, got2.BEDPE) {
		t.Errorf("BEDPE synthesis is not deterministic")
	}

	// Empty BED derives nothing (dependent cells SKIP).
	empty, err := deriveInputs(Inputs{}, out, "")
	if err != nil {
		t.Fatalf("deriveInputs(empty): %v", err)
	}
	if empty.BED4 != "" || empty.BEDPE != "" || empty.Window != "" {
		t.Errorf("empty BED should derive nothing, got %+v", empty)
	}
}

// TestDeriveInputs_SampleRename verifies the deterministic one-line sample-rename
// map is written whenever a VCF is present (feeding `bcftools reheader -s`), and
// that it is absent with no VCF.
func TestDeriveInputs_SampleRename(t *testing.T) {
	dir := t.TempDir()
	// A non-existent VCF path is fine: deriveInputs only needs in.VCF to be
	// non-empty to decide the rename map is relevant (it does not read the VCF).
	got, err := deriveInputs(Inputs{VCF: filepath.Join(dir, "x.vcf.gz")}, dir, "")
	if err != nil {
		t.Fatalf("deriveInputs: %v", err)
	}
	if got.SampleRename == "" {
		t.Fatalf("sample-rename map not derived from a present VCF")
	}
	if content := readAll(t, got.SampleRename); content != "RB_SAMPLE\n" {
		t.Errorf("sample-rename content = %q, want %q", content, "RB_SAMPLE\n")
	}

	none, err := deriveInputs(Inputs{}, dir, "")
	if err != nil {
		t.Fatalf("deriveInputs(empty): %v", err)
	}
	if none.SampleRename != "" {
		t.Errorf("no VCF should derive no sample-rename map, got %q", none.SampleRename)
	}
}

// TestDeriveInputs_BAMTransforms verifies the prerequisite BAM transforms are
// produced when a real samtools is available: a name-collated BAM (fixmate
// input) and a name-sort|fixmate -m|coord-sort BAM (markdup input). Both must be
// real BAMs the same samtools can decode. It SKIPs when no samtools is on PATH,
// and asserts that an empty BAM (or empty samtools) derives neither.
func TestDeriveInputs_BAMTransforms(t *testing.T) {
	samtoolsBin := locateTestSamtools(t)
	if samtoolsBin == "" {
		t.Skip("no samtools binary found (PATH or bin/{upstream,ours}); skipping BAM-transform derivation test")
	}
	dir := t.TempDir()
	// Build a small coord-sorted BAM from SAM via the same samtools.
	sam := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n" +
		"r1\t99\tchr1\t100\t60\t5M\t=\t200\t105\tACGTA\tIIIII\n" +
		"r1\t147\tchr1\t200\t40\t5M\t=\t100\t-105\tTGCAT\tIIIII\n" +
		"r2\t99\tchr1\t300\t50\t4M\t=\t400\t104\tACGT\tIIII\n" +
		"r2\t147\tchr1\t400\t55\t4M\t=\t300\t-104\tTGCA\tIIII\n"
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatal(err)
	}
	bamPath := filepath.Join(dir, "in.bam")
	if out, err := exec.Command(samtoolsBin, "view", "-b", "-o", bamPath, samPath).CombinedOutput(); err != nil {
		t.Fatalf("building input BAM: %v\n%s", err, out)
	}

	work := filepath.Join(dir, "work")
	if err := os.MkdirAll(work, 0o755); err != nil {
		t.Fatal(err)
	}
	got, err := deriveInputs(Inputs{BAM: bamPath}, work, samtoolsBin)
	if err != nil {
		t.Fatalf("deriveInputs: %v", err)
	}
	if got.NameBAM == "" {
		t.Fatalf("name-collated BAM not derived")
	}
	if got.FixmateBAM == "" {
		t.Fatalf("markdup-ready (fixmate'd) BAM not derived")
	}
	// Both derived files must be real BAMs the same samtools can decode.
	for _, p := range []string{got.NameBAM, got.FixmateBAM} {
		if out, err := exec.Command(samtoolsBin, "quickcheck", p).CombinedOutput(); err != nil {
			t.Errorf("derived BAM %s failed quickcheck: %v\n%s", p, err, out)
		}
	}
	// The markdup-ready BAM must carry the ms tag that fixmate -m added.
	out, err := exec.Command(samtoolsBin, "view", got.FixmateBAM).Output()
	if err != nil {
		t.Fatalf("view fixmate'd BAM: %v", err)
	}
	if !strings.Contains(string(out), "ms:i:") {
		t.Errorf("markdup-ready BAM missing ms tag from `fixmate -m`:\n%s", out)
	}

	// An empty BAM (or empty samtools) derives neither transform.
	none, err := deriveInputs(Inputs{}, work, samtoolsBin)
	if err != nil {
		t.Fatalf("deriveInputs(empty): %v", err)
	}
	if none.NameBAM != "" || none.FixmateBAM != "" {
		t.Errorf("no BAM should derive no BAM transforms, got %+v", none)
	}
	noBin, err := deriveInputs(Inputs{BAM: bamPath}, work, "")
	if err != nil {
		t.Fatalf("deriveInputs(no samtools): %v", err)
	}
	if noBin.NameBAM != "" || noBin.FixmateBAM != "" {
		t.Errorf("no samtools should derive no BAM transforms, got %+v", noBin)
	}
}

// locateTestSamtools resolves a samtools binary for the BAM-transform test: it
// tries PATH first, then walks up from the test's cwd looking for
// bin/upstream/samtools then bin/ours/samtools (the repo's built binaries).
// Returns "" when none is found, so the caller can SKIP.
func locateTestSamtools(t *testing.T) string {
	t.Helper()
	if p, err := exec.LookPath("samtools"); err == nil {
		return p
	}
	dir, err := os.Getwd()
	if err != nil {
		return ""
	}
	for {
		for _, rel := range []string{"bin/upstream/samtools", "bin/ours/samtools"} {
			cand := filepath.Join(dir, rel)
			if fileExists(cand) {
				return cand
			}
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return ""
		}
		dir = parent
	}
}

func readAll(t *testing.T, path string) string {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return string(b)
}

func nonEmptyLines(s string) []string {
	var out []string
	for _, l := range strings.Split(s, "\n") {
		if l != "" {
			out = append(out, l)
		}
	}
	return out
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

// TestDeriveInputs_FastqPlain verifies that deriveInputs decompresses a gzipped
// Fastq1 into a plain FASTQ (so prinseq, which cannot read gzip, gets an input
// it can parse), that the decompressed bytes round-trip, and that an absent
// Fastq1 derives nothing (the prinseq cells then SKIP).
func TestDeriveInputs_FastqPlain(t *testing.T) {
	dir := t.TempDir()
	plainBytes := "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nHHHHHHHH\n"

	// Write a gzipped R1 (the realbench tier hands prinseq the bgzipped R1).
	gzPath := filepath.Join(dir, "R1.fq.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	gw := gzip.NewWriter(f)
	if _, err := gw.Write([]byte(plainBytes)); err != nil {
		t.Fatal(err)
	}
	if err := gw.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()

	work := filepath.Join(dir, "work")
	if err := os.MkdirAll(work, 0o755); err != nil {
		t.Fatal(err)
	}
	got, err := deriveInputs(Inputs{Fastq1: gzPath}, work, "")
	if err != nil {
		t.Fatalf("deriveInputs: %v", err)
	}
	if got.FastqPlain == "" {
		t.Fatalf("FastqPlain not derived from gzipped Fastq1")
	}
	if out := readAll(t, got.FastqPlain); out != plainBytes {
		t.Errorf("decompressed FASTQ mismatch:\n got %q\nwant %q", out, plainBytes)
	}

	// No Fastq1 -> no plain FASTQ (dependent cells SKIP).
	none, err := deriveInputs(Inputs{}, work, "")
	if err != nil {
		t.Fatalf("deriveInputs(empty): %v", err)
	}
	if none.FastqPlain != "" {
		t.Errorf("empty inputs should derive no plain FASTQ, got %q", none.FastqPlain)
	}
}

// TestPrinseqCellArgWiring asserts the prinseq cells consume the DECOMPRESSED
// plain FASTQ (not the gzipped {fastq1}) and require NeedFastqPlain, so both
// ours and upstream (neither of which reads gzip) get a comparable input.
func TestPrinseqCellArgWiring(t *testing.T) {
	contains := func(args []string, want string) bool {
		for _, a := range args {
			if a == want {
				return true
			}
		}
		return false
	}
	cells := trimmerCells()
	var stats, filter *CellSpec
	for i := range cells {
		switch cells[i].Name {
		case "prinseq_stats":
			stats = &cells[i]
		case "prinseq_filter":
			filter = &cells[i]
		}
	}
	if stats == nil || filter == nil {
		t.Fatal("prinseq_stats/prinseq_filter cells not found")
	}
	for _, c := range []*CellSpec{stats, filter} {
		if !contains(c.OurArgs, phFastqPlain) {
			t.Errorf("%s must feed the plain FASTQ placeholder %q, got %v", c.Name, phFastqPlain, c.OurArgs)
		}
		if contains(c.OurArgs, phFastq1) {
			t.Errorf("%s must NOT feed the gzipped {fastq1}, got %v", c.Name, c.OurArgs)
		}
		if c.Need&NeedFastqPlain == 0 {
			t.Errorf("%s must require NeedFastqPlain", c.Name)
		}
	}
	// The stats cell must avoid the nondeterministic stats_tag group (and the
	// -stats_all superset that pulls it in): upstream's stats_tag midseq value
	// is an unsorted Perl-hash join and would DIFF run-to-run.
	if contains(stats.OurArgs, "-stats_all") {
		t.Errorf("prinseq_stats must not use -stats_all (pulls in nondeterministic stats_tag), got %v", stats.OurArgs)
	}
	if contains(stats.OurArgs, "-stats_tag") {
		t.Errorf("prinseq_stats must not use -stats_tag (nondeterministic upstream), got %v", stats.OurArgs)
	}
	if !contains(stats.OurArgs, "-stats_info") || !contains(stats.OurArgs, "-stats_len") {
		t.Errorf("prinseq_stats must request the deterministic stats subset, got %v", stats.OurArgs)
	}
}
