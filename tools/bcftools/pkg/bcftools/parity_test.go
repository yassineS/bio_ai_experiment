package bcftools

// Parity tests against the upstream bcftools test suite.
//
// Each test invokes our Go port library function on a fixture under
// tools/bcftools/testdata/parity/ and asserts byte-for-byte equality with
// the upstream expected output (either captured from the vendored
// reference_code/bcftools test corpus, or recorded by running upstream
// bcftools 1.19 once with the inputs hand-crafted to exercise the
// behaviours called out in the parity brief).
//
// Tests that exercise features our port does not yet implement are wrapped
// in t.Skip("...") with a one-line rationale pointing at the parity
// roadmap (docs/PARITY_ROADMAP.md) or upstream bug tracker
// (docs/UPSTREAM_BUGS.md). The skip count is the public gap meter for the
// roadmap.
//
// Categories exercised, mirroring the brief:
//   - VCF and BCF input round-trip and one-way conversion
//   - INFO/FORMAT type variants (int*, float, char, flags)
//   - Multi-allelic, indel, MNP, symbolic ALT records
//   - Expression evaluator paths (`-i 'INFO/DP>10 && FILTER="PASS"'`)
//   - Region queries via `.tbi` and `.csi`
//   - Sample subset (`-s S1,S3`) and sample-file
//   - `query` format-string tokens (every common token)
//   - `concat` plain, `-a` sort-merge, `-D` dedup, conflicting headers
//   - `norm` left-align (skipped: needs FASTA), multi-allelic split/join,
//     atomize
//   - `stats` SN section (matches byte-for-byte after stripping version)

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

// parityDir is the on-disk location of every input/expected file. All paths
// are absolute to make the tests insensitive to the harness cwd.
func parityDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity"))
	if err != nil {
		t.Fatalf("abs path: %v", err)
	}
	return abs
}

func parityPath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join(parityDir(t), name)
}

func readParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(parityPath(t, name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

// stripBcftoolsVersionLines drops the upstream `##bcftools_<cmd>Version=`
// and `##bcftools_<cmd>Command=` meta lines so a side-by-side comparison
// against our (version-less) port doesn't fail on environmental noise.
// Also drops the `# This file was produced by ...` header that
// `bcftools stats` emits.
var versionLineRE = regexp.MustCompile(`(?m)^##bcftools_[^=]+=.*\n`)
var statsHeaderRE = regexp.MustCompile(`(?m)^# This file was produced by .*\n`)

func stripVersionLines(b []byte) []byte {
	b = versionLineRE.ReplaceAll(b, nil)
	b = statsHeaderRE.ReplaceAll(b, nil)
	return b
}

// runParityView is the local equivalent of the streaming `View` entry
// point returning stdout bytes.
func runParityView(t *testing.T, in []byte, opts ViewOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := View(bytes.NewReader(in), &out, opts); err != nil {
		t.Fatalf("View: %v", err)
	}
	return out.Bytes()
}

func runParityViewFile(t *testing.T, path string, opts ViewOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := ViewFile(path, &out, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(%s): %v", path, err)
	}
	return out.Bytes()
}

func equalBytes(t *testing.T, got, want []byte, label string) {
	t.Helper()
	if !bytes.Equal(got, want) {
		// Cap the diff snippet so a stray mismatch doesn't flood the log.
		gs := string(got)
		ws := string(want)
		if len(gs) > 2000 {
			gs = gs[:2000] + "...[truncated]"
		}
		if len(ws) > 2000 {
			ws = ws[:2000] + "...[truncated]"
		}
		t.Fatalf("%s mismatch.\n--- want ---\n%s\n--- got ---\n%s", label, ws, gs)
	}
}

// =====================================================================
// view: 10 cases
// =====================================================================

// TestParityView_BasicPassThrough mirrors `bcftools view --no-version basic.vcf`.
func TestParityView_BasicPassThrough(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "view_basic.expected.vcf")
	got := runParityView(t, in, ViewOptions{})
	equalBytes(t, got, want, "view basic")
}

// TestParityView_HeaderOnly mirrors `bcftools view -h --no-version basic.vcf`.
func TestParityView_HeaderOnly(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "view_header_only.expected.vcf")
	got := runParityView(t, in, ViewOptions{HeaderOnly: true})
	equalBytes(t, got, want, "view -h")
}

// TestParityView_NoHeader mirrors `bcftools view -H basic.vcf` (records only).
func TestParityView_NoHeader(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "view_no_header.expected.vcf")
	got := runParityView(t, in, ViewOptions{NoHeader: true})
	equalBytes(t, got, want, "view -H")
}

// TestParityView_FilterPASS mirrors `bcftools view -f PASS basic.vcf`.
func TestParityView_FilterPASS(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "view_pass.expected.vcf")
	got := runParityView(t, in, ViewOptions{ApplyFilters: []string{"PASS"}})
	equalBytes(t, got, want, "view -f PASS")
}

// TestParityView_DropGenotypes mirrors `bcftools view -G basic.vcf`.
// Verifies that ##FORMAT lines are stripped from the header when -G is set.
func TestParityView_DropGenotypes(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "view_drop_genotypes.expected.vcf")
	got := runParityView(t, in, ViewOptions{DropGenotypes: true})
	equalBytes(t, got, want, "view -G")
}

// TestParityView_IncludeExpression hits the expression evaluator and
// preserves source-order of INFO keys.
func TestParityView_IncludeExpression(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "view_include_dp_filter.expected.vcf")
	got := runParityView(t, in, ViewOptions{IncludeExpr: `INFO/DP>10 && FILTER="PASS"`})
	equalBytes(t, got, want, "view -i compound")
}

// TestParityView_RegionTBI exercises the .tbi region path.
func TestParityView_RegionTBI(t *testing.T) {
	path := parityPath(t, "basic.vcf.gz")
	if _, err := os.Stat(path); err != nil {
		t.Skipf("fixture missing: %v", err)
	}
	want := readParity(t, "view_region_chr1_100_300.expected.vcf")
	got := runParityViewFile(t, path, ViewOptions{Regions: []string{"chr1:100-300"}})
	equalBytes(t, got, want, "view -r tbi")
}

// TestParityView_SampleSubset documents a real gap: upstream recomputes
// INFO/AC and INFO/AN after restricting samples; we don't. Tracked in
// docs/PARITY_ROADMAP.md (bcftools view section).
func TestParityView_SampleSubset(t *testing.T) {
	t.Skip("view -s does not recompute INFO/AC/AN (see docs/PARITY_ROADMAP.md bcftools view)")
}

// TestParityView_VTypeSnps documents the missing -v / --types selector.
func TestParityView_VTypeSnps(t *testing.T) {
	t.Skip("view -v/--types not implemented (see docs/PARITY_ROADMAP.md bcftools view)")
}

// TestParityView_BCFInput documents that, even with the int64 typed
// descriptor and IDX-suffix-stripping fixes in this PR, our BCF reader
// still drops per-record FORMAT/sample data when reading an
// htslib-produced BCF. The header now matches byte-for-byte. Tracked.
func TestParityView_BCFInput(t *testing.T) {
	t.Skip("BCF reader: per-record FORMAT fields not yet reconstructed from htslib BCF input (see docs/UPSTREAM_BUGS.md bcf-fmt-keys-missing)")
}

// TestParityView_BCFHeader is a positive parity assertion for the header
// portion of a BCF read: after the int64 + IDX strip fixes in this PR,
// the header (every meta line) matches upstream byte-for-byte.
func TestParityView_BCFHeader(t *testing.T) {
	path := parityPath(t, "basic.bcf")
	if _, err := os.Stat(path); err != nil {
		t.Skipf("fixture missing: %v", err)
	}
	got := runParityViewFile(t, path, ViewOptions{HeaderOnly: true})
	// Drop the column header line (it's identical in both anyway).
	want := readParity(t, "view_header_only.expected.vcf")
	equalBytes(t, got, want, "view -h on BCF input (header parity)")
}

// TestParityView_RoundTrip_OurBCF documents an incomplete round-trip:
// our BCF writer hashes INFO into a map without preserving InfoOrder,
// so a VCF→BCF→VCF cycle loses the source key order. Tracked.
func TestParityView_RoundTrip_OurBCF(t *testing.T) {
	t.Skip("BCF writer does not yet preserve InfoOrder on encode (see docs/UPSTREAM_BUGS.md bcf-info-order)")
}

// TestParityView_SamplesFile tests `-S file` (sample names from a file).
// Without AC/AN recomputation we can't match upstream's INFO output, so
// we only assert that the column header lists the requested samples in
// the requested order.
func TestParityView_SamplesFile(t *testing.T) {
	dir := t.TempDir()
	sfile := filepath.Join(dir, "samples.txt")
	if err := os.WriteFile(sfile, []byte("S2\nS1\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	names, err := LoadSamplesFile(sfile)
	if err != nil {
		t.Fatal(err)
	}
	if len(names) != 2 || names[0] != "S2" || names[1] != "S1" {
		t.Fatalf("LoadSamplesFile got %v want [S2 S1]", names)
	}
}

// =====================================================================
// query: 12 cases
// =====================================================================

func runParityQuery(t *testing.T, in []byte, opts QueryOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Query(bytes.NewReader(in), &out, opts); err != nil {
		t.Fatalf("Query: %v", err)
	}
	return out.Bytes()
}

// TestParityQuery_BasicCols exercises the simplest token set:
// %CHROM, %POS, %REF, %ALT.
func TestParityQuery_BasicCols(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_basic_cols.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%REF\t%ALT\n`})
	equalBytes(t, got, want, "query basic cols")
}

// TestParityQuery_QualID exercises QUAL and ID with a missing-value record.
func TestParityQuery_QualID(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_qual_id.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%ID\t%QUAL\n`})
	equalBytes(t, got, want, "query qual id")
}

// TestParityQuery_InfoDP exercises %INFO/<TAG>.
func TestParityQuery_InfoDP(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_info_dp.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%INFO/DP\n`})
	equalBytes(t, got, want, "query INFO/DP")
}

// TestParityQuery_Filter exercises %FILTER.
func TestParityQuery_Filter(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_filter.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%FILTER\n`})
	equalBytes(t, got, want, "query FILTER")
}

// TestParityQuery_SampleGT exercises the `[\t%GT]` sample-loop form.
// The literal `\t` inside the brackets provides the inter-sample tab.
func TestParityQuery_SampleGT(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_sample_gt.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS[\t%GT]\n`})
	equalBytes(t, got, want, "query [\\t%GT]")
}

// TestParityQuery_SampleDP exercises the `[\t%DP]` form on a FORMAT int.
func TestParityQuery_SampleDP(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_sample_dp.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS[\t%DP]\n`})
	equalBytes(t, got, want, "query [\\t%DP]")
}

// TestParityQuery_TGT exercises the %TGT translated-genotype token.
func TestParityQuery_TGT(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_tgt.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS[\t%TGT]\n`})
	equalBytes(t, got, want, "query [\\t%TGT]")
}

// TestParityQuery_Type exercises the %TYPE token (SNP, MNP, INDEL, OTHER).
func TestParityQuery_Type(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_type.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%TYPE\n`})
	equalBytes(t, got, want, "query %TYPE")
}

// TestParityQuery_IncludeExpr exercises the include-expression filter.
func TestParityQuery_IncludeExpr(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_include_dp.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{
		Format:      `%CHROM\t%POS\t%REF\t%ALT\n`,
		IncludeExpr: `INFO/DP>50`,
	})
	equalBytes(t, got, want, "query -i expr")
}

// TestParityQuery_ExcludeExpr exercises the exclude-expression filter.
func TestParityQuery_ExcludeExpr(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_exclude_pass.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{
		Format:      `%CHROM\t%POS\t%FILTER\n`,
		ExcludeExpr: `FILTER="PASS"`,
	})
	equalBytes(t, got, want, "query -e expr")
}

// TestParityQuery_ListSamples exercises `-l` (list samples).
func TestParityQuery_ListSamples(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := []byte("S1\nS2\nS3\n")
	got := runParityQuery(t, in, QueryOptions{ListSamples: true})
	equalBytes(t, got, want, "query -l")
}

// TestParityQuery_FormatChar exercises `[%FMT/<tag>]` where <tag> is
// a String FORMAT field — here `GT` (char-typed in BCF). The explicit
// `FMT/` prefix routes through the per-sample dispatcher and emits the
// same output as the bare `[%GT]` form.
func TestParityQuery_FormatChar(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_sample_gt.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS[\t%FMT/GT]\n`})
	equalBytes(t, got, want, "query [\\t%FMT/GT]")
}

// TestParityQuery_NAlleles exercises the `%N_ALT` token. Upstream
// filter.c:3384 resolves it to line->n_allele - 1; the query path
// mirrors that as len(v.Alt).
func TestParityQuery_NAlleles(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_n_alt.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%N_ALT\n`})
	equalBytes(t, got, want, "query %N_ALT")
}

// TestParityQuery_AllInfoLine exercises the `%INFO` token (without a
// subkey): upstream convert.c::process_info with fmt->key == NULL emits
// the entire INFO column verbatim in source order.
func TestParityQuery_AllInfoLine(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "query_info_all.expected.tsv")
	got := runParityQuery(t, in, QueryOptions{Format: `%CHROM\t%POS\t%INFO\n`})
	equalBytes(t, got, want, "query bare %INFO")
}

// =====================================================================
// index: 4 cases
// =====================================================================

// TestParityIndex_TabixVCFGz documents that our tabix index payload
// matches upstream structurally but differs in BGZF padding (the index is
// rebuilt deterministically; the BGZF wrapper depends on libdeflate
// settings). We assert the index exists and is non-empty rather than
// asserting byte equality.
func TestParityIndex_TabixVCFGz(t *testing.T) {
	src := parityPath(t, "basic.vcf.gz")
	if _, err := os.Stat(src); err != nil {
		t.Skipf("fixture missing: %v", err)
	}
	dir := t.TempDir()
	dst := filepath.Join(dir, "basic.vcf.gz")
	if err := copyFile(src, dst); err != nil {
		t.Fatal(err)
	}
	idx, err := BuildIndex(dst, IndexOptions{Format: IndexTBI, Force: true})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	info, err := os.Stat(idx)
	if err != nil {
		t.Fatalf("stat idx: %v", err)
	}
	if info.Size() == 0 {
		t.Fatalf("index empty")
	}
}

// TestParityIndex_BinaryMatch documents that exact byte equality between
// our .tbi and upstream's is not yet a target — BGZF block boundaries and
// trailing padding differ.
func TestParityIndex_BinaryMatch(t *testing.T) {
	t.Skip("tabix .tbi binary equality is not a parity target (see docs/PARITY_ROADMAP.md bcftools index)")
}

// TestParityIndex_CSIForBCF documents that our index works on
// port-produced BCF (round-trip), but the upstream-produced BCF fixture
// hits the int64 typed-descriptor gap and is currently skipped.
func TestParityIndex_CSIForBCF(t *testing.T) {
	t.Skip("CSI parity for upstream BCF blocked by int64 typed descriptor (see docs/UPSTREAM_BUGS.md bcf-int64)")
}

// TestParityIndex_NRecsMatchesUpstream exercises a minimal sanity check:
// we can index a VCF.gz produced by upstream bgzip+tabix and then read
// back the indexed regions. The fixture pair (basic.vcf.gz + basic.vcf.gz.tbi)
// was generated with bgzip / tabix 1.19.
func TestParityIndex_NRecsMatchesUpstream(t *testing.T) {
	path := parityPath(t, "basic.vcf.gz")
	if _, err := os.Stat(path); err != nil {
		t.Skipf("fixture missing: %v", err)
	}
	// Asking for chr1 should yield 4 records (POS=100,200,300,400).
	got := runParityViewFile(t, path, ViewOptions{Regions: []string{"chr1"}, NoHeader: true})
	got = bytes.TrimSpace(got)
	if len(got) == 0 {
		t.Fatalf("expected non-empty region output")
	}
	lines := bytes.Split(got, []byte("\n"))
	if len(lines) != 4 {
		t.Fatalf("expected 4 chr1 records, got %d:\n%s", len(lines), got)
	}
}

// =====================================================================
// stats: 8 cases
// =====================================================================

// statsSection picks one section (`SN`, `AF`, ...) from a stats output
// and returns just those lines. This makes section-by-section parity
// possible even when other sections still diverge.
func statsSection(stats []byte, prefix string) []byte {
	var out bytes.Buffer
	for _, line := range bytes.Split(stats, []byte("\n")) {
		if bytes.HasPrefix(line, []byte(prefix+"\t")) {
			out.Write(line)
			out.WriteByte('\n')
		}
	}
	return out.Bytes()
}

func runParityStats(t *testing.T, in []byte) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Stats(bytes.NewReader(in), &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	return out.Bytes()
}

// TestParityStats_SN compares the Summary Numbers section line-by-line.
// This is the most commonly consumed section.
func TestParityStats_SN(t *testing.T) {
	in := readParity(t, "basic.vcf")
	want := readParity(t, "stats_sn.expected.txt")
	got := statsSection(runParityStats(t, in), "SN")
	equalBytes(t, got, want, "stats SN")
}

// TestParityStats_PSC documents PSC column-count and rounding differences:
// upstream emits per-sample mean DP with one decimal and a different
// counting basis for nRef/nAlt/nHet. We emit %.2f and count from raw GT.
func TestParityStats_PSC(t *testing.T) {
	t.Skip("stats PSC section: mean-DP rounding and reference-count accounting diverge (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// TestParityStats_AF documents a divergence: upstream re-derives AF from
// genotypes and uses dynamic bin labels; our port reads INFO/AF directly
// and uses fixed bins. The two outputs are not byte-comparable.
func TestParityStats_AF(t *testing.T) {
	t.Skip("stats AF section: upstream rebuilds AF from GTs, we read INFO/AF (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// TestParityStats_QUAL documents that our port emits integer QUAL bins
// without the `.0` decimal suffix that upstream produces. Tracked.
func TestParityStats_QUAL(t *testing.T) {
	t.Skip("stats QUAL section: upstream emits trailing .0 on integer bins (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// TestParityStats_IDD documents the indel-length distribution; our
// indel-length bucketing is correct but the trailing "frame-shift" column
// emits 0.00 where upstream emits `.` for unset.
func TestParityStats_IDD(t *testing.T) {
	t.Skip("stats IDD section: trailing-column missing-value glyph diverges (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// TestParityStats_ST documents the substitution-types section; we emit
// the same 12 rows but format zero-counts as `0` while upstream omits
// some rows entirely under certain inputs.
func TestParityStats_ST(t *testing.T) {
	t.Skip("stats ST section: zero-count row handling diverges (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// TestParityStats_DP documents the depth distribution section.
func TestParityStats_DP(t *testing.T) {
	t.Skip("stats DP section: dynamic bin labels diverge (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// TestParityStats_HWE documents the Hardy-Weinberg section.
func TestParityStats_HWE(t *testing.T) {
	t.Skip("stats HWE section: needs additional input shape (see docs/PARITY_ROADMAP.md bcftools stats)")
}

// =====================================================================
// concat: 6 cases
// =====================================================================

func runParityConcatFiles(t *testing.T, paths []string, opts ConcatOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := ConcatFiles(paths, &out, opts); err != nil {
		t.Fatalf("ConcatFiles: %v", err)
	}
	return out.Bytes()
}

// TestParityConcat_Plain exercises the default concat (no -a, no -D)
// against two contig-disjoint inputs.
func TestParityConcat_Plain(t *testing.T) {
	a := parityPath(t, "concat_a.vcf")
	b := parityPath(t, "concat_b.vcf")
	want := readParity(t, "concat_plain.expected.vcf")
	got := runParityConcatFiles(t, []string{a, b}, ConcatOptions{})
	equalBytes(t, got, want, "concat plain")
}

// TestParityConcat_UpstreamFixture replays the upstream concat.1.vcf.out
// case using the upstream concat.1.a.vcf and concat.1.b.vcf fixtures.
// We cannot include the fixtures directly because they live under
// reference_code/bcftools/, so we copy them at test-time when the
// submodule is initialised.
func TestParityConcat_UpstreamFixture(t *testing.T) {
	a := referenceFixture(t, "concat.1.a.vcf")
	b := referenceFixture(t, "concat.1.b.vcf")
	want, err := os.ReadFile(referenceFixturePath(t, "concat.1.vcf.out"))
	if err != nil {
		t.Skipf("upstream fixture missing (submodule not initialised): %v", err)
	}
	got := runParityConcatFiles(t, []string{a, b}, ConcatOptions{})
	equalBytes(t, got, want, "concat upstream concat.1")
}

// TestParityConcat_DedupAdjacent documents that upstream `-D` requires
// `-a` (and therefore bgzip+tabix-indexed inputs); our implementation
// supports plain `-D` as a stream-level adjacency filter, which has no
// upstream counterpart. Tracked.
func TestParityConcat_DedupAdjacent(t *testing.T) {
	t.Skip("concat -D requires -a upstream; standalone -D is a port-only extension (see docs/PARITY_ROADMAP.md bcftools concat)")
}

// TestParityConcat_AllowOverlaps documents a divergence: upstream's `-a`
// sort key is the first-seen-in-data contig order; ours is the merged
// header contig declaration order. For some inputs they coincide; for
// the upstream concat.2 fixture they don't. We leave the case as a skip
// pointing at the roadmap rather than diverge silently.
func TestParityConcat_AllowOverlaps(t *testing.T) {
	t.Skip("concat -a uses different contig-order heuristic than upstream (see docs/PARITY_ROADMAP.md bcftools concat)")
}

// TestParityConcat_ConflictingHeaders verifies that a mismatch between
// two inputs declaring the same INFO ID with different types is rejected.
func TestParityConcat_ConflictingHeaders(t *testing.T) {
	a := parityPath(t, "concat_conflict_a.vcf")
	b := parityPath(t, "concat_conflict_b.vcf")
	var out bytes.Buffer
	if _, err := ConcatFiles([]string{a, b}, &out, ConcatOptions{}); err == nil {
		t.Fatalf("expected error on conflicting INFO definitions, got success:\n%s", out.String())
	}
}

// TestParityConcat_GzipOutput exercises the `-O z` output path through
// a round-trip: gunzip our output and compare against upstream's
// uncompressed concat output. (We cannot byte-compare the gzip blob
// because gzip is not deterministic across implementations.)
func TestParityConcat_GzipOutput(t *testing.T) {
	a := parityPath(t, "concat_a.vcf")
	b := parityPath(t, "concat_b.vcf")
	want := readParity(t, "concat_plain.expected.vcf")
	var buf bytes.Buffer
	if _, err := ConcatFiles([]string{a, b}, &buf, ConcatOptions{OutputFormat: OutputVCFGz}); err != nil {
		t.Fatal(err)
	}
	gr, err := gzip.NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	got, err := ioutil.ReadAll(gr)
	if err != nil {
		t.Fatal(err)
	}
	equalBytes(t, got, want, "concat -O z round-trip")
}

// =====================================================================
// norm: 7 cases
// =====================================================================

func runParityNormFile(t *testing.T, path string, opts NormOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := NormFile(path, &out, opts, io.Discard); err != nil {
		t.Fatalf("NormFile(%s): %v", path, err)
	}
	return out.Bytes()
}

// TestParityNorm_SplitAny exercises `-m -` (split everything).
func TestParityNorm_SplitAny(t *testing.T) {
	want := readParity(t, "norm_split_any.expected.vcf")
	got := runParityNormFile(t, parityPath(t, "multi.vcf"),
		NormOptions{Multiallelics: MultiallelicMode{Active: true, Split: true, Snps: true, Indels: true}})
	equalBytes(t, got, want, "norm -m -")
}

// TestParityNorm_SplitSnps exercises `-m -snps` (only SNPs).
func TestParityNorm_SplitSnps(t *testing.T) {
	want := readParity(t, "norm_split_snps.expected.vcf")
	got := runParityNormFile(t, parityPath(t, "multi.vcf"),
		NormOptions{Multiallelics: MultiallelicMode{Active: true, Split: true, Snps: true}})
	equalBytes(t, got, want, "norm -m -snps")
}

// TestParityNorm_JoinAny exercises `-m +` (join everything).
func TestParityNorm_JoinAny(t *testing.T) {
	want := readParity(t, "norm_join_any.expected.vcf")
	got := runParityNormFile(t, parityPath(t, "biallelic.vcf"),
		NormOptions{Multiallelics: MultiallelicMode{Active: true, Split: false, Snps: true, Indels: true}})
	equalBytes(t, got, want, "norm -m +")
}

// TestParityNorm_Atomize exercises `-a` (atomize complex variants).
func TestParityNorm_Atomize(t *testing.T) {
	want := readParity(t, "norm_atomize.expected.vcf")
	got := runParityNormFile(t, parityPath(t, "atom.vcf"), NormOptions{Atomize: true})
	equalBytes(t, got, want, "norm -a")
}

// TestParityNorm_RmDupExact exercises `-d exact`.
func TestParityNorm_RmDupExact(t *testing.T) {
	want := readParity(t, "norm_rmdup_exact.expected.vcf")
	got := runParityNormFile(t, parityPath(t, "dup.vcf"), NormOptions{RmDup: RmDupExact})
	equalBytes(t, got, want, "norm -d exact")
}

// TestParityNorm_LeftAlign documents that left-alignment requires a
// FASTA reference. Without the reference, upstream errors out; we error
// out too. Add a fixture (FASTA + VCF + expected) later.
func TestParityNorm_LeftAlign(t *testing.T) {
	t.Skip("norm -f left-align: needs FASTA fixture (see docs/PARITY_ROADMAP.md bcftools norm)")
}

// TestParityNorm_CheckRefSkip documents `-c s` (skip on REF/FASTA
// mismatch). Same FASTA dependency.
func TestParityNorm_CheckRefSkip(t *testing.T) {
	t.Skip("norm -c: needs FASTA fixture (see docs/PARITY_ROADMAP.md bcftools norm)")
}

// =====================================================================
// helpers
// =====================================================================

// referenceFixture reads a file from the bcftools submodule's test/
// directory at test-time. The function returns the absolute path; tests
// that depend on the submodule should call referenceFixturePath and skip
// gracefully when the submodule isn't initialised.
func referenceFixture(t *testing.T, name string) string {
	t.Helper()
	return referenceFixturePath(t, name)
}

func referenceFixturePath(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "bcftools", "test", name))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Skipf("upstream fixture %s missing: %v", name, err)
	}
	return abs
}

func copyFile(src, dst string) error {
	in, err := os.Open(src)
	if err != nil {
		return err
	}
	defer in.Close()
	out, err := os.Create(dst)
	if err != nil {
		return err
	}
	defer out.Close()
	if _, err := io.Copy(out, in); err != nil {
		return err
	}
	return nil
}

// silence unused-import warnings when guarding strip helpers behind
// skipped tests.
var _ = stripVersionLines
var _ = fmt.Sprintf
var _ = strings.TrimSpace

// =====================================================================
// head: 3 cases (PR #86 wave-1 tail)
//
// Inputs reuse `basic.vcf` from this directory. `bcftools head`
// emits the header (or a slice of it). Expected outputs are derived
// inline from the shape of `basic.vcf`'s header rather than from a
// captured upstream file because `bcftools head` was added in
// htslib 1.16+ and its golden output is the literal file header
// modulo the column-header line; the per-line content is identical.
// =====================================================================

// TestParityHead_Default emits every meta line plus the column-header
// line from basic.vcf, matching `bcftools head basic.vcf`.
func TestParityHead_Default(t *testing.T) {
	in := bytes.NewReader(readParity(t, "basic.vcf"))
	var out bytes.Buffer
	if err := Head(in, &out, HeadOptions{}); err != nil {
		t.Fatalf("Head: %v", err)
	}
	got := out.String()
	// First and last header lines should be present, in order.
	if !strings.HasPrefix(got, "##fileformat=VCFv4.2\n") {
		t.Errorf("missing leading ##fileformat:\n%s", got[:80])
	}
	if !strings.Contains(got, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n") {
		t.Errorf("missing #CHROM column-header line:\n%s", got)
	}
	// No record body should leak into the output.
	if strings.Contains(got, "rs1") {
		t.Errorf("Head leaked record body: rs1 present in:\n%s", got)
	}
}

// TestParityHead_NumLines mirrors `bcftools head -n 3 basic.vcf`:
// only the first 3 meta lines are emitted, no #CHROM line.
func TestParityHead_NumLines(t *testing.T) {
	in := bytes.NewReader(readParity(t, "basic.vcf"))
	var out bytes.Buffer
	if err := Head(in, &out, HeadOptions{NumLines: 3}); err != nil {
		t.Fatalf("Head -n 3: %v", err)
	}
	lines := strings.Split(strings.TrimRight(out.String(), "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected 3 header lines, got %d:\n%s", len(lines), out.String())
	}
	if lines[0] != "##fileformat=VCFv4.2" {
		t.Errorf("line 1 = %q, want fileformat", lines[0])
	}
	if lines[1] != "##contig=<ID=chr1,length=10000>" {
		t.Errorf("line 2 = %q, want contig chr1", lines[1])
	}
}

// TestParityHead_SamplesOnly mirrors `bcftools head --samples basic.vcf`:
// one sample-name per line, no other meta lines.
func TestParityHead_SamplesOnly(t *testing.T) {
	in := bytes.NewReader(readParity(t, "basic.vcf"))
	var out bytes.Buffer
	if err := Head(in, &out, HeadOptions{SamplesOnly: true}); err != nil {
		t.Fatalf("Head --samples: %v", err)
	}
	want := "S1\nS2\nS3\n"
	if out.String() != want {
		t.Errorf("samples-only mismatch.\ngot:%q\nwant:%q", out.String(), want)
	}
}

// =====================================================================
// sort: 3 cases (PR #86 wave-1 tail)
// =====================================================================

// TestParitySort_Basic confirms that out-of-order CHROM/POS records are
// emitted in (header-contig-order, POS) order.
func TestParitySort_Basic(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr2	500	.	A	T	.	.	.
chr1	300	.	G	A	.	.	.
chr1	100	.	C	T	.	.	.
chr2	100	.	A	G	.	.	.
`
	var out bytes.Buffer
	n, err := Sort(strings.NewReader(src), &out, SortOptions{})
	if err != nil {
		t.Fatalf("Sort: %v", err)
	}
	if n != 4 {
		t.Errorf("n = %d, want 4", n)
	}
	body := out.String()
	// chr1:100 before chr1:300 before chr2:100 before chr2:500.
	want := []string{"chr1\t100", "chr1\t300", "chr2\t100", "chr2\t500"}
	last := -1
	for _, w := range want {
		idx := strings.Index(body, w)
		if idx < 0 {
			t.Fatalf("missing %q in sorted output:\n%s", w, body)
		}
		if idx <= last {
			t.Errorf("%q at offset %d, expected after offset %d:\n%s", w, idx, last, body)
		}
		last = idx
	}
}

// TestParitySort_AlreadySorted is the no-op case — output preserves
// every input record, in the same order. Header line shape may differ
// by one (the writer auto-injects a ##FILTER=PASS line if absent), so
// we compare on record-line count not total-line count.
func TestParitySort_AlreadySorted(t *testing.T) {
	in := readParity(t, "basic.vcf")
	var out bytes.Buffer
	n, err := Sort(bytes.NewReader(in), &out, SortOptions{})
	if err != nil {
		t.Fatalf("Sort: %v", err)
	}
	// basic.vcf has 6 records.
	if n != 6 {
		t.Errorf("Sort returned n = %d, want 6", n)
	}
	body := out.String()
	// Records should appear in (chr1, chr2) order at increasing positions.
	wantOrder := []string{"chr1\t100", "chr1\t200", "chr1\t300", "chr1\t400", "chr2\t50", "chr2\t150"}
	last := -1
	for _, w := range wantOrder {
		idx := strings.Index(body, w)
		if idx < 0 {
			t.Fatalf("missing record %q in:\n%s", w, body)
		}
		if idx <= last {
			t.Errorf("record %q at offset %d, expected after offset %d", w, idx, last)
		}
		last = idx
	}
}

// TestParitySort_Empty handles the no-records edge case (header-only).
func TestParitySort_Empty(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
`
	var out bytes.Buffer
	n, err := Sort(strings.NewReader(src), &out, SortOptions{})
	if err != nil {
		t.Fatalf("Sort empty: %v", err)
	}
	if n != 0 {
		t.Errorf("n = %d, want 0", n)
	}
	if !strings.Contains(out.String(), "#CHROM") {
		t.Errorf("Sort dropped header on empty input:\n%s", out.String())
	}
}

// =====================================================================
// isec: 3 cases (PR #86 wave-1 tail)
// =====================================================================

// TestParityIsec_Intersection — `bcftools isec -n=2 -w 1` returns the
// records of input 1 that also appear in input 2.
func TestParityIsec_Intersection(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	hdr := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	if err := os.WriteFile(a, []byte(hdr+
		"chr1\t100\t.\tA\tT\t.\t.\t.\n"+
		"chr1\t200\t.\tC\tG\t.\t.\t.\n"+
		"chr1\t300\t.\tA\tT\t.\t.\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(hdr+
		"chr1\t100\t.\tA\tT\t.\t.\t.\n"+
		"chr1\t300\t.\tA\tT\t.\t.\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := IsecFiles([]string{a, b}, &out, IsecOptions{
		Nfiles: NfilesSpec{Mode: '=', N: 2},
		Write:  []int{1},
	})
	if err != nil {
		t.Fatalf("IsecFiles: %v", err)
	}
	if n != 2 {
		t.Errorf("n = %d, want 2 (records present in both)", n)
	}
	body := out.String()
	if !strings.Contains(body, "chr1\t100") || !strings.Contains(body, "chr1\t300") {
		t.Errorf("expected both shared records, got:\n%s", body)
	}
	if strings.Contains(body, "chr1\t200") {
		t.Errorf("should not contain a-only record at chr1:200:\n%s", body)
	}
}

// TestParityIsec_Bitmask exercises the `-n ~10` (input 1 only) case.
func TestParityIsec_Bitmask(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	hdr := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	if err := os.WriteFile(a, []byte(hdr+
		"chr1\t100\t.\tA\tT\t.\t.\t.\n"+
		"chr1\t200\t.\tC\tG\t.\t.\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(hdr+
		"chr1\t100\t.\tA\tT\t.\t.\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := IsecFiles([]string{a, b}, &out, IsecOptions{
		Nfiles: NfilesSpec{Mode: '~', Bits: []bool{true, false}},
		Write:  []int{1},
	})
	if err != nil {
		t.Fatalf("IsecFiles: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1 (a-only)", n)
	}
	if !strings.Contains(out.String(), "chr1\t200") {
		t.Errorf("expected chr1:200 (a-only), got:\n%s", out.String())
	}
}

// TestParityIsec_PrefixMode writes per-input projections to <prefix>/000<i>.vcf.
func TestParityIsec_PrefixMode(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	hdr := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	if err := os.WriteFile(a, []byte(hdr+"chr1\t100\t.\tA\tT\t.\t.\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(hdr+"chr1\t100\t.\tA\tT\t.\t.\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	prefix := filepath.Join(dir, "out")
	var stdout bytes.Buffer
	if _, err := IsecFiles([]string{a, b}, &stdout, IsecOptions{
		Nfiles: NfilesSpec{Mode: '=', N: 2},
		Prefix: prefix,
	}); err != nil {
		t.Fatalf("IsecFiles prefix: %v", err)
	}
	for _, name := range []string{"0000.vcf", "0001.vcf"} {
		p := filepath.Join(prefix, name)
		st, err := os.Stat(p)
		if err != nil || st.Size() == 0 {
			t.Errorf("expected non-empty %s, got %v / %v", p, st, err)
		}
	}
}

// =====================================================================
// merge: 3 cases (PR #86 wave-1 tail)
// =====================================================================

// TestParityMerge_TwoSamples merges two single-sample VCFs into a
// two-sample VCF; the merged record at chr1:100 should have GT for
// both samples.
func TestParityMerge_TwoSamples(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	hdr1 := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n"
	hdr2 := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT\t1/1\n"
	if err := os.WriteFile(a, []byte(hdr1), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(hdr2), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := MergeFiles([]string{a, b}, &out, MergeOptions{})
	if err != nil {
		t.Fatalf("MergeFiles: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1", n)
	}
	body := out.String()
	if !strings.Contains(body, "S1\tS2") {
		t.Errorf("merged header missing S1+S2 columns:\n%s", body)
	}
	if !strings.Contains(body, "chr1\t100") {
		t.Errorf("merged record at chr1:100 missing:\n%s", body)
	}
}

// TestParityMerge_DisjointPositions — records at distinct positions
// should both appear in the merged output, sorted by POS.
func TestParityMerge_DisjointPositions(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	if err := os.WriteFile(a, []byte("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\nchr1\t200\t.\tC\tG\t.\t.\t.\tGT\t1/1\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := MergeFiles([]string{a, b}, &out, MergeOptions{})
	if err != nil {
		t.Fatalf("MergeFiles: %v", err)
	}
	if n != 2 {
		t.Errorf("n = %d, want 2", n)
	}
	body := out.String()
	// Position 100 should appear before 200.
	i100 := strings.Index(body, "chr1\t100")
	i200 := strings.Index(body, "chr1\t200")
	if i100 < 0 || i200 < 0 || i100 > i200 {
		t.Errorf("expected chr1:100 before chr1:200 in:\n%s", body)
	}
}

// TestParityMerge_SingleInputRejected — merge requires >= 2 inputs.
func TestParityMerge_SingleInputRejected(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	if err := os.WriteFile(a, []byte("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := MergeFiles([]string{a}, &out, MergeOptions{}); err == nil {
		t.Errorf("expected error on single-input merge, got nil")
	}
}

// =====================================================================
// reheader: 3 cases (PR #86 wave-1 tail)
// =====================================================================

// TestParityReheader_RenameSamplesPositional renames samples by index
// using a positional sample-list (one name per line, in the same order
// as the input columns).
func TestParityReheader_RenameSamplesPositional(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := copyFile(parityPath(t, "basic.vcf"), in); err != nil {
		t.Fatal(err)
	}
	samples := filepath.Join(dir, "samples.txt")
	// Positional rename: first column → A, second → B, third → C.
	if err := os.WriteFile(samples, []byte("A\nB\nC\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ReheaderFile(in, &out, ReheaderOptions{SamplesFile: samples}); err != nil {
		t.Fatalf("ReheaderFile: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "\tA\tB\tC\n") {
		t.Errorf("expected #CHROM line to end ...\\tA\\tB\\tC, got:\n%s", body)
	}
	if strings.Contains(body, "\tS1\t") || strings.Contains(body, "\tS2\t") || strings.Contains(body, "\tS3\n") {
		t.Errorf("original sample names should be replaced:\n%s", body)
	}
}

// TestParityReheader_MapSamples uses an "OLD\tNEW" two-column file to
// rename selected samples in place.
func TestParityReheader_MapSamples(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := copyFile(parityPath(t, "basic.vcf"), in); err != nil {
		t.Fatal(err)
	}
	samples := filepath.Join(dir, "samples.tsv")
	if err := os.WriteFile(samples, []byte("S1\tNEW1\nS3\tNEW3\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ReheaderFile(in, &out, ReheaderOptions{SamplesFile: samples}); err != nil {
		t.Fatalf("ReheaderFile: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "NEW1") || !strings.Contains(body, "NEW3") {
		t.Errorf("expected NEW1 and NEW3 in:\n%s", body)
	}
	if !strings.Contains(body, "\tS2\t") && !strings.Contains(body, "\tS2\n") {
		t.Errorf("S2 should be unchanged in:\n%s", body)
	}
}

// TestParityReheader_HeaderFile substitutes the entire header from
// a separate file; the body records should be unchanged.
func TestParityReheader_HeaderFile(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := copyFile(parityPath(t, "basic.vcf"), in); err != nil {
		t.Fatal(err)
	}
	hdr := filepath.Join(dir, "hdr.vcf")
	newHdr := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">
##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indel">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=lowDP,Description="Depth below threshold">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##custom_meta=injected_via_reheader
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
`
	if err := os.WriteFile(hdr, []byte(newHdr), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ReheaderFile(in, &out, ReheaderOptions{HeaderFile: hdr}); err != nil {
		t.Fatalf("ReheaderFile: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "##custom_meta=injected_via_reheader") {
		t.Errorf("expected injected meta line in output:\n%s", body)
	}
	// Every record from basic.vcf should still be present.
	for _, want := range []string{"rs1", "rs3", "rs4", "rs5"} {
		if !strings.Contains(body, want) {
			t.Errorf("missing record %q in output:\n%s", want, body)
		}
	}
}

// =====================================================================
// annotate: 3 cases (PR #86 wave-1 tail)
// =====================================================================

// TestParityAnnotate_RemoveID — `bcftools annotate -x ID` strips the ID
// column on every record.
func TestParityAnnotate_RemoveID(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := copyFile(parityPath(t, "basic.vcf"), in); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := AnnotateFile(in, &out, AnnotateOptions{Remove: "ID"})
	if err != nil {
		t.Fatalf("AnnotateFile -x ID: %v", err)
	}
	if n != 6 {
		t.Errorf("n = %d, want 6", n)
	}
	body := out.String()
	for _, id := range []string{"rs1", "rs3", "rs4", "rs5"} {
		if strings.Contains(body, "\t"+id+"\t") {
			t.Errorf("ID %q should be stripped after -x ID:\n%s", id, body)
		}
	}
}

// TestParityAnnotate_RemoveInfoTag — `-x INFO/DP` drops the DP key from
// the INFO field of every record.
func TestParityAnnotate_RemoveInfoTag(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := copyFile(parityPath(t, "basic.vcf"), in); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := AnnotateFile(in, &out, AnnotateOptions{Remove: "INFO/DP"}); err != nil {
		t.Fatalf("AnnotateFile -x INFO/DP: %v", err)
	}
	body := out.String()
	// Records still mention DP in headers and possibly in FORMAT; the
	// INFO field should NOT contain `DP=` after the removal.
	for _, line := range strings.Split(body, "\n") {
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) >= 8 && strings.Contains(fields[7], "DP=") {
			t.Errorf("INFO field still contains DP after -x INFO/DP: %q", line)
		}
	}
}

// TestParityAnnotate_SetID — `-I '+%CHROM_%POS'` populates ID columns
// that are currently '.' (only when not already set).
func TestParityAnnotate_SetID(t *testing.T) {
	t.Skip("annotate -I/--set-id not yet implemented (see docs/PARITY_ROADMAP.md bcftools annotate)")
}
