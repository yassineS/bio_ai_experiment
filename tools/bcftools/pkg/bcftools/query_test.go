package bcftools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
	"github.com/yassineS/bio_ai_experiment/tools/tabix/pkg/tabix"
)

const queryVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=200>
##contig=<ID=chr2,length=300>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##FILTER=<ID=q10,Description="below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	rs1	A	T	30	PASS	DP=50;AF=0.25;H2	GT:DP	0/0:10	0/1:20	1/1:30
chr1	200	.	C	G	10	q10	DP=5;AF=0.05	GT:DP	0/0:5	0/0:5	0/1:5
chr2	50	rs2	G	A,C	50	PASS	DP=80;AF=0.5	GT:DP	1/2:40	0/1:25	0/0:15
chr2	60	rs3	AT	A	40	PASS	DP=70	GT:DP	0/0:1	0/1:2	1/1:3
chr2	70	rs4	AT	GC	40	PASS	DP=60	GT:DP	0/0:1	0/1:2	1/1:3
`

func runQuery(t *testing.T, input string, opts QueryOptions) string {
	t.Helper()
	var out bytes.Buffer
	if _, err := Query(strings.NewReader(input), &out, opts); err != nil {
		t.Fatalf("Query: %v", err)
	}
	return out.String()
}

func TestQueryBasicFormat(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\t%REF\t%ALT\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 5 {
		t.Fatalf("expected 5 lines, got %d:\n%s", len(lines), out)
	}
	if lines[0] != "chr1\t100\tA\tT" {
		t.Errorf("first line wrong: %q", lines[0])
	}
	if lines[2] != "chr2\t50\tG\tA,C" {
		t.Errorf("third line wrong: %q", lines[2])
	}
}

func TestQueryInfoTag(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%INFO/DP\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if lines[0] != "chr1\t50" || lines[1] != "chr1\t5" {
		t.Errorf("INFO/DP not extracted: %v", lines)
	}
}

func TestQueryInfoMissing(t *testing.T) {
	// rs3/rs4 don't carry AF — that should print "." not an empty cell.
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\t%INFO/AF\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if !strings.HasSuffix(lines[3], "\t.") {
		t.Errorf("missing INFO/AF should emit '.': %q", lines[3])
	}
}

func TestQueryInfoFlag(t *testing.T) {
	// H2 is a flag (no value). When present it should render as "1".
	out := runQuery(t, queryVCF, QueryOptions{Format: `%INFO/H2\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if lines[0] != "1" {
		t.Errorf("flag INFO should render as 1, got %q", lines[0])
	}
	if lines[1] != "." {
		t.Errorf("absent flag should render as '.', got %q", lines[1])
	}
}

func TestQuerySampleRepeated(t *testing.T) {
	// Upstream bcftools repeats the inner [...] body verbatim with no
	// auto-inserted separator between samples. The literal in front of the
	// sample placeholder is what produces the inter-sample gap (or its
	// absence).
	out := runQuery(t, queryVCF, QueryOptions{Format: `[%GT]\t[%DP]\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	want := "0/00/11/1\t102030"
	if lines[0] != want {
		t.Errorf("got %q want %q", lines[0], want)
	}
}

func TestQuerySampleNarrowing(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `[%GT]\n`, Samples: []string{"S1", "S3"}})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if lines[0] != "0/01/1" {
		t.Errorf("first row: got %q want %q", lines[0], "0/01/1")
	}
}

func TestQueryType(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\t%TYPE\n`})
	want := map[string]string{
		"chr1\t100": "SNP",
		"chr1\t200": "SNP",
		"chr2\t50":  "SNP",
		"chr2\t60":  "INDEL",
		"chr2\t70":  "MNP",
	}
	for _, line := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		fields := strings.SplitN(line, "\t", 3)
		key := fields[0] + "\t" + fields[1]
		if want[key] != fields[2] {
			t.Errorf("%s: got %q want %q", key, fields[2], want[key])
		}
	}
}

func TestQueryTGT(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `[%TGT]\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// rs2: REF=G, ALT=A,C with GT 1/2,0/1,0/0 -> A/C, G/A, G/G,
	// concatenated with no separator (upstream behaviour).
	if lines[2] != "A/CG/AG/G" {
		t.Errorf("rs2 TGT: got %q", lines[2])
	}
}

func TestQueryHeader(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\t%REF\n`, PrintHeader: true})
	first := strings.SplitN(out, "\n", 2)[0]
	if !strings.HasPrefix(first, "# ") {
		t.Errorf("header row missing '# ' prefix: %q", first)
	}
	if !strings.Contains(first, "CHROM") || !strings.Contains(first, "POS") || !strings.Contains(first, "REF") {
		t.Errorf("header row missing token names: %q", first)
	}
}

func TestQueryHeaderWithSamples(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t[%GT]\n`, PrintHeader: true})
	first := strings.SplitN(out, "\n", 2)[0]
	// Per-sample columns should carry the sample name.
	for _, want := range []string{"GT:S1", "GT:S2", "GT:S3"} {
		if !strings.Contains(first, want) {
			t.Errorf("header missing %q: %q", want, first)
		}
	}
}

func TestQueryIncludeFilter(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\n`, IncludeExpr: `INFO/DP>10`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// rs1 (DP=50), rs2 (DP=80), rs3 (DP=70), rs4 (DP=60) pass; q10 row (DP=5) drops.
	if len(lines) != 4 {
		t.Fatalf("expected 4 lines, got %d:\n%s", len(lines), out)
	}
	for _, line := range lines {
		if strings.HasPrefix(line, "chr1\t200") {
			t.Errorf("DP=5 row should be excluded: %s", line)
		}
	}
}

func TestQueryExcludeFilter(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\n`, ExcludeExpr: `FILTER="q10"`})
	if strings.Contains(out, "chr1\t200") {
		t.Errorf("q10 row should be excluded:\n%s", out)
	}
}

func TestQueryApplyFilters(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%POS\n`, ApplyFilters: []string{"PASS"}})
	if strings.Contains(out, "chr1\t200") {
		t.Errorf("non-PASS row should be excluded:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t100") {
		t.Errorf("PASS row missing:\n%s", out)
	}
}

func TestQueryListSamplesFromFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "x.vcf")
	if err := os.WriteFile(path, []byte(queryVCF), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := QueryFile(path, &out, QueryOptions{ListSamples: true}, io.Discard); err != nil {
		t.Fatal(err)
	}
	got := strings.TrimRight(out.String(), "\n")
	want := "S1\nS2\nS3"
	if got != want {
		t.Errorf("list-samples: got %q want %q", got, want)
	}
}

func TestQueryFileTargets(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "x.vcf")
	if err := os.WriteFile(path, []byte(queryVCF), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	opts := QueryOptions{Format: `%CHROM\t%POS\n`, Targets: []string{"chr1"}}
	if _, err := QueryFile(path, &out, opts, io.Discard); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "chr2") {
		t.Errorf("chr2 records should be filtered: %s", out.String())
	}
}

func TestQueryFileRegionsFallback(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "x.vcf")
	if err := os.WriteFile(path, []byte(queryVCF), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	var stderr bytes.Buffer
	opts := QueryOptions{Format: `%CHROM\t%POS\n`, Regions: []string{"chr2"}}
	if _, err := QueryFile(path, &out, opts, &stderr); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "chr1") {
		t.Errorf("chr1 records should be filtered: %s", out.String())
	}
	if !strings.Contains(stderr.String(), "no .tbi/.csi index found") {
		t.Errorf("expected fallback warning, got: %s", stderr.String())
	}
}

func TestQueryEscapes(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\\%POS\nfoo`})
	// Two records (first two lines, ignoring rest): chr1\100\nfoo... etc.
	if !strings.HasPrefix(out, `chr1\100`) {
		t.Errorf("escaped backslash not preserved: %q", out[:20])
	}
	if !strings.Contains(out, "\nfoochr1") {
		t.Errorf("literal foo missing between records: %q", out)
	}
}

func TestParseFormatStringErrors(t *testing.T) {
	cases := []string{
		``,
		`%`,
		`[%CHROM`,
		`]oops`,
		`[[%X]`,
		`%CHROM\`,
	}
	for _, c := range cases {
		if _, err := ParseFormatString(c); err == nil {
			t.Errorf("expected error for %q", c)
		}
	}
}

func TestQueryEmptyFormatErrors(t *testing.T) {
	var out bytes.Buffer
	if _, err := Query(strings.NewReader(queryVCF), &out, QueryOptions{Format: ""}); err == nil {
		t.Error("expected error for empty format string")
	}
}

func TestQueryBadIncludeExpression(t *testing.T) {
	var out bytes.Buffer
	if _, err := Query(strings.NewReader(queryVCF), &out, QueryOptions{Format: `%POS\n`, IncludeExpr: `nope (`}); err == nil {
		t.Error("expected error for malformed expression")
	}
}

func TestQueryUnknownPlaceholder(t *testing.T) {
	// Unknown placeholder is allowed by the tokenizer but should render as ".".
	out := runQuery(t, queryVCF, QueryOptions{Format: `%CHROM\t%NOSUCH\n`})
	if !strings.Contains(out, "chr1\t.\n") {
		t.Errorf("unknown placeholder should render '.': %q", out)
	}
}

func TestQueryListSamplesViaStream(t *testing.T) {
	var out bytes.Buffer
	if _, err := Query(strings.NewReader(queryVCF), &out, QueryOptions{ListSamples: true}); err != nil {
		t.Fatal(err)
	}
	want := "S1\nS2\nS3\n"
	if out.String() != want {
		t.Errorf("ListSamples streaming: got %q want %q", out.String(), want)
	}
}

func TestQueryFmtPrefix(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `[%FMT/GT]\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// No auto separator between samples (upstream-compatible behaviour).
	if lines[0] != "0/00/11/1" {
		t.Errorf("FMT/GT: got %q", lines[0])
	}
}

func TestQueryQualMissing(t *testing.T) {
	// Construct a tiny VCF with QUAL = "." and confirm formatting.
	const minimalVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	T	.	PASS	.
`
	out := runQuery(t, minimalVCF, QueryOptions{Format: `%QUAL\t%ID\t%FILTER\n`})
	if strings.TrimRight(out, "\n") != ".\t.\tPASS" {
		t.Errorf("missing QUAL/ID handling: got %q", out)
	}
}

func TestQueryTypeEdgeCases(t *testing.T) {
	// Symbolic ALT and split alleles get TYPE=OTHER per upstream.
	const svVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	<DEL>	.	PASS	.
chr1	2	.	A	*	.	PASS	.
`
	out := runQuery(t, svVCF, QueryOptions{Format: `%POS\t%TYPE\n`})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if lines[0] != "1\tOTHER" || lines[1] != "2\tOTHER" {
		t.Errorf("symbolic/star ALT TYPE: %v", lines)
	}
}

// TestQueryBCFStream feeds a BCF stream into Query and confirms the BCF
// dispatch path renders correctly.
func TestQueryBCFStream(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	rsB	A	T	30	PASS	DP=42	GT	0/1
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatalf("View(-O b): %v", err)
	}
	f.Close()
	var out bytes.Buffer
	if _, err := QueryFile(bcfPath, &out, QueryOptions{Format: `%CHROM\t%POS\t%INFO/DP\t[%GT]\n`}, io.Discard); err != nil {
		t.Fatalf("QueryFile(bcf): %v", err)
	}
	got := strings.TrimRight(out.String(), "\n")
	if got != "chr1\t100\t42\t0/1" {
		t.Errorf("BCF stream query: got %q", got)
	}
}

// TestQueryBCFListSamples exercises the list-samples path against a BCF input.
func TestQueryBCFListSamples(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	1	.	A	T	.	PASS	DP=1	GT	0/0	0/1
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	var out bytes.Buffer
	if _, err := QueryFile(bcfPath, &out, QueryOptions{ListSamples: true}, io.Discard); err != nil {
		t.Fatal(err)
	}
	if strings.TrimRight(out.String(), "\n") != "A\nB" {
		t.Errorf("list-samples (BCF): got %q", out.String())
	}
}

// TestQueryBCFRegions exercises the CSI region path for BCF.
func TestQueryBCFRegions(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr2	200	c	G	A	30	PASS	DP=30
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := QueryFile(bcfPath, &out, QueryOptions{Format: `%CHROM\t%POS\n`, Regions: []string{"chr2"}}, io.Discard); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "chr1") {
		t.Errorf("chr1 should be filtered:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "chr2\t200") {
		t.Errorf("chr2:200 missing:\n%s", out.String())
	}
}

// TestQueryVCFRegions exercises the TBI-backed region path for bgzipped VCF.
func TestQueryVCFRegions(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr2	200	c	G	A	30	PASS	DP=30
`
	dir := t.TempDir()
	gzPath := filepath.Join(dir, "x.vcf.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(vcfText)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()
	cfg, _ := tabix.PresetConfig(tabix.PresetVCF)
	idx, err := tabix.Build(gzPath, cfg)
	if err != nil {
		t.Fatalf("tabix.Build: %v", err)
	}
	if err := idx.WriteFile(gzPath + ".tbi"); err != nil {
		t.Fatalf("tabix.WriteFile: %v", err)
	}
	var out bytes.Buffer
	if _, err := QueryFile(gzPath, &out, QueryOptions{Format: `%CHROM\t%POS\n`, Regions: []string{"chr1:50-200"}}, io.Discard); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "chr2") {
		t.Errorf("chr2 should not appear:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "chr1\t100") {
		t.Errorf("chr1:100 missing:\n%s", out.String())
	}
}

func TestQueryFilterSamplesByName(t *testing.T) {
	got := filterSamplesByName([]string{"S2", "MISSING", "S1"}, []string{"S1", "S2", "S3"})
	if len(got) != 2 || got[0] != "S2" || got[1] != "S1" {
		t.Errorf("got %v", got)
	}
}

func TestQuerySampleFieldOutOfRange(t *testing.T) {
	v := &vcf.Variant{}
	got := sampleField(v, 99, "GT")
	if got != "." {
		t.Errorf("out-of-range sampleField: got %q", got)
	}
}

// TestQueryBCFRegionsHeader exercises the header / sample-filter branches of
// the BCF region path.
func TestQueryBCFRegionsHeader(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	a	A	T	30	PASS	DP=10	GT	0/1	1/1
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	opts := QueryOptions{
		Format:      `%CHROM\t%POS\t[%GT]\n`,
		PrintHeader: true,
		Regions:     []string{"chr1"},
		Samples:     []string{"B", "A"},
	}
	if _, err := QueryFile(bcfPath, &out, opts, io.Discard); err != nil {
		t.Fatal(err)
	}
	got := out.String()
	if !strings.Contains(got, "GT:B") || !strings.Contains(got, "GT:A") {
		t.Errorf("header missing sample columns:\n%s", got)
	}
	// The sample filter should reorder GTs as B then A; no auto-separator
	// between repeated samples in the inner [...] body.
	if !strings.Contains(got, "1/10/1") {
		t.Errorf("sample reorder failed:\n%s", got)
	}
}

// TestQueryBCFRegionsIncludeExpr drives the include-expression path under
// CSI regions.
func TestQueryBCFRegionsIncludeExpr(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=5
chr1	200	b	C	G	30	PASS	DP=50
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	opts := QueryOptions{Format: `%POS\n`, Regions: []string{"chr1"}, IncludeExpr: `INFO/DP>10`}
	if _, err := QueryFile(bcfPath, &out, opts, io.Discard); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "100") {
		t.Errorf("DP=5 record should be excluded:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "200") {
		t.Errorf("DP=50 record missing:\n%s", out.String())
	}
}

// TestQueryVCFRegionsHeader exercises the TBI region path with -H/-s.
func TestQueryVCFRegionsHeader(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	a	A	T	30	PASS	.	GT	0/0	1/1
chr1	200	b	C	G	30	PASS	.	GT	0/1	1/1
`
	dir := t.TempDir()
	gzPath := filepath.Join(dir, "x.vcf.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(vcfText)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()
	cfg, _ := tabix.PresetConfig(tabix.PresetVCF)
	idx, err := tabix.Build(gzPath, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := idx.WriteFile(gzPath + ".tbi"); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	opts := QueryOptions{
		Format:      `%POS\t[%GT]\n`,
		PrintHeader: true,
		Regions:     []string{"chr1:50-150"},
		Samples:     []string{"B"},
	}
	if _, err := QueryFile(gzPath, &out, opts, io.Discard); err != nil {
		t.Fatal(err)
	}
	got := out.String()
	if !strings.Contains(got, "GT:B") {
		t.Errorf("header should have GT:B column:\n%s", got)
	}
	// Region cut: only chr1:100 (POS=100) should appear.
	if strings.Contains(got, "200") {
		t.Errorf("out-of-region record leaked:\n%s", got)
	}
	if !strings.Contains(got, "100\t1/1") {
		t.Errorf("region match missing:\n%s", got)
	}
}

// TestQueryBCFStreamHeader exercises header + filter + sample-narrowing on a
// plain BCF stream (no CSI index).
func TestQueryBCFStreamHeader(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	a	A	T	30	PASS	DP=10	GT	0/0	1/1
chr1	200	b	C	G	30	PASS	DP=50	GT	0/1	1/1
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	// Open through QueryFile (no index) → streaming dispatch + filter branches.
	var out bytes.Buffer
	opts := QueryOptions{
		Format:      `%POS\t[%GT]\n`,
		PrintHeader: true,
		IncludeExpr: `INFO/DP>20`,
		Samples:     []string{"B"},
	}
	if _, err := QueryFile(bcfPath, &out, opts, io.Discard); err != nil {
		t.Fatal(err)
	}
	got := out.String()
	if !strings.Contains(got, "GT:B") {
		t.Errorf("header missing GT:B:\n%s", got)
	}
	if strings.Contains(got, "\n100\t") {
		t.Errorf("DP=10 record should be excluded:\n%s", got)
	}
	if !strings.Contains(got, "200\t1/1") {
		t.Errorf("DP=50 record missing:\n%s", got)
	}
}

// TestQueryListSamplesBCFStream covers list-samples on a stream wrapped in
// QueryFile (path-based, BCF).
func TestQueryListSamplesBCFStream(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	T	.	PASS	DP=1
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	var out bytes.Buffer
	// list-samples on a no-sample BCF should produce empty output.
	if _, err := QueryFile(bcfPath, &out, QueryOptions{ListSamples: true}, io.Discard); err != nil {
		t.Fatal(err)
	}
	if out.String() != "" {
		t.Errorf("expected empty list-samples output, got %q", out.String())
	}
}

// TestQueryPlaceholderEdgeCases covers the rarely-hit fall-through branches
// in formatPlaceholder.
func TestQueryPlaceholderEdgeCases(t *testing.T) {
	const text = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	T	30.5	.	.	GT	0/1
`
	// QUAL=30.5 hits the float-format branch; missing ID/FILTER hit the "."
	// branches; %GT outside [] returns "."; %FMT/GT outside [] returns ".".
	out := runQuery(t, text, QueryOptions{Format: `%QUAL\t%ID\t%FILTER\t%GT\t%FMT/GT\n`})
	if strings.TrimRight(out, "\n") != "30.5\t.\t.\t.\t." {
		t.Errorf("placeholder edge cases: got %q", out)
	}
}

// TestQuerySampleLiterals covers the inner-literal branch of [%...] groups.
func TestQuerySampleLiterals(t *testing.T) {
	out := runQuery(t, queryVCF, QueryOptions{Format: `[%GT=%DP]\n`, PrintHeader: true})
	// Header row has inner literals between placeholders; data rows have
	// "0/0=10\t0/1=20\t1/1=30" etc.
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if !strings.Contains(lines[1], "0/0=10") {
		t.Errorf("inner literal not emitted: %q", lines[1])
	}
	if !strings.Contains(lines[0], "=") {
		t.Errorf("header should preserve inner literal: %q", lines[0])
	}
}

// TestQueryTGTOutOfRangeIndex covers the "tok parse fails" path of TGT
// translation.
func TestQueryTGTOutOfRangeIndex(t *testing.T) {
	const text = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	T	.	PASS	.	GT	9/0
`
	out := runQuery(t, text, QueryOptions{Format: `[%TGT]\n`})
	// "9" is out of range (only REF + 1 ALT = indices 0..1) so the tokeniser
	// emits the raw index back.
	if strings.TrimRight(out, "\n") != "9/A" {
		t.Errorf("out-of-range TGT: got %q", out)
	}
}

// TestQueryAltMissing covers the no-ALT (len==0) branch.
func TestQueryAltMissing(t *testing.T) {
	// Build a Variant with zero ALT (e.g. gVCF reference block uses ".").
	// Since vcf.Read uses "." as a single-element split, we test directly:
	v := &vcf.Variant{Chrom: "x", Pos: 1, Ref: "A"}
	if got := formatPlaceholder("ALT", v, -1); got != "." {
		t.Errorf("empty ALT: got %q", got)
	}
}

// TestQueryTranslatedGenotypeMissing exercises the missing-token branch of
// translatedGenotype.
func TestQueryTranslatedGenotypeMissing(t *testing.T) {
	const minimalVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	T	.	PASS	.	GT	./.
`
	out := runQuery(t, minimalVCF, QueryOptions{Format: `[%TGT]\n`})
	// "./." should translate to "./." (missing tokens preserved).
	if strings.TrimRight(out, "\n") != "./." {
		t.Errorf("missing TGT: got %q", out)
	}
}

func TestQueryTGTPhased(t *testing.T) {
	const phasedVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	T	.	PASS	.	GT	0|1
`
	out := runQuery(t, phasedVCF, QueryOptions{Format: `[%TGT]\n`})
	if strings.TrimRight(out, "\n") != "A|T" {
		t.Errorf("phased TGT: got %q", out)
	}
}
