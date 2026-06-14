package bcftools

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

const sampleVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=200>
##contig=<ID=chr2,length=300>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FILTER=<ID=q10,Description="below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	rs1	A	T	30	PASS	DP=50;AF=0.25	GT:DP	0/0:10	0/1:20	1/1:30
chr1	200	.	C	G	10	q10	DP=15;AF=0.05	GT:DP	0/0:5	0/0:5	0/1:5
chr2	50	rs2	G	A	50	PASS	DP=80;AF=0.5	GT:DP	1/1:40	0/1:25	0/0:15
`

func runView(t *testing.T, input string, opts ViewOptions) string {
	t.Helper()
	var out bytes.Buffer
	if _, err := View(strings.NewReader(input), &out, opts); err != nil {
		t.Fatalf("View: %v", err)
	}
	return out.String()
}

// recordsOf counts data lines (non-header, non-empty) in a VCF string.
func recordsOf(out string) int {
	n := 0
	for _, line := range strings.Split(out, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		n++
	}
	return n
}

func TestViewPassThrough(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{})
	if !strings.Contains(out, "rs1") || !strings.Contains(out, "rs2") {
		t.Fatalf("expected both rs1 and rs2 in output:\n%s", out)
	}
	if !strings.HasPrefix(out, "##fileformat=") {
		t.Fatalf("expected header at top:\n%s", out)
	}
}

func TestViewHeaderOnly(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{HeaderOnly: true})
	if strings.Contains(out, "rs1") {
		t.Fatalf("header-only output should not have records:\n%s", out)
	}
	if !strings.Contains(out, "#CHROM") {
		t.Fatalf("missing column header:\n%s", out)
	}
}

func TestViewNoHeader(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{NoHeader: true})
	if strings.Contains(out, "##fileformat") {
		t.Fatalf("expected no header lines:\n%s", out)
	}
	if !strings.Contains(out, "rs1") {
		t.Fatalf("records should still be there")
	}
}

func TestViewDropGenotypes(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{DropGenotypes: true})
	for _, name := range []string{"S1", "S2", "S3", "GT:DP"} {
		if strings.Contains(out, name) {
			t.Errorf("dropped output should not contain %q:\n%s", name, out)
		}
	}
}

func TestViewApplyFilterPASS(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{ApplyFilters: []string{"PASS"}})
	if strings.Contains(out, "\tq10\t") {
		t.Fatalf("PASS filter should drop q10 record:\n%s", out)
	}
	if recordsOf(out) != 2 {
		t.Fatalf("expected 2 records, got %d:\n%s", recordsOf(out), out)
	}
}

func TestViewSamplesRestrict(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{Samples: []string{"S1", "S3"}})
	if strings.Contains(out, "\tS2\t") || strings.Contains(out, "\tS2\n") {
		t.Fatalf("S2 should not appear:\n%s", out)
	}
	if !strings.Contains(out, "\tS1\tS3") {
		t.Fatalf("S1, S3 should appear:\n%s", out)
	}
}

func TestViewIncludeExpression(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{IncludeExpr: `INFO/DP>30 && FILTER="PASS"`})
	// rs1 has DP=50, rs2 has DP=80, q10 record has DP=15 → only rs1 and rs2.
	if recordsOf(out) != 2 {
		t.Fatalf("expected 2 records, got %d:\n%s", recordsOf(out), out)
	}
	if !strings.Contains(out, "rs1") || !strings.Contains(out, "rs2") {
		t.Fatalf("expected rs1 and rs2 in output:\n%s", out)
	}
}

func TestViewExcludeExpression(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{ExcludeExpr: `FILTER="q10"`})
	if strings.Contains(out, "\tq10\t") {
		t.Fatalf("q10 record should be excluded:\n%s", out)
	}
}

func TestViewMinAC(t *testing.T) {
	// AC counts: rs1 = (0,0,0,1,1,1) → 3, q10 row = (0,0,0,0,0,1) → 1, rs2 = (1,1,0,1,0,0) → 3
	out := runView(t, sampleVCF, ViewOptions{MinAlleleCount: 2})
	if recordsOf(out) != 2 {
		t.Fatalf("expected 2 records, got %d:\n%s", recordsOf(out), out)
	}
}

func TestViewMaxAC(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{MaxAlleleCount: 1})
	if recordsOf(out) != 1 {
		t.Fatalf("expected 1 record (the q10 row), got %d:\n%s", recordsOf(out), out)
	}
	if strings.Contains(out, "rs1\t") {
		t.Fatalf("rs1 (AC=3) should be filtered out:\n%s", out)
	}
}

func TestViewMinAF(t *testing.T) {
	// rs1 AF = 3/6 = 0.5, q10 row AF = 1/6 ≈ 0.167, rs2 AF = 3/6 = 0.5
	out := runView(t, sampleVCF, ViewOptions{MinAlleleFreq: 0.3})
	if recordsOf(out) != 2 {
		t.Fatalf("expected 2 records, got %d:\n%s", recordsOf(out), out)
	}
}

func TestViewMaxAF(t *testing.T) {
	out := runView(t, sampleVCF, ViewOptions{MaxAlleleFreq: 0.3})
	if recordsOf(out) != 1 {
		t.Fatalf("expected 1 record (the q10 row), got %d:\n%s", recordsOf(out), out)
	}
	if strings.Contains(out, "rs1\t") {
		t.Fatalf("rs1 (AF=0.5) should be filtered out:\n%s", out)
	}
}

func TestViewTargetsPostFilter(t *testing.T) {
	// Write the VCF to disk so we exercise ViewFile (which routes -t via the
	// streaming path with targets evaluated as a post-filter).
	dir := t.TempDir()
	path := filepath.Join(dir, "x.vcf")
	if err := os.WriteFile(path, []byte(sampleVCF), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ViewFile(path, &out, ViewOptions{Targets: []string{"chr1:90-150"}}, io.Discard); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "rs1") {
		t.Fatalf("rs1 should be kept:\n%s", out.String())
	}
	if strings.Contains(out.String(), "rs2") {
		t.Fatalf("rs2 should be excluded:\n%s", out.String())
	}
}

// TestViewRegionsWithTabix exercises the indexed-region fast path. We need a
// bgzipped VCF with a sibling .tbi for this to apply.
func TestViewRegionsWithTabix(t *testing.T) {
	dir := t.TempDir()
	bgzPath := filepath.Join(dir, "x.vcf.gz")
	f, err := os.Create(bgzPath)
	if err != nil {
		t.Fatal(err)
	}
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(sampleVCF)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()

	cfg, _ := tabix.PresetConfig(tabix.PresetVCF)
	idx, err := tabix.Build(bgzPath, cfg)
	if err != nil {
		t.Fatalf("tabix.Build: %v", err)
	}
	if err := idx.WriteFile(bgzPath + ".tbi"); err != nil {
		t.Fatalf("tabix.WriteFile: %v", err)
	}

	var out bytes.Buffer
	if _, err := ViewFile(bgzPath, &out, ViewOptions{Regions: []string{"chr1:90-150"}}, io.Discard); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "rs1") {
		t.Fatalf("rs1 should be in region output:\n%s", out.String())
	}
	if strings.Contains(out.String(), "rs2") {
		t.Fatalf("rs2 (chr2) should not be in chr1 region:\n%s", out.String())
	}
}

func TestParseOutputFormat(t *testing.T) {
	cases := []struct {
		in   string
		want OutputFormat
	}{
		{"", OutputVCF},
		{"v", OutputVCF},
		{"V", OutputVCF},
		{"z", OutputVCFGz},
		{"b", OutputBCF},
		{"u", OutputBCFUncompressed},
	}
	for _, tc := range cases {
		got, err := ParseOutputFormat(tc.in)
		if err != nil {
			t.Fatalf("%q: %v", tc.in, err)
		}
		if got != tc.want {
			t.Errorf("%q: got %d, want %d", tc.in, got, tc.want)
		}
	}
	if _, err := ParseOutputFormat("xyz"); err == nil {
		t.Fatal("expected error for unknown format")
	}
}

func TestViewVCFGzOutput(t *testing.T) {
	var out bytes.Buffer
	if _, err := View(strings.NewReader(sampleVCF), &out, ViewOptions{OutputFormat: OutputVCFGz}); err != nil {
		t.Fatal(err)
	}
	gr, err := gzip.NewReader(&out)
	if err != nil {
		t.Fatalf("output is not valid gzip: %v", err)
	}
	defer gr.Close()
	data, err := io.ReadAll(gr)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(data), "rs1") {
		t.Fatalf("missing record in gzip output:\n%s", data)
	}
}

// TestViewBCFOutputRoundTrip writes VCF to BCF and back, checking that the
// records survive a -O b → -O v round-trip via the on-disk fast path.
func TestViewBCFOutputRoundTrip(t *testing.T) {
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(sampleVCF), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatalf("View(-O b): %v", err)
	}
	f.Close()

	var back bytes.Buffer
	if _, err := ViewFile(bcfPath, &back, ViewOptions{}, io.Discard); err != nil {
		t.Fatalf("ViewFile(bcf→vcf): %v", err)
	}
	got := back.String()
	for _, want := range []string{"rs1", "rs2", "DP=50", "GT:DP", "S1\tS2\tS3"} {
		if !strings.Contains(got, want) {
			t.Errorf("BCF round-trip missing %q in:\n%s", want, got)
		}
	}
}

// TestViewBCFUncompressedRoundTrip exercises -O u (uncompressed BCF).
func TestViewBCFUncompressedRoundTrip(t *testing.T) {
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.ubcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(sampleVCF), f, ViewOptions{OutputFormat: OutputBCFUncompressed}); err != nil {
		t.Fatalf("View(-O u): %v", err)
	}
	f.Close()

	var back bytes.Buffer
	if _, err := ViewFile(bcfPath, &back, ViewOptions{}, io.Discard); err != nil {
		t.Fatalf("ViewFile(ubcf→vcf): %v", err)
	}
	if !strings.Contains(back.String(), "rs1") {
		t.Fatalf("ubcf round-trip missing rs1:\n%s", back.String())
	}
}

func TestViewBCFInput(t *testing.T) {
	// Build a minimal BCF stream in memory and feed it to View.
	stream := buildBCFFixture(t)
	var out bytes.Buffer
	if _, err := View(bytes.NewReader(stream), &out, ViewOptions{}); err != nil {
		t.Fatalf("View(BCF): %v", err)
	}
	if !strings.Contains(out.String(), "rs7") {
		t.Fatalf("expected rs7 in BCF→VCF output:\n%s", out.String())
	}
}

// buildBCFFixture constructs a one-record BCF byte stream that the BCF
// decoder accepts and ToVariant prints as "chr1 100 rs7 A T 30 PASS DP=50".
func buildBCFFixture(t *testing.T) []byte {
	t.Helper()
	const text = "##fileformat=VCFv4.3\n" +
		"##contig=<ID=chr1,length=200>\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	full := text + "\x00"
	var buf bytes.Buffer
	buf.Write(bcf.Magic[:])
	binary.Write(&buf, binary.LittleEndian, uint32(len(full)))
	buf.WriteString(full)

	// Build the shared portion by hand.
	var shared bytes.Buffer
	binary.Write(&shared, binary.LittleEndian, int32(0))  // chrom id 0 = chr1
	binary.Write(&shared, binary.LittleEndian, int32(99)) // pos
	binary.Write(&shared, binary.LittleEndian, int32(1))  // rlen
	binary.Write(&shared, binary.LittleEndian, math.Float32bits(30))
	binary.Write(&shared, binary.LittleEndian, uint16(1)) // n_info
	binary.Write(&shared, binary.LittleEndian, uint16(2)) // n_allele
	binary.Write(&shared, binary.LittleEndian, uint32(0)) // n_sample=0, n_fmt=0
	shared.Write(bcf.EncodeTypedString("rs7"))
	shared.Write(bcf.EncodeTypedString("A"))
	shared.Write(bcf.EncodeTypedString("T"))
	shared.Write(bcf.EncodeTypedInt8(0)) // FILTER = PASS (dict idx 0)
	shared.Write(bcf.EncodeTypedInt8(2)) // INFO key 2 = DP
	shared.Write(bcf.EncodeTypedInt8(50))

	binary.Write(&buf, binary.LittleEndian, uint32(shared.Len()))
	binary.Write(&buf, binary.LittleEndian, uint32(0))
	buf.Write(shared.Bytes())
	return buf.Bytes()
}

func TestLoadSamplesFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "samples.txt")
	content := "# comment\n\nS1\nS2\tignored\n# another\nS3 trailing\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	got, err := LoadSamplesFile(path)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"S1", "S2", "S3"}
	if len(got) != len(want) {
		t.Fatalf("got %v want %v", got, want)
	}
	for i := range got {
		if got[i] != want[i] {
			t.Errorf("got[%d]=%q want %q", i, got[i], want[i])
		}
	}
	if _, err := LoadSamplesFile(filepath.Join(dir, "missing")); err == nil {
		t.Fatal("expected error for missing file")
	}
}

func TestLoadRegionsFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "regions.bed")
	content := "# header\nchr1\t99\t200\nchr2\t10\t20\nchr3\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	got, err := LoadRegionsFile(path)
	if err != nil {
		t.Fatal(err)
	}
	// BED 99-200 is 0-based half-open -> 1-based inclusive 100-200.
	want := []string{"chr1:100-200", "chr2:11-20", "chr3"}
	if len(got) != len(want) {
		t.Fatalf("got %v want %v", got, want)
	}
	for i := range got {
		if got[i] != want[i] {
			t.Errorf("got[%d]=%q want %q", i, got[i], want[i])
		}
	}
}

func TestLoadRegionsFileBad(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "bad.bed")
	if err := os.WriteFile(path, []byte("chr1\toops\t10\n"), 0644); err != nil {
		t.Fatal(err)
	}
	if _, err := LoadRegionsFile(path); err == nil {
		t.Fatal("expected error for non-numeric start")
	}
}

func TestParseRegionsErrors(t *testing.T) {
	if _, err := parseRegions([]string{"chr1:abc"}); err == nil {
		t.Fatal("expected error")
	}
	if _, err := parseRegions([]string{"chr1:1-xyz"}); err == nil {
		t.Fatal("expected error")
	}
}

func TestParseRegionsShapes(t *testing.T) {
	cases := []struct {
		in        string
		wantChrom string
		wantBeg   int
		wantEnd   int
	}{
		{"chr1", "chr1", 1, 1 << 30},
		{"chr1:100", "chr1", 100, 100},
		{"chr1:50-200", "chr1", 50, 200},
		{"chr1:50-", "chr1", 50, 1 << 30},
	}
	for _, tc := range cases {
		got, err := parseRegions([]string{tc.in})
		if err != nil {
			t.Fatal(err)
		}
		if got[0].chrom != tc.wantChrom || got[0].beg != tc.wantBeg || got[0].end != tc.wantEnd {
			t.Errorf("%q: got %+v", tc.in, got[0])
		}
	}
}

func TestComputeACFromInfo(t *testing.T) {
	v := mkVariant()
	v.Samples = nil
	v.Info["AC"] = "3,5"
	v.Info["AN"] = "10"
	ac, an := computeAC(v)
	if ac != 8 || an != 10 {
		t.Errorf("got ac=%d an=%d", ac, an)
	}
}

func TestRestrictSamplesUnknown(t *testing.T) {
	v := mkVariant()
	v.Samples = []vcf.Sample{{Name: "S1"}, {Name: "S2"}}
	restrictSamples(v, []string{"S2", "S3"})
	if len(v.Samples) != 1 || v.Samples[0].Name != "S2" {
		t.Errorf("got %+v", v.Samples)
	}
}

// privateVCF mirrors testdata/parity/view/private.vcf: subset {S1,S2}.
//
//	priv1   S1=0/1 S2=0/0 S3=0/0 -> acSub=1 acFull=1 -> private
//	shared  S1=0/1 S2=0/0 S3=0/1 -> acSub=1 acFull=2 -> not private
//	priv2   S1=0/0 S2=1/1 S3=0/0 -> acSub=2 acFull=2 -> private
//	noalt   all 0/0              -> acSub=0          -> not private
//	outonly S1=0/0 S2=0/0 S3=0/1 -> acSub=0 acFull=1 -> not private
const privateVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	priv1	A	T	30	PASS	DP=10	GT	0/1	0/0	0/0
chr1	200	shared	C	G	30	PASS	DP=10	GT	0/1	0/0	0/1
chr1	300	priv2	G	A	30	PASS	DP=10	GT	0/0	1/1	0/0
chr1	400	noalt	T	C	30	PASS	DP=10	GT	0/0	0/0	0/0
chr1	500	outonly	A	G	30	PASS	DP=10	GT	0/0	0/0	0/1
`

// recordIDs returns the ID column of every data record in a VCF string,
// preserving order.
func recordIDs(out string) []string {
	var ids []string
	for _, line := range strings.Split(out, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) >= 3 {
			ids = append(ids, fields[2])
		}
	}
	return ids
}

func equalStrings(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// TestViewPrivate exercises -x/--private and -X/--exclude-private record
// selection against the embedded private fixture.
func TestViewPrivate(t *testing.T) {
	cases := []struct {
		name string
		opts ViewOptions
		want []string
	}{
		{
			name: "private subset S1,S2",
			opts: ViewOptions{Samples: []string{"S1", "S2"}, Private: true},
			want: []string{"priv1", "priv2"},
		},
		{
			name: "exclude-private subset S1,S2",
			opts: ViewOptions{Samples: []string{"S1", "S2"}, ExcludePrivate: true},
			want: []string{"shared", "noalt", "outonly"},
		},
		{
			name: "private subset S3",
			opts: ViewOptions{Samples: []string{"S3"}, Private: true},
			want: []string{"outonly"},
		},
		{
			name: "exclude-private subset S3",
			opts: ViewOptions{Samples: []string{"S3"}, ExcludePrivate: true},
			want: []string{"priv1", "shared", "priv2", "noalt"},
		},
		{
			// Without a sample subset the private filter is a no-op: every
			// record is kept (matches upstream gating on n_samples).
			name: "private without subset is no-op",
			opts: ViewOptions{Private: true},
			want: []string{"priv1", "shared", "priv2", "noalt", "outonly"},
		},
		{
			// ExcludePrivate is likewise a no-op without a subset: nothing is
			// dropped (upstream gates the test on n_samples > 0).
			name: "exclude-private without subset is no-op",
			opts: ViewOptions{ExcludePrivate: true},
			want: []string{"priv1", "shared", "priv2", "noalt", "outonly"},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			out := runView(t, privateVCF, tc.opts)
			got := recordIDs(out)
			if !equalStrings(got, tc.want) {
				t.Fatalf("got IDs %v, want %v\noutput:\n%s", got, tc.want, out)
			}
		})
	}
}

// TestViewPrivateNoGT verifies the private filter is a no-op when the records
// carry no GT field (upstream evaluates it only after recomputing AC from
// genotypes).
func TestViewPrivateNoGT(t *testing.T) {
	const noGT = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr1	200	b	C	G	30	PASS	DP=10
`
	out := runView(t, noGT, ViewOptions{Samples: []string{"S1"}, Private: true})
	if got := recordIDs(out); !equalStrings(got, []string{"a", "b"}) {
		t.Fatalf("private filter should be a no-op without GT; got %v\n%s", got, out)
	}
}

// dataRecordsStripINFO returns each VCF data line with the INFO column (field
// index 7) blanked. Upstream `bcftools view -s ... -x/-X` recomputes INFO
// AC/AN after subsetting (a separate, documented gap — see
// TestParityView_SampleSubset); blanking INFO lets us assert byte-for-byte
// parity on the record *selection* and every other column the private filter
// governs.
func dataRecordsStripINFO(vcfText string) []string {
	var out []string
	for _, line := range strings.Split(vcfText, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) > 7 {
			fields[7] = "."
		}
		out = append(out, strings.Join(fields, "\t"))
	}
	return out
}

func TestRecomputeACAN(t *testing.T) {
	// Two kept samples, multiallelic ALT=T,G: S1=1/2 (one each), S2=0/1.
	v := &vcf.Variant{
		Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"T", "G"},
		Info:      map[string]string{"DP": "99"},
		InfoOrder: []string{"DP"},
		Samples: []vcf.Sample{
			{Name: "S1", Data: map[string]string{"GT": "1/2"}},
			{Name: "S2", Data: map[string]string{"GT": "0/1"}},
		},
	}
	recomputeACAN(v)
	if v.Info["AC"] != "2,1" {
		t.Errorf("AC = %q, want 2,1", v.Info["AC"])
	}
	if v.Info["AN"] != "4" {
		t.Errorf("AN = %q, want 4", v.Info["AN"])
	}
	// AC and AN appended after the pre-existing DP, in that order.
	if got := strings.Join(v.InfoOrder, ","); got != "DP,AC,AN" {
		t.Errorf("InfoOrder = %q, want DP,AC,AN", got)
	}
}

func TestPassesTypeFilter(t *testing.T) {
	snp := &vcf.Variant{Ref: "A", Alt: []string{"T"}}
	ins := &vcf.Variant{Ref: "A", Alt: []string{"AT"}}
	mnp := &vcf.Variant{Ref: "AC", Alt: []string{"GT"}}

	if !(ViewOptions{IncludeTypes: []string{"snps"}}).passesTypeFilter(snp) {
		t.Error("-v snps should keep a SNP")
	}
	if (ViewOptions{IncludeTypes: []string{"snps"}}).passesTypeFilter(ins) {
		t.Error("-v snps should drop an indel")
	}
	if !(ViewOptions{IncludeTypes: []string{"indels"}}).passesTypeFilter(ins) {
		t.Error("-v indels should keep an indel")
	}
	if !(ViewOptions{IncludeTypes: []string{"mnps"}}).passesTypeFilter(mnp) {
		t.Error("-v mnps should keep an MNP")
	}
	if (ViewOptions{ExcludeTypes: []string{"snps"}}).passesTypeFilter(snp) {
		t.Error("-V snps should drop a SNP")
	}
	if !(ViewOptions{ExcludeTypes: []string{"snps"}}).passesTypeFilter(ins) {
		t.Error("-V snps should keep an indel")
	}
}
