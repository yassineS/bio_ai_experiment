package vcftools

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// TestParseFilterList covers comma-separated FILTER-name parsing.
func TestParseFilterList(t *testing.T) {
	got := parseFilterList(" q10 , s50 , ")
	if len(got) != 2 {
		t.Errorf("got %d, want 2: %v", len(got), got)
	}
	if _, ok := got["q10"]; !ok {
		t.Error("missing q10")
	}
	if _, ok := got["s50"]; !ok {
		t.Error("missing s50")
	}
	if parseFilterList("") != nil {
		t.Error("expected nil for empty input")
	}
}

// TestPassRemoveFilteredNames tests the per-variant FILTER-name drop.
func TestPassRemoveFilteredNames(t *testing.T) {
	set := parseFilterList("q10,s50")
	cases := []struct {
		filter []string
		want   bool
	}{
		{nil, true},              // no filter -> kept
		{[]string{"PASS"}, true}, // PASS -> kept
		{[]string{"q10"}, false}, // hit -> dropped
		{[]string{"s50"}, false}, // hit -> dropped
		{[]string{"q10", "s50"}, false},
		{[]string{"other"}, true},
	}
	for _, tc := range cases {
		v := &vcf.Variant{Filter: tc.filter}
		got := passRemoveFilteredNames(v, set)
		if got != tc.want {
			t.Errorf("filter=%v: got %v want %v", tc.filter, got, tc.want)
		}
	}
}

// TestPassKeepFilteredNames tests --keep-filtered inverse behaviour.
func TestPassKeepFilteredNames(t *testing.T) {
	set := parseFilterList("q10")
	cases := []struct {
		filter []string
		want   bool
	}{
		{[]string{"PASS"}, false}, // PASS doesn't list q10 -> dropped
		{[]string{"q10"}, true},   // matches -> kept
		{[]string{"q10", "s50"}, true},
		{[]string{"s50"}, false},
		{nil, false},
	}
	for _, tc := range cases {
		v := &vcf.Variant{Filter: tc.filter}
		got := passKeepFilteredNames(v, set)
		if got != tc.want {
			t.Errorf("filter=%v: got %v want %v", tc.filter, got, tc.want)
		}
	}
}

// TestFilterRecodeInfo covers --keep-INFO / --remove-INFO restriction.
func TestFilterRecodeInfo(t *testing.T) {
	in := map[string]string{"AF": "0.5", "DP": "10", "AA": "T"}

	// Empty sets -> unchanged.
	out := filterRecodeInfo(in, nil, nil)
	if len(out) != 3 {
		t.Errorf("empty sets: got %d, want 3", len(out))
	}

	// keep-INFO restricts to listed tags.
	keep := parseInfoTagList("AF,DP")
	out = filterRecodeInfo(in, keep, nil)
	if len(out) != 2 || out["AF"] != "0.5" || out["DP"] != "10" {
		t.Errorf("keep-INFO AF,DP: got %v", out)
	}

	// remove-INFO strips listed tags.
	rem := parseInfoTagList("DP")
	out = filterRecodeInfo(in, nil, rem)
	if len(out) != 2 || out["AF"] != "0.5" || out["AA"] != "T" {
		t.Errorf("remove-INFO DP: got %v", out)
	}

	// Both compose: keep-INFO AF,DP then strip DP -> only AF survives.
	out = filterRecodeInfo(in, keep, rem)
	if len(out) != 1 || out["AF"] != "0.5" {
		t.Errorf("keep+remove: got %v", out)
	}
}

// TestRun_RemoveFiltered_Integration: --remove-filtered drops sites that
// list any of the named FILTERs.
func TestRun_RemoveFiltered_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n" +
		"1\t200\t.\tA\tG\t.\tq10\t.\tGT\t0/1\n" +
		"1\t300\t.\tA\tG\t.\ts50\t.\tGT\t1/1\n" +
		"1\t400\t.\tA\tG\t.\tq10;s50\t.\tGT\t0/1\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Recode: true, RemoveFiltered: "q10"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode: %v", err)
	}
	body := string(data)
	if strings.Contains(body, "\t200\t") {
		t.Errorf("expected 200 (q10) dropped; got:\n%s", body)
	}
	if strings.Contains(body, "\t400\t") {
		t.Errorf("expected 400 (q10;s50) dropped; got:\n%s", body)
	}
	if !strings.Contains(body, "\t100\t") || !strings.Contains(body, "\t300\t") {
		t.Errorf("expected 100 and 300 kept; got:\n%s", body)
	}
}

// TestRun_KeepFiltered_Integration: --keep-filtered keeps only sites that
// list at least one named FILTER.
func TestRun_KeepFiltered_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n" +
		"1\t200\t.\tA\tG\t.\tq10\t.\tGT\t0/1\n" +
		"1\t300\t.\tA\tG\t.\ts50\t.\tGT\t1/1\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Recode: true, KeepFiltered: "q10"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".recode.vcf")
	body := string(data)
	if !strings.Contains(body, "\t200\t") {
		t.Errorf("expected 200 (q10) kept; got:\n%s", body)
	}
	if strings.Contains(body, "\t100\t") || strings.Contains(body, "\t300\t") {
		t.Errorf("expected only 200 kept; got:\n%s", body)
	}
}

// TestRun_GetINFO_Integration extracts AF and DP into .INFO.
func TestRun_GetINFO_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tAF=0.5;DP=10\tGT\t0/0\n" +
		"1\t200\t.\tA\tC\t.\tPASS\tAF=0.25\tGT\t0/1\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, GetINFO: "AF,DP"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".INFO")
	if err != nil {
		t.Fatalf("read .INFO: %v", err)
	}
	want := "CHROM\tPOS\tREF\tALT\tAF\tDP\n" +
		"1\t100\tA\tG\t0.5\t10\n" +
		"1\t200\tA\tC\t0.25\t.\n"
	if string(data) != want {
		t.Errorf("got:\n%s\nwant:\n%s", data, want)
	}
}

// TestRun_KeepINFO_Integration restricts INFO in recoded output.
func TestRun_KeepINFO_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=AF,Number=1,Type=Float,Description=\"\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tAF=0.5;DP=10;AA=T\tGT\t0/0\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Recode: true, KeepINFO: "AF"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".recode.vcf")
	body := string(data)
	// Find the data row
	for _, line := range strings.Split(body, "\n") {
		if strings.HasPrefix(line, "1\t") {
			fields := strings.Split(line, "\t")
			if len(fields) < 8 {
				t.Fatalf("unexpected line %q", line)
			}
			info := fields[7]
			// Only AF should survive.
			if info != "AF=0.5" {
				t.Errorf("INFO = %q, want %q", info, "AF=0.5")
			}
			return
		}
	}
	t.Fatalf("data row not found in:\n%s", body)
}

// TestRun_RemoveINFO_Integration strips a tag from recoded output.
func TestRun_RemoveINFO_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=AF,Number=1,Type=Float,Description=\"\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tAF=0.5;DP=10\tGT\t0/0\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	// Without --keep-INFO or --recode-INFO-all, vcftools normally strips all
	// INFO. Combining --recode-INFO-all is not how vcftools' --remove-INFO is
	// most often used, but our implementation treats --remove-INFO alone as
	// "preserve everything except the listed tags" by default, matching
	// upstream behaviour.
	params := &Params{OutPrefix: prefix, Recode: true, RemoveINFO: "DP"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".recode.vcf")
	body := string(data)
	for _, line := range strings.Split(body, "\n") {
		if strings.HasPrefix(line, "1\t") {
			fields := strings.Split(line, "\t")
			info := fields[7]
			if info != "AF=0.5" {
				t.Errorf("INFO = %q, want %q", info, "AF=0.5")
			}
			return
		}
	}
	t.Fatalf("data row not found")
}
