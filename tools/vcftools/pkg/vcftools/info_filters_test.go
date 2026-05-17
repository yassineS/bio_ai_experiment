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

// TestFilterRecodeInfo covers --recode-INFO / --remove-INFO restriction
// of the recoded INFO column. The function is now exclusively driven by
// `Params.RecodeINFO` (and the historical recode-column `RemoveINFO`); the
// site-filter `--keep-INFO` flow is covered by TestPassKeepINFOSite below.
func TestFilterRecodeInfo(t *testing.T) {
	in := map[string]string{"AF": "0.5", "DP": "10", "AA": "T"}

	// Empty sets -> unchanged.
	out := filterRecodeInfo(in, nil, nil)
	if len(out) != 3 {
		t.Errorf("empty sets: got %d, want 3", len(out))
	}

	// recode-INFO restricts to listed tags.
	keep := parseInfoTagList("AF,DP")
	out = filterRecodeInfo(in, keep, nil)
	if len(out) != 2 || out["AF"] != "0.5" || out["DP"] != "10" {
		t.Errorf("recode-INFO AF,DP: got %v", out)
	}

	// remove-INFO strips listed tags.
	rem := parseInfoTagList("DP")
	out = filterRecodeInfo(in, nil, rem)
	if len(out) != 2 || out["AF"] != "0.5" || out["AA"] != "T" {
		t.Errorf("remove-INFO DP: got %v", out)
	}

	// Both compose: recode-INFO AF,DP then strip DP -> only AF survives.
	out = filterRecodeInfo(in, keep, rem)
	if len(out) != 1 || out["AF"] != "0.5" {
		t.Errorf("recode+remove: got %v", out)
	}
}

// TestPassKeepINFOSite covers the new upstream-aligned `--keep-INFO TAG`
// site filter (entry_filters.cpp:1033-1063). A site passes when any of
// the named INFO Flag-type tags is present; multiple tags compose via OR.
func TestPassKeepINFOSite(t *testing.T) {
	flags := parseInfoTagList("FLAG_A,FLAG_B")
	cases := []struct {
		name string
		info map[string]string
		want bool
	}{
		{"flagA present (bare)", map[string]string{"FLAG_A": ""}, true},
		{"flagA present (=1)", map[string]string{"FLAG_A": "1"}, true},
		{"flagB present", map[string]string{"FLAG_B": "", "DP": "10"}, true},
		{"both present", map[string]string{"FLAG_A": "", "FLAG_B": ""}, true},
		{"neither present", map[string]string{"DP": "10"}, false},
		{"flag explicitly zero", map[string]string{"FLAG_A": "0"}, false},
		{"empty info", map[string]string{}, false},
	}
	for _, tc := range cases {
		v := &vcf.Variant{Info: tc.info}
		got := passKeepINFOSite(v, flags)
		if got != tc.want {
			t.Errorf("%s: got %v want %v", tc.name, got, tc.want)
		}
	}

	// Empty flag set is a no-op (always pass).
	if !passKeepINFOSite(&vcf.Variant{}, nil) {
		t.Error("empty flag set should always pass")
	}
}

// TestLookupInfoMeta confirms the header scanner extracts the Type field
// from `##INFO=<...>` declarations, including those with quoted
// Description strings that embed commas.
func TestLookupInfoMeta(t *testing.T) {
	h := &vcf.Header{MetaInfo: []string{
		`##INFO=<ID=FLAG_A,Number=0,Type=Flag,Description="A flag, with comma">`,
		`##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">`,
	}}
	if meta, ok := lookupInfoMeta(h, "FLAG_A"); !ok || meta.Type != "Flag" {
		t.Errorf("FLAG_A: got (%+v, %v)", meta, ok)
	}
	if meta, ok := lookupInfoMeta(h, "DP"); !ok || meta.Type != "Integer" {
		t.Errorf("DP: got (%+v, %v)", meta, ok)
	}
	if _, ok := lookupInfoMeta(h, "MISSING"); ok {
		t.Error("MISSING: expected not-found")
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

// TestRun_RecodeINFO_Integration restricts INFO in recoded output via
// the upstream-canonical `--recode-INFO TAG` flag (parameters.cpp:319,
// recode_INFO_to_keep). This test previously used `KeepINFO` because
// the port had misaligned `--keep-INFO` with the recode-column selector
// semantic; wave-17 separates the two flags.
func TestRun_RecodeINFO_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=AF,Number=1,Type=Float,Description=\"\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tAF=0.5;DP=10;AA=T\tGT\t0/0\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Recode: true, RecodeINFO: "AF"}
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

// TestRun_KeepINFO_SiteFilter_Integration covers the wave-17 fix: the
// port's `--keep-INFO TAG` flag now matches upstream's SITE FILTER
// semantic (entry_filters.cpp:1033-1063). Only sites where the named
// INFO Flag is present should survive.
func TestRun_KeepINFO_SiteFilter_Integration(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=FLAG_A,Number=0,Type=Flag,Description=\"\">\n" +
		"##INFO=<ID=FLAG_B,Number=0,Type=Flag,Description=\"\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tFLAG_A;DP=10\tGT\t0/0\n" +
		"1\t200\t.\tA\tC\t.\tPASS\tFLAG_B;DP=20\tGT\t0/1\n" +
		"1\t300\t.\tA\tT\t.\tPASS\tFLAG_A;FLAG_B;DP=30\tGT\t1/1\n" +
		"1\t400\t.\tA\tG\t.\tPASS\tDP=40\tGT\t0/1\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	// --keep-INFO FLAG_A --recode --recode-INFO-all: only sites with
	// FLAG_A present should pass (100, 300).
	params := &Params{OutPrefix: prefix, Recode: true, RecodeInfoAll: true, KeepINFO: "FLAG_A"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".recode.vcf")
	body := string(data)
	for _, pos := range []string{"\t100\t", "\t300\t"} {
		if !strings.Contains(body, pos) {
			t.Errorf("expected site %q kept; got:\n%s", pos, body)
		}
	}
	for _, pos := range []string{"\t200\t", "\t400\t"} {
		if strings.Contains(body, pos) {
			t.Errorf("expected site %q dropped; got:\n%s", pos, body)
		}
	}
}

// TestRun_KeepINFO_SiteFilter_OR confirms multiple --keep-INFO tags
// compose via OR (entry_filters.cpp:1049-1062 loops over flags_to_keep
// and sets `keep=true` on the first present tag).
func TestRun_KeepINFO_SiteFilter_OR(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=FLAG_A,Number=0,Type=Flag,Description=\"\">\n" +
		"##INFO=<ID=FLAG_B,Number=0,Type=Flag,Description=\"\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tFLAG_A;DP=10\tGT\t0/0\n" +
		"1\t200\t.\tA\tC\t.\tPASS\tFLAG_B;DP=20\tGT\t0/1\n" +
		"1\t400\t.\tA\tG\t.\tPASS\tDP=40\tGT\t0/1\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Recode: true, RecodeInfoAll: true, KeepINFO: "FLAG_A,FLAG_B"}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".recode.vcf")
	body := string(data)
	for _, pos := range []string{"\t100\t", "\t200\t"} {
		if !strings.Contains(body, pos) {
			t.Errorf("expected site %q kept; got:\n%s", pos, body)
		}
	}
	if strings.Contains(body, "\t400\t") {
		t.Errorf("expected site 400 dropped; got:\n%s", body)
	}
}

// TestRun_KeepINFO_SiteFilter_NonFlagType confirms the port errors out
// when --keep-INFO names an INFO key that is not declared as Type=Flag
// in the header. Mirrors upstream entry_filters.cpp:1053-1054.
func TestRun_KeepINFO_SiteFilter_NonFlagType(t *testing.T) {
	vcfText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"1\t100\t.\tA\tG\t.\tPASS\tDP=10\tGT\t0/0\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Recode: true, KeepINFO: "DP"}
	err := Run(strings.NewReader(vcfText), params)
	if err == nil {
		t.Fatal("expected --keep-INFO on non-Flag type to error; got nil")
	}
	if !strings.Contains(err.Error(), "non flag type") {
		t.Errorf("error %q: want substring %q", err.Error(), "non flag type")
	}
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
