package bcftools

import (
	"bytes"
	"strings"
	"testing"
)

const filterHeader = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`

// TestVCFFilterIncludeAndExcludeSoftTag verifies that records failing
// the -i / -e expression have the FILTER column set to the soft-filter
// name AND that passing records are emitted unchanged.
func TestVCFFilterIncludeAndExcludeSoftTag(t *testing.T) {
	input := filterHeader +
		"chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/1\t1/1\n" +
		"chr1\t200\t.\tC\tG\t.\tPASS\tDP=5\tGT\t0/0\t0/1\n" +
		"chr1\t300\t.\tG\tA\t.\tPASS\tDP=20\tGT\t1/1\t0/0\n"

	cases := []struct {
		name      string
		opts      VCFFilterOptions
		wantLines []string
	}{
		{
			name: "include keeps DP>=10 and soft-filters the rest",
			opts: VCFFilterOptions{
				IncludeExpr: "INFO/DP >= 10",
				SoftFilter:  "LowDP",
			},
			wantLines: []string{
				"chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/1\t1/1",
				"chr1\t200\t.\tC\tG\t.\tLowDP\tDP=5\tGT\t0/0\t0/1",
				"chr1\t300\t.\tG\tA\t.\tPASS\tDP=20\tGT\t1/1\t0/0",
			},
		},
		{
			name: "exclude soft-filters DP<10",
			opts: VCFFilterOptions{
				ExcludeExpr: "INFO/DP < 10",
				SoftFilter:  "LowDP",
			},
			wantLines: []string{
				"chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/1\t1/1",
				"chr1\t200\t.\tC\tG\t.\tLowDP\tDP=5\tGT\t0/0\t0/1",
				"chr1\t300\t.\tG\tA\t.\tPASS\tDP=20\tGT\t1/1\t0/0",
			},
		},
		{
			name: "mode x resets FILTER to PASS on passing sites",
			opts: VCFFilterOptions{
				ExcludeExpr: "INFO/DP < 10",
				SoftFilter:  "LowDP",
				Mode:        FilterModeReset,
			},
			wantLines: []string{
				"chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/1\t1/1",
				"chr1\t200\t.\tC\tG\t.\tLowDP\tDP=5\tGT\t0/0\t0/1",
				"chr1\t300\t.\tG\tA\t.\tPASS\tDP=20\tGT\t1/1\t0/0",
			},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			opts := tc.opts
			opts.NoVersion = true
			n, err := VCFFilter(strings.NewReader(input), &out, opts)
			if err != nil {
				t.Fatalf("VCFFilter: %v", err)
			}
			if n != 3 {
				t.Fatalf("wrote %d records, want 3", n)
			}
			gotLines := dataLines(out.String())
			if len(gotLines) != len(tc.wantLines) {
				t.Fatalf("line count: got %d, want %d\nOUT:\n%s", len(gotLines), len(tc.wantLines), out.String())
			}
			for i, w := range tc.wantLines {
				if gotLines[i] != w {
					t.Errorf("line %d:\n got  %q\n want %q", i, gotLines[i], w)
				}
			}
		})
	}
}

// TestVCFFilterModeAddPreservesExistingFilters verifies that `-m +`
// appends rather than replaces.
func TestVCFFilterModeAddPreservesExistingFilters(t *testing.T) {
	input := filterHeader +
		"chr1\t100\t.\tA\tT\t.\tOldFail\tDP=5\tGT\t0/0\t0/0\n"
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(input), &out, VCFFilterOptions{
		ExcludeExpr: "INFO/DP < 10",
		SoftFilter:  "LowDP",
		Mode:        FilterModeAdd,
		NoVersion:   true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	lines := dataLines(out.String())
	if len(lines) != 1 {
		t.Fatalf("want 1 record, got %d", len(lines))
	}
	// Expect "OldFail;LowDP" — preserved + appended.
	if !strings.Contains(lines[0], "OldFail;LowDP") {
		t.Errorf("expected merged FILTER=OldFail;LowDP, got line %q", lines[0])
	}
}

// TestVCFFilterSetGTsMissing verifies that -S . rewrites GTs of failing
// records and preserves the phase separator.
func TestVCFFilterSetGTsMissing(t *testing.T) {
	input := filterHeader +
		"chr1\t100\t.\tA\tT\t.\tPASS\tDP=5\tGT\t0|1\t1/1\n"
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(input), &out, VCFFilterOptions{
		ExcludeExpr: "INFO/DP < 10",
		SoftFilter:  "LowDP",
		SetGTs:      SetGTsMissing,
		NoVersion:   true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	lines := dataLines(out.String())
	if len(lines) != 1 {
		t.Fatalf("want 1 record, got %d", len(lines))
	}
	// Sample 1 phased -> "./." with "|"; sample 2 unphased -> "./." with "/".
	if !strings.Contains(lines[0], ".|.") {
		t.Errorf("expected sample1 GT=.|., got %q", lines[0])
	}
	if !strings.Contains(lines[0], "./.") {
		t.Errorf("expected sample2 GT=./., got %q", lines[0])
	}
}

// TestVCFFilterSetGTsRef verifies that -S 0 emits "0/0" for failing
// records.
func TestVCFFilterSetGTsRef(t *testing.T) {
	input := filterHeader +
		"chr1\t100\t.\tA\tT\t.\tPASS\tDP=5\tGT\t0|1\t1/1\n"
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(input), &out, VCFFilterOptions{
		ExcludeExpr: "INFO/DP < 10",
		SoftFilter:  "LowDP",
		SetGTs:      SetGTsRef,
		NoVersion:   true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	lines := dataLines(out.String())
	if !strings.Contains(lines[0], "0|0") || !strings.Contains(lines[0], "0/0") {
		t.Errorf("expected reset GTs 0|0 and 0/0, got %q", lines[0])
	}
}

// TestVCFFilterSnpGap tags SNPs within INT bp of an indel.
func TestVCFFilterSnpGap(t *testing.T) {
	input := filterHeader +
		"chr1\t100\t.\tA\tT\t.\tPASS\tDP=20\tGT\t0/1\t0/1\n" +
		"chr1\t105\t.\tCT\tC\t.\tPASS\tDP=20\tGT\t0/1\t0/1\n" +
		"chr1\t200\t.\tG\tA\t.\tPASS\tDP=20\tGT\t0/1\t0/1\n"
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(input), &out, VCFFilterOptions{
		SnpGap:     10,
		SoftFilter: "SnpGap",
		NoVersion:  true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	lines := dataLines(out.String())
	if len(lines) != 3 {
		t.Fatalf("want 3 records, got %d", len(lines))
	}
	// SNP @100 is within 10 bp of indel @105 -> SnpGap.
	if !strings.Contains(lines[0], "\tSnpGap\t") {
		t.Errorf("SNP @100: expected SnpGap tag, got %q", lines[0])
	}
	// Indel itself is unaffected by SnpGap (only SNPs are tagged).
	if strings.Contains(lines[1], "\tSnpGap\t") {
		t.Errorf("indel @105: expected no SnpGap tag, got %q", lines[1])
	}
	// SNP @200 is far away -> PASS.
	if !strings.Contains(lines[2], "\tPASS\t") {
		t.Errorf("SNP @200: expected PASS, got %q", lines[2])
	}
}

// TestVCFFilterAutoUniqueName verifies the upstream "+" sentinel auto-
// picks a name like "Filter1".
func TestVCFFilterAutoUniqueName(t *testing.T) {
	input := filterHeader +
		"chr1\t100\t.\tA\tT\t.\tPASS\tDP=5\tGT\t0/0\t0/0\n"
	var out bytes.Buffer
	_, err := VCFFilter(strings.NewReader(input), &out, VCFFilterOptions{
		ExcludeExpr: "INFO/DP < 10",
		SoftFilter:  "+",
		NoVersion:   true,
	})
	if err != nil {
		t.Fatalf("VCFFilter: %v", err)
	}
	if !strings.Contains(out.String(), "##FILTER=<ID=Filter1") {
		t.Errorf("expected auto-named ##FILTER=<ID=Filter1 header line, got:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "\tFilter1\t") {
		t.Errorf("expected record tagged Filter1, got:\n%s", out.String())
	}
}

func TestParseFilterMode(t *testing.T) {
	cases := []struct {
		in   string
		want FilterMode
		err  bool
	}{
		{"", FilterModeReplace, false},
		{"+", FilterModeAdd, false},
		{"x", FilterModeReset, false},
		{"+x", FilterModeAdd | FilterModeReset, false},
		{"q", 0, true},
	}
	for _, tc := range cases {
		got, err := ParseFilterMode(tc.in)
		if tc.err && err == nil {
			t.Errorf("ParseFilterMode(%q) expected error", tc.in)
		}
		if !tc.err && got != tc.want {
			t.Errorf("ParseFilterMode(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}
}

func TestParseSetGTsMode(t *testing.T) {
	cases := []struct {
		in   string
		want SetGTsMode
		err  bool
	}{
		{"", SetGTsOff, false},
		{".", SetGTsMissing, false},
		{"0", SetGTsRef, false},
		{"1", 0, true},
	}
	for _, tc := range cases {
		got, err := ParseSetGTsMode(tc.in)
		if tc.err && err == nil {
			t.Errorf("ParseSetGTsMode(%q) expected error", tc.in)
		}
		if !tc.err && got != tc.want {
			t.Errorf("ParseSetGTsMode(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}
}

// dataLines strips the header and blank lines, returning the data
// records as a slice (one per row).
func dataLines(s string) []string {
	var out []string
	for _, l := range strings.Split(s, "\n") {
		if l == "" || strings.HasPrefix(l, "#") {
			continue
		}
		out = append(out, l)
	}
	return out
}
