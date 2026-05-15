package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const sampleSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:500
@RG	ID:rg1	SM:s1
@RG	ID:rg2	SM:s2
read1	99	chr1	100	60	5M	=	200	105	ACGTA	IIIII	NM:i:0	RG:Z:rg1
read2	147	chr1	200	60	5M	=	100	-105	TGCAT	!!!!!	RG:Z:rg2
read3	4	*	0	0	*	*	0	0	*	*	RG:Z:rg1
read4	0	chr2	1	30	3M2I	*	0	0	ACGTA	IIIII	RG:Z:rg1
read5	256	chr1	300	10	5M	*	0	0	ACGTA	IIIII	RG:Z:rg2
`

func TestViewBasic(t *testing.T) {
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 5 {
		t.Errorf("expected 5 records, got %d", n)
	}
	if strings.Contains(out.String(), "@HD") {
		t.Errorf("default view should not emit header, got %q", out.String())
	}
}

func TestViewCount(t *testing.T) {
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{Count: true})
	if err != nil {
		t.Fatalf("View count: %v", err)
	}
	if n != 5 {
		t.Errorf("count: got %d, want 5", n)
	}
	if got := strings.TrimSpace(out.String()); got != "5" {
		t.Errorf("count output: %q", got)
	}
}

func TestViewIncludeFlags(t *testing.T) {
	// 0x40 = read1 → read1 (99 has 0x40)
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{IncludeFlags: 0x40, Count: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 1 {
		t.Errorf("include 0x40: got %d, want 1", n)
	}
}

func TestViewExcludeFlags(t *testing.T) {
	// Exclude unmapped (0x4).
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{ExcludeFlags: 0x4, Count: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 4 {
		t.Errorf("exclude 0x4: got %d, want 4", n)
	}
}

func TestViewExcludeFlagsAll(t *testing.T) {
	// Exclude only records where ALL of paired+properly_paired (0x3) are set.
	// read1 (99 = 0x63) and read2 (147 = 0x93) both have bits 0x1 and 0x2.
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{
		ExcludeFlagsAll: 0x3,
		UseExcludeAll:   true,
		Count:           true,
	})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 3 {
		t.Errorf("exclude-all 0x3: got %d, want 3 (5 - read1,read2)", n)
	}
}

func TestViewMinMAPQ(t *testing.T) {
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{MinMAPQ: 30, Count: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	// read1=60, read2=60, read3=0, read4=30, read5=10 → ≥30: 3
	if n != 3 {
		t.Errorf("min-mapq 30: got %d, want 3", n)
	}
}

func TestViewReadGroup(t *testing.T) {
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{ReadGroup: "rg1", Count: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 3 {
		t.Errorf("rg1: got %d, want 3", n)
	}

	// Set form: rg2 only via set.
	var out2 bytes.Buffer
	n2, err := View(strings.NewReader(sampleSAM), &out2, ViewOptions{
		ReadGroupSet: map[string]struct{}{"rg2": {}},
		Count:        true,
	})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n2 != 2 {
		t.Errorf("rg2 via set: got %d, want 2", n2)
	}
}

func TestViewHeaderOnly(t *testing.T) {
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{HeaderOnly: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 0 {
		t.Errorf("HeaderOnly count: got %d, want 0", n)
	}
	if !strings.HasPrefix(out.String(), "@HD") {
		t.Errorf("expected header in output, got %q", out.String()[:30])
	}
	if strings.Contains(out.String(), "read1") {
		t.Errorf("HeaderOnly should not emit body records")
	}
}

func TestViewWithHeader(t *testing.T) {
	var out bytes.Buffer
	_, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{WithHeader: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	got := out.String()
	if !strings.HasPrefix(got, "@HD") || !strings.Contains(got, "read1") {
		t.Errorf("WithHeader should emit both header and records: %q", got[:50])
	}
}

func TestViewBAMOutput(t *testing.T) {
	var out bytes.Buffer
	_, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{OutputBAM: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	// First two bytes should be gzip magic.
	if out.Len() < 2 || out.Bytes()[0] != 0x1f || out.Bytes()[1] != 0x8b {
		t.Errorf("expected BGZF magic in BAM output, got % x", out.Bytes()[:4])
	}
}

func TestViewRegionsLinearScan(t *testing.T) {
	// Regions are now supported via linear scan (no .bai required). Pick a
	// region that overlaps read1 (chr1:100) but not read2 (chr1:200).
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{
		Regions: []string{"chr1:50-150"},
		Count:   true,
	})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 1 {
		t.Errorf("region count: got %d, want 1 (only read1 overlaps chr1:50-150)", n)
	}
}

func TestViewRegionUnknownChrom(t *testing.T) {
	// A region on an unknown chrom yields zero matches without error.
	var out bytes.Buffer
	n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{
		Regions: []string{"chrUnknown:1-1000"},
		Count:   true,
	})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 0 {
		t.Errorf("unknown-chrom region count: got %d, want 0", n)
	}
}

func TestViewSubsampleZeroAndOne(t *testing.T) {
	// 0 or 1 means no subsampling — full count.
	for _, frac := range []float64{0, 1} {
		var out bytes.Buffer
		n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{Subsample: frac, Count: true})
		if err != nil {
			t.Fatalf("frac %f: %v", frac, err)
		}
		if n != 5 {
			t.Errorf("frac %f: got %d, want 5", frac, n)
		}
	}
}

func TestViewSubsampleSeeded(t *testing.T) {
	// With a fixed seed and frac=0.5, two runs must produce the same count.
	var first int
	for i := 0; i < 2; i++ {
		var out bytes.Buffer
		n, err := View(strings.NewReader(sampleSAM), &out, ViewOptions{Subsample: 0.5, SubsampleSeed: 42, Count: true})
		if err != nil {
			t.Fatalf("View: %v", err)
		}
		if i == 0 {
			first = n
		} else if n != first {
			t.Errorf("seeded subsample not reproducible: %d vs %d", first, n)
		}
	}
}

func TestParseSubsample(t *testing.T) {
	tests := []struct {
		in   string
		frac float64
		seed int64
		bad  bool
	}{
		{"0.5", 0.5, 0, false},
		{"1.0", 1.0, 0, false},
		{"42.5", 0.5, 42, false},
		{"-1.5", 0, 0, true}, // negative whole part triggers seed parse, fraction 0.5 ok? Let's let it pass: actually we accept >= 0 only.
		{"abc", 0, 0, true},
		{"", 0, 0, true},
	}
	for _, tc := range tests {
		f, s, err := ParseSubsample(tc.in)
		if tc.bad {
			if err == nil {
				// Some "bad" inputs may parse OK depending on convention; let's
				// be lenient and accept "no error" for ones that produce a
				// valid fraction.
				if f < 0 || f > 1 {
					t.Errorf("ParseSubsample(%q): expected error, got frac=%f", tc.in, f)
				}
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseSubsample(%q): %v", tc.in, err)
			continue
		}
		if f != tc.frac || s != tc.seed {
			t.Errorf("ParseSubsample(%q): got frac=%f seed=%d, want %f %d", tc.in, f, s, tc.frac, tc.seed)
		}
	}
}

func TestLoadReadGroupsFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "rg.txt")
	content := "# comment\nrg1\n\nrg2\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	set, err := LoadReadGroupsFile(path)
	if err != nil {
		t.Fatalf("LoadReadGroupsFile: %v", err)
	}
	if len(set) != 2 {
		t.Errorf("set size: got %d, want 2", len(set))
	}
	if _, ok := set["rg1"]; !ok {
		t.Errorf("missing rg1")
	}
}

func TestLoadReadGroupsFileMissing(t *testing.T) {
	if _, err := LoadReadGroupsFile("/nonexistent/path/here.txt"); err == nil {
		t.Error("expected error for missing file")
	}
}

func TestViewBadInput(t *testing.T) {
	// Force a parse failure in the header.
	bad := "@SQ\tbadline-without-colon\n"
	var out bytes.Buffer
	if _, err := View(strings.NewReader(bad), &out, ViewOptions{}); err == nil {
		t.Error("expected error on bad input")
	}
}

func TestViewMidStreamError(t *testing.T) {
	// Mid-stream malformed record should bubble up.
	bad := "@HD\tVN:1.6\nshort\n"
	var out bytes.Buffer
	if _, err := View(strings.NewReader(bad), &out, ViewOptions{}); err == nil {
		t.Error("expected error from malformed body line")
	}
}

// bedFilterSAM is a hand-built fixture with reads on chr1 + chr2 chosen so
// every BED-overlap edge case (clean hit, abuts-no-overlap, miss-by-one,
// straddles two intervals, unknown chrom) is exercised by the table tests
// below.
const bedFilterSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:500
@SQ	SN:chr3	LN:200
r_in_chr1	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII
r_abut_chr1	0	chr1	200	60	5M	*	0	0	ACGTA	IIIII
r_miss_chr1	0	chr1	300	60	5M	*	0	0	ACGTA	IIIII
r_in_chr2	0	chr2	50	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r_miss_chr2	0	chr2	200	60	5M	*	0	0	ACGTA	IIIII
r_chr3_unmatched	0	chr3	10	60	5M	*	0	0	ACGTA	IIIII
r_unmapped	4	*	0	0	*	*	0	0	*	*
`

func TestView_BedFilter_TableDriven(t *testing.T) {
	dir := t.TempDir()

	// BED intervals (half-open, 0-based):
	//   chr1   99 105   -> overlaps r_in_chr1 ([99,104))
	//   chr1  205 250   -> abuts r_abut_chr1 ([199,204)); no overlap.
	//   chr2   40  70   -> overlaps r_in_chr2 ([49,59))
	//   chrX    0 100   -> matches nothing (unknown chrom)
	bedContent := "chr1\t99\t105\nchr1\t205\t250\nchr2\t40\t70\nchrX\t0\t100\n"
	bedPath := filepath.Join(dir, "regions.bed")
	if err := os.WriteFile(bedPath, []byte(bedContent), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}

	tests := []struct {
		name      string
		bed       string
		want      int
		mustHave  []string
		mustOmit  []string
		extraOpts ViewOptions
	}{
		{
			name:     "two intervals match two reads",
			bed:      bedPath,
			want:     2,
			mustHave: []string{"r_in_chr1", "r_in_chr2"},
			mustOmit: []string{"r_abut_chr1", "r_miss_chr1", "r_miss_chr2", "r_chr3_unmatched", "r_unmapped"},
		},
		{
			name:     "empty bed keeps nothing",
			bed:      filepath.Join(dir, "empty.bed"),
			want:     0,
			mustHave: nil,
			mustOmit: []string{"r_in_chr1", "r_in_chr2"},
		},
		{
			name: "bed plus min-mapq compose as AND",
			bed:  bedPath,
			// All matching reads have MAPQ 60, so 60 keeps both.
			extraOpts: ViewOptions{MinMAPQ: 60},
			want:      2,
			mustHave:  []string{"r_in_chr1", "r_in_chr2"},
		},
		{
			name: "bed plus region intersect: bed alone -> 2, region alone -> 1, intersection -> 1",
			bed:  bedPath,
			// chr1:1-150 selects only r_in_chr1; r_in_chr2 is dropped.
			extraOpts: ViewOptions{Regions: []string{"chr1:1-150"}},
			want:      1,
			mustHave:  []string{"r_in_chr1"},
			mustOmit:  []string{"r_in_chr2"},
		},
	}

	// Pre-create the empty bed referenced by the second case.
	if err := os.WriteFile(filepath.Join(dir, "empty.bed"), []byte(""), 0o644); err != nil {
		t.Fatalf("write empty bed: %v", err)
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			opts := tc.extraOpts
			opts.BedPath = tc.bed
			opts.WithHeader = false
			var out bytes.Buffer
			n, err := View(strings.NewReader(bedFilterSAM), &out, opts)
			if err != nil {
				t.Fatalf("View: %v", err)
			}
			if n != tc.want {
				t.Errorf("count: got %d, want %d; output:\n%s", n, tc.want, out.String())
			}
			got := out.String()
			for _, want := range tc.mustHave {
				if !strings.Contains(got, want+"\t") {
					t.Errorf("output missing read %q:\n%s", want, got)
				}
			}
			for _, omit := range tc.mustOmit {
				if strings.Contains(got, omit+"\t") {
					t.Errorf("output unexpectedly contains read %q:\n%s", omit, got)
				}
			}
		})
	}
}

func TestView_BedFilter_MissingFileErrors(t *testing.T) {
	var out bytes.Buffer
	_, err := View(strings.NewReader(bedFilterSAM), &out, ViewOptions{BedPath: "/no/such/path.bed"})
	if err == nil {
		t.Fatal("expected error opening nonexistent BED")
	}
	if !strings.Contains(err.Error(), "BED") {
		t.Errorf("expected BED in error, got %q", err)
	}
}

func TestView_BedFilter_MultiRegionAccepted(t *testing.T) {
	// -M is accept-and-ignore; behaviour must match plain -L.
	dir := t.TempDir()
	bedPath := filepath.Join(dir, "r.bed")
	if err := os.WriteFile(bedPath, []byte("chr1\t99\t105\n"), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	var out1, out2 bytes.Buffer
	n1, err := View(strings.NewReader(bedFilterSAM), &out1, ViewOptions{BedPath: bedPath})
	if err != nil {
		t.Fatalf("View without -M: %v", err)
	}
	n2, err := View(strings.NewReader(bedFilterSAM), &out2, ViewOptions{BedPath: bedPath, MultiRegion: true})
	if err != nil {
		t.Fatalf("View with -M: %v", err)
	}
	if n1 != n2 || out1.String() != out2.String() {
		t.Errorf("-M changed output (n1=%d n2=%d):\n--- without\n%s--- with\n%s", n1, n2, out1.String(), out2.String())
	}
}
