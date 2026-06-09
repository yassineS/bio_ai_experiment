package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

// statsUSRSection extracts the data rows of the USR (user-tstv) section from
// a full stats report, dropping the `# USR` comment/header line so the golden
// comparison only covers the per-bin rows (matching the expected fixtures,
// which were captured from upstream `bcftools stats` with grep '^USR').
func statsUSRSection(stats []byte) []byte {
	var out bytes.Buffer
	for _, line := range bytes.Split(stats, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("USR:")) {
			out.Write(line)
			out.WriteByte('\n')
		}
	}
	return out.Bytes()
}

func runUserTSTVStats(t *testing.T, in []byte, spec UserTSTVSpec) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Stats(bytes.NewReader(in), &out, StatsOptions{UserTSTV: []UserTSTVSpec{spec}}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	return out.Bytes()
}

// TestParseUserTSTV covers the TAG[:min:max:n] grammar including defaults,
// an explicit value index, partial binning specs, and error cases.
func TestParseUserTSTV(t *testing.T) {
	tests := []struct {
		name    string
		in      string
		want    UserTSTVSpec
		wantErr bool
	}{
		{"tag only defaults", "DP", UserTSTVSpec{Tag: "DP", Idx: 0, Min: 0, Max: 1, NBins: 100}, false},
		{"full spec", "BQB:0:1:10", UserTSTVSpec{Tag: "BQB", Idx: 0, Min: 0, Max: 1, NBins: 10}, false},
		{"int range", "MQ:0:100:5", UserTSTVSpec{Tag: "MQ", Idx: 0, Min: 0, Max: 100, NBins: 5}, false},
		{"indexed tag", "PV4[1]:0:1:4", UserTSTVSpec{Tag: "PV4", Idx: 1, Min: 0, Max: 1, NBins: 4}, false},
		{"indexed tag no binning", "PV4[2]", UserTSTVSpec{Tag: "PV4", Idx: 2, Min: 0, Max: 1, NBins: 100}, false},
		{"partial min only", "DP:5", UserTSTVSpec{Tag: "DP", Idx: 0, Min: 5, Max: 1, NBins: 100}, false},
		{"negative min", "DP:-2:2:4", UserTSTVSpec{Tag: "DP", Idx: 0, Min: -2, Max: 2, NBins: 4}, false},
		{"empty", "", UserTSTVSpec{}, true},
		{"bad min", "DP:x:1:4", UserTSTVSpec{}, true},
		{"bad n", "DP:0:1:x", UserTSTVSpec{}, true},
		{"zero bins", "DP:0:1:0", UserTSTVSpec{}, true},
		{"negative index", "PV4[-1]", UserTSTVSpec{}, true},
		{"too many fields", "DP:0:1:4:9", UserTSTVSpec{}, true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseUserTSTV(tt.in)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("ParseUserTSTV(%q) = %+v, want error", tt.in, got)
				}
				return
			}
			if err != nil {
				t.Fatalf("ParseUserTSTV(%q) unexpected error: %v", tt.in, err)
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Fatalf("ParseUserTSTV(%q) = %+v, want %+v", tt.in, got, tt.want)
			}
		})
	}
}

// TestUserTSTVBin checks the binning rule against the upstream
// (val-min)/(max-min)*(nbins-1) formula, including the clamped endpoints.
func TestUserTSTVBin(t *testing.T) {
	spec := UserTSTVSpec{Min: 0, Max: 1, NBins: 10}
	tests := []struct {
		val  float64
		want int
	}{
		{-0.5, 0}, // <= min clamps to 0
		{0, 0},    // == min
		{0.05, 0}, // 0.05 * 9 = 0.45 -> 0
		{0.15, 1}, // 0.15 * 9 = 1.35 -> 1
		{0.95, 8}, // 0.95 * 9 = 8.55 -> 8
		{1.0, 9},  // >= max clamps to last bin
		{1.5, 9},  // > max clamps to last bin
	}
	for _, tt := range tests {
		if got := userTSTVBin(spec, tt.val); got != tt.want {
			t.Errorf("userTSTVBin(%v) = %d, want %d", tt.val, got, tt.want)
		}
	}
}

// TestParityStats_UserTSTV_Float validates the Float-typed (BQB) USR section
// byte-for-byte against the upstream-captured golden fixture.
func TestParityStats_UserTSTV_Float(t *testing.T) {
	in := readStatsFixture(t, "user_tstv.vcf")
	want := readStatsFixture(t, "user_tstv_bqb.expected.txt")
	got := statsUSRSection(runUserTSTVStats(t, in, UserTSTVSpec{Tag: "BQB", Min: 0, Max: 1, NBins: 10}))
	equalBytes(t, got, want, "stats -u BQB (float)")
}

// TestParityStats_UserTSTV_Int validates the Integer-typed (MQ) USR section
// byte-for-byte against the upstream-captured golden fixture (note the int
// bin labels print without a decimal exponent).
func TestParityStats_UserTSTV_Int(t *testing.T) {
	in := readStatsFixture(t, "user_tstv.vcf")
	want := readStatsFixture(t, "user_tstv_mq.expected.txt")
	got := statsUSRSection(runUserTSTVStats(t, in, UserTSTVSpec{Tag: "MQ", Min: 0, Max: 100, NBins: 5}))
	equalBytes(t, got, want, "stats -u MQ (int)")
}

// TestStatsUserTSTV_MissingTag confirms that records without the requested
// INFO tag are silently skipped (no panic, empty section).
func TestStatsUserTSTV_MissingTag(t *testing.T) {
	in := []byte("##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=AB,Number=1,Type=Float,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t50\t.\t.\n")
	got := statsUSRSection(runUserTSTVStats(t, in, UserTSTVSpec{Tag: "AB", Min: 0, Max: 1, NBins: 4}))
	if len(got) != 0 {
		t.Fatalf("expected empty USR section for missing tag, got:\n%s", got)
	}
}

// readStatsFixture reads a file from testdata/parity/stats, skipping the test
// when the fixture (or upstream-derived golden) is unavailable.
func readStatsFixture(t *testing.T, name string) []byte {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "stats", name))
	if err != nil {
		t.Fatalf("abs path: %v", err)
	}
	data, err := os.ReadFile(abs)
	if err != nil {
		t.Skipf("stats fixture %s unavailable: %v", name, err)
	}
	return data
}
