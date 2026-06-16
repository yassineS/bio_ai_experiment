package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// statsRegionFixtureSAM covers two references with reads at separated
// positions so a positional region restricts the count in a way the live
// parity check can pin: a1/a2/a3 on chr1 (10, 50, 90), b1/b2 on chr2 (10, 60).
const statsRegionFixtureSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:100
@SQ	SN:chr2	LN:100
a1	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
a2	0	chr1	50	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
a3	0	chr1	90	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
b1	0	chr2	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
b2	0	chr2	60	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
`

// buildSortedIndexedBAM writes samText to a temp dir and uses the upstream
// samtools binary to produce a coordinate-sorted, indexed BAM, returning its
// path. The upstream binary is the live oracle harness's binary (built or
// symlinked); a failure is fatal, never skipped.
func buildSortedIndexedBAM(t *testing.T, bin, samText string) string {
	t.Helper()
	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(samText), 0o600); err != nil {
		t.Fatalf("write sam: %v", err)
	}
	bamPath := filepath.Join(dir, "in.sorted.bam")
	sortCmd := exec.Command(bin, "sort", "-o", bamPath, samPath)
	if out, err := sortCmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools sort: %v\n%s", err, out)
	}
	idxCmd := exec.Command(bin, "index", bamPath)
	if out, err := idxCmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools index: %v\n%s", err, out)
	}
	return bamPath
}

// stripStatsComments drops the leading "# ..." header lines (the version
// banner and the verbatim command line, which differ by construction between
// upstream and the port) so the comparison covers every data section.
func stripStatsComments(blob string) string {
	var b strings.Builder
	for _, line := range strings.Split(blob, "\n") {
		if strings.HasPrefix(line, "#") {
			continue
		}
		b.WriteString(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// TestStatsPositionalRegionsUpstreamParity is the LIVE parity check for
// command-line positional region arguments. It builds a sorted+indexed BAM,
// runs upstream `samtools stats in.bam <region...>` and the Go port over the
// SAME input, and asserts every data section (all sections except the
// command-line/version header comments) is byte-identical. The upstream binary
// is built on demand and a build failure is fatal — never skipped.
func TestStatsPositionalRegionsUpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)
	bamPath := buildSortedIndexedBAM(t, bin, statsRegionFixtureSAM)

	cases := []struct {
		name    string
		regions []string
	}{
		{"single_closed", []string{"chr1:1-60"}},
		{"whole_chrom", []string{"chr1"}},
		{"two_regions", []string{"chr1:1-60", "chr2"}},
		{"partial_overlap_clip", []string{"chr1:15-55"}},
		{"open_ended", []string{"chr2:30-"}},
		{"empty_region", []string{"chr1:1-9"}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Upstream.
			upArgs := append([]string{"stats", bamPath}, tc.regions...)
			upCmd := exec.Command(bin, upArgs...)
			var upOut, upErr bytes.Buffer
			upCmd.Stdout = &upOut
			upCmd.Stderr = &upErr
			if err := upCmd.Run(); err != nil {
				t.Fatalf("upstream stats %v: %v\n%s", tc.regions, err, upErr.String())
			}

			// Go port: stream the same BAM through Stats with the same regions.
			in, err := os.Open(bamPath)
			if err != nil {
				t.Fatalf("open bam: %v", err)
			}
			defer in.Close()
			var goOut bytes.Buffer
			if err := Stats(in, &goOut, StatsOptions{Regions: tc.regions}); err != nil {
				t.Fatalf("Stats: %v", err)
			}

			want := stripStatsComments(upOut.String())
			got := stripStatsComments(goOut.String())
			if got != want {
				t.Fatalf("positional-region stats differ for %v\n--- upstream ---\n%s\n--- go ---\n%s",
					tc.regions, want, got)
			}
		})
	}
}

// TestStatsPositionalRegionsWithCovThreshold pins parity for `-g` combined
// with positional regions (upstream accepts -g without -t when positional
// regions provide the target set via replicate_regions).
func TestStatsPositionalRegionsWithCovThreshold(t *testing.T) {
	bin := upstreamSamtools(t)
	bamPath := buildSortedIndexedBAM(t, bin, statsRegionFixtureSAM)

	upCmd := exec.Command(bin, "stats", "-g", "0", bamPath, "chr1:1-60")
	var upOut, upErr bytes.Buffer
	upCmd.Stdout = &upOut
	upCmd.Stderr = &upErr
	if err := upCmd.Run(); err != nil {
		t.Fatalf("upstream stats -g: %v\n%s", err, upErr.String())
	}

	in, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open bam: %v", err)
	}
	defer in.Close()
	var goOut bytes.Buffer
	if err := Stats(in, &goOut, StatsOptions{Regions: []string{"chr1:1-60"}, CovThreshold: 0}); err != nil {
		t.Fatalf("Stats: %v", err)
	}

	want := stripStatsComments(upOut.String())
	got := stripStatsComments(goOut.String())
	if got != want {
		t.Fatalf("positional-region -g stats differ\n--- upstream ---\n%s\n--- go ---\n%s", want, got)
	}
}
