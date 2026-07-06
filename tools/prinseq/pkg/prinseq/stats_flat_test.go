package prinseq

// Tests for the flat `-stats_*` reporting path (CollectFlatStats) and a
// gzip-input round-trip check for the filter path.
//
// The stats tests prefer the live upstream prinseq-lite.pl oracle: they run
// `perl reference_code/prinseq/prinseq-lite.pl -fastq <in> -stats_<group>` and
// compare its STDOUT byte-for-byte against CollectFlatStats. When perl or the
// submodule is unavailable they fall back to asserting the TSV shape (group,
// key, value ordering) so the test still guards the output contract.

import (
	"bytes"
	"compress/gzip"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// statsFixtureFastq is a small deterministic FASTQ corpus exercising every
// stats group: varying lengths (assembly/length stats), an N-containing read
// (stats_ns), and an exact duplicate pair (stats_dupl).
const statsFixtureFastq = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n" +
	"@read2\nACGTNNACGTAC\n+\nIIIIIIIIIIII\n" +
	"@read3\nGGCCGGCCGGCCGGCCGGCC\n+\nHHHHHHHHHHHHHHHHHHHH\n" +
	"@read4\nGGCCGGCCGGCCGGCCGGCC\n+\nHHHHHHHHHHHHHHHHHHHH\n" +
	"@read5\nACGT\n+\nIIII\n"

// runUpstreamStats runs the upstream oracle for the given stats flag and
// returns its STDOUT, or ("", false) when the oracle is unavailable.
func runUpstreamStats(t *testing.T, input []byte, flag string) (string, bool) {
	t.Helper()
	pl := upstreamPrinseqPath()
	if pl == "" {
		return "", false
	}
	if _, err := exec.LookPath("perl"); err != nil {
		return "", false
	}
	dir := t.TempDir()
	inPath := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(inPath, input, 0o644); err != nil {
		t.Fatalf("write upstream input: %v", err)
	}
	cmd := exec.Command("perl", pl, "-fastq", inPath, flag)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	_ = cmd.Run() // prinseq may exit non-zero; STDOUT is what we compare.
	return stdout.String(), true
}

// TestFlatStatsUpstreamParity compares CollectFlatStats to the live oracle for
// each -stats_* group (and -stats_all). stats_tag is excluded from the
// byte-exact comparison because upstream's getTagFrequency depends on Perl's
// randomised hash iteration order and is nondeterministic on tiny inputs; our
// deterministic output is the fix-on-port improvement.
func TestFlatStatsUpstreamParity(t *testing.T) {
	input := []byte(statsFixtureFastq)
	cases := []struct {
		flag   string
		groups StatsGroups
	}{
		{"-stats_info", StatsGroups{Info: true}},
		{"-stats_len", StatsGroups{Len: true}},
		{"-stats_ns", StatsGroups{Ns: true}},
		{"-stats_dupl", StatsGroups{Dupl: true}},
		{"-stats_dinuc", StatsGroups{Dinuc: true}},
		{"-stats_assembly", StatsGroups{Assembly: true}},
	}
	for _, tc := range cases {
		t.Run(tc.flag, func(t *testing.T) {
			want, ok := runUpstreamStats(t, input, tc.flag)
			if !ok {
				t.Skip("upstream prinseq-lite.pl / perl unavailable")
			}
			lines, err := CollectFlatStats(bytes.NewReader(input), true, tc.groups)
			if err != nil {
				t.Fatalf("CollectFlatStats: %v", err)
			}
			got := strings.Join(lines, "\n") + "\n"
			if got != want {
				t.Errorf("%s mismatch:\n--- got ---\n%s\n--- want ---\n%s", tc.flag, got, want)
			}
		})
	}
}

// TestFlatStatsAllUpstreamParity compares the full -stats_all output (minus the
// stats_tag lines, which are nondeterministic upstream) against the oracle.
func TestFlatStatsAllUpstreamParity(t *testing.T) {
	input := []byte(statsFixtureFastq)
	want, ok := runUpstreamStats(t, input, "-stats_all")
	if !ok {
		t.Skip("upstream prinseq-lite.pl / perl unavailable")
	}
	lines, err := CollectFlatStats(bytes.NewReader(input), true, StatsGroupsAll())
	if err != nil {
		t.Fatalf("CollectFlatStats: %v", err)
	}
	got := strings.Join(lines, "\n") + "\n"
	if dropTagLines(got) != dropTagLines(want) {
		t.Errorf("-stats_all mismatch (excluding stats_tag):\n--- got ---\n%s\n--- want ---\n%s", got, want)
	}
}

// dropTagLines removes stats_tag rows (nondeterministic upstream) from a stats
// TSV blob.
func dropTagLines(s string) string {
	var out []string
	for _, l := range strings.Split(s, "\n") {
		if strings.HasPrefix(l, "stats_tag\t") {
			continue
		}
		out = append(out, l)
	}
	return strings.Join(out, "\n")
}

// TestFlatStatsShape asserts the TSV contract independent of the oracle: every
// emitted line is "stats_<group>\t<key>\t<value>", groups appear in sorted
// order, keys within a group are sorted, and -stats_all surfaces every group.
func TestFlatStatsShape(t *testing.T) {
	input := []byte(statsFixtureFastq)
	lines, err := CollectFlatStats(bytes.NewReader(input), true, StatsGroupsAll())
	if err != nil {
		t.Fatalf("CollectFlatStats: %v", err)
	}
	if len(lines) == 0 {
		t.Fatal("no stats lines emitted")
	}
	seenGroups := map[string]bool{}
	var lastGroup string
	lastKey := map[string]string{}
	for _, l := range lines {
		f := strings.Split(l, "\t")
		if len(f) != 3 {
			t.Fatalf("line %q has %d tab fields, want 3", l, len(f))
		}
		group, key := f[0], f[1]
		if !strings.HasPrefix(group, "stats_") {
			t.Errorf("group %q does not start with stats_", group)
		}
		if lastGroup != "" && group != lastGroup && group < lastGroup {
			t.Errorf("groups not sorted: %q after %q", group, lastGroup)
		}
		if prev, ok := lastKey[group]; ok && key < prev {
			t.Errorf("keys not sorted within %s: %q after %q", group, key, prev)
		}
		lastKey[group] = key
		lastGroup = group
		seenGroups[group] = true
	}
	for _, want := range []string{
		"stats_info", "stats_len", "stats_dupl", "stats_tag",
		"stats_dinuc", "stats_ns", "stats_assembly",
	} {
		if !seenGroups[want] {
			t.Errorf("-stats_all missing group %q", want)
		}
	}
}

// TestFlatStatsIndividualGroups checks that -stats_len and -stats_info emit only
// their own group's lines with the expected keys.
func TestFlatStatsIndividualGroups(t *testing.T) {
	input := []byte(statsFixtureFastq)

	infoLines, err := CollectFlatStats(bytes.NewReader(input), true, StatsGroups{Info: true})
	if err != nil {
		t.Fatalf("CollectFlatStats(info): %v", err)
	}
	wantInfo := []string{"stats_info\tbases\t", "stats_info\treads\t"}
	if len(infoLines) != len(wantInfo) {
		t.Fatalf("stats_info: got %d lines %v, want %d", len(infoLines), infoLines, len(wantInfo))
	}
	for i, prefix := range wantInfo {
		if !strings.HasPrefix(infoLines[i], prefix) {
			t.Errorf("stats_info line %d = %q, want prefix %q", i, infoLines[i], prefix)
		}
	}

	lenLines, err := CollectFlatStats(bytes.NewReader(input), true, StatsGroups{Len: true})
	if err != nil {
		t.Fatalf("CollectFlatStats(len): %v", err)
	}
	wantLenKeys := []string{"max", "mean", "median", "min", "mode", "modeval", "range", "stddev"}
	if len(lenLines) != len(wantLenKeys) {
		t.Fatalf("stats_len: got %d lines %v, want %d", len(lenLines), lenLines, len(wantLenKeys))
	}
	for i, key := range wantLenKeys {
		want := "stats_len\t" + key + "\t"
		if !strings.HasPrefix(lenLines[i], want) {
			t.Errorf("stats_len line %d = %q, want prefix %q", i, lenLines[i], want)
		}
	}
}

// TestFilterGzipInputRoundTrips verifies that Filter over a gzipped FASTQ (read
// through iohelper the same way the CLI does) yields the same records as the
// plain FASTQ. Filter itself takes an io.Reader, so this drives the gzip decode
// explicitly and asserts the two outputs match.
func TestFilterGzipInputRoundTrips(t *testing.T) {
	plain := []byte(statsFixtureFastq)

	// Compress the fixture.
	var gzbuf bytes.Buffer
	gw := gzip.NewWriter(&gzbuf)
	if _, err := gw.Write(plain); err != nil {
		t.Fatal(err)
	}
	if err := gw.Close(); err != nil {
		t.Fatal(err)
	}

	opts := FilterOptions{MinLen: 5, QualType: "sanger"}

	var plainOut bytes.Buffer
	if err := Filter(bytes.NewReader(plain), &plainOut, true, opts); err != nil {
		t.Fatalf("Filter(plain): %v", err)
	}

	gr, err := gzip.NewReader(bytes.NewReader(gzbuf.Bytes()))
	if err != nil {
		t.Fatal(err)
	}
	var gzOut bytes.Buffer
	if err := Filter(gr, &gzOut, true, opts); err != nil {
		t.Fatalf("Filter(gzip): %v", err)
	}
	gr.Close()

	if plainOut.String() != gzOut.String() {
		t.Errorf("gzip and plain inputs produced different filter output:\n--- plain ---\n%s\n--- gzip ---\n%s",
			plainOut.String(), gzOut.String())
	}
	// Sanity: the -min_len 5 filter drops the 4bp read5, keeping 4 records.
	if got := strings.Count(plainOut.String(), "@read"); got != 4 {
		t.Errorf("expected 4 kept records after -min_len 5, got %d", got)
	}
}

// TestFlatStatsDupl_RevcompStableSort is a regression guard for the reverse-
// complement dedup phase. Upstream prinseq-lite.pl sorts stably; our port must
// use sort.SliceStable in checkForDupl's revcomp expansion (graphdata.go),
// otherwise tied equal-key entries are ordered differently from Perl and the
// exactrevcomp duplicate count diverges (real chr20 R1 gave 5656 vs upstream's
// 3780 before the fix). This input is 6 distinct 20-mers each paired with its
// reverse complement as a separate read, so every pair is an exact-revcomp
// duplicate: exactrevcomp and total must be 6, deterministically.
func TestFlatStatsDupl_RevcompStableSort(t *testing.T) {
	const fq = `@r0
ACGTACGTACGTAAACCCGT
+
IIIIIIIIIIIIIIIIIIII
@r1
ACGGGTTTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@r2
TTGGCCAATTGGCCAATTGG
+
IIIIIIIIIIIIIIIIIIII
@r3
CCAATTGGCCAATTGGCCAA
+
IIIIIIIIIIIIIIIIIIII
@r4
GATCGATCGATCTTAACCGG
+
IIIIIIIIIIIIIIIIIIII
@r5
CCGGTTAAGATCGATCGATC
+
IIIIIIIIIIIIIIIIIIII
@r6
ACACACACGTGTGTGTAAAC
+
IIIIIIIIIIIIIIIIIIII
@r7
GTTTACACACACGTGTGTGT
+
IIIIIIIIIIIIIIIIIIII
@r8
GGGGTTTTCCCCAAAAGTAC
+
IIIIIIIIIIIIIIIIIIII
@r9
GTACTTTTGGGGAAAACCCC
+
IIIIIIIIIIIIIIIIIIII
@r10
TACGTACGTACGACGTTTGA
+
IIIIIIIIIIIIIIIIIIII
@r11
TCAAACGTCGTACGTACGTA
+
IIIIIIIIIIIIIIIIIIII
`
	lines, err := CollectFlatStats(strings.NewReader(fq), true, StatsGroups{Dupl: true})
	if err != nil {
		t.Fatalf("CollectFlatStats: %v", err)
	}
	got := strings.Join(lines, "\n")
	for _, want := range []string{"stats_dupl\texactrevcomp\t6", "stats_dupl\ttotal\t6"} {
		if !strings.Contains(got, want) {
			t.Errorf("missing %q in stats_dupl output:\n%s", want, got)
		}
	}
	// Determinism: the unstable-sort bug could vary run-to-run; assert stable.
	for i := 0; i < 5; i++ {
		l2, _ := CollectFlatStats(strings.NewReader(fq), true, StatsGroups{Dupl: true})
		if strings.Join(l2, "\n") != got {
			t.Fatalf("stats_dupl not deterministic across runs")
		}
	}
}
