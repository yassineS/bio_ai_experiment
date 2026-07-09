package bedcluster

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// Live-upstream parity tests for `bedtools cluster`: they run the real upstream
// binary from the vendored submodule and compare against this port.
//
// `bedtools cluster` (and our port) assigns cluster IDs after a per-chromosome
// std::sort by start alone. For records that share a start coordinate, the
// resulting intra-chromosome ORDER is C++-stdlib defined: our pkg/cppsort (a
// libstdc++ introsort port) reproduces the libstdc++ oracle (the CI/container
// upstream) exactly, whereas a libc++ oracle (e.g. the local arm64-macOS
// bedtools) may order equal-start records differently. That reordering does not
// change the CLUSTERING (which records overlap is coordinate-defined, not
// order-defined) — only the row order within a cluster and, with it, the exact
// line-by-line output.
//
// So these tests assert the order-INDEPENDENT invariants that hold on ANY
// conforming oracle: the SET of records per cluster is identical, and every
// input record is assigned to exactly one cluster. The exact intra-cluster row
// order is stdlib-tie-dependent and intentionally not asserted here (the
// canonical ordered layout is locked by parity_test.go against the libstdc++
// oracle). They t.Skip only when the upstream binary is genuinely absent.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
)

func repoRoot(t *testing.T) string {
	t.Helper()
	_, here, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(here)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root (go.mod) above %s", here)
		}
		dir = parent
	}
}

func upstreamBedtools(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsOnce.Do(func() {
		root := repoRoot(t)
		bin := filepath.Join(root, "reference_code", "bedtools", "bin", "bedtools")
		if _, err := os.Stat(bin); err == nil {
			upstreamBedtoolsPath = bin
		}
	})
	if upstreamBedtoolsPath == "" {
		t.Skipf("upstream bedtools binary not found at reference_code/bedtools/bin/bedtools " +
			"(run `git submodule update --init reference_code/bedtools` and build it)")
	}
	return upstreamBedtoolsPath
}

func fixtureAbs(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs path: %v", err)
	}
	return abs
}

// parseClusters maps each output line to (record-without-cluster-id ->
// clusterID) and (clusterID -> set of records). The trailing column is the
// cluster ID.
func parseClusters(t *testing.T, out []byte) (map[string]string, map[string]map[string]int) {
	t.Helper()
	recToID := make(map[string]string)
	idToRecs := make(map[string]map[string]int)
	for _, line := range strings.Split(strings.TrimRight(string(out), "\n"), "\n") {
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 2 {
			t.Fatalf("unexpected cluster output line: %q", line)
		}
		id := fields[len(fields)-1]
		rec := strings.Join(fields[:len(fields)-1], "\t")
		recToID[rec] = id
		if idToRecs[id] == nil {
			idToRecs[id] = make(map[string]int)
		}
		idToRecs[id][rec]++
	}
	return recToID, idToRecs
}

// clusterSets normalises the id->records map into a set-of-sets so two cluster
// partitions can be compared independently of the numeric IDs and intra-cluster
// row order. Each cluster becomes a sorted, newline-joined key.
func clusterSets(idToRecs map[string]map[string]int) map[string]bool {
	out := make(map[string]bool, len(idToRecs))
	for _, recs := range idToRecs {
		lines := make([]string, 0, len(recs))
		for rec, c := range recs {
			for i := 0; i < c; i++ {
				lines = append(lines, rec)
			}
		}
		sortStrings(lines)
		out[strings.Join(lines, "\n")] = true
	}
	return out
}

func sortStrings(s []string) {
	// small slices; insertion sort keeps the test dependency-free
	for i := 1; i < len(s); i++ {
		for j := i; j > 0 && s[j-1] > s[j]; j-- {
			s[j-1], s[j] = s[j], s[j-1]
		}
	}
}

func assertClusterInvariants(t *testing.T, fixture string, stranded bool) {
	t.Helper()
	bin := upstreamBedtools(t)
	in := fixtureAbs(t, fixture)

	args := []string{"cluster", "-i", in}
	if stranded {
		args = append(args, "-s")
	}
	want, err := exec.Command(bin, args...).Output()
	if err != nil {
		t.Fatalf("upstream cluster (stranded=%v): %v", stranded, err)
	}

	data, err := os.ReadFile(in)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	var got bytes.Buffer
	if _, err := Cluster(bytes.NewReader(data), &got, Options{StrandSpec: stranded}); err != nil {
		t.Fatalf("Cluster (stranded=%v): %v", stranded, err)
	}

	wantRecToID, wantIDToRecs := parseClusters(t, want)
	gotRecToID, gotIDToRecs := parseClusters(t, got.Bytes())

	// Every input record must be present in both outputs exactly once (cluster
	// never drops or duplicates a record).
	if len(wantRecToID) != len(gotRecToID) {
		t.Fatalf("record-count mismatch (stranded=%v): upstream=%d ours=%d",
			stranded, len(wantRecToID), len(gotRecToID))
	}
	for rec := range wantRecToID {
		if _, ok := gotRecToID[rec]; !ok {
			t.Errorf("record %q assigned by upstream but missing from ours (stranded=%v)", rec, stranded)
		}
	}

	// The partition into clusters (the SET of records per cluster) must be
	// identical, independent of the numeric ID labels and intra-cluster row
	// order — the order-independent invariant that holds on any oracle.
	wantSets := clusterSets(wantIDToRecs)
	gotSets := clusterSets(gotIDToRecs)
	if len(wantSets) != len(gotSets) {
		t.Fatalf("cluster-count mismatch (stranded=%v): upstream=%d ours=%d\nupstream:\n%s\nours:\n%s",
			stranded, len(wantSets), len(gotSets), want, got.Bytes())
	}
	for key := range wantSets {
		if !gotSets[key] {
			t.Errorf("cluster {%s} present upstream but not in ours (stranded=%v)", key, stranded)
		}
	}
}

// The equal-start-tie fixture stresses the per-chromosome start-sort tie order
// that distinguishes libstdc++ (our pkg/cppsort target and the CI oracle) from
// libc++ (the local arm64 binary). These assert the order-independent cluster
// invariants that must hold on either oracle.
func TestLiveParity_ClusterTie_Unstranded(t *testing.T) {
	assertClusterInvariants(t, "equal_start_ties.bed", false)
}

func TestLiveParity_ClusterTie_Stranded(t *testing.T) {
	assertClusterInvariants(t, "equal_start_ties.bed", true)
}
