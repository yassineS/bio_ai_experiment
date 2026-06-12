package samtools

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
)

// cramRegionFixture encodes the coordinate-sorted, multi-read mpileup.1.sam
// fixture into a CRAM v3 file whose slices are capped (seqs_per_slice=50) to
// force several containers, then builds a sibling .crai with the Go
// CreateCRAI. It returns the CRAM path, the .crai path, and the reference
// FASTA path. The upstream samtools binary is the encoder, so the wiring is
// exercised against a real samtools-written CRAM; a build failure is a hard
// error (never a skip), per the project's testing rules.
func cramRegionFixture(t *testing.T) (cramPath, craiPath, refPath string) {
	t.Helper()
	bin := upstreamSamtoolsBinary(t)
	root := repoRootForTest(t)
	dat := filepath.Join(root, "reference_code", "samtools", "test", "dat")
	srcSAM := filepath.Join(dat, "mpileup.1.sam")
	refPath = filepath.Join(dat, "mpileup.ref.fa")
	for _, p := range []string{srcSAM, refPath} {
		if _, err := os.Stat(p); err != nil {
			t.Fatalf("fixture %s missing: %v", p, err)
		}
	}
	dir := t.TempDir()
	cramPath = filepath.Join(dir, "in.cram")
	cmd := exec.Command(bin, "view", "-C", "-T", refPath,
		"--output-fmt-option", "seqs_per_slice=50", "-o", cramPath, srcSAM)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream samtools view -C: %v\n%s", err, out)
	}
	craiPath = cramPath + ".crai"
	if err := cram.CreateCRAI(cramPath, craiPath); err != nil {
		t.Fatalf("CreateCRAI: %v", err)
	}
	// Confirm the fixture really has multiple containers, otherwise the
	// multi-container-spanning region case below would prove nothing.
	idx, err := cram.OpenCRAI(craiPath)
	if err != nil {
		t.Fatalf("OpenCRAI: %v", err)
	}
	distinct := map[int64]struct{}{}
	for _, e := range idx.Entries {
		distinct[e.ContainerOffset] = struct{}{}
	}
	if len(distinct) < 2 {
		t.Fatalf("expected a multi-container CRAM, got %d containers", len(distinct))
	}
	t.Logf("CRAM fixture spans %d containers (%d crai entries)", len(distinct), len(idx.Entries))
	return cramPath, craiPath, refPath
}

// upstreamViewRegion runs the upstream samtools `view <region> file.cram -T
// ref` and returns its body record lines (no header), sorted for an
// order-independent comparison.
func upstreamViewRegion(t *testing.T, bin, cramPath, refPath, region string) []string {
	t.Helper()
	cmd := exec.Command(bin, "view", "-T", refPath, cramPath, region)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream samtools view %s: %v", region, err)
	}
	return sortedBodyLines(string(out))
}

// upstreamViewRegionKeepMDNM is upstreamViewRegion but keeps the
// reference-derived MD:Z and NM:i aux tags, so a caller can assert byte-for-
// byte parity including the regenerated tags.
func upstreamViewRegionKeepMDNM(t *testing.T, bin, cramPath, refPath, region string) []string {
	t.Helper()
	cmd := exec.Command(bin, "view", "-T", refPath, cramPath, region)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream samtools view %s: %v", region, err)
	}
	return sortedBodyLinesKeepMDNM(string(out))
}

// sortedBodyLinesKeepMDNM is sortedBodyLines without the MD/NM strip: it
// splits SAM text into non-header, non-empty lines verbatim and sorts them.
// It is used by the parity test that asserts the MD:Z and NM:i tags the
// CRAM decoder now regenerates match upstream exactly.
func sortedBodyLinesKeepMDNM(s string) []string {
	var out []string
	for _, line := range strings.Split(s, "\n") {
		if line == "" || strings.HasPrefix(line, "@") {
			continue
		}
		out = append(out, line)
	}
	sort.Strings(out)
	return out
}

// sortedBodyLines splits SAM text into non-header, non-empty lines, strips
// the auto-regenerated MD/NM aux tags, and sorts the result. Sorting makes
// the comparison independent of intra-region record ordering (both
// implementations emit file order, but the test should not depend on it).
//
// The MD/NM strip isolates the region-query behaviour under test from a
// known, orthogonal decoder difference: upstream `samtools view` of a
// reference-backed CRAM recomputes the MD:Z and NM:i tags from the reference,
// whereas this repo's CRAM decoder does not yet regenerate them. Every other
// field — flag, position, CIGAR, SEQ, QUAL, and all remaining tags — is
// compared verbatim, so the region's record set and its reconstructed content
// are still pinned exactly against upstream.
func sortedBodyLines(s string) []string {
	var out []string
	for _, line := range strings.Split(s, "\n") {
		if line == "" || strings.HasPrefix(line, "@") {
			continue
		}
		out = append(out, stripMDNM(line))
	}
	sort.Strings(out)
	return out
}

// stripMDNM removes any MD:Z: and NM:i: aux fields from a tab-separated SAM
// record line.
func stripMDNM(line string) string {
	fields := strings.Split(line, "\t")
	kept := fields[:0]
	for _, f := range fields {
		if strings.HasPrefix(f, "MD:Z:") || strings.HasPrefix(f, "NM:i:") {
			continue
		}
		kept = append(kept, f)
	}
	return strings.Join(kept, "\t")
}

// TestViewCRAMIndexedRegionUpstreamParity is the LIVE parity check for the
// .crai-indexed CRAM region query wired into ViewFile. It builds (or reuses)
// the vendored upstream samtools binary, encodes a multi-container CRAM, and
// compares `samtools view <region> in.cram -T ref` (upstream) against
// ViewFile record-for-record over several regions: an empty region, a
// single-container region, a region spanning multiple containers, and the
// whole reference. Per the project rules it never skips on a build failure.
func TestViewCRAMIndexedRegionUpstreamParity(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	bin := upstreamSamtoolsBinary(t)
	cramPath, _, refPath := cramRegionFixture(t)

	regions := []string{
		"17:1-200",       // first container.
		"17:1000-2000",   // spans several containers.
		"17:4000-4200",   // tail of the reference.
		"17",             // whole reference (open-ended, every container).
		"17:99000-99999", // empty: no overlapping reads.
	}
	for _, region := range regions {
		t.Run(region, func(t *testing.T) {
			want := upstreamViewRegion(t, bin, cramPath, refPath, region)

			var buf bytes.Buffer
			n, err := ViewFile(cramPath, &buf, ViewOptions{
				Regions:    []string{region},
				Reference:  refPath,
				WithHeader: false,
			}, io.Discard)
			if err != nil {
				t.Fatalf("ViewFile(%s): %v", region, err)
			}
			got := sortedBodyLines(buf.String())
			if n != len(got) {
				t.Fatalf("region %s: ViewFile reported %d records but wrote %d lines", region, n, len(got))
			}
			if len(got) != len(want) {
				t.Fatalf("region %s: got %d records, upstream %d", region, len(got), len(want))
			}
			for i := range want {
				if got[i] != want[i] {
					t.Fatalf("region %s record %d mismatch:\n got=%q\nwant=%q", region, i, got[i], want[i])
				}
			}
		})
	}
}

// TestViewCRAMNoIndexFallback proves that when no .crai is present, a CRAM
// region query falls back to the streaming linear-scan path, emits a warning,
// and still produces exactly the same records the indexed seek path would.
// (Upstream samtools cannot itself answer a region query on an unindexed
// CRAM, so the indexed result is the oracle here.)
func TestViewCRAMNoIndexFallback(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	cramPath, craiPath, refPath := cramRegionFixture(t)
	region := "17:1000-2000"

	// Indexed result (the oracle) while the .crai still exists.
	var indexedBuf bytes.Buffer
	if _, err := ViewFile(cramPath, &indexedBuf, ViewOptions{
		Regions:   []string{region},
		Reference: refPath,
	}, io.Discard); err != nil {
		t.Fatalf("ViewFile(indexed): %v", err)
	}
	want := sortedBodyLines(indexedBuf.String())
	if len(want) == 0 {
		t.Fatalf("indexed query returned no records for %s", region)
	}

	// Remove the .crai so the streaming fallback path is taken.
	if err := os.Remove(craiPath); err != nil {
		t.Fatalf("remove crai: %v", err)
	}
	var warn bytes.Buffer
	var buf bytes.Buffer
	if _, err := ViewFile(cramPath, &buf, ViewOptions{
		Regions:   []string{region},
		Reference: refPath,
	}, &warn); err != nil {
		t.Fatalf("ViewFile(fallback): %v", err)
	}
	got := sortedBodyLines(buf.String())
	if len(got) != len(want) {
		t.Fatalf("fallback: got %d records, indexed oracle %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("fallback record %d mismatch:\n got=%q\nwant=%q", i, got[i], want[i])
		}
	}
	if !strings.Contains(warn.String(), "no CRAM index") {
		t.Errorf("expected a no-index warning, got %q", warn.String())
	}
}

// TestViewCRAMRemoteIndexedQuery proves the .crai-indexed CRAM region query
// works against a CRAM served over HTTP with Range support: ViewFile fetches
// the sibling .crai via hfile, opens the CRAM as a ranged-GET-backed seekable
// handle, and seeks to the relevant containers — never reading the whole
// object. The remote result must equal the local result.
func TestViewCRAMRemoteIndexedQuery(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	cramPath, craiPath, refPath := cramRegionFixture(t)
	cramBytes, err := os.ReadFile(cramPath)
	if err != nil {
		t.Fatalf("read cram: %v", err)
	}
	craiBytes, err := os.ReadFile(craiPath)
	if err != nil {
		t.Fatalf("read crai: %v", err)
	}

	mux := http.NewServeMux()
	serve := func(name string, body []byte) {
		mux.HandleFunc("/"+name, func(w http.ResponseWriter, r *http.Request) {
			http.ServeContent(w, r, name, time.Unix(0, 0), bytes.NewReader(body))
		})
	}
	serve("in.cram", cramBytes)
	serve("in.cram.crai", craiBytes)
	srv := httptest.NewServer(mux)
	defer srv.Close()

	region := "17:1000-2000"
	opts := ViewOptions{Regions: []string{region}, Reference: refPath}

	var localOut bytes.Buffer
	if _, err := ViewFile(cramPath, &localOut, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(local): %v", err)
	}
	var remoteOut bytes.Buffer
	if _, err := ViewFile(srv.URL+"/in.cram", &remoteOut, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(remote): %v", err)
	}
	local := sortedBodyLines(localOut.String())
	remote := sortedBodyLines(remoteOut.String())
	if len(remote) == 0 {
		t.Fatalf("remote query returned no records (expected a non-empty region)")
	}
	if len(local) != len(remote) {
		t.Fatalf("remote/local count differ: local %d, remote %d", len(local), len(remote))
	}
	for i := range local {
		if local[i] != remote[i] {
			t.Fatalf("remote/local record %d differ:\n local=%q\nremote=%q", i, local[i], remote[i])
		}
	}
}

// TestViewCRAMMDNMUpstreamParity is the LIVE parity check for the
// reference-derived MD:Z and NM:i aux tags the CRAM decoder regenerates. It
// asserts the Go ViewFile output matches `samtools view -T ref file.cram`
// byte-for-byte INCLUDING the MD and NM tags — the strip that the other CRAM
// region tests apply is deliberately absent here. The whole-reference query
// exercises every container; the sub-region exercises the indexed seek path.
// Per the project rules it never skips on a build failure.
func TestViewCRAMMDNMUpstreamParity(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	bin := upstreamSamtoolsBinary(t)
	cramPath, _, refPath := cramRegionFixture(t)

	// Sanity: upstream must actually emit MD/NM for this reference-backed
	// CRAM, otherwise the assertion would be vacuous.
	whole := upstreamViewRegionKeepMDNM(t, bin, cramPath, refPath, "17")
	sawTags := false
	for _, line := range whole {
		if strings.Contains(line, "\tMD:Z:") && strings.Contains(line, "\tNM:i:") {
			sawTags = true
			break
		}
	}
	if !sawTags {
		t.Fatalf("upstream emitted no MD/NM tags for the fixture — the parity check would be vacuous")
	}

	regions := []string{
		"17",           // whole reference: every container.
		"17:1000-2000", // a sub-region spanning several containers.
	}
	for _, region := range regions {
		t.Run(region, func(t *testing.T) {
			want := upstreamViewRegionKeepMDNM(t, bin, cramPath, refPath, region)

			var buf bytes.Buffer
			if _, err := ViewFile(cramPath, &buf, ViewOptions{
				Regions:    []string{region},
				Reference:  refPath,
				WithHeader: false,
			}, io.Discard); err != nil {
				t.Fatalf("ViewFile(%s): %v", region, err)
			}
			got := sortedBodyLinesKeepMDNM(buf.String())
			if len(got) != len(want) {
				t.Fatalf("region %s: got %d records, upstream %d", region, len(got), len(want))
			}
			for i := range want {
				if got[i] != want[i] {
					t.Fatalf("region %s record %d mismatch (MD/NM kept):\n got=%q\nwant=%q", region, i, got[i], want[i])
				}
			}
		})
	}
}
