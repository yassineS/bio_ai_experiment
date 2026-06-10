package bcftools

// Live upstream-binary parity for `bcftools view -x/--private` and
// `-X/--exclude-private`.
//
// Unlike the table-driven unit tests (view_test.go), this test runs the
// *actual upstream C binary* on the same fixture and compares its record
// selection against our Go port in-process. No committed golden/snapshot
// file is involved: the expected output is produced live by the binary the
// shared upstreamBcftools helper (upstream_test.go) locates or builds under
// reference_code/bcftools.
//
// By project rule the test must always run remotely; it therefore never
// t.Skip. If the upstream binary genuinely cannot be located or built it
// fails loudly via t.Fatalf so the gap is visible in CI.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// runUpstreamView runs the upstream bcftools view on path with the given
// extra flags and returns stdout.
func runUpstreamView(t *testing.T, bin, path string, extraFlags ...string) []byte {
	t.Helper()
	args := append([]string{"view", "--no-version"}, extraFlags...)
	args = append(args, path)
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools %v failed: %v\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestView_PrivateUpstreamParity runs the upstream C binary and our Go port
// on the SAME fixture and asserts the record SELECTION (and every column the
// private filter governs) matches in-process — no committed snapshot.
//
// Upstream additionally recomputes INFO/AC/AN after sample subsetting (a
// separate, pre-existing gap tracked by TestParityView_SampleSubset). The
// private filter does not govern the INFO column, so we blank INFO (field 7)
// on both sides via dataRecordsStripINFO before comparing; every other
// column — CHROM/POS/ID/REF/ALT/QUAL/FILTER/FORMAT and the retained samples'
// genotypes — must match byte-for-byte.
func TestView_PrivateUpstreamParity(t *testing.T) {
	bin := upstreamBcftools(t)

	fixture := parityPath(t, filepath.Join("view", "private.vcf"))
	in, readErr := os.ReadFile(fixture)
	if readErr != nil {
		t.Fatalf("read fixture %s: %v", fixture, readErr)
	}

	cases := []struct {
		name string
		flag string
		opts ViewOptions
	}{
		{"private", "-x", ViewOptions{Samples: []string{"S1", "S2"}, Private: true}},
		{"exclude-private", "-X", ViewOptions{Samples: []string{"S1", "S2"}, ExcludePrivate: true}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			upstream := runUpstreamView(t, bin, fixture, "-s", "S1,S2", tc.flag)
			got := runParityView(t, in, tc.opts)

			wantRecs := dataRecordsStripINFO(string(upstream))
			gotRecs := dataRecordsStripINFO(string(got))
			if !equalStrings(gotRecs, wantRecs) {
				t.Fatalf("record selection mismatch vs live upstream bcftools.\nflag: %s\nwant: %v\ngot:  %v",
					tc.flag, wantRecs, gotRecs)
			}
		})
	}
}
