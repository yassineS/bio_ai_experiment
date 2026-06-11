package cram

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// upstreamSamtoolsCramV21 returns the live upstream samtools binary used
// by the CRAM v2.1 record-counter parity test. It delegates to the
// shared upstreamSamtoolsCram builder (one build per process, t.Fatalf
// on failure, never t.Skip) — the build machinery is identical; only the
// caller intent differs, so the v2.1 test reads as a self-contained unit
// while still sharing the single sync.Once build.
func upstreamSamtoolsCramV21(t *testing.T) string {
	t.Helper()
	return upstreamSamtoolsCram(t)
}

// TestV21RecordCounterParity proves the CRAM v2.1 slice-header record
// counter — read as ITF-8 (32-bit) for major version 2, distinct from
// the LTF-8 (64-bit) v3+ form — decodes correctly end to end against a
// real samtools-written v2.1 file.
//
// A fixture with >= 2^28 records (the only point at which ITF-8 and LTF-8
// diverge on the wire) is impractical to build, so this test pins the
// realistic path two ways: (1) it confirms every slice header of a real
// v2.1 file parses through the ITF-8 branch with a monotonic,
// non-negative record counter matching the running record total, and
// (2) it asserts our full decode of the v2.1 file is byte-for-byte equal
// to `samtools view`. The dedicated unit test
// TestParseSliceHeaderRecordCounterWidth covers the >= 2^28 divergence
// directly. Together they close the documented v2.1 counter gap.
func TestV21RecordCounterParity(t *testing.T) {
	samtools := upstreamSamtoolsCramV21(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.sam")
	if _, err := os.Stat(srcSAM); err != nil {
		t.Fatalf("source SAM fixture missing: %v", err)
	}
	cramPath := filepath.Join(t.TempDir(), "v21.cram")
	cmd := exec.Command(samtools, "view", "-C",
		"--output-fmt-option", "version=2.1",
		"-o", cramPath, srcSAM)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -C version=2.1: %v\n%s", err, out)
	}

	// Confirm the file really is v2.1 (major version 2) so the slice
	// headers genuinely exercise the ITF-8 counter branch.
	rd, err := Open(cramPath)
	if err != nil {
		t.Fatalf("Open %s: %v", cramPath, err)
	}
	if got := rd.FileDefinition().Major; got != 2 {
		rd.Close()
		t.Fatalf("expected CRAM major version 2, got %d", got)
	}
	conts, err := rd.Containers()
	if err != nil {
		rd.Close()
		t.Fatalf("Containers: %v", err)
	}
	rd.Close()

	// The record counter must advance monotonically: each slice's counter
	// is the running total of records in all preceding slices. A
	// width-misparse would either go negative or jump non-monotonically.
	var running int64
	sawSlice := false
	for _, c := range conts {
		if c.Major != 2 {
			t.Fatalf("container %d reports major %d, want 2", c.Index, c.Major)
		}
		dc, derr := ParseDataContainer(c)
		if derr != nil {
			continue
		}
		for _, sl := range dc.Slices {
			sawSlice = true
			if sl.Header.RecordCounter < 0 {
				t.Fatalf("slice record counter is negative: %d", sl.Header.RecordCounter)
			}
			if sl.Header.RecordCounter != running {
				t.Fatalf("slice record counter = %d, want running total %d",
					sl.Header.RecordCounter, running)
			}
			running += int64(sl.Header.NumRecords)
		}
	}
	if !sawSlice {
		t.Fatalf("v2.1 fixture produced no slices to check the record counter")
	}

	// Full-decode parity: our reader vs live `samtools view`, byte-for-byte.
	want := samtoolsViewRecords(t, samtools, cramPath)
	got := ourViewRecords(t, cramPath)
	if len(got) != len(want) {
		t.Fatalf("decoded %d records, samtools decoded %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("record %d mismatch:\n got=%q\nwant=%q", i, got[i], want[i])
		}
	}
}
