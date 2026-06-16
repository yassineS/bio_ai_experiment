package samtools

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// This file exercises the `-@/--threads` block-parallel BGZF *input* decode
// added to the read-heavy subcommands (view, flagstat, idxstats, stats, depth).
// The contract under test is the cross-thread byte-identity invariant: parallel
// BGZF decode is deterministic, so every subcommand's output must be
// byte-for-byte identical for any thread count. Only decode throughput changes.

// bigBAM builds a BGZF-wrapped BAM stream from bigSAM(n) by running View with
// BAM output. It is the multi-block fixture the input-threading tests decode at
// several thread counts. The BAM body fills many BGZF blocks so the parallel
// reader's worker pool and ordering collector are genuinely exercised.
func bigBAM(t *testing.T, n int) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := View(bytes.NewReader(bigSAM(n)), &out, ViewOptions{OutputBAM: true, WithHeader: true}); err != nil {
		t.Fatalf("build BAM fixture: %v", err)
	}
	return out.Bytes()
}

// TestInputThreads_View_ByteIdentical asserts View's text-SAM output decoded
// from a BAM input is byte-identical for -@ {1,2,4,8}.
func TestInputThreads_View_ByteIdentical(t *testing.T) {
	bam := bigBAM(t, 20000)
	viewSAM := func(threads int) []byte {
		var out bytes.Buffer
		if _, err := View(bytes.NewReader(bam), &out, ViewOptions{WithHeader: true, Threads: threads}); err != nil {
			t.Fatalf("View -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := viewSAM(1)
	for _, threads := range []int{2, 4, 8} {
		got := viewSAM(threads)
		if !bytes.Equal(base, got) {
			t.Fatalf("view -@ %d output differs from -@ 1 (%d vs %d bytes)", threads, len(got), len(base))
		}
	}
}

// TestInputThreads_Flagstat_ByteIdentical asserts the flagstat report is
// byte-identical across thread counts when decoding a BAM input.
func TestInputThreads_Flagstat_ByteIdentical(t *testing.T) {
	bam := bigBAM(t, 20000)
	report := func(threads int) []byte {
		var out bytes.Buffer
		if err := FlagstatThreaded(bytes.NewReader(bam), &out, threads); err != nil {
			t.Fatalf("FlagstatThreaded -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := report(1)
	for _, threads := range []int{2, 4, 8} {
		got := report(threads)
		if !bytes.Equal(base, got) {
			t.Fatalf("flagstat -@ %d report differs from -@ 1:\n--- base ---\n%s\n--- got ---\n%s", threads, base, got)
		}
	}
}

// TestInputThreads_Stats_ByteIdentical asserts the stats report is
// byte-identical across thread counts. Stats reads from a raw (still BGZF)
// stream so the parallel reader engages.
func TestInputThreads_Stats_ByteIdentical(t *testing.T) {
	bam := bigBAM(t, 20000)
	report := func(threads int) []byte {
		var out bytes.Buffer
		if err := Stats(bytes.NewReader(bam), &out, StatsOptions{Threads: threads}); err != nil {
			t.Fatalf("Stats -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := report(1)
	for _, threads := range []int{2, 4, 8} {
		got := report(threads)
		if !bytes.Equal(base, got) {
			t.Fatalf("stats -@ %d report differs from -@ 1 (%d vs %d bytes)", threads, len(got), len(base))
		}
	}
}

// TestInputThreads_Depth_ByteIdentical asserts the per-position depth output is
// byte-identical across thread counts.
func TestInputThreads_Depth_ByteIdentical(t *testing.T) {
	bam := bigBAM(t, 8000)
	depth := func(threads int) []byte {
		var out bytes.Buffer
		if err := Depth([]io.Reader{bytes.NewReader(bam)}, &out, DepthOptions{Threads: threads}); err != nil {
			t.Fatalf("Depth -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := depth(1)
	for _, threads := range []int{2, 4, 8} {
		got := depth(threads)
		if !bytes.Equal(base, got) {
			t.Fatalf("depth -@ %d output differs from -@ 1 (%d vs %d bytes)", threads, len(got), len(base))
		}
	}
}

// TestInputThreads_Idxstats_ByteIdentical asserts the no-index idxstats scan is
// byte-identical across thread counts. The BAM is written to a temp file with no
// sibling .bai so the parallel fallback scan is what runs.
func TestInputThreads_Idxstats_ByteIdentical(t *testing.T) {
	bam := bigBAM(t, 20000)
	dir := t.TempDir()
	path := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(path, bam, 0o644); err != nil {
		t.Fatalf("write BAM: %v", err)
	}
	report := func(threads int) []byte {
		var out bytes.Buffer
		if err := IdxstatsFile(path, &out, threads); err != nil {
			t.Fatalf("IdxstatsFile -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := report(1)
	for _, threads := range []int{2, 4, 8} {
		got := report(threads)
		if !bytes.Equal(base, got) {
			t.Fatalf("idxstats -@ %d output differs from -@ 1:\n%s\nvs\n%s", threads, base, got)
		}
	}
}

// TestInputThreads_RecordsIdenticalToSerial decodes the BAM fixture through the
// parallel reader at several thread counts and confirms the exact record stream
// (one QNAME/FLAG/RNAME/POS key per record) matches the single-threaded decode.
func TestInputThreads_RecordsIdenticalToSerial(t *testing.T) {
	bam := bigBAM(t, 20000)
	base := recordKeys(t, decompressBGZF(t, bam)) // serial decode of the BAM body
	viewKeys := func(threads int) []string {
		var out bytes.Buffer
		if _, err := View(bytes.NewReader(bam), &out, ViewOptions{OutputBAM: true, WithHeader: true, Threads: threads}); err != nil {
			t.Fatalf("View -@ %d: %v", threads, err)
		}
		return recordKeys(t, decompressBGZF(t, out.Bytes()))
	}
	for _, threads := range []int{1, 2, 4, 8} {
		got := viewKeys(threads)
		if len(got) != len(base) {
			t.Fatalf("threads=%d record count %d != serial %d", threads, len(got), len(base))
		}
		for i := range base {
			if got[i] != base[i] {
				t.Fatalf("threads=%d record %d mismatch: %q vs %q", threads, i, got[i], base[i])
			}
		}
	}
}

// BenchmarkInputThreads_FlagstatBAM measures BAM input-decode throughput at 1 vs
// several BGZF decode workers, demonstrating the parallel read path's speedup.
// Run with:
//
//	go test -run=^$ -bench=BenchmarkInputThreads_FlagstatBAM ./tools/samtools/...
func BenchmarkInputThreads_FlagstatBAM(b *testing.B) {
	var sb bytes.Buffer
	if _, err := View(bytes.NewReader(bigSAM(200000)), &sb, ViewOptions{OutputBAM: true, WithHeader: true}); err != nil {
		b.Fatalf("build BAM fixture: %v", err)
	}
	bam := sb.Bytes()
	for _, threads := range []int{1, 2, 4, 8} {
		b.Run(fmt.Sprintf("threads=%d", threads), func(b *testing.B) {
			b.SetBytes(int64(len(bam)))
			for i := 0; i < b.N; i++ {
				var out bytes.Buffer
				if err := FlagstatThreaded(bytes.NewReader(bam), &out, threads); err != nil {
					b.Fatalf("FlagstatThreaded: %v", err)
				}
			}
		})
	}
}
