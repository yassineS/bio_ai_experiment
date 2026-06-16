package fastp

// Live-parity tests for fastp's *multi-threaded* output splitting
// (-s/-S/-d combined with -w/--thread).
//
// Upstream fastp distributes the input across worker threads in fixed packs of
// PACK_SIZE (256) reads, assigning pack i to thread i % thread, and each thread
// owns a disjoint set of split files (start index = threadId, stride = thread).
// The result is fully deterministic for a fixed thread count, so we can compare
// every split file byte-for-byte against the upstream binary -- not just for
// -w 1 (which the parity_tail_test.go cases already cover) but for -w 2, 3 and
// 4. These tests close the "multi-thread --split file-boundary distribution"
// parity gap recorded in PROJECT_STATUS / docs/PARITY_ROADMAP.
//
// The upstream binary is built once via the shared upstreamFastp sync.Once
// helper in parity_tail_test.go; a missing/unbuildable binary t.Fatalf's with
// the exact init/build hint (env-guard policy). Any per-file byte mismatch is a
// hard failure.
//
// Inputs are generated at runtime (large enough to span many packs) rather than
// checked in, keeping testdata small. Generation is deterministic so the
// comparison is reproducible.

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// genSplitFASTQ writes n deterministic 40bp all-passing records to path. Read
// IDs are sequential ("@readK") so a per-file diff pinpoints any boundary
// mismatch. The bases cycle so the records are not all identical.
func genSplitFASTQ(t *testing.T, path string, n int) {
	t.Helper()
	var b strings.Builder
	bases := []byte("ACGT")
	for i := 0; i < n; i++ {
		b.WriteString(fmt.Sprintf("@read%d\n", i))
		seq := make([]byte, 40)
		for j := range seq {
			seq[j] = bases[(i+j)%4]
		}
		b.Write(seq)
		b.WriteString("\n+\n")
		b.WriteString(strings.Repeat("I", 40))
		b.WriteByte('\n')
	}
	if err := os.WriteFile(path, []byte(b.String()), 0o644); err != nil {
		t.Fatalf("write %s: %v", path, err)
	}
}

// TestParity_Fastp_SplitByNumber_MultiThread compares the upstream and Go
// per-file split output for -s N across several worker-thread counts. With all
// reads passing, byFileNumber rollover is keyed on input pack counts, so this
// directly validates the pack-to-thread round-robin distribution and the
// thread-strided file assignment.
func TestParity_Fastp_SplitByNumber_MultiThread(t *testing.T) {
	bin := ensureUpstream(t)

	in := filepath.Join(t.TempDir(), "in.fq")
	genSplitFASTQ(t, in, 6000) // 24 packs of 256 (last pack 96)

	for _, w := range []int{2, 3, 4} {
		w := w
		t.Run(fmt.Sprintf("w%d", w), func(t *testing.T) {
			upDir := t.TempDir()
			args := append([]string{
				"-i", in, "-o", filepath.Join(upDir, "out.fq"),
				"-s", "4", "-w", fmt.Sprintf("%d", w),
				"-j", filepath.Join(upDir, "up.json"), "-h", filepath.Join(upDir, "up.html"),
			}, commonDisableFlags...)
			runUpstream(t, bin, args)

			goDir := t.TempDir()
			opts := permissiveOpts()
			opts.SplitNumber = 4
			opts.Threads = w
			runGoFastpSESplit(t, in, filepath.Join(goDir, "out.fq"), opts)

			compareSplitDirs(t, upDir, goDir, "out.fq")
		})
	}
}

// TestParity_Fastp_SplitByLines_MultiThread compares the upstream and Go
// per-file split output for -S L across several worker-thread counts.
// byFileLines rollover is keyed on the per-thread *passed* read count (here
// equal to input, since all reads pass) and is uncapped, so this validates the
// no-cap multi-thread rollover and the trailing-file generation.
func TestParity_Fastp_SplitByLines_MultiThread(t *testing.T) {
	bin := ensureUpstream(t)

	in := filepath.Join(t.TempDir(), "in.fq")
	genSplitFASTQ(t, in, 6000)

	for _, w := range []int{2, 3, 4} {
		w := w
		t.Run(fmt.Sprintf("w%d", w), func(t *testing.T) {
			upDir := t.TempDir()
			args := append([]string{
				"-i", in, "-o", filepath.Join(upDir, "out.fq"),
				"-S", "4000", "-w", fmt.Sprintf("%d", w),
				"-j", filepath.Join(upDir, "up.json"), "-h", filepath.Join(upDir, "up.html"),
			}, commonDisableFlags...)
			runUpstream(t, bin, args)

			goDir := t.TempDir()
			opts := permissiveOpts()
			opts.SplitByLines = 4000
			opts.Threads = w
			runGoFastpSESplit(t, in, filepath.Join(goDir, "out.fq"), opts)

			compareSplitDirs(t, upDir, goDir, "out.fq")
		})
	}
}

// TestParity_Fastp_SplitByLines_MultiThread_Filtered exercises the byFileLines
// passed-count rollover where passed != input: a length filter drops ~half the
// reads. Because the drop is purely length-based (deterministic and identical
// between upstream and the Go port), the surviving set and per-file boundaries
// must still match byte-for-byte. This is the path where the rollover count
// (readPassed) genuinely diverges from the input pack count.
func TestParity_Fastp_SplitByLines_MultiThread_Filtered(t *testing.T) {
	bin := ensureUpstream(t)

	// Half the reads are 40bp (kept), half are 20bp (dropped by -l 30).
	in := filepath.Join(t.TempDir(), "in.fq")
	genMixedLenFASTQ(t, in, 4000)

	for _, w := range []int{1, 2, 3} {
		w := w
		t.Run(fmt.Sprintf("w%d", w), func(t *testing.T) {
			upDir := t.TempDir()
			args := []string{
				"-i", in, "-o", filepath.Join(upDir, "out.fq"),
				"-S", "2000", "-w", fmt.Sprintf("%d", w),
				"-l", "30",
				"--disable_quality_filtering", "--disable_adapter_trimming",
				"-j", filepath.Join(upDir, "up.json"), "-h", filepath.Join(upDir, "up.html"),
			}
			runUpstream(t, bin, args)

			goDir := t.TempDir()
			opts := permissiveOpts()
			opts.MinLength = 30
			opts.LengthRequired = 30
			opts.SplitByLines = 2000
			opts.Threads = w
			runGoFastpSESplit(t, in, filepath.Join(goDir, "out.fq"), opts)

			compareSplitDirs(t, upDir, goDir, "out.fq")
		})
	}
}

// genMixedLenFASTQ writes n records alternating 40bp (kept by -l 30) and 20bp
// (dropped), all high quality. The drop pattern is deterministic so the
// surviving set is identical between upstream and the Go port.
func genMixedLenFASTQ(t *testing.T, path string, n int) {
	t.Helper()
	var b strings.Builder
	bases := []byte("ACGT")
	for i := 0; i < n; i++ {
		length := 40
		if i%2 == 1 {
			length = 20
		}
		b.WriteString(fmt.Sprintf("@read%d\n", i))
		seq := make([]byte, length)
		for j := range seq {
			seq[j] = bases[(i+j)%4]
		}
		b.Write(seq)
		b.WriteString("\n+\n")
		b.WriteString(strings.Repeat("I", length))
		b.WriteByte('\n')
	}
	if err := os.WriteFile(path, []byte(b.String()), 0o644); err != nil {
		t.Fatalf("write %s: %v", path, err)
	}
}

// TestParity_Fastp_SplitPE_MultiThread compares paired-end split output across
// worker-thread counts, verifying both the R1 and R2 numbered files match
// upstream byte-for-byte and stay aligned across the rollover.
func TestParity_Fastp_SplitPE_MultiThread(t *testing.T) {
	bin := ensureUpstream(t)

	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	genSplitFASTQ(t, r1, 3000)
	genSplitFASTQ(t, r2, 3000)

	for _, w := range []int{2, 3, 4} {
		w := w
		t.Run(fmt.Sprintf("w%d", w), func(t *testing.T) {
			upDir := t.TempDir()
			args := append([]string{
				"-i", r1, "-I", r2,
				"-o", filepath.Join(upDir, "out1.fq"), "-O", filepath.Join(upDir, "out2.fq"),
				"-s", "4", "-w", fmt.Sprintf("%d", w),
				"-j", filepath.Join(upDir, "up.json"), "-h", filepath.Join(upDir, "up.html"),
			}, commonDisableFlags...)
			runUpstream(t, bin, args)

			goDir := t.TempDir()
			opts := permissiveOpts()
			opts.SplitNumber = 4
			opts.Threads = w
			runGoFastpPESplit(t, r1, r2,
				filepath.Join(goDir, "out1.fq"), filepath.Join(goDir, "out2.fq"), opts)

			compareSplitDirs(t, upDir, goDir, "out1.fq")
			compareSplitDirs(t, upDir, goDir, "out2.fq")
		})
	}
}

// runGoFastpPESplit drives ProcessPairedEndSplit against the given base paths.
func runGoFastpPESplit(t *testing.T, r1, r2, outBase1, outBase2 string, opts ProcessOptions) {
	t.Helper()
	f1, err := os.Open(r1)
	if err != nil {
		t.Fatalf("open %s: %v", r1, err)
	}
	defer f1.Close()
	f2, err := os.Open(r2)
	if err != nil {
		t.Fatalf("open %s: %v", r2, err)
	}
	defer f2.Close()
	if _, err := ProcessPairedEndSplit(f1, f2, outBase1, outBase2, defaultEncoding(), opts); err != nil {
		t.Fatalf("ProcessPairedEndSplit: %v", err)
	}
}
