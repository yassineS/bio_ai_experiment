package samtools

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// This file exercises the real `-@/--threads` parallelism added to the
// BAM-writing subcommands (view, sort, markdup). The contract under test is
// that turning on multiple BGZF compression workers never changes the data:
// the decompressed BAM body is byte-identical regardless of thread count, and
// it decodes to the same records as the upstream C `samtools` binary. The
// compressed bytes themselves may differ (block boundaries depend on buffering
// timing), so every comparison is made on the *decoded* BAM body, mirroring
// the approach the bgzip MultiWriter parity test takes.

// decompressBGZF fully inflates a BGZF stream and returns the underlying
// plaintext (here, a BAM body: magic + header + records). Two BGZF streams that
// decode to the same bytes carry the same data even when their block framing
// differs, which is exactly the invariant multi-threaded compression must keep.
func decompressBGZF(t *testing.T, b []byte) []byte {
	t.Helper()
	r, err := bgzip.NewReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("open BGZF reader: %v", err)
	}
	out, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("inflate BGZF: %v", err)
	}
	return out
}

// bigSAM builds a synthetic SAM stream with n records spanning enough bytes to
// fill many BGZF blocks (MaxBlockSize is ~64 KiB), so the worker pool is
// genuinely exercised — a single-block fixture would not distinguish serial
// from parallel framing.
func bigSAM(n int) []byte {
	var sb strings.Builder
	sb.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	sb.WriteString("@SQ\tSN:chr1\tLN:100000000\n")
	seq := strings.Repeat("ACGT", 25) // 100 bp
	qual := strings.Repeat("I", 100)
	for i := 0; i < n; i++ {
		pos := (i*37)%99999000 + 1
		fmt.Fprintf(&sb, "read%06d\t0\tchr1\t%d\t60\t100M\t*\t0\t0\t%s\t%s\tNM:i:0\n",
			i, pos, seq, qual)
	}
	return []byte(sb.String())
}

// viewBAM streams in (SAM bytes) through View producing BAM with the given
// thread count, returning the raw BGZF output bytes.
func viewBAM(t *testing.T, in []byte, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := View(bytes.NewReader(in), &out, ViewOptions{OutputBAM: true, Threads: threads}); err != nil {
		t.Fatalf("View -@ %d: %v", threads, err)
	}
	return out.Bytes()
}

// TestThreads_View_BAMByteIdenticalAcrossThreadCounts asserts that View's BAM
// output decodes to identical bytes for -@ 1 and several higher thread counts.
func TestThreads_View_BAMByteIdenticalAcrossThreadCounts(t *testing.T) {
	in := bigSAM(20000)
	base := decompressBGZF(t, viewBAM(t, in, 1))
	for _, threads := range []int{2, 3, 4, 8} {
		got := decompressBGZF(t, viewBAM(t, in, threads))
		if !bytes.Equal(base, got) {
			t.Fatalf("view -@ %d decoded BAM differs from -@ 1 (%d vs %d bytes)",
				threads, len(got), len(base))
		}
	}
}

// TestThreads_Sort_BAMByteIdenticalAcrossThreadCounts asserts the same
// invariant for sort, which also compresses temporary shards in parallel.
func TestThreads_Sort_BAMByteIdenticalAcrossThreadCounts(t *testing.T) {
	in := bigSAM(20000)
	sortBAM := func(threads int) []byte {
		var out bytes.Buffer
		// A small MaxMem forces the external-merge path so the parallel shard
		// writer is exercised, not just the final-output writer.
		opts := SortOptions{Order: SortCoordinate, Threads: threads, MaxMemBytes: 256 * 1024}
		if err := Sort(bytes.NewReader(in), &out, opts); err != nil {
			t.Fatalf("Sort -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := decompressBGZF(t, sortBAM(1))
	for _, threads := range []int{2, 4, 8} {
		got := decompressBGZF(t, sortBAM(threads))
		if !bytes.Equal(base, got) {
			t.Fatalf("sort -@ %d decoded BAM differs from -@ 1 (%d vs %d bytes)",
				threads, len(got), len(base))
		}
	}
}

// TestThreads_Markdup_BAMByteIdenticalAcrossThreadCounts asserts the invariant
// for markdup's BAM output.
func TestThreads_Markdup_BAMByteIdenticalAcrossThreadCounts(t *testing.T) {
	in := bigSAM(20000)
	markdupBAM := func(threads int) []byte {
		var out bytes.Buffer
		opener := func() (io.ReadCloser, error) {
			return io.NopCloser(bytes.NewReader(in)), nil
		}
		if _, err := Markdup(opener, &out, MarkdupOptions{Threads: threads}); err != nil {
			t.Fatalf("Markdup -@ %d: %v", threads, err)
		}
		return out.Bytes()
	}
	base := decompressBGZF(t, markdupBAM(1))
	for _, threads := range []int{2, 4, 8} {
		got := decompressBGZF(t, markdupBAM(threads))
		if !bytes.Equal(base, got) {
			t.Fatalf("markdup -@ %d decoded BAM differs from -@ 1 (%d vs %d bytes)",
				threads, len(got), len(base))
		}
	}
}

// recordKeys decodes a BAM stream and returns one "QNAME\tFLAG\tRNAME\tPOS"
// key per record, for comparing record sets against the upstream binary
// independent of header key ordering.
func recordKeys(t *testing.T, bamBytes []byte) []string {
	t.Helper()
	r, err := sam.NewReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatalf("open BAM reader: %v", err)
	}
	var keys []string
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("read BAM record: %v", err)
		}
		keys = append(keys, fmt.Sprintf("%s\t%d\t%s\t%d", rec.QName, rec.Flag, rec.RName, rec.Pos))
	}
	return keys
}

// TestThreads_View_UpstreamParity confirms that our multi-threaded BAM output
// decodes to the same record set the upstream C `samtools` produces from the
// same SAM. The upstream binary is built live from the vendored submodule
// (t.Fatalf, never t.Skip).
func TestThreads_View_UpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)
	in := bigSAM(5000)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, in, 0o644); err != nil {
		t.Fatalf("write SAM: %v", err)
	}

	// Upstream: SAM -> BAM with several threads.
	upBAM, err := run(dir, bin, "view", "-@", "4", "-b", "-o", filepath.Join(dir, "up.bam"), samPath)
	if err != nil {
		t.Fatalf("upstream samtools view: %v\n%s", err, upBAM)
	}
	upBytes, err := os.ReadFile(filepath.Join(dir, "up.bam"))
	if err != nil {
		t.Fatalf("read upstream BAM: %v", err)
	}
	wantKeys := recordKeys(t, upBytes)

	gotKeys := recordKeys(t, viewBAM(t, in, 4))
	if len(gotKeys) != len(wantKeys) {
		t.Fatalf("record count mismatch: ours %d, upstream %d", len(gotKeys), len(wantKeys))
	}
	for i := range wantKeys {
		if gotKeys[i] != wantKeys[i] {
			t.Fatalf("record %d mismatch: ours %q upstream %q", i, gotKeys[i], wantKeys[i])
		}
	}
}

// TestThreads_Sort_UpstreamParity confirms our multi-threaded sort produces the
// same coordinate-sorted record set as upstream samtools sort.
func TestThreads_Sort_UpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)
	in := bigSAM(5000)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, in, 0o644); err != nil {
		t.Fatalf("write SAM: %v", err)
	}
	upOut := filepath.Join(dir, "up.sorted.bam")
	if out, err := run(dir, bin, "sort", "-@", "4", "-o", upOut, samPath); err != nil {
		t.Fatalf("upstream samtools sort: %v\n%s", err, out)
	}
	upBytes, err := os.ReadFile(upOut)
	if err != nil {
		t.Fatalf("read upstream sorted BAM: %v", err)
	}
	wantKeys := recordKeys(t, upBytes)

	var ourOut bytes.Buffer
	if err := Sort(bytes.NewReader(in), &ourOut, SortOptions{Order: SortCoordinate, Threads: 4, MaxMemBytes: 256 * 1024}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	gotKeys := recordKeys(t, ourOut.Bytes())
	if len(gotKeys) != len(wantKeys) {
		t.Fatalf("record count mismatch: ours %d, upstream %d", len(gotKeys), len(wantKeys))
	}
	for i := range wantKeys {
		if gotKeys[i] != wantKeys[i] {
			t.Fatalf("sorted record %d mismatch: ours %q upstream %q", i, gotKeys[i], wantKeys[i])
		}
	}
}

// BenchmarkThreads_ViewBAM measures BAM-output throughput at 1 vs 4 BGZF
// workers, demonstrating the speedup of the parallel path. Run with
// `go test -run=^$ -bench=BenchmarkThreads_ViewBAM ./tools/samtools/...`.
func BenchmarkThreads_ViewBAM(b *testing.B) {
	in := bigSAM(50000)
	for _, threads := range []int{1, 2, 4, 8} {
		b.Run(fmt.Sprintf("threads=%d", threads), func(b *testing.B) {
			b.SetBytes(int64(len(in)))
			for i := 0; i < b.N; i++ {
				var out bytes.Buffer
				if _, err := View(bytes.NewReader(in), &out, ViewOptions{OutputBAM: true, Threads: threads}); err != nil {
					b.Fatalf("View: %v", err)
				}
			}
		})
	}
}
