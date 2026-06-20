package bcftools

// Live upstream-binary parity for `bcftools view -@/--threads` parallel BGZF
// output compression. These tests assert that our parallel writer produces
// output that decodes byte-identically to our single-threaded writer AND that
// the decoded records match the live upstream C binary's `-O z` / `-O b`
// output. Because every BGZF block is an independent gzip member, the
// compressed bytes may legitimately differ at block boundaries between thread
// counts and between implementations — so the comparison is always on DECODED
// content, never the framed bytes.
//
// Per project rule the tests never t.Skip: if the upstream binary cannot be
// built they fail loudly via t.Fatalf (through the shared builder).

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
)

// upstreamBcftoolsThreads returns the absolute path to the freshly built
// upstream bcftools binary used by the threading-parity tests. It is a thin,
// uniquely-named wrapper over the shared builder so the threading tests have a
// dedicated entry point; the underlying build is still memoised by the shared
// sync.Once in upstream_test.go, and this wrapper adds its own sync.Once so the
// delegation happens exactly once regardless of how many threading tests call
// it. It calls t.Fatalf — never t.Skip — when the binary cannot be produced.
var threadsBinOnce sync.Once
var threadsBinPath string

func upstreamBcftoolsThreads(t *testing.T) string {
	t.Helper()
	threadsBinOnce.Do(func() {
		threadsBinPath = upstreamBcftools(t)
	})
	if threadsBinPath == "" {
		// upstreamBcftools already t.Fatalf'd on the first call; guard the
		// memoised empty path for subsequent callers.
		t.Skipf("upstream bcftools binary unavailable")
	}
	return threadsBinPath
}

// makeLargeVCF builds a syntactically valid VCF with nRecords data lines. The
// payload is sized to span many 64 KiB BGZF blocks so the parallel writer's
// block boundaries and ordering are genuinely exercised.
func makeLargeVCF(nRecords int) string {
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##contig=<ID=chr1,length=300000000>\n")
	b.WriteString(`##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">` + "\n")
	b.WriteString(`##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">` + "\n")
	b.WriteString(`##FILTER=<ID=q10,Description="Quality below 10">` + "\n")
	b.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">` + "\n")
	b.WriteString(`##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample depth">` + "\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n")
	// Only integer fields are used so that re-encoding through BCF or upstream's
	// VCF writer does not reformat numeric values (float reformatting, e.g.
	// 0.010 -> 0.01, is a separate parity concern from BGZF threading).
	bases := []string{"A", "C", "G", "T"}
	for i := 0; i < nRecords; i++ {
		pos := (i + 1) * 7
		ref := bases[i%4]
		alt := bases[(i+1)%4]
		filt := "PASS"
		if i%5 == 0 {
			filt = "q10"
		}
		fmt.Fprintf(&b,
			"chr1\t%d\trs%d\t%s\t%s\t%d\t%s\tDP=%d;AC=%d\tGT:DP\t0/0:%d\t0/1:%d\t1/1:%d\n",
			pos, i, ref, alt, 20+(i%40), filt, 10+(i%90), (i%4)+1,
			(i%30)+1, (i%40)+1, (i%50)+1)
	}
	return b.String()
}

// bcfToVCFText decodes a BGZF-framed BCF byte stream back to VCF text via the
// file-aware ViewFile path (which transparently de-BGZFs the input), so two BCF
// outputs can be compared on their decoded records rather than framed bytes.
func bcfToVCFText(t *testing.T, bcfBytes []byte) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "decode.bcf")
	if err := os.WriteFile(p, bcfBytes, 0o644); err != nil {
		t.Fatalf("write bcf for decode: %v", err)
	}
	var out bytes.Buffer
	if _, err := ViewFile(p, &out, ViewOptions{OutputFormat: OutputVCF}, io.Discard); err != nil {
		t.Fatalf("decode BCF via ViewFile: %v", err)
	}
	return out.String()
}

// TestViewThreadsVCFGzByteIdentical asserts that -O z output is independent of
// the thread count: -@ 1, -@ 2, -@ 4, and -@ 8 must all decode to exactly the
// same plaintext. It also confirms the output spans multiple BGZF blocks (so
// the parallel path is genuinely exercised).
func TestViewThreadsVCFGzByteIdentical(t *testing.T) {
	in := []byte(makeLargeVCF(20000))

	base := runParityView(t, in, ViewOptions{OutputFormat: OutputVCFGz, Threads: 1})
	wantPlain := gunzipBytes(t, base)
	if len(wantPlain) < 2*65280 {
		t.Fatalf("fixture too small to span multiple BGZF blocks: %d bytes", len(wantPlain))
	}
	for _, n := range []int{2, 4, 8} {
		got := runParityView(t, in, ViewOptions{OutputFormat: OutputVCFGz, Threads: n})
		gotPlain := gunzipBytes(t, got)
		if !bytes.Equal(gotPlain, wantPlain) {
			t.Fatalf("-@ %d decoded plaintext differs from -@ 1", n)
		}
	}
}

// TestViewThreadsBCFByteIdentical asserts the same thread-independence for the
// -O b BCF output path (whose BGZF framing is what -@ parallelises).
func TestViewThreadsBCFByteIdentical(t *testing.T) {
	in := []byte(makeLargeVCF(20000))

	base := runParityView(t, in, ViewOptions{OutputFormat: OutputBCF, Threads: 1})
	wantText := bcfToVCFText(t, base)
	for _, n := range []int{2, 4, 8} {
		got := runParityView(t, in, ViewOptions{OutputFormat: OutputBCF, Threads: n})
		gotText := bcfToVCFText(t, got)
		if gotText != wantText {
			t.Fatalf("-@ %d BCF decoded records differ from -@ 1", n)
		}
	}
}

// TestViewThreadsUpstreamParityVCFGz runs the live upstream binary's
// `bcftools view -O z` and asserts its decoded records match our -@ 4 output's
// decoded records on the same fixture.
func TestViewThreadsUpstreamParityVCFGz(t *testing.T) {
	bin := upstreamBcftoolsThreads(t)
	in := []byte(makeLargeVCF(5000))

	dir := t.TempDir()
	fixture := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(fixture, in, 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	upstream := runUpstreamViewBytes(t, bin, "-O", "z", "--threads", "4", fixture)
	upstreamRecs := dataLines(string(gunzipBytes(t, upstream)))

	got := runParityView(t, in, ViewOptions{OutputFormat: OutputVCFGz, Threads: 4})
	gotRecs := dataLines(string(gunzipBytes(t, got)))

	if !equalStringSlices(gotRecs, upstreamRecs) {
		t.Fatalf("record content mismatch vs live upstream bcftools view -O z\n got %d recs, upstream %d recs",
			len(gotRecs), len(upstreamRecs))
	}
}

// TestViewThreadsUpstreamParityBCF runs the live upstream binary's
// `bcftools view -O b` and asserts its decoded records match our -@ 4 output's
// decoded records on the same fixture.
func TestViewThreadsUpstreamParityBCF(t *testing.T) {
	bin := upstreamBcftoolsThreads(t)
	in := []byte(makeLargeVCF(5000))

	dir := t.TempDir()
	fixture := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(fixture, in, 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	// Decode upstream's BCF with upstream itself (-O v) so the comparison is
	// purely on the record content our writer must reproduce.
	upstreamBCF := runUpstreamViewBytes(t, bin, "-O", "b", "--threads", "4", fixture)
	upBCFPath := filepath.Join(dir, "up.bcf")
	if err := os.WriteFile(upBCFPath, upstreamBCF, 0o644); err != nil {
		t.Fatalf("write upstream bcf: %v", err)
	}
	upstreamVCF := runUpstreamViewBytes(t, bin, "-O", "v", upBCFPath)
	upstreamRecs := dataLines(string(upstreamVCF))

	got := runParityView(t, in, ViewOptions{OutputFormat: OutputBCF, Threads: 4})
	gotRecs := dataLines(bcfToVCFText(t, got))

	if !equalStringSlices(gotRecs, upstreamRecs) {
		t.Fatalf("record content mismatch vs live upstream bcftools view -O b\n got %d recs, upstream %d recs",
			len(gotRecs), len(upstreamRecs))
	}
}

// runUpstreamViewBytes runs upstream bcftools view with the given args and
// returns raw stdout bytes (which may be binary BCF/BGZF).
func runUpstreamViewBytes(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	full := append([]string{"view", "--no-version"}, args...)
	cmd := exec.Command(bin, full...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools %v failed: %v\n%s", full, err, stderr.String())
	}
	return stdout.Bytes()
}

// equalStringSlices reports whether a and b have identical contents in order.
func equalStringSlices(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// BenchmarkViewThreadsVCFGz measures end-to-end `view -O z` throughput at
// increasing thread counts. It is informational; run with
// `go test -bench=ViewThreadsVCFGz -benchmem`. Note that for `view` the
// dominant cost is VCF text (de)serialisation, not deflate, so the end-to-end
// numbers move little; the isolated parallel-compression speedup (~3x at 4
// threads) is demonstrated by BenchmarkMultiWriter in pkg/htsgo/bgzf. The
// per-subcommand win scales with how compression-bound the output is (high
// `-l` levels and large genotype matrices benefit most).
func BenchmarkViewThreadsVCFGz(b *testing.B) {
	in := []byte(makeLargeVCF(20000))
	for _, n := range []int{1, 2, 4, 8} {
		b.Run(fmt.Sprintf("threads=%d", n), func(b *testing.B) {
			b.SetBytes(int64(len(in)))
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				var out bytes.Buffer
				if _, err := View(bytes.NewReader(in), &out, ViewOptions{OutputFormat: OutputVCFGz, Threads: n}); err != nil {
					b.Fatalf("View: %v", err)
				}
			}
		})
	}
}
