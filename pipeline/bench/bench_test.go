// Package bench holds the performance side of the pipeline: testing.B harnesses
// that time OUR tool binary against the vendored UPSTREAM binary on a medium or
// large fixture and report the ratio.
//
// Run with:
//
//	PIPELINE_SCALE=medium go test -bench=. ./pipeline/bench
//	PIPELINE_SCALE=large  go test -bench=. -benchtime=3x ./pipeline/bench
//
// The fixtures are generated once (and cached) on first use via the shared
// fixtures package, so a bench run reuses whatever the parity runner produced.
// Each benchmark times both sides and logs "ours/upstream" so the wins (or
// regressions) are visible; later agents add more heavy ops following this
// exact pattern (see benchPair).
//
// These are integration benchmarks: they shell out to the built binaries
// rather than calling library code, so they measure the same thing a user sees
// on the command line.
package bench

import (
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

var (
	manOnce sync.Once
	man     *fixtures.Manifest
	manErr  error
)

// benchManifest generates (or reuses) the fixture set at PIPELINE_SCALE,
// defaulting to medium for benchmarks (heavier than the small parity default).
func benchManifest(tb testing.TB) *fixtures.Manifest {
	tb.Helper()
	manOnce.Do(func() {
		scale, err := fixtures.ParseScale(os.Getenv("PIPELINE_SCALE"))
		if err != nil {
			manErr = err
			return
		}
		if os.Getenv("PIPELINE_SCALE") == "" {
			scale = fixtures.Medium
		}
		man, manErr = fixtures.Generate(fixtures.Options{Scale: scale})
	})
	if manErr != nil {
		tb.Fatalf("fixtures: %v\n(hint: ensure upstream binaries exist under reference_code/; "+
			"in a worktree, symlink them from the main checkout)", manErr)
	}
	return man
}

// cacheDir returns the shared binary cache used by both the runner and benches.
func cacheDir(tb testing.TB) string {
	tb.Helper()
	root, err := upstream.RepoRoot()
	if err != nil {
		tb.Fatal(err)
	}
	return filepath.Join(root, "pipeline", ".fixtures", "bin")
}

// benchPair times OUR binary and the UPSTREAM binary over b.N iterations of the
// same command and reports each side's mean wall-clock plus the ours/upstream
// ratio as custom metrics. This is the reusable pattern every heavy bench
// should follow.
func benchPair(b *testing.B, ourTool, upstreamKey string, ourArgs, upstreamArgs []string) {
	ourBin, err := upstream.OurBinary(ourTool, cacheDir(b))
	if err != nil {
		b.Fatal(err)
	}
	upBin, err := upstream.Binary(upstreamKey)
	if err != nil {
		b.Fatal(err)
	}

	run := func(bin string, args []string) time.Duration {
		cmd := exec.Command(bin, args...)
		cmd.Stdout = nil // discard
		cmd.Stderr = nil
		start := time.Now()
		if err := cmd.Run(); err != nil {
			b.Fatalf("%s %v: %v", bin, args, err)
		}
		return time.Since(start)
	}

	var ourTotal, upTotal time.Duration
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ourTotal += run(ourBin, ourArgs)
		upTotal += run(upBin, upstreamArgs)
	}
	b.StopTimer()

	n := float64(b.N)
	ourMs := float64(ourTotal.Milliseconds()) / n
	upMs := float64(upTotal.Milliseconds()) / n
	b.ReportMetric(ourMs, "ours_ms/op")
	b.ReportMetric(upMs, "upstream_ms/op")
	if upMs > 0 {
		b.ReportMetric(ourMs/upMs, "ratio_ours/up")
	}
}

// BenchmarkSamtoolsViewBAMtoBAM times a full BAM read + re-encode to BAM (the
// canonical heavy samtools path).
func BenchmarkSamtoolsViewBAMtoBAM(b *testing.B) {
	m := benchManifest(b)
	out := filepath.Join(b.TempDir(), "out.bam")
	benchPair(b, "samtools", "samtools",
		[]string{"view", "-b", "-o", out, m.Path("bam")},
		[]string{"view", "-b", "-o", out, m.Path("bam")},
	)
}

// BenchmarkBcftoolsViewFilter times an INFO/QUAL filter over the whole VCF.
func BenchmarkBcftoolsViewFilter(b *testing.B) {
	m := benchManifest(b)
	args := []string{"view", "-H", "-i", "QUAL>30 && INFO/DP>20", m.Path("vcf_plain")}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBedtoolsIntersect times the write-A-write-B intersection, the
// heaviest interval-output path.
func BenchmarkBedtoolsIntersect(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-a", m.Path("bed"), "-b", m.Path("bed12"), "-wa", "-wb"}
	upArgs := append([]string{"intersect"}, ourArgs...)
	benchPair(b, "bedintersect", "bedtools", ourArgs, upArgs)
}
