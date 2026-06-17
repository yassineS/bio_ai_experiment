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
// regressions) are visible.
//
// These are integration benchmarks: they shell out to the built binaries
// rather than calling library code, so they measure the same thing a user sees
// on the command line. They MEASURE; they do not assert parity — a bench fails
// only on a hard error (a non-zero exit / setup failure), never on a perf
// delta. Each our-side command still produces real output so the timing is
// meaningful.
//
// The set covers the genuinely heavy operations across the toolset:
//
//	samtools : view BAM->BAM (+ -@ threaded), view BAM->CRAM, sort, index,
//	           flagstat, stats, depth, mpileup, markdup
//	bcftools : view filter, query, norm, stats, concat, call, mpileup, +fill-tags
//	bedtools : intersect, coverage, merge, sort, genomecov, map, closest,
//	           makewindows
//	vcftools : a heavy stat pass (--freq / --site-pi / --relatedness2)
//	htslib   : bgzip compress + decompress, tabix build
//	QC       : fastp and seqtk seq throughput on the FASTQ fixture
//
// Heavy threaded variants are gated to the `large` tier (see requireLarge) so a
// plain `go test -bench` at the medium default finishes in reasonable time.
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

// requireLarge skips a benchmark unless PIPELINE_SCALE=large. The very heavy
// threaded re-encode variants only become interesting at the large tier; at the
// medium default they would just lengthen a routine `go test -bench` run.
func requireLarge(b *testing.B) {
	b.Helper()
	if fixtures.Scale(os.Getenv("PIPELINE_SCALE")) != fixtures.Large {
		b.Skipf("heavy variant: set PIPELINE_SCALE=large to run %s", b.Name())
	}
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

// prep runs an untimed setup command (resolving our or upstream binaries) and
// fails the benchmark on a non-zero exit. It is used to materialise inputs that
// one of the timed commands then consumes (e.g. an mpileup VCF for `call`).
func prep(b *testing.B, bin string, args ...string) {
	b.Helper()
	cmd := exec.Command(bin, args...)
	if out, err := cmd.CombinedOutput(); err != nil {
		b.Fatalf("setup %s %v: %v\n%s", bin, args, err, out)
	}
}

// ---------------------------------------------------------------------------
// samtools
// ---------------------------------------------------------------------------

// BenchmarkSamtoolsViewBAMtoBAM times a full BAM read + re-encode to BAM (the
// canonical heavy samtools path).
func BenchmarkSamtoolsViewBAMtoBAM(b *testing.B) {
	m := benchManifest(b)
	out := filepath.Join(b.TempDir(), "out.bam")
	args := []string{"view", "-b", "-o", out, m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsViewBAMtoBAMThreaded re-encodes BAM->BAM with 4 worker
// threads. Gated to the large tier (the threaded win only shows at scale).
func BenchmarkSamtoolsViewBAMtoBAMThreaded(b *testing.B) {
	requireLarge(b)
	m := benchManifest(b)
	out := filepath.Join(b.TempDir(), "out.bam")
	args := []string{"view", "-@", "4", "-b", "-o", out, m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsViewBAMtoCRAM times BAM->CRAM re-encoding (reference
// compression, the heaviest encode path).
func BenchmarkSamtoolsViewBAMtoCRAM(b *testing.B) {
	m := benchManifest(b)
	out := filepath.Join(b.TempDir(), "out.cram")
	args := []string{"view", "-C", "-T", m.Path("fasta"), "-o", out, m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsViewBAMtoCRAMThreaded is the threaded BAM->CRAM encode,
// gated to the large tier.
func BenchmarkSamtoolsViewBAMtoCRAMThreaded(b *testing.B) {
	requireLarge(b)
	m := benchManifest(b)
	out := filepath.Join(b.TempDir(), "out.cram")
	args := []string{"view", "-@", "4", "-C", "-T", m.Path("fasta"), "-o", out, m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsSort times a coordinate sort + BAM re-encode.
func BenchmarkSamtoolsSort(b *testing.B) {
	m := benchManifest(b)
	out := filepath.Join(b.TempDir(), "sorted.bam")
	args := []string{"sort", "-o", out, m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsIndex times BAI index construction over the sorted BAM. Each
// side indexes a private copy so the two never race on the same .bai.
func BenchmarkSamtoolsIndex(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	ourCopy := copyFixture(b, m.Path("bam"), filepath.Join(dir, "ours.bam"))
	upCopy := copyFixture(b, m.Path("bam"), filepath.Join(dir, "up.bam"))
	benchPair(b, "samtools", "samtools",
		[]string{"index", ourCopy},
		[]string{"index", upCopy},
	)
}

// BenchmarkSamtoolsFlagstat times the flag-count pass over the whole BAM.
func BenchmarkSamtoolsFlagstat(b *testing.B) {
	m := benchManifest(b)
	args := []string{"flagstat", m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsStats times the full `stats` summary pass.
func BenchmarkSamtoolsStats(b *testing.B) {
	m := benchManifest(b)
	args := []string{"stats", m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsDepth times the per-position depth scan over the BAM.
func BenchmarkSamtoolsDepth(b *testing.B) {
	m := benchManifest(b)
	args := []string{"depth", m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsMpileup times the per-base pileup against the reference.
func BenchmarkSamtoolsMpileup(b *testing.B) {
	m := benchManifest(b)
	args := []string{"mpileup", "-f", m.Path("fasta"), m.Path("bam")}
	benchPair(b, "samtools", "samtools", args, args)
}

// BenchmarkSamtoolsMarkdup times duplicate marking (two-pass) on the
// coordinate-sorted BAM, writing BAM to a private temp file per side.
func BenchmarkSamtoolsMarkdup(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	benchPair(b, "samtools", "samtools",
		[]string{"markdup", m.Path("bam"), filepath.Join(dir, "ours.bam")},
		[]string{"markdup", m.Path("bam"), filepath.Join(dir, "up.bam")},
	)
}

// ---------------------------------------------------------------------------
// bcftools
// ---------------------------------------------------------------------------

// BenchmarkBcftoolsViewFilter times an INFO/QUAL filter over the whole VCF.
func BenchmarkBcftoolsViewFilter(b *testing.B) {
	m := benchManifest(b)
	args := []string{"view", "-H", "-i", "QUAL>30 && INFO/DP>20", m.Path("vcf_plain")}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBcftoolsQuery times a formatted per-record extraction with a filter.
func BenchmarkBcftoolsQuery(b *testing.B) {
	m := benchManifest(b)
	args := []string{"query", "-i", "QUAL>20", "-f", `%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n`, m.Path("vcf_plain")}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBcftoolsNorm times left-alignment + multiallelic splitting against
// the reference.
func BenchmarkBcftoolsNorm(b *testing.B) {
	m := benchManifest(b)
	args := []string{"norm", "-f", m.Path("fasta"), "-m-", "-Ov", m.Path("vcf")}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBcftoolsStats times the full stats pass over the VCF.
func BenchmarkBcftoolsStats(b *testing.B) {
	m := benchManifest(b)
	args := []string{"stats", m.Path("vcf_plain")}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBcftoolsConcat times concatenating the bgzipped VCF with itself to
// uncompressed VCF (the decode+re-emit path). -a (allow overlaps) is required
// because the two inputs share coordinate ranges; upstream errors otherwise.
func BenchmarkBcftoolsConcat(b *testing.B) {
	m := benchManifest(b)
	args := []string{"concat", "-a", "-Ov", m.Path("vcf"), m.Path("vcf")}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBcftoolsMpileup times genotype-likelihood mpileup over the BAM,
// writing uncompressed BCF to a private temp file per side.
func BenchmarkBcftoolsMpileup(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	benchPair(b, "bcftools", "bcftools",
		[]string{"mpileup", "-f", m.Path("fasta"), "-Ou", "-o", filepath.Join(dir, "ours.bcf"), m.Path("bam")},
		[]string{"mpileup", "-f", m.Path("fasta"), "-Ou", "-o", filepath.Join(dir, "up.bcf"), m.Path("bam")},
	)
}

// BenchmarkBcftoolsCall times variant calling. The genotype-likelihood input is
// produced once (untimed) via upstream mpileup so both sides consume an
// identical VCF; the bench measures only the `call` step.
func BenchmarkBcftoolsCall(b *testing.B) {
	m := benchManifest(b)
	upBin, err := upstream.Binary("bcftools")
	if err != nil {
		b.Fatal(err)
	}
	pile := filepath.Join(b.TempDir(), "pileup.vcf")
	prep(b, upBin, "mpileup", "-f", m.Path("fasta"), "-Ov", "-o", pile, m.Path("bam"))
	args := []string{"call", "-mv", "-Ov", pile}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// BenchmarkBcftoolsFillTags times the +fill-tags plugin recomputing INFO tags
// over the multi-sample VCF. Upstream loads the plugin as a shared object from
// BCFTOOLS_PLUGINS (shipped under reference_code/bcftools/plugins); our port has
// it built in and ignores the variable. Setting it in the inherited process
// environment is harmless to our side and required for upstream.
func BenchmarkBcftoolsFillTags(b *testing.B) {
	m := benchManifest(b)
	upBin, err := upstream.Binary("bcftools")
	if err != nil {
		b.Fatal(err)
	}
	// The plugins live next to the *real* binary; in a worktree the bcftools
	// path is a symlink into the main checkout, so resolve it before locating
	// the plugins/ dir.
	if real, rerr := filepath.EvalSymlinks(upBin); rerr == nil {
		upBin = real
	}
	b.Setenv("BCFTOOLS_PLUGINS", filepath.Join(filepath.Dir(upBin), "plugins"))
	args := []string{"+fill-tags", m.Path("vcf_multi"), "-Ov"}
	benchPair(b, "bcftools", "bcftools", args, args)
}

// ---------------------------------------------------------------------------
// bedtools (each bed* tool is its own binary; only upstream takes the subcommand)
// ---------------------------------------------------------------------------

// BenchmarkBedtoolsIntersect times the write-A-write-B intersection, the
// heaviest interval-output path.
func BenchmarkBedtoolsIntersect(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-a", m.Path("bed"), "-b", m.Path("bed12"), "-wa", "-wb"}
	upArgs := append([]string{"intersect"}, ourArgs...)
	benchPair(b, "bedintersect", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsCoverage times per-A coverage computation against B.
func BenchmarkBedtoolsCoverage(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-a", m.Path("bed"), "-b", m.Path("bed12")}
	upArgs := append([]string{"coverage"}, ourArgs...)
	benchPair(b, "bedcoverage", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsMerge times merging the overlapping/adjacent intervals.
func BenchmarkBedtoolsMerge(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-i", m.Path("bed")}
	upArgs := append([]string{"merge"}, ourArgs...)
	benchPair(b, "bedmerge", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsSort times the lexicographic interval sort.
func BenchmarkBedtoolsSort(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-i", m.Path("bed")}
	upArgs := append([]string{"sort"}, ourArgs...)
	benchPair(b, "bedsort", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsGenomecov times the per-genome coverage histogram.
func BenchmarkBedtoolsGenomecov(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-i", m.Path("bed"), "-g", m.Path("genome")}
	upArgs := append([]string{"genomecov"}, ourArgs...)
	benchPair(b, "bedgenomecov", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsMap times mapping/aggregating B values onto A.
func BenchmarkBedtoolsMap(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-a", m.Path("bed"), "-b", m.Path("bed12"), "-c", "5", "-o", "mean"}
	upArgs := append([]string{"map"}, ourArgs...)
	benchPair(b, "bedmap", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsClosest times the nearest-feature search of A against B.
func BenchmarkBedtoolsClosest(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-a", m.Path("bed"), "-b", m.Path("bed12")}
	upArgs := append([]string{"closest"}, ourArgs...)
	benchPair(b, "bedclosest", "bedtools", ourArgs, upArgs)
}

// BenchmarkBedtoolsMakewindows times tiling the genome into fixed windows.
func BenchmarkBedtoolsMakewindows(b *testing.B) {
	m := benchManifest(b)
	ourArgs := []string{"-g", m.Path("genome"), "-w", "1000"}
	upArgs := append([]string{"makewindows"}, ourArgs...)
	benchPair(b, "bedmakewindows", "bedtools", ourArgs, upArgs)
}

// ---------------------------------------------------------------------------
// vcftools (writes to a --out prefix; each side gets its own temp prefix)
// ---------------------------------------------------------------------------

// BenchmarkVcftoolsFreq times the allele-frequency stat pass.
func BenchmarkVcftoolsFreq(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	benchPair(b, "vcftools", "vcftools",
		[]string{"--vcf", m.Path("vcf_plain"), "--freq", "--out", filepath.Join(dir, "ours")},
		[]string{"--vcf", m.Path("vcf_plain"), "--freq", "--out", filepath.Join(dir, "up")},
	)
}

// BenchmarkVcftoolsSitePi times the per-site nucleotide-diversity stat pass.
func BenchmarkVcftoolsSitePi(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	benchPair(b, "vcftools", "vcftools",
		[]string{"--vcf", m.Path("vcf_plain"), "--site-pi", "--out", filepath.Join(dir, "ours")},
		[]string{"--vcf", m.Path("vcf_plain"), "--site-pi", "--out", filepath.Join(dir, "up")},
	)
}

// BenchmarkVcftoolsRelatedness2 times the O(N^2-per-site) KING relatedness pass
// over the multi-sample VCF (the heaviest vcftools stat path here).
func BenchmarkVcftoolsRelatedness2(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	benchPair(b, "vcftools", "vcftools",
		[]string{"--vcf", m.Path("vcf_multi_plain"), "--relatedness2", "--out", filepath.Join(dir, "ours")},
		[]string{"--vcf", m.Path("vcf_multi_plain"), "--relatedness2", "--out", filepath.Join(dir, "up")},
	)
}

// ---------------------------------------------------------------------------
// htslib: bgzip / tabix
// ---------------------------------------------------------------------------

// BenchmarkBgzipCompress times BGZF compression of the plain VCF to a private
// temp file per side (binary output never goes to stdout).
func BenchmarkBgzipCompress(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	ourSrc := copyFixture(b, m.Path("vcf_plain"), filepath.Join(dir, "ours.vcf"))
	upSrc := copyFixture(b, m.Path("vcf_plain"), filepath.Join(dir, "up.vcf"))
	// -kf keeps the source and overwrites any stale .gz, so b.N iterations
	// recompress the same input deterministically.
	benchPair(b, "bgzip", "bgzip",
		[]string{"-kf", ourSrc},
		[]string{"-kf", upSrc},
	)
}

// BenchmarkBgzipDecompress times BGZF decompression of the bgzipped VCF to a
// private temp file per side.
func BenchmarkBgzipDecompress(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	ourSrc := copyFixture(b, m.Path("vcf"), filepath.Join(dir, "ours.vcf.gz"))
	upSrc := copyFixture(b, m.Path("vcf"), filepath.Join(dir, "up.vcf.gz"))
	benchPair(b, "bgzip", "bgzip",
		[]string{"-dkf", ourSrc},
		[]string{"-dkf", upSrc},
	)
}

// BenchmarkTabixBuild times tabix VCF index construction over a private copy of
// the bgzipped VCF per side (the .tbi is written next to the input).
func BenchmarkTabixBuild(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	ourSrc := copyFixture(b, m.Path("vcf"), filepath.Join(dir, "ours.vcf.gz"))
	upSrc := copyFixture(b, m.Path("vcf"), filepath.Join(dir, "up.vcf.gz"))
	benchPair(b, "tabix", "tabix",
		[]string{"-p", "vcf", "-f", ourSrc},
		[]string{"-p", "vcf", "-f", upSrc},
	)
}

// ---------------------------------------------------------------------------
// QC throughput: fastp / seqtk
// ---------------------------------------------------------------------------

// BenchmarkFastp times default QC/adapter trimming over the single-end FASTQ
// fixture, writing trimmed reads + reports to a private temp dir per side.
func BenchmarkFastp(b *testing.B) {
	m := benchManifest(b)
	dir := b.TempDir()
	benchPair(b, "fastp", "fastp",
		[]string{"-i", m.Path("fastq"), "-o", filepath.Join(dir, "ours.fq"),
			"--json", filepath.Join(dir, "ours.json"), "--html", filepath.Join(dir, "ours.html")},
		[]string{"-i", m.Path("fastq"), "-o", filepath.Join(dir, "up.fq"),
			"--json", filepath.Join(dir, "up.json"), "--html", filepath.Join(dir, "up.html")},
	)
}

// BenchmarkSeqtkSeq times FASTQ->FASTA conversion throughput over the FASTQ
// fixture (output discarded to stdout, the canonical seqtk seq path).
func BenchmarkSeqtkSeq(b *testing.B) {
	m := benchManifest(b)
	args := []string{"seq", "-A", m.Path("fastq")}
	benchPair(b, "seqtk", "seqtk", args, args)
}

// copyFixture copies a fixture file to dst (so b.N iterations of an in-place /
// alongside-output command operate on a private, deterministic input) and
// returns dst.
func copyFixture(b *testing.B, src, dst string) string {
	b.Helper()
	data, err := os.ReadFile(src)
	if err != nil {
		b.Fatalf("read fixture %s: %v", src, err)
	}
	if err := os.WriteFile(dst, data, 0o644); err != nil {
		b.Fatalf("write %s: %v", dst, err)
	}
	return dst
}
