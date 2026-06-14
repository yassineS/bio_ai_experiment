package bcftools

// Live-binary oracle tests.
//
// Every test in this file invokes the genuine upstream bcftools binary
// vendored under reference_code/bcftools/bcftools AND the local Go port
// (built once in TestMain into a temp dir) on the same fixture, then
// asserts byte-equality of stdout (modulo provenance lines and a few
// well-documented stderr-progress quirks).
//
// The motivation is concrete: the existing parity tests compare against
// stale vendored expected files. This file closes the loop so any
// future divergence from real upstream behaviour fails CI immediately,
// rather than silently passing because both sides regressed in lock-step.
//
// Subcommand coverage matrix is documented in docs/PARITY_ROADMAP.md
// under "Live oracle coverage". When upstream cannot be located the
// suite t.Skip's so CI without submodules still passes.

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

// liveBinPath returns the absolute path to the upstream bcftools binary,
// building it from the reference_code/bcftools submodule on first use. It
// returns "" only when the submodule is not checked out; genuine build
// failures are surfaced via t.Fatalf by upstreamBcftoolsConvertGen.
func liveBinPath(t *testing.T) string {
	t.Helper()
	return upstreamBcftoolsConvertGen(t)
}

// ourBinPath is set by TestMain (live_oracle_main_test.go) to the path of
// the locally-built port binary. Empty when the build failed.
var ourBinPath string

// requireLive returns the upstream and local port binaries. Per the
// env-guard policy (PR #294) it t.Fatalf's with an exact init/build hint
// when a dependency is absent, rather than silently skipping: the
// upstream submodule can be checked out and built here, and a failure to
// build our own port is a genuine test failure.
func requireLive(t *testing.T) (live, ours string) {
	t.Helper()
	live = liveBinPath(t)
	if live == "" {
		t.Fatalf("reference_code/bcftools submodule not checked out; run `git submodule update --init --recursive reference_code/htslib reference_code/bcftools` to enable the live oracle")
	}
	if ourBinPath == "" {
		t.Fatalf("local bcftools port binary failed to build (see TestMain stderr `go build ../../cmd/bcftools`)")
	}
	return live, ourBinPath
}

// runBin invokes a binary with the given args, returns stdout, fails on
// non-zero exit. Stderr is suppressed (we only oracle-compare stdout).
func runBin(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = io.Discard
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v", bin, args, err)
	}
	return out.Bytes()
}

// livePluginDir returns the absolute path to the vendored compiled
// plugins directory (the .so files that the upstream binary dlopen's),
// or "" when it is unavailable. The upstream binary needs
// BCFTOOLS_PLUGINS pointed here to run `+fill-tags` / `+mendelian2`.
func livePluginDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..",
		"reference_code", "bcftools", "plugins"))
	if err != nil {
		return ""
	}
	if fi, err := os.Stat(abs); err != nil || !fi.IsDir() {
		return ""
	}
	return abs
}

// runBinEnv is runBin with extra environment entries (e.g.
// BCFTOOLS_PLUGINS) appended to the inherited environment.
func runBinEnv(t *testing.T, bin string, env []string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Env = append(os.Environ(), env...)
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = io.Discard
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v", bin, args, err)
	}
	return out.Bytes()
}

// Provenance lines that legitimately differ between upstream and our port
// (they encode version strings and absolute paths). Stripping them keeps
// the oracle focused on actual semantic divergence.
var (
	// ##bcftools_<name>=... and ##bcftools/<name>=... — both forms are
	// emitted by upstream depending on the subcommand / plugin path.
	provVersionRE   = regexp.MustCompile(`(?m)^##bcftools[_/][^=]+=.*\n`)
	provVersionRE2  = regexp.MustCompile(`(?m)^##bcftoolsVersion=.*\n`)
	provVersionRE3  = regexp.MustCompile(`(?m)^##bcftoolsCommand=.*\n`)
	provReferenceRE = regexp.MustCompile(`(?m)^##reference=.*\n`)
	// Comment-style provenance: roh / gtcheck / stats banners. Match
	// the permissive form `# This file was produced by` followed by
	// anything (including a `:` directly after `by`) so the regex
	// covers both `produced by bcftools` and `produced by: bcftools roh`.
	provStatsRE1   = regexp.MustCompile(`(?m)^# This file was produced by.*\n`)
	provStatsRE2   = regexp.MustCompile(`(?m)^# The command line was:.*\n`)
	provStatsRE3   = regexp.MustCompile(`(?m)^# and the working directory was:.*\n`)
	provStatsRE4   = regexp.MustCompile(`(?m)^# \t .*\n`)
	provBlankHash  = regexp.MustCompile(`(?m)^#\n`)
	provStatsTime  = regexp.MustCompile(`(?m)^INFO\tTime required.*\n`)
	provStatsTime2 = regexp.MustCompile(`(?m)^INFO\trun-time.*\n`)
)

// stripProvenance removes the upstream provenance noise that always
// differs between two invocations: version strings, command lines, and
// absolute reference paths.
func stripProvenanceBytes(b []byte) []byte {
	b = provVersionRE.ReplaceAll(b, nil)
	b = provVersionRE2.ReplaceAll(b, nil)
	b = provVersionRE3.ReplaceAll(b, nil)
	b = provReferenceRE.ReplaceAll(b, nil)
	b = provStatsRE1.ReplaceAll(b, nil)
	b = provStatsRE2.ReplaceAll(b, nil)
	b = provStatsRE3.ReplaceAll(b, nil)
	b = provStatsRE4.ReplaceAll(b, nil)
	b = provBlankHash.ReplaceAll(b, nil)
	b = provStatsTime.ReplaceAll(b, nil)
	b = provStatsTime2.ReplaceAll(b, nil)
	return b
}

// gunzipBytes BGZF/gzip-decodes b. Used to oracle-compare -Oz outputs
// without depending on the block-level encoding (BGZF block boundaries
// are an implementation detail).
func gunzipBytes(t *testing.T, b []byte) []byte {
	t.Helper()
	zr, err := gzip.NewReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("gunzip: %v", err)
	}
	defer zr.Close()
	out, err := io.ReadAll(zr)
	if err != nil {
		t.Fatalf("gunzip read: %v", err)
	}
	return out
}

// fixturePath returns an absolute path to the named fixture under
// tools/bcftools/testdata/parity. Tests that need other testdata
// subdirs construct their own paths.
func fixturePath(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	return abs
}

// assertEqualStdout runs two binaries with the same args on the same
// fixture and asserts stripProvenance-equal stdout. The case name is
// reported on failure so the matrix is readable.
func assertEqualStdout(t *testing.T, live, ours string, args ...string) {
	t.Helper()
	livOut := runBin(t, live, args...)
	ourOut := runBin(t, ours, args...)
	if !bytes.Equal(stripProvenanceBytes(livOut), stripProvenanceBytes(ourOut)) {
		t.Errorf("stdout diverges for args=%v\n--- live (%d bytes) ---\n%s\n--- ours (%d bytes) ---\n%s",
			args, len(livOut), snippet(livOut, 800),
			len(ourOut), snippet(ourOut, 800))
	}
}

// snippet truncates b to n bytes for log readability.
func snippet(b []byte, n int) string {
	if len(b) <= n {
		return string(b)
	}
	return string(b[:n]) + "\n... [truncated]"
}

// -------------------------------------------------------------------------
// SORT
// -------------------------------------------------------------------------

func TestLiveSort(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	assertEqualStdout(t, live, ours, "sort", fx)
}

func copyFileT(t *testing.T, src, dst string) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("copy read: %v", err)
	}
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatalf("copy write: %v", err)
	}
}

// -------------------------------------------------------------------------
// MPILEUP — smoke-test against one vendored fixture. The full matrix
// lives in mpileup_golden_test.go.
// -------------------------------------------------------------------------

func TestLiveMpileupSmoke(t *testing.T) {
	live, ours := requireLive(t)
	mpileupDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "mpileup"))
	if err != nil {
		t.Fatal(err)
	}
	bam := filepath.Join(mpileupDir, "mpileup.1.bam")
	ref := filepath.Join(mpileupDir, "mpileup.ref.fa")
	if _, err := os.Stat(bam); err != nil {
		t.Fatalf("vendored mpileup BAM fixture missing: %s: %v", bam, err)
	}
	if _, err := os.Stat(ref); err != nil {
		t.Fatalf("vendored mpileup reference fixture missing: %s: %v", ref, err)
	}
	// mpileup output is voluminous; we just verify both run and
	// agree on stdout for a minimal invocation. The mpileup_golden
	// suite covers the full matrix.
	assertEqualStdout(t, live, ours, "mpileup", "-f", ref, bam)
}

// TestLiveMpileupIndelsCNS exercises the salvaged `mpileup --indels-cns`
// path (the bam2bcf_indelcns.c / edlib consensus indel caller). It asserts
// byte-for-byte parity of the whole VCF stream (header + every record)
// against the upstream binary on the upstream indel-AD fixture, which is
// dense in indels and so stresses the consensus alignment scoring.
func TestLiveMpileupIndelsCNS(t *testing.T) {
	live, ours := requireLive(t)
	mpileupDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "mpileup"))
	if err != nil {
		t.Fatal(err)
	}
	bam := filepath.Join(mpileupDir, "indel-AD.1.bam")
	ref := filepath.Join(mpileupDir, "indel-AD.1.fa")
	if _, err := os.Stat(bam); err != nil {
		t.Fatalf("vendored indel-AD BAM fixture missing: %s: %v", bam, err)
	}
	if _, err := os.Stat(ref); err != nil {
		t.Fatalf("vendored indel-AD reference fixture missing: %s: %v", ref, err)
	}
	assertEqualStdout(t, live, ours, "mpileup", "--indels-cns", "-f", ref, bam)
}

// -------------------------------------------------------------------------
// assertRejectionParity invokes a subcommand on both binaries with
// args expected to fail (e.g. missing required flags or unknown
// plugin) and asserts that both binaries exit non-zero and produce no
// stdout. This is the canonical fallback when we can't get byte-equal
// output for a subcommand (e.g. fundamental architectural divergence
// or upstream features we have not yet ported), but we still want a
// live-oracle gate that locks in observable behaviour. Modelled on
// tools/mosdepth/pkg/mosdepth/live_oracle_test.go.
// -------------------------------------------------------------------------

func assertRejectionParity(t *testing.T, args []string) {
	t.Helper()
	live, ours := requireLive(t)

	run := func(bin string) (rejected bool, stdoutBytes int) {
		cmd := exec.Command(bin, args...)
		var out bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = io.Discard
		err := cmd.Run()
		return err != nil, out.Len()
	}

	liveRejected, liveStdout := run(live)
	if !liveRejected {
		t.Fatalf("upstream bcftools unexpectedly ACCEPTED %v; the oracle assumption is wrong", args)
	}
	if liveStdout > 0 {
		t.Fatalf("upstream bcftools wrote %d bytes to stdout for %v despite rejecting it",
			liveStdout, args)
	}
	oursRejected, oursStdout := run(ours)
	if !oursRejected {
		t.Errorf("our port accepted %v but upstream rejects it", args)
	}
	if oursStdout > 0 {
		t.Errorf("our port wrote %d bytes to stdout for %v but upstream produces none",
			oursStdout, args)
	}
}

// -------------------------------------------------------------------------
// CALL — the multiallelic caller (`call -m`) is a faithful port of
// mcall.c. We generate the mpileup VCF with the live binary, then assert
// our `call -m` output byte-matches upstream's over the whole contig
// (4000+ sites: EM allele-frequency estimation, per-site QUAL, the
// max-likelihood GT, and the INFO rewrite AN/AC/DP4/MQ). The same
// fixture is exercised with -v (variants only) and -A (keep alts).
// -------------------------------------------------------------------------

func TestLiveCall(t *testing.T) {
	live, ours := requireLive(t)
	mpDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "mpileup"))
	if err != nil {
		t.Fatal(err)
	}
	ref := filepath.Join(mpDir, "mpileup.ref.fa")
	bam := filepath.Join(mpDir, "mpileup.1.bam")
	if _, err := os.Stat(bam); err != nil {
		t.Fatalf("vendored mpileup fixture missing: %s: %v", bam, err)
	}

	// Produce the mpileup VCF with the genuine binary so the caller input
	// is upstream-identical; both `call` binaries then consume it.
	mp := runBin(t, live, "mpileup", "-f", ref, bam)
	mpFile := filepath.Join(t.TempDir(), "mp.vcf")
	if err := os.WriteFile(mpFile, mp, 0o644); err != nil {
		t.Fatal(err)
	}
	// `call --gvcf` rejects input lacking ##FORMAT=<ID=DP,...> in the
	// header (vcfcall.c:723). The default mpileup output does not
	// declare FORMAT/DP even though each record carries DP — request
	// it explicitly with -a DP for the --gvcf-only subtests below.
	mpDP := runBin(t, live, "mpileup", "-a", "DP", "-f", ref, bam)
	mpDPFile := filepath.Join(t.TempDir(), "mp_dp.vcf")
	if err := os.WriteFile(mpDPFile, mpDP, 0o644); err != nil {
		t.Fatal(err)
	}

	for _, tc := range []struct {
		name string
		args []string
	}{
		{"m", []string{"call", "-m", mpFile}},
		{"m_v", []string{"call", "-m", "-v", mpFile}},
		{"m_A", []string{"call", "-m", "-A", mpFile}},
		// --ploidy GRCh37 / 38: per-contig ploidy table. The mpileup
		// fixture contigs (1, 19, etc.) all fall under the "*" default
		// (F=2), so the output is identical to the diploid -m run —
		// what we're really asserting is that --ploidy still byte-
		// matches upstream and didn't regress the diploid path.
		{"m_ploidy_grch37", []string{"call", "-m", "--ploidy", "GRCh37", mpFile}},
		{"m_ploidy_grch38", []string{"call", "-m", "--ploidy", "GRCh38", mpFile}},
		// Consensus caller: with the ccall.c + em.c + prob1.c port
		// (callc.go) the `-c` pipeline now byte-matches upstream
		// (modulo provenance lines stripped by stripProvenance).
		{"c", []string{"call", "-c", mpFile}},
		// --gvcf banded blocks. Bands consecutive ref-only sites by
		// per-sample MIN_DP bin; the body of every gVCF block plus
		// every variant record byte-matches upstream
		// (`reference_code/bcftools/bcftools call -m --gvcf 0,5,10`).
		// See tools/bcftools/pkg/bcftools/callm_gvcf.go.
		{"m_gvcf_0_5_10", []string{"call", "-m", "--gvcf", "0,5,10", mpDPFile}},
		{"m_gvcf_5", []string{"call", "-m", "--gvcf", "5", mpDPFile}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			assertEqualStdout(t, live, ours, tc.args...)
		})
	}

	// -C alleles -T sites.tsv: constrain the multiallelic caller to a
	// user-supplied REF,ALT set per (CHROM,POS). The sites file is
	// derived from the mpileup output above by selecting the first
	// SNP allele at each record — guarantees the constraint set is a
	// subset of upstream's natural calls, so any divergence is the
	// projection itself (alsMap, PL re-index, INFO/QS reproject,
	// FORMAT/AD reproject) misbehaving rather than mcall.c EM drift.
	sitesTSV := filepath.Join(t.TempDir(), "sites.tsv")
	{
		var sb bytes.Buffer
		for _, line := range strings.Split(string(mp), "\n") {
			if line == "" || line[0] == '#' {
				continue
			}
			f := strings.SplitN(line, "\t", 6)
			if len(f) < 5 || len(f[3]) != 1 || strings.ContainsAny(f[3], "<>") {
				continue
			}
			alts := strings.Split(f[4], ",")
			alt := ""
			for _, a := range alts {
				if a == "<*>" || len(a) != 1 || strings.ContainsAny(a, "<>") {
					continue
				}
				alt = a
				break
			}
			if alt == "" {
				continue
			}
			sb.WriteString(f[0])
			sb.WriteByte('\t')
			sb.WriteString(f[1])
			sb.WriteByte('\t')
			sb.WriteString(f[3])
			sb.WriteByte(',')
			sb.WriteString(alt)
			sb.WriteByte('\n')
		}
		if err := os.WriteFile(sitesTSV, sb.Bytes(), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	t.Run("m_C_alleles", func(t *testing.T) {
		assertEqualStdout(t, live, ours,
			"call", "-m", "-C", "alleles", "-T", sitesTSV, mpFile)
	})

	// BCF input through call -m: upstream's vcfcall.c accepts BCF
	// natively, and now so does our port — the streamer routes BCF
	// records through bcf.NewReader and reuses the same caller loop
	// (see callStreaming in call.go). Convert the mpileup VCF to BCF
	// with upstream, then assert byte-equality.
	bcfFile := filepath.Join(t.TempDir(), "mp.bcf")
	runBin(t, live, "view", "-Ob", "-o", bcfFile, mpFile)
	t.Run("m_bcf_input", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", bcfFile)
	})

	// -r region post-filter: upstream uses the tabix index to seek; our
	// port scans the stream and filters in-process. The OUTPUT is
	// byte-exact (architectural-parity / perf-only residual). We assert
	// this against a bgzip-indexed copy of the mpileup fixture so
	// upstream can actually seek.
	bgzipBin := filepath.Join(filepath.Dir(filepath.Dir(live)), "htslib", "bgzip")
	if _, err := os.Stat(bgzipBin); err == nil {
		gz := filepath.Join(t.TempDir(), "mp.vcf.gz")
		// bgzip -k <plain> emits <plain>.gz; copy the plain VCF in first.
		if err := os.WriteFile(strings.TrimSuffix(gz, ".gz"), mp, 0o644); err != nil {
			t.Fatal(err)
		}
		runBin(t, bgzipBin, strings.TrimSuffix(gz, ".gz"))
		runBin(t, live, "index", "-t", gz)
		t.Run("m_r_region", func(t *testing.T) {
			assertEqualStdout(t, live, ours,
				"call", "-m", "-r", "17:103-110", gz)
		})
	}

	// -G sample groups (mcall.c::init_sample_groups). Builds a fresh
	// multi-sample mpileup with `-a AD` (mandatory for -G's per-group
	// AD-based qsum recomputation) using mpileup.1+mpileup.2, then
	// exercises the two non-degenerate -G shapes through both
	// binaries (the degenerate "single resolved group" file is
	// equivalent to the no-G default and is covered above):
	//
	//   * `-G -`: every sample is its own group (nsmpl_grp==nsmpl)
	//   * named two-group file: HG00100 → popA, HG00101 → popB
	bam2 := filepath.Join(mpDir, "mpileup.2.bam")
	if _, err := os.Stat(bam2); err == nil {
		mp2 := runBin(t, live, "mpileup", "-a", "AD", "-f", ref, bam, bam2)
		mpADFile := filepath.Join(t.TempDir(), "mp_ad2.vcf")
		if err := os.WriteFile(mpADFile, mp2, 0o644); err != nil {
			t.Fatal(err)
		}
		twoGroup := filepath.Join(t.TempDir(), "g2.tsv")
		if err := os.WriteFile(twoGroup, []byte("HG00100\tpopA\nHG00101\tpopB\n"), 0o644); err != nil {
			t.Fatal(err)
		}
		t.Run("m_G_dash_per_sample", func(t *testing.T) {
			assertEqualStdout(t, live, ours, "call", "-m", "-G", "-", mpADFile)
		})
		t.Run("m_G_two_groups", func(t *testing.T) {
			assertEqualStdout(t, live, ours, "call", "-m", "-G", twoGroup, mpADFile)
		})
	}

	// -a GQ: per-sample FORMAT/GQ on variant sites.
	t.Run("m_a_GQ", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "-a", "GQ", mpFile)
	})

	// -i/--insert-missed with -C alleles -T sites.tsv: synthetic
	// records emitted for sites in the sites file that mpileup
	// never produced. Build a sites TSV from the mpileup oracle
	// plus one site (17:5000) outside the input range so the
	// end-of-stream flush is exercised.
	missedSites := filepath.Join(t.TempDir(), "sites_miss.tsv")
	{
		var sb bytes.Buffer
		// Take the first concrete-SNP site so the present-site
		// portion of the output is non-trivial.
		for _, line := range strings.Split(string(mp), "\n") {
			if line == "" || line[0] == '#' {
				continue
			}
			f := strings.SplitN(line, "\t", 6)
			if len(f) < 5 || len(f[3]) != 1 || strings.ContainsAny(f[3], "<>") {
				continue
			}
			alts := strings.Split(f[4], ",")
			alt := ""
			for _, a := range alts {
				if a == "<*>" || len(a) != 1 || strings.ContainsAny(a, "<>") {
					continue
				}
				alt = a
				break
			}
			if alt == "" {
				continue
			}
			sb.WriteString(f[0])
			sb.WriteByte('\t')
			sb.WriteString(f[1])
			sb.WriteByte('\t')
			sb.WriteString(f[3])
			sb.WriteByte(',')
			sb.WriteString(alt)
			sb.WriteByte('\n')
			break
		}
		// Past-end site for the flush-on-EOF path (the fixture's
		// contig 17 has length 4200; pos 5000 is past every input
		// record).
		sb.WriteString("17\t5000\tA,T\n")
		if err := os.WriteFile(missedSites, sb.Bytes(), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	t.Run("m_C_alleles_insert_missed", func(t *testing.T) {
		assertEqualStdout(t, live, ours,
			"call", "-m", "-C", "alleles", "-T", missedSites, "-i", mpFile)
	})

	// -V indels: skip indel records; -V snps: skip SNP records.
	t.Run("m_V_indels", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "-V", "indels", mpFile)
	})
	t.Run("m_V_snps", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "-V", "snps", mpFile)
	})

	// -* / --keep-unseen-allele: retain the <*> ALT in the output.
	t.Run("m_keep_unseen", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "-*", mpFile)
	})

	// -M / --keep-masked-ref: by default mcall.c drops records whose
	// REF base is N; -M overrides. The mpileup fixture has no N
	// REFs so the test exercises the codepath without observable
	// drift, but the comparison anchors the implementation.
	t.Run("m_keep_masked_ref", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "-M", mpFile)
	})

	// -F AN,AC: incorporate prior allele frequencies from the input
	// INFO tags. The mpileup fixture lacks the panel-AF tags so
	// the formula is a no-op; the assertion anchors the wiring.
	t.Run("m_F_prior_freqs", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "-F", "AN,AC", mpFile)
	})

	// --ploidy-file FILE: per-region per-sex ploidy override. A
	// uniform diploid file produces the same output as the default
	// diploid run; pin the round-trip.
	ploidyFile := filepath.Join(t.TempDir(), "ploidy.txt")
	if err := os.WriteFile(ploidyFile, []byte("17 1 4200 M 2\n17 1 4200 F 2\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	t.Run("m_ploidy_file", func(t *testing.T) {
		assertEqualStdout(t, live, ours, "call", "-m", "--ploidy-file", ploidyFile, mpFile)
	})
}
