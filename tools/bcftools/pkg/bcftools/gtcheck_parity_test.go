package bcftools

// Live parity tests for bcftools gtcheck.
//
// Unlike the recorded-golden tests elsewhere in this package, these build
// the upstream `bcftools` binary from the reference_code/bcftools
// submodule and run it directly, comparing its `.tsv` output to the Go
// port byte-for-byte. There is intentionally NO t.Skip: if the submodule
// or its build environment is unavailable the test fails (t.Fatalf), so a
// broken parity loop is loud rather than silent.
//
// Two upstream output lines are NOT reproducible and are stripped before
// comparison on both sides:
//   - the "# This file was produced by bcftools ..." provenance block
//     (depends on argv + working directory);
//   - the "INFO\tTime required to process one record ..." timing line.
// Everything else (the INFO stats block, the DCv2 table and the DS block)
// is compared verbatim.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"testing"
)

// upstreamGtcheckOnce guards the one-time upstream build.
var upstreamGtcheckOnce sync.Once

// upstreamGtcheckPath is the built binary path (empty on failure).
var upstreamGtcheckPath string

// upstreamGtcheckErr records any build failure for reporting.
var upstreamGtcheckErr error

// upstreamBcftoolsGtcheck builds (once) and returns the path to the
// upstream bcftools binary from reference_code/bcftools.
func upstreamBcftoolsGtcheck(t *testing.T) string {
	t.Helper()
	upstreamGtcheckOnce.Do(func() {
		repo, err := repoRootGtcheck()
		if err != nil {
			upstreamGtcheckErr = err
			return
		}
		// Search for an already-built upstream binary across candidate
		// repo roots: the (possibly worktree) module root first, then the
		// canonical checkout that a git worktree shares its submodules
		// from (strip the .claude/worktrees/<name> suffix).
		for _, cand := range candidateRepoRoots(repo) {
			bin := filepath.Join(cand, "reference_code", "bcftools", "bcftools")
			if fi, statErr := os.Stat(bin); statErr == nil && !fi.IsDir() {
				upstreamGtcheckPath = bin
				return
			}
		}
		// No prebuilt binary found: build in this repo's submodule tree.
		bcfDir := filepath.Join(repo, "reference_code", "bcftools")
		bin := filepath.Join(bcfDir, "bcftools")
		// Build htslib then bcftools.
		htsDir := filepath.Join(repo, "reference_code", "htslib")
		for _, step := range []struct {
			dir  string
			name string
			args []string
		}{
			{htsDir, "autoreconf", []string{"-i"}},
			{htsDir, "sh", []string{"-c", "./configure"}},
			{htsDir, "make", []string{"-j4"}},
			{bcfDir, "sh", []string{"-c", "autoheader; autoconf; ./configure"}},
			{bcfDir, "make", []string{"-j4"}},
		} {
			cmd := exec.Command(step.name, step.args...)
			cmd.Dir = step.dir
			if out, runErr := cmd.CombinedOutput(); runErr != nil {
				// autoreconf/autoheader may be unnecessary if configure
				// already exists; only the make steps are fatal.
				if step.name == "make" {
					upstreamGtcheckErr = &buildError{dir: step.dir, out: out, err: runErr}
					return
				}
			}
		}
		if fi, statErr := os.Stat(bin); statErr == nil && !fi.IsDir() {
			upstreamGtcheckPath = bin
			return
		}
		upstreamGtcheckErr = &buildError{dir: bcfDir, out: nil, err: os.ErrNotExist}
	})
	if upstreamGtcheckErr != nil {
		t.Fatalf("upstream bcftools build failed: %v", upstreamGtcheckErr)
	}
	if upstreamGtcheckPath == "" {
		t.Fatalf("upstream bcftools binary not found after build")
	}
	return upstreamGtcheckPath
}

// candidateRepoRoots returns repo roots to probe for a prebuilt upstream
// binary. Git worktrees live under <main>/.claude/worktrees/<name> and
// share the main checkout's submodules, so the main checkout is included.
func candidateRepoRoots(repo string) []string {
	roots := []string{repo}
	const marker = "/.claude/worktrees/"
	if i := strings.Index(repo, marker); i >= 0 {
		roots = append(roots, repo[:i])
	}
	return roots
}

// repoRootGtcheck walks up from the working directory until it finds the
// module's go.mod.
func repoRootGtcheck() (string, error) {
	dir, err := os.Getwd()
	if err != nil {
		return "", err
	}
	for {
		if _, statErr := os.Stat(filepath.Join(dir, "go.mod")); statErr == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", os.ErrNotExist
		}
		dir = parent
	}
}

// stripNonReproducible removes the provenance header block and the timing
// line from gtcheck output so the rest can be compared byte-for-byte.
func stripNonReproducible(s string) string {
	var b strings.Builder
	for _, line := range strings.Split(s, "\n") {
		switch {
		case strings.HasPrefix(line, "# This file was produced by"):
			continue
		case strings.HasPrefix(line, "# \t"):
			continue
		case strings.HasPrefix(line, "# and the working directory"):
			continue
		case line == "#":
			continue
		case strings.HasPrefix(line, "INFO\tTime required"):
			continue
		}
		b.WriteString(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// runUpstreamGtcheck runs the upstream binary with args on stdin data.
func runUpstreamGtcheck(t *testing.T, bin string, args []string, dir string) string {
	t.Helper()
	cmd := exec.Command(bin, append([]string{"gtcheck"}, args...)...)
	cmd.Dir = dir
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream gtcheck %v failed: %v\nstderr:\n%s", args, err, errBuf.String())
	}
	return out.String()
}

// runGoGtcheck runs the Go port with the matching options.
func runGoGtcheck(t *testing.T, queryPath string, opts GtcheckOptions) string {
	t.Helper()
	var out bytes.Buffer
	if _, err := GtcheckFile(queryPath, &out, opts); err != nil {
		t.Fatalf("Go GtcheckFile failed: %v", err)
	}
	return out.String()
}

// parityFixtureVCF is a 6-sample, multi-site biallelic VCF carrying both
// GT and PL plus AC,AN, exercising every scoring path.
const parityFixtureVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	1000	.	A	T	.	.	AC=4;AN=12	GT:PL	0/0:0,30,255	0/0:0,40,200	1/1:255,30,0	0/1:30,0,30	1/1:200,20,0	0/0:0,50,255
chr1	2000	.	C	G	.	.	AC=6;AN=12	GT:PL	0/1:30,0,30	1/1:255,30,0	0/1:25,0,25	1/1:180,20,0	0/0:0,30,250	1/1:255,40,0
chr1	3000	.	G	A	.	.	AC=3;AN=12	GT:PL	1/1:255,30,0	1/1:255,30,0	0/0:0,30,255	0/0:0,40,200	0/1:30,0,30	0/0:0,50,255
chr1	4000	.	T	C	.	.	AC=5;AN=12	GT:PL	0/0:0,30,255	0/1:30,0,30	1/1:255,30,0	1/1:200,20,0	0/0:0,30,250	1/1:255,40,0
chr1	5000	.	A	G	.	.	AC=2;AN=12	GT:PL	0/0:0,30,255	0/0:0,40,200	0/1:30,0,30	0/0:0,40,200	1/1:200,20,0	0/0:0,50,255
chr1	6000	.	C	T	.	.	AC=7;AN=12	GT:PL	1/1:255,30,0	0/1:30,0,30	1/1:255,30,0	0/1:25,0,25	1/1:200,20,0	0/0:0,50,255
chr1	7000	.	G	C	.	.	AC=4;AN=12	GT:PL	0/1:30,0,30	0/0:0,40,200	1/1:255,30,0	0/0:0,40,200	0/1:30,0,30	1/1:255,40,0
chr1	8000	.	T	A	.	.	AC=6;AN=12	GT:PL	1/1:255,30,0	1/1:255,30,0	0/1:25,0,25	1/1:200,20,0	0/0:0,30,250	0/1:30,0,30
`

// writeParityFixture writes the fixture into a temp dir, returning the
// dir and the file path.
func writeParityFixture(t *testing.T) (dir, path string) {
	t.Helper()
	dir = t.TempDir()
	path = filepath.Join(dir, "query.vcf")
	if err := os.WriteFile(path, []byte(parityFixtureVCF), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	return dir, path
}

// TestParityGtcheck_AllModes compares the Go port against the live
// upstream binary across every scoring and output mode on a shared
// multi-sample fixture, asserting byte-for-byte equality of the
// non-provenance output.
func TestParityGtcheck_AllModes(t *testing.T) {
	bin := upstreamBcftoolsGtcheck(t)
	dir, path := writeParityFixture(t)

	cases := []struct {
		name string
		args []string
		opts GtcheckOptions
	}{
		{"default-auto", nil, GtcheckOptions{}},
		{"use-GT", []string{"-u", "GT"}, GtcheckOptions{UseTag: "GT"}},
		{"use-PL", []string{"-u", "PL"}, GtcheckOptions{UseTag: "PL"}},
		{"use-GT-PL", []string{"-u", "GT,PL"}, GtcheckOptions{UseTag: "GT,PL"}},
		{"use-PL-GT", []string{"-u", "PL,GT"}, GtcheckOptions{UseTag: "PL,GT"}},
		{"E0-integer", []string{"-u", "GT", "-E", "0"}, GtcheckOptions{UseTag: "GT", ErrorProbabilityZero: true}},
		{"E0-PL-integer", []string{"-u", "PL", "-E", "0"}, GtcheckOptions{UseTag: "PL", ErrorProbabilityZero: true}},
		{"E20", []string{"-u", "GT", "-E", "20"}, GtcheckOptions{UseTag: "GT", ErrorProbability: 20}},
		{"E60", []string{"-u", "GT", "-E", "60"}, GtcheckOptions{UseTag: "GT", ErrorProbability: 60}},
		{"no-HWE", []string{"-u", "GT", "--no-HWE-prob"}, GtcheckOptions{UseTag: "GT", NoHWEProb: true}},
		{"no-HWE-PL", []string{"-u", "PL", "--no-HWE-prob"}, GtcheckOptions{UseTag: "PL", NoHWEProb: true}},
		{"nmatches-2", []string{"--n-matches", "2"}, GtcheckOptions{NMatches: 2}},
		{"nmatches-neg2", []string{"--n-matches", "-2"}, GtcheckOptions{NMatches: -2}},
		{"nmatches-GT-3", []string{"-u", "GT", "--n-matches", "3"}, GtcheckOptions{UseTag: "GT", NMatches: 3}},
		{"pairs", []string{"-p", "S3,S1,S1,S2,S2,S6"}, GtcheckOptions{PairsSpec: "S3,S1,S1,S2,S2,S6"}},
		// -i/-e filter expressions (the AC INFO field varies per site in the
		// fixture, so these drop a subset of sites and change the scoring).
		{"include-AC", []string{"-i", "AC>3"}, GtcheckOptions{IncludeExpr: "AC>3"}},
		{"include-qry-AC", []string{"-i", "qry:AC>3"}, GtcheckOptions{IncludeExpr: "qry:AC>3"}},
		{"include-GT-AC", []string{"-u", "GT", "-i", "AC>=4"}, GtcheckOptions{UseTag: "GT", IncludeExpr: "AC>=4"}},
		// Use the qry:-prefixed -e form: upstream's bare `-e EXPR` has a quirk
		// where strtol(EXPR) clobbers the error-probability to 0 (integer
		// scoring); the prefixed form avoids it on both sides. See
		// docs/UPSTREAM_BUGS.md.
		{"exclude-qry-AC", []string{"-e", "qry:AC<4"}, GtcheckOptions{ExcludeExpr: "qry:AC<4"}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			up := stripNonReproducible(runUpstreamGtcheck(t, bin, append(tc.args, path), dir))
			tc.opts.OutputType = "t"
			got := stripNonReproducible(runGoGtcheck(t, path, tc.opts))
			if up != got {
				t.Fatalf("parity mismatch [%s]\n--- upstream ---\n%s\n--- go ---\n%s", tc.name, up, got)
			}
		})
	}
}

// TestParityGtcheck_DistinctiveSites compares the --distinctive-sites
// block (and the surrounding DCv2 table) against upstream for several NUM
// thresholds. The hts_lrand48 tie-break is reproduced, so the DS rows
// match byte-for-byte.
func TestParityGtcheck_DistinctiveSites(t *testing.T) {
	bin := upstreamBcftoolsGtcheck(t)
	dir, path := writeParityFixture(t)

	pairs := "S1,S2,S2,S3,S1,S3,S4,S5,S1,S6,S3,S6"
	for _, num := range []string{"1", "0.5", "2", "3", "0.25"} {
		t.Run("ds-"+num, func(t *testing.T) {
			args := []string{"-p", pairs, "--distinctive-sites", num, path}
			up := stripNonReproducible(runUpstreamGtcheck(t, bin, args, dir))

			dsNum := mustParseFloat(t, num)
			got := stripNonReproducible(runGoGtcheck(t, path, GtcheckOptions{
				PairsSpec:           pairs,
				DistinctiveSites:    dsNum,
				HasDistinctiveSites: true,
				OutputType:          "t",
			}))
			if up != got {
				t.Fatalf("distinctive parity mismatch [num=%s]\n--- upstream ---\n%s\n--- go ---\n%s", num, up, got)
			}
		})
	}
}

// TestParityGtcheck_Monoallelic checks the no-ALT monoallelic skip and
// its --keep-refs override against upstream.
func TestParityGtcheck_Monoallelic(t *testing.T) {
	bin := upstreamBcftoolsGtcheck(t)
	const mono = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
chr1	100	.	A	.	.	.	.	GT	0/0	0/0	0/0
chr1	200	.	C	G	.	.	.	GT	0/0	0/1	1/1
chr1	300	.	G	A	.	.	.	GT	1/1	0/0	0/1
`
	dir := t.TempDir()
	path := filepath.Join(dir, "mono.vcf")
	if err := os.WriteFile(path, []byte(mono), 0o644); err != nil {
		t.Fatalf("write mono fixture: %v", err)
	}

	for _, tc := range []struct {
		name string
		args []string
		opts GtcheckOptions
	}{
		{"skip", []string{"-u", "GT"}, GtcheckOptions{UseTag: "GT"}},
		{"keep-refs", []string{"-u", "GT", "--keep-refs"}, GtcheckOptions{UseTag: "GT", KeepRefs: true}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			up := stripNonReproducible(runUpstreamGtcheck(t, bin, append(tc.args, path), dir))
			tc.opts.OutputType = "t"
			got := stripNonReproducible(runGoGtcheck(t, path, tc.opts))
			if up != got {
				t.Fatalf("monoallelic parity mismatch [%s]\n--- upstream ---\n%s\n--- go ---\n%s", tc.name, up, got)
			}
		})
	}
}

func mustParseFloat(t *testing.T, s string) float64 {
	t.Helper()
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		t.Fatalf("parse %q: %v", s, err)
	}
	return f
}
