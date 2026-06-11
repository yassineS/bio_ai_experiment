package samtools

import (
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
	"time"
)

// repoRootForTest returns the absolute path to the repository root,
// derived from this test file's own location (tools/samtools/pkg/samtools
// is four directories below the root). Deriving it from the source path
// rather than the working directory keeps the helper robust regardless of
// where `go test` is invoked from.
func repoRootForTest(t *testing.T) string {
	t.Helper()
	_, thisFile, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed; cannot locate repo root")
	}
	// thisFile == <root>/tools/samtools/pkg/samtools/consensus_upstream_test.go
	root := filepath.Clean(filepath.Join(filepath.Dir(thisFile), "..", "..", "..", ".."))
	if _, err := os.Stat(filepath.Join(root, "go.mod")); err != nil {
		t.Fatalf("repo root %q has no go.mod: %v", root, err)
	}
	return root
}

// upstreamSamtoolsOnce memoises the (one-time, expensive) submodule
// init + htslib/samtools build across all tests in the package so that
// multiple upstream-parity tests share a single built binary.
var upstreamSamtoolsOnce struct {
	sync.Once
	path string
	err  error
}

// upstreamSamtoolsBinary locates, or builds from the vendored
// submodules, the upstream samtools binary and returns its absolute
// path. The build is performed at most once per test binary. It never
// calls t.Skip: per the project's testing rules the upstream parity
// check must actually execute, so a genuine inability to produce the
// binary is a hard failure (t.Fatalf), not a skip.
//
// Build steps mirror the docs and the CI upstream-parity job:
//  1. `git submodule update --init --recursive reference_code/htslib
//     reference_code/samtools` (only if the sources are missing), with
//     exponential backoff retries for transient network failures. The
//     --recursive flag is required because htslib carries its own nested
//     htscodecs submodule, without which its ./configure aborts.
//  2. `autoreconf -i && ./configure && make -j` in reference_code/htslib.
//  3. `make -j` in reference_code/samtools (its own ./configure is
//     skipped: it is optional and can clobber htslib's config.mk).
//
// An already-built reference_code/samtools/samtools is reused as-is.
func upstreamSamtoolsBinary(t *testing.T) string {
	t.Helper()
	root := repoRootForTest(t)
	upstreamSamtoolsOnce.Do(func() {
		upstreamSamtoolsOnce.path, upstreamSamtoolsOnce.err = buildUpstreamSamtools(root)
	})
	if upstreamSamtoolsOnce.err != nil {
		t.Fatalf("could not obtain upstream samtools binary: %v", upstreamSamtoolsOnce.err)
	}
	return upstreamSamtoolsOnce.path
}

// buildUpstreamSamtools performs the locate-or-build work for
// upstreamSamtoolsBinary, returning the binary path or an error.
func buildUpstreamSamtools(root string) (string, error) {
	samtoolsDir := filepath.Join(root, "reference_code", "samtools")
	htslibDir := filepath.Join(root, "reference_code", "htslib")
	bin := filepath.Join(samtoolsDir, "samtools")

	// Fast path: reuse an already-built binary.
	if fi, err := os.Stat(bin); err == nil && !fi.IsDir() {
		return bin, nil
	}

	// Ensure the submodule sources are present. --recursive pulls
	// htslib's nested htscodecs submodule, which its ./configure needs.
	htscodecsC := filepath.Join(htslibDir, "htscodecs", "htscodecs", "htscodecs.c")
	if _, err := os.Stat(filepath.Join(samtoolsDir, "bam_consensus.c")); err != nil {
		if err := runWithRetry(root, "git", "submodule", "update", "--init", "--recursive",
			"reference_code/htslib", "reference_code/samtools"); err != nil {
			return "", err
		}
	} else if _, err := os.Stat(htscodecsC); err != nil {
		// Sources present but the nested htscodecs submodule is missing
		// (e.g. a non-recursive prior init): pull it now.
		if err := runWithRetry(htslibDir, "git", "submodule", "update", "--init", "--recursive"); err != nil {
			return "", err
		}
	}

	// Build htslib: autoreconf -i (a git checkout ships no ./configure)
	// then ./configure && make. We always re-run ./configure so a stale
	// or error-stamped config.mk (samtools' own configure can overwrite
	// htslib's during bundled-htslib detection) is regenerated cleanly.
	if _, err := os.Stat(filepath.Join(htslibDir, "configure")); err != nil {
		if err := runCmd(htslibDir, "autoreconf", "-i"); err != nil {
			return "", err
		}
	}
	if err := runCmd(htslibDir, "./configure"); err != nil {
		return "", err
	}
	if err := runCmd(htslibDir, "make", "-j"); err != nil {
		return "", err
	}

	// Build samtools with a plain `make`. Per samtools' INSTALL, running
	// its own ./configure is optional and only diagnoses build problems;
	// the Makefile builds against the bundled ../htslib we just built.
	// We deliberately skip samtools' ./configure because it regenerates
	// (and can corrupt) htslib's config.mk.
	if err := runCmd(samtoolsDir, "make", "-j"); err != nil {
		return "", err
	}

	if fi, err := os.Stat(bin); err != nil || fi.IsDir() {
		return "", err
	}
	return bin, nil
}

// runCmd runs name+args in dir, returning a descriptive error (with
// captured combined output) on failure.
func runCmd(dir, name string, args ...string) error {
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	out, err := cmd.CombinedOutput()
	if err != nil {
		return &cmdError{name: name, args: args, dir: dir, out: out, err: err}
	}
	return nil
}

// runWithRetry runs name+args in dir, retrying up to four times with
// exponential backoff (2/4/8/16s) to absorb transient network failures
// (e.g. a flaky `git submodule` clone).
func runWithRetry(dir, name string, args ...string) error {
	backoff := []time.Duration{2 * time.Second, 4 * time.Second, 8 * time.Second, 16 * time.Second}
	var lastErr error
	for attempt := 0; attempt <= len(backoff); attempt++ {
		if attempt > 0 {
			time.Sleep(backoff[attempt-1])
		}
		if err := runCmd(dir, name, args...); err != nil {
			lastErr = err
			continue
		}
		return nil
	}
	return lastErr
}

// runUpstreamConsensus invokes the upstream samtools `consensus`
// subcommand on samPath and returns its stdout. samtools writes
// progress/warnings to stderr, which is discarded.
func runUpstreamConsensus(t *testing.T, bin, samPath string, args ...string) string {
	t.Helper()
	full := append([]string{"consensus"}, args...)
	full = append(full, samPath)
	cmd := exec.Command(bin, full...)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream samtools %v failed: %v", full, err)
	}
	return string(out)
}

// indelConsensusSAM is an indel-rich fixture used by the live indel-parity
// test. It mixes:
//   - a 3bp insertion (5M3I5M) shared by three reads, one with a SNP inside
//     the inserted run so the insertion columns exercise both homozygous and
//     heterozygous calls;
//   - a 2bp deletion (3M2D5M) on a strand-mixed pair;
//   - a read (a6, 8M) that terminates exactly at the reference base preceding
//     a downstream insertion — it must NOT pad that insertion column, which
//     is the membership rule ported from consensus_pileup.c::get_next_base;
//   - a 1bp insertion (4M1I4M) carrying a heterozygous inserted base.
//
// Coordinates are deliberately chosen so the insertion-column depth differs
// between the engine's "spanning" set and a naive "all reads in the base
// column" set, which is exactly what regressed before the fix.
const indelConsensusSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:20
a1	0	chr1	1	60	5M3I5M	*	0	0	ACGTAGGGCATGC	IIIIIIIIIIIII
a2	16	chr1	1	60	5M3I5M	*	0	0	ACGTAGGGCATGC	IIIIIIIIIIIII
a3	0	chr1	1	60	5M3I5M	*	0	0	ACGTAGTGCATGC	IIIIIIIIIIIII
a7	0	chr1	1	60	10M	*	0	0	ACGTAGCATG	II#III+III
a4	0	chr1	3	60	3M2D5M	*	0	0	GTACATGC	IIIIIIII
a5	16	chr1	3	60	3M2D5M	*	0	0	GTACATGC	IIIIIIII
a6	0	chr1	8	60	8M	*	0	0	CATGCTTA	IIIIIIII
b1	0	chr1	12	60	4M1I4M	*	0	0	TTAAGCGTA	HHHHHHHHH
b2	0	chr1	12	60	4M1I4M	*	0	0	TTAATCGTA	HHHHHHHHH
`

// upstreamSamtoolsConsensusIndel is the LIVE indel-parity check. It builds
// (or reuses) the vendored upstream samtools binary and compares its
// `consensus` output byte-for-byte against the Go port over an indel-rich
// fixture, sweeping the simple and bayesian modes, the FASTA/FASTQ/pileup
// formats, and the indel-relevant flags (--mark-ins, --show-del, -A, -a).
// Per the project's testing rules it never calls t.Skip on a build failure:
// an inability to produce the upstream binary is a hard failure.
func upstreamSamtoolsConsensusIndel(t *testing.T) {
	t.Helper()
	bin := upstreamSamtoolsBinary(t)

	samPath := filepath.Join(t.TempDir(), "indel.sam")
	if err := os.WriteFile(samPath, []byte(indelConsensusSAM), 0o600); err != nil {
		t.Fatalf("write indel fixture: %v", err)
	}

	modes := []struct {
		cli string
		mod ConsensusMode
	}{
		{"simple", ConsensusModeSimple},
		{"bayesian", ConsensusModeBayesian},
	}
	formats := []struct {
		cli string
		fmt ConsensusFormat
	}{
		{"fasta", ConsensusFASTA},
		{"fastq", ConsensusFASTQ},
		{"pileup", ConsensusPileup},
	}
	// Each variant pairs the extra upstream CLI flags with the matching Go
	// ConsensusOptions mutation, so the two invocations are equivalent.
	variants := []struct {
		name    string
		cliArgs []string
		apply   func(*ConsensusOptions)
	}{
		{"plain", nil, func(*ConsensusOptions) {}},
		{"mark-ins", []string{"--mark-ins"}, func(o *ConsensusOptions) { o.MarkIns = true }},
		{"show-del", []string{"--show-del", "yes"}, func(o *ConsensusOptions) { o.ShowDel = true }},
		{"ambig", []string{"-A"}, func(o *ConsensusOptions) { o.AmbigCodes = true }},
		{"all-pos", []string{"-a"}, func(o *ConsensusOptions) { o.AllPositions = true }},
		{"no-show-ins", []string{"--show-ins", "no"}, func(o *ConsensusOptions) { o.NoShowIns = true }},
		// --use-qual and --min-BQ route the insertion-column pads through
		// the simple gap bucket; these variants pin the gap-quality
		// weighting (score[16] += 8*q) and the pad min-qual gate so the
		// inserted-vs-pad call fraction matches upstream.
		{"use-qual", []string{"-q"}, func(o *ConsensusOptions) { o.UseQual = true }},
		{"min-bq", []string{"--min-BQ", "20"}, func(o *ConsensusOptions) { o.MinBaseQ = 20 }},
	}

	for _, m := range modes {
		for _, f := range formats {
			for _, v := range variants {
				name := m.cli + "_" + f.cli + "_" + v.name
				t.Run(name, func(t *testing.T) {
					args := append([]string{"-m", m.cli, "-f", f.cli}, v.cliArgs...)
					up := runUpstreamConsensus(t, bin, samPath, args...)

					opts := ConsensusOptions{Mode: m.mod, Format: f.fmt}
					v.apply(&opts)
					got := runConsensusOnSAM(t, indelConsensusSAM, opts)

					if got != up {
						t.Fatalf("indel parity mismatch (%s):\n--- upstream ---\n%s\n--- go ---\n%s",
							name, up, got)
					}
				})
			}
		}
	}
}

// TestConsensus_IndelUpstreamParity runs the live indel-parity sweep. It is
// skipped only in -short mode (which omits the expensive upstream build).
func TestConsensus_IndelUpstreamParity(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	upstreamSamtoolsConsensusIndel(t)
}

// allPosConsensusSAM is a fixture exercising the -a/--all-positions pileup
// placeholder-row emission. It is engineered to contain every flavour of
// position the placeholder mechanism must reproduce:
//
//   - a LEADING zero-coverage gap (positions 1-2 are uncovered);
//   - a covered run (positions 3-7);
//   - a DELETION-only run inside the reads (the 3bp deletion at positions
//     8-10 from `5M3D5M` at pos 3), which upstream emits as placeholder
//     `N 0 * *` rows AND duplicates at the leading deletion positions
//     (last_pos is not advanced for a suppressed '*' column);
//   - an INTERNAL zero-coverage gap between two covered blocks on chr1;
//   - a TRAILING zero-coverage gap (the contig is longer than any read);
//   - a wholly UNCOVERED contig (chr2) that only -aa fills.
const allPosConsensusSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:28
@SQ	SN:chr2	LN:6
d1	0	chr1	3	60	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII
d2	0	chr1	3	60	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII
d3	16	chr1	3	60	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII
e1	0	chr1	20	60	4M	*	0	0	TTGG	IIII
e2	0	chr1	20	60	4M	*	0	0	TTGG	IIII
`

// refSkipConsensusSAM exercises the -a pileup ref-skip (CIGAR N / intron)
// edge case. Three spliced reads carry a 10bp N intron (`5M10N5M` and a
// strand-mixed `4M10N5M`), so chr1:6-15 is covered ONLY by ref-skip events.
// Upstream keeps those reads in the column, counts them in depth, and renders
// each as '.', so an intron row is `depth 3 / N / ... / !!!` — NOT a
// zero-coverage placeholder. The exon bases bordering the intron (positions 5
// and 16) additionally carry upstream's p->ref_skip flag, which only the
// bayesian/Gap5 caller honours (excluding them from its depth), so the
// simple- and bayesian-mode boundaries differ; this fixture pins both.
const refSkipConsensusSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:24
s1	0	chr1	1	60	5M10N5M	*	0	0	ACGTAGGGCA	IIIIIIIIII
s2	0	chr1	1	60	5M10N5M	*	0	0	ACGTAGGGCA	IIIIIIIIII
s3	16	chr1	2	60	4M10N5M	*	0	0	CGTAGGGCA	IIIIIIIII
`

// delInsRunConsensusSAM exercises the interaction of a suppressed deletion
// RUN with an insertion at its trailing edge under -a. Three reads carry a
// 5bp deletion (`3M5D3M` at pos 8 -> del at chr1:11-15); a fourth read
// (`8M2I3M`) spans the whole run with real bases AND inserts 2bp immediately
// after position 15. With --show-del off, the nth==0 deletion rows for
// 11-15 are all suppressed (3 dels outvote 1 base), so c->last_pos never
// advances through the run. At position 15 the insertion columns (nth>0) are
// also pad-suppressed, and upstream RE-RUNS its gap fill at every nth
// (bam_consensus.c:2227), emitting the leading `11..14` placeholder block
// once per nth before last_pos finally advances. This is the per-nth gap-fill
// duplication and the insertion-row last_pos advance the Go port must mirror.
const delInsRunConsensusSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:24
g1	0	chr1	8	60	3M5D3M	*	0	0	AAACCC	IIIIII
g2	0	chr1	8	60	3M5D3M	*	0	0	AAACCC	IIIIII
g3	0	chr1	8	60	3M5D3M	*	0	0	AAACCC	IIIIII
g4	0	chr1	8	60	8M2I3M	*	0	0	AAAGGGGGTTCCC	IIIIIIIIIIIII
`

// upstreamSamtoolsConsensusAll is the LIVE parity check for the
// -a/--all-positions pileup placeholder rows. It builds (or reuses) the
// vendored upstream samtools binary and compares its `consensus -a/-aa
// --format pileup` output byte-for-byte against the Go port over three
// fixtures:
//   - allPosConsensusSAM: leading/internal/trailing zero-coverage gaps, a
//     deletion-only run (upstream's duplicate placeholder rows), and an
//     entirely uncovered contig (-aa);
//   - refSkipConsensusSAM: a spliced read with a long CIGAR N intron, so
//     intron positions are ref-skip columns (depth>0, '.' bases), not gaps;
//   - delInsRunConsensusSAM: a suppressed deletion run whose trailing edge
//     carries an insertion, exercising the per-nth gap-fill duplication and
//     the insertion-row last_pos advance.
//
// It sweeps the simple and bayesian modes and the --show-del on/off setting.
// Per the project's testing rules it never calls t.Skip on a build failure:
// an inability to produce the upstream binary is a hard failure.
func upstreamSamtoolsConsensusAll(t *testing.T) {
	t.Helper()
	bin := upstreamSamtoolsBinary(t)

	fixtures := []struct {
		name string
		sam  string
	}{
		{"gaps", allPosConsensusSAM},
		{"refskip", refSkipConsensusSAM},
		{"delins-run", delInsRunConsensusSAM},
	}

	modes := []struct {
		cli string
		mod ConsensusMode
	}{
		{"simple", ConsensusModeSimple},
		{"bayesian", ConsensusModeBayesian},
	}
	variants := []struct {
		name    string
		cliArgs []string
		apply   func(*ConsensusOptions)
	}{
		{"a", []string{"-a"}, func(o *ConsensusOptions) { o.AllPositions = true }},
		{"aa", []string{"-aa"}, func(o *ConsensusOptions) {
			o.AllPositions = true
			o.AllContigs = true
		}},
		{"a_show-del", []string{"-a", "--show-del", "yes"}, func(o *ConsensusOptions) {
			o.AllPositions = true
			o.ShowDel = true
		}},
		{"aa_show-del", []string{"-aa", "--show-del", "yes"}, func(o *ConsensusOptions) {
			o.AllPositions = true
			o.AllContigs = true
			o.ShowDel = true
		}},
	}

	for _, fx := range fixtures {
		samPath := filepath.Join(t.TempDir(), fx.name+".sam")
		if err := os.WriteFile(samPath, []byte(fx.sam), 0o600); err != nil {
			t.Fatalf("write %s fixture: %v", fx.name, err)
		}
		for _, m := range modes {
			for _, v := range variants {
				name := fx.name + "_" + m.cli + "_" + v.name
				t.Run(name, func(t *testing.T) {
					args := append([]string{"-m", m.cli, "-f", "pileup"}, v.cliArgs...)
					up := runUpstreamConsensus(t, bin, samPath, args...)

					opts := ConsensusOptions{Mode: m.mod, Format: ConsensusPileup}
					v.apply(&opts)
					got := runConsensusOnSAM(t, fx.sam, opts)

					if got != up {
						t.Fatalf("all-positions pileup parity mismatch (%s):\n--- upstream ---\n%s\n--- go ---\n%s",
							name, up, got)
					}
				})
			}
		}
	}
}

// TestConsensus_AllPositionsUpstreamParity runs the live -a/--all-positions
// pileup placeholder-row parity sweep. It is skipped only in -short mode
// (which omits the expensive upstream build).
func TestConsensus_AllPositionsUpstreamParity(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	upstreamSamtoolsConsensusAll(t)
}

// cmdError is a build/command failure carrying the captured output so
// the test log shows exactly what went wrong.
type cmdError struct {
	name string
	args []string
	dir  string
	out  []byte
	err  error
}

func (e *cmdError) Error() string {
	return strings.TrimSpace(
		e.name + " " + strings.Join(e.args, " ") + " (in " + e.dir + ") failed: " +
			e.err.Error() + "\n--- output ---\n" + string(e.out))
}

// upstreamBugSAM mixes heterozygous and homozygous positions: four reads
// over chr1:1-3 give a het flank (A/C at pos1 and pos3) and a homozygous
// middle (all-G at pos2). Because it contains BOTH kinds of position, any
// genuine --het-only filter (like ours) visibly changes the output by
// dropping the homozygous column — which is exactly what lets the test
// below distinguish "filter applied" from "flag ignored".
const upstreamBugSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:3
r1	0	chr1	1	60	3M	*	0	0	AGA	III
r2	0	chr1	1	60	3M	*	0	0	AGA	III
r3	0	chr1	1	60	3M	*	0	0	CGC	III
r4	0	chr1	1	60	3M	*	0	0	CGC	III
`

// TestConsensus_HetOnlyUpstreamBug is the LIVE check that documents the
// upstream dead-option bug and our intentional, correct divergence from
// it. It builds (or reuses) the vendored upstream samtools binary and
// proves two facts by execution:
//
//  1. UPSTREAM IGNORES --het-only (the bug). samtools parses --het-only
//     into opts.het_only (bam_consensus.c) but never reads that variable
//     anywhere in the calling/output path, so its consensus output is
//     byte-for-byte IDENTICAL with and without the flag — even though the
//     fixture contains a homozygous position the flag's name says should
//     be suppressed. Documented in docs/UPSTREAM_BUGS.md.
//  2. OUR Go --het-only DIFFERS from our own output without the flag: it
//     suppresses the homozygous column (rendered 'N' in FASTA, omitted in
//     pileup), which is the behaviour the flag implies. This is the
//     intentional divergence — we fix upstream's dead option.
//
// No committed golden/snapshot file is used: upstream and Go outputs are
// produced live and compared directly. If the upstream binary cannot be
// built, the test fails hard (it never skips), per the project rules.
func TestConsensus_HetOnlyUpstreamBug(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	bin := upstreamSamtoolsBinary(t)

	// Write the fixture to a tmp SAM file for the upstream binary.
	samPath := filepath.Join(t.TempDir(), "bug.sam")
	if err := os.WriteFile(samPath, []byte(upstreamBugSAM), 0o600); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	cases := []struct {
		name      string
		cliFormat string
		goFormat  ConsensusFormat
	}{
		{"fasta", "fasta", ConsensusFASTA},
		{"pileup", "pileup", ConsensusPileup},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Fact 1: upstream --het-only is a no-op (the bug). Output
			// is identical with and without the flag, despite the
			// fixture's homozygous pos2 that a real filter would drop.
			// We run with -A/--ambig so heterozygous positions are
			// actually emitted as IUPAC codes; that makes the homozygous
			// vs het distinction visible and guarantees a genuine filter
			// would change the output.
			upNoFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-A", "-f", c.cliFormat)
			upWithFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-A", "-f", c.cliFormat, "--het-only")
			if upNoFlag != upWithFlag {
				t.Fatalf("upstream --het-only changed %s output, expected the "+
					"documented no-op (bug):\n--- without ---\n%s\n--- with ---\n%s",
					c.name, upNoFlag, upWithFlag)
			}

			// Fact 2: OUR --het-only differs from our own no-flag output
			// — we suppress the homozygous column. This is the
			// intentional divergence that fixes the upstream dead option.
			goNoFlag := runConsensusOnSAM(t, upstreamBugSAM,
				ConsensusOptions{Format: c.goFormat, AmbigCodes: true})
			goWithFlag := runConsensusOnSAM(t, upstreamBugSAM,
				ConsensusOptions{Format: c.goFormat, AmbigCodes: true, HetOnly: true})
			if goNoFlag == goWithFlag {
				t.Fatalf("Go --het-only %s output did not change (expected the "+
					"homozygous column to be suppressed):\n--- without ---\n%s",
					c.name, goNoFlag)
			}
			// And our no-flag output must still match upstream's no-flag
			// output (we only diverge WHEN the flag is set), confirming we
			// are not accidentally changing the baseline consensus.
			if goNoFlag != upNoFlag {
				t.Fatalf("Go baseline (no --het-only) %s output differs from "+
					"upstream:\n--- go ---\n%s\n--- upstream ---\n%s",
					c.name, goNoFlag, upNoFlag)
			}
		})
	}
}
