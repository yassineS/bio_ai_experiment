package fastp

// Validated-parity tests for the fastp Go port against the upstream C++
// reference implementation (OpenGene/fastp 1.0.1).
//
// Process:
//   - The upstream binary is built from reference_code/fastp (a git
//     submodule) via `cd reference_code/fastp && make` — see
//     tools/PARITY_VALIDATION.md for the build instructions and the
//     library dependencies (libisal, libdeflate, libpthread).
//   - The static FASTQ fixtures live under tools/fastp/testdata/parity/
//     and are deterministic (regenerable via testdata/parity/generate.py).
//   - Each test runs the upstream binary AND the Go port on the same
//     fixture with the same flags, then compares either the output FASTQ
//     bytes (for trimming/filter tests) or the JSON counter subset (for
//     tests where the FASTQ output is identical and we want to also
//     verify the report).
//
// All tests are skipped automatically when the upstream binary is not
// available, so `go test ./tools/fastp/...` still passes on systems that
// have not built the C++ reference (e.g. a fresh CI worker that has not
// run `git submodule update --init reference_code/fastp` yet).
//
// Bugs surfaced by this audit and fixed inline (in this PR):
//   - UMI tag format: Go was always emitting ":UMI_<umi>" regardless of
//     whether --umi_prefix was set. Upstream emits ":<umi>" when no
//     prefix is set and ":<prefix>_<umi>" otherwise. Fixed in
//     tools/fastp/pkg/fastp/umi.go.
//   - Low-complexity formula: Go used "unique 2-mers / total 2-mers"
//     which classifies AT-repeats as low-complexity. Upstream defines
//     complexity as "fraction of adjacent positions where seq[i] !=
//     seq[i+1]". An AT-repeat is now (correctly) considered HIGH
//     complexity. Fixed in tools/fastp/pkg/fastp/fastp.go.
//   - low_complexity_reads counter was missing from the JSON report.
//     Added.
//
// Documented divergences (test t.Skip + entry in
// docs/UPSTREAM_BUGS.md / tools/PARITY_VALIDATION.md):
//   - poly-G trimming: upstream allows 1 mismatch per 8 bases scanned
//     (max 5 total). Our Go port does a strict consecutive-G count.
//   - sliding-window cut_right/cut_tail/cut_front: upstream's algorithm
//     keeps the leading high-quality bases of the offending window and
//     advances past 'N's; we cut the whole window. Off-by-1..2 bp.
//   - adapter auto-detect (SE & PE): different kmer/overlap heuristics
//     and a different default sample size.

import (
	"bytes"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// ensureUpstream returns the upstream fastp binary, building it from the
// reference_code/fastp submodule on first use. Per the env-guard policy
// (PR #294) it t.Fatalf's with an exact init/build hint when the binary
// cannot be located or built, rather than silently skipping.
func ensureUpstream(t *testing.T) string {
	t.Helper()
	p, err := upstreamFastp(t)
	if err != nil {
		t.Fatalf("upstream fastp binary unavailable; run `git submodule update --init reference_code/fastp && make -C reference_code/fastp`: %v", err)
	}
	return p
}

// parityInput returns the absolute path of the fixture named name under
// tools/fastp/testdata/parity/.
func parityInput(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("resolve parity fixture %q: %v", name, err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Fatalf("parity fixture %q missing: %v", name, err)
	}
	return abs
}

// runUpstreamSE invokes the upstream binary with `-i <in> -o <out>
// --json <json> -h <html>` plus extraFlags, in dir. Returns the
// generated output FASTQ bytes and the parsed JSON report.
func runUpstreamSE(t *testing.T, bin, in, dir string, extraFlags []string) ([]byte, map[string]any) {
	t.Helper()
	out := filepath.Join(dir, "up_out.fq")
	jsonPath := filepath.Join(dir, "up_report.json")
	htmlPath := filepath.Join(dir, "up_report.html")
	args := []string{"-i", in, "-o", out, "--json", jsonPath, "-h", htmlPath}
	args = append(args, extraFlags...)
	cmd := exec.Command(bin, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream fastp %v failed: %v", args, err)
	}
	return readFile(t, out), readJSON(t, jsonPath)
}

// runUpstreamPE invokes the upstream binary for paired-end with the
// `-i <r1> -I <r2> -o <or1> -O <or2>` skeleton + extraFlags. Returns
// (out R1, out R2, parsed JSON).
func runUpstreamPE(t *testing.T, bin, r1, r2, dir string, extraFlags []string) ([]byte, []byte, map[string]any) {
	t.Helper()
	or1 := filepath.Join(dir, "up_r1.fq")
	or2 := filepath.Join(dir, "up_r2.fq")
	jsonPath := filepath.Join(dir, "up_report.json")
	htmlPath := filepath.Join(dir, "up_report.html")
	args := []string{"-i", r1, "-I", r2, "-o", or1, "-O", or2, "--json", jsonPath, "-h", htmlPath}
	args = append(args, extraFlags...)
	cmd := exec.Command(bin, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream fastp PE %v failed: %v", args, err)
	}
	return readFile(t, or1), readFile(t, or2), readJSON(t, jsonPath)
}

// runGoFastpSE runs the Go binary in-process via the library API so we
// don't depend on the CLI being built. Uses ProcessSingleEnd to keep
// determinism (no threading, no time-of-day calls).
func runGoFastpSE(t *testing.T, in string, opts ProcessOptions) ([]byte, *ProcessStats) {
	t.Helper()
	f, err := os.Open(in)
	if err != nil {
		t.Fatalf("open %s: %v", in, err)
	}
	defer f.Close()
	var out bytes.Buffer
	stats, err := ProcessSingleEnd(f, &out, defaultEncoding(), opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	return out.Bytes(), stats
}

// runGoFastpPE drives ProcessPairedEnd in the same in-process manner.
func runGoFastpPE(t *testing.T, r1, r2 string, opts ProcessOptions) ([]byte, []byte, *ProcessStats) {
	t.Helper()
	in1, err := os.Open(r1)
	if err != nil {
		t.Fatalf("open %s: %v", r1, err)
	}
	defer in1.Close()
	in2, err := os.Open(r2)
	if err != nil {
		t.Fatalf("open %s: %v", r2, err)
	}
	defer in2.Close()
	var o1, o2 bytes.Buffer
	stats, err := ProcessPairedEnd(in1, in2, &o1, &o2, defaultEncoding(), opts)
	if err != nil {
		t.Fatalf("ProcessPairedEnd: %v", err)
	}
	return o1.Bytes(), o2.Bytes(), stats
}

func readFile(t *testing.T, p string) []byte {
	t.Helper()
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read %s: %v", p, err)
	}
	return b
}

func readJSON(t *testing.T, p string) map[string]any {
	t.Helper()
	b := readFile(t, p)
	var m map[string]any
	if err := json.Unmarshal(b, &m); err != nil {
		t.Fatalf("parse %s: %v", p, err)
	}
	return m
}

// mustEqualBytes fails the test if got != want with a small inline diff.
func mustEqualBytes(t *testing.T, label string, got, want []byte) {
	t.Helper()
	if bytes.Equal(got, want) {
		return
	}
	// Find the first differing line so the failure message stays small.
	gotLines := strings.Split(string(got), "\n")
	wantLines := strings.Split(string(want), "\n")
	maxLen := len(gotLines)
	if len(wantLines) > maxLen {
		maxLen = len(wantLines)
	}
	for i := 0; i < maxLen; i++ {
		var g, w string
		if i < len(gotLines) {
			g = gotLines[i]
		}
		if i < len(wantLines) {
			w = wantLines[i]
		}
		if g != w {
			t.Fatalf("%s mismatch at line %d:\n  want: %q\n  got:  %q\n  (want %d lines, got %d lines)",
				label, i+1, w, g, len(wantLines), len(gotLines))
		}
	}
	t.Fatalf("%s: bytes unequal but no per-line diff", label)
}

// assertCounters checks that the named counters in the Go JSON report
// match the upstream JSON report.  Path is a JSON path like
// "summary.before_filtering.total_reads".  We fail with a helpful
// message on mismatch.  Numbers are compared as float64 (which is how
// encoding/json decodes them) to avoid int/float coercion noise.
func assertCounters(t *testing.T, label string, upJSON, goJSON map[string]any, paths ...string) {
	t.Helper()
	for _, p := range paths {
		up, upOK := lookupJSON(upJSON, p)
		got, goOK := lookupJSON(goJSON, p)
		if !upOK || !goOK {
			t.Fatalf("%s: path %s missing (upstream=%v go=%v)", label, p, upOK, goOK)
		}
		upF, _ := toFloat(up)
		goF, _ := toFloat(got)
		if upF != goF {
			t.Fatalf("%s: %s mismatch: upstream=%v go=%v", label, p, up, got)
		}
	}
}

func lookupJSON(m map[string]any, path string) (any, bool) {
	parts := strings.Split(path, ".")
	cur := any(m)
	for _, p := range parts {
		mp, ok := cur.(map[string]any)
		if !ok {
			return nil, false
		}
		v, ok := mp[p]
		if !ok {
			return nil, false
		}
		cur = v
	}
	return cur, true
}

func toFloat(v any) (float64, bool) {
	switch x := v.(type) {
	case float64:
		return x, true
	case int:
		return float64(x), true
	case int64:
		return float64(x), true
	}
	return 0, false
}

// defaultEncoding returns the same quality encoding the CLI uses by
// default (Phred33).
func defaultEncoding() fastq.QualityEncoding {
	return fastq.Phred33
}

// ----------------------------------------------------------------------
// Case 1 — SE basic: defaults across the board, no filtering should fire.
//
// Upstream flags: (none)
// Validates: byte parity of out.fq, JSON before/after counters.
func TestParity_Fastp_Case01_SEBasic(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_clean.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir, nil)

	opts := DefaultProcessOptions()
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case01 SE clean FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case01 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"summary.before_filtering.total_bases",
		"summary.after_filtering.total_reads",
		"summary.after_filtering.total_bases",
		"filtering_result.passed_filter_reads",
		"filtering_result.low_quality_reads",
		"filtering_result.too_many_N_reads",
		"filtering_result.too_short_reads",
		"filtering_result.too_long_reads",
	)
}

// Case 2 — SE custom adapter: explicit -a AGATCGGAAGAGC.
//
// Upstream flags: -a AGATCGGAAGAGC
// Validates: byte parity of out.fq + adapter-trimmed counters.
func TestParity_Fastp_Case02_SEAdapter(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_adapter.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-a", "AGATCGGAAGAGC"})

	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case02 SE adapter FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case02 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// Case 3 — SE length filter (--length_required 50). After adapter-trim
// many adapter reads fall below 50bp and should be filtered.
//
// Upstream flags: -a AGATCGGAAGAGC -l 50
// Validates: too_short_reads parity + passed_filter_reads parity.
func TestParity_Fastp_Case03_SELengthFilter(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_adapter.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-a", "AGATCGGAAGAGC", "-l", "50"})

	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.MinLength = 50
	opts.LengthRequired = 50
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case03 SE length filter FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case03 counters", upJSON, goJSON,
		"filtering_result.passed_filter_reads",
		"filtering_result.too_short_reads",
	)
}

// Case 4 — SE max length (--length_limit 60). On 100bp clean reads
// every read is "too long" so passed=0.
//
// Upstream flags: --length_limit 60
// Validates: too_long_reads parity.
func TestParity_Fastp_Case04_SEMaxLength(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_clean.fq")
	dir := t.TempDir()

	_, upJSON := runUpstreamSE(t, bin, in, dir, []string{"--length_limit", "60"})

	opts := DefaultProcessOptions()
	opts.MaxLength = 60
	opts.LengthLimit = 60
	_, goStats := runGoFastpSE(t, in, opts)

	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case04 counters", upJSON, goJSON,
		"filtering_result.passed_filter_reads",
		"filtering_result.too_long_reads",
	)
}

// Case 5 — SE N filter (-n 2). Several reads have N counts above 2.
//
// Upstream flags: -n 2
// Validates: byte parity + too_many_N_reads counter.
func TestParity_Fastp_Case05_SENFilter(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_ns.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-n", "2"})

	opts := DefaultProcessOptions()
	opts.MaxNCount = 2
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case05 SE N filter FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case05 counters", upJSON, goJSON,
		"filtering_result.passed_filter_reads",
		"filtering_result.too_many_N_reads",
	)
}

// Case 6 — SE Q20/Q30 filter (-q 30 -u 20). 30% of bases must be Q30+.
//
// Upstream flags: -q 30 -u 20
// Validates: low_quality_reads + passed_filter_reads parity.
func TestParity_Fastp_Case06_SEQualityFilter(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_qmix.fq")
	dir := t.TempDir()

	_, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-q", "30", "-u", "20"})

	opts := DefaultProcessOptions()
	opts.QualThreshold = 30
	// Upstream's -u 20 means "20% of bases may be below threshold". Our
	// QualPercent is the inverse: "percent that must MEET the threshold"
	// — i.e. 100 - upstream's -u value.
	opts.QualPercent = 80
	_, goStats := runGoFastpSE(t, in, opts)

	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case06 counters", upJSON, goJSON,
		"filtering_result.passed_filter_reads",
		"filtering_result.low_quality_reads",
	)
}

// Case 7 — SE low-complexity filter (-y -Y 30). Validates the
// upstream-aligned complexity definition (fraction of adjacent
// differences).
//
// Upstream flags: -y -Y 30
// Validates: low_complexity_reads counter + byte parity.
//
// Note: the homopolymer fixture has 5 zero-complexity reads which both
// upstream and (now-correctly) Go reject.
func TestParity_Fastp_Case07_SELowComplexity(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_homopolymer.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-y", "-Y", "30"})

	opts := DefaultProcessOptions()
	opts.LowComplexity = true
	opts.ComplexityThreshold = 0.30
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case07 SE low-complexity FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case07 counters", upJSON, goJSON,
		"filtering_result.passed_filter_reads",
		"filtering_result.low_complexity_reads",
	)
}

// Case 8 — SE UMI extraction (--umi --umi_loc read1 --umi_len 8).
// Validates the upstream-aligned name suffix ":<umi>" (no implicit
// UMI_ prefix when --umi_prefix is unset).
//
// Upstream flags: -U --umi_loc read1 --umi_len 8
func TestParity_Fastp_Case08_SEUMI(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_umi.fq")
	dir := t.TempDir()

	upFq, _ := runUpstreamSE(t, bin, in, dir, []string{"-U", "--umi_loc", "read1", "--umi_len", "8"})

	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = "read1"
	opts.UMILen = 8
	goFq, _ := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case08 SE UMI FASTQ", goFq, upFq)
}

// Case 9 — SE UMI with --umi_prefix UMI: tag format becomes ":UMI_<umi>".
//
// Upstream flags: -U --umi_loc read1 --umi_len 8 --umi_prefix UMI
func TestParity_Fastp_Case09_SEUMIPrefix(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_umi.fq")
	dir := t.TempDir()

	upFq, _ := runUpstreamSE(t, bin, in, dir, []string{"-U", "--umi_loc", "read1", "--umi_len", "8", "--umi_prefix", "UMI"})

	opts := DefaultProcessOptions()
	opts.UMI = true
	opts.UMILoc = "read1"
	opts.UMILen = 8
	opts.UMIPrefix = "UMI"
	goFq, _ := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case09 SE UMI prefix FASTQ", goFq, upFq)
}

// Case 10 — PE basic: defaults across the board.
//
// Upstream flags: (none)
// Validates: byte parity of out R1 + R2, before/after counters.
func TestParity_Fastp_Case10_PEBasic(t *testing.T) {
	bin := ensureUpstream(t)
	r1 := parityInput(t, "pe_r1.fq")
	r2 := parityInput(t, "pe_r2.fq")
	dir := t.TempDir()

	upR1, upR2, upJSON := runUpstreamPE(t, bin, r1, r2, dir, nil)

	opts := DefaultProcessOptions()
	goR1, goR2, goStats := runGoFastpPE(t, r1, r2, opts)

	mustEqualBytes(t, "case10 PE R1", goR1, upR1)
	mustEqualBytes(t, "case10 PE R2", goR2, upR2)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case10 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"summary.before_filtering.total_bases",
		"summary.after_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// Case 11 — PE with explicit adapter (-a + --adapter_sequence_r2).
//
// Upstream flags: -a AGATCGGAAGAGC --adapter_sequence_r2 AGATCGGAAGAGC
func TestParity_Fastp_Case11_PEAdapter(t *testing.T) {
	bin := ensureUpstream(t)
	r1 := parityInput(t, "pe_r1.fq")
	r2 := parityInput(t, "pe_r2.fq")
	dir := t.TempDir()

	_, _, upJSON := runUpstreamPE(t, bin, r1, r2, dir,
		[]string{"-a", "AGATCGGAAGAGC", "--adapter_sequence_r2", "AGATCGGAAGAGC"})

	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	_, _, goStats := runGoFastpPE(t, r1, r2, opts)

	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case11 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"summary.before_filtering.total_bases",
		"filtering_result.passed_filter_reads",
	)
}

// Case 11b — PE with --detect_adapter_for_pe (overlap-based adapter
// detection between mates). Validates byte parity of R1+R2 and the
// summary counters.  This works because the PE adapter-detection path
// uses overlap analysis (which we DO implement) rather than the
// SE kmer-tree path.
//
// Upstream flags: --detect_adapter_for_pe
func TestParity_Fastp_Case11b_PEDetectAdapter(t *testing.T) {
	bin := ensureUpstream(t)
	r1 := parityInput(t, "pe_r1.fq")
	r2 := parityInput(t, "pe_r2.fq")
	dir := t.TempDir()

	upR1, upR2, upJSON := runUpstreamPE(t, bin, r1, r2, dir, []string{"--detect_adapter_for_pe"})

	opts := DefaultProcessOptions()
	opts.DetectAdapterPE = true
	goR1, goR2, goStats := runGoFastpPE(t, r1, r2, opts)

	mustEqualBytes(t, "case11b PE R1 (detect_adapter_for_pe)", goR1, upR1)
	mustEqualBytes(t, "case11b PE R2 (detect_adapter_for_pe)", goR2, upR2)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case11b counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// Case 12 — SE poly-G trimming (-g flag, upstream's NovaSeq-style trim).
// Our Go port now runs the verbatim PolyX::trimPolyG algorithm from
// reference_code/fastp/src/polyx.cpp:16-42 (1 mismatch per 8 bases scanned,
// capped at 5; anchors on the last-G seen).
func TestParity_Fastp_Case12_SEPolyG(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_clean.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-g"})

	opts := DefaultProcessOptions()
	opts.TrimPolyG = true
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case12 SE poly-G FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case12 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// Case 13 — SE sliding-window cut_right (--cut_right with default window=4
// and quality=20). Our Go port now mirrors the upstream cut_right walk
// from reference_code/fastp/src/filter.cpp:144-178, including the
// high-Q-prefix preservation inside the offending window.
func TestParity_Fastp_Case13_SECutRight(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_lqtail.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir,
		[]string{"--cut_right", "--cut_right_window_size", "4", "--cut_right_mean_quality", "20"})

	opts := DefaultProcessOptions()
	opts.CutRight = true
	opts.CutWindowSize = 4
	opts.CutMeanQuality = 20
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case13 SE cut_right FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case13 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// Case 14 — SE sliding-window cut_front + cut_tail. Tests both ends of
// upstream's symmetric trimAndCut algorithm (filter.cpp:111-209),
// including the trailing-N skip step at filter.cpp:138-139 / 206-207.
func TestParity_Fastp_Case14_SECutFrontTail(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_lqtail.fq")
	dir := t.TempDir()

	upFq, upJSON := runUpstreamSE(t, bin, in, dir,
		[]string{"--cut_front", "--cut_tail",
			"--cut_front_window_size", "4", "--cut_front_mean_quality", "20",
			"--cut_tail_window_size", "4", "--cut_tail_mean_quality", "20"})

	opts := DefaultProcessOptions()
	opts.CutFront = true
	opts.CutTail = true
	opts.CutWindowSize = 4
	opts.CutMeanQuality = 20
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case14 SE cut_front+cut_tail FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	assertCounters(t, "case14 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// Case 15 — SE adapter auto-detection (default upstream flags: -a is
// "auto" so upstream runs Evaluator::evalAdapterAndReadNum). Our Go
// port now uses the verbatim upstream algorithm (a kmer overlap-tree
// over the first ~10000 reads) instead of a simple known-adapter
// substring search. The se_adapter.fq fixture has 20 reads, which is
// below upstream's gate at evaluator.cpp:344 (records >= 10000), so
// both upstream and our port return "" from the evaluator and the
// FASTQ output is bit-identical to case 1 (no adapter trimming
// applied). This case validates the SE auto-detect code path runs
// to completion and produces byte-parity output.
func TestParity_Fastp_Case15_SEAutoDetect(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_adapter.fq")
	dir := t.TempDir()

	// Upstream's default: -a defaults to "auto" so the evaluator runs
	// unconditionally. No extra flags needed.
	upFq, upJSON := runUpstreamSE(t, bin, in, dir, nil)

	opts := DefaultProcessOptions()
	// Mirror upstream's default: enable SE adapter detection.
	opts.DetectAdapterSE = true
	goFq, goStats := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "case15 SE adapter auto-detect FASTQ", goFq, upFq)
	goJSON := jsonFromStats(t, goStats)
	// adapter_cutting.* counters are only present when upstream detects
	// an adapter; for this small fixture detection returns "" so we
	// only check the always-present summary counters.
	assertCounters(t, "case15 counters", upJSON, goJSON,
		"summary.before_filtering.total_reads",
		"filtering_result.passed_filter_reads",
	)
}

// jsonFromStats serialises a ProcessStats through WriteJSONReport so the
// counter-comparison logic operates on the same shape that the CLI
// emits. Internally writes to a temp file because WriteJSONReport
// targets a path.
func jsonFromStats(t *testing.T, stats *ProcessStats) map[string]any {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "go_report.json")
	if err := WriteJSONReport(p, stats); err != nil {
		t.Fatalf("WriteJSONReport: %v", err)
	}
	return readJSON(t, p)
}
