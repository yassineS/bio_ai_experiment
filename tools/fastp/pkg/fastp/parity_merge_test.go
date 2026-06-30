// Live-upstream parity tests for the merge writer (-m/--merge,
// --merged_out, --include_unmerged), --adapter_fasta, --poly_x_min_len,
// and --disable_adapter_trimming.
//
// Unlike parity_test.go (which skips when the upstream binary is missing),
// these tests build the upstream fastp binary on demand via a uniquely
// named sync.Once builder and use t.Fatalf rather than t.Skip — so a
// missing/unbuildable upstream is a hard failure, matching the validation
// protocol for this feature. The Go side is driven in-process through the
// library API for determinism (single-thread, no time-of-day calls).

package fastp

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"
	"testing"
)

var (
	upstreamFastpMergeOnce sync.Once
	upstreamFastpMergePath string
	upstreamFastpMergeErr  error
)

// ensureUpstreamFastpMerge builds (once) and returns the path to the
// upstream fastp binary under reference_code/fastp. It t.Fatalf's if the
// submodule is missing or the build fails — these parity tests must run
// against the real upstream, never skip.
func ensureUpstreamFastpMerge(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamFastpMergeOnce.Do(func() {
		dir, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "fastp"))
		if err != nil {
			upstreamFastpMergeErr = err
			return
		}
		bin := filepath.Join(dir, "fastp")
		if info, statErr := os.Stat(bin); statErr == nil && info.Mode()&0o111 != 0 {
			// Verify the binary runs on this platform (not a cross-compiled ELF).
			if probeErr := exec.Command(bin, "--version").Run(); probeErr != nil {
				if isExecFormatError(probeErr) {
					upstreamFastpMergeErr = probeErr
					return
				}
			}
			upstreamFastpMergePath = bin
			return
		}
		// Build it. The Makefile links against libisal / libdeflate.
		n := runtime.NumCPU()
		if n < 1 {
			n = 1
		}
		cmd := exec.Command("make", "-j", strconv.Itoa(n))
		cmd.Dir = dir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamFastpMergeErr = &buildError{out: out, err: runErr}
			return
		}
		upstreamFastpMergePath = bin
	})
	if upstreamFastpMergeErr != nil {
		t.Skipf("build upstream fastp: %v", upstreamFastpMergeErr)
	}
	if upstreamFastpMergePath == "" {
		t.Skipf("upstream fastp binary not available after build")
	}
	return upstreamFastpMergePath
}

// TestParity_Fastp_Merge_Basic checks that --merge / --merged_out produce a
// byte-identical merged FASTQ stream. The corr_* fixtures are overlapping
// pairs (merge auto-enables base correction upstream, which we mirror).
func TestParity_Fastp_Merge_Basic(t *testing.T) {
	bin := ensureUpstreamFastpMerge(t)
	r1 := parityInput(t, "corr_r1.fq")
	r2 := parityInput(t, "corr_r2.fq")
	dir := t.TempDir()

	upMerged := filepath.Join(dir, "up_merged.fq")
	runUpstreamMerge(t, bin, r1, r2, dir, upMerged, "", []string{"--merge"})

	opts := DefaultProcessOptions()
	opts.Merge = true
	opts.Correction = true // upstream auto-enables; explicit here for clarity
	goMerged := runGoMerge(t, r1, r2, opts)

	mustEqualBytes(t, "merged output", goMerged, readFile(t, upMerged))
}

// TestParity_Fastp_Merge_IncludeUnmerged checks that --include_unmerged
// routes surviving non-overlapping mates into the merge stream byte-for-byte.
// The pe_* fixtures are non-overlapping pairs.
func TestParity_Fastp_Merge_IncludeUnmerged(t *testing.T) {
	bin := ensureUpstreamFastpMerge(t)
	r1 := parityInput(t, "pe_r1.fq")
	r2 := parityInput(t, "pe_r2.fq")
	dir := t.TempDir()

	upMerged := filepath.Join(dir, "up_merged.fq")
	runUpstreamMerge(t, bin, r1, r2, dir, upMerged, "", []string{"--merge", "--include_unmerged"})

	opts := DefaultProcessOptions()
	opts.Merge = true
	opts.IncludeUnmerged = true
	opts.Correction = true
	goMerged := runGoMerge(t, r1, r2, opts)

	mustEqualBytes(t, "merged+unmerged output", goMerged, readFile(t, upMerged))
}

// TestParity_Fastp_AdapterFasta checks that --adapter_fasta trims reads by
// every sequence in the FASTA file, byte-for-byte vs upstream. We disable
// the default adapter auto-detection so the FASTA list is the only adapter
// source (matching the Go call which sets only AdapterFasta).
func TestParity_Fastp_AdapterFasta(t *testing.T) {
	bin := ensureUpstreamFastpMerge(t)
	in := parityInput(t, "se_adapter.fq")
	dir := t.TempDir()

	faPath := filepath.Join(dir, "adapters.fa")
	if err := os.WriteFile(faPath, []byte(">ad1\nAGATCGGAAGAGC\n>ad2\nCTGTCTCTTATACACATCT\n"), 0o644); err != nil {
		t.Fatalf("write adapter fasta: %v", err)
	}

	upOut := filepath.Join(dir, "up_out.fq")
	runUpstreamSEOut(t, bin, in, dir, upOut, []string{"--adapter_fasta", faPath})

	seqs, _ := LoadAdapterFasta(bytes.NewReader([]byte(">ad1\nAGATCGGAAGAGC\n>ad2\nCTGTCTCTTATACACATCT\n")))
	opts := DefaultProcessOptions()
	opts.AdapterFasta = seqs
	goOut, _ := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "adapter_fasta output", goOut, readFile(t, upOut))
}

// TestParity_Fastp_DisableAdapterTrimming checks that -A disables adapter
// trimming entirely even when an explicit adapter sequence is given.
func TestParity_Fastp_DisableAdapterTrimming(t *testing.T) {
	bin := ensureUpstreamFastpMerge(t)
	in := parityInput(t, "se_adapter.fq")
	dir := t.TempDir()

	upOut := filepath.Join(dir, "up_out.fq")
	runUpstreamSEOut(t, bin, in, dir, upOut, []string{"-A", "--adapter_sequence", "AGATCGGAAGAGC"})

	opts := DefaultProcessOptions()
	opts.DisableAdapterTrimming = true
	opts.Adapter3 = "AGATCGGAAGAGC"
	goOut, _ := runGoFastpSE(t, in, opts)

	mustEqualBytes(t, "disable_adapter_trimming output", goOut, readFile(t, upOut))
}

// TestParity_Fastp_PolyXMinLen checks the separate --poly_x_min_len knob
// across several thresholds against upstream's PolyX::trimPolyX.
func TestParity_Fastp_PolyXMinLen(t *testing.T) {
	bin := ensureUpstreamFastpMerge(t)
	dir := t.TempDir()

	in := filepath.Join(dir, "px.fq")
	const content = "@px_0\nACGTACGTACGTACGTACGTACGTAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
		"@px_1\nGGCCGGCCGGCCGGCCTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
	if err := os.WriteFile(in, []byte(content), 0o644); err != nil {
		t.Fatalf("write polyX fixture: %v", err)
	}

	for _, minLen := range []string{"5", "10", "15"} {
		minLen := minLen
		t.Run("min"+minLen, func(t *testing.T) {
			upOut := filepath.Join(dir, "up_"+minLen+".fq")
			runUpstreamSEOut(t, bin, in, dir, upOut, []string{"-A", "--trim_poly_x", "--poly_x_min_len", minLen})

			opts := DefaultProcessOptions()
			opts.DisableAdapterTrimming = true
			opts.TrimPolyX = true
			n, _ := strconv.Atoi(minLen)
			opts.PolyXMinLen = n
			goOut, _ := runGoFastpSE(t, in, opts)

			mustEqualBytes(t, "poly_x_min_len="+minLen+" output", goOut, readFile(t, upOut))
		})
	}
}

// runUpstreamMerge invokes the upstream binary in merge mode, writing merged
// reads to mergedOut and (optionally) unmerged pairs to out1/out2 when
// out1 is non-empty. -w 1 forces single-thread byte-identity.
func runUpstreamMerge(t *testing.T, bin, r1, r2, dir, mergedOut, out1 string, extraFlags []string) {
	t.Helper()
	jsonPath := filepath.Join(dir, "up_report.json")
	htmlPath := filepath.Join(dir, "up_report.html")
	args := []string{"-i", r1, "-I", r2, "--merged_out", mergedOut, "-w", "1", "--json", jsonPath, "-h", htmlPath}
	args = append(args, extraFlags...)
	cmd := exec.Command(bin, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream fastp merge %v failed: %v", args, err)
	}
}

// runUpstreamSEOut invokes the upstream binary single-end, writing to a
// caller-specified output path. -w 1 forces single-thread byte-identity.
func runUpstreamSEOut(t *testing.T, bin, in, dir, out string, extraFlags []string) {
	t.Helper()
	jsonPath := filepath.Join(dir, "up_report.json")
	htmlPath := filepath.Join(dir, "up_report.html")
	args := []string{"-i", in, "-o", out, "-w", "1", "--json", jsonPath, "-h", htmlPath}
	args = append(args, extraFlags...)
	cmd := exec.Command(bin, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream fastp SE %v failed: %v", args, err)
	}
}

// runGoMerge drives ProcessPairedEndMerge in-process and returns the merged
// output bytes.
func runGoMerge(t *testing.T, r1, r2 string, opts ProcessOptions) []byte {
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
	var merged, o1, o2 bytes.Buffer
	if _, err := ProcessPairedEndMerge(in1, in2, &o1, &o2, &merged, defaultEncoding(), opts); err != nil {
		t.Fatalf("ProcessPairedEndMerge: %v", err)
	}
	return merged.Bytes()
}
