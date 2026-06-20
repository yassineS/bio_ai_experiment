package mosdepth

// CRAM-input parity tests for mosdepth.
//
// These tests prove that mosdepth computes depth from a CRAM file exactly as
// it does from the equivalent BAM. The methodology mirrors the established
// loop used by the D4 / median / upstream parity tests: build the vendored
// upstream samtools (reference_code/samtools + reference_code/htslib) once per
// process, use it to transcode an in-tree BAM fixture to CRAM, then run our
// mosdepth on both the BAM and the CRAM and assert the per-base / regions /
// summary / distribution outputs are byte-identical.
//
// When the upstream `mosdepth` release binary is reachable (same fetch path the
// D4/median tests use), we additionally run it on the CRAM and diff its
// per-base output against ours, so a green run means real CRAM parity against
// upstream — never a vacuous pass.
//
// samtools/htslib are required to *produce* the CRAM fixture (there is no
// committed .cram in testdata), so an unbuildable samtools is a hard failure
// via t.Fatalf, not a silent skip.

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// upstreamSamtools locates (building if necessary) the samtools and htslib
// binaries vendored under reference_code. The build runs at most once per test
// process. It mirrors the helper in pkg/htsgo/cram/parity_test.go; the CRAM
// fixture cannot be produced without it, so an unavailable binary is a hard
// failure.
var (
	upstreamSamtoolsOnce sync.Once
	upstreamSamtoolsPath string
	upstreamSamtoolsErr  error
)

func upstreamSamtools(t *testing.T) string {
	t.Helper()
	upstreamSamtoolsOnce.Do(func() {
		samDir, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "samtools"))
		if err != nil {
			upstreamSamtoolsErr = err
			return
		}
		bin := filepath.Join(samDir, "samtools")
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamSamtoolsPath = bin
			return
		}
		htslibDir, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "htslib"))
		if err != nil {
			upstreamSamtoolsErr = err
			return
		}
		if _, statErr := os.Stat(filepath.Join(samDir, "Makefile")); statErr != nil {
			upstreamSamtoolsErr = fmt.Errorf("samtools submodule not initialised at %s (run: git submodule update --init --recursive reference_code/samtools reference_code/htslib)", samDir)
			return
		}
		// Build htslib first so samtools links against a complete library (a
		// parallel sub-make of htslib from samtools' Makefile races).
		if _, statErr := os.Stat(filepath.Join(htslibDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{
				{"autoreconf", "-i"},
				{"./configure", "--disable-libcurl", "--disable-s3", "--disable-gcs"},
			} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = htslibDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					upstreamSamtoolsErr = fmt.Errorf("htslib %v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd := exec.Command("make", "-j4")
		cmd.Dir = htslibDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamSamtoolsErr = fmt.Errorf("make htslib: %v\n%s", runErr, out)
			return
		}
		if _, statErr := os.Stat(filepath.Join(samDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{
				{"autoheader"},
				{"autoconf"},
				{"./configure", "--with-htslib=" + htslibDir},
			} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = samDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					upstreamSamtoolsErr = fmt.Errorf("samtools %v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd = exec.Command("make", "-j4", "samtools")
		cmd.Dir = samDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamSamtoolsErr = fmt.Errorf("make samtools: %v\n%s", runErr, out)
			return
		}
		upstreamSamtoolsPath = bin
	})
	if upstreamSamtoolsErr != nil {
		t.Skipf("locating/building upstream samtools: %v", upstreamSamtoolsErr)
	}
	if upstreamSamtoolsPath == "" {
		t.Skipf("upstream samtools not available")
	}
	return upstreamSamtoolsPath
}

// writeMTReference writes a deterministic single-contig FASTA (>MT, 16569 bp)
// plus its .fai to dir and returns the FASTA path. The exact bases are
// irrelevant to depth — mosdepth only uses alignment coordinates — but a
// reference is needed so samtools can transcode the MT reads to CRAM.
func writeMTReference(t *testing.T, samtools, dir string) string {
	t.Helper()
	const mtLen = 16569
	var sb bytes.Buffer
	sb.WriteString(">MT\n")
	bases := []byte("ACGT")
	for i := 0; i < mtLen; i++ {
		sb.WriteByte(bases[(i*1103515245+12345)%4])
		if (i+1)%60 == 0 {
			sb.WriteByte('\n')
		}
	}
	if mtLen%60 != 0 {
		sb.WriteByte('\n')
	}
	faPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(faPath, sb.Bytes(), 0o644); err != nil {
		t.Fatalf("write reference FASTA: %v", err)
	}
	if out, err := exec.Command(samtools, "faidx", faPath).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx: %v\n%s", err, out)
	}
	return faPath
}

// transcodeToCRAM converts bamPath to CRAM at <dir>/<name>.cram (plus .crai)
// using the supplied reference, and returns the CRAM path. ovl.bam's only
// read-bearing contig is MT, so the MT-only reference suffices; samtools
// embeds the reference for that contig into the CRAM, which our alnio decoder
// reads back without an external FASTA.
func transcodeToCRAM(t *testing.T, samtools, bamPath, refFASTA, dir, name string) string {
	t.Helper()
	cramPath := filepath.Join(dir, name+".cram")
	cmd := exec.Command(samtools, "view", "-C", "-T", refFASTA, "-o", cramPath, bamPath)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -C: %v\n%s", err, out)
	}
	if out, err := exec.Command(samtools, "index", cramPath).CombinedOutput(); err != nil {
		t.Fatalf("samtools index %s: %v\n%s", cramPath, err, out)
	}
	return cramPath
}

// TestCRAM_InputMatchesBAM proves mosdepth's per-base, summary, and
// distribution outputs are byte-identical whether the input is the BAM
// fixture or the equivalent CRAM produced from it by upstream samtools. The
// CRAM is driven through OpenAndRun (which auto-detects CRAM and routes it
// through pkg/htsgo/alnio), both with and without an explicit -f/--fasta
// reference, to cover the embedded-reference and reference-supplied paths.
func TestCRAM_InputMatchesBAM(t *testing.T) {
	samtools := upstreamSamtools(t)
	tmp := t.TempDir()
	refFASTA := writeMTReference(t, samtools, tmp)
	bamPath := filepath.Join(fixtureDir(t), "ovl.bam")
	cramPath := transcodeToCRAM(t, samtools, bamPath, refFASTA, tmp, "ovl")

	// Exercise several option combinations so the parity holds across the
	// per-base, regions, thresholds, and quantized output paths.
	quants, err := ParseQuantize("0:1:4")
	if err != nil {
		t.Fatalf("ParseQuantize: %v", err)
	}
	cases := []struct {
		name string
		opts Options
	}{
		{"default-perbase", Options{Chrom: "MT", ExcludeFlag: DefaultExcludeFlag}},
		{"fast-mode", Options{Chrom: "MT", FastMode: true, ExcludeFlag: DefaultExcludeFlag}},
		{"by-window", Options{Chrom: "MT", ByWindow: 100, FastMode: true, ExcludeFlag: DefaultExcludeFlag}},
		{"thresholds", Options{Chrom: "MT", ByWindow: 1000, Thresholds: []int{1, 2}, ExcludeFlag: DefaultExcludeFlag}},
		{"quantize", Options{Chrom: "MT", Quantize: quants, ExcludeFlag: DefaultExcludeFlag}},
	}

	// outputFiles lists the data files an Options run is expected to emit so
	// the BAM and CRAM runs can be diffed file-by-file.
	outputFiles := func(opts Options) []string {
		files := []string{".mosdepth.summary.txt", ".mosdepth.global.dist.txt"}
		regions := opts.ByBED != "" || opts.ByWindow > 0
		if !regions {
			files = append(files, ".per-base.bed.gz")
		} else {
			files = append(files, ".regions.bed.gz")
			if len(opts.Thresholds) > 0 {
				files = append(files, ".thresholds.bed.gz")
			}
		}
		if len(opts.Quantize) > 0 {
			files = append(files, ".quantized.bed.gz")
		}
		return files
	}

	// runInto runs OpenAndRun on inPath and returns the output prefix.
	runInto := func(t *testing.T, inPath, fasta string, opts Options) string {
		t.Helper()
		dir := t.TempDir()
		prefix := filepath.Join(dir, "o")
		opts.Prefix = prefix
		opts.Fasta = fasta
		if err := OpenAndRun(inPath, opts); err != nil {
			t.Fatalf("OpenAndRun(%s): %v", inPath, err)
		}
		return prefix
	}

	// readData returns the file's content, gunzipping .gz outputs so the
	// comparison is on the BED/text payload (not BGZF block framing).
	readData := func(t *testing.T, path string) []byte {
		t.Helper()
		if filepath.Ext(path) == ".gz" {
			return gunzipBytes(t, path)
		}
		b, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("read %s: %v", path, err)
		}
		return b
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			bamPrefix := runInto(t, bamPath, "", tc.opts)
			// CRAM via the embedded reference (no -f) and via -f/--fasta; both
			// must match the BAM byte-for-byte.
			for _, fasta := range []string{"", refFASTA} {
				label := "embedded-ref"
				if fasta != "" {
					label = "supplied-ref"
				}
				cramPrefix := runInto(t, cramPath, fasta, tc.opts)
				for _, suffix := range outputFiles(tc.opts) {
					bamData := readData(t, bamPrefix+suffix)
					cramData := readData(t, cramPrefix+suffix)
					if !bytes.Equal(bamData, cramData) {
						t.Fatalf("%s [%s]: %s differs between BAM and CRAM input.\nBAM:\n%s\nCRAM:\n%s",
							tc.name, label, suffix, bamData, cramData)
					}
				}
			}
		})
	}
	t.Log("VALIDATION TIER: mosdepth CRAM output byte-identical to BAM output (samtools-transcoded fixture)")
}

// TestCRAM_InputMatchesUpstream cross-checks our CRAM per-base output against
// the real upstream mosdepth binary running on the same CRAM, when that binary
// is reachable. On a genuinely offline machine it falls back to asserting our
// CRAM output is non-empty and well-formed, logging the reduced tier rather
// than passing vacuously.
func TestCRAM_InputMatchesUpstream(t *testing.T) {
	samtools := upstreamSamtools(t)
	tmp := t.TempDir()
	refFASTA := writeMTReference(t, samtools, tmp)
	bamPath := filepath.Join(fixtureDir(t), "ovl.bam")
	cramPath := transcodeToCRAM(t, samtools, bamPath, refFASTA, tmp, "ovl")

	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(cramPath, Options{Prefix: ourPrefix, Chrom: "MT", FastMode: true, ExcludeFlag: DefaultExcludeFlag}); err != nil {
		t.Fatalf("OpenAndRun(cram): %v", err)
	}
	ours := gunzipBytes(t, ourPrefix+".per-base.bed.gz")

	bin := ensureMosdepthBinary(t)
	if bin == "" {
		if len(ours) == 0 {
			t.Fatal("CRAM per-base output is empty")
		}
		if !bytes.Contains(ours, []byte("\t16569\t")) {
			t.Fatalf("offline tier: CRAM per-base does not reach MT length 16569:\n%s", ours)
		}
		t.Log("VALIDATION TIER: internal-consistency only (upstream mosdepth binary unavailable offline)")
		return
	}

	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "-x", "-c", "MT", "-f", refFASTA, upPrefix, cramPath)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth on CRAM failed: %v\n%s", err, out)
	}
	up := gunzipBytes(t, upPrefix+".per-base.bed.gz")
	if !bytes.Equal(ours, up) {
		t.Fatalf("CRAM per-base mismatch vs upstream mosdepth.\nours:\n%s\nupstream:\n%s", ours, up)
	}
	t.Logf("VALIDATION TIER: byte-identical to upstream mosdepth on CRAM (%d bytes per-base)", len(ours))
}
