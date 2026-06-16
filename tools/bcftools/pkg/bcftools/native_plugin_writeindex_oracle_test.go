package bcftools

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// Live-oracle parity tests for the -W/--write-index support added to the native
// plugins contrast, isecGT, mendelian2, scatter and split. Each test drives BOTH
// the real upstream bcftools and OUR port through their CLIs with the same
// upstream-accepted argv, then validates the produced indexes against upstream.
//
// Index files (.csi / .tbi) are BGZF-wrapped, and our DEFLATE backend
// (klauspost) frames blocks differently from upstream zlib, so the raw on-disk
// index bytes can never be byte-identical even for identical index content (the
// data file's framing differs too). The meaningful, rigorous comparison is
// therefore: the DECOMPRESSED index content our plugin writes for a given output
// file must equal what upstream `bcftools index` produces for that SAME file.
// readIndexContent decodes the index, and indexOracle re-indexes our produced
// data file with the real upstream binary and diffs the decoded bytes. Combined
// with the data-content parity (decoded with `bcftools view`), this validates
// the -W feature end-to-end against upstream 1.23.1.

// readIndexContent returns the BGZF-decompressed bytes of a .csi/.tbi index
// file, i.e. the index structure independent of the gzip framing.
func readIndexContent(t *testing.T, path string) []byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip %s: %v", path, err)
	}
	gr.Multistream(true)
	data, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return data
}

// upstreamIndexContent runs the real upstream `bcftools index` over dataFile
// (CSI by default, TBI when wantTBI) and returns the decompressed index bytes.
// It is the oracle: our plugin's index must match this for the identical file.
func upstreamIndexContent(t *testing.T, bin, dataFile string, wantTBI bool) []byte {
	t.Helper()
	dir := t.TempDir()
	local := filepath.Join(dir, filepath.Base(dataFile))
	src, err := os.ReadFile(dataFile)
	if err != nil {
		t.Fatalf("read %s: %v", dataFile, err)
	}
	if err := os.WriteFile(local, src, 0o644); err != nil {
		t.Fatalf("write %s: %v", local, err)
	}
	argv := []string{"index", "-f"}
	suffix := ".csi"
	if wantTBI {
		argv = append(argv, "-t")
		suffix = ".tbi"
	}
	argv = append(argv, local)
	cmd := exec.Command(bin, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream index %v: %v\n%s", argv, err, out)
	}
	return readIndexContent(t, local+suffix)
}

// dataAndIndexFiles partitions a directory's regular files into data files
// (.vcf.gz / .bcf) and index files (.csi / .tbi), returning their base names.
func dataAndIndexFiles(t *testing.T, dir string) (data, index []string) {
	t.Helper()
	ents, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range ents {
		if e.IsDir() {
			continue
		}
		name := e.Name()
		switch {
		case strings.HasSuffix(name, ".csi"), strings.HasSuffix(name, ".tbi"):
			index = append(index, name)
		case strings.HasSuffix(name, ".vcf.gz"), strings.HasSuffix(name, ".bcf"):
			data = append(data, name)
		}
	}
	sort.Strings(data)
	sort.Strings(index)
	return data, index
}

// assertWriteIndexParity runs upstream and ours into separate directories with
// the given argv (a function of the output directory), then:
//   - asserts the produced file sets are identical (data + index names);
//   - for each data file, asserts our decoded content matches upstream's
//     (provenance stripped);
//   - for each index, asserts our decoded index bytes match what the upstream
//     `bcftools index` binary produces for OUR data file — a true byte-level
//     index oracle that is robust to the unavoidable BGZF framing difference.
func assertWriteIndexParity(t *testing.T, bin string, wantTBI bool, mkArgv func(dir string) []string) {
	t.Helper()
	upDir := t.TempDir()
	ourDir := t.TempDir()
	_ = runUpstreamPluginCmd(t, bin, mkArgv(upDir))
	_ = runOursPluginCmd(t, mkArgv(ourDir))

	upData, upIdx := dataAndIndexFiles(t, upDir)
	ourData, ourIdx := dataAndIndexFiles(t, ourDir)

	if strings.Join(upData, ",") != strings.Join(ourData, ",") {
		t.Fatalf("data file set differs: upstream %v vs ours %v", upData, ourData)
	}
	if strings.Join(upIdx, ",") != strings.Join(ourIdx, ",") {
		t.Fatalf("index file set differs: upstream %v vs ours %v", upIdx, ourIdx)
	}
	if len(ourIdx) == 0 {
		t.Fatalf("no index files produced (ours): %v", ourData)
	}

	// Data content parity (decoded, provenance-stripped).
	for _, name := range ourData {
		up := canonicalizeFile(t, bin, filepath.Join(upDir, name))
		our := canonicalizeFile(t, bin, filepath.Join(ourDir, name))
		if !bytes.Equal(stripProvenanceBytes(up), stripProvenanceBytes(our)) {
			t.Fatalf("data file %q diverges\n--- upstream ---\n%s\n--- ours ---\n%s",
				name, snippet(up, 1500), snippet(our, 1500))
		}
	}

	// Index byte parity: our index decoded == upstream `bcftools index` over our
	// produced data file decoded.
	for _, name := range ourIdx {
		dataName := strings.TrimSuffix(strings.TrimSuffix(name, ".csi"), ".tbi")
		dataPath := filepath.Join(ourDir, dataName)
		if _, err := os.Stat(dataPath); err != nil {
			t.Fatalf("index %q has no sibling data file %q", name, dataName)
		}
		ours := readIndexContent(t, filepath.Join(ourDir, name))
		want := upstreamIndexContent(t, bin, dataPath, wantTBI)
		if !bytes.Equal(ours, want) {
			t.Fatalf("index %q content diverges from upstream bcftools index (len ours=%d want=%d)",
				name, len(ours), len(want))
		}
	}
}

// TestNativePluginContrastWriteIndex checks contrast's -W output indexing for
// VCF.gz/BCF CSI and VCF.gz TBI, plus the non-indexable / stdout error cases.
func TestNativePluginContrastWriteIndex(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	// VCF.gz cases use gt_plugins.vcf (its phased-missing GT survives the VCF
	// text round-trip). BCF cases use split_multi.vcf, which carries no
	// phased-missing genotype — our BCF encoder canonicalises ".|." to "./."
	// (a documented encoder behaviour independent of -W), so a fixture without
	// one keeps the BCF data-content comparison about the plugin, not the
	// encoder.
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "split_multi.vcf")
	base := []string{"+contrast", "-0", "S1,S2", "-1", "S3,S4"}

	cases := []struct {
		name    string
		fix     string
		oFile   string
		wflag   string
		oType   string
		wantTBI bool
	}{
		{"vcfgz_csi", gt, "out.vcf.gz", "-W", "-Oz", false},
		{"vcfgz_csi_explicit", gt, "out.vcf.gz", "-W=csi", "-Oz", false},
		{"vcfgz_tbi", gt, "out.vcf.gz", "-W=tbi", "-Oz", true},
		{"bcf_csi", multi, "out.bcf", "-W", "-Ob", false},
		{"bcf_tbi_fallback", multi, "out.bcf", "-W=tbi", "-Ob", false}, // upstream writes CSI for BCF
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertWriteIndexParity(t, bin, tc.wantTBI, func(dir string) []string {
				argv := append([]string{}, base...)
				argv = append(argv, tc.oType, "-o", filepath.Join(dir, tc.oFile), tc.wflag, tc.fix)
				return argv
			})
		})
	}

	// Non-indexable output (-Ov) and stdout (-W with no -o) must error on BOTH
	// binaries with a non-zero exit.
	for _, tc := range []struct {
		name string
		args []string
	}{
		{"plain_vcf_errors", []string{"-Ov", "-o", filepath.Join(t.TempDir(), "out.vcf"), "-W", gt}},
		{"stdout_errors", []string{"-Oz", "-W", gt}},
	} {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			argv := append(append([]string{}, base...), tc.args...)
			assertPluginWriteIndexErrors(t, bin, argv)
		})
	}
}

// TestNativePluginIsecGTWriteIndex checks isecGT's -W output indexing.
func TestNativePluginIsecGTWriteIndex(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	aPlain := parityFixture(t, "isecgt_a.vcf")
	bPlain := parityFixture(t, "isecgt_b.vcf")
	// Upstream's synced reader needs indexed bgzipped inputs; our port reads the
	// same bgzipped inputs transparently, so both binaries use the identical
	// indexed input pair (the output and its index are independent of the input
	// container anyway).
	prep := t.TempDir()
	aGz := bgzipAndIndex(t, bin, aPlain, filepath.Join(prep, "a.vcf.gz"))
	bGz := bgzipAndIndex(t, bin, bPlain, filepath.Join(prep, "b.vcf.gz"))

	for _, tc := range []struct {
		name    string
		oFile   string
		wflag   string
		oType   string
		wantTBI bool
	}{
		{"vcfgz_csi", "out.vcf.gz", "-W", "-Oz", false},
		{"vcfgz_tbi", "out.vcf.gz", "-W=tbi", "-Oz", true},
		{"bcf_csi", "out.bcf", "-W", "-Ob", false},
	} {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertWriteIndexParity(t, bin, tc.wantTBI, func(dir string) []string {
				out := filepath.Join(dir, tc.oFile)
				return []string{"+isecGT", aGz, bGz, tc.oType, "-o", out, tc.wflag}
			})
		})
	}
}

// TestNativePluginSplitWriteIndex checks split's per-output -W indexing (one
// index file PER sample output), plus the non-indexable error.
func TestNativePluginSplitWriteIndex(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "split_multi.vcf")
	for _, tc := range []struct {
		name    string
		fix     string
		wflag   string
		oType   string
		wantTBI bool
	}{
		{"vcfgz_csi", gt, "-W", "-Oz", false},
		{"vcfgz_tbi", gt, "-W=tbi", "-Oz", true},
		{"bcf_csi", multi, "-W", "-Ob", false},
	} {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertWriteIndexParity(t, bin, tc.wantTBI, func(dir string) []string {
				return []string{"+split", "-o", dir, tc.oType, tc.wflag, tc.fix}
			})
		})
	}
	t.Run("plain_vcf_errors", func(t *testing.T) {
		assertPluginWriteIndexErrors(t, bin, []string{"+split", "-o", t.TempDir(), "-W", gt})
	})
}

// TestNativePluginScatterWriteIndex checks scatter's per-output -W indexing
// (one index file PER scattered chunk).
func TestNativePluginScatterWriteIndex(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "split_multi.vcf")
	for _, tc := range []struct {
		name    string
		fix     string
		wflag   string
		oType   string
		wantTBI bool
	}{
		{"chunks_vcfgz_csi", gt, "-W", "-Oz", false},
		{"chunks_vcfgz_tbi", gt, "-W=tbi", "-Oz", true},
		{"chunks_bcf_csi", multi, "-W", "-Ob", false},
	} {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertWriteIndexParity(t, bin, tc.wantTBI, func(dir string) []string {
				return []string{"+scatter", "-o", dir, tc.oType, "-n", "2", tc.wflag, tc.fix}
			})
		})
	}
}

// TestNativePluginMendelian2WriteIndex checks mendelian2's -W indexing for a
// VCF-emitting mode, plus that the text-counts mode (-m c) writes no index and
// that -W to stdout errors.
func TestNativePluginMendelian2WriteIndex(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "split_multi.vcf")
	base := []string{"+mendelian2", "-p", "S1,S2,S3", "-m", "a"}
	for _, tc := range []struct {
		name    string
		fix     string
		oFile   string
		wflag   string
		oType   string
		wantTBI bool
	}{
		{"vcfgz_csi", gt, "out.vcf.gz", "-W", "-Oz", false},
		{"vcfgz_tbi", gt, "out.vcf.gz", "-W=tbi", "-Oz", true},
		{"bcf_csi", multi, "out.bcf", "-W", "-Ob", false},
	} {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertWriteIndexParity(t, bin, tc.wantTBI, func(dir string) []string {
				argv := append([]string{}, base...)
				return append(argv, tc.oType, "-o", filepath.Join(dir, tc.oFile), tc.wflag, tc.fix)
			})
		})
	}

	// -m c (text counts) emits no VCF, so -W produces no index on either binary.
	t.Run("counts_mode_no_index", func(t *testing.T) {
		upDir := t.TempDir()
		ourDir := t.TempDir()
		_ = runUpstreamPluginCmd(t, bin, []string{"+mendelian2", "-p", "S1,S2,S3", "-m", "c", "-Oz", "-o", filepath.Join(upDir, "o.txt"), "-W", gt})
		_ = runOursPluginCmd(t, []string{"+mendelian2", "-p", "S1,S2,S3", "-m", "c", "-Oz", "-o", filepath.Join(ourDir, "o.txt"), "-W", gt})
		_, upIdx := dataAndIndexFiles(t, upDir)
		_, ourIdx := dataAndIndexFiles(t, ourDir)
		if len(upIdx) != 0 || len(ourIdx) != 0 {
			t.Fatalf("counts mode should produce no index: upstream %v ours %v", upIdx, ourIdx)
		}
	})

	t.Run("stdout_errors", func(t *testing.T) {
		assertPluginWriteIndexErrors(t, bin, []string{"+mendelian2", "-p", "S1,S2,S3", "-m", "a", "-Oz", "-W", gt})
	})
}

// assertPluginWriteIndexErrors runs the given argv on BOTH binaries and asserts
// both exit non-zero (the non-indexable-output / stdout-index error cases).
func assertPluginWriteIndexErrors(t *testing.T, bin string, argv []string) {
	t.Helper()
	run := func(b string) error {
		cmd := exec.Command(b, argv...)
		cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
		cmd.Stdout = io.Discard
		cmd.Stderr = io.Discard
		return cmd.Run()
	}
	if run(bin) == nil {
		t.Fatalf("upstream unexpectedly succeeded for %v", argv)
	}
	if run(ourBinPath) == nil {
		t.Fatalf("ours unexpectedly succeeded for %v", argv)
	}
}
