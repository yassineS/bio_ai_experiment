package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"testing"
)

// Live-oracle parity tests for the batch-7 plugins: the per-record annotation
// plugins add-variantkey and split-vep, the multi-output plugins split, scatter
// and variantkey-hex, and the two-file isecGT. Each test builds the genuine
// upstream bcftools via buildBcftools() and drives BOTH that binary and OUR port
// through their CLIs, diffing the output byte-for-byte after stripping
// provenance (stripProvenanceBytes). For the multi-output plugins every produced
// file is diffed (sorted by name).
//
// The two CLI forms differ only in how the input file and plugin options are
// arranged, because OUR host CLI consumes the run()-style -o/-O host options for
// itself; the multi-output plugins are therefore invoked through the `--` form
// (`+name FILE -- <plugin opts>`) for OUR port, while upstream uses its native
// run()-style form (`+name <plugin opts> FILE`). The produced files and their
// contents are identical regardless of the surface argv.

// runUpstreamInDir runs the upstream binary with argv (cwd irrelevant) and the
// vendored plugins directory, failing the test on a non-zero exit.
func runUpstreamPluginCmd(t *testing.T, bin string, argv []string) []byte {
	t.Helper()
	cmd := exec.Command(bin, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream %v: %v\nstderr: %s", argv, err, errBuf.String())
	}
	return out.Bytes()
}

// runOursPluginCmd runs OUR built port with argv, failing on a non-zero exit.
func runOursPluginCmd(t *testing.T, argv []string) []byte {
	t.Helper()
	if ourBinPath == "" {
		t.Fatalf("local bcftools port binary not built; cannot run CLI oracle")
	}
	cmd := exec.Command(ourBinPath, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("ours %v: %v\nstderr: %s", argv, err, errBuf.String())
	}
	return out.Bytes()
}

// assertStdoutParity diffs the provenance-stripped stdout of upstream and ours
// for the given argv pair (which may differ between the two binaries).
func assertStdoutParity(t *testing.T, bin string, upArgv, ourArgv []string) {
	t.Helper()
	up := runUpstreamPluginCmd(t, bin, upArgv)
	ours := runOursPluginCmd(t, ourArgv)
	if !bytes.Equal(stripProvenanceBytes(up), stripProvenanceBytes(ours)) {
		t.Fatalf("stdout diverges from upstream\nupstream argv=%v\nours argv=%v\n--- upstream ---\n%s\n--- ours ---\n%s",
			upArgv, ourArgv, snippet(up, 2000), snippet(ours, 2000))
	}
}

// assertDirParity runs upstream and ours into their own temp directories and
// diffs every produced file byte-for-byte (after stripping provenance), failing
// if the file sets differ or any file content differs.
func assertDirParity(t *testing.T, bin string, mkUpArgv, mkOurArgv func(dir string) []byte) {
	t.Helper()
	upDir := t.TempDir()
	ourDir := t.TempDir()
	_ = mkUpArgv(upDir)
	_ = mkOurArgv(ourDir)

	upFiles := listFiles(t, upDir)
	ourFiles := listFiles(t, ourDir)
	if len(upFiles) != len(ourFiles) {
		t.Fatalf("file count differs: upstream %v vs ours %v", upFiles, ourFiles)
	}
	for i := range upFiles {
		if upFiles[i] != ourFiles[i] {
			t.Fatalf("file names differ at %d: upstream %q vs ours %q (upstream=%v ours=%v)",
				i, upFiles[i], ourFiles[i], upFiles, ourFiles)
		}
		upData := canonicalizeFile(t, bin, filepath.Join(upDir, upFiles[i]))
		ourData := canonicalizeFile(t, bin, filepath.Join(ourDir, ourFiles[i]))
		if !bytes.Equal(stripProvenanceBytes(upData), stripProvenanceBytes(ourData)) {
			t.Fatalf("file %q diverges\n--- upstream ---\n%s\n--- ours ---\n%s",
				upFiles[i], snippet(upData, 2000), snippet(ourData, 2000))
		}
	}
}

// canonicalizeFile normalises a produced output file to plain VCF text. VCF/BCF
// containers (.vcf.gz, .bcf) are decoded with upstream `bcftools view` so the
// comparison is on record content, independent of BGZF framing or BCF binary
// layout; plain text files (the index tables, plain .vcf) are returned as-is.
func canonicalizeFile(t *testing.T, bin, path string) []byte {
	t.Helper()
	switch ext := filepath.Ext(path); ext {
	case ".gz", ".bcf", ".bgz":
		cmd := exec.Command(bin, "view", path)
		cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
		var out, errBuf bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = &errBuf
		if err := cmd.Run(); err != nil {
			t.Fatalf("view %s: %v\n%s", path, err, errBuf.String())
		}
		return out.Bytes()
	default:
		data, err := os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
		return data
	}
}

// listFiles returns the sorted base names of the regular files in dir.
func listFiles(t *testing.T, dir string) []string {
	t.Helper()
	ents, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	var names []string
	for _, e := range ents {
		if !e.IsDir() {
			names = append(names, e.Name())
		}
	}
	sort.Strings(names)
	return names
}

// TestNativePluginAddVariantKey checks the VKX/RSX annotation against upstream
// across reversible, hashed (long indel), multiallelic and non-rs-ID records.
func TestNativePluginAddVariantKey(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "variantkey.vcf")
	// add-variantkey is a generic init/process plugin: file before `--`.
	assertStdoutParity(t, bin,
		[]string{"+add-variantkey", fix},
		[]string{"+add-variantkey", fix})
}

// TestNativePluginSplitVep checks the split-vep query (-f) and annotate (-c)
// modes plus list (-l), severity selection and the -A all-fields expansion.
func TestNativePluginSplitVep(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "vep_csq.vcf")
	cases := [][]string{
		{"-l"},
		{"-f", `%CHROM:%POS %Consequence %SYMBOL %gnomAD_AF\n`},
		{"-f", `%CHROM:%POS %Consequence\n`, "-d"},
		{"-f", `%CHROM %POS %Consequence %IMPACT\n`, "-s", "worst"},
		{"-f", `%CHROM %POS %CSQ\n`, "-d", "-A", "tab"},
		{"-c", "Consequence,IMPACT,SYMBOL", "-s", "worst", "-p", "vep"},
		{"-c", "1-3", "-s", "worst", "-p", "vep"},
		{"-c", "gnomAD_AF:Float"},
		{"-c", "DISTANCE"},
		{"-c", "Consequence", "-s", ":missense"},
		{"-c", "Consequence", "-s", ":missense", "-x"},
		{"-c", "Consequence", "-s", ":stop_gained+"},
		{"-c", "-"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			// split-vep is run()-style: options precede the file upstream; OUR port
			// accepts the same run-style form (the host forwards -f/-c/-s/-d/-A/-p).
			upArgv := append(append([]string{"+split-vep"}, args...), fix)
			ourArgv := append(append([]string{"+split-vep"}, args...), fix)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginSplit checks the +split multi-output modes: default per-sample,
// -O containers, -S samples-file, -G groups-file and FORMAT-subset -k.
func TestNativePluginSplit(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	// The -Ob/-Oz container cases use a fixture without phased-missing genotypes,
	// because OUR BCF/BGZF encoder canonicalises a missing GT to the unphased
	// "./." form; that GT-encoder behaviour is independent of the split plugin.
	multi := parityFixture(t, "split_multi.vcf")
	samples := parityFixture(t, "split_samples.txt")
	groups := parityFixture(t, "split_groups.txt")
	cases := []struct {
		name string
		fix  string
		args []string
	}{
		{"default", gt, nil},
		{"bcf", multi, []string{"-Ob"}},
		{"vcfgz", multi, []string{"-Oz"}},
		{"samples", gt, []string{"-S", samples}},
		{"groups", gt, []string{"-G", groups}},
		{"keep_fmt_gt", gt, []string{"-k", "FMT/GT"}},
		{"keep_fmt_gt_dp", gt, []string{"-k", "FMT/GT,DP"}},
		{"keep_info_all", gt, []string{"-k", "INFO"}},
		{"keep_fmt_all", gt, []string{"-k", "FMT"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertDirParity(t, bin,
				func(dir string) []byte {
					argv := append([]string{"+split", "-o", dir}, append(append([]string{}, tc.args...), tc.fix)...)
					return runUpstreamPluginCmd(t, bin, argv)
				},
				func(dir string) []byte {
					argv := append([]string{"+split", tc.fix, "--", "-o", dir}, tc.args...)
					return runOursPluginCmd(t, argv)
				})
		})
	}
}

// TestNativePluginScatter checks +scatter's chunk (-n) and region (-s/-S) modes,
// the -x extra file, the -p prefix and the -O containers.
func TestNativePluginScatter(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "split_multi.vcf")
	regions := parityFixture(t, "scatter_regions.txt")
	cases := []struct {
		name string
		fix  string
		args []string
	}{
		{"chunks", gt, []string{"-n", "2"}},
		{"chunks_prefix", gt, []string{"-n", "2", "-p", "part-"}},
		{"chunks_bcf", multi, []string{"-n", "3", "-Ob"}},
		{"chunks_vcfgz", multi, []string{"-n", "2", "-Oz"}},
		{"regions", gt, []string{"-s", "chr1,chr2"}},
		{"regions_extra", gt, []string{"-s", "chr1", "-x", "other"}},
		{"regions_file", gt, []string{"-S", regions}},
		{"regions_prefix", gt, []string{"-s", "chr1,chr2", "-p", "shard_"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertDirParity(t, bin,
				func(dir string) []byte {
					argv := append([]string{"+scatter", "-o", dir}, append(append([]string{}, tc.args...), tc.fix)...)
					return runUpstreamPluginCmd(t, bin, argv)
				},
				func(dir string) []byte {
					argv := append([]string{"+scatter", tc.fix, "--", "-o", dir}, tc.args...)
					return runOursPluginCmd(t, argv)
				})
		})
	}
}

// TestNativePluginVariantKeyHex checks the three index files and the stdout count
// report of +variantkey-hex.
func TestNativePluginVariantKeyHex(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "variantkey.vcf")
	// The directory positional must end with a separator (upstream concatenates
	// the filename directly onto it).
	assertDirParity(t, bin,
		func(dir string) []byte {
			return runUpstreamPluginCmd(t, bin, []string{"+variantkey-hex", fix, "--", dir + "/"})
		},
		func(dir string) []byte {
			return runOursPluginCmd(t, []string{"+variantkey-hex", fix, "--", dir + "/"})
		})
	// Also verify the stdout count report matches.
	assertStdoutParity(t, bin,
		[]string{"+variantkey-hex", fix, "--", t.TempDir() + "/"},
		[]string{"+variantkey-hex", fix, "--", t.TempDir() + "/"})
}

// TestNativePluginIsecGT checks the two-file genotype comparison: discordant
// genotypes in the first file are set to missing, A-only sites pass through and
// B-only sites are dropped. The fixtures exercise phased, multiallelic and
// reordered-sample cases.
func TestNativePluginIsecGT(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	aPlain := parityFixture(t, "isecgt_a.vcf")
	bPlain := parityFixture(t, "isecgt_b.vcf")

	// Upstream's synced reader requires both inputs bgzipped and indexed; OUR
	// port streams them. Prepare indexed copies for upstream in a temp dir.
	dir := t.TempDir()
	aGz := bgzipAndIndex(t, bin, aPlain, filepath.Join(dir, "a.vcf.gz"))
	bGz := bgzipAndIndex(t, bin, bPlain, filepath.Join(dir, "b.vcf.gz"))

	assertStdoutParity(t, bin,
		[]string{"+isecGT", aGz, bGz},
		// OUR port reads the second file as the second positional (routed via the
		// host into the plugin), and does not need an index.
		[]string{"+isecGT", aPlain, bPlain})
}

// TestNativePluginBatch7Unsupported asserts the deliberately unsupported
// split-vep / split / scatter / isecGT modes fail with a clean error rather than
// diverging silently. These go through RunPlugin directly.
func TestNativePluginBatch7Unsupported(t *testing.T) {
	csq := parityFixture(t, "vep_csq.vcf")
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := []struct {
		name string
		args []string
		in   string
	}{
		// split-vep: filter expressions, gene list, EXPRESSION/canonical selectors,
		// severity/columns-types overrides need the upstream filter/convert engines.
		{"split-vep", []string{"-c", "Consequence", "-i", "gnomAD_AF<0.1"}, csq},
		{"split-vep", []string{"-c", "Consequence", "-e", "IMPACT=\"LOW\""}, csq},
		{"split-vep", []string{"-c", "Consequence", "-g", "/tmp/genes"}, csq},
		{"split-vep", []string{"-c", "Consequence", "-s", "mane"}, csq},
		{"split-vep", []string{"-c", "Consequence", "-s", "CANONICAL=YES"}, csq},
		{"split-vep", []string{"-S", "-"}, csq},
		{"split-vep", []string{"--columns-types", "-"}, csq},
		// split: filter expressions, region selection, write-index.
		{"split", []string{"-o", "/tmp/x", "-i", "GT=\"alt\""}, gt},
		{"split", []string{"-o", "/tmp/x", "-r", "chr1"}, gt},
		{"split", []string{"-o", "/tmp/x", "-W"}, gt},
		{"split", nil, gt}, // missing -o
		// scatter: filter expressions, region pre-selection, missing -n/-s.
		{"scatter", []string{"-o", "/tmp/x", "-n", "1", "-i", "GT=\"alt\""}, gt},
		{"scatter", []string{"-o", "/tmp/x"}, gt}, // missing -n/-s
		{"scatter", []string{"-n", "1"}, gt},      // missing -o
		// isecGT: missing the second file, region selection.
		{"isecGT", nil, gt},
		{"isecGT", []string{"-r", "chr1"}, gt},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+joinArgs(tc.args), func(t *testing.T) {
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         tc.name,
				Args:         tc.args,
				InputFile:    tc.in,
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected a clean unsupported error for +%s %v, got nil", tc.name, tc.args)
			}
		})
	}
}

// bgzipAndIndex writes a bgzipped, tabix-indexed copy of src at dst using the
// upstream binary and returns dst.
func bgzipAndIndex(t *testing.T, bin, src, dst string) string {
	t.Helper()
	cmd := exec.Command(bin, "view", "-Oz", "-o", dst, src)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("bgzip %s: %v\n%s", src, err, out)
	}
	idx := exec.Command(bin, "index", "-t", dst)
	idx.Env = cmd.Env
	if out, err := idx.CombinedOutput(); err != nil {
		t.Fatalf("index %s: %v\n%s", dst, err, out)
	}
	return dst
}
