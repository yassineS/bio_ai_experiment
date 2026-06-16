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
// The run()-style multi-output plugins split and scatter are now driven with the
// SAME upstream form on BOTH binaries (`+name -o DIR <plugin opts> FILE`): OUR
// host forwards -o/-O and every other option to the plugin, which parses them
// itself, so the earlier `--`-rewritten asymmetry is gone. variantkey-hex is a
// generic init/process plugin upstream (its run() is an init() reading argv[1]),
// so it keeps the `+name FILE -- DIR/` form on both sides. The produced files and
// their contents are identical regardless of the surface argv.

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

// TestNativePluginSplitVepFilter checks the -i/-e filter expressions, which
// upstream evaluates against the expanded per-transcript CSQ subfields it
// registers as INFO tags on the output header (filter_init on hdr_out). It
// covers regex (~), string equality, a numeric subfield, the -d per-transcript
// vs collapsed semantics, the auto-registration of a filter-only subfield (no
// -c), the -f text mode, and combination with -s/-c/-p.
func TestNativePluginSplitVepFilter(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "vep_csq.vcf")
	cases := [][]string{
		// Regex / string equality on a CSQ subfield, include and exclude.
		{"-c", "Consequence", "-i", `Consequence~"missense"`},
		{"-c", "IMPACT", "-e", `IMPACT="LOW"`},
		// Numeric subfield, with and without an explicit -c (auto-registration).
		{"-c", "gnomAD_AF", "-i", "gnomAD_AF<0.01"},
		{"-i", "gnomAD_AF<0.01"},
		{"-i", "gnomAD_AF>0.4"},
		// -d per-transcript: the filter selects which CSQ entries survive.
		{"-c", "Consequence", "-d", "-i", `Consequence~"missense"`},
		{"-c", "gnomAD_AF", "-d", "-i", "gnomAD_AF<0.01"},
		{"-c", "Consequence", "-d", "-e", `IMPACT="MODIFIER"`},
		// Array-OR semantics on the collapsed (non -d) record: IMPACT=MODERATE,MODIFIER.
		{"-c", "IMPACT", "-i", `IMPACT="MODIFIER"`},
		// Combine with severity / multiple columns / prefix.
		{"-c", "Consequence,IMPACT", "-s", "worst", "-i", `IMPACT="HIGH"`},
		{"-c", "Consequence,IMPACT", "-s", "worst", "-p", "vep", "-i", `vepIMPACT="HIGH"`},
		{"-c", "Consequence,IMPACT", "-e", "gnomAD_AF>0.4"},
		// Conjunction across two subfields.
		{"-i", `Consequence~"missense"&&gnomAD_AF<0.01`},
		// drop-sites (-x) / keep-sites (-X) interplay with the filter.
		{"-c", "Consequence", "-i", `Consequence~"missense"`, "-x"},
		{"-c", "Consequence", "-i", `Consequence~"missense"`, "-X"},
		{"-c", "IMPACT", "-e", `IMPACT="LOW"`, "-x"},
		// Text (-f) mode with the filter.
		{"-f", `%CHROM:%POS %Consequence %gnomAD_AF\n`, "-i", "gnomAD_AF<0.01"},
		{"-f", `%CHROM:%POS %Consequence\n`, "-d", "-i", `Consequence~"missense"`},
		{"-f", `%CHROM:%POS %Consequence %IMPACT\n`, "-e", `IMPACT="LOW"`},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			upArgv := append(append([]string{"+split-vep"}, args...), fix)
			ourArgv := append(append([]string{"+split-vep"}, args...), fix)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginSplitVepSelectors checks the transcript-selection surface that
// the native port previously rejected: the EXPRESSION selectors (primary => CANONICAL=YES,
// pick => PICK=1, mane => MANE_SELECT!=""), the arbitrary <FIELD><OP><VALUE> forms
// (=, !=, ~, !~), and the PRN qualifier (:worst rewrites the printed Consequence
// to its single worst term). It drives both binaries with a VEP fixture that
// carries CANONICAL/MANE_SELECT/PICK subfields and diffs byte-for-byte.
func TestNativePluginSplitVepSelectors(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "vep_select.vcf")
	cases := [][]string{
		// primary / pick / mane keyword selectors, text and annotate.
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", "primary"},
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", "pick"},
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", "mane"},
		{"-c", "Consequence,SYMBOL", "-s", "primary", "-p", "vep"},
		{"-c", "Consequence,Feature", "-s", "pick"},
		// EXPRESSION selectors with each operator, per-transcript (-d) text.
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", `CANONICAL=YES`, "-d"},
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", `CANONICAL!=YES`, "-d"},
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", `Feature~"ENST[13]"`, "-d"},
		{"-f", `%CHROM:%POS %Consequence %Feature\n`, "-s", `Feature!~"ENST[13]"`, "-d"},
		// EXPRESSION with -p prefix: the FIELD is prefix-resolved.
		{"-c", "Consequence", "-s", `CANONICAL=YES`, "-p", "vep"},
		// EXPRESSION-only with no -c/-f: drop_sites defaults to 1 (keep matching sites).
		{"-s", "mane"},
		{"-s", "pick"},
		// PRN :worst rewrites the worst Consequence term, in -c, -f and -d.
		{"-c", "Consequence", "-s", "all:any:worst"},
		{"-f", `%CHROM:%POS %Consequence\n`, "-s", "all:any:worst", "-d"},
		{"-c", "Consequence,IMPACT", "-s", "worst:any:worst"},
		// PRN :worst combined with a severity filter.
		{"-c", "Consequence", "-s", "all:missense+:worst"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			upArgv := append(append([]string{"+split-vep"}, args...), fix)
			ourArgv := append(append([]string{"+split-vep"}, args...), fix)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginSplitVepGeneList checks the -g/--gene-list restriction and
// prioritisation machinery: restrict mode (only listed-gene transcripts survive),
// prioritise mode (the leading "+" keeps all but moves listed genes first), the
// --gene-list-fields override, and the no-match case. Both binaries read the same
// gene-list file written into a temp dir.
func TestNativePluginSplitVepGeneList(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "vep_select.vcf")
	dir := t.TempDir()
	genes := filepath.Join(dir, "genes.txt")
	if err := os.WriteFile(genes, []byte("BRCA2\nEGFR\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	geneGene := filepath.Join(dir, "genes_gene.txt")
	if err := os.WriteFile(geneGene, []byte("ENSG2\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	geneNone := filepath.Join(dir, "genes_none.txt")
	if err := os.WriteFile(geneNone, []byte("NONEXISTENT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	cases := [][]string{
		// Restrict: only BRCA2/EGFR transcripts survive (per-transcript text).
		{"-f", `%CHROM:%POS %SYMBOL %Consequence\n`, "-g", genes, "-d"},
		// Prioritise: all kept, BRCA2/EGFR moved to the front.
		{"-f", `%CHROM:%POS %SYMBOL %Consequence\n`, "-g", "+" + genes, "-d"},
		// Restrict, annotate (-c): sites kept, annotations limited to listed genes.
		{"-c", "SYMBOL,Consequence", "-g", genes},
		// --gene-list-fields override: match the Gene subfield instead of SYMBOL.
		{"-f", `%CHROM:%POS %Gene %Consequence\n`, "-g", geneGene, "--gene-list-fields", "Gene", "-d"},
		// No gene matches: restrict mode yields empty annotations.
		{"-c", "SYMBOL", "-g", geneNone},
		{"-f", `%CHROM:%POS %SYMBOL\n`, "-g", geneNone, "-d"},
		// Gene list combined with worst-transcript selection.
		{"-c", "SYMBOL,Consequence", "-s", "worst", "-g", genes},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			upArgv := append(append([]string{"+split-vep"}, args...), fix)
			ourArgv := append(append([]string{"+split-vep"}, args...), fix)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginSplitVepSeverity checks the -S/--severity file override: a
// custom scale re-orders the worst-transcript selection and the severity-range
// filter (term[+|-]). Both binaries read the same scale file.
func TestNativePluginSplitVepSeverity(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "vep_select.vcf")
	dir := t.TempDir()
	// Custom scale making synonymous the most severe, missense the least.
	scale := filepath.Join(dir, "severity.txt")
	if err := os.WriteFile(scale, []byte("# custom\nmissense\nstop_gained\nsynonymous\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	cases := [][]string{
		{"-c", "Consequence", "-s", "worst", "-S", scale},
		{"-f", `%CHROM:%POS %Consequence\n`, "-s", "worst", "-S", scale, "-d"},
		{"-f", `%CHROM:%POS %Consequence\n`, "-s", ":synonymous+", "-S", scale, "-d"},
		{"-c", "Consequence", "-s", ":missense-", "-S", scale},
		// PRN :worst rewrite under the custom scale.
		{"-c", "Consequence", "-s", "all:any:worst", "-S", scale},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			upArgv := append(append([]string{"+split-vep"}, args...), fix)
			ourArgv := append(append([]string{"+split-vep"}, args...), fix)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginSplitVepColumnsTypes checks the --columns-types FILE override:
// the regex-matched type table replaces the built-in presets and drives both the
// emitted ##INFO header Type and the numeric re-parsing of the column values.
// Both binaries read the same types file.
func TestNativePluginSplitVepColumnsTypes(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fix := parityFixture(t, "vep_select.vcf")
	dir := t.TempDir()
	// Override DISTANCE to Float and gnomAD_AF to Integer (both flip from the
	// presets), plus a regex rule.
	ct := filepath.Join(dir, "ctypes.txt")
	if err := os.WriteFile(ct, []byte("# custom\nDISTANCE Float\ngnomAD_AF Integer\n.*_AF Integer\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	cases := [][]string{
		{"-c", "DISTANCE,gnomAD_AF", "--columns-types", ct},
		{"-c", "DISTANCE", "--columns-types", ct, "-f", `%CHROM:%POS %DISTANCE\n`},
		{"-c", "gnomAD_AF", "--columns-types", ct, "-i", "gnomAD_AF<1"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
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
		// Per-output -i/-e filter: upstream compiles the expression against each
		// subset header and tests the already-subsetted record, so a FORMAT
		// expression sees only that file's samples.
		{"filter_site_include", gt, []string{"-i", "QUAL>10"}},
		{"filter_site_exclude", gt, []string{"-e", "QUAL<30"}},
		{"filter_fmt_include", gt, []string{"-i", `GT="het"`}},
		{"filter_fmt_dp", gt, []string{"-i", "FMT/DP>10"}},
		{"filter_groups", gt, []string{"-G", groups, "-i", `GT="het"`}},
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
					// Same upstream run()-style form for OUR port: the host
					// forwards -o/-O and split parses them itself.
					argv := append([]string{"+split", "-o", dir}, append(append([]string{}, tc.args...), tc.fix)...)
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
		// -i/-e are accepted but applied to nothing, exactly as upstream scatter.c
		// (which never calls filter_init/filter_test). The output must be identical
		// to the same run without the filter.
		{"filter_noop_include", gt, []string{"-n", "2", "-i", "QUAL>10"}},
		{"filter_noop_exclude", gt, []string{"-n", "2", "-e", "QUAL<30"}},
		{"filter_noop_fmt", gt, []string{"-n", "2", "-i", `GT="het"`}},
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
					// Same upstream run()-style form for OUR port: the host
					// forwards -o/-O/-n/-s and scatter parses them itself.
					argv := append([]string{"+scatter", "-o", dir}, append(append([]string{}, tc.args...), tc.fix)...)
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
		// split-vep: the -i/-e filter expressions, the gene list, the
		// EXPRESSION/primary/pick/mane transcript selectors, the PRN :worst
		// qualifier, and the -S/--columns-types file overrides are all now
		// supported and parity-checked in the dedicated TestNativePluginSplitVep*
		// tests. The cases below remain genuine clean errors.
		//
		// A -g/--gene-list file that cannot be read is a clean error.
		{"split-vep", []string{"-c", "Consequence", "-g", "/tmp/no_such_gene_list_file"}, csq},
		// A -i/-e expression that references a tag which is neither a CSQ subfield
		// nor declared in the header is a clean error, matching upstream's
		// "the tag ... is not defined in the VCF header or in INFO/CSQ".
		{"split-vep", []string{"-i", "NoSuchField>1"}, csq},
		// An EXPRESSION/mane selector whose FIELD is absent from INFO/CSQ errors
		// (this fixture has no CANONICAL/MANE_SELECT subfields).
		{"split-vep", []string{"-c", "Consequence", "-s", "mane"}, csq},
		{"split-vep", []string{"-c", "Consequence", "-s", "CANONICAL=YES"}, csq},
		// -S - and --columns-types - print the default table to stderr and exit
		// non-zero, exactly as upstream's pre-init checks do.
		{"split-vep", []string{"-S", "-"}, csq},
		{"split-vep", []string{"--columns-types", "-"}, csq},
		// An unknown consequence term in the severity filter is a clean error.
		{"split-vep", []string{"-c", "Consequence", "-s", ":nosuchterm"}, csq},
		// A PRN qualifier other than all/worst is a clean error.
		{"split-vep", []string{"-c", "Consequence", "-s", "all:any:bogus"}, csq},
		// split: -W/--write-index is now supported and parity-checked in
		// TestNativePluginSplitWriteIndex (the per-output -i/-e filter is
		// parity-checked in TestNativePluginSplit; -r/-R/-t/-T region/target
		// selection in TestNativePluginRegionTarget).
		{"split", nil, gt}, // missing -o
		// scatter: region pre-selection and missing -n/-s remain unsupported. (-i/-e
		// are accepted but applied to nothing, matching upstream's dead option; the
		// "only one of -i or -e" guard and the no-op behaviour are covered by
		// TestNativePluginScatter. Giving both -i and -e is still an error.)
		{"scatter", []string{"-o", "/tmp/x", "-n", "1", "-i", "QUAL>1", "-e", "QUAL<1"}, gt},
		{"scatter", []string{"-o", "/tmp/x"}, gt}, // missing -n/-s
		{"scatter", []string{"-n", "1"}, gt},      // missing -o
		// isecGT: missing the second file. (-r/-R/-t/-T region/target selection is
		// now supported and parity-checked in TestNativePluginRegionTarget.)
		{"isecGT", nil, gt},
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
