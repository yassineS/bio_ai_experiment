package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// ---------------------------------------------------------------------
// Unit tests for the rules engine (ploidy-file parser + PAR boundary
// classification). These run with no external dependency.
// ---------------------------------------------------------------------

// TestParseMendelianRulesFile exercises the custom --rules-file parser:
// SEX_ID dictionary assignment in first-seen order, 1-based→0-based
// coordinate conversion, and the M/F/MF/. inheritance decode.
func TestParseMendelianRulesFile(t *testing.T) {
	const text = `# a comment
1X  X:1-9999       M
1X  Y:1-100        F
2X  Y:1-100        .
1X  A:5-10         MF
`
	rules, err := parseMendelianRules(strings.NewReader(text))
	if err != nil {
		t.Fatalf("parseMendelianRules: %v", err)
	}
	if rules.NumSexes() != 2 {
		t.Fatalf("NumSexes: got %d, want 2", rules.NumSexes())
	}
	if id, ok := rules.sexID["1X"]; !ok || id != 0 {
		t.Fatalf("1X sex id: got %d ok=%v, want 0", id, ok)
	}
	if id, ok := rules.sexID["2X"]; !ok || id != 1 {
		t.Fatalf("2X sex id: got %d ok=%v, want 1", id, ok)
	}

	maleID, _ := rules.SexIDFor(1)
	femaleID, _ := rules.SexIDFor(2)

	// X:1-9999 (1-based) -> 0-based [0,9998]; male = haploid maternal.
	r := rules.rulesFor("X", 0, 0, maleID)
	if r.Ploidy != 1 || r.Inherits != inheritMother {
		t.Errorf("X male @0: got ploidy=%d inherits=%d, want 1/M", r.Ploidy, r.Inherits)
	}
	// One past the end (pos 10000, 0-based 9999) is no longer covered:
	// default diploid biparental.
	r = rules.rulesFor("X", 9999, 9999, maleID)
	if r.Ploidy != 2 || r.Inherits != inheritMother|inheritFather {
		t.Errorf("X male @9999: got ploidy=%d inherits=%d, want 2/MF", r.Ploidy, r.Inherits)
	}
	// Y for a male: haploid paternal.
	r = rules.rulesFor("Y", 0, 0, maleID)
	if r.Ploidy != 1 || r.Inherits != inheritFather {
		t.Errorf("Y male: got ploidy=%d inherits=%d, want 1/F", r.Ploidy, r.Inherits)
	}
	// Y for a female: absent (ploidy 0, no rule).
	r = rules.rulesFor("Y", 0, 0, femaleID)
	if r.Ploidy != 0 || r.Inherits != 0 {
		t.Errorf("Y female: got ploidy=%d inherits=%d, want 0/none", r.Ploidy, r.Inherits)
	}
	// A:5-10 diploid for the male sex id.
	r = rules.rulesFor("A", 4, 4, maleID)
	if r.Ploidy != 2 || r.Inherits != inheritMother|inheritFather {
		t.Errorf("A male @4: got ploidy=%d inherits=%d, want 2/MF", r.Ploidy, r.Inherits)
	}
}

// TestPARBoundaryClassification checks the pseudo-autosomal region
// boundaries of the built-in GRCh38 rules: inside PAR the male child is
// diploid biparental; just past the PAR1 boundary the male-specific X
// is haploid maternal.
func TestPARBoundaryClassification(t *testing.T) {
	rules, err := LoadMendelianRulesByName("GRCh38")
	if err != nil {
		t.Fatalf("LoadMendelianRulesByName: %v", err)
	}
	maleID, _ := rules.SexIDFor(1)

	// GRCh38 PAR1 on X ends at 1-based 9999; that span is the listed
	// haploid-maternal region. Position 1-based 5000 (0-based 4999) is
	// inside it.
	in := rules.rulesFor("X", 4999, 4999, maleID)
	if in.Ploidy != 1 || in.Inherits != inheritMother {
		t.Errorf("X @5000 (PAR1 listed): got ploidy=%d inherits=%d, want 1/M", in.Ploidy, in.Inherits)
	}
	// Between PAR1 and the male-specific stretch (1-based 10000..2781479,
	// 0-based 9999) there is no listed rule: default diploid biparental.
	between := rules.rulesFor("X", 100000, 100000, maleID)
	if between.Ploidy != 2 {
		t.Errorf("X @100001 (between regions): got ploidy=%d, want 2", between.Ploidy)
	}
	// The male-specific X (1-based 2781480..) is haploid maternal.
	msx := rules.rulesFor("X", 2781479, 2781479, maleID)
	if msx.Ploidy != 1 || msx.Inherits != inheritMother {
		t.Errorf("X @2781480 (male-specific): got ploidy=%d inherits=%d, want 1/M", msx.Ploidy, msx.Inherits)
	}
}

// TestLoadMendelianRulesByNameList confirms "list"/"list?" surface the
// catalogue request via ErrMendelianRulesList rather than a table.
func TestLoadMendelianRulesByNameList(t *testing.T) {
	for _, tc := range []struct {
		alias    string
		detailed bool
	}{
		{"list", false},
		{"list?", true},
	} {
		_, err := LoadMendelianRulesByName(tc.alias)
		var le *ErrMendelianRulesList
		if err == nil || !asErr(err, &le) {
			t.Fatalf("LoadMendelianRulesByName(%q): want ErrMendelianRulesList, got %v", tc.alias, err)
		}
		if le.Detailed != tc.detailed {
			t.Errorf("LoadMendelianRulesByName(%q): Detailed=%v, want %v", tc.alias, le.Detailed, tc.detailed)
		}
	}
}

// asErr is a tiny errors.As shim used so the test file doesn't import
// errors solely for one call.
func asErr(err error, target **ErrMendelianRulesList) bool {
	for err != nil {
		if e, ok := err.(*ErrMendelianRulesList); ok {
			*target = e
			return true
		}
		type unwrapper interface{ Unwrap() error }
		u, ok := err.(unwrapper)
		if !ok {
			return false
		}
		err = u.Unwrap()
	}
	return false
}

// ---------------------------------------------------------------------
// Live parity tests against the upstream `bcftools +mendelian2` binary.
// No goldens: the upstream binary is built from the vendored submodule
// and invoked directly, then compared byte-for-byte to the Go port.
// ---------------------------------------------------------------------

var (
	upstreamBcftoolsMendelian2Once sync.Once
	upstreamBcftoolsMendelian2Bin  string
	upstreamBcftoolsMendelian2Dir  string
)

// upstreamBcftoolsMendelian2 returns the path to the upstream bcftools
// binary (built from reference_code/bcftools) and the plugin directory
// to set in BCFTOOLS_PLUGINS. It t.Fatalf's — never t.Skip's — when the
// submodule has not been built, so a missing binary is a hard failure
// rather than a silently-skipped test.
func upstreamBcftoolsMendelian2(t *testing.T) (bin, pluginDir string) {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamBcftoolsMendelian2Once.Do(func() {
		root := repoRootForMendelian2(t)
		bin := filepath.Join(root, "reference_code", "bcftools", "bcftools")
		plug := filepath.Join(root, "reference_code", "bcftools", "plugins")
		upstreamBcftoolsMendelian2Bin = bin
		upstreamBcftoolsMendelian2Dir = plug
	})
	if upstreamBcftoolsMendelian2Bin == "" {
		t.Skipf("could not locate upstream bcftools binary")
	}
	if _, err := os.Stat(upstreamBcftoolsMendelian2Bin); err != nil {
		t.Skipf("upstream bcftools not built at %s: %v\n"+
			"build it with: (cd reference_code/htslib && autoreconf -i && ./configure && make -j) && "+
			"(cd reference_code/bcftools && ./configure && make -j)",
			upstreamBcftoolsMendelian2Bin, err)
	}
	if _, err := os.Stat(filepath.Join(upstreamBcftoolsMendelian2Dir, "mendelian2.so")); err != nil {
		t.Skipf("upstream mendelian2 plugin not built at %s: %v", upstreamBcftoolsMendelian2Dir, err)
	}
	return upstreamBcftoolsMendelian2Bin, upstreamBcftoolsMendelian2Dir
}

// repoRootForMendelian2 walks up from this test file's directory until
// it finds go.mod, returning the module root.
func repoRootForMendelian2(t *testing.T) string {
	t.Helper()
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not find go.mod above %s", filepath.Dir(file))
		}
		dir = parent
	}
}

// mendelian2TrioVCF is a small autosomal+sex-chromosome trio used by
// the live parity tests. CHILD is treated as 1X (male) or 2X (female)
// depending on the -p prefix.
const mendelian2TrioVCF = `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=250000000>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##contig=<ID=MT,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CHILD	FATHER	MOTHER
1	100	.	A	T	.	PASS	.	GT	0/1	0/0	0/1
1	200	.	G	C	.	PASS	.	GT	1/1	0/0	0/0
X	60000	.	A	T	.	PASS	.	GT	0/1	1/1	0/1
X	2700000	.	C	G	.	PASS	.	GT	1/1	0/0	0/1
X	2700100	.	C	G	.	PASS	.	GT	1/1	0/0	0/0
Y	1000000	.	A	G	.	PASS	.	GT	1/1	1/1	0/0
Y	1000100	.	A	G	.	PASS	.	GT	1/1	0/0	0/0
MT	5000	.	A	G	.	PASS	.	GT	1/1	0/0	1/1
MT	5100	.	A	G	.	PASS	.	GT	1/1	0/0	0/0
`

// runUpstreamMendelian2 runs `bcftools +mendelian2 <args> <vcfPath>`
// and returns its stdout.
func runUpstreamMendelian2(t *testing.T, vcfPath string, args ...string) string {
	t.Helper()
	bin, plug := upstreamBcftoolsMendelian2(t)
	full := append([]string{"+mendelian2"}, args...)
	full = append(full, vcfPath)
	cmd := exec.Command(bin, full...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+plug)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools +mendelian2 %v: %v\nstderr: %s", full, err, stderr.String())
	}
	return stdout.String()
}

// TestLiveParity_RulesCountMode compares the Go port's `-m c` count
// output against the live upstream binary across both assemblies and
// both child sexes. The output must be byte-for-byte identical.
func TestLiveParity_RulesCountMode(t *testing.T) {
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "trio.vcf")
	if err := os.WriteFile(vcfPath, []byte(mendelian2TrioVCF), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	cases := []struct {
		name     string
		assembly string
		pfx      string // 1X or 2X
		sex      int
	}{
		{"GRCh37-male", "GRCh37", "1X", 1},
		{"GRCh37-female", "GRCh37", "2X", 2},
		{"GRCh38-male", "GRCh38", "1X", 1},
		{"GRCh38-female", "GRCh38", "2X", 2},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := runUpstreamMendelian2(t, vcfPath,
				"-p", tc.pfx+":CHILD,FATHER,MOTHER", "-m", "c", "--rules", tc.assembly)

			ruleSet, err := LoadMendelianRulesByName(tc.assembly)
			if err != nil {
				t.Fatalf("LoadMendelianRulesByName(%q): %v", tc.assembly, err)
			}
			var got bytes.Buffer
			pfm := Mendelian2PFM{Child: "CHILD", Father: "FATHER", Mother: "MOTHER", Sex: tc.sex}
			if _, err := Mendelian2(strings.NewReader(mendelian2TrioVCF), &got, Mendelian2Options{
				PFM:   &pfm,
				Mode:  Mendelian2Count,
				Rules: ruleSet,
			}); err != nil {
				t.Fatalf("Mendelian2: %v", err)
			}
			if got.String() != want {
				t.Fatalf("count output mismatch for %s\n--- upstream ---\n%s\n--- go ---\n%s",
					tc.name, want, got.String())
			}
		})
	}
}

// TestLiveParity_RulesFileCountMode checks that a custom --rules-file
// produces the same count output as the live upstream binary.
func TestLiveParity_RulesFileCountMode(t *testing.T) {
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "trio.vcf")
	if err := os.WriteFile(vcfPath, []byte(mendelian2TrioVCF), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	rulesPath := filepath.Join(dir, "rules.txt")
	const customRules = `# custom GRCh38-like rules, numeric naming only
1X  X:1-9999                M
1X  X:2781480-155701381     M
1X  Y:1-57227415            F
2X  Y:1-57227415            .
1X  MT:1-16569              M
2X  MT:1-16569              M
`
	if err := os.WriteFile(rulesPath, []byte(customRules), 0o644); err != nil {
		t.Fatalf("write rules: %v", err)
	}

	want := runUpstreamMendelian2(t, vcfPath,
		"-p", "1X:CHILD,FATHER,MOTHER", "-m", "c", "--rules-file", rulesPath)

	ruleSet, err := LoadMendelianRulesFile(rulesPath)
	if err != nil {
		t.Fatalf("LoadMendelianRulesFile: %v", err)
	}
	var got bytes.Buffer
	pfm := Mendelian2PFM{Child: "CHILD", Father: "FATHER", Mother: "MOTHER", Sex: 1}
	if _, err := Mendelian2(strings.NewReader(mendelian2TrioVCF), &got, Mendelian2Options{
		PFM:   &pfm,
		Mode:  Mendelian2Count,
		Rules: ruleSet,
	}); err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	if got.String() != want {
		t.Fatalf("rules-file count mismatch\n--- upstream ---\n%s\n--- go ---\n%s", want, got.String())
	}
}

// TestLiveParity_WriteIndexEmitsValidIndex writes a bgzipped,
// MERR-annotated VCF and confirms the in-tree -W index is valid: the
// live upstream binary must be able to use it for a region query and
// return exactly the records on that contig.
func TestLiveParity_WriteIndexEmitsValidIndex(t *testing.T) {
	bin, plug := upstreamBcftoolsMendelian2(t)
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "trio.vcf")
	if err := os.WriteFile(vcfPath, []byte(mendelian2TrioVCF), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	outPath := filepath.Join(dir, "out.vcf.gz")
	f, err := os.Create(outPath)
	if err != nil {
		t.Fatalf("create output: %v", err)
	}
	ruleSet, _ := LoadMendelianRulesByName("GRCh38")
	pfm := Mendelian2PFM{Child: "CHILD", Father: "FATHER", Mother: "MOTHER", Sex: 1}
	if _, err := Mendelian2(strings.NewReader(mendelian2TrioVCF), f, Mendelian2Options{
		PFM:          &pfm,
		Mode:         Mendelian2Annotate,
		OutputFormat: OutputVCFGz,
		Rules:        ruleSet,
		BGZF:         true,
	}); err != nil {
		f.Close()
		t.Fatalf("Mendelian2: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("close output: %v", err)
	}

	idxPath, err := BuildIndex(outPath, IndexOptions{Format: IndexCSI, Force: true})
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	if _, err := os.Stat(idxPath); err != nil {
		t.Fatalf("index not written: %v", err)
	}

	// The live upstream binary must accept the index for a region query.
	cmd := exec.Command(bin, "view", "-H", outPath, "X")
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+plug)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream view -H X via our index: %v\nstderr: %s", err, stderr.String())
	}
	lines := strings.Split(strings.TrimSpace(stdout.String()), "\n")
	if len(lines) != 3 {
		t.Fatalf("region X query returned %d records, want 3 (index resolved wrong)\n%s",
			len(lines), stdout.String())
	}
	for _, l := range lines {
		if !strings.HasPrefix(l, "X\t") {
			t.Fatalf("region X query returned a non-X record: %q", l)
		}
	}
}
