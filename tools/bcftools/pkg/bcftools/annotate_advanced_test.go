package bcftools

// Live-parity and unit tests for the advanced `bcftools annotate` flags:
// --set-id macros, --merge-logic, --min-overlap, --pair-logic,
// --single-overlaps and --rename-annots.
//
// The parity tests build the vendored upstream bcftools (and its htslib)
// once via upstreamBcftoolsAnnotate and compare its stdout byte-for-byte with
// our Go port over fixtures crafted to exercise each flag. The only
// normalisation applied is dropping the upstream-only
// `##bcftools_annotate{Version,Command}=` provenance lines (environmental
// noise our version-less port never emits).
//
// These tests never t.Skip: when the upstream toolchain cannot be built they
// t.Fatalf with the remediation (init the submodules + a C toolchain), so a
// silent gap can never masquerade as a pass.

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sync"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// upstreamAnnotateOnce guards the one-time upstream build.
var upstreamAnnotateOnce sync.Once

// upstreamAnnotateBin is the path to the freshly built upstream bcftools, or
// "" if the build failed (with upstreamAnnotateErr describing why).
var (
	upstreamAnnotateBin string
	upstreamAnnotateErr error
)

// repoRootForAnnotate walks up from the test's working directory to the repo
// root (the directory that holds go.mod).
func repoRootForAnnotate(t *testing.T) string {
	t.Helper()
	dir, err := os.Getwd()
	if err != nil {
		t.Fatalf("getwd: %v", err)
	}
	for i := 0; i < 12; i++ {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		dir = filepath.Dir(dir)
	}
	t.Fatalf("could not locate repo root (go.mod) from working directory")
	return ""
}

// upstreamBcftoolsAnnotate builds htslib + bcftools from reference_code once
// and returns the path to the bcftools binary. It t.Fatalf's (never skips) on
// any failure so a missing toolchain surfaces as a hard error.
func upstreamBcftoolsAnnotate(t *testing.T) string {
	t.Helper()
	root := repoRootForAnnotate(t)
	upstreamAnnotateOnce.Do(func() {
		htslib := filepath.Join(root, "reference_code", "htslib")
		bcftools := filepath.Join(root, "reference_code", "bcftools")
		if _, err := os.Stat(filepath.Join(bcftools, "vcfannotate.c")); err != nil {
			upstreamAnnotateErr = errf("reference_code/bcftools not initialised: %v; run `git submodule update --init --recursive reference_code/htslib reference_code/bcftools`", err)
			return
		}
		bin := filepath.Join(bcftools, "bcftools")
		// Reuse a previously built binary when present (the build is slow).
		if fi, err := os.Stat(bin); err == nil && fi.Mode()&0o111 != 0 {
			upstreamAnnotateBin = bin
			return
		}
		// htslib: autoreconf + configure if config.mk is missing.
		if _, err := os.Stat(filepath.Join(htslib, "config.mk")); err != nil {
			if out, err := runIn(htslib, "autoreconf", "-i"); err != nil {
				upstreamAnnotateErr = errf("htslib autoreconf failed: %v\n%s", err, out)
				return
			}
			if out, err := runIn(htslib, "./configure", "--disable-bz2", "--disable-lzma", "--disable-libcurl"); err != nil {
				upstreamAnnotateErr = errf("htslib configure failed: %v\n%s", err, out)
				return
			}
		}
		if out, err := runIn(htslib, "make", "-j4"); err != nil {
			upstreamAnnotateErr = errf("htslib make failed: %v\n%s", err, out)
			return
		}
		// bcftools: configure if config.mk missing.
		if _, err := os.Stat(filepath.Join(bcftools, "config.mk")); err != nil {
			_, _ = runIn(bcftools, "autoheader")
			if out, err := runIn(bcftools, "autoconf"); err != nil {
				upstreamAnnotateErr = errf("bcftools autoconf failed: %v\n%s", err, out)
				return
			}
			if out, err := runIn(bcftools, "./configure", "--with-htslib=../htslib"); err != nil {
				upstreamAnnotateErr = errf("bcftools configure failed: %v\n%s", err, out)
				return
			}
		}
		if out, err := runIn(bcftools, "make", "-j4"); err != nil {
			upstreamAnnotateErr = errf("bcftools make failed: %v\n%s", err, out)
			return
		}
		upstreamAnnotateBin = bin
	})
	if upstreamAnnotateErr != nil {
		t.Fatalf("upstream bcftools unavailable: %v", upstreamAnnotateErr)
	}
	if upstreamAnnotateBin == "" {
		t.Fatalf("upstream bcftools binary not produced")
	}
	return upstreamAnnotateBin
}

// runIn runs a command in dir and returns its combined output.
func runIn(dir, name string, args ...string) ([]byte, error) {
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	return cmd.CombinedOutput()
}

// errf wraps fmt.Errorf for the build-failure messages.
func errf(format string, a ...interface{}) error {
	return fmt.Errorf(format, a...)
}

// annotateVersionRE drops the upstream provenance lines.
var annotateVersionRE = regexp.MustCompile(`(?m)^##bcftools_annotate(Version|Command)=.*\n`)

func stripAnnotateProvenance(b []byte) []byte {
	return annotateVersionRE.ReplaceAll(b, nil)
}

// upstreamHtslibBin returns the directory holding the freshly built bgzip and
// tabix helpers, so fixtures can be bgzipped + indexed exactly as upstream
// expects.
func upstreamHtslibBin(t *testing.T) string {
	t.Helper()
	upstreamBcftoolsAnnotate(t) // ensure htslib is built
	return filepath.Join(repoRootForAnnotate(t), "reference_code", "htslib")
}

// writeFile writes content to a file under dir and returns its path.
func writeAnnFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

// bgzipIndex bgzips src in place (-> src.gz) and tabixes it with the given
// args, using the freshly built htslib helpers.
func bgzipIndex(t *testing.T, htslibDir, src string, tabixArgs ...string) string {
	t.Helper()
	bgzip := filepath.Join(htslibDir, "bgzip")
	tabix := filepath.Join(htslibDir, "tabix")
	if out, err := exec.Command(bgzip, "-f", src).CombinedOutput(); err != nil {
		t.Fatalf("bgzip %s: %v\n%s", src, err, out)
	}
	gz := src + ".gz"
	args := append(append([]string{}, tabixArgs...), gz)
	if out, err := exec.Command(tabix, args...).CombinedOutput(); err != nil {
		t.Fatalf("tabix %v: %v\n%s", args, err, out)
	}
	return gz
}

// runUpstreamAnnotate runs `bcftools annotate <args>` and returns stdout with
// the provenance lines stripped.
func runUpstreamAnnotate(t *testing.T, args ...string) []byte {
	t.Helper()
	bin := upstreamBcftoolsAnnotate(t)
	cmd := exec.Command(bin, append([]string{"annotate"}, args...)...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream annotate %v failed: %v\n%s", args, err, stderr.String())
	}
	return stripAnnotateProvenance(stdout.Bytes())
}

// runGoAnnotate runs the Go port via AnnotateFile and returns stdout.
func runGoAnnotate(t *testing.T, inputPath string, opts AnnotateOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := AnnotateFile(inputPath, &out, opts); err != nil {
		t.Fatalf("Go AnnotateFile(%s): %v", inputPath, err)
	}
	return out.Bytes()
}

// ---------------------------------------------------------------------------
// Live parity tests
// ---------------------------------------------------------------------------

const setIDInputVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##FILTER=<ID=q10,Description="q10">
##INFO=<ID=AC,Number=A,Type=Integer,Description="ac">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	50	PASS	AC=3
1	200	rs9	C	T,TT	.	q10	AC=1,2
1	300	.	G	.	10	.	.
1	400	.	AT	A	.	.	.
1	500	.	A	AT	.	.	.
`

func TestAnnotateParity_SetIDMacros(t *testing.T) {
	dir := t.TempDir()
	in := writeAnnFile(t, dir, "in.vcf", setIDInputVCF)
	formats := []string{
		`+%CHROM\_%POS`,
		`%CHROM:%POS:%REF:%FIRST_ALT`,
		`%CHROM-%INFO/AC`,
		`%ALT|%QUAL|%FILTER|%TYPE|%END`,
		`%CHROM\_%POS0\_%END0`,
		`set_%TYPE`,
	}
	for _, f := range formats {
		f := f
		t.Run(f, func(t *testing.T) {
			want := runUpstreamAnnotate(t, "--set-id", f, in)
			got := runGoAnnotate(t, in, AnnotateOptions{SetID: f})
			if !bytes.Equal(got, want) {
				t.Fatalf("--set-id %q mismatch.\n--- want ---\n%s\n--- got ---\n%s", f, want, got)
			}
		})
	}
}

const rangeMainVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	130	.	A	G	.	.	.
1	160	.	C	T	.	.	.
1	550	.	C	T	.	.	.
`

const rangeAnnTAB = `1	100	150	geneA	5
1	120	200	geneB	8
1	500	600	geneC	2
`

const rangeHdrLines = `##INFO=<ID=GENE,Number=.,Type=String,Description="g">
##INFO=<ID=SCORE,Number=.,Type=Integer,Description="s">
`

// setupRangeFixture writes the range main VCF, the bgzipped+tabixed
// annotation table and the header-lines file, returning their paths.
func setupRangeFixture(t *testing.T) (mainVCF, annGz, hdr string) {
	t.Helper()
	htslibDir := upstreamHtslibBin(t)
	dir := t.TempDir()
	mainVCF = writeAnnFile(t, dir, "main.vcf", rangeMainVCF)
	hdr = writeAnnFile(t, dir, "hdr.txt", rangeHdrLines)
	tab := writeAnnFile(t, dir, "ann.tab", rangeAnnTAB)
	annGz = bgzipIndex(t, htslibDir, tab, "-s1", "-b2", "-e3")
	return mainVCF, annGz, hdr
}

func TestAnnotateParity_MergeLogic(t *testing.T) {
	mainVCF, annGz, hdr := setupRangeFixture(t)
	cases := []struct {
		name       string
		mergeFlags []string // pairs of (--merge-logic, VALUE)
		mergeOpt   string
		columns    string
	}{
		{"first_default", nil, "", "CHROM,BEG,END,GENE,SCORE"},
		{"unique_sum", []string{"--merge-logic", "GENE:unique", "--merge-logic", "SCORE:sum"}, "GENE:unique,SCORE:sum", "CHROM,BEG,END,GENE,SCORE"},
		{"append", []string{"--merge-logic", "GENE:append"}, "GENE:append", "CHROM,BEG,END,GENE"},
		{"min", []string{"--merge-logic", "SCORE:min"}, "SCORE:min", "CHROM,BEG,END,-,SCORE"},
		{"max", []string{"--merge-logic", "SCORE:max"}, "SCORE:max", "CHROM,BEG,END,-,SCORE"},
		{"avg", []string{"--merge-logic", "SCORE:avg"}, "SCORE:avg", "CHROM,BEG,END,-,SCORE"},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			upArgs := []string{"-a", annGz, "-h", hdr, "-c", tc.columns}
			upArgs = append(upArgs, tc.mergeFlags...)
			upArgs = append(upArgs, mainVCF)
			want := runUpstreamAnnotate(t, upArgs...)
			got := runGoAnnotate(t, mainVCF, AnnotateOptions{
				Annotations: annGz,
				HeaderLines: hdr,
				Columns:     tc.columns,
				MergeLogic:  tc.mergeOpt,
			})
			if !bytes.Equal(got, want) {
				t.Fatalf("%s mismatch.\n--- want ---\n%s\n--- got ---\n%s", tc.name, want, got)
			}
		})
	}
}

func TestAnnotateParity_MinOverlap(t *testing.T) {
	mainVCF, annGz, hdr := setupRangeFixture(t)
	for _, mo := range []string{"0:1", "0.5", "0.01:0.5", ":1"} {
		mo := mo
		t.Run(mo, func(t *testing.T) {
			want := runUpstreamAnnotate(t, "-a", annGz, "-h", hdr,
				"-c", "CHROM,BEG,END,GENE", "--merge-logic", "GENE:unique",
				"--min-overlap", mo, mainVCF)
			got := runGoAnnotate(t, mainVCF, AnnotateOptions{
				Annotations: annGz,
				HeaderLines: hdr,
				Columns:     "CHROM,BEG,END,GENE",
				MergeLogic:  "GENE:unique",
				MinOverlap:  mo,
			})
			if !bytes.Equal(got, want) {
				t.Fatalf("--min-overlap %s mismatch.\n--- want ---\n%s\n--- got ---\n%s", mo, want, got)
			}
		})
	}
}

func TestAnnotateParity_SingleOverlaps(t *testing.T) {
	mainVCF, annGz, hdr := setupRangeFixture(t)
	want := runUpstreamAnnotate(t, "-a", annGz, "-h", hdr,
		"-c", "CHROM,BEG,END,GENE", "--single-overlaps", mainVCF)
	got := runGoAnnotate(t, mainVCF, AnnotateOptions{
		Annotations:    annGz,
		HeaderLines:    hdr,
		Columns:        "CHROM,BEG,END,GENE",
		SingleOverlaps: true,
	})
	if !bytes.Equal(got, want) {
		t.Fatalf("--single-overlaps mismatch.\n--- want ---\n%s\n--- got ---\n%s", want, got)
	}
}

const pairAnnVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##INFO=<ID=TAG,Number=1,Type=Integer,Description="t">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	annA	A	G	.	.	TAG=11
1	150	annR	A	.	.	.	TAG=99
1	200	annB	C	T	.	.	TAG=22
1	300	annC	G	GA	.	.	TAG=33
1	400	annD	G	T,A	.	.	TAG=44
`

const pairQueryVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	.	.
1	150	.	A	C	.	.	.
1	200	.	C	A	.	.	.
1	300	.	G	GA	.	.	.
1	400	.	G	A	.	.	.
`

func TestAnnotateParity_PairLogic(t *testing.T) {
	htslibDir := upstreamHtslibBin(t)
	dir := t.TempDir()
	annVCF := writeAnnFile(t, dir, "ann.vcf", pairAnnVCF)
	annGz := bgzipIndex(t, htslibDir, annVCF, "-p", "vcf")
	qVCF := writeAnnFile(t, dir, "q.vcf", pairQueryVCF)
	qGz := bgzipIndex(t, htslibDir, qVCF, "-p", "vcf")

	for _, mode := range []string{"some", "exact", "all", "any", "snps", "indels", "both", "id"} {
		mode := mode
		t.Run(mode, func(t *testing.T) {
			want := runUpstreamAnnotate(t, "-a", annGz, "-c", "INFO/TAG", "--pair-logic", mode, qGz)
			got := runGoAnnotate(t, qGz, AnnotateOptions{
				Annotations: annGz,
				Columns:     "INFO/TAG",
				PairLogic:   mode,
			})
			if !bytes.Equal(got, want) {
				t.Fatalf("--pair-logic %s mismatch.\n--- want ---\n%s\n--- got ---\n%s", mode, want, got)
			}
		})
	}
}

const renameInputVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FILTER=<ID=LowQ,Description="lq">
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	100	.	A	G	.	LowQ	DP=7	GT:GQ	0/1:30
`

const renameMap = `INFO/DP	DEPTH
FORMAT/GQ	GENOQUAL
FILTER/LowQ	LOWQUAL
`

func TestAnnotateParity_RenameAnnots(t *testing.T) {
	dir := t.TempDir()
	in := writeAnnFile(t, dir, "in.vcf", renameInputVCF)
	mapFile := writeAnnFile(t, dir, "map.txt", renameMap)
	want := runUpstreamAnnotate(t, "--rename-annots", mapFile, in)
	got := runGoAnnotate(t, in, AnnotateOptions{RenameAnnots: mapFile})
	if !bytes.Equal(got, want) {
		t.Fatalf("--rename-annots mismatch.\n--- want ---\n%s\n--- got ---\n%s", want, got)
	}
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

// mkVariant builds a vcf.Variant from compact fields for unit tests.
func mkAnnVariant(chrom string, pos int, id, ref string, alt []string, qual float64, filter []string, info map[string]string, order []string) *vcf.Variant {
	return &vcf.Variant{Chrom: chrom, Pos: pos, ID: id, Ref: ref, Alt: alt, Qual: qual, Filter: filter, Info: info, InfoOrder: order}
}

func TestSetIDMacroExpansion(t *testing.T) {
	v := mkAnnVariant("chr1", 1000, "rsX", "AC", []string{"A", "ACT"}, 37, []string{"PASS"},
		map[string]string{"AF": "0.5", "FLG": ""}, []string{"AF", "FLG"})
	cases := []struct {
		format string
		want   string
	}{
		{`%CHROM`, "chr1"},
		{`%POS`, "1000"},
		{`%POS0`, "999"},
		{`%REF`, "AC"},
		{`%ALT`, "A,ACT"},
		{`%FIRST_ALT`, "A"},
		{`%QUAL`, "37"},
		{`%FILTER`, "PASS"},
		{`%ID`, "rsX"},
		{`%END`, "1001"},  // pos + len(REF) - 1
		{`%END0`, "1000"}, // END - 1
		{`%INFO/AF`, "0.5"},
		{`%AF`, "0.5"},     // bare tag -> INFO
		{`%INFO/FLG`, "1"}, // flag present -> 1
		{`%INFO/MISSING`, "."},
		{`%CHROM\_%POS`, "chr1_1000"},
		{`%CHROM:%POS:%REF:%FIRST_ALT`, "chr1:1000:AC:A"},
		{`literal-%CHROM-text`, "literal-chr1-text"},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.format, func(t *testing.T) {
			prog, err := ParseSetID(tc.format)
			if err != nil {
				t.Fatalf("ParseSetID(%q): %v", tc.format, err)
			}
			got, ok := prog.expand(v)
			if !ok {
				t.Fatalf("ParseSetID(%q): empty expansion", tc.format)
			}
			if got != tc.want {
				t.Fatalf("ParseSetID(%q): got %q want %q", tc.format, got, tc.want)
			}
		})
	}
}

func TestSetIDOnlyIfEmpty(t *testing.T) {
	prog, err := ParseSetID(`+%CHROM`)
	if err != nil {
		t.Fatalf("ParseSetID: %v", err)
	}
	if !prog.onlyIfEmpty {
		t.Fatalf("expected onlyIfEmpty for +-prefixed format")
	}
	prog2, _ := ParseSetID(`%CHROM`)
	if prog2.onlyIfEmpty {
		t.Fatalf("did not expect onlyIfEmpty for plain format")
	}
}

func TestSetIDTypeMacro(t *testing.T) {
	cases := []struct {
		ref  string
		alt  []string
		want string
	}{
		{"A", []string{"G"}, "SNP"},
		{"AT", []string{"A"}, "INDEL"},
		{"A", []string{"AT"}, "INDEL"},
		{"AT", []string{"GC"}, "MNP"},
		{"A", []string{"<DEL>"}, "OTHER"},
		{"A", []string{"*"}, "OVERLAP"},
		{"C", []string{"T", "TT"}, "SNP,OTHER"},
		{"G", []string{"."}, "REF"},
		{"A", []string{"a"}, "SNP"}, // case-sensitive single-base compare
		{"A", []string{"A"}, "REF"},
	}
	prog, _ := ParseSetID(`%TYPE`)
	for _, tc := range cases {
		tc := tc
		t.Run(tc.ref+">"+joinAlt(tc.alt), func(t *testing.T) {
			v := mkAnnVariant("1", 1, ".", tc.ref, tc.alt, -1, nil, nil, nil)
			got, _ := prog.expand(v)
			if got != tc.want {
				t.Fatalf("TYPE(%s>%v): got %q want %q", tc.ref, tc.alt, got, tc.want)
			}
		})
	}
}

func joinAlt(a []string) string {
	out := ""
	for i, s := range a {
		if i > 0 {
			out += ","
		}
		out += s
	}
	return out
}

func TestPairLogicTruthTable(t *testing.T) {
	// Each row: query REF/ALT vs annotation REF/ALT at the same site, plus the
	// expected match result per mode. Verified against upstream behaviour.
	type row struct {
		qRef, qAlt   string
		aRef, aAlt   string
		some, exact  bool
		all, snps    bool
		indels, both bool
	}
	rows := []row{
		// identical SNP: matches everywhere (shared allele).
		{"A", "G", "A", "G", true, true, true, true, true, true},
		// different SNP alleles, shared REF: only "all"/"any" and (no shared
		// allele) so SNP<->SNP needs SNPS score (snps mode is SNP_REF, not
		// SNP<->SNP) -> matches only under all/both.
		{"C", "A", "C", "T", false, false, true, false, false, true},
		// query SNP vs annotation REF-only: snps mode pairs (SNP_REF).
		{"A", "G", "A", ".", false, false, true, true, false, true},
		// query indel vs annotation REF-only: indels mode pairs.
		{"A", "AT", "A", ".", false, false, true, false, true, true},
		// identical indel: shared allele -> matches everywhere.
		{"G", "GA", "G", "GA", true, true, true, true, true, true},
	}
	modes := []struct {
		name  string
		logic PairLogic
		pick  func(r row) bool
	}{
		{"some", PairSome, func(r row) bool { return r.some }},
		{"exact", PairExact, func(r row) bool { return r.exact }},
		{"all", PairAll, func(r row) bool { return r.all }},
		{"snps", PairSNPs, func(r row) bool { return r.snps }},
		{"indels", PairIndels, func(r row) bool { return r.indels }},
		{"both", PairBoth, func(r row) bool { return r.both }},
	}
	for _, r := range rows {
		r := r
		q := mkAnnVariant("1", 100, ".", r.qRef, splitAltField(r.qAlt), -1, nil, nil, nil)
		a := mkAnnVariant("1", 100, ".", r.aRef, splitAltField(r.aAlt), -1, nil, nil, nil)
		for _, m := range modes {
			m := m
			want := m.pick(r)
			got := pairLogicMatches(q, a, m.logic)
			if got != want {
				t.Errorf("pair %s: %s>%s vs %s>%s got %v want %v",
					m.name, r.qRef, r.qAlt, r.aRef, r.aAlt, got, want)
			}
		}
	}
}

// splitAltField turns "." into a no-ALT record and otherwise splits on commas.
func splitAltField(s string) []string {
	if s == "." || s == "" {
		return nil
	}
	var out []string
	cur := ""
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			out = append(out, cur)
			cur = ""
			continue
		}
		cur += string(s[i])
	}
	return append(out, cur)
}

func TestParseMinOverlap(t *testing.T) {
	cases := []struct {
		in       string
		ann, vcf float64
		wantErr  bool
	}{
		{"", 0, 0, false},
		{"0.5", 0.5, 0, false},
		{"0.3:0.7", 0.3, 0.7, false},
		{":1", 0, 1, false},
		{"1:0", 1, 0, false},
		{"2", 0, 0, true},
		{"abc", 0, 0, true},
		{"0.5:9", 0, 0, true},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.in, func(t *testing.T) {
			mo, err := ParseMinOverlap(tc.in)
			if tc.wantErr {
				if err == nil {
					t.Fatalf("ParseMinOverlap(%q): expected error", tc.in)
				}
				return
			}
			if err != nil {
				t.Fatalf("ParseMinOverlap(%q): %v", tc.in, err)
			}
			if mo.Ann != tc.ann || mo.Vcf != tc.vcf {
				t.Fatalf("ParseMinOverlap(%q): got {%v %v} want {%v %v}", tc.in, mo.Ann, mo.Vcf, tc.ann, tc.vcf)
			}
		})
	}
}

func TestParseMergeLogic(t *testing.T) {
	m, err := ParseMergeLogicSpec("GENE:unique,SCORE:sum,DP:append-missing")
	if err != nil {
		t.Fatalf("ParseMergeLogicSpec: %v", err)
	}
	if m["GENE"] != MergeUnique || m["SCORE"] != MergeSum || m["DP"] != MergeAppendMissing {
		t.Fatalf("unexpected merge map: %+v", m)
	}
	if _, err := ParseMergeLogicSpec("BAD:nonsense"); err == nil {
		t.Fatalf("expected error for unknown logic")
	}
}

func TestParsePairLogic(t *testing.T) {
	cases := map[string]PairLogic{
		"some": PairSome, "exact": PairExact, "none": PairExact,
		"all": PairAll, "any": PairAll, "snps": PairSNPs,
		"indels": PairIndels, "both": PairBoth, "id": PairID,
	}
	for in, want := range cases {
		got, err := ParsePairLogic(in)
		if err != nil || got != want {
			t.Fatalf("ParsePairLogic(%q): got %v err %v want %v", in, got, err, want)
		}
	}
	if _, err := ParsePairLogic("bogus"); err == nil {
		t.Fatalf("expected error for bogus pair-logic")
	}
}
