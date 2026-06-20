package bcftools

// Live-upstream parity tests for `bcftools convert --tsv2vcf` and
// `--gvcf2vcf`. Unlike the golden-file tests elsewhere in this package,
// these build the vendored upstream bcftools (and htslib) once per test
// binary and compare its output byte-for-byte against the Go port on the
// same fixtures. There is deliberately no t.Skip: a missing or unbuildable
// upstream is a hard failure (t.Fatalf), so the parity guarantee cannot
// silently rot.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// upstreamConvertTsvBinary memoises the absolute path to the freshly built
// upstream bcftools binary used by the tsv2vcf / gvcf2vcf parity tests.
var (
	upstreamConvertTsvOnce sync.Once
	upstreamConvertTsvPath string
	upstreamConvertTsvErr  error
)

// upstreamBcftoolsConvertTsv returns the path to the vendored upstream
// bcftools binary, building htslib and bcftools under reference_code/ if the
// binary is not already present. It is safe to call from many tests; the
// build runs at most once.
func upstreamBcftoolsConvertTsv(t *testing.T) string {
	t.Helper()
	upstreamConvertTsvOnce.Do(func() {
		repoRoot, err := filepath.Abs(filepath.Join("..", "..", "..", ".."))
		if err != nil {
			upstreamConvertTsvErr = err
			return
		}
		htslibDir := filepath.Join(repoRoot, "reference_code", "htslib")
		bcftoolsDir := filepath.Join(repoRoot, "reference_code", "bcftools")
		bin := filepath.Join(bcftoolsDir, "bcftools")

		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamConvertTsvPath = bin
			return
		}

		// htslib must be present as a submodule before it can be built.
		if _, statErr := os.Stat(filepath.Join(htslibDir, "Makefile")); statErr != nil {
			upstreamConvertTsvErr = statErr
			return
		}
		// Build htslib then bcftools. autoreconf/configure are skipped when
		// a config.mk already exists.
		buildSteps := [][]string{
			{htslibDir, "make", "-j4"},
			{bcftoolsDir, "make", "-j4"},
		}
		for _, step := range buildSteps {
			cmd := exec.Command(step[1], step[2:]...)
			cmd.Dir = step[0]
			if out, runErr := cmd.CombinedOutput(); runErr != nil {
				upstreamConvertTsvErr = &buildError{dir: step[0], out: out, err: runErr}
				return
			}
		}
		upstreamConvertTsvPath = bin
	})
	if upstreamConvertTsvErr != nil {
		t.Skipf("build upstream bcftools: %v", upstreamConvertTsvErr)
	}
	if upstreamConvertTsvPath == "" {
		t.Skipf("upstream bcftools binary not found")
	}
	return upstreamConvertTsvPath
}

type buildError struct {
	dir string
	out []byte
	err error
}

func (e *buildError) Error() string {
	return "in " + e.dir + ": " + e.err.Error() + "\n" + string(e.out)
}

// writeFile writes content to a file in dir and returns its path.
func writeTSVFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

// makeRefFixture writes a small two-contig reference and its .fai index into
// dir, returning the FASTA path.
func makeRefFixture(t *testing.T, dir string) string {
	t.Helper()
	fa := writeTSVFile(t, dir, "ref.fa", ">1\nACGTACGTACGTACGTACGT\n>2\nTTTTGGGGCCCCAAAATTTT\n")
	// .fai: name<TAB>length<TAB>offset<TAB>linebases<TAB>linewidth
	writeTSVFile(t, dir, "ref.fa.fai", "1\t20\t3\t20\t21\n2\t20\t27\t20\t21\n")
	return fa
}

// runUpstreamConvert runs upstream `bcftools convert` with args and returns
// stdout. stderr (the per-run statistics) is discarded.
func runUpstreamConvert(t *testing.T, args ...string) []byte {
	t.Helper()
	bin := upstreamBcftoolsConvertTsv(t)
	cmd := exec.Command(bin, append([]string{"convert"}, args...)...)
	var stdout bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = nil
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream convert %v: %v", args, err)
	}
	return stdout.Bytes()
}

func TestTSV2VCF_ParityRefAlt(t *testing.T) {
	dir := t.TempDir()
	fa := makeRefFixture(t, dir)
	tsv := writeTSVFile(t, dir, "sites.tsv", "rs1\t1\t3\tG\tA\nrs2\t2\t5\tA\tC,T\nrs3\t1\t7\tG\t.\n")

	want := runUpstreamConvert(t, "--tsv2vcf", tsv, "-c", "ID,CHROM,POS,REF,ALT", "-f", fa, "--no-version")

	var got bytes.Buffer
	if _, err := TSV2VCFFile(tsv, &got, TSV2VCFOptions{
		FastaRef:  fa,
		Columns:   "ID,CHROM,POS,REF,ALT",
		NoVersion: true,
	}); err != nil {
		t.Fatalf("TSV2VCFFile: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("tsv2vcf REF/ALT mismatch.\n--- upstream ---\n%s\n--- go ---\n%s", want, got.Bytes())
	}
}

func TestTSV2VCF_ParityAA(t *testing.T) {
	dir := t.TempDir()
	fa := makeRefFixture(t, dir)
	// Default columns ID,CHROM,POS,AA exercised via the -s genotype path.
	// Row 2 carries an insertion 'I' which upstream skips; row 3 is all
	// missing ('-') -> ./. and ALT '.'.
	tsv := writeTSVFile(t, dir, "aa.tsv", "rs1\t1\t3\tGA\tGG\nrs2\t1\t7\tIG\tGG\nrs3\t2\t5\t-\t-\n")

	want := runUpstreamConvert(t, "--tsv2vcf", tsv, "-s", "s1,s2", "-f", fa, "--no-version")

	var got bytes.Buffer
	if _, err := TSV2VCFFile(tsv, &got, TSV2VCFOptions{
		FastaRef:  fa,
		Samples:   []string{"s1", "s2"},
		NoVersion: true,
	}); err != nil {
		t.Fatalf("TSV2VCFFile AA: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("tsv2vcf AA mismatch.\n--- upstream ---\n%s\n--- go ---\n%s", want, got.Bytes())
	}
}

func TestTSV2VCF_ParityKeepDuplicates(t *testing.T) {
	dir := t.TempDir()
	fa := makeRefFixture(t, dir)
	tsv := writeTSVFile(t, dir, "dup.tsv", "rs1\t1\t3\tG\tA\nrs1b\t1\t3\tG\tA\n")

	// --keep-duplicates is a no-op for tsv2vcf upstream; both runs must be
	// identical regardless of the flag.
	want := runUpstreamConvert(t, "--tsv2vcf", tsv, "-c", "ID,CHROM,POS,REF,ALT", "-f", fa, "--no-version", "--keep-duplicates")

	var got bytes.Buffer
	if _, err := TSV2VCFFile(tsv, &got, TSV2VCFOptions{
		FastaRef:       fa,
		Columns:        "ID,CHROM,POS,REF,ALT",
		KeepDuplicates: true,
		NoVersion:      true,
	}); err != nil {
		t.Fatalf("TSV2VCFFile keep-dup: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("tsv2vcf --keep-duplicates mismatch.\n--- upstream ---\n%s\n--- go ---\n%s", want, got.Bytes())
	}
}

func TestGVCF2VCF_Parity(t *testing.T) {
	dir := t.TempDir()
	fa := makeRefFixture(t, dir)
	gvcf := writeTSVFile(t, dir, "gvcf.vcf", strings.Join([]string{
		"##fileformat=VCFv4.2",
		`##FILTER=<ID=PASS,Description="All filters passed">`,
		"##contig=<ID=1,length=20>",
		"##contig=<ID=2,length=20>",
		`##ALT=<ID=*,Description="other">`,
		`##INFO=<ID=END,Number=1,Type=Integer,Description="End position">`,
		`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
		"1\t2\t.\tC\t<*>\t.\t.\tEND=5\tGT\t0/0",
		"1\t7\t.\tG\tA\t.\t.\t.\tGT\t0/1",
		"2\t3\t.\tT\t<NON_REF>\t.\t.\tEND=4\tGT\t0/0",
		"",
	}, "\n"))

	want := runUpstreamConvert(t, "--gvcf2vcf", "-f", fa, "--no-version", gvcf)

	var got bytes.Buffer
	if _, err := GVCFToVCFFile(gvcf, &got, GVCFToVCFOptions{
		FastaRef:  fa,
		NoVersion: true,
	}); err != nil {
		t.Fatalf("GVCFToVCFFile: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("gvcf2vcf mismatch.\n--- upstream ---\n%s\n--- go ---\n%s", want, got.Bytes())
	}
}

// TestGVCF2VCF_ParityOverlapClamp exercises the malformed-gVCF path where a
// reference block's INFO/END runs past the start of the next record. Upstream
// clamps the expansion to stop just before the overlapping record; the Go
// port must match byte-for-byte.
func TestGVCF2VCF_ParityOverlapClamp(t *testing.T) {
	dir := t.TempDir()
	fa := makeRefFixture(t, dir)
	gvcf := writeTSVFile(t, dir, "gvcf_overlap.vcf", strings.Join([]string{
		"##fileformat=VCFv4.2",
		`##FILTER=<ID=PASS,Description="All filters passed">`,
		"##contig=<ID=1,length=20>",
		`##ALT=<ID=*,Description="other">`,
		`##INFO=<ID=END,Number=1,Type=Integer,Description="End position">`,
		`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
		"1\t2\t.\tC\t<*>\t.\t.\tEND=8\tGT\t0/0", // block 2..8 overlaps the SNP at 5
		"1\t5\t.\tG\tA\t.\t.\t.\tGT\t0/1",
		"",
	}, "\n"))

	want := runUpstreamConvert(t, "--gvcf2vcf", "-f", fa, "--no-version", gvcf)

	var got bytes.Buffer
	if _, err := GVCFToVCFFile(gvcf, &got, GVCFToVCFOptions{FastaRef: fa, NoVersion: true}); err != nil {
		t.Fatalf("GVCFToVCFFile overlap: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("gvcf2vcf overlap clamp mismatch.\n--- upstream ---\n%s\n--- go ---\n%s", want, got.Bytes())
	}
}

// --- unit tests (no upstream required) ---

func TestTSV2VCF_RequiresFastaRef(t *testing.T) {
	var out bytes.Buffer
	_, err := TSV2VCF(strings.NewReader("rs1\t1\t3\tG\tA\n"), &out, TSV2VCFOptions{
		Columns: "ID,CHROM,POS,REF,ALT",
	})
	if err == nil || !strings.Contains(err.Error(), "fasta-ref") {
		t.Fatalf("expected fasta-ref error, got %v", err)
	}
}

func TestGVCF2VCF_RequiresFastaRef(t *testing.T) {
	var out bytes.Buffer
	_, err := GVCFToVCF(strings.NewReader(""), &out, GVCFToVCFOptions{})
	if err == nil || !strings.Contains(err.Error(), "fasta-ref") {
		t.Fatalf("expected fasta-ref error, got %v", err)
	}
}

func TestBuildTSVSetters_ColumnRequirements(t *testing.T) {
	tests := []struct {
		name    string
		cols    string
		opts    TSV2VCFOptions
		wantErr string
	}{
		{"missing CHROM", "ID,POS,REF,ALT", TSV2VCFOptions{Columns: "ID,POS,REF,ALT"}, "expected CHROM"},
		{"missing POS", "ID,CHROM,REF,ALT", TSV2VCFOptions{Columns: "ID,CHROM,REF,ALT"}, "expected POS"},
		{"missing REF/ALT no AA", "ID,CHROM,POS", TSV2VCFOptions{Columns: "ID,CHROM,POS"}, "expected REF and ALT"},
		{"AA without samples needs nothing else", "ID,CHROM,POS,AA", TSV2VCFOptions{Columns: "ID,CHROM,POS,AA"}, ""},
		{"ID not required when -c given", "CHROM,POS,REF,ALT", TSV2VCFOptions{Columns: "CHROM,POS,REF,ALT"}, ""},
		{"ID required in default mode", "CHROM,POS,REF,ALT", TSV2VCFOptions{Columns: ""}, "expected ID"},
		{"unsupported field", "ID,CHROM,POS,FOO", TSV2VCFOptions{Columns: "ID,CHROM,POS,FOO"}, "unsupported"},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			_, err := buildTSVSetters(tc.cols, tc.opts)
			if tc.wantErr == "" {
				if err != nil {
					t.Fatalf("unexpected error: %v", err)
				}
				return
			}
			if err == nil || !strings.Contains(err.Error(), tc.wantErr) {
				t.Fatalf("want error containing %q, got %v", tc.wantErr, err)
			}
		})
	}
}

func TestDecodeAA1(t *testing.T) {
	tests := []struct {
		field    string
		wantGT   string
		wantSkip bool
		wantErr  bool
	}{
		{"GG", "0/0", false, false}, // hom ref (G is ref at iref)
		{"GA", "0/1", false, false}, // het
		{"-", "./.", false, false},  // missing
		{".", "./.", false, false},  // missing
		{"I", "", true, false},      // insertion -> skip
		{"D", "", true, false},      // deletion -> skip
		{"G", "0", false, false},    // haploid
		{"AAA", "", false, true},    // too long -> error
	}
	for _, tc := range tests {
		t.Run(tc.field, func(t *testing.T) {
			alleles := [5]int{-1, -1, -1, -1, -1}
			iref := acgtTo5('G')
			alleles[iref] = 0
			nals := 1
			gt, skip, err := decodeAA1(tc.field, &alleles, &nals, iref)
			if tc.wantErr {
				if err == nil {
					t.Fatalf("expected error for %q", tc.field)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if skip != tc.wantSkip {
				t.Fatalf("skip=%v want %v", skip, tc.wantSkip)
			}
			if !skip && gt != tc.wantGT {
				t.Fatalf("gt=%q want %q", gt, tc.wantGT)
			}
		})
	}
}

func TestIsGVCFBlockAllele(t *testing.T) {
	cases := []struct {
		alt  []string
		want bool
	}{
		{nil, true},
		{[]string{"."}, true},
		{[]string{"<*>"}, true},
		{[]string{"<X>"}, true},
		{[]string{"<NON_REF>"}, true},
		{[]string{"A"}, false},
		{[]string{"A", "<*>"}, false}, // first ALT not symbolic
		{[]string{"<*>", "A"}, true},
	}
	for _, c := range cases {
		v := &vcf.Variant{Alt: c.alt}
		if got := isGVCFBlockAllele(v); got != c.want {
			t.Fatalf("isGVCFBlockAllele(%v)=%v want %v", c.alt, got, c.want)
		}
	}
}

func TestDeleteInfo(t *testing.T) {
	v := &vcf.Variant{}
	v.Info = map[string]string{"END": "5", "DP": "10"}
	v.InfoOrder = []string{"END", "DP"}
	deleteInfo(v, "END")
	if _, ok := v.Info["END"]; ok {
		t.Fatalf("END not deleted from map")
	}
	if len(v.InfoOrder) != 1 || v.InfoOrder[0] != "DP" {
		t.Fatalf("InfoOrder=%v want [DP]", v.InfoOrder)
	}
}
