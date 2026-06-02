package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const isecHdr = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
`

func isecVCF(records ...string) string {
	return isecHdr + strings.Join(records, "\n") + "\n"
}

func writeIsecInputs(t *testing.T, contents []string) []string {
	t.Helper()
	dir := t.TempDir()
	out := make([]string, len(contents))
	for i, c := range contents {
		p := filepath.Join(dir, "in"+string(rune('0'+i))+".vcf")
		if err := os.WriteFile(p, []byte(c), 0644); err != nil {
			t.Fatal(err)
		}
		out[i] = p
	}
	return out
}

// Hand-computed: input A has rs1, rs2; input B has rs1, rs3.
// Intersection (-n =2) yields rs1 (chr1:100 A→T) only. Without -p/-w the
// default output is upstream's tuple "list of sites" form
// (CHROM\tPOS\tREF\tALT\tBITS).
func TestIsecIntersection(t *testing.T) {
	a := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1",
	)
	b := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t300\trs3\tA\tT\t.\tPASS\t.\tGT\t0/1",
	)
	paths := writeIsecInputs(t, []string{a, b})
	var stdout bytes.Buffer
	n, err := IsecFiles(paths, &stdout, IsecOptions{
		Nfiles: NfilesSpec{Mode: '=', N: 2},
	})
	if err != nil {
		t.Fatal(err)
	}
	if n != 1 {
		t.Errorf("intersection size = %d, want 1", n)
	}
	want := "chr1\t100\tA\tT\t11\n"
	if stdout.String() != want {
		t.Errorf("tuple output mismatch:\ngot:  %q\nwant: %q", stdout.String(), want)
	}
}

// Hand-computed: union with -n+1 (without -n the two-input default is
// upstream's OP_VENN error). Tuple output lists every site, marked with
// per-input membership bits.
func TestIsecUnion(t *testing.T) {
	a := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1",
	)
	b := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t300\trs3\tA\tT\t.\tPASS\t.\tGT\t0/1",
	)
	paths := writeIsecInputs(t, []string{a, b})
	var stdout bytes.Buffer
	n, err := IsecFiles(paths, &stdout, IsecOptions{
		Nfiles: NfilesSpec{Mode: '+', N: 1},
	})
	if err != nil {
		t.Fatal(err)
	}
	if n != 3 {
		t.Errorf("union size = %d, want 3", n)
	}
	for _, want := range []string{"chr1\t100\tA\tT\t11", "chr1\t200\tC\tG\t10", "chr1\t300\tA\tT\t01"} {
		if !strings.Contains(stdout.String(), want) {
			t.Errorf("missing tuple line %q in:\n%s", want, stdout.String())
		}
	}
}

// Hand-computed: complement (records in A only) ⇒ -n ~10 — rs2 only, with
// tuple bits = 10.
func TestIsecBitmaskAOnly(t *testing.T) {
	a := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1",
	)
	b := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t300\trs3\tA\tT\t.\tPASS\t.\tGT\t0/1",
	)
	paths := writeIsecInputs(t, []string{a, b})
	var stdout bytes.Buffer
	n, err := IsecFiles(paths, &stdout, IsecOptions{
		Nfiles: NfilesSpec{Mode: '~', Bits: []bool{true, false}},
	})
	if err != nil {
		t.Fatal(err)
	}
	if n != 1 {
		t.Errorf("A-only = %d, want 1 (rs2)", n)
	}
	want := "chr1\t200\tC\tG\t10\n"
	if stdout.String() != want {
		t.Errorf("tuple output mismatch:\ngot:  %q\nwant: %q", stdout.String(), want)
	}
}

// Tuple "list of sites" default: upstream emits an advisory stderr note
// and the same tuple body. We capture the stderr writer to confirm parity.
func TestIsecTupleDefaultStderrNote(t *testing.T) {
	a := isecVCF("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	b := isecVCF("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	paths := writeIsecInputs(t, []string{a, b})
	var stdout, stderr bytes.Buffer
	if _, err := IsecFiles(paths, &stdout, IsecOptions{
		Nfiles: NfilesSpec{Mode: '=', N: 2},
		Stderr: &stderr,
	}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(stderr.String(), "printing list of sites") {
		t.Errorf("missing upstream stderr note; got:\n%s", stderr.String())
	}
}

// Two inputs, no -n: upstream switches to OP_VENN and requires -p, erroring
// with the exact text we mirror.
func TestIsecVennRequiresPrefix(t *testing.T) {
	a := isecVCF("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	b := isecVCF("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	paths := writeIsecInputs(t, []string{a, b})
	_, err := IsecFiles(paths, &bytes.Buffer{}, IsecOptions{})
	if err == nil || !strings.Contains(err.Error(), "Expected the -p option") {
		t.Errorf("want OP_VENN error 'Expected the -p option', got %v", err)
	}
}

// CollapseSome: REF must match AND at least one ALT must intersect. With
// fixture A=chr1:100 A→T, chr1:200 C→G,T; B=chr1:100 A→C, chr1:200 C→G:
// only chr1:200 collapses to a shared row (G in common). chr1:100 has
// disjoint ALTs, so it splits into two single-reader rows.
func TestIsecCollapseSome(t *testing.T) {
	a := isecVCF(
		"chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\t.\tC\tG,T\t.\tPASS\t.\tGT\t0/1",
	)
	b := isecVCF(
		"chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\t.\tC\tG\t.\tPASS\t.\tGT\t0/1",
	)
	paths := writeIsecInputs(t, []string{a, b})
	var stdout bytes.Buffer
	n, err := IsecFiles(paths, &stdout, IsecOptions{
		Collapse: CollapseSome,
		Nfiles:   NfilesSpec{Mode: '=', N: 2},
	})
	if err != nil {
		t.Fatal(err)
	}
	// =2 keeps only the chr1:200 cluster; chr1:100 splits into two
	// non-intersecting singletons whose bits are 10 and 01.
	if n != 1 {
		t.Errorf("CollapseSome =2 size = %d, want 1", n)
	}
	want := "chr1\t200\tC\tG,T\t11\n"
	if stdout.String() != want {
		t.Errorf("CollapseSome tuple output mismatch:\ngot:  %q\nwant: %q", stdout.String(), want)
	}
}

// Same fixture under -n+1: confirms the non-intersecting chr1:100 rows are
// split rather than merged.
func TestIsecCollapseSomeSplitDisjointAlts(t *testing.T) {
	a := isecVCF(
		"chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
	)
	b := isecVCF(
		"chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
	)
	paths := writeIsecInputs(t, []string{a, b})
	var stdout bytes.Buffer
	n, err := IsecFiles(paths, &stdout, IsecOptions{
		Collapse: CollapseSome,
		Nfiles:   NfilesSpec{Mode: '+', N: 1},
	})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("CollapseSome split count = %d, want 2", n)
	}
	out := stdout.String()
	if !strings.Contains(out, "chr1\t100\tA\tT\t10") || !strings.Contains(out, "chr1\t100\tA\tC\t01") {
		t.Errorf("CollapseSome split output missing one of the rows:\n%s", out)
	}
}

func TestIsecPrefixMode(t *testing.T) {
	a := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1",
	)
	b := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",
	)
	paths := writeIsecInputs(t, []string{a, b})
	prefixDir := t.TempDir()
	var stdout bytes.Buffer
	if _, err := IsecFiles(paths, &stdout, IsecOptions{
		Nfiles: NfilesSpec{Mode: '+', N: 1},
		Prefix: prefixDir,
	}); err != nil {
		t.Fatal(err)
	}
	for _, name := range []string{"0000.vcf", "0001.vcf", "sites.txt", "README.txt"} {
		_, err := os.Stat(filepath.Join(prefixDir, name))
		if err != nil {
			t.Errorf("missing prefix file %s: %v", name, err)
		}
	}
	// 0000.vcf should hold all of A's records (rs1, rs2); 0001.vcf only rs1.
	a0, _ := os.ReadFile(filepath.Join(prefixDir, "0000.vcf"))
	a1, _ := os.ReadFile(filepath.Join(prefixDir, "0001.vcf"))
	if !bytes.Contains(a0, []byte("rs1")) || !bytes.Contains(a0, []byte("rs2")) {
		t.Errorf("0000.vcf missing records:\n%s", a0)
	}
	if !bytes.Contains(a1, []byte("rs1")) || bytes.Contains(a1, []byte("rs2")) {
		t.Errorf("0001.vcf shape wrong:\n%s", a1)
	}
}

func TestParseNfilesSpec(t *testing.T) {
	type want struct {
		mode byte
		n    int
		bits []bool
		ok   bool
	}
	cases := []struct {
		in   string
		want want
	}{
		{"", want{ok: true}},
		{"+2", want{mode: '+', n: 2, ok: true}},
		{"=3", want{mode: '=', n: 3, ok: true}},
		{"~10", want{mode: '~', bits: []bool{true, false}, ok: true}},
		{"3", want{mode: '+', n: 3, ok: true}},
		{"foo", want{ok: false}},
		{"~1x", want{ok: false}},
	}
	for _, c := range cases {
		got, err := ParseNfilesSpec(c.in)
		if !c.want.ok {
			if err == nil {
				t.Errorf("%q should error", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("%q err: %v", c.in, err)
			continue
		}
		if got.Mode != c.want.mode {
			t.Errorf("%q mode = %q want %q", c.in, got.Mode, c.want.mode)
		}
		if got.N != c.want.n {
			t.Errorf("%q n = %d want %d", c.in, got.N, c.want.n)
		}
		if c.want.bits != nil {
			if len(got.Bits) != len(c.want.bits) {
				t.Errorf("%q bits len mismatch", c.in)
				continue
			}
			for i := range got.Bits {
				if got.Bits[i] != c.want.bits[i] {
					t.Errorf("%q bits[%d] = %v want %v", c.in, i, got.Bits[i], c.want.bits[i])
				}
			}
		}
	}
}

func TestParseCollapseMode(t *testing.T) {
	for _, in := range []string{"none", "snps", "indels", "both", "all", "some", "id", ""} {
		if _, err := ParseCollapseMode(in); err != nil {
			t.Errorf("ParseCollapseMode(%q): %v", in, err)
		}
	}
	if _, err := ParseCollapseMode("xx"); err == nil {
		t.Error("expected error for unknown mode")
	}
}

func TestIsecNeedsTwoInputs(t *testing.T) {
	if _, err := IsecFiles([]string{"only-one.vcf"}, &bytes.Buffer{}, IsecOptions{}); err == nil {
		t.Error("expected error for single input")
	}
}

func TestIsecBitmaskWrongLen(t *testing.T) {
	a := isecVCF("chr1\t1\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	b := isecVCF("chr1\t1\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	paths := writeIsecInputs(t, []string{a, b})
	if _, err := IsecFiles(paths, &bytes.Buffer{}, IsecOptions{
		Nfiles: NfilesSpec{Mode: '~', Bits: []bool{true}}, // wrong length
	}); err == nil {
		t.Error("expected error for bitmask length mismatch")
	}
}
