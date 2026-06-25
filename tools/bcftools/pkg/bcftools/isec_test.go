package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TestIsecStreamingMatchesInMemory pins the streaming k-way merge byte-identical
// to the whole-corpus in-memory path. It runs the same streamable inputs through
// IsecFiles (which takes isecStreaming) and through Isec on pre-read groups,
// across several forced flush thresholds — including 1 (flush every position) —
// so a batch boundary lands between, and inside, runs of intra-position
// duplicates that must stay paired across inputs. The inputs span two contigs,
// same-position different-allele keys, private records on both sides, and an
// unpaired duplicate occurrence.
func TestIsecStreamingMatchesInMemory(t *testing.T) {
	hdr := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
`
	mk := func(recs ...string) string { return hdr + strings.Join(recs, "\n") + "\n" }
	a := mk(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1",  // shared
		"chr1\t100\trs1b\tA\tC\t.\tPASS\t.\tGT\t0/1", // same pos, different ALT key
		"chr1\t200\trsA\tG\tC\t.\tPASS\t.\tGT\t0/1",  // private to A
		"chr1\t300\trs2\tA\tT\t.\tPASS\t.\tGT\t0/1",  // dup-key occ 1, shared
		"chr1\t300\trs2\tA\tT\t.\tPASS\t.\tGT\t0/1",  // dup-key occ 2, no B partner
		"chr2\t50\trsC\tA\tT\t.\tPASS\t.\tGT\t0/1",   // shared on chr2
		"chr2\t60\trsD\tA\tT\t.\tPASS\t.\tGT\t0/1",   // private to A
	)
	b := mk(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1", // shared
		"chr1\t180\trsB\tT\tA\t.\tPASS\t.\tGT\t0/1", // private to B
		"chr1\t300\trs2\tA\tT\t.\tPASS\t.\tGT\t0/1", // pairs with A's first 300
		"chr2\t50\trsC\tA\tT\t.\tPASS\t.\tGT\t0/1",  // shared on chr2
	)
	paths := writeIsecInputs(t, []string{a, b})

	// In-memory reference: read both files whole and run Isec once.
	readGroups := func() ([]*vcf.Header, [][]*vcf.Variant) {
		hs := make([]*vcf.Header, len(paths))
		gs := make([][]*vcf.Variant, len(paths))
		for i, p := range paths {
			in, err := iohelper.OpenReader(p)
			if err != nil {
				t.Fatalf("open %s: %v", p, err)
			}
			h, recs, err := readAllVariants(in)
			_ = in.Close()
			if err != nil {
				t.Fatalf("read %s: %v", p, err)
			}
			hs[i], gs[i] = h, recs
		}
		return hs, gs
	}

	optsCases := []struct {
		name   string
		mk     func(prefix string) IsecOptions
		prefix bool
	}{
		{"venn_prefix", func(p string) IsecOptions { return IsecOptions{Prefix: p} }, true},
		{"intersect_stdout", func(string) IsecOptions { return IsecOptions{Nfiles: NfilesSpec{Mode: '=', N: 2}} }, false},
		{"union_stdout", func(string) IsecOptions { return IsecOptions{} }, false},
		{"collapse_all_prefix", func(p string) IsecOptions { return IsecOptions{Prefix: p, Collapse: CollapseAll} }, true},
	}

	defaultFlush := isecStreamFlushRecords
	defer func() { isecStreamFlushRecords = defaultFlush }()

	for _, oc := range optsCases {
		for _, flush := range []int{1, 2, 3, defaultFlush} {
			t.Run(oc.name+"/flush="+strconv.Itoa(flush), func(t *testing.T) {
				isecStreamFlushRecords = flush

				// Streaming path via IsecFiles.
				streamDir := t.TempDir()
				var streamOut bytes.Buffer
				nStream, err := IsecFiles(paths, &streamOut, oc.mk(streamDir))
				if err != nil {
					t.Fatalf("streaming IsecFiles: %v", err)
				}

				// In-memory path via Isec on pre-read groups.
				memDir := t.TempDir()
				var memOut bytes.Buffer
				hs, gs := readGroups()
				nMem, err := Isec(hs, gs, &memOut, oc.mk(memDir))
				if err != nil {
					t.Fatalf("in-memory Isec: %v", err)
				}

				if nStream != nMem {
					t.Fatalf("kept count: streaming %d != in-memory %d", nStream, nMem)
				}
				if !bytes.Equal(streamOut.Bytes(), memOut.Bytes()) {
					t.Fatalf("stdout differs:\n--- streaming ---\n%s\n--- in-memory ---\n%s",
						streamOut.String(), memOut.String())
				}
				if oc.prefix {
					// Compare every emitted file except README.txt, which embeds
					// the (differing) prefix-dir paths.
					ents, err := os.ReadDir(streamDir)
					if err != nil {
						t.Fatal(err)
					}
					for _, e := range ents {
						if e.Name() == "README.txt" {
							continue
						}
						s, _ := os.ReadFile(filepath.Join(streamDir, e.Name()))
						m, err := os.ReadFile(filepath.Join(memDir, e.Name()))
						if err != nil {
							t.Fatalf("in-memory missing %s: %v", e.Name(), err)
						}
						if !bytes.Equal(s, m) {
							t.Fatalf("%s differs:\n--- streaming ---\n%s\n--- in-memory ---\n%s",
								e.Name(), s, m)
						}
					}
				}
			})
		}
	}
}

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
// Intersection (-n =2) yields rs1 only. Default writes from input A.
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
	if !strings.Contains(stdout.String(), "rs1") {
		t.Errorf("expected rs1 in stdout dump:\n%s", stdout.String())
	}
	if strings.Contains(stdout.String(), "rs2") || strings.Contains(stdout.String(), "rs3") {
		t.Errorf("non-intersection records leaked:\n%s", stdout.String())
	}
}

// Hand-computed: union (no -n constraint) = rs1, rs2, rs3.
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
	n, err := IsecFiles(paths, &stdout, IsecOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if n != 3 {
		t.Errorf("union size = %d, want 3", n)
	}
}

// Hand-computed: complement (records in A only) ⇒ -n ~10
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
	if !strings.Contains(stdout.String(), "rs2") {
		t.Errorf("rs2 missing from stdout:\n%s", stdout.String())
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

// TestIsecVennMode verifies the two-input, no-constraint Venn decomposition:
// records private to input 1 -> 0000, private to input 2 -> 0001, shared
// records -> 0002 (input-1 copy) and 0003 (input-2 copy), and one sites.txt
// line per matched occurrence (including intra-position duplicates paired
// across inputs).
func TestIsecVennMode(t *testing.T) {
	a := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1", // shared
		"chr1\t150\trsA\tG\tC\t.\tPASS\t.\tGT\t0/1", // private to A
		"chr1\t200\trs2\tA\tT\t.\tPASS\t.\tGT\t0/1", // duplicate key, shared
		"chr1\t200\trs3\tA\tT\t.\tPASS\t.\tGT\t0/1", // duplicate key, shared
	)
	b := isecVCF(
		"chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1", // shared
		"chr1\t180\trsB\tT\tA\t.\tPASS\t.\tGT\t0/1", // private to B
		"chr1\t200\trs2\tA\tT\t.\tPASS\t.\tGT\t0/1", // pairs with A's first 200
	)
	paths := writeIsecInputs(t, []string{a, b})
	prefixDir := t.TempDir()
	if _, err := IsecFiles(paths, &bytes.Buffer{}, IsecOptions{Prefix: prefixDir}); err != nil {
		t.Fatal(err)
	}
	read := func(name string) string {
		b, err := os.ReadFile(filepath.Join(prefixDir, name))
		if err != nil {
			t.Fatalf("read %s: %v", name, err)
		}
		return string(b)
	}
	for _, name := range []string{"0000.vcf", "0001.vcf", "0002.vcf", "0003.vcf", "sites.txt", "README.txt"} {
		if _, err := os.Stat(filepath.Join(prefixDir, name)); err != nil {
			t.Errorf("missing %s: %v", name, err)
		}
	}
	if f0 := read("0000.vcf"); !strings.Contains(f0, "rsA") || strings.Contains(f0, "rs1") {
		t.Errorf("0000.vcf (private to A) wrong:\n%s", f0)
	}
	if f1 := read("0001.vcf"); !strings.Contains(f1, "rsB") || strings.Contains(f1, "rs1") {
		t.Errorf("0001.vcf (private to B) wrong:\n%s", f1)
	}
	// rs1 + the two 200 records: A's first 200 pairs with B's 200 (shared);
	// A's second 200 has no B partner (private to A -> 0000).
	if f2 := read("0002.vcf"); !strings.Contains(f2, "rs1") || !strings.Contains(f2, "rs2") {
		t.Errorf("0002.vcf (shared, A copy) wrong:\n%s", f2)
	}
	if f3 := read("0003.vcf"); !strings.Contains(f3, "rs1") {
		t.Errorf("0003.vcf (shared, B copy) wrong:\n%s", f3)
	}
	// The unpaired second chr1:200 record from A is private to A.
	if f0 := read("0000.vcf"); !strings.Contains(f0, "rs3") {
		t.Errorf("0000.vcf should hold the unpaired duplicate rs3:\n%s", f0)
	}
	// README uses upstream's "for records private to" / "shared by both" legend.
	if rm := read("README.txt"); !strings.Contains(rm, "for records private to") ||
		!strings.Contains(rm, "shared by both") || !strings.Contains(rm, "produced by vcfisec") {
		t.Errorf("README.txt legend wrong:\n%s", rm)
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
