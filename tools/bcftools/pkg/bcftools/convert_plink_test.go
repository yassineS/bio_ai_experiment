package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// plinkTestVCF is a 3-sample, 3-variant fixture exercising hom-REF, het,
// hom-ALT, missing, and a non-biallelic site (the third record) so the
// skip path is covered.
const plinkTestVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##contig=<ID=X,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100	rs1	A	G	.	.	.	GT	0/0	0/1	1/1
X	200	.	C	T	.	.	.	GT	1/1	./.	0|1
1	300	rs3	A	G,C	.	.	.	GT	0/1	0/2	1/2
`

// writePlinkInput writes the fixture VCF to a temp file and returns its
// path.
func writePlinkInput(t *testing.T, content string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write input: %v", err)
	}
	return p
}

func readFileStr(t *testing.T, name string) string {
	t.Helper()
	b, err := os.ReadFile(name)
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return string(b)
}

// --- .map / .ped text exporter ---------------------------------------------

func TestVCFToPlink_PedMap(t *testing.T) {
	in := writePlinkInput(t, plinkTestVCF)
	prefix := filepath.Join(t.TempDir(), "out")
	var stderr bytes.Buffer
	n, err := VCFToPlink(in, PlinkConvertOptions{Prefix: prefix}, &stderr)
	if err != nil {
		t.Fatalf("VCFToPlink: %v", err)
	}
	// Two biallelic records kept (the multi-allelic third is skipped).
	if n != 2 {
		t.Fatalf("records written = %d, want 2", n)
	}

	// .map: CHROM SNP_ID 0 BP. The X contig maps to chromosome code 23,
	// and the missing-ID record uses CHROM:POS.
	wantMap := "1\trs1\t0\t100\n" +
		"23\tX:200\t0\t200\n"
	if got := readFileStr(t, prefix+".map"); got != wantMap {
		t.Fatalf(".map mismatch:\n got: %q\nwant: %q", got, wantMap)
	}

	// .ped: one line per sample, 6 mandatory cols then the two allele
	// letters per variant. S2's missing X genotype is "0 0".
	wantPed := "S1 S1 0 0 0 -9 A A T T\n" +
		"S2 S2 0 0 0 -9 A G 0 0\n" +
		"S3 S3 0 0 0 -9 G G C T\n"
	if got := readFileStr(t, prefix+".ped"); got != wantPed {
		t.Fatalf(".ped mismatch:\n got: %q\nwant: %q", got, wantPed)
	}

	if !strings.Contains(stderr.String(), "non-biallelic records are skipped") {
		t.Errorf("expected multi-allelic warning, got: %q", stderr.String())
	}
}

// roundTripPed parses a .ped/.map pair back into per-(sample,variant)
// allele pairs keyed by SNP id, so the text round-trip can be asserted
// against the source genotypes.
func roundTripPed(t *testing.T, pedPath, mapPath string) (snpIDs []string, byVarSample map[string][]string) {
	t.Helper()
	for _, line := range strings.Split(strings.TrimRight(readFileStr(t, mapPath), "\n"), "\n") {
		f := strings.Fields(line)
		if len(f) != 4 {
			t.Fatalf("bad .map line: %q", line)
		}
		snpIDs = append(snpIDs, f[1])
	}
	byVarSample = make(map[string][]string)
	for _, line := range strings.Split(strings.TrimRight(readFileStr(t, pedPath), "\n"), "\n") {
		f := strings.Fields(line)
		if len(f) < 6 {
			t.Fatalf("bad .ped line: %q", line)
		}
		alleles := f[6:]
		if len(alleles) != 2*len(snpIDs) {
			t.Fatalf(".ped allele count = %d, want %d", len(alleles), 2*len(snpIDs))
		}
		for i, id := range snpIDs {
			pair := alleles[2*i] + alleles[2*i+1]
			byVarSample[id] = append(byVarSample[id], pair)
		}
	}
	return snpIDs, byVarSample
}

func TestVCFToPlink_PedMapRoundTrip(t *testing.T) {
	in := writePlinkInput(t, plinkTestVCF)
	prefix := filepath.Join(t.TempDir(), "out")
	if _, err := VCFToPlink(in, PlinkConvertOptions{Prefix: prefix}, nil); err != nil {
		t.Fatalf("VCFToPlink: %v", err)
	}
	ids, byVar := roundTripPed(t, prefix+".ped", prefix+".map")

	if len(ids) != 2 || ids[0] != "rs1" || ids[1] != "X:200" {
		t.Fatalf("snp ids = %v, want [rs1 X:200]", ids)
	}
	// rs1: S1 hom-REF(AA), S2 het(AG), S3 hom-ALT(GG).
	if got := byVar["rs1"]; got[0] != "AA" || got[1] != "AG" || got[2] != "GG" {
		t.Errorf("rs1 alleles = %v, want [AA AG GG]", got)
	}
	// X:200: S1 hom-ALT(TT), S2 missing(00), S3 phased het(CT).
	if got := byVar["X:200"]; got[0] != "TT" || got[1] != "00" || got[2] != "CT" {
		t.Errorf("X:200 alleles = %v, want [TT 00 CT]", got)
	}
}

// --- .tped / .tfam transposed text exporter --------------------------------

func TestVCFToPlinkTransposed(t *testing.T) {
	in := writePlinkInput(t, plinkTestVCF)
	prefix := filepath.Join(t.TempDir(), "out")
	n, err := VCFToPlinkTransposed(in, PlinkConvertOptions{Prefix: prefix}, nil)
	if err != nil {
		t.Fatalf("VCFToPlinkTransposed: %v", err)
	}
	if n != 2 {
		t.Fatalf("records written = %d, want 2", n)
	}

	// .tped: CHROM SNP_ID 0 BP then the two alleles for each sample.
	wantTped := "1 rs1 0 100 A A A G G G\n" +
		"23 X:200 0 200 T T 0 0 C T\n"
	if got := readFileStr(t, prefix+".tped"); got != wantTped {
		t.Fatalf(".tped mismatch:\n got: %q\nwant: %q", got, wantTped)
	}

	// .tfam: same 6 columns as the .ped prefix, one line per sample.
	wantTfam := "S1 S1 0 0 0 -9\nS2 S2 0 0 0 -9\nS3 S3 0 0 0 -9\n"
	if got := readFileStr(t, prefix+".tfam"); got != wantTfam {
		t.Fatalf(".tfam mismatch:\n got: %q\nwant: %q", got, wantTfam)
	}
}

// --- .bed / .bim / .fam binary exporter ------------------------------------

func TestVCFToPlinkBinary_BedBytes(t *testing.T) {
	in := writePlinkInput(t, plinkTestVCF)
	prefix := filepath.Join(t.TempDir(), "out")
	n, err := VCFToPlinkBinary(in, PlinkConvertOptions{Prefix: prefix}, nil)
	if err != nil {
		t.Fatalf("VCFToPlinkBinary: %v", err)
	}
	if n != 2 {
		t.Fatalf("records written = %d, want 2", n)
	}

	// Hand-encode the expected .bed body. 3 samples => ceil(3/4)=1 byte
	// per variant. Codes (A1=ALT, A2=REF; on-disk = #ALT alleles):
	//   00 = hom-REF, 10 = het, 11 = hom-ALT, 01 = missing.
	// Bits are little-endian within the byte: sample0=bits0-1, etc.
	//
	// rs1: S1 0/0->00, S2 0/1->10, S3 1/1->11
	//   byte = (00) | (10<<2) | (11<<4) = 0b00111000 = 0x38
	// X:200: S1 1/1->11, S2 ./.->01, S3 0|1->10
	//   byte = (11) | (01<<2) | (10<<4) = 0b00100111 = 0x27
	wantBed := []byte{0x6c, 0x1b, 0x01, 0x38, 0x27}
	gotBed, err := os.ReadFile(prefix + ".bed")
	if err != nil {
		t.Fatalf("read .bed: %v", err)
	}
	if !bytes.Equal(gotBed, wantBed) {
		t.Fatalf(".bed mismatch:\n got: % x\nwant: % x", gotBed, wantBed)
	}

	// .bim: CHROM SNP_ID 0 BP A1 A2 with A1=ALT, A2=REF.
	wantBim := "1\trs1\t0\t100\tG\tA\n" +
		"23\tX:200\t0\t200\tT\tC\n"
	if got := readFileStr(t, prefix+".bim"); got != wantBim {
		t.Fatalf(".bim mismatch:\n got: %q\nwant: %q", got, wantBim)
	}

	// .fam: the 6 .ped columns, one line per sample.
	wantFam := "S1 S1 0 0 0 -9\nS2 S2 0 0 0 -9\nS3 S3 0 0 0 -9\n"
	if got := readFileStr(t, prefix+".fam"); got != wantFam {
		t.Fatalf(".fam mismatch:\n got: %q\nwant: %q", got, wantFam)
	}
}

// TestVCFToPlinkBinary_Padding exercises the per-variant padding for a
// sample count that is not a multiple of 4 (5 samples => 2 bytes/variant,
// the high 6 bits of the second byte are zero padding).
func TestVCFToPlinkBinary_Padding(t *testing.T) {
	const vcf5 = `##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C	D	E
1	10	v1	A	T	.	.	.	GT	0/0	0/1	1/1	./.	0/1
`
	in := writePlinkInput(t, vcf5)
	prefix := filepath.Join(t.TempDir(), "out")
	if _, err := VCFToPlinkBinary(in, PlinkConvertOptions{Prefix: prefix}, nil); err != nil {
		t.Fatalf("VCFToPlinkBinary: %v", err)
	}
	// 5 samples => 2 bytes per variant.
	// byte0: A 0/0->00, B 0/1->10, C 1/1->11, D ./.->01
	//   = (00) | (10<<2) | (11<<4) | (01<<6) = 0b01111000 = 0x78
	// byte1: E 0/1->10 in bits0-1, rest padding 0
	//   = 0b00000010 = 0x02
	wantBed := []byte{0x6c, 0x1b, 0x01, 0x78, 0x02}
	gotBed, err := os.ReadFile(prefix + ".bed")
	if err != nil {
		t.Fatalf("read .bed: %v", err)
	}
	if !bytes.Equal(gotBed, wantBed) {
		t.Fatalf(".bed padding mismatch:\n got: % x\nwant: % x", gotBed, wantBed)
	}
}

// --- explicit comma-separated filenames ------------------------------------

func TestVCFToPlink_ExplicitNames(t *testing.T) {
	in := writePlinkInput(t, plinkTestVCF)
	dir := t.TempDir()
	ped := filepath.Join(dir, "geno.ped")
	mp := filepath.Join(dir, "geno.map")
	if _, err := VCFToPlink(in, PlinkConvertOptions{Prefix: ped + "," + mp}, nil); err != nil {
		t.Fatalf("VCFToPlink explicit names: %v", err)
	}
	if _, err := os.Stat(ped); err != nil {
		t.Errorf("expected %s to exist: %v", ped, err)
	}
	if _, err := os.Stat(mp); err != nil {
		t.Errorf("expected %s to exist: %v", mp, err)
	}
}

func TestPlinkChromCodes(t *testing.T) {
	cases := []struct{ in, want string }{
		{"1", "1"},
		{"22", "22"},
		{"chr7", "7"},
		{"X", "23"},
		{"chrX", "23"},
		{"Y", "24"},
		{"XY", "25"},
		{"MT", "26"},
		{"M", "26"},
		{"chrM", "26"},
		{"GL000192.1", "GL000192.1"},
		{"HLA-A", "HLA-A"},
	}
	for _, c := range cases {
		if got := plinkChrom(c.in); got != c.want {
			t.Errorf("plinkChrom(%q) = %q, want %q", c.in, got, c.want)
		}
	}
}

func TestBedCode(t *testing.T) {
	cases := []struct {
		gt   string
		want byte
	}{
		{"0/0", 0b00},
		{"0|0", 0b00},
		{"0/1", 0b10},
		{"1/0", 0b10},
		{"1/1", 0b11},
		{"./.", 0b01},
		{".", 0b01},
		{"0", 0b00},
		{"1", 0b11},
		{"0/.", 0b01},
	}
	for _, c := range cases {
		if got := bedCode(c.gt); got != c.want {
			t.Errorf("bedCode(%q) = %02b, want %02b", c.gt, got, c.want)
		}
	}
}
