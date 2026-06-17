package main

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// writeBGZF writes text as a bgzipped file under dir and returns its path.
func writeBGZF(t *testing.T, dir, name, text string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	defer f.Close()
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(text)); err != nil {
		t.Fatalf("write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	return path
}

func writePlain(t *testing.T, dir, name, text string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(text), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	return path
}

const sampleVCF = "##fileformat=VCFv4.2\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\t.\tA\tT\t.\t.\t.\n" +
	"chr1\t150\t.\tC\tG\t.\t.\t.\n" +
	"chr1\t200\t.\tG\tA\t.\t.\t.\n" +
	"chr2\t300\t.\tA\tG\t.\t.\t.\n"

// indexSample builds a sample bgzipped VCF plus .tbi and returns the data path.
func indexSample(t *testing.T, dir string) string {
	t.Helper()
	gz := writeBGZF(t, dir, "in.vcf.gz", sampleVCF)
	cfg, err := tabix.PresetConfig(tabix.PresetVCF)
	if err != nil {
		t.Fatalf("preset: %v", err)
	}
	idx, err := tabix.Build(gz, cfg)
	if err != nil {
		t.Fatalf("build: %v", err)
	}
	if err := idx.WriteFile(gz + ".tbi"); err != nil {
		t.Fatalf("write tbi: %v", err)
	}
	return gz
}

func TestRunTargetsStrictFilter(t *testing.T) {
	dir := t.TempDir()
	gz := indexSample(t, dir)
	targets := writePlain(t, dir, "t.txt", "chr1\t150\t200\n")

	var out, errb bytes.Buffer
	code := run([]string{"-p", "vcf", "-T", targets, gz}, nil, &out, &errb)
	if code != 0 {
		t.Fatalf("run exit=%d stderr=%s", code, errb.String())
	}
	want := "chr1\t150\t.\tC\tG\t.\t.\t.\nchr1\t200\t.\tG\tA\t.\t.\t.\n"
	if out.String() != want {
		t.Errorf("targets output:\n got %q\nwant %q", out.String(), want)
	}
}

func TestRunTargetsWithRegion(t *testing.T) {
	dir := t.TempDir()
	gz := indexSample(t, dir)
	targets := writePlain(t, dir, "t.txt", "chr1\t140\t160\n")

	var out, errb bytes.Buffer
	code := run([]string{"-p", "vcf", "-T", targets, gz, "chr1:100-200"}, nil, &out, &errb)
	if code != 0 {
		t.Fatalf("run exit=%d stderr=%s", code, errb.String())
	}
	want := "chr1\t150\t.\tC\tG\t.\t.\t.\n"
	if out.String() != want {
		t.Errorf("targets+region output:\n got %q\nwant %q", out.String(), want)
	}
}

func TestRunReheader(t *testing.T) {
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "in.vcf.gz", sampleVCF)
	newHdr := writePlain(t, dir, "hdr.txt",
		"##fileformat=VCFv4.2\n##NEW=yes\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	var out, errb bytes.Buffer
	code := run([]string{"-p", "vcf", "-r", newHdr, gz}, nil, &out, &errb)
	if code != 0 {
		t.Fatalf("run exit=%d stderr=%s", code, errb.String())
	}
	// Decompress the emitted stream and check the header was replaced.
	br, err := bgzip.NewReader(bytes.NewReader(out.Bytes()))
	if err != nil {
		t.Fatalf("bgzf reader: %v", err)
	}
	var dec bytes.Buffer
	if _, err := dec.ReadFrom(br); err != nil {
		t.Fatalf("decode: %v", err)
	}
	body := "chr1\t100\t.\tA\tT\t.\t.\t.\nchr1\t150\t.\tC\tG\t.\t.\t.\n" +
		"chr1\t200\t.\tG\tA\t.\t.\t.\nchr2\t300\t.\tA\tG\t.\t.\t.\n"
	want := "##fileformat=VCFv4.2\n##NEW=yes\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + body
	if dec.String() != want {
		t.Errorf("reheader output:\n got %q\nwant %q", dec.String(), want)
	}
}

func TestRunReheaderMissingHeaderFile(t *testing.T) {
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "in.vcf.gz", sampleVCF)
	var out, errb bytes.Buffer
	code := run([]string{"-p", "vcf", "-r", filepath.Join(dir, "nope.txt"), gz}, nil, &out, &errb)
	if code == 0 {
		t.Fatalf("expected non-zero exit for missing header file")
	}
}

func TestRunMissingDataFileArg(t *testing.T) {
	var out, errb bytes.Buffer
	if code := run(nil, nil, &out, &errb); code == 0 {
		t.Fatalf("expected error when no data file is given")
	}
}

// TestSortRegionsRegidx pins the htslib regidx region-iteration order used by
// tabix -R: chromosomes in first-appearance order, and within a chromosome by
// start ascending then end DESCENDING (a stable sort). The end-descending
// tie-break — two regions sharing a start emit the longer one first — is the
// non-obvious part that makes -R byte-exact against the upstream binary.
func TestSortRegionsRegidx(t *testing.T) {
	in := []region{
		{"chr1", 100, 200},
		{"chr1", 100, 500}, // same start, larger end -> must come BEFORE the 100-200 region
		{"chr1", 50, 60},
		{"chr2", 10, 20},
		{"chr1", 100, 500}, // exact duplicate -> stable, keeps input order after the first
		{"chr2", 5, 9},
	}
	got := sortRegionsRegidx(in)
	want := []region{
		{"chr1", 50, 60},
		{"chr1", 100, 500},
		{"chr1", 100, 500},
		{"chr1", 100, 200},
		{"chr2", 5, 9}, // chr2 in first-appearance order (after chr1), start asc
		{"chr2", 10, 20},
	}
	if len(got) != len(want) {
		t.Fatalf("len = %d, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("region %d = %+v, want %+v", i, got[i], want[i])
		}
	}
}
