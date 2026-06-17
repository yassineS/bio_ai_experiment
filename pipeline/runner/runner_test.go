package runner

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
)

// fakeManifest returns a manifest with the new fixture kinds populated.
func fakeManifest() *fixtures.Manifest {
	return &fixtures.Manifest{
		Files: map[string]string{
			"fastq":           "/fix/reads.fastq",
			"fastq1":          "/fix/r1.fastq",
			"fastq2":          "/fix/r2.fastq",
			"gff":             "/fix/a.gff3",
			"vcf_multi_plain": "/fix/multi.vcf",
			"bed":             "/fix/i.bed",
		},
	}
}

// TestResolvePlaceholders_NewKinds covers the FASTQ/GFF/multi-VCF tokens and
// the per-invocation {out} prefix.
func TestResolvePlaceholders_NewKinds(t *testing.T) {
	m := fakeManifest()
	got, err := resolvePlaceholders([]string{"-i", "{fastq}", "--in2", "{fastq2}", "-g", "{gff}", "-o", "{out}.fastq"}, m, "/tmp/run/out")
	if err != nil {
		t.Fatalf("resolve: %v", err)
	}
	want := []string{"-i", "/fix/reads.fastq", "--in2", "/fix/r2.fastq", "-g", "/fix/a.gff3", "-o", "/tmp/run/out.fastq"}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("arg %d = %q, want %q", i, got[i], want[i])
		}
	}
}

// TestResolvePlaceholders_OutWithoutPrefix errors when {out} is used but no
// prefix was provided (an entry that forgot to declare OutputFiles).
func TestResolvePlaceholders_OutWithoutPrefix(t *testing.T) {
	if _, err := resolvePlaceholders([]string{"-o", "{out}.x"}, fakeManifest(), ""); err == nil {
		t.Fatal("expected error for {out} without a prefix")
	}
}

// TestResolvePlaceholders_MissingFixture errors clearly for an unknown kind.
func TestResolvePlaceholders_MissingFixture(t *testing.T) {
	if _, err := resolvePlaceholders([]string{"{cram}"}, fakeManifest(), ""); err == nil {
		t.Fatal("expected error for fixture absent from the manifest")
	}
}
