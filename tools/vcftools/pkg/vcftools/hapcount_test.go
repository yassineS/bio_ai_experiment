package vcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestParity_Hapcount verifies --hapcount produces byte-for-byte output
// matching upstream vcftools. The fixture exercises:
//
//   - Multiple chromosomes (chr 1, chr 2) with multiple bins each.
//   - A "sentinel" chromosome (chr 3) at the end of the VCF whose only
//     role is to TRIGGER the chr 2 → chr 3 transition so that chr 2's
//     bins are emitted. This is required because upstream's final-flush
//     path is read-after-free UB and (in observed behaviour) drops the
//     last chromosome from the output. See hapcount.go for the
//     full explanation.
//   - The "wrong-bin overwrite" bug: chr 1 site at pos 1500 falls into
//     bin (1000,2000] AFTER a flush of bin (100,500] at site 5
//     (pos 550 — outside both bins). Upstream's prev_bin_idx-shift bug
//     then causes the chr-1 transition flush to write site-6's data
//     into bin 0's slots, so the output shows bin 0 with N_SNP=1
//     (site 6's count) rather than N_SNP=4 (the actual count for
//     sites 1-4).
func TestParity_Hapcount(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_fixture.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Hapcount_MissingGT verifies that missing-allele genotypes
// (`.|.`) are accepted as phased and carry through to the haplotype
// vectors as `-1` (upstream lines 1361-1366). The fixture's site at
// pos 200 has S2 = `.|.`; the two unique missing-allele haplotypes
// contribute to N_UNIQ_HAPS in the bin's row.
func TestParity_Hapcount_MissingGT(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_missing.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount_missing.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount (missing) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Hapcount_DropsUnphased verifies that the implicit
// `phased_only` filter (upstream parameters.cpp:248) drops `0/1`-style
// sites before they reach the hapcount accumulator. The fixture has
// three chr-1 sites of which the middle one is `0/1`; the output
// row for chr 1 bin 0 reflects only the two phased sites.
func TestParity_Hapcount_DropsUnphased(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_unphased.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount_unphased.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount (unphased) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHapcount_OverlappingBED rejects a BED file where two bins on the
// same chromosome overlap. Upstream variant_file_output.cpp:1212-1214
// errors with "BED file must be non-overlapping." (exit code 33). We
// surface a fatal Run-time error rather than producing garbage.
func TestHapcount_OverlappingBED(t *testing.T) {
	tmp := t.TempDir()
	bedPath := filepath.Join(tmp, "overlap.bed")
	if err := os.WriteFile(bedPath, []byte("#hdr\n1\t100\t300\n1\t200\t400\n"), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "hapcount_fixture.vcf"))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	err = Run(in, &Params{
		OutPrefix:   filepath.Join(tmp, "out"),
		HapcountBED: bedPath,
		MinAlleles:  2,
	})
	if err == nil {
		t.Fatal("expected Run() to error on overlapping BED, got nil")
	}
	if !strings.Contains(err.Error(), "non-overlapping") {
		t.Errorf("want error mentioning 'non-overlapping', got: %v", err)
	}
}

// TestHapcount_MissingBED reports a clear error when the BED file
// cannot be opened. Mirrors upstream variant_file_output.cpp:1178.
func TestHapcount_MissingBED(t *testing.T) {
	tmp := t.TempDir()
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "hapcount_fixture.vcf"))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	err = Run(in, &Params{
		OutPrefix:   filepath.Join(tmp, "out"),
		HapcountBED: filepath.Join(tmp, "does-not-exist.bed"),
		MinAlleles:  2,
	})
	if err == nil {
		t.Fatal("expected Run() to error on missing BED file, got nil")
	}
	if !strings.Contains(err.Error(), "Could not open BED file") {
		t.Errorf("want error 'Could not open BED file', got: %v", err)
	}
}

// TestHapcount_EmptyVCF_StillWritesHeader pins that the output file is
// created (header-only) even when no input variant ever touches a BED
// bin — useful as a defensive guarantee for downstream pipelines.
func TestHapcount_EmptyVCF_StillWritesHeader(t *testing.T) {
	tmp := t.TempDir()
	vcfPath := filepath.Join(tmp, "empty.vcf")
	if err := os.WriteFile(vcfPath, []byte("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	in, err := os.Open(vcfPath)
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	prefix := filepath.Join(tmp, "out")
	if err := Run(in, &Params{
		OutPrefix:   prefix,
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := readFileBytes(t, prefix+".hapcount")
	want := []byte("#CHROM\tBIN_START\tBIN_END\tN_SNP\tN_UNIQ_HAPS\tN_GROUPS\t{MULTIPLICITY:FREQ}\n")
	if !bytes.Equal(got, want) {
		t.Errorf("header-only output mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParams_TempDirAccepted is a smoke check that --temp is accepted
// without altering output. It mirrors upstream parameters.cpp:341
// (DIR is stored as the spill-file base path); we don't spill so the
// flag is observably a no-op.
func TestParams_TempDirAccepted(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvBurden: true,
		TempDir:    "/some/path/we/never/touch",
		MinAlleles: 2,
	})
	// The TempDir field is silent — assert the iburden output still
	// matches the no-temp baseline.
	got := readFileBytes(t, prefix+".iburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden.expected.iburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".iburden (with --temp) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}
