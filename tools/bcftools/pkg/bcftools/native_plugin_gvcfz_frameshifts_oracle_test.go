package bcftools

import (
	"os"
	"path/filepath"
	"testing"
)

// CLI-to-CLI live-oracle parity tests for the natively-ported gvcfz and
// frameshifts plugins. Both run the genuine upstream bcftools 1.23.1 binary
// (BCFTOOLS_PLUGINS pointed at the vendored .so directory) and OUR port through
// their CLIs with the SAME upstream-accepted argv, comparing stdout byte-for-byte
// after stripping provenance. These cases run (not skip) whenever the upstream
// submodule + plugin .so directory are present.

// TestNativePluginGvcfz exercises the gvcfz block-resizing state machine against
// upstream: the multi-group example from the man page, the single-block alt
// grouping (with and without -a), a custom non-PASS FILTER label, the RGQ
// gq-key path, and the -i record-level pre-filter.
func TestNativePluginGvcfz(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gvcfz := parityFixture(t, "gvcfz.vcf")
	rgq := parityFixture(t, "gvcfz_rgq.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{gvcfz, []string{"-g", `PASS:GQ>60 & DP<20; PASS:GQ>40 & DP<15; Flt1:GQ>20; Flt2:-`}},
		{gvcfz, []string{"-g", `PASS:GT!="alt"`}},
		{gvcfz, []string{"-a", "-g", `PASS:GT!="alt"`}},
		{gvcfz, []string{"-g", `Lo:GQ<60; Hi:-`}},
		{gvcfz, []string{"-g", `Flt1:GT!="alt"; PASS:GQ>1`}},
		{gvcfz, []string{"-i", "GQ>40", "-g", `PASS:-`}},
		{gvcfz, []string{"-e", "GQ>40", "-g", `PASS:-`}},
		{rgq, []string{"-g", `PASS:RGQ>10`}},
		{rgq, []string{"-g", `Lo:RGQ<30; Hi:-`}},
	}
	for _, tc := range cases {
		t.Run(filepath.Base(tc.fixture)+" "+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "gvcfz", tc.args...)
		})
	}
}

// TestNativePluginFrameshifts exercises the frameshifts exon-overlap + OOF
// annotation against upstream, with both a plain BED and a bgzipped+tabixed BED
// exon file (the two synced-reader region-cursor paths). The corrected
// in-frame/out-of-frame computation (--fix-oof) is OUR-only and is covered by
// the binary-free unit tests, not the oracle.
func TestNativePluginFrameshifts(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "frameshifts.vcf")

	// Plain BED exon file (in-memory regs path).
	plainBED := parityFixture(t, "frameshifts_exons.bed")

	// bgzipped + tabixed BED exon file (tabix-index cursor path). Built in a temp
	// dir from the same source so the committed fixture stays a plain BED. The
	// htslib bgzip/tabix binaries come from the vendored reference_code/htslib
	// directory directly (no need to rebuild bcftools from source, which the
	// worktree's unpopulated submodule cannot do).
	htslibDir := vendoredHtslibDir(t)
	tmp := t.TempDir()
	src, rerr := os.ReadFile(plainBED)
	if rerr != nil {
		t.Fatalf("read exon bed: %v", rerr)
	}
	bedCopy := filepath.Join(tmp, "exons.bed")
	if werr := os.WriteFile(bedCopy, src, 0o644); werr != nil {
		t.Fatalf("write exon bed: %v", werr)
	}
	gzBED := bgzipIndex(t, htslibDir, bedCopy, "-p", "bed")

	for _, exon := range []string{plainBED, gzBED} {
		exon := exon
		t.Run(filepath.Base(exon), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "frameshifts", "-e", exon)
		})
	}
}

// vendoredHtslibDir returns the absolute path of the vendored reference_code/htslib
// directory holding the bgzip/tabix binaries the oracle uses to build a
// bgzipped+tabixed exon file. It fails (not skips) when bgzip/tabix are missing.
func vendoredHtslibDir(t *testing.T) string {
	t.Helper()
	dir, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "htslib"))
	if err != nil {
		t.Fatalf("htslib dir abs: %v", err)
	}
	for _, bin := range []string{"bgzip", "tabix"} {
		if _, serr := os.Stat(filepath.Join(dir, bin)); serr != nil {
			t.Fatalf("vendored htslib %s missing: %v", bin, serr)
		}
	}
	return dir
}
