package bcftools

import (
	"testing"
)

// Live-oracle parity tests for the LD / annotation / soft-filter / rand /
// default-maxAF / keep-sites modes of the native prune plugin — the modes that
// were previously rejected as "not supported". Each case drives BOTH the
// genuine upstream bcftools (built via buildBcftools, BCFTOOLS_PLUGINS pointed
// at the vendored .so dir) and OUR port through their CLIs with the SAME
// run()-style argv and diffs stdout byte-for-byte after stripping provenance.
// prune is run()-style: `+prune <opts> FILE`.
//
// The LD math (r2 / Lewontin's D' / Ragsdale's RD) is exercised on a
// multi-sample diploid fixture (prune_ld.vcf) so the genotype-correlation is
// non-trivial; rand selection and randomize-missing are pinned with an explicit
// --random-seed so the drand48 draw order is reproducible byte-for-byte.

// TestNativePluginPruneLD covers -m R2=/LD=/RD= thresholding (hard drop) and
// -a r2/LD/RD/count annotation across bp and site-count windows.
func TestNativePluginPruneLD(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	ld := parityFixture(t, "prune_ld.vcf")
	cases := [][]string{
		// LD thresholding (hard drop), bp + site-count windows.
		{"-m", "R2=0.5", "-w", "100bp"},
		{"-m", "R2=0.3", "-w", "1000bp"},
		{"-m", "R2=0.6", "-w", "3"},
		{"-m", "LD=0.4", "-w", "200bp"},
		{"-m", "LD=0.5", "-w", "5"},
		{"-m", "RD=0.01", "-w", "100bp"},
		{"-m", "RD=0.005", "-w", "1000bp"},
		{"-m", "0.5", "-w", "100bp"}, // bare number == r2
		// Annotation: positions + values, single and combined measures.
		{"-a", "r2", "-w", "1"},
		{"-a", "LD", "-w", "1"},
		{"-a", "RD", "-w", "1"},
		{"-a", "r2,LD,RD", "-w", "1", "-f", "."},
		{"-a", "r2,LD,RD", "-w", "1000bp"},
		{"-a", "r2", "-w", "2"},
		{"-a", "count", "-w", "100bp"},
		{"-a", "r2,LD,RD,count", "-w", "2"},
		// HD is an alias for RD.
		{"-m", "HD=0.01", "-w", "100bp"},
		{"-a", "HD", "-w", "1"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, ld, "prune", args...)
		})
	}
}

// TestNativePluginPruneSoftFilter covers -f LABEL soft filtering with -m
// thresholds (the FILTER column is set instead of the record being dropped),
// including the ##FILTER header line.
func TestNativePluginPruneSoftFilter(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	ld := parityFixture(t, "prune_ld.vcf")
	cases := [][]string{
		{"-m", "R2=0.5", "-f", "MAX_R2", "-w", "100bp"},
		{"-m", "R2=0.4", "-f", "MAX_R2", "-w", "1000bp"},
		{"-m", "LD=0.5", "-f", "MAX_LD", "-w", "200bp"},
		{"-m", "RD=0.005", "-f", "MAX_RD", "-w", "100bp"},
		{"-m", "R2=0.4", "-f", "MAX_R2", "-w", "5"}, // site-count window header
		// Annotate AND soft-filter together.
		{"-m", "0.4", "-f", "MAX_R2", "-a", "r2", "-w", "1000bp"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, ld, "prune", args...)
		})
	}
}

// TestNativePluginPruneRand covers the "rand" nsites selection mode with a
// fixed --random-seed, asserting byte-parity of the deterministic drand48 draw
// order across bp and site-count windows.
func TestNativePluginPruneRand(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	ld := parityFixture(t, "prune_ld.vcf")
	cases := [][]string{
		{"-n", "1", "-N", "rand", "--random-seed", "42", "-w", "100bp"},
		{"-n", "2", "-N", "rand", "--random-seed", "1", "-w", "1000bp"},
		{"-n", "1", "-N", "rand", "--random-seed", "7", "-w", "3"},
		{"-n", "2", "-N", "rand", "--random-seed", "12345", "-w", "5"},
		{"-n", "1", "-N", "rand", "--random-seed", "0", "-w", "200bp"}, // default seed
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, ld, "prune", args...)
		})
	}
}

// TestNativePluginPruneDefaultMaxAF covers the default maxAF selection without
// an explicit --AF-tag: the allele frequency is computed from INFO/AC+AN (or
// the genotypes) the way bcf_calc_ac does, including the upstream af=alt/ref
// quirk, and ties are resolved with the stable-sort order matching glibc qsort.
func TestNativePluginPruneDefaultMaxAF(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	ld := parityFixture(t, "prune_ld.vcf")
	plain := parityFixture(t, "prune.vcf")
	tie := parityFixture(t, "prune_tie.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{ld, []string{"-n", "1", "-w", "100bp"}},
		{ld, []string{"-n", "2", "-w", "1000bp"}},
		{ld, []string{"-n", "1", "-N", "maxAF", "-w", "1Mb"}},
		{plain, []string{"-n", "1", "-w", "100bp"}},
		{plain, []string{"-n", "2", "-w", "1000bp"}},
		{plain, []string{"-n", "1", "-i", "QUAL>40", "-w", "100bp"}},
		// All-equal and partial-tie AF: exercises the stable tie order.
		{tie, []string{"-n", "1", "-w", "1000bp"}},
		{tie, []string{"-n", "2", "-w", "1000bp"}},
		{tie, []string{"-n", "3", "-w", "1000bp"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "prune", tc.args...)
		})
	}
}

// TestNativePluginPruneKeepSites covers -k/--keep-sites with -m count= and the
// LD modes: records failing the -i/-e expression are left in place (and exempt
// from cluster counting/removal) instead of being discarded.
func TestNativePluginPruneKeepSites(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	ld := parityFixture(t, "prune_ld.vcf")
	cases := [][]string{
		{"-m", "count=2", "-k", "-i", "QUAL>40", "-w", "100bp"},
		{"-m", "count=2", "-k", "-e", "QUAL>40", "-w", "1000bp"},
		{"-m", "count=1", "-k", "-i", `INFO/AC>5`, "-w", "200bp"},
		{"-m", "R2=0.5", "-k", "-i", "QUAL>40", "-w", "100bp"},
		{"-a", "count", "-k", "-i", "QUAL>40", "-w", "100bp"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, ld, "prune", args...)
		})
	}
}

// TestNativePluginPruneRandMissing covers --randomize-missing: missing
// genotypes are replaced by a deterministic drand48 draw against the site
// allele frequency before the LD computation. A fixed --random-seed pins the
// draw order for byte-parity.
func TestNativePluginPruneRandMissing(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	miss := parityFixture(t, "prune_missing.vcf")
	cases := [][]string{
		{"-a", "r2,LD,RD", "-w", "1", "--randomize-missing", "--random-seed", "7"},
		{"-a", "r2,LD,RD", "-w", "1", "--randomize-missing", "--random-seed", "99"},
		{"-m", "R2=0.5", "-w", "100bp", "--randomize-missing", "--random-seed", "3"},
		{"-a", "r2,LD,RD", "-w", "1"}, // no rand-missing: missing GTs break the pair
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, miss, "prune", args...)
		})
	}
}
