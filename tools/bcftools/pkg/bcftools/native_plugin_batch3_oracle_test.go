package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the batch-3 FASTA-backed native plugins
// (fill-from-fasta, fixref) and the unsupported vrfs entry. Each builds the
// genuine upstream bcftools via buildBcftools() and drives BOTH that binary and
// OUR port through their CLIs with the SAME upstream-accepted argv, diffing
// stdout byte-for-byte after stripping provenance — exactly the batch-2 harness
// (pluginCLIArgs / runUpstreamPlugin / runOursPlugin / assertPluginParity /
// stripProvenanceBytes), kept strictly CLI-to-CLI.
//
// Both fill-from-fasta and fixref are generic init/process plugins (their
// upstream .c exports init/process/destroy, not run), so the invocation form is
// `+name FILE -- <plugin opts>`; pluginCLIArgs builds that automatically. The
// FASTA reference is the testdata/parity/fixref_ref.fa fixture (two contigs,
// indexed by a sibling .fai); the VCF fixtures carry REF=N, strand-swapped, and
// ambiguous-pair (A/T, C/G) sites so the fixref flip/swap/top conversions are
// genuinely exercised.

// fastaRefFixture returns the absolute path to the shared two-contig FASTA
// reference used by the batch-3 plugins.
func fastaRefFixture(t *testing.T) string {
	t.Helper()
	return parityFixture(t, "fixref_ref.fa")
}

func TestNativePluginFillFromFasta(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "fillfromfasta.vcf")
	fa := fastaRefFixture(t)

	// A header-lines file declaring a fresh INFO tag, exercised via -h.
	hdrFile := filepath.Join(t.TempDir(), "p3.hdr")
	if err := os.WriteFile(hdrFile, []byte("##INFO=<ID=P3,Number=1,Type=String,Description=\"mask\">\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	cases := [][]string{
		{"-c", "REF", "-f", fa},               // fix the REF allele from the FASTA
		{"-c", "REF", "-f", fa, "-N"},         // with non-ACGTN replacement
		{"-c", "AA", "-f", fa},                // fill a declared String INFO tag
		{"-c", "RN", "-f", fa},                // fill a declared Integer INFO tag (single base)
		{"-c", "P3", "-f", fa, "-h", hdrFile}, // append header lines then fill
		{"-c", "INFO/AA", "-f", fa},           // INFO/ prefix on the column name
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "fill-from-fasta", args...)
		})
	}
}

func TestNativePluginFixref(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "fixref.vcf")
	fa := fastaRefFixture(t)

	// stats mode produces no stdout (output suppressed); the conversion modes
	// emit a VCF whose REF/ALT/GT and INFO/FIXREF annotation are parity-checked.
	cases := [][]string{
		{"-f", fa},                                 // default: stats (no stdout)
		{"-f", fa, "-m", "stats"},                  // explicit stats
		{"-f", fa, "-m", "ref-alt"},                // swap/flip REF/ALT, leave GTs
		{"-f", fa, "-m", "swap"},                   // swap only, leave GTs
		{"-f", fa, "-m", "flip"},                   // flip/swap + GT for non-ambiguous
		{"-f", fa, "-m", "flip", "-d"},             // discard the unresolved
		{"-f", fa, "-m", "flip-all"},               // flip/swap all SNPs
		{"-f", fa, "-m", "top"},                    // TOP -> fwd with ambiguous walking
		{"-f", fa, "-m", "top", "-d"},              // top + discard
		{"-f", fa, "-m", "ref-alt", "-t", "MYFIX"}, // custom tag name
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "fixref", args...)
		})
	}
}

// TestNativePluginBatch3Unsupported asserts the deliberately unsupported path
// (vrfs) and the fixref id-mode misuse (-m id with no -i file) fail cleanly from
// Init rather than diverging. The fixref id mode itself is supported and is
// parity-checked by TestNativePluginFixrefUseID.
func TestNativePluginBatch3Unsupported(t *testing.T) {
	fixture := parityFixture(t, "fixref.vcf")
	fa := fastaRefFixture(t)
	cases := []struct {
		name string
		args []string
	}{
		{"fixref", []string{"-f", fa, "-m", "id"}}, // -m id without -i is an error
		{"vrfs", []string{"-f", fa, "-a", "bams.txt", "-s", "sites.txt"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+joinArgs(tc.args), func(t *testing.T) {
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         tc.name,
				Args:         tc.args,
				InputFile:    fixture,
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected an unsupported error for +%s %v, got nil", tc.name, tc.args)
			}
		})
	}
}
