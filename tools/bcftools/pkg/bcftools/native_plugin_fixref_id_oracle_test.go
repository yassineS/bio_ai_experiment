package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"testing"
)

// Live-oracle parity tests for the fixref `id`/--use-id mode (MODE_USE_ID).
//
// The mode determines the correct REF allele from a separate dbSNP VCF keyed by
// the ID (rsID) column rather than from strand convention. Both the genuine
// upstream bcftools (BCFTOOLS_PLUGINS -> vendored fixref.so) and OUR port are
// driven through their CLIs with the SAME upstream-accepted argv; BOTH the
// corrected VCF on stdout (provenance-stripped) AND the stderr stats summary are
// compared byte-for-byte. The dbSNP fixture is bgzipped and tabix-indexed (the
// upstream synced reader requires the index for its region restriction); our
// port streams it through iohelper's transparent BGZF decoding.

// runFixrefBoth runs bin with argv (and BCFTOOLS_PLUGINS set) and returns its
// stdout and stderr separately. A non-zero exit is a hard failure.
func runFixrefBoth(t *testing.T, bin string, argv ...string) (stdout, stderr []byte) {
	t.Helper()
	cmd := exec.Command(bin, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v\nstderr: %s", bin, argv, err, errBuf.String())
	}
	return out.Bytes(), errBuf.Bytes()
}

// TestNativePluginFixrefUseID drives the upstream binary and our port with the
// same `+fixref FILE -- -f REF -i dbsnp.vcf.gz [...]` argv and asserts BOTH the
// stdout VCF (provenance-stripped) and the stderr stats summary match. The
// cases exercise: a clean match (none), an ALT->REF swap with GT flip, an
// rsID not present in dbSNP (unresolved -> "skip"), an input ID="." (unresolved
// -> "skip"), the --discard variant (unresolved sites dropped), and a custom
// -t tag name.
func TestNativePluginFixrefUseID(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "fixref_id.vcf")
	fa := fastaRefFixture(t)
	dbsnp := parityFixture(t, "fixref_dbsnp.vcf.gz")
	// The tabix index must exist next to the bgzipped dbSNP file for upstream.
	parityFixture(t, "fixref_dbsnp.vcf.gz.tbi")

	cases := [][]string{
		{"-f", fa, "-i", dbsnp},                   // implies -m id
		{"-f", fa, "-m", "id", "-i", dbsnp},       // explicit -m id
		{"-f", fa, "-i", dbsnp, "-d"},             // discard the unresolved
		{"-f", fa, "-i", dbsnp, "-t", "MYID"},     // custom tag name
		{"-f", fa, "-m", "id", "-i", dbsnp, "-d"}, // explicit + discard
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			argv := pluginCLIArgs("fixref", fixture, args)
			upOut, upErr := runFixrefBoth(t, bin, argv...)
			ourOut, ourErr := runFixrefBoth(t, ourBinPath, argv...)
			if !bytes.Equal(stripProvenanceBytes(upOut), stripProvenanceBytes(ourOut)) {
				t.Fatalf("fixref id stdout diverges (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					argv, snippet(upOut, 1200), snippet(ourOut, 1200))
			}
			if !bytes.Equal(upErr, ourErr) {
				t.Fatalf("fixref id stderr summary diverges (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					argv, snippet(upErr, 2000), snippet(ourErr, 2000))
			}
		})
	}
}
